"""
This module contains functions to compute the Cross Spectral Density (CSD) of surface pressure data and the far-field FW-H data from an HDF5 file using the Welch method.
The computation is parallelized over blocks of nodes to improve performance. 
The main function `CSD_surface_data` the input and segments the data into blocks for parallel processing, 
while the worker function `compute_CSD_block` handles the CSD computation for each block of nodes.
"""
import os
from datetime import datetime
import numpy as np
import h5py
import scipy
import scipy.signal as signal
from scipy.signal import welch, coherence
import multiprocessing
from multiprocessing import Pool, cpu_count, Value, Lock
from antares import Reader, Base, Zone, Instant, Writer
from .utils import next_greater_power_of_2


def CSD_surface_data(farfield_dir:str, mic_file:str, mic_num:int, surface_pressure_data: str, var: str, dt: float, reload: bool = True, block_size: int = 1000, nchunk: int = 3):
    """
    The CSD_surface_data function computes the Cross Spectral Density (CSD) and Magnitude Squared Coherence between surface pressure data and a reference far-field microphone signal. 
    The computation is performed in blocks to optimize memory usage and speed, and it can be parallelized across multiple CPU cores.
    
    Args: 
        farfield_dir (str): parent path to the HDF5 file containing the far-field microphone data.
        mic_file (str): Name of the HDF5 file containing the far-field microphone data.
        mic_num (int): The microphone number to be used as the reference signal.
        surface_pressure_data (str): Path to the HDF5 file containing the surface pressure data.
        var (str): The variable name in the HDF5 file for surface pressure data.
        dt (float): Time-step for the data.
        reload (bool, optional): If True, recompute CSD even if output file exists. Defaults to True.
        block_size (int, optional): Number of nodes to process in each block. Defaults to 1000.
        nchunk (int, optional): Number of chunks for FFT computation. Defaults to 3.
    
    Returns:
        csd_file_path (str): Path to the HDF5 file where the computed CSD and coherence data are saved.
    """
    
    header = "Performing Magnitude Squared Coherence on surface pressure data"
    print(f"\n{header:.^80}\n")
    # The sampling frequency
    fs = 1/dt
    # Computing the CSD in blocks
    csd_file_path = surface_pressure_data.replace('.hdf5', '_csd.hdf5')
    farfield_data = os.path.join(farfield_dir, mic_file)
    
    if os.path.exists(csd_file_path) and not reload:
        print("----> CSD file already exist, Reload is set to False.")
        print(f"      CSD File Path: {csd_file_path} ")
        with h5py.File(csd_file_path, 'r') as h5f:
            date_extracted = h5f.attrs.get('Date Computed', 'Unknown')
        print("      The CSD data was previously computed on: {0:s}".format(date_extracted))
    else:
        print("----> Computing CSD for surface pressure data.")
        # loading the the reference pressure signal
        print("      The reference signal is loaded from \n        {0:s}".format(farfield_data))
        print("      The reference microphone number is mic {0:d}".format(mic_num))
        with h5py.File(farfield_data, 'r') as h5f:
            p_ref = h5f['Microphone_data/mic_{}/mic_{}_pressure'.format(mic_num, mic_num)][:].flatten()
        # Load the pressure data from the HDF5 file.
        with h5py.File(surface_pressure_data, 'r') as h5f:
            data = h5f[var][:]
            # Ensure time is the first dimension.
            if np.shape(data)[0] > np.shape(data)[1]:
                data = np.swapaxes(data, 0, 1)
                
        nt, nx = data.shape
        print("      The surface data is %d (nodes) x %d (timesteps)" % (nx, nt))
        print("      The reference data is %d (timesteps)" % (len(p_ref)))
        n_ref = p_ref.shape[0]
        # Align lengths by truncation
        n_timesteps = min(nt, n_ref)
        if n_timesteps != nt or n_timesteps != n_ref:
            data  = data[:n_timesteps, :]
            p_ref = p_ref[:n_timesteps]
            print(f"        After alignment: {n_timesteps} timesteps each")
    
        # Determine number of blocks.
        nblocks = int(np.ceil(nx / block_size))
        num_processes = multiprocessing.cpu_count()
        print("\n----> Processing %d blocks of up to %d nodes each using %d cores." % (nblocks, block_size, num_processes))
        
        # Create list of arguments for each block.
        block_args = []
        for b in range(nblocks):
            start = b * block_size
            end = min((b + 1) * block_size, nx)
            block_data = data[:, start:end]
            block_args.append((block_data, p_ref, dt, nchunk, b))
        
        # Process each block in parallel.
        with Pool(processes=num_processes) as pool:
            results = pool.map(compute_CSD_block, block_args)
        
        # Extract frequency (assumed identical across blocks) and combine PSD blocks.
        f_csd = results[0][2]
        f_coh = results[0][3]
        print("\n      The csd frequency shape is: ", f_csd.shape)
        print("      The coherence frequency shape is: ", f_coh.shape)
        # Concatenate along axis 0 so that the final shape is (n_nodes, n_freq).
        csd_all = np.concatenate([res[0] for res in results], axis=0)
        coh_all = np.concatenate([res[1] for res in results], axis=0)
        
        # Save the computed frequency and PSD data.
        print("      Saving CSD data to file:", csd_file_path)
        with h5py.File(csd_file_path, 'w') as h5f_out:
            h5f_out.create_dataset('frequency_csd', data=f_csd)
            h5f_out.create_dataset('pressure_csd', data=csd_all)
            h5f_out.create_dataset('frequency_coh', data=f_coh)
            h5f_out.create_dataset('pressure_coh', data=coh_all)
            h5f_out.attrs['Date Computed'] = datetime.now().strftime("%Y-%m-%d")
            h5f_out.attrs['Reference Microphone'] = mic_num
            h5f_out.attrs['Reference Signal'] = 'Microphone_data/mic_{mic_num}/mic_{mic_num}_pressure'
            h5f_out.attrs['Reference Signal Path'] = farfield_data
            h5f_out.attrs['Reference Signal Sampling Frequency'] = fs
            h5f_out.attrs['Date Computed'] = datetime.now().strftime("%Y-%m-%d")
            h5f_out.attrs['Nodes'] = nx
            h5f_out.attrs['Time Steps'] = nt
    
    header = "       CSD Computation Complete"
    print(f"\n{header:.^80}\n")
    return csd_file_path


def compute_CSD_block(args):
    """
    This is the worker function that computes the Cross Spectral Density (CSD) and Magnitude Squared Coherence between 
    surface pressure data and a reference far-field microphone signal for a block of nodes.
    It uses Welch's method to compute the CSD and coherence for each node in the block, and returns the results along with the corresponding frequency arrays.
    
    Parallelization is achieved by dividing the surface pressure data into blocks, and each block is processed independently in a separate process, similar to the PSD computation. 
    
    Parameters:
        args: tuple containing (block_data, p_ref, dt, nchunk, block_number)
            - block_data: 2D numpy array of shape (nt, n_nodes_in_block)
            - p_ref: 1D numpy array of reference pressure data
            - dt: time-step (float)
            - nchunk: number of chunks to divide the data for the FFT
            - block_number: integer identifier for the current block

    Returns:
        csd_block: 2D numpy array of shape (n_nodes_in_block, n_freq) containing the CSD values
        coh_block: 2D numpy array of shape (n_nodes_in_block, n_freq) containing the coherence values
        f_csd: 1D numpy array of frequencies corresponding to the CSD values
        f_coh: 1D numpy array of frequencies corresponding to the coherence values
    """
    block_data, p_ref, dt, nchunk, block_number = args
    if block_number % 100 == 0:  
        print(f"      Computing CSD for block {block_number}")
        
    n_nodes = block_data.shape[1]
    n_time = block_data.shape[0]
    csd_list = []
    coh_list = []
    f_csd, f_coh = None, None
    
    for i in range(n_nodes):
        lensg = n_time
        nperseg = int(lensg / nchunk)
        nfft = next_greater_power_of_2(nperseg)   
        # Compute the CSDfor the i-th node time-series using Welch's method.
        fcsd, csd = signal.csd(block_data[:, i],p_ref,fs=1.0/dt,window='hann',nperseg=nperseg,nfft=nfft)
        fcoh, coh = signal.coherence(block_data[:, i],p_ref, fs=1.0/dt,window='hann', nperseg=nperseg, nfft=nfft)
        csd = 10*np.log10(np.abs(csd)/2e-5**2)
        
        if f_csd is None:
            f_csd = fcsd
        if f_coh is None:
            f_coh = fcoh
        # Store CSD as a row (node) in the list.
        csd_list.append(csd)
        coh_list.append(coh)
        
    # Create an array with shape (n_nodes_in_block, n_freq)
    csd_block = np.array(csd_list)
    coh_block = np.array(coh_list)
    return csd_block, coh_block, f_csd, f_coh


def source_csd(output_path: str, surface_mesh: str, data: str,data_csd: str,
            freq_select: list = [500, 2000], 
            half_window: int = 0) -> None:
    """
    This function outputs the surface pressure CSD data for visualization in ParaView. 
    It reads the surface mesh and the corresponding pressure data (mean, RMS, min) and PSD data from HDF5 files, 
    computes band-averaged or band-summed CSD values for selected frequencies, and saves the results in a new HDF5 file.
    
    Args: 
        output_path (str): Directory where the output HDF5 file will be saved.
        mesh_fileDir (str): Directory containing the surface mesh file.
        mesh_fileName (str): Name of the surface mesh file.
        surface_mesh (str): Path to the surface mesh file.
        data (str): Path to the HDF5 file containing time-domain pressure statistics.
        data_csd (str): Path to the HDF5 file containing CSD and coherence spectra.
        freq_select (list, optional): List of frequencies for which to compute band-averaged CSD values. Defaults to [500, 2000].
        half_window (int, optional): Half-width of the frequency window for averaging. Defaults to 0 (no averaging).
    """
    # ——— Load surface mesh geometry ———
    print(f"\n{'Performing Surface Source Localization based on CSD':.^80}\n")
    print('----> Loading the Airfoil Surface Mesh')
    reader = Reader('hdf_antares')
    reader['filename'] = surface_mesh
    base = reader.read()
    base.show()
    x, y, z = base[0][0]["x"], base[0][0]["y"], base[0][0]["z"]
    num_nodes = len(x)

    # ——— Load time-domain pressure stats ———
    print(f'\n----> Loading the Pressure Data: {data}')
    with h5py.File(data, 'r') as h5f:
        p_mean = h5f['mean_pressure'][:]
        p_rms  = h5f['rms_pressure'][:]
        p_min  = h5f['min_pressure'][:]

    # ——— Load CSD & coherence spectra ———
    print(f'\n----> Loading the CSD of the Pressure Data: {data_csd}')
    with h5py.File(data_csd, 'r') as h5f:
        p_hat     = h5f['pressure_csd'][:]
        p_hat_coh = h5f['pressure_coh'][:]
        freq      = h5f['frequency_csd'][:]
        freq_coh  = h5f['frequency_coh'][:]

    assert p_hat.shape[0] == num_nodes, \
        'Node count mismatch between mesh and pressure data'
    print(f"      The pressure data is {p_hat.shape[0]} nodes × {p_hat.shape[1]} frequency bins")

    # ——— Prepare animated base ———
    print('\n----> Saving the surface FFT data for visualization')
    animated_base = Base()
    animated_base['0'] = Zone()
    z0 = animated_base[0]
    z0.shared["x"], z0.shared["y"], z0.shared["z"] = x, y, z
    z0.shared.connectivity = base[0][0].connectivity
    z0[str(0)] = Instant()
    record = z0[str(0)]

    # ——— Inline frequency-slicing & storage ———
    nt, n_f  = p_hat.shape
    _, n_fc = p_hat_coh.shape

    for f0 in freq_select:
        # find nearest-bin indices
        i_orig = np.argmin(np.abs(freq     - f0))
        i_coh  = np.argmin(np.abs(freq_coh - f0))

        # compute/clamp window bounds
        i0, i1 = max(0, i_orig-half_window), min(n_f,   i_orig+half_window+1)
        j0, j1 = max(0, i_coh -half_window), min(n_fc,  i_coh +half_window+1)

        # average over the small band
        band_orig = p_hat[:, i0:i1].mean(axis=1)
        band_coh  = p_hat_coh[:, j0:j1].mean(axis=1)

        label = f"{int(round(f0)):02d}"
        record[f"frequency_{label}_Hz_csd"] = band_orig
        record[f"frequency_{label}_Hz_coh"] = band_coh

    # ——— Write out full dataset ———
    record['Pressure_Mean'] = p_mean
    record['Pressure_RMS']  = p_rms
    record['Pressure_Min']  = p_min

    writer = Writer('hdf_antares')
    writer['filename'] = os.path.join(output_path, 'Surface_csd')
    writer['base']     = animated_base
    writer['dtype']    = 'float32'
    writer.dump()
    del animated_base  # free memory

    outfile = os.path.join(output_path, 'Surface_csd.hdf5')
    print(f'\n----> Saving the output surface psd data as: {outfile}')
    print(f"\n{'Surface CSD Complete!':.^80}\n")
