"""
This module contains functions to compute the Power Spectral Density (PSD) of surface pressure data from an HDF5 file using the Welch method.
The computation is parallelized over blocks of nodes to improve performance. 
The main function `PSD_surface_data` prepares the input and segments the data into blocks for parallel processing, 
while the worker function `compute_PSD_block` handles the PSD computation for each block of nodes.
"""
import scipy
import scipy.signal as signal
from scipy.signal import welch, coherence
import multiprocessing
from multiprocessing import Pool, cpu_count, Value, Lock

def PSD_surface_data(surface_pressure_data: str, var: str, dt: float, reload: bool = True, block_size: int = 1000, nchunk: int = 3):
    """
    This is the main function that computes the Power Spectral Density (PSD) of surface pressure data from an HDF5 file via the Welch method. 
    It parallelizes the computation over blocks of nodes to improve performance.
    
    Calculates the PSD of surface pressure fluctuations for all nodes on an airfoil grid.
    The PSD is computed via Welch's method and parallelized over blocks of nodes.
    
    Args:
        surface_pressure_data (str): Path to the HDF5 file containing surface pressure data.
        var (str): Variable name in the HDF5 file to extract.
        dt (float): Time-step between successive samples.
        reload (bool): If True, forces recalculation even if the output file exists.
        block_size (int): Number of nodes per block for multiprocessing.
        nchunk (int): Number of chunks for dividing the time-series data when computing FFT.
        
    Returns:
        str: Path to the output HDF5 file containing the computed PSD data.
    """
    header = "Performing PSD on surface pressure data"
    print(f"\n{header:.^80}\n")
    
    psd_file_path = surface_pressure_data.replace('.hdf5', '_psd.hdf5')
    if os.path.exists(psd_file_path) and not reload:
        print(f"{psd_file_path} already exists. Skipping PSD computation.")
    else:
        print("---->Computing PSD for surface pressure data.")
        # Load the pressure data from the HDF5 file.
        with h5py.File(surface_pressure_data, 'r') as h5f:
            data = h5f[var][:]
            # Ensure time is the first dimension.
            if np.shape(data)[0] > np.shape(data)[1]:
                data = np.swapaxes(data, 0, 1)
                
        nt, nx = data.shape
        print("The surface data is %d (nodes) x %d (timesteps)" % (nx, nt))
        
        # Determine number of blocks.
        nblocks = int(np.ceil(nx / block_size))
        num_processes = multiprocessing.cpu_count()
        print("Processing %d blocks of up to %d nodes each using %d cores." % (nblocks, block_size, num_processes))
        
        # Create list of arguments for each block.
        block_args = []
        for b in range(nblocks):
            start = b * block_size
            end = min((b + 1) * block_size, nx)
            block_data = data[:, start:end]
            block_args.append((block_data, dt, nchunk, b))
        
        # Process each block in parallel.
        with Pool(processes=num_processes) as pool:
            results = pool.map(compute_PSD_block, block_args)
        
        # Extract frequency (assumed identical across blocks) and combine PSD blocks.
        freq = results[0][1]
        # Concatenate along axis 0 so that the final shape is (n_nodes, n_freq).
        psd_all = np.concatenate([res[0] for res in results], axis=0)
        
        # Save the computed frequency and PSD data.
        print("Saving PSD data to file:", psd_file_path)
        with h5py.File(psd_file_path, 'w') as h5f_out:
            h5f_out.create_dataset('frequency', data=freq)
            h5f_out.create_dataset('pressure_psd', data=psd_all)
            h5f_out.attrs['Date Computed'] = datetime.now().strftime("%Y-%m-%d")
            h5f_out.attrs['Reference Signal'] = var
            h5f_out.attrs['Reference Signal Path'] = surface_pressure_data
            h5f_out.attrs['Reference Signal Sampling Frequency'] = 1/dt
            h5f_out.attrs['Nodes'] = nx
            h5f_out.attrs['Time Steps'] = nt
    
    header = "PSD Computation Complete"
    print(f"\n{header:.^80}\n")
    return psd_file_path


def compute_PSD_block(args):
    """
    This is the WORKER function that computes the PSD for a block of nodes. It is designed to be called in parallel using multiprocessing.
    Computes the PSD for a block of nodes for the specific worker process. This function is called by the multiprocessing Pool.
    The PSD is computed using Welch's method for each node in the block over the time-series data.
    
    args: tuple containing (block_data, dt, nchunk, block_number)
        - block_data: 2D numpy array of shape (nt, n_nodes_in_block)
        - dt: time-step (float)
        - nchunk: number of chunks to divide the data for the FFT
        - block_number: integer identifier for the current block
          
    Returns:
        - psd_block  (n_nodes_in_block, n_freq): 2D numpy array of PSD values with shape
        - freq (identical for all nodes in the block): 1D numpy array of frequency values 
    """
    block_data, dt, nchunk, block_number = args
    if block_number % 100 == 0:  
        print(f"    Computing PSD for block {block_number}")
        
    n_nodes = block_data.shape[1]
    n_time = block_data.shape[0]
    psd_list = []
    freq = None
    
    for i in range(n_nodes):
        lensg = n_time
        nperseg = int(lensg / nchunk)
        nfft = next_greater_power_of_2(nperseg)   
        # Compute the PSD for the i-th node time-series using Welch's method.
        f, psd = welch(block_data[:, i],fs=1.0/dt,window='hann',nperseg=nperseg,nfft=nfft,scaling='density')
        if freq is None:
            freq = f
        # Store PSD as a row (node) in the list.
        psd_list.append(psd)
        
    # Create an array with shape (n_nodes_in_block, n_freq)
    psd_block = np.array(psd_list)
    return psd_block, freq


def source_psd(output_path: str,mesh_fileDir: str,mesh_fileName: str, surface_mesh: str, data: str, data_psd: str,
    freq_select: list = (500, 2000),
    df: float = 50.0,                  # <-- band width in Hz
    band_stat: str = "mean",           # "mean" (band-average) or "sum" (band-integrated)
    pref: float = 2e-5,                # acoustic reference pressure [Pa]
    )-> None:
    """
    This function outputs the surface pressure PSD data for visualization in ParaView. 
    It reads the surface mesh and the corresponding pressure data (mean, RMS, min) and PSD data from HDF5 files, 
    computes band-averaged or band-summed PSD values for selected frequencies, and saves the results in a new HDF5 file.
    
    Args: 
        output_path (str): Directory where the output HDF5 file will be saved.
        mesh_fileDir (str): Directory containing the surface mesh file.
        mesh_fileName (str): Name of the surface mesh file.
        surface_mesh (str): Path to the surface mesh file.
        data (str): Path to the HDF5 file containing pressure statistics (mean, RMS, min).
        data_psd (str): Path to the HDF5 file containing pressure PSD data.
        freq_select (list): List of frequencies for which to compute band-averaged or band-summed PSD values.
        df (float): Bandwidth in Hz for averaging/summing around each selected frequency.
        band_stat (str): Method for aggregating PSD values within the band; either "mean" or "sum".
        pref (float): Acoustic reference pressure in Pascals for converting PSD to dB.
    """
    
    
    base_mesh = os.path.join(mesh_fileDir, mesh_fileName)
    text = 'Performing Surface Source Localization based on PSD'
    print(f'\n{text:.^80}\n')

    # Loading the surface mesh
    print('----> Loading the Airfoil Surface Mesh')
    reader = Reader('hdf_antares')
    reader['filename'] = surface_mesh
    base = reader.read()
    base.show()

    # Extract mesh coordinates
    x, y, z = base[0][0]["x"], base[0][0]["y"], base[0][0]["z"]
    num_nodes = len(x)

    # Loading the Surface pressure statistics
    print(f'\n----> Loading the Pressure Data: {data}')
    with h5py.File(data, 'r') as h5f:
        p_mean = h5f['mean_pressure'][:]
        p_rms  = h5f['rms_pressure'][:]
        p_min  = h5f['min_pressure'][:]

    # Loading the PSD data
    print(f'\n----> Loading the PSD of the Pressure Data: {data_psd}')
    with h5py.File(data_psd, 'r') as h5f:
        p_hat = h5f['pressure_psd'][:]   # shape: (num_nodes, nfreq)
        freq  = h5f['frequency'][:]      # shape: (nfreq,)

    assert p_hat.shape[0] == num_nodes, 'The number of nodes in the pressure data and the mesh do not match'
    print('The pressure data is %d (nodes) x %d (frequency bins)' % (p_hat.shape[0], p_hat.shape[1]))

    # Helper: band indices around f0 using df
    def _band_indices(freq_array: np.ndarray, f0: float, df_hz: float):
        half = 0.5 * df_hz
        lo, hi = f0 - half, f0 + half
        idx = np.where((freq_array >= lo) & (freq_array <= hi))[0]
        if idx.size == 0:
            # fallback: at least include the nearest bin
            idx = np.array([int(np.argmin(np.abs(freq_array - f0)))], dtype=int)
        return idx

    # Create output base
    print('\n----> Saving the surface fft/psd data for visualization')
    animated_base = Base()
    animated_base['0'] = Zone()
    animated_base[0].shared["x"], animated_base[0].shared["y"], animated_base[0].shared["z"] = x, y, z
    animated_base[0].shared.connectivity = base[0][0].connectivity
    animated_base[0][str(0)] = Instant()

    # Save band-averaged PSD for selected frequencies
    for f0 in freq_select:
        idx_band = _band_indices(freq, float(f0), float(df))

        # IMPORTANT: average/sum in linear units
        # p_hat is assumed real-valued PSD already; keep real part for safety
        Spp_band = np.real(p_hat[:, idx_band])

        if band_stat.lower() == "mean":
            Spp_agg = np.mean(Spp_band, axis=1)
            tag = f"df{int(df)}Hz_mean"
        elif band_stat.lower() == "sum":
            Spp_agg = np.sum(Spp_band, axis=1)
            tag = f"df{int(df)}Hz_sum"
        else:
            raise ValueError("band_stat must be 'mean' or 'sum'")

        # Convert to dB re pref^2 (guard against zeros/negatives)
        Spp_agg_safe = np.maximum(Spp_agg, np.finfo(float).tiny)
        dB_agg = 10.0 * np.log10(Spp_agg_safe / (pref ** 2))

        # Store results in both dB and linear
        f0_int = int(round(float(f0)))
        animated_base[0][str(0)][f'frequency_{f0_int}_Hz_P_band_{tag}_dB']  = dB_agg
        animated_base[0][str(0)][f'frequency_{f0_int}_Hz_P_band_{tag}_Spp'] = Spp_agg
    # Store pressure statistics
    animated_base[0][str(0)]['Pressure_Mean'] = p_mean
    animated_base[0][str(0)]['Pressure_RMS']  = p_rms
    animated_base[0][str(0)]['Pressure_Min']  = p_min
    animated_base[0][str(0)].attrs['band_stat'] = band_stat
    animated_base[0][str(0)].attrs['band_width_Hz'] = df

    # Write output
    w = Writer('hdf_antares')
    w['filename'] = os.path.join(output_path, 'Surface_psd')
    w['base'] = animated_base
    w['dtype'] = 'float32'
    w.dump()
    del animated_base

    print(f'\n----> Saving the output surface psd data as: {os.path.join(output_path, "Surface_psd.hdf5")}')
    text = 'Surface PSD Complete!'
    print(f'\n{text:.^80}\n')
