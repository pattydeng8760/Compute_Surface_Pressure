"""


"""

import os
import glob
import shutil
from datetime import datetime
import numpy as np
import h5py
import scipy
from scipy.interpolate import griddata, interp1d
from scipy.ndimage import gaussian_filter
from scipy.spatial import KDTree, cKDTree
from antares import Reader, Writer, Base, Zone, Instant, Treatment, Family
from .utils import replace_zeros_vectorized


def extract_files(sol_dir:str, working_dir:str, option:int=1, nskip:int=1, max_file:int=5000, reload:bool=False):
    """ Copying the surface pressure datafrom the AVBP surface FWH solution directory to the working directory, 
        Required to copy locally to avoid I/O issues and corruption of the data inside the source directory
    Args:
        sol_dir (str): path to the AVBP surface FWH solution directory
        working_dir (str): path to the working directory
        option (int, optional): extract options. 1 = sequential, 2 = skip every n files. Defaults to 1.
        nskip (int, optional): skip option. Defaults to 1.
        max_file (int, optional): maximum number of files to extract. Defaults to 5000.
        reload (bool, optional): if True, reload the files even if they already exist. Defaults to False.
    Returns:
        ntime (int): number of time steps extracted
        output (str): path to the extracted files in the working directory
    """
    text = 'Extracting the Transient Files'
    print(f'\n{text:.^80}\n')  
    output = os.path.join(working_dir,'FWH_Data')
    # An array with all the files in the source directory
    arr,arr_dest = [], []
    arr += sorted(glob.glob('{0:s}/{1:s}*.h5'.format(sol_dir,'FWH_Airfoil_')))
    print('----> Beginning the extraction of the transient files from the source directory')
    if os.path.exists(output):
        print('      There are {0:d} files in the source directory'.format(len(arr)))
        arr_dest += sorted(glob.glob('{0:s}/*.h5'.format(output)))
        print('      There are {0:d} files in the destination directory'.format(len(arr_dest)))
        print('      After extraction there should be {0:d} files in the destination directory'.format(int(np.floor(len(arr))/nskip)))
    # Checking if all the files from the source directory are already extracted
    if int(len(arr_dest)) == int(np.floor(len(arr))/nskip) or (reload == False and os.path.exists(output)):
        text = '\n----> All the files are already extracted to destination directory, Extraction is skipped'
        print(f'{text}')
    else:
        # Removing the directory if exists
        if os.path.exists(output):
            shutil.rmtree(output, ignore_errors=True)
        os.makedirs(output, exist_ok = False)
        if option == 1:
            # Only take n files
            arr = os.listdir(sol_dir)
            arr = list(arr)
            arr.sort()
            nameonly = np.array([])
            for i in range(0,np.min([int(max_file),int(len(arr))])):
                source = os.path.join(sol_dir,arr[i])
                destination = output
                shutil.copy(source,destination)
                print('    Copying file: %s' %os.path.split(source)[1]) if np.mod(i,100) == 0 else None
                print('    Copying file: %s' %os.path.split(source)[1]) if i == np.min([int(max_file),int(len(arr))])-1 else None
            arr2 = os.listdir(output)
            arr2 = list(arr2)
            arr2.sort()
        elif option == 2:
            # Skip every n files
            arr = os.listdir(sol_dir)
            arr = list(arr)
            arr.sort()
            nameonly = np.array([]) 
            for i in range(0,int(len(arr)/nskip)):
                source = os.path.join(sol_dir,arr[int(i*nskip-1)])
                destination = output
                shutil.copy(source,destination)
                print('    Copying file %s' %os.path.split(source)[1]) if np.mod(i,100) == 0 else None
    list_files=[]
    list_files+=sorted(glob.glob('{0:s}/*.h5'.format(output)))
    ntime = np.shape(list_files)[0]
    print('\n----> The files are copied to: %s' %output)
    text = 'File Extraction Complete'
    print(f'\n{text:.^80}\n')  
    output = os.path.abspath(output)
    return ntime, output


def extract_surface(mesh_file:str,input_surface:list, working_dir:str,reload:bool=False):
    """ Extracting the airfoil surface mesh from the main LES mesh
    Args:
        mesh_file (str): path to the full mesh file
        input_surface (list): list of surface names to extract from the full mesh
        working_dir (str): current working directory
        reload (bool, optional): if True, reload the surface mesh even if it already exists. Defaults to False.
    Returns:
        airfoil_mesh (str): path to the extracted airfoil surface mesh
        nodes (int): number of nodes on the airfoil surface mesh
    """
    text = 'Extracting Surface Mesh'
    print(f'\n{text:.^80}\n')
    airfoil_mesh = os.path.join(working_dir,'Airfoil_Surface_Mesh.h5')
    if os.path.exists(airfoil_mesh) == True and reload == False:
        # Loading the mesh
        r = Reader('hdf_antares')
        r['filename'] = airfoil_mesh
        mesh = r.read()
        mesh.show()
        nodes = mesh[0][0]['x'].shape[0]
    else:
        ## Surface Extraction of airfoil 
        print('----> Beginning Surface Extraction')
        ## Loading the mesh
        print('      Loading the Main LES Mesh File')
        ## Loading the Main LES Mesh File
        r = Reader('hdf_avbp')
        r['filename'] = mesh_file
        base  = r.read() # b is the Base object of the Antares API
        print('\n----> The complete mesh:')
        base.show()
        airfoil_base = Family()
        for surf_name in input_surface:
            try:
                airfoil_base[surf_name] = base.families['Patches'][surf_name]
                print('      Surface %s extracted' %surf_name)
            except: 
                print('      Warning: Surface %s not found in the mesh file.' %surf_name)
        base.families['SKIN'] = airfoil_base
        skin_base = base[base.families['SKIN']]
        # The data for the original mesh
        print('\n----> The Original Mesh Surface:')
        skin_base.show()
        ## Merging the extracted base objects to the same zone
        print('\n----> Merging the Base Objects')
        myt = Treatment('merge')
        myt['base'] = skin_base
        myt['duplicates_detection'] = False
        myt['tolerance_decimals'] = 13
        # Writing the extraced mesh
        print('----> Writing the Mesh File')
        merged = myt.execute()
        writer = Writer('hdf_antares')
        writer['base'] = merged
        writer['filename'] = airfoil_mesh.replace('.h5','')
        writer.dump()
        # The data for the extracted and merged mesh
        print('----> The Post Extraced Mesh Surface:')
        merged.show()
        if len(input_surface)==1:
            nodes = merged[surf_name][0].shape[0]
        else: 
            nodes = merged['0000'][0].shape[0]
    print('\n----> Output Statistics')
    print('      The Extracted surface mesh is saved in: %s' %airfoil_mesh)
    print('      The number of nodes in on the airofoil surface is: %d nodes' %nodes)
    text = 'Surface Extraction Complete!'
    print(f'\n{text:.^80}\n')  
    airfoil_mesh = os.path.abspath(airfoil_mesh)
    return airfoil_mesh


def extract_data(working_dir, data_dir, airfoil_mesh, dtype='float64', reload:bool=False):
    """
    This function copies the surface pressure data from the FWH files in the specified FWH source directory and saves it as an equivalent HDF5 file in the working directory.
    
    Args:
        working_dir (str): Path to the working directory where the extracted data will be saved.
        data_dir (str): Path to the FWH source directory containing the pressure data files.
        airfoil_mesh (str): Path to the airfoil surface mesh file.
        dtype (str, optional): Data type for the pressure data. Defaults to 'float64'.
        reload (bool, optional): If True, forces re-extraction of data even if it already exists. Defaults to False.
    Returns: 
        surface_pressure_data (str): Path to the extracted surface pressure data HDF5 file.
        dt (float): Time step between successive samples.
    """
    
    text = 'Beginning Data Extraction'
    print(f'\n{text:.^80}\n') 
    if os.path.exists(os.path.join(working_dir,'pressure_airfoil.hdf5')) == True and reload == False:
        print('----> The pressure data is already extracted at: {0:s}'.format(working_dir,'pressure_airfoil.hdf5'))
        print('\n----> Loading the pressure data')
        surface_pressure_data = os.path.join(working_dir,'pressure_airfoil.hdf5')
        with h5py.File(surface_pressure_data, 'r') as f:
            dt = f.attrs['dt']
            last_extracted = f.attrs['Extracted Date']
            print('      The pressure data was last extracted on: %s' %last_extracted)
    else:  
        print('\n----> Extacting the pressure data from FWH files')
        # The directory information
        l=sorted(glob.glob(os.path.join(data_dir,'FWH_Airfoil_0000*.h5')))
        # The number of files (timesteps)
        nb_files=len(l)
        print('      The number of time steps is %d\n' %nb_files)
        # Extract the number of nodal points from the mesh
        r = Reader('hdf_antares')
        r['filename'] = airfoil_mesh
        base  = r.read() # b is the Base object of the Antares API
        nb_points = int(base[0][0].shape[0])
        # Pre allocating space for the pressure data array of size nb_points (nodes) x nb_files (time steps)
        data = np.zeros((nb_points,nb_files), dtype=dtype)
        data_time = np.zeros((nb_files), dtype=dtype)
        print('----> The surface data will be extracted to a %d (nodes) x %d (timestep), array' %(nb_points,nb_files))
        # Running the loop to extract the pressure data and export as hdf5 file
        for it,filename in enumerate(l):
                print('      Extracting file %s ...' %os.path.split(filename)[1]) if np.mod(it,100) == 0 else None
                with h5py.File(filename, 'r') as f:
                        # The pressure array from FWH data into 1D 
                        press = f['frame_data/pressure'][()]
                        if np.any(press == 0):
                            print("Warning: Zero values detected in pressure data.")
                        # Check if the number of pressure points match the number of nodal points
                        if len(press) != nb_points:
                                raise ValueError(f"Number of pressure points in {filename} ({len(press)}) does not match expected number of nodal points ({nb_points})")
                        # Saving the pressure data into the array
                        time = f['frame_data/time'][()]
                        data[:,it]=press.astype(dtype)
                        data_time[it] = time[0]
                print('      The last file extracted is %s' %os.path.split(filename)[1]) if it == nb_files-1 else None
        # Saving the output as hdf5 file
        dt = np.mean(np.diff(data_time))
        data = replace_zeros_vectorized(data)
        mean_pressure = np.mean(data, axis=1)  # Mean across timesteps
        rms_pressure = np.sqrt(np.mean((data - mean_pressure[:, None]) ** 2, axis=1))  # RMS using mean-subtracted values
        
        surface_pressure_data = os.path.join(working_dir,'pressure_airfoil.hdf5')
        surface_pressure_data = os.path.abspath(surface_pressure_data)
        with h5py.File(surface_pressure_data, 'w') as f:
            f.create_dataset('pressure', data=data, dtype=dtype)
            f.create_dataset("mean_pressure", data=mean_pressure, dtype=dtype)
            f.create_dataset("rms_pressure", data=rms_pressure, dtype=dtype)
            f.create_dataset("min_pressure", data=np.min(data, axis=1), dtype=dtype)
            f.create_dataset("max_pressure", data=np.max(data, axis=1), dtype=dtype)
            f.attrs['dt'] = dt
            f.attrs['Extracted Date'] = datetime.now().strftime("%Y-%m-%d")
            f.attrs['Source Path'] = data_dir
            f.attrs['Mesh Path'] = airfoil_mesh
        # Path to the surface pressure data
    print('\n----> Statistics of the extracted pressure data')
    print('      The pressure data is saved in: %s' %surface_pressure_data)
    print('      The pressure data storage size is {0:2.4f} MB' .format(os.path.getsize(surface_pressure_data)/1e6))
    print('      The time step is dt = {0:5.6e}' .format(dt))
    text = 'Data Extraction Complete!'
    print(f'\n{text:.^80}\n')  
    return surface_pressure_data, dt


def extract_surface_line( airfoil_mesh: str, 
    surface_pressure_data: str,
    AoA: int = 10,
    Uinf: int = 30,
    orientation: list = [1, 0, 0],
    cut_loc_percent: float = 95,
    z_loc_percent: float = 10,
    cut_length_percent: float = 10,
    data_size: int = 1000,
    airfoil_file : str = None,
    camber_file : str = None
    ):
    """
    Extract a 1D line of surface nodes (nearest-neighbor sampled from the full surface mesh) and save the corresponding pressure time series to an HDF5 file.
    The line can be either spanwise (x-normal, xi) or chordwise (z-normal, zeta) based on the orientation vector. 
    Note that the data is sampled from the original surface mesh, not the cut mesh, to ensure proper node correspondence.
    The spatial resolution is also controlled by the `data_size` parameter, which defines the number of points along the extracted line which are upsampled from the original mesh.
    
    The output surface line can be used to calculate spanwise or chordwise two-point coherence statistics, space-time coherence, or spatial-temporal correlations of the surface pressure data.
    
    Args:
        airfoil_mesh (str): Path to the airfoil surface mesh file.
        surface_pressure_data (str): Path to the HDF5 file containing surface pressure data.
        AoA (int): Angle of attack in degrees. Default is 10.
        Uinf (int): Free-stream velocity. Default is 30.
        orientation (list): Orientation vector for the cut plane. Default is [1, 0, 0] for spanwise.
        cut_loc_percent (float): Location of the cut along the chord as a percentage. Default is 95.
        z_loc_percent (float): Location of the cut along the span as a percentage. Default is 10.
        cut_length_percent (float): Length of the cut along the chord as a percentage. Default is 10.
        data_size (int): Number of points along the extracted line. Default is 1000.
        airfoil_file (str): Path to the airfoil coordinate file for chordwise extraction. Default is None.
        camber_file (str): Path to the airfoil camber line file for chordwise extraction. Default is None.
    
    Returns:
        str: Path to the output HDF5 file containing the extracted surface line data.
    """

    assert orientation in [[1, 0, 0], [0, 0, 1]], \
        "The orientation vector must be [1,0,0] for spanwise, or [0,0,1] for chordwise"

    text = "Beginning Surface Line Extraction"
    print(f"\n{text:.^80}\n")
    # Load airfoil surface mesh
    r = Reader("hdf_antares")
    r["filename"] = os.path.join(airfoil_mesh)
    print("----> Loading the airfoil surface mesh")
    mesh = r.read()
    mesh.compute_cell_volume()
    mesh.cell_to_node()
    mesh = mesh[:, :, [("x", "node"), ("y", "node"), ("z", "node"), ("cell_volume", "node")]]
    mesh.show()
    
    # Load surface pressure data
    var = "pressure"
    print("\n----> Open results file at: \n{0:s}".format(surface_pressure_data))
    with h5py.File(os.path.join(surface_pressure_data), "r") as fin:
        data = {var: fin[f"/{var}"][:]}
        dt = fin.attrs["dt"]
    print("      The data shape is: {0}".format(data[var].shape))

    # Create a new base for cut treatment
    b = Base()
    b["0"] = Zone()
    b[0].shared["cell_volume"] = mesh[0][0]["cell_volume"]
    b[0].shared["x"] = mesh[0][0]["x"]
    b[0].shared["y"] = mesh[0][0]["y"]
    b[0].shared["z"] = mesh[0][0]["z"]

    # Add a few instants (Antares cut expects instants present!)
    n_inst = min(10, int(data[var].shape[1]))
    for i in range(n_inst):
        instant_name = f"snapshot_{i}"
        b[0][instant_name] = Instant()
        b[0][instant_name][var] = np.asarray(data[var][:, i], dtype=np.float64).ravel(order="F")

    # Define cut plane
    # chord length computed from x-extent corrected by AoA (your original logic)
    x_nodes = mesh[0][0]["x"]
    chord_length = (x_nodes.max() - x_nodes.min()) / np.cos(AoA * np.pi / 180.0)
    cut_loc = x_nodes.min() + cut_loc_percent * chord_length / 100.0
    print("\n----> Applying the surface cut")
    t = Treatment("cut")
    t["base"] = b
    t["expert"] = "tri"
    t["type"] = "plane"
    
    # origin definition (kept as you had it)
    x_o, y_o, z_o = cut_loc, 0.0, -0.1034 - chord_length * z_loc_percent / 100.0
    t["origin"] = [float(x_o), float(y_o), float(z_o)]
    t["with_boundaries"] = True
    x_n, y_n, z_n = orientation
    t["normal"] = [int(x_n), int(y_n), int(z_n)]
    b_cut = t.execute()
    b_cut.show()
    print("      The cut origin is: (x,y,z) = (%f,%f,%f) with a normal vector of (x,y,z) = (%d,%d,%d).\n"
        % (x_o, y_o, z_o, x_n, y_n, z_n))
    
    # Build clean source/dest arrays (N,3) float64 with NaN/Inf guards
    # The Source is the original mesh nodes from the mesh file
    x_src = np.asarray(b[0].shared["x"], dtype=np.float64).ravel()
    y_src = np.asarray(b[0].shared["y"], dtype=np.float64).ravel()
    z_src = np.asarray(b[0].shared["z"], dtype=np.float64).ravel()
    source = np.column_stack((x_src, y_src, z_src))  # (N,3)
    # Guard against non-finite values in the source mesh
    if source.ndim != 2 or source.shape[1] != 3:
        raise ValueError(f"KDTree source must be (N,3). Got {source.shape}.")
    finite_mask = np.isfinite(source).all(axis=1)
    if not np.all(finite_mask):
        bad = np.where(~finite_mask)[0]
        print(f"WARNING: {bad.size} non-finite coordinate rows in source (showing first 10): {bad[:10]}")
        source_valid = source[finite_mask]
        orig_idx = np.nonzero(finite_mask)[0]  # mapping back to original indices
    else:
        source_valid = source
        orig_idx = None
        
    # Dest = cut mesh nodes (only used for spanwise sampling bounds in your logic)
    x_dst = np.asarray(b_cut[0].shared["x"], dtype=np.float64).ravel()
    y_dst = np.asarray(b_cut[0].shared["y"], dtype=np.float64).ravel()
    z_dst = np.asarray(b_cut[0].shared["z"], dtype=np.float64).ravel()
    dest = np.column_stack((x_dst, y_dst, z_dst))
    
    # Build sample points along requested line
    sample = np.zeros((data_size, 3), dtype=np.float64)
    
    # The Logic for Spanwise aligned line, this requies a cut normal to the x-axis (x_n=1, streamwise-normal) to create a 2-D plane, then extracting the surface points only from the 2D plane
    if x_n == 1:
        # Spanwise line at fixed x=cut_loc, fixed y=max(dest_y), varying z
        print("      Applying cut along the suction span at {0:.4f} % chord length from Leading Edge".format(
                cut_loc_percent
            )
        )
        # guard: need at least a handful of dest points
        if dest.shape[0] < 10:
            raise ValueError(f"Cut mesh has too few points ({dest.shape[0]}) to define a spanwise line safely.")
        
        z_sorted = np.sort(dest[:, 2])
        z_lo = z_sorted[3]
        z_hi = z_sorted[-3]
        length_arr = np.linspace(z_lo, z_hi, data_size)
        
        y_const = float(np.max(dest[:, 1]))
        for i in range(data_size):
            sample[i, :] = [float(cut_loc), y_const, float(length_arr[i])]
            
    # The Logic for Chordwise aligned line, this requies a cut normal to the z-axis (z_n=1, spanwise-normal) to create a 2-D plane, then extracting the surface points only from the 2D plane
    elif z_n == 1:
        # NOTE: requies a physical airfoil coordinate file and camber line file to identify the suction surface! 
        # Chordwise line at fixed z=z_o, using airfoil coordinate file for (x,y) curve
        print("      Applying cut along the suction chord at {0:.2f} % chord length from the tip".format(z_loc_percent))
        
        # You used AoA-specific coord files; preserve that convention
        airfoil_file = f"/project/rrg-plavoie/denggua1/Coordinates/airfoil_{AoA}_AOA.txt" if airfoil_file is not None else airfoil_file
        camber_file = f"/project/rrg-plavoie/denggua1/Coordinates/airfoil_camber_{AoA}_AOA.txt" if camber_file is not None else camber_file
        if not os.path.exists(airfoil_file):
            raise FileNotFoundError(f"Airfoil coordinate file not found: {airfoil_file}")
        if not os.path.exists(camber_file):
            raise FileNotFoundError(f"Camber coordinate file not found: {camber_file}")
        
        airfoil = np.loadtxt(airfoil_file, dtype=float, delimiter=",")
        _ = np.loadtxt(camber_file, dtype=float, delimiter=",")  # kept for parity; not used below
        
        dummy_x = airfoil[0:100, 0]
        dummy_y = airfoil[0:100, 1]
        airfoil_x = dummy_x[::-1]
        airfoil_y = dummy_y[::-1]
        
        # Create an interpolation function for the airfoil coordinates to obtain the y-coordinate for any x-coordinate along the chord
        f = scipy.interpolate.interp1d(airfoil_x, airfoil_y, fill_value="extrapolate")
        
        x0p = max(0.0, cut_loc_percent - cut_length_percent)
        x1p = min(100.0, cut_loc_percent + cut_length_percent)
        # Map percent -> index in [0, len(airfoil_x)-1]
        i0 = int(round((x0p / 100.0) * (len(airfoil_x) - 1)))
        i1 = int(round((x1p / 100.0) * (len(airfoil_x) - 1)))
        i0, i1 = min(i0, i1), max(i0, i1)
        
        airfoil_x_interp = np.linspace(airfoil_x[i0], airfoil_x[i1], data_size)
        airfoil_y_interp = f(airfoil_x_interp)
        
        for i in range(data_size):
            sample[i, :] = [float(airfoil_x_interp[i]), float(airfoil_y_interp[i]), float(z_o)]
            
    # Sanity check for sample points
    if not np.isfinite(sample).all():
        raise ValueError("Sample points contain NaN/Inf; check cut definition / interpolation inputs.")
    
    
    # Nearest-neighbor query using cKDTree (vectorized) to find the closest nodes in the original mesh to the sample points along the cut line
    print("\n---> Building cKDTree and querying nearest nodes")
    tree = cKDTree(source_valid, leafsize=64, compact_nodes=True, balanced_tree=True)
    dist, idx = tree.query(sample, k=1)
    idx = idx.astype(np.int64)
    # Map back to original indices if we removed non-finite rows
    index = orig_idx[idx] if orig_idx is not None else idx
    
    # Extract time series at selected nodes (vectorized)
    print("\n---> Extracting Time Series to the Surface Line")
    x_line = source[index, 0]
    y_line = source[index, 1]
    z_line = source[index, 2]
    # data[var] is (n_nodes, n_timesteps) in your extraction
    # Pull all nodes for all timesteps at once: (data_size, n_timesteps)
    time_series = np.asarray(data[var][index, :], dtype=np.float64)
    print("      The timeseries size %d nodes x %d timesteps\n" % (time_series.shape[0], time_series.shape[1]))
    
    
    # Save output
    print("---> Saving the Surface Line Data to HDF5 file")
    orient = "span" if x_n == 1 else "chord"
    outputname = f"B_{AoA}AOA_U{Uinf}_LES_{cut_loc_percent}c_{z_loc_percent}z_{orient}_Surface_Line.hdf5"

    with h5py.File(outputname, "w") as file:
        grid_data = file.create_group("Grid Data")
        grid_data.create_dataset("x", data=x_line.astype(np.float64).reshape(1, -1))
        grid_data.create_dataset("y", data=y_line.astype(np.float64).reshape(1, -1))
        grid_data.create_dataset("z", data=z_line.astype(np.float64).reshape(1, -1))
        file.attrs["cut_loc_percent"] = float(cut_loc_percent)
        file.attrs["z_loc_percent"] = float(z_loc_percent)
        file.attrs["cut_length_percent"] = float(cut_length_percent)
        file.attrs["origin"] = [float(x_o), float(y_o), float(z_o)]
        file.attrs["normal"] = [int(x_n), int(y_n), int(z_n)]
        file.attrs["dt"] = float(dt)
        time_data = file.create_group("Timeseries Data")
        time_data.create_dataset("time_series", data=time_series.astype(np.float64))

    print("      The surface line data is saved in:\n %s" % outputname)
    text = "Surface Line Extraction Complete"
    print(f"\n{text:.^80}\n")
    return outputname
