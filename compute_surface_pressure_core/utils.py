"""
This module contains utility functions for the compute_surface_pressure_core package.
"""

def timer(func):
    """ Decorator to time the function func to track the time taken for the function to run"""
    def inner(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        elapsed = end - start
        print('The total compute time is: {0:1.0f} s'.format(elapsed))
        return elapsed
    return inner


def next_greater_power_of_2(x):
    return 2**(x-1).bit_length()


def fft_surface_data(surface_pressure_data:str, var:str, dt:float, weight='default', nOvlp=128, nDFT=256, window='default', method='fast', reload=False):
    """
    This function performs a Fast Fourier Transform (FFT) on surface pressure data stored in an HDF5 file.
    The compute is done in blocks to manage memory usage, and the results are saved in a new HDF5 file. 
    NOTE: It is recommended to use welch method for PSD and CSD calculations instead of this function, as it is more efficient and robust. But it will output almost the same result as the welch method.
    
    Args:
        surface_pressure_data (str): Path to the HDF5 file containing surface pressure data.
        var (str): Variable name in the HDF5 file to extract.
        dt (float): Time-step between successive samples.
        weight (str or array-like, optional): Weighting function for the FFT. Defaults to 'default'.
        nOvlp (int, optional): Number of overlapping points between blocks. Defaults to 128.
        nDFT (int, optional): Number of points for the DFT. Defaults to 256.
        window (str or array-like, optional): Window function for the FFT. Defaults to 'default'.
        method (str, optional): Method for FFT computation ('fast' or 'lowRAM'). Defaults to 'fast'.
        reload (bool, optional): If True, forces recalculation even if the output file exists. Defaults to False.
    
    Returns:
        str: Path to the output HDF5 file containing the FFT results.
    """
    
    text = 'Performing FFT on surface pressure data'
    print(f'\n{text:.^80}\n')
    fft_file_path = surface_pressure_data.replace('.hdf5', '_fft.hdf5')
    if os.path.exists(fft_file_path) == True and reload == False:  
        pass 
    else:
        # Load the pressure data
        with h5py.File(surface_pressure_data , 'r') as h5f:
        # Extracting the data from the file as a numpy array
            data = h5f[var][:]  # Read the flow field data
            # Check for the chape o
            if np.shape(data)[0] > np.shape(data)[1]:
                data = np.swapaxes(data, 0, 1) # Swap the axes to have the time as the first dimension
        nt = np.shape(data)[0]
        nx = np.shape(data)[1] 
        print('The surface data is %d (nodes) x %d (timestep)' %(nx,nt))
        # SPOD parser
        [weight, window, nOvlp, nDFT, nBlks] = spod_parser(nt, nx, window, weight, nOvlp, nDFT, method)
        print('Calculating temporal DFT'              )
        print('--------------------------------------')
        # calculate time-averaged result
        x_mean  = np.mean(data,axis=0)
        # obtain frequency axis
        f     = np.arange(0,int(np.ceil(nDFT/2)+1))
        f     = f/dt/nDFT
        nFreq = f.shape[0]   
        # initialize all DFT result in frequency domain
        if method == 'fast':
            Q_hat = np.zeros((nx,nFreq,nBlks),dtype = complex) # RAM demanding here       
            
        elif method == 'lowRAM':
            Q_hat = h5py.File(os.path.join(save_path,'Q_hat.h5'), 'w')
            Q_hat.create_dataset('Q_hat', shape=(nx,nFreq,nBlks), chunks=True, dtype = complex, compression="gzip")
        # initialize block data in time domain
        Q_blk = np.zeros((nDFT,nx))
        Q_blk_hat = np.zeros((nx,nFreq),dtype = complex)
        # loop over each block
        for iBlk in range(nBlks):
            # get time index for present block
            it_end   = min(iBlk*(nDFT-nOvlp)+nDFT, nt)
            it_start = it_end-nDFT
            print('block {0:d} / {1:d} ({2:d} : {3:d})'.format(iBlk+1, nBlks, it_start+1, it_end))
            # subtract time-averaged results from the block
            Q_blk = data[it_start:it_end,:] - x_mean # column-wise broadcasting
            # add window function to block
            Q_blk = Q_blk.T * window # row-wise broadcasting
            # Fourier transform on block
            Q_blk_hat = 1/np.mean(window)/nDFT*fft(Q_blk)       
            Q_blk_hat = Q_blk_hat[:,0:nFreq]           
            # correct Fourier coefficients for one-sided spectrum
            Q_blk_hat[:,1:(nFreq-1)] *= 2
            # save block result to the whole domain result 
            if method == 'fast':
                Q_hat[:,:,iBlk] = Q_blk_hat
                
            elif method == 'lowRAM':
                Q_hat['Q_hat'][:,:,iBlk] = Q_blk_hat
            # remove vars to release RAM
        del data, Q_blk, Q_blk_hat
        print('--------------------------------------')
        print('Calculating FFT'                       )
        print('--------------------------------------')
        # # initialize output vars
        # if method == 'fast':
        #     L = np.zeros((nFreq,nBlks))
        #     P = np.zeros((nFreq,nx,nBlks),dtype = complex) # RAM demanding here
            
        # elif method == 'lowRAM':
        #     h5f = h5py.File(sol_file, 'w')
        #     h5f.create_dataset('L', shape=(nFreq,nBlks), compression="gzip")
        #     h5f.create_dataset('P', shape=(nFreq,nx,nBlks), chunks=True, dtype = complex, compression="gzip")
        #     h5f.create_dataset('f', data=f, compression="gzip")     
        # Initialize data for the FFT of the pressure data
        data_fft = np.zeros((nx, nFreq), dtype=complex)
        # loop over each frequency
        for iFreq in range(nFreq):
            print('Frequency {0:d} / {1:d} (f = {2:3.3f})'.format(iFreq+1,nFreq,f[iFreq]))
            if method == 'fast':
                Q_hat_f = Q_hat[:,iFreq,:]
                data_fft[:,iFreq] = np.mean(np.abs(Q_hat_f),axis=1)
            elif method == 'lowRAM':
                Q_hat_f = Q_hat['Q_hat'][:,iFreq,:]
                data_fft[:,iFreq] = np.mean(np.abs(Q_hat_f),axis=1)
        print('After FFT, the data is %d (nodes) x %d (frequency bins)' %(nx,int(len(f))))        
        # Saving the Fourier Transform output and frequency information as an hdf5 file
        with h5py.File(fft_file_path, 'w') as hdf:
            hdf.create_dataset('pressure_fft', data=data_fft, dtype=complex)
            hdf.create_dataset('frequency', data=f, dtype=float)
    # Printing information regaring the post-fft data
    print('\nThe Fourier Transform of the pressure data is saved in {0:s}' .format(fft_file_path ))
    print('The Fourier Transform of the pressure data storage size is {0:2.4f} MB' .format(os.path.getsize(fft_file_path )/1e6))
    text = 'FFT Complete!'
    print(f'\n{text:.^80}\n')
    return fft_file_path


def spod_parser(nt, nx, window, weight, nOvlp, nDFT, method):
    '''
    Purpose: determine data structure/shape for SPOD
    
    Parameters
    ----------
    nt     : int; number of time snapshots
    nx     : int; number of grid point * number of variable
    window : expect 1D numpy array, float; specified window function values
    weight : expect 1D numpy array; specified weight function
    nOvlp  : expect int; specified number of overlap
    nDFT   : expect int; specified number of DFT points (expect to be same as weight.shape[0])
    method : expect string; specified running mode of SPOD
    
    Returns
    -------
    weight : 1D numpy array; calculated/specified weight function
    window : 1D numpy array, float; window function values
    nOvlp  : int; calculated/specified number of overlap
    nDFT   : int; calculated/specified number of DFT points
    nBlks  : int; calculated/specified number of blocks
    '''
    # check SPOD running method
    try:
        # user specified method
        if method not in ['fast', 'lowRAM']:
            print('WARNING: user specified method not supported')
            raise ValueError
        else:
            print('Using user specified method...')
    except:        
        # default method
        method = 'lowRAM'
        print('Using default low RAM method...')
    # check specified weight function value
    try:
        # user specified weight
        nweight = weight.shape[0]
        if nweight != nx:
            print('WARNING: weight does not match with data')
            raise ValueError
        else:
            wgt_name = 'user specified'
            print('Using user specified weight...')
    except:        
        # default weight
        weight   = np.ones(nx)
        wgt_name = 'unity'
        print('Using default weight...')
    # calculate or specify window function value
    try:
        # user sepcified window
        nWinLen  = window.shape[0]
        win_name = 'user specified'
        nDFT     = nWinLen # use window shape to over-write nDFT (if specified)
        print('Using user specified nDFT from window length...')         
        print('Using user specified window function...')  
        
    except:
        # default window with specified/default nDFT
        try:
            # user specified nDFT
            nDFT  = int(nDFT)
            nDFT  = int(2**(np.floor(np.log2(nDFT)))) # round-up to 2**n type int
            print('Using user specified nDFT ...')             
                
        except:
            # default nDFT
            nDFT  = 2**(np.floor(np.log2(nt/10)))
            nDFT  = int(nDFT)
            print('Using default nDFT...')
            
        window   = hammwin(nDFT)
        win_name = 'Hamming'
        print('Using default Hamming window...') 
    # calculate or specify nOvlp
    try:
        # user specified nOvlp
        nOvlp = int(nOvlp)
        
        # test feasibility
        if nOvlp > nDFT-1:
            print('WARNING: nOvlp too large')
            raise ValueError
        else:
            print('Using user specified nOvlp...')
    except:            
        # default nOvlp
        nOvlp = int(np.floor(nDFT/2))
        print('Using default nOvlp...')
    # calculate nBlks from nOvlp and nDFT    
    nBlks = int(np.floor((nt-nOvlp)/(nDFT-nOvlp)))
    # test feasibility
    if (nDFT < 4) or (nBlks < 2):
        raise ValueError('User sepcified window and nOvlp leads to wrong nDFT and nBlk.')
    print('--------------------------------------')
    print('SPOD parameters summary:'              )
    print('--------------------------------------')
    print('number of DFT points :{0:d}'.format(int(nDFT)))
    print('number of blocks is  :{0:d}'.format(int(nBlks)))
    print('number of overlap percent is :{0:d}'.format(nOvlp))
    print('Window function      :{0:s}'.format(win_name))
    print('Weight function      :{0:s}'.format(wgt_name))
    print('Running method       :{0:s}'.format(method))
    
    return weight, window, nOvlp, nDFT, nBlks

def hammwin(N):
    '''
    Purpose: standard Hamming window
    
    Parameters
    ----------
    N : int; window lengh

    Returns
    -------
    window : 1D numpy array; containing window function values
             n = nDFT
    '''
    
    window = np.arange(0, N)
    window = 0.54 - 0.46*np.cos(2*np.pi*window/(N-1))
    window = np.array(window)

    return window


def replace_zeros_vectorized(data):
    """ Efficiently replaces zeros in a 2D array with the nearest nonzero neighbor using forward & backward filling.
        Prints the node index and replaced value for each zero detected.
    """
    
    mask = data == 0  # Boolean mask of zero entries
    
    if not np.any(mask):  # If no zeros are found, return early
        print("No zeros detected. Skipping replacement.")
        return data  

    print("Replacing zeros with nearest neighbor values...")

    # Forward fill: Copy last nonzero value forward
    for i in range(1, data.shape[1]):
        update_mask = mask[:, i] & (data[:, i - 1] != 0)  # Only replace where needed
        if np.any(update_mask):  # Print details of replacement
            replaced_nodes = np.where(update_mask)[0]  # Indices of affected nodes
            replaced_values = data[replaced_nodes, i - 1]
            for node, value in zip(replaced_nodes, replaced_values):
                print(f"Node {node}: Zero replaced with forward value {value}")
        data[:, i] = np.where(update_mask, data[:, i - 1], data[:, i])

    # Backward fill: Copy first nonzero value backward
    for i in range(data.shape[1] - 2, -1, -1):
        update_mask = mask[:, i] & (data[:, i + 1] != 0)
        if np.any(update_mask):  # Print details of replacement
            replaced_nodes = np.where(update_mask)[0]
            replaced_values = data[replaced_nodes, i + 1]
            for node, value in zip(replaced_nodes, replaced_values):
                print(f"Node {node}: Zero replaced with backward value {value}")
        data[:, i] = np.where(update_mask, data[:, i + 1], data[:, i])

    return data
