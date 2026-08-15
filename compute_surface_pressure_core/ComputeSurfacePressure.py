import os
import sys
import time
import logging
import glob
import argparse
from types import SimpleNamespace
from datetime import datetime
import numpy as np
import scipy
import scipy.signal as signal
import h5py

from .extractor import extract_files, extract_surface, extract_data, extract_surface_line
from .surface_psd import PSD_surface_data, source_psd
from .surface_csd import CSD_surface_data, source_csd
from .utils import fft_surface_data, timer

def parse_arguments(argv=None):
    """
    This is the CLI parser for the compute_surface_pressure_core module. It defines the command-line arguments that can be passed to the script and returns the parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Compute surface pressure post-processing from FWH data."
    )
    # Manditory arguments
    parser.add_argument("--sol-dir", type=str, required=True, help="Solution directory for the FWH data from AVBP.")
    parser.add_argument("--freq-select",type=float, nargs="+",default=[1000, 1500, 2000, 2500, 3000], help="Comma-separated list of frequencies to select for analysis. Default is 1000 to 3000 Hz in steps of 500 Hz as list.")
    parser.add_argument("--mesh-file", type=str, required=True, help="AVBP mesh file.")
    parser.add_argument("--surface-patches", "-sp",type=str, nargs="+",default=["Airfoil_Surface"], help="List of surface patches to include in the analysis.")
    parser.add_argument("--working-dir", type=str, default=os.getcwd(), help="Destination directory for the extracted files. Default is current working directory.")
    parser.add_argument("--option", type=int, default=1, help="Option for extraction. Default is 1.")
    parser.add_argument("--nskip", type=int, default=1, help="Number of files to skip during extraction. Default is 1.")
    parser.add_argument("--max-file", type=int, default=5000, help="Maximum number of files to process. Default is 5000.")
    parser.add_argument("--reload", action="store_true", help="Flag to re-extract all data. If set, will re-extract data even if it already exists.") 
    parser.add_argument("--var", type=str, default="pressure", help="Variable to process. Default is 'pressure'.")
    
    # Compute specific workflow
    parser.add_argument("--csd", action="store_true", help="Flag to compute CSD. If set, will compute the cross-spectral density.")
    parser.add_argument("--psd", action="store_true", help="Flag to compute PSD. If set, will compute the power spectral density.")
    parser.add_argument("--fft", action="store_true", help="Flag to compute FFT. If set, will compute the fast Fourier transform.")
    parser.add_argument("--surf_line", action="store_true", help="Flag to extract surface line data. If set, will extract surface line data.")
    parser.add_argument("--block-size", type=int, default=1000, help="Block size for PSD/CSD computation. Default is 1000.")
    parser.add_argument("--nchunk", type=int, default=4, help="Number of chunks for PSD/CSD computation. Default is 4.")
    
    # Arguments required for PSD 
    parser.add_argument("--band-stat", type=str, default="mean", choices=["mean", "sum"], help="Statistic to use for averaging PSD over frequency bands. Default is 'mean'.")
    parser.add_argument("--df", type=float, default=50.0, help="Frequency resolution for PSD computation if averaging frequencies. Default is 50 Hz.")
    
    # Arguments required for CSD
    parser.add_argument("--mic-dir", type=str,  help="Directory for the FWH computed microphone data. Should be a hdf5 file containing microphone data.")
    parser.add_argument("--mic-file", type=str, help="Name of the microphone data file (hdf5 format).")
    parser.add_argument("--mic-num", type=int, default=1, help="Microphone reference for CSD computation. Should be an integer corresponding to the microphone number.")
    
    # Surface line extraction parameters
    parser.add_argument("--AoA", type=float, default=10.0, help="Angle of attack for surface line extraction. Default is 10 degrees.")
    parser.add_argument("--Uinf", type=float, default=30.0, help="Freestream velocity for surface line extraction. Default is 30 m/s.")
    parser.add_argument("--orientation", type=float, nargs=3, default=[0, 0, 1], help="Orientation vector for surface line extraction. Default is [0, 0, 1], the other option is [1,0,0]")
    # NOTE for spanwise cut all three arguments below are required, for streamwise cut only cut_loc_percent is required!
    parser.add_argument("--cut-loc-percent", type=float, default=70.0, help="Streamwise chord location of the spanwise cut as a percent chord for surface line extraction. Default is 70%. i.e., a spanwise cut along 70 percent chord along the whole span")
    parser.add_argument("--z-loc-percent", type=float, default=70.0, help="Spanwise center location of the spanwise cut as a percent span for surface line extraction. Default is 70%. i.e., a streamwise cut along 70 percent span")
    parser.add_argument("--cut-length-percent", type=float, default=10, help="Length of the spanwise cut as a percent span for surface line extraction. Default is 10%. i.e., a streamwise cut centered at 70 percent span and 70 percent chord with a length of 2*10 percent chord ")
    parser.add_argument("--airfoil-file", type=str, default=None, help="Optional airfoil file in x,y,z stored in .txt for surface line extraction. If provided, will use this airfoil file for extraction.")
    parser.add_argument("--camber-file", type=str, default=None, help="Optional camber line file in x,y,z stored in .txt for surface line extraction. If provided, will use this camber line file for extraction.")
    parser.add_argument("--data-size", type=int, default=1000, help="Number of data points to extract for surface line extraction. Default is 1000.")
    
    return parser.parse_args(argv)


class SurfacePressure():
    def __init__(self, args):
        self.args = args
        logfile = os.path.join(args.working_dir,'log_Surface_Pressure_'+datetime.now().strftime("%Y%m%d_%H%M")+'.txt')
        sys.stdout = open(logfile, "w", buffering=1)
        self.start_time = time.time()
        text = "Beginning Compute Surface Pressure Post-Processing"
        print(f'\n{text:=^100}\n')
        self.sol_dir = args.sol_dir
        self.freq_select = args.freq_select
        self.mesh_file = args.mesh_file
        self.surface_patches = args.surface_patches
        self.working_dir = args.working_dir
        self.option = args.option
        self.nskip = args.nskip
        self.max_file = args.max_file
        self.reload = args.reload
        self.var = args.var
        
        # PSD specific arguments
        self.psd = SimpleNamespace()
        self.psd.run = args.psd
        if self.psd.run:
            self.psd.band_stat = args.band_stat
            self.psd.df = args.df
            self.psd.block_size = args.block_size
            self.psd.nchunk = args.nchunk
        
        # CSD specific arguments
        self.csd = SimpleNamespace()
        self.csd.run = args.csd
        if self.csd.run:
            self.csd.mic_dir = args.mic_dir
            self.csd.mic_file = args.mic_file
            self.csd.mic_num = args.mic_num
            self.csd.block_size = args.block_size
            self.csd.nchunk = args.nchunk
        
        # FFT specific arguments
        self.fft = SimpleNamespace()
        self.fft.run = args.fft
        
        # Surface line extraction specific arguments
        self.surf_line = SimpleNamespace()
        self.surf_line.run = args.surf_line
        if self.surf_line.run:
            self.surf_line.data_size = args.data_size
            self.surf_line.AoA = args.AoA
            self.surf_line.Uinf = args.Uinf
            self.surf_line.orientation = args.orientation
            self.surf_line.cut_loc_percent = args.cut_loc_percent
            self.surf_line.z_loc_percent = args.z_loc_percent
            self.surf_line.cut_length_percent = args.cut_length_percent
            self.surf_line.airfoil_file = args.airfoil_file
            self.surf_line.camber_file = args.camber_file
            
        self._print_args()
    
    def _print_args(self):
        """Print all CLI / input arguments stored in self.args."""
        text = " Compute Surface Pressure Input Arguments "
        print(f'\n{text:.^80}\n')
        for name, value in vars(self.args).items():
            print(f"{name:20s}: {value}")
        text = " End of Input Arguments "
        print(f'\n{text:.^80}\n')
    
    def prepare_inputs(self):
        # Extracting the FWH data files and surface mesh
        try:
            self.ntime, self.FWH_data_dir = extract_files(self.sol_dir, self.working_dir, self.option, self.nskip, self.max_file, reload=self.reload)
        except: 
            logging.error("Error in extracting files. Please check the FWH solution directory and file names.")
            raise
        
        # Extracting the surface mesh from the AVBP mesh file
        try:
            self.airfoil_mesh = extract_surface(self.mesh_file, self.surface_patches, self.working_dir, reload=self.reload)
        except:
            logging.error("Error in extracting surface mesh. Please check the mesh file and surface patch names.")
            raise
        
        # Extracting the surface pressure data
        try:
            self.surface_pressure_data, self.dt = extract_data(self.working_dir,self.FWH_data_dir,self.airfoil_mesh,dtype='float64',reload=self.reload)
        except:
            logging.error("Error in extracting surface pressure data. Please check the FWH data files and surface mesh.")
            raise
    
    def run_compute(self):
        if self.psd.run:
            # The source PSD
            surface_pressure_psd_data = PSD_surface_data(self.surface_pressure_data, 
                self.var, 
                self.dt, 
                reload=self.reload, 
                block_size=self.psd.block_size,
                nchunk = self.psd.nchunk
            )
            
            source_psd(self.working_dir, 
                self.airfoil_mesh, 
                self.surface_pressure_data, 
                surface_pressure_psd_data, 
                freq_select = self.freq_select, 
                df=self.psd.df, 
                band_stat=self.psd.band_stat,
                pref=2e-5
            )
        
        if self.csd.run:
            # The source coherence
            surface_pressure_csd_data = CSD_surface_data(self.csd.mic_dir, 
                self.csd.mic_file, 
                self.csd.mic_num, 
                self.surface_pressure_data, 
                self.var, 
                self.dt, 
                reload=self.reload, 
                block_size=self.csd.block_size,
                nchunk=self.csd.nchunk
            )
            
            source_csd(self.working_dir, 
                self.airfoil_mesh, 
                self.surface_pressure_data, 
                surface_pressure_csd_data, 
                freq_select = self.freq_select
            )
        
        if self.fft.run:
            surface_pressure_fft_data = fft_surface_data(self.surface_pressure_data, 
                self.var, 
                self.dt, 
                weight='default',
                nOvlp=128,
                nDFT=256,
                window='default',
                method='fast', 
                reload=self.reload
            )
            
            source_fft(self.working_dir, 
                self.mesh_file, 
                self.airfoil_mesh, 
                self.surface_pressure_data, 
                surface_pressure_fft_data, 
                freq_select = self.freq_select
            )
        
        if self.surf_line.run:
            time_series = extract_surface_line(self.airfoil_mesh, self.surface_pressure_data, 
                AoA=self.surf_line.AoA,
                Uinf=self.surf_line.Uinf, 
                orientation=self.surf_line.orientation, 
                cut_loc_percent=self.surf_line.cut_loc_percent, 
                z_loc_percent=self.surf_line.z_loc_percent, 
                cut_length_percent=self.surf_line.cut_length_percent, 
                data_size=self.surf_line.data_size, 
                airfoil_file=self.surf_line.airfoil_file, 
                camber_file=self.surf_line.camber_file
            )
        self.end_time = time.time()
        print(f"\nTotal execution time: {self.end_time - self.start_time:.2f} s\n")
        text = "Complete Surface Pressure Post-Processing"
        print(f'\n{text:=^100}\n')

def main(argv=None):
    """
    Main entry point for surface-pressure post-processing.
    """
    args = parse_arguments(argv)
    surface_pressure = SurfacePressure(args)
    surface_pressure.prepare_inputs()
    surface_pressure.run_compute()


if __name__ == "__main__":
    main()