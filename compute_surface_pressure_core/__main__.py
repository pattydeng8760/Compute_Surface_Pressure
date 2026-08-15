"""
Module runner for the compute_surface_pressure_core package.
 
Allows the package to be invoked directly as a script:
 
    python -m compute_surface_pressure_core --sol-dir ... --mesh-file ... --psd
 
This is the internal (core-module) entry point. For use outside of the core
package directory, see the standalone `run_compute_surface_pressure.py`
script and `run_compute_surface_pressure.sh` SLURM wrapper distributed
alongside this package.
"""
 
import sys
 
from .ComputeSurfacePressure import main

if __name__ == "__main__":
    main(sys.argv[1:])