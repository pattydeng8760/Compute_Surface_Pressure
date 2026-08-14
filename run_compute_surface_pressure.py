#!/usr/bin/env python
"""
run_compute_surface_pressure.py
================================

External, top-level entry point for the `compute_surface_pressure_core`
package. This script is meant to live OUTSIDE of the core package directory
(e.g. one level up, or anywhere on the cluster) so it can be called directly
from the command line or from a SLURM batch script
(`run_compute_surface_pressure.sh`) without needing to know Python packaging
details.

It simply forwards all CLI arguments to `compute_surface_pressure_core.main`.

Usage
-----
    python run_compute_surface_pressure.py \\
        --sol-dir /path/to/AVBP/FWH/solution \\
        --mesh-file /path/to/mesh.mesh.h5 \\
        --surface-patches Airfoil_Surface \\
        --working-dir /path/to/working/dir \\
        --psd --csd \\
        --mic-dir /path/to/mic/dir --mic-file mic_data.hdf5 --mic-num 3 \\
        --freq-select 500 1000 2000

"""

import os
import sys

# Convenience fallback: if compute_surface_pressure_core isn't installed,
# allow running this script directly from the repository root without
# requiring an editable install.
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)

from compute_surface_pressure_core import main


if __name__ == "__main__":
    main(sys.argv[1:])