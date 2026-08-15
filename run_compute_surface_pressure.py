#!/usr/bin/env python3
"""Top-level runner for :mod:`compute_surface_pressure_core`.

There are two supported ways to use this file:

1. Edit ``CONFIG`` below and run::

       python run_compute_surface_pressure.py

2. Pass command-line arguments. If any CLI arguments are supplied, they are
   forwarded directly to ``compute_surface_pressure_core.main`` and ``CONFIG``
   is ignored::

       python run_compute_surface_pressure.py --sol-dir ... --mesh-file ... --psd

The companion ``run_compute_surface_pressure.sh`` uses the second mode.
"""

from __future__ import annotations

import os
import sys
from typing import Any


# =============================================================================
# USER CONFIGURATION
# =============================================================================
# Used only when this script is executed WITHOUT command-line arguments.
# Paths can be absolute or relative to the directory from which the script is
# launched.
CONFIG: dict[str, Any] = {
    # Required inputs
    "sol_dir": "/scratch/denggua1/Bombardier_LES/B_10AOA_LES/RUN_Fine_Jul25/FWH_Airfoil/FWH_Data_TTG",
    "mesh_file": "/project/rrg-plavoie/denggua1/BBDB_10AOA/MESH_Fine_Jul25/Bombardier_10AOA_U30_Combine_Fine.mesh.h5",

    # General settings
    "surface_patches": [
        "Airfoil_Surface",
        "Airfoil_Trailing_Edge",
        "Airfoil_Side_LE",
        "Airfoil_Side_Mid",
        "Airfoil_Side_TE",
    ],
    "working_dir": "./surface_pressure_output",
    "freq_select": [500, 1000, 1500, 2000, 2500, 3000],
    "option": 1,
    "nskip": 1,
    "max_file": 5000,
    "reload": False,
    "var": "pressure",

    # Workflows
    "psd": False,
    "csd": False,
    "fft": True,
    "surf_line": False,

    # PSD/CSD settings
    "block_size": 1000,
    "nchunk": 4,
    "band_stat": "mean",
    "df": 50.0,

    # CSD settings; required only when csd=True
    "mic_dir": "/scratch/denggua1/Bombardier_LES/B_10AOA_LES/RUN_Fine_Jul25/FWH_Airfoil/FWHpostpro_TTG",
    "mic_file": 'B_10AOA_U30_LES_Solid_Airfoil_Mic_Data.h5',
    "mic_num": 68,

    # Surface-line settings; used only when surf_line=True
    "AoA": 10.0,
    "Uinf": 30.0,
    "orientation": [0, 0, 1],
    "cut_loc_percent": 70.0,
    "z_loc_percent": 70.0,
    "cut_length_percent": 10.0,
    "airfoil_file": None,
    "camber_file": None,
    "data_size": 1000,
}


# Boolean options are emitted only when True. None-valued settings are skipped.
_BOOL_OPTIONS = {"reload", "psd", "csd", "fft", "surf_line"}


def _flag(name: str) -> str:
    """Convert a CONFIG key to its command-line spelling."""
    return "--" + name.replace("_", "-")


def build_argv(config: dict[str, Any]) -> list[str]:
    """Convert ``CONFIG`` into the argument list expected by the core parser."""
    argv: list[str] = []

    for key, value in config.items():
        if key in _BOOL_OPTIONS:
            if value:
                argv.append(_flag(key))
            continue

        if value is None:
            continue

        argv.append(_flag(key))
        if isinstance(value, (list, tuple)):
            argv.extend(str(item) for item in value)
        else:
            argv.append(str(value))

    return argv


def run(argv: list[str] | None = None) -> None:
    """Run the surface-pressure workflow.

    Explicit command-line arguments take precedence over ``CONFIG``. This makes
    the same entry point convenient for interactive use and for SLURM scripts.
    """
    if argv is None:
        argv = sys.argv[1:]

    if not argv:
        argv = build_argv(CONFIG)
        print("No CLI arguments supplied; using CONFIG in run_compute_surface_pressure.py.")

    # Keep the repository root importable when the package is not installed.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)

    from compute_surface_pressure_core import main

    main(argv)


if __name__ == "__main__":
    run()
