# Surface Pressure Core — AVBP FWH Post-Processing

Post-processing tools for computing surface-pressure spectral quantities
(PSD, CSD/coherence, FFT) and surface-line extractions from AVBP FWH
(Ffowcs Williams–Hawkings) solution data, for aeroacoustic surface source
localization on an airfoil.

The core logic lives in the `compute_surface_pressure_core` package. Two
standalone entry points (`run_compute_surface_pressure.py` and
`run_compute_surface_pressure.sh`) sit **outside** that package directory so
the workflow can be launched without caring about Python packaging details —
e.g. from a SLURM batch job.

```
.
├── compute_surface_pressure_core/   # the core package (importable module)
│   ├── __init__.py
│   ├── __main__.py                  # `python -m compute_surface_pressure_core ...`
│   ├── ComputeSurfacePressure.py    # CLI parser, SurfacePressure workflow class, main()
│   ├── extractor.py                 # file/mesh/pressure extraction from AVBP output
│   ├── surface_psd.py               # Welch PSD + ParaView-ready source export
│   ├── surface_csd.py               # Welch CSD/coherence + ParaView-ready source export
│   └── utils.py                     # FFT/SPOD helpers, timers, misc numeric helpers
├── run_compute_surface_pressure.py  # external Python entry point
├── run_compute_surface_pressure.sh  # external SLURM/shell entry point
└── README.md
```

## Requirements

- Python ≥ 3.9
- `numpy`, `scipy`, `h5py` (installed automatically via `pyproject.toml`)
- [`antares`](https://cerfacs.fr/antares/) — the Cerfacs CFD post-processing
  API used for mesh I/O (`Reader`, `Writer`, `Base`, `Zone`, `Instant`,
  `Treatment`, `Family`). This is **not** on PyPI and is not installed
  automatically — it must already be available in your environment
  (commonly provided as a module on HPC clusters, e.g. `module load antares`).

## Installation

From the repository root:

```bash
pip install -e .
```

This makes `compute_surface_pressure_core` importable from anywhere in that
Python environment, which is what both external entry-point scripts rely on.
(`run_compute_surface_pressure.py` also falls back to adding its own
directory to `sys.path`, so it will still work even without the `pip
install -e .` step, as long as it stays next to the
`compute_surface_pressure_core/` folder.)

## Usage

### As a script (outside the core package directory)

```bash
python run_compute_surface_pressure.py \
    --sol-dir /path/to/AVBP/FWH/solution \
    --mesh-file /path/to/mesh.mesh.h5 \
    --surface-patches Airfoil_Surface \
    --working-dir ./surface_pressure_output \
    --psd \
    --freq-select 500 1000 1500 2000 2500 3000
```

Or on a SLURM cluster, edit the variables at the top of
`run_compute_surface_pressure.sh` (or export them as environment variables
beforehand) and submit:

```bash
sbatch run_compute_surface_pressure.sh
```

### As a module

```bash
python -m compute_surface_pressure_core --sol-dir ... --mesh-file ... --psd
```

### As a library

```python
from compute_surface_pressure_core import parse_arguments, SurfacePressure

args = parse_arguments([
    "--sol-dir", "/path/to/solution",
    "--mesh-file", "/path/to/mesh.mesh.h5",
    "--psd", "--csd",
    "--mic-dir", "/path/to/mics", "--mic-file", "mic_data.hdf5", "--mic-num", "3",
])
sp = SurfacePressure(args)
sp.prepare_inputs()   # copies FWH files, extracts surface mesh, assembles pressure time series
sp.run_compute()      # runs whichever of --psd/--csd/--fft/--surf_line were requested
```

## Workflow options

| Flag | What it does |
|---|---|
| `--psd` | Welch PSD of surface pressure per node, band-averaged/summed at `--freq-select` frequencies, exported for ParaView via `source_psd`. |
| `--csd` | Welch CSD + magnitude-squared coherence between surface pressure and a reference far-field microphone (`--mic-dir/--mic-file/--mic-num`), exported via `source_csd`. |
| `--fft` | Computes the FFT/SPOD-style spectrum via `fft_surface_data`. **Not yet fully wired up** — see [Known gaps](#known-gaps) below. |
| `--surf_line` | Extracts a 1D spanwise or chordwise line of surface nodes (nearest-neighbor sampled) and its pressure time series, for two-point coherence/correlation analysis. |

## What was fixed in this pass

The package had several bugs that would have prevented it from importing or
running at all; these are now fixed:

- **Missing/incomplete imports.** Every module (`extractor.py`,
  `surface_psd.py`, `surface_csd.py`, `utils.py`,
  `ComputeSurfacePressure.py`) is now self-contained: each imports
  everything it actually uses (`os`, `glob`, `shutil`, `h5py`, `numpy`,
  `datetime`, `scipy.signal`, etc.) at module scope, and the two
  multiprocessing worker functions (`compute_PSD_block`,
  `compute_CSD_block`) additionally re-import their few dependencies
  locally so they remain safe to pickle/run under the `spawn` start method
  (macOS/Windows), not just `fork` (Linux default).
  `antares` is imported at the top of `extractor.py` (it's needed by nearly
  every function there) but lazily inside `source_psd`/`source_csd` in the
  spectral modules, since those are the only functions there that touch it.
- **`__init__.py`** referenced a non-existent `Surface_Pressure_Functions`
  module — removed, and the public API re-exported explicitly (`__all__`).
- **`argparse` syntax error**: a stray extra comma in the `--freq-select`
  argument definition, and `returns parser.parse_args(argv)` instead of
  `return ...`.
- **Literal `%` characters in three `--help` strings** (e.g. "Default is
  70%.") — these crash argparse's help formatter, since `%` is treated as a
  string-format character. Reworded to "70 percent".
- **`SurfacePressure` class**: `_print_args(self)` was being called as
  `self._print_args(self)` (passing `self` twice) and referenced a
  non-existent `self.args`; `run_compute` didn't take `self` as a
  parameter; `prepare_inputs` called `extract_surface` with arguments in
  the wrong order and referenced undefined local names (`FWH_data_dir`,
  `airfoil_mesh`) instead of `self.FWH_data_dir`/`self.airfoil_mesh`. All
  fixed.
- **`extractor.py`**: `extract_files` compared `Option` (capitalized,
  undefined) against the lowercase `option` parameter — always raised
  `NameError`. Fixed to use `option` consistently. Also, the chordwise
  surface-line branch had its airfoil/camber file default-path logic
  inverted (`if airfoil_file is not None else airfoil_file` meant a
  user-supplied path was silently thrown away and replaced by the
  hard-coded cluster path, while `None` inputs raised `FileNotFoundError`);
  this is now `is None`, so user-supplied paths are respected and only
  missing ones fall back to the default.
- **`utils.py`**: the `lowRAM` branch of `fft_surface_data` referenced an
  undefined `save_path`; it now derives an output directory from
  `surface_pressure_data`'s own path.

## Known gaps

- **`--fft` is not fully wired up.** `run_compute` calls
  `fft_surface_data`, but there is no `source_fft` counterpart to
  `source_psd`/`source_csd` for exporting the FFT result to a
  ParaView-ready surface file anywhere in the original codebase. Calling
  `--fft` will raise a clear `NotImplementedError` rather than silently
  doing nothing or fabricating an export routine — add a `surface_fft.py`
  module with a `source_fft(...)` function (mirroring `source_psd`) if you
  need this path.
- `extract_data`/`extract_surface` assume a specific AVBP/FWH directory and
  file-naming convention (`FWH_Airfoil_*.h5`, `frame_data/pressure`,
  `frame_data/time`) — adjust if your solver output differs.
- `extract_surface_line`'s chordwise-cut branch depends on user-supplied
  airfoil/camber coordinate `.txt` files if you don't pass `--airfoil-file`
  / `--camber-file` explicitly; the hard-coded cluster fallback paths in
  `extractor.py` are almost certainly specific to the original author's
  project storage and should be replaced with your own paths (or always
  pass the flags explicitly).