# Surface Pressure Core — AVBP FWH Post-Processing

A Python post-processing toolkit for extracting, processing, and analyzing unsteady surface-pressure data from AVBP Ffowcs Williams–Hawkings (FWH) output.

The package is intended primarily for aeroacoustic analysis and surface-source localization on airfoils and wings. It provides utilities for:

- extracting transient FWH surface-pressure files;
- extracting selected surface patches from an AVBP mesh;
- assembling nodal surface-pressure time histories into HDF5 datasets;
- computing surface-pressure Power Spectral Density (PSD);
- computing surface-to-microphone Cross-Spectral Density (CSD) and coherence;
- exporting frequency-dependent surface quantities to ANTARES HDF5 files for visualization in ParaView;
- computing FFT-based quantities;
- extracting spanwise or chordwise surface-pressure lines for coherence, correlation, and space-time analysis.

The computational logic is contained in the `compute_surface_pressure_core` package. Two external run scripts are kept outside the package:

- `run_compute_surface_pressure.py` — user-facing Python runner;
- `run_compute_surface_pressure.sh` — shell/SLURM launcher.

This separation allows the package internals to remain reusable while simulation-specific settings can be defined in the external run scripts.

---

# 1. Repository Structure

```text
.
├── compute_surface_pressure_core/
│   ├── __init__.py
│   ├── __main__.py
│   ├── ComputeSurfacePressure.py
│   ├── extractor.py
│   ├── surface_psd.py
│   ├── surface_csd.py
│   └── utils.py
├── run_compute_surface_pressure.py
├── run_compute_surface_pressure.sh
└── README.md
```

## Core files

### `ComputeSurfacePressure.py`
Contains the CLI parser, `SurfacePressure` workflow class, input preparation, processing dispatch, and `main()` entry point.

```text
parse arguments
      ↓
SurfacePressure(args)
      ↓
prepare_inputs()
      ↓
run_compute()
```

### `extractor.py`
Provides:

```python
extract_files(...)
extract_surface(...)
extract_data(...)
extract_surface_line(...)
```

These functions copy FWH files, extract AVBP surface patches, assemble pressure time histories, and optionally extract one-dimensional surface lines.

### `surface_psd.py`
Provides:

```python
PSD_surface_data(...)
compute_PSD_block(...)
source_psd(...)
```

It computes nodal surface-pressure PSDs with Welch's method and exports selected spectral quantities back onto the surface mesh.

### `surface_csd.py`
Provides CSD/coherence calculations between the surface-pressure field and a reference microphone.

### `utils.py`
Contains common numerical helpers including FFT/SPOD utilities, timers, window functions, and helper routines.

### `__init__.py`
Defines the package public API.

### `__main__.py`
Allows:

```bash
python -m compute_surface_pressure_core ...
```

### `run_compute_surface_pressure.py`
Top-level user runner. It lives outside the core package and calls `compute_surface_pressure_core.main(...)`.

### `run_compute_surface_pressure.sh`
Shell/SLURM wrapper for launching the Python runner on HPC systems.

---

# 2. Requirements

## Python
Python 3.8 or newer is recommended. On HPC systems, compatibility with the installed ANTARES version may determine the exact Python version.

## Python packages

```text
numpy
scipy
h5py
```

The standard library modules used include `argparse`, `datetime`, `glob`, `logging`, `multiprocessing`, `os`, and `shutil`.

## ANTARES
The package requires CERFACS ANTARES for AVBP mesh I/O, surface extraction, and output writing.

Objects used include:

```python
Reader
Writer
Base
Zone
Instant
Treatment
Family
```

ANTARES is normally supplied through the CFD/HPC environment rather than PyPI.

---

# 3. Expected AVBP FWH Data

Transient files are expected to follow a convention similar to:

```text
FWH_Airfoil_00000001.h5
FWH_Airfoil_00000002.h5
...
```

Each file is expected to contain:

```text
frame_data/
├── pressure
└── time
```

The extractor accesses:

```python
f["frame_data/pressure"]
f["frame_data/time"]
```

If the AVBP output uses a different naming convention or HDF5 hierarchy, modify `extractor.py`.

---

# 4. Typical Input and Output Layout

A typical case may look like:

```text
Case/
├── FWH/
│   ├── FWH_Airfoil_00000001.h5
│   ├── ...
│   └── FWH_Airfoil_00004925.h5
├── mesh.mesh.h5
└── Compute_Surface_Pressure/
    ├── compute_surface_pressure_core/
    ├── run_compute_surface_pressure.py
    ├── run_compute_surface_pressure.sh
    └── README.md
```

With:

```text
--working-dir ./surface_pressure_output
```

a typical output directory becomes:

```text
surface_pressure_output/
├── FWH_Data/
├── Airfoil_Surface_Mesh.h5
├── pressure_airfoil.hdf5
├── pressure_airfoil_psd.hdf5
├── pressure_airfoil_csd.hdf5
├── Surface_psd.hdf5
└── ...
```

---

# 5. Installation and Import

The simplest workflow does not require package installation as long as `run_compute_surface_pressure.py` remains beside `compute_surface_pressure_core/`.

Run:

```bash
python run_compute_surface_pressure.py
```

If packaging metadata are available, an editable installation can be used:

```bash
pip install -e .
```

Verify the import with:

```bash
python -c "import compute_surface_pressure_core; print(compute_surface_pressure_core.__file__)"
```

---

# 6. Running the Code

## External runner

```bash
python run_compute_surface_pressure.py     --sol-dir /path/to/FWH     --mesh-file /path/to/mesh.mesh.h5     --surface-patches Airfoil_Surface     --working-dir ./surface_pressure_output     --psd     --freq-select 500 1000 1500 2000 2500 3000
```

## SLURM

```bash
sbatch run_compute_surface_pressure.sh
```

## Package module

```bash
python -m compute_surface_pressure_core     --sol-dir /path/to/FWH     --mesh-file /path/to/mesh.mesh.h5     --psd
```

## Programmatic use

```python
from compute_surface_pressure_core import parse_arguments, SurfacePressure

args = parse_arguments([
    "--sol-dir", "/path/to/FWH",
    "--mesh-file", "/path/to/mesh.mesh.h5",
    "--surface-patches", "Airfoil_Surface",
    "--working-dir", "./surface_pressure_output",
    "--psd",
    "--freq-select", "500", "1000", "1500", "2000",
])

surface_pressure = SurfacePressure(args)
surface_pressure.prepare_inputs()
surface_pressure.run_compute()
```

---

# 7. Main Workflow

```text
Original AVBP FWH files
          │
          ▼
    extract_files()
          │
          ▼
working_dir/FWH_Data/
          │
AVBP mesh ─────► extract_surface()
                    │
                    ▼
           Airfoil_Surface_Mesh.h5
                    │
FWH_Data/ ──────────┤
                    ▼
               extract_data()
                    │
                    ▼
           pressure_airfoil.hdf5
                    │
        ┌───────────┼─────────────┐
        │           │             │
        ▼           ▼             ▼
       PSD         CSD           FFT
        │           │
        ▼           ▼
 source_psd()   source_csd()
        │           │
        ▼           ▼
ParaView-ready surface quantities
```

Surface-line extraction uses the extracted mesh and pressure database:

```text
Airfoil_Surface_Mesh.h5
          +
pressure_airfoil.hdf5
          │
          ▼
extract_surface_line()
          │
          ▼
1-D surface-pressure time histories
```

---

# 8. Main CLI Arguments

| Argument | Default | Description |
|---|---:|---|
| `--sol-dir` | required | Directory containing AVBP FWH files. |
| `--mesh-file` | required | AVBP mesh file. |
| `--working-dir` | current directory | Output/intermediate directory. |
| `--surface-patches`, `-sp` | `Airfoil_Surface` | One or more AVBP patch-family names. |
| `--option` | `1` | File extraction strategy. |
| `--nskip` | `1` | File skipping interval. |
| `--max-file` | `5000` | Maximum number of FWH files. |
| `--reload` | false | Force regeneration instead of reuse. |
| `--var` | `pressure` | Variable to process. |
| `--psd` | false | Run PSD workflow. |
| `--csd` | false | Run CSD/coherence workflow. |
| `--fft` | false | Run FFT workflow. |
| `--surf_line` | false | Extract a surface line. |
| `--block-size` | `1000` | Number of surface nodes per processing block. |
| `--nchunk` | `4` | Number of Welch/FFT chunks. |
| `--df` | `50` | PSD source-map band width in Hz. |
| `--band-stat` | `mean` | `mean` or `sum` over selected PSD band. |
| `--freq-select` | case-dependent | Selected frequencies for exported maps. |

CSD additionally uses:

```text
--mic-dir
--mic-file
--mic-num
```

Surface-line extraction additionally uses:

```text
--AoA
--Uinf
--orientation
--cut-loc-percent
--z-loc-percent
--cut-length-percent
--airfoil-file
--camber-file
--data-size
```

---

# 9. FWH File Extraction

`extract_files()` copies the requested FWH files into:

```text
<working_dir>/FWH_Data/
```

This reduces repeated I/O against the original simulation directory.

### Option 1

```bash
--option 1
```

Processes files sequentially up to `--max-file`.

### Option 2

```bash
--option 2 --nskip N
```

Processes approximately every `N`th file.

Skipping snapshots changes the effective temporal sampling. Ensure the spectral calculation uses the correct physical sample spacing.

---

# 10. Surface Mesh Extraction

The function signature is:

```python
extract_surface(
    mesh_file,
    input_surface,
    working_dir,
    reload=False,
)
```

The argument order matters.

Correct:

```python
airfoil_mesh = extract_surface(
    "/path/to/mesh.mesh.h5",
    ["Airfoil_Surface"],
    "./surface_pressure_output",
)
```

The function returns the path to:

```text
Airfoil_Surface_Mesh.h5
```

Downstream routines should use that returned path directly:

```python
reader["filename"] = airfoil_mesh
```

Do not prepend `working_dir` again, or duplicated paths can result.

---

# 11. Surface-Pressure Database

`extract_data()` assembles the FWH snapshots into:

```text
pressure_airfoil.hdf5
```

A typical file contains:

```text
pressure
mean_pressure
rms_pressure
min_pressure
max_pressure
```

and attributes including:

```text
dt
Extracted Date
Source Path
Mesh Path
```

The pressure dataset is organized as:

```text
(n_nodes, n_time)
```

---

# 12. PSD Workflow

Enable with:

```bash
--psd
```

The workflow is:

```text
pressure_airfoil.hdf5
        ↓
PSD_surface_data()
        ↓
pressure_airfoil_psd.hdf5
        ↓
source_psd()
        ↓
Surface_psd.hdf5
```

## Welch calculation

For each node:

```python
f, Spp = scipy.signal.welch(
    pressure,
    fs=1.0 / dt,
    window="hann",
    nperseg=nperseg,
    nfft=nfft,
    scaling="density",
)
```

where approximately:

```python
nperseg = n_time / nchunk
```

and:

```python
nfft = next_greater_power_of_2(nperseg)
```

For pressure in Pa, the PSD has units of approximately:

```text
Pa²/Hz
```

## Parallelization

Nodes are divided into blocks controlled by:

```bash
--block-size 1000
```

The worker pool processes spatial blocks in parallel using Python `multiprocessing`.

---

# 13. PSD Output

The full PSD is stored in:

```text
pressure_airfoil_psd.hdf5
```

with:

```text
frequency
pressure_psd
```

where:

```text
pressure_psd.shape = (n_nodes, n_frequency_bins)
```

Metadata include the computation date, reference signal, sampling frequency, node count, and time-step count.

---

# 14. Reusing Existing PSD Results

If `pressure_airfoil_psd.hdf5` already exists and `reload=False`, the file should be reused.

Read an existing file with:

```python
with h5py.File(psd_file_path, "r") as h5f:
    ...
```

Never inspect an existing HDF5 file with:

```python
h5py.File(psd_file_path, "w")
```

because `"w"` truncates the existing file.

A valid cached PSD file should contain:

```python
"frequency" in h5f
"pressure_psd" in h5f
```

If the file exists but these datasets are missing, the PSD should be recomputed.

The global reload value should be propagated into the PSD calculation:

```python
surface_pressure_psd_data = PSD_surface_data(
    self.surface_pressure_data,
    self.var,
    self.dt,
    reload=self.reload,
    block_size=self.psd.block_size,
    nchunk=self.psd.nchunk,
)
```

---

# 15. Frequency Selection and Band Averaging

`--freq-select` controls which frequencies are mapped back onto the surface.

Example:

```bash
--freq-select 500 1000 1500 2000 2500 3000
```

The full PSD is computed first; selected bands are extracted afterward by `source_psd()`.

For center frequency `f0` and bandwidth `df`:

```text
f0 - df/2 ≤ f ≤ f0 + df/2
```

Example:

```bash
--freq-select 2000 --df 100
```

selects approximately 1950–2050 Hz.

Two statistics are supported:

```bash
--band-stat mean
```

and:

```bash
--band-stat sum
```

The logarithmic level is computed using:

```python
10.0 * np.log10(Spp / pref**2)
```

with:

```python
pref = 2e-5
```

by default.

Because PSD has units of Pa²/Hz while a frequency-integrated mean-square pressure has units of Pa², interpret `mean` and `sum` according to the chosen spectral discretization.

---

# 16. Surface PSD Export

`source_psd()` maps selected PSD bands back onto the extracted surface mesh and writes:

```text
Surface_psd.hdf5
```

Fields may include:

```text
frequency_1000_Hz_P_band_df50Hz_mean_dB
frequency_1000_Hz_P_band_df50Hz_mean_Spp
Pressure_Mean
Pressure_RMS
Pressure_Min
```

The output can be visualized through ANTARES/ParaView-compatible workflows.

---

# 17. CSD and Coherence

Enable with:

```bash
--csd
```

Example:

```bash
python run_compute_surface_pressure.py     --sol-dir /path/to/FWH     --mesh-file /path/to/mesh.mesh.h5     --working-dir ./surface_pressure_output     --csd     --mic-dir /path/to/microphone/results     --mic-file microphone_data.hdf5     --mic-num 3
```

CSD/coherence measures the frequency-domain statistical relationship between each surface node and the selected reference microphone.

This differs from PSD:

- PSD indicates where surface-pressure fluctuations are strong.
- CSD/coherence indicates which surface fluctuations are statistically related to the reference acoustic signal.

A large hydrodynamic surface PSD does not necessarily imply strong far-field acoustic coupling.

---

# 18. FFT Mode

Enable with:

```bash
--fft
```

The calculation uses:

```python
fft_surface_data(...)
```

However, the current codebase does not provide a complete `source_fft(...)` export routine analogous to `source_psd()` and `source_csd()`.

FFT-to-surface visualization should therefore be considered only partially integrated.

---

# 19. Surface-Line Extraction

Enable with:

```bash
--surf_line
```

`extract_surface_line()` generates a one-dimensional line of mesh nodes and extracts their complete pressure time histories.

Potential applications include:

- two-point coherence;
- spanwise coherence length;
- chordwise coherence;
- space-time correlation;
- convection velocity;
- hydrodynamic wavenumber analysis.

Supported orientations are:

### Spanwise

```bash
--orientation 1 0 0
```

### Chordwise

```bash
--orientation 0 0 1
```

The routine constructs requested sampling positions, builds a `scipy.spatial.cKDTree` from the original mesh nodes, finds the nearest physical mesh nodes, and extracts the corresponding time histories.

Chordwise extraction may require:

```bash
--airfoil-file /path/to/airfoil_coordinates.txt --camber-file /path/to/camber_coordinates.txt
```

Explicit paths are recommended because any hard-coded fallback paths are project-specific.

---

# 20. Example PSD Run

```bash
python run_compute_surface_pressure.py     --sol-dir /scratch/user/Case/FWH     --mesh-file /scratch/user/Case/mesh.mesh.h5     --surface-patches Airfoil_Surface     --working-dir ./surface_pressure_output     --psd     --freq-select 500 1000 1500 2000 2500 3000     --df 50     --band-stat mean     --block-size 1000     --nchunk 4
```

---

# 21. Example PSD + CSD Run

```bash
python run_compute_surface_pressure.py     --sol-dir /path/to/FWH     --mesh-file /path/to/mesh.mesh.h5     --surface-patches Airfoil_Surface     --working-dir ./surface_pressure_output     --psd     --csd     --mic-dir /path/to/microphone/results     --mic-file microphone_data.hdf5     --mic-num 3     --freq-select 1000 1500 2000 3000 5000     --block-size 1000     --nchunk 4
```

---

# 22. Example Surface-Line Runs

## Spanwise

```bash
python run_compute_surface_pressure.py     --sol-dir /path/to/FWH     --mesh-file /path/to/mesh.mesh.h5     --surface-patches Airfoil_Surface     --working-dir ./surface_pressure_output     --surf_line     --AoA 10     --Uinf 30     --orientation 1 0 0     --cut-loc-percent 70     --z-loc-percent 70     --cut-length-percent 10     --data-size 1000
```

## Chordwise

```bash
python run_compute_surface_pressure.py     --sol-dir /path/to/FWH     --mesh-file /path/to/mesh.mesh.h5     --surface-patches Airfoil_Surface     --working-dir ./surface_pressure_output     --surf_line     --AoA 10     --Uinf 30     --orientation 0 0 1     --z-loc-percent 70     --cut-loc-percent 50     --cut-length-percent 40     --data-size 1000     --airfoil-file /path/to/airfoil_10_AOA.txt     --camber-file /path/to/airfoil_camber_10_AOA.txt
```

---

# 23. Example SLURM Script

```bash
#!/bin/bash
#SBATCH --job-name=surface_psd
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=0
#SBATCH --time=04:00:00
#SBATCH --output=surface_psd_%j.out
#SBATCH --error=surface_psd_%j.err

module purge

# Load project-specific modules/environment here.
source /path/to/python/environment/bin/activate

export OMP_NUM_THREADS=1

python run_compute_surface_pressure.py     --sol-dir /path/to/FWH     --mesh-file /path/to/mesh.mesh.h5     --surface-patches Airfoil_Surface     --working-dir ./surface_pressure_output     --psd     --freq-select 1000 1500 2000 2500 3000     --block-size 1000     --nchunk 4
```

`OMP_NUM_THREADS=1` can help prevent oversubscription when Python multiprocessing is already being used.

---

# 24. HPC Memory Considerations

The PSD workflow currently loads the complete pressure matrix:

```python
data = h5f[var][:]
```

For double precision, approximate raw storage is:

```text
8 × n_nodes × n_time bytes
```

Example:

```text
500,000 × 5,000 × 8 bytes ≈ 20 GB
```

before accounting for worker copies, PSD arrays, FFT workspace, and Python overhead.

Large cases may therefore require substantial memory.

Potential future improvements include chunked HDF5 reads, shared memory, memory-mapped arrays, node-block streaming, and explicit worker limits.

---

# 25. Multiprocessing on SLURM

The current implementation may use:

```python
multiprocessing.cpu_count()
```

This can report more CPUs than allocated by the scheduler.

A future improvement should expose an `--nproc` argument or infer the process count from:

```text
SLURM_CPUS_PER_TASK
```

For current runs, ensure the scheduler allocation is consistent with the worker count.

---

# 26. File-Path Conventions

Functions should return usable paths, and downstream routines should use those paths directly.

For example:

```python
self.airfoil_mesh = extract_surface(...)
reader["filename"] = self.airfoil_mesh
```

Avoid:

```python
os.path.join(self.working_dir, self.airfoil_mesh)
```

if `self.airfoil_mesh` already contains the working directory.

Similarly, output files should be written explicitly into the intended directory:

```python
surface_pressure_data = os.path.join(
    working_dir,
    "pressure_airfoil.hdf5",
)
```

---

# 27. HDF5 File Modes

### Read-only

```python
h5py.File(filename, "r")
```

Use for existing results.

### Write

```python
h5py.File(filename, "w")
```

Creates a new file or truncates an existing one.

### Append/read-write

```python
h5py.File(filename, "a")
```

Opens an existing file for modification or creates it if necessary.

Never use `"w"` to inspect an existing cache file.

---

# 28. Restartability

The intended behavior is:

```text
First run:
copy FWH files
extract surface mesh
assemble pressure data
compute PSD
export surface PSD

Later run:
reuse FWH files
reuse surface mesh
reuse pressure database
reuse full PSD
rerun only the selected source-map export
```

This is particularly useful when changing only:

```text
--freq-select
--df
--band-stat
```

because the full PSD does not need to be recomputed.

---

# 29. Debugging and Validation

## Syntax check

```bash
python -m compileall     compute_surface_pressure_core     run_compute_surface_pressure.py
```

## Import check

```bash
python - <<'PY'
import compute_surface_pressure_core
from compute_surface_pressure_core import SurfacePressure, parse_arguments, main
print("Package import successful")
PY
```

## Argument parser check

```bash
python run_compute_surface_pressure.py --help
```

## Inspect an HDF5 PSD file

```bash
python - <<'PY'
import h5py

filename = "./surface_pressure_output/pressure_airfoil_psd.hdf5"

with h5py.File(filename, "r") as h5f:
    print("Datasets:", list(h5f.keys()))
    print("Attributes:")
    for key, value in h5f.attrs.items():
        print(key, "=", value)
PY
```

A valid PSD file should contain:

```text
frequency
pressure_psd
```

---

# 30. Common Errors

## `ImportError: cannot import name 'main'`

Ensure `ComputeSurfacePressure.py` defines:

```python
def main(argv=None):
    args = parse_arguments(argv)
    surface_pressure = SurfacePressure(args)
    surface_pressure.prepare_inputs()
    surface_pressure.run_compute()
```

## `TypeError: expected str, bytes or os.PathLike object, not list`

Usually caused by incorrect `extract_surface()` argument order.

Correct:

```python
extract_surface(mesh_file, surface_patches, working_dir)
```

## Duplicated mesh path

Example:

```text
./surface_pressure_output/./surface_pressure_output/Airfoil_Surface_Mesh.h5
```

Use the returned `airfoil_mesh` path directly rather than joining it to `working_dir` again.

## `KeyError: object 'pressure_psd' doesn't exist`

A previously computed PSD file may have been accidentally opened with mode `"w"` and truncated.

Use `"r"` when reusing existing files.

## `NameError: Option is not defined`

Use the lowercase function argument consistently:

```python
if option == 1:
```

not:

```python
if Option == 1:
```

## Surface patch not found

The patch name must exactly match an AVBP family under the mesh's patch families.

## Surface-pressure node mismatch

The number of pressure values in each FWH snapshot must match the number of nodes in the extracted surface mesh.

Possible causes include an incorrect mesh, wrong patch selection, inconsistent FWH generation mesh, or corrupted output.

---

# 31. Physical Interpretation

`Pressure_Mean` is the time-average surface pressure.

`Pressure_RMS` measures total unsteady pressure amplitude.

`Spp(f)` describes the frequency distribution of surface-pressure fluctuations.

CSD measures frequency-dependent statistical coupling between a surface node and a reference microphone.

Magnitude-squared coherence normalizes the cross-spectrum and measures the degree of linear frequency-domain correlation.

These quantities should not be interpreted interchangeably.

Large wall-pressure PSD does not automatically imply efficient far-field acoustic radiation. Hydrodynamic pressure fluctuations from turbulent boundary layers, shear layers, and vortical structures can dominate the surface field while coupling only weakly to the acoustic far field.

---

# 32. Recommended Validation Checks

Before using results quantitatively, verify:

- the FWH sampling interval `dt`;
- Nyquist frequency;
- total record duration;
- Welch segment length;
- FFT-bin spacing;
- number of statistical averages;
- spectral convergence;
- mesh-node consistency;
- microphone sampling rate;
- synchronization of microphone and surface signals;
- selected surface patches.

The Nyquist frequency is:

```text
f_N = 1 / (2 dt)
```

The full-record frequency scale is approximately:

```text
1 / T
```

although the actual Welch resolution is governed by the segment length and FFT size.

---

# 33. Known Limitations

- FFT surface export is not fully implemented because a complete `source_fft()` routine is absent.
- File naming currently assumes `FWH_Airfoil_*.h5`.
- HDF5 paths currently assume `frame_data/pressure` and `frame_data/time`.
- Surface pressure is assumed to map directly to the extracted surface mesh nodes.
- Large cases can require substantial RAM.
- Chordwise line extraction contains geometry-specific assumptions.
- Hard-coded fallback coordinate paths are not portable.
- The package does not automatically establish physical or statistical convergence.

---

# 34. Development Notes

Important fixes made while consolidating the codebase include:

- invalid `argparse` syntax;
- `returns` changed to `return`;
- addition of the missing `main()` entry point;
- correction of `extract_surface()` argument order;
- correction of class-method `self` handling;
- storage of `self.args`;
- replacement of undefined locals by `self.FWH_data_dir` and `self.airfoil_mesh`;
- `Option` changed to `option`;
- correction of duplicated working-directory paths;
- explicit writing of output HDF5 files to the configured working directory;
- prevention of accidental PSD-cache truncation by opening existing files with `"r"` rather than `"w"`;
- propagation of `self.reload` to the PSD workflow;
- correction of airfoil/camber fallback-path logic;
- correction of undefined utility output paths.

Because this originated as research post-processing software rather than a general-purpose library, new geometries and workflows should be validated carefully.

---

# 35. Potential Future Improvements

Potential improvements include:

- expose multiprocessing count through `--nproc`;
- detect `SLURM_CPUS_PER_TASK`;
- stream pressure data by spatial blocks;
- use shared-memory arrays;
- add HDF5 compression and chunking;
- implement `source_fft()`;
- separate extraction and spectral reload flags;
- expose Welch window and overlap settings;
- save all runtime arguments into output metadata;
- save actual spectral resolution;
- validate microphone and surface sampling consistency automatically;
- remove hard-coded geometry/project paths;
- add unit and integration tests;
- improve structured logging and exception reporting;
- support alternate surface variables and naming conventions;
- implement physically explicit band integration using frequency-bin widths.

---

# 36. Quick Reference

Typical PSD run:

```bash
python run_compute_surface_pressure.py     --sol-dir /path/to/FWH     --mesh-file /path/to/mesh.mesh.h5     --surface-patches Airfoil_Surface     --working-dir ./surface_pressure_output     --psd     --freq-select 1000 1500 2000 2500 3000     --df 50     --band-stat mean     --block-size 1000     --nchunk 4
```

Major outputs:

```text
surface_pressure_output/
├── FWH_Data/
├── Airfoil_Surface_Mesh.h5
├── pressure_airfoil.hdf5
├── pressure_airfoil_psd.hdf5
└── Surface_psd.hdf5
```

For later changes to visualization frequencies, retain `pressure_airfoil_psd.hdf5` and use `reload=False` so only the surface-map export is repeated.
