# DSMC_erlend

A Direct Simulation Monte Carlo (DSMC) solver for gas flow around an axis-aligned
box obstacle. The solver uses two spatial coordinates and three velocity components,
with MPI domain decomposition for parallel runs.

Particles move through a Cartesian grid, undergo binary elastic collisions, and
reflect diffusely from the obstacle surface. Sampled particle moments provide
number density, bulk velocity, and translational temperature fields.

## Requirements

- A C++20 compiler.
- CMake 3.20 or newer.
- An MPI implementation with C++ build support, such as Open MPI or MPICH.
- Python 3, NumPy, Matplotlib, and SciPy for plotting.

## Build

Run these commands from the project root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

The executable is `build/dsmc`. If CMake cannot find MPI, ensure its development
headers and libraries are installed. You can also specify the MPI compiler:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build --parallel
```

Use a fresh build directory when changing compilers.

## Run

Run with one MPI process:

```bash
mpiexec -n 1 ./build/dsmc
```

Or run with four processes:

```bash
mpiexec -n 4 ./build/dsmc
```

Set the angle of attack in degrees without rebuilding:

```bash
mpirun -np 4 ./build/dsmc --angle-of-attack 15
```

`--aoa` is an alias for `--angle-of-attack`. Accepted angles range from -180 to
180 degrees. At zero degrees the freestream flows toward positive x; positive
angles tilt it toward negative y. The box stays aligned with the grid. Use
`./build/dsmc --help` to see the options.

MPI automatically divides the grid into a two-dimensional process layout.
Progress messages report the time step, total particle count, and collision
count every 100 steps. Each process uses a time-based random seed, so repeated
runs need not produce identical results.

## Simulation settings

Time, gas, and sampling settings are defined near the top of
`src/Initialization.cpp`. Domain bounds, grid dimensions, and obstacle geometry
are defined in `src/Geometry.cpp`. Edit these settings and rebuild to change the
simulation. The flow angle can also be set with `--angle-of-attack`; other
parameters currently require editing the source.

| Setting | Default | Meaning |
| --- | --- | --- |
| `Dt` | `1e-5` s | Time step |
| `NStep` | `100000` | Number of time steps |
| `FirstSamplingStep` | `10000` | Sampling starts when `Step > FirstSamplingStep` |
| `SamplingPeriod` | `1` | Sample every this many steps after sampling starts |
| `NX`, `NY` | `100`, `100` | Global grid dimensions |
| `X1`, `X2`, `Y1`, `Y2` | `-20`, `20`, `-20`, `20` m | Domain bounds |
| `DZ` | `0.1` m | Thickness used to calculate cell volume |
| `NPCFree` | `10` | Target freestream simulation particles per cell |
| `MaFree`, `NFree`, `TFree` | `4`, `1e18` m⁻³, `200` K | Freestream Mach number, number density, and temperature |
| `BoxX`, `BoxY` | `0`, `0` m | Box centre |
| `BoxWidth`, `BoxHeight` | `2`, `2` m | Box dimensions along x and y |
| `AttackAngle` | `0` degrees | Default flow angle, overridden by `--angle-of-attack` |
| `Tw` | `300` K | Wall temperature |

Array capacities are defined in `include/Constants.h`. The solver retains shared
global state, declared in `include/SimulationState.h` and `include/Geometry.h`,
with definitions in `src/Initialization.cpp` and `src/Geometry.cpp`.

## Results

At the end of a run, rank 0 writes the gathered sampled field to:

```text
results/DSMC2D_Field.plt
```

The destination is the source project's `results/` directory, recorded by CMake
at configuration time. It does not depend on the working directory used to
launch the executable. The solver creates the directory if needed and overwrites
the field file on each completed run. Copy or rename the file to preserve a run.

`results/DSMC2D_Field.json` records the box geometry and angle of attack. Keep
this file alongside the field when archiving a run, using the same basename
(for example, `run1.plt` and `run1.json`). Field plots read this metadata
automatically. The Knudsen number reported during setup uses `BoxHeight` as
the reference length.

The file contains two header lines followed by these columns:

| Column | Quantity | Unit |
| --- | --- | --- |
| `X`, `Y` | Cell-centre coordinates | m |
| `N` | Number density | m⁻³ |
| `Ux`, `Uy` | Mean velocity components | m/s |
| `Tx`, `Ty` | Directional translational temperatures | K |
| `T` | Mean translational temperature over all three velocity components | K |

Field output is currently written only at the end of the run. `PrintPeriod` does
not enable intermediate field files in the current time loop.

## Plotting

Create a Python environment and install the plotting dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r postprocessing/requirements.txt
```

Display density, temperature, pressure, and velocity fields:

```bash
python postprocessing/plot_field.py
```

Display temperature and pressure along the grid column closest to a chosen
`x` coordinate:

```bash
python postprocessing/plot_x.py --x -5
```

Both scripts default to the project's `results/DSMC2D_Field.plt`. To choose a
different run and save PNG files instead of displaying figures:

```bash
python postprocessing/plot_field.py results/run1/DSMC2D_Field.plt --output-dir results/run1/plots
python postprocessing/plot_x.py results/run1/DSMC2D_Field.plt --x -5 --output-dir results/run1/slices
```

On a machine without a graphical display, prefix a plotting command with
`MPLBACKEND=Agg` and use `--output-dir`.

The field plot masks the obstacle interior using the saved JSON metadata without
drawing a box overlay. If metadata
is absent, it assumes the default 6 m × 6 m box at the origin. Old circular-obstacle
results should be plotted with the earlier plotting script; generate a new run
for the box geometry.
Use `--help` on either script to see its options.

## Tests

The build includes checks for face and corner intersections, grazing trajectories,
initial particle exclusion, flow-angle calculations, and outgoing wall reflections:

```bash
ctest --test-dir build --output-on-failure
```

## Project layout

```text
DSMC_erlend/
├── CMakeLists.txt
├── include/          # Data types, constants, module interfaces, shared state
├── src/              # Entry point and solver implementations
├── rng/              # Random number generators
├── postprocessing/   # Python plotting scripts and requirements
├── tests/            # Geometry and wall-reflection checks
├── results/          # Generated field data and plots
├── build/            # Generated build files and executable
└── .vscode/          # IntelliSense configuration
```

| Module | Responsibility |
| --- | --- |
| `Simulation` | MPI lifecycle and time loop |
| `Initialization` | Time/gas/sampling parameters, derived gas properties, and initial conditions |
| `Boundary` | Particle injection, physical boundaries, motion, and wall reflection |
| `Geometry` | Domain/grid/obstacle settings, cell dimensions, indexing, and obstacle exclusion |
| `Collision` | Binary collision selection and elastic scattering |
| `Sampling` | Accumulate moments and write field output |
| `MPIManager` | Domain decomposition, particle exchange, and field gathering |

The mathematical Gamma function used during setup is a private helper in
`src/Initialization.cpp`.

## VS Code

CMake generates `build/compile_commands.json`. The checked-in IntelliSense
configuration uses it to find the same headers, MPI include paths, and compiler
definitions as the build. Configure the project with CMake before opening source
files. If stale include errors remain, run **Developer: Reload Window** from the
Command Palette.

## Version control

Commit source files, plotting scripts, dependency lists, documentation, and build
configuration. `build/`, generated files under `results/`, `.venv/`, Python
caches, and `.save` backups are ignored. `results/.gitkeep` retains an empty
results directory in a fresh checkout.
