# DSMC_erlend

A Direct Simulation Monte Carlo (DSMC) solver for gas flow around an axis-aligned
box obstacle. The solver uses two spatial coordinates and three velocity components,
with MPI domain decomposition for parallel runs.

Particles move through a Cartesian grid, undergo binary elastic collisions, and
reflect diffusely or specularly from the obstacle surface. Sampled particle moments provide
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

Choose the wall interaction near the top of `src/Boundary.cpp`, then rebuild:

```cpp
WallModel WallInteraction = WallModel::Specular;
```

Use `WallModel::Diffuse` to select diffuse reflection. Diffuse reflection samples outgoing
velocities at the wall temperature `Tw`. Specular reflection from the stationary
box reverses the normal velocity and preserves tangential velocity and kinetic
energy; it does not use `Tw`.

For Cercignani–Lampis–Lord (CLL) reflection, select and configure it in
`src/Boundary.cpp`, then rebuild:

```cpp
WallModel WallInteraction = WallModel::CLL;
double CLLNormalAccommodation = 0.8;
double CLLTangentialAccommodation = 0.8;
```

Both accommodation coefficients range from 0 to 1 and use SPARTA's `acc_n` and
`acc_t` convention. `(0, 0)` gives specular reflection, while `(1, 1)` gives the
fully accommodated diffuse distribution at `Tw`. Intermediate values retain
part of the incident velocity distribution. The tangential coefficient is an
energy accommodation parameter; if using a tangential momentum coefficient
`sigma_t` in [0, 1], convert it with `acc_t = sigma_t * (2 - sigma_t)`. The solver models
translational motion of a monatomic gas, so no rotational or vibrational
accommodation is included. See [SPARTA's wall-model documentation](https://sparta.github.io/doc/surf_collide.html).

To save both cases, run with `WallModel::Specular`, then change the setting to
`WallModel::Diffuse`, rebuild, and run the same `dsmc` executable again.

MPI automatically divides the grid into a two-dimensional process layout.
Rank 0 prints setup information and progress messages. Progress reports include
the time step, total particle count, and collision count every 100 steps.
Errors identify the process where they occur. Each process uses a time-based random seed, so repeated
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
| `MaFree`, `NFree`, `TFree` | `4`, `1e17` m⁻³, `200` K | Freestream Mach number, number density, and temperature |
| `BoxX`, `BoxY` | `0`, `0` m | Box centre |
| `BoxWidth`, `BoxHeight` | `4`, `4` m | Box dimensions along x and y |
| `AttackAngle` | `0` degrees | Default flow angle, overridden by `--angle-of-attack` |
| `Tw` | `300` K | Wall temperature |

Array capacities are defined in `include/Constants.h`. The solver retains shared
global state, declared in `include/SimulationState.h` and `include/Geometry.h`,
with definitions in `src/Initialization.cpp` and `src/Geometry.cpp`.

## Results

At the end of each run, rank 0 saves the field and metadata according to the
wall model selected in `Boundary.cpp`:

```text
results/
├── specular/
│   ├── DSMC2D_Field.plt
│   └── DSMC2D_Field.json
├── diffuse/
│   ├── DSMC2D_Field.plt
│   └── DSMC2D_Field.json
└── cll/
    ├── DSMC2D_Field.plt
    └── DSMC2D_Field.json
```

A run replaces only the files for its selected wall model. The solver creates
the folders automatically. Results
always go into the source project's `results/` folder, regardless of where the
executable is launched.

The JSON file records geometry, wall model, flow angle, and simulation settings.
Keep it alongside the field for plotting. The Knudsen number printed during setup
uses `BoxHeight` as the reference length.

CLL metadata also records both accommodation coefficients.

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

Save density, temperature, pressure, and velocity plots to `results/plots/`:

```bash
python postprocessing/plot_field.py
```

Save temperature and pressure plots along the grid column closest to a chosen
`x` coordinate:

```bash
python postprocessing/plot_x.py --x -5
```

After completing both cases, compare them with:

```bash
python postprocessing/compare_runs.py
```

The script reads the field from each model folder and overlays temperature,
pressure, density, and speed. Plots are saved in `results/plots/comparison/`.
If a CLL result exists in `results/cll/`, it is included automatically with its
accommodation coefficients in the legend.
Change `PROFILE_AXIS` and `PROFILE_POSITION` near the top of `compare_runs.py`
to choose the line. The default is along x at y = 0; use `PROFILE_AXIS = "y"`
and `PROFILE_POSITION = -5.0` for a vertical profile at x = -5.

For single-field plots, select `WALL_MODEL = "specular"`, `"diffuse"`, or `"cll"` near
the top of `postprocessing/plot_cli.py`. No path arguments are needed.

The comparison also reads `sparta_y_x_-5_specular.plt`,
`sparta_y_x_-5_diffuse.plt`, and `sparta_y_x_-5_cll.plt` from the project root,
with separate SPARTA labels for each wall model. The exports must have the same eight field columns
listed above, in SI units, and two header lines. For this y-profile at x = -5 m,
set `PROFILE_AXIS = "y"` and `PROFILE_POSITION = -5` in `compare_runs.py`.
Set an entry in `SPARTA_FILES` to `None` to omit that reference, or edit its path
for another file. Replace a file with a new SPARTA export to update the comparison.

SPARTA values are plotted without density rescaling. Its export has no simulation
metadata, so match freestream density, geometry, angle, wall model, and sampling
settings yourself before interpreting agreement with the DSMC results.

Plots use a noninteractive backend, so saving also works when desktop display
access is unavailable. You can override the default output directory:

```bash
python postprocessing/plot_field.py --output-dir results/plots
```

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

With the Python plotting dependencies installed, check profile extraction and
comparison validation:

```bash
python -m unittest discover -s tests -p 'test_*.py'
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
