# Z-Pinch Stability Simulator: MHD Plasma Analysis

![C++17](https://img.shields.io/badge/Language-C++17-blue)  ![Python](https://img.shields.io/badge/Visualization-Python-green)
![Plasma Physics](https://img.shields.io/badge/Physics-Plasma_Physics-red)

---

## Description

This project implements a **Z-Pinch plasma stability simulator** for **MHD analysis**. It supports **equilibrium calculations, linear stability analysis, and nonlinear evolution** of plasma instabilities such as **kink and sausage modes**.

Key features:

* Equilibrium force-balance solver
* Linear stability analysis for kink and sausage modes
* Nonlinear time evolution of plasma instabilities
* Full analysis pipeline combining equilibrium, linear, and nonlinear simulations
* Python visualization scripts for field evolution and mode structures
* Radial analysis for detailed plasma diagnostics

---

## Project Structure

```
ZPinchStabilitySim/
├── include/                 # Header files
│   ├── diagnostics.hh       # Diagnostics and data analysis
│   ├── equilibrium_solver.hh# Equilibrium solver
│   ├── linear_stability.hh  # Linear stability analyzer
│   ├── nonlinear_evolution.hh # Nonlinear evolution solver
│   └── zpinch_params.hh     # Plasma and simulation parameters
├── src/                     # Source files
│   ├── main.cpp             # Main program with mode selection
│   ├── diagnostics.cpp      # Diagnostics implementation
│   ├── equilibrium_solver.cpp # Equilibrium solver implementation
│   ├── linear_stability.cpp # Linear stability implementation
│   └── nonlinear_evolution.cpp # Nonlinear evolution implementation
├── pyscripts/               # Python visualization scripts
│   ├── animation_2d.py      # 2D plasma evolution animation
│   ├── radial_animation.py  # Radial profile evolution animation
│   ├── equilibrium_2d.py    # Equilibrium visualization
│   ├── instability_growth.py # Instability growth analysis
│   ├── mode_structure.py    # Mode structure plotting
│   └── stability_diagram.py # Linear stability diagrams
├── data/                    # Output CSV and snapshot files
├── CMakeLists.txt           # Build configuration
├── LICENSE                  # MIT License
└── README.md                # Project documentation
```

---

## Simulation Modes

* `equilibrium`  → Run equilibrium analysis only
* `linear`       → Linear stability analysis (kink/sausage modes)
* `nonlinear`    → Nonlinear evolution
* `full`         → Full pipeline (equilibrium + linear + nonlinear)
* `help`         → Show usage instructions

**Mode Options:**

* `--no-kink` → Disable kink mode
* `--sausage` → Enable sausage mode

**Default behavior:**

* Kink: ENABLED
* Sausage: DISABLED

---

## Build Instructions

### Prerequisites

* C++17 compatible compiler
* CMake 3.12+
* Python 3.6+ with **NumPy, Matplotlib**

### Build

```bash
mkdir build
cd build
cmake ..
cmake --build . --config Release
```

Executable: `zpinch_sim` (Linux/macOS) or `zpinch_sim.exe` (Windows)

---

## Usage Examples

### Equilibrium Analysis

```bash
./zpinch_sim equilibrium
```

### Linear Stability Analysis

```bash
./zpinch_sim linear --sausage
```

### Nonlinear Evolution

```bash
./zpinch_sim nonlinear --no-kink --sausage
```

### Full Analysis Pipeline

```bash
./zpinch_sim full --sausage
```

---

## Output Files

* `data/equilibrium_profiles.csv`
* `data/stability_diagram.csv` (linear analysis)
* `data/time_history.csv`
* `data/mode_structures.csv` (linear analysis)
* Snapshots: `data/snap/output_state_*.csv`
* Reports: `data/full_analysis_report.txt`

### Generated Visualization Files

* `zpinch_evolution.gif` (2D animation)
* `radial_evolution.gif` (Radial profiles animation)
* `radial_profile_comparison.png`
* `radial_growth_analysis.png`
* `pressure_magnetic_profiles.png`
* `time_evolution_plots.png`
* `zpinch_static_comparison.png`

---

## Python Visualization

* `animation_2d.py` → 2D plasma animation
* `radial_animation.py` → Radial profile evolution
* `equilibrium_2d.py` → Equilibrium plots
* `instability_growth.py` → Growth rates
* `mode_structure.py` → Linear eigenmodes
* `stability_diagram.py` → Stability diagram

**Quick Start:**

```bash
1. ./zpinch_sim nonlinear
2. python radial_analysis.py
3. View generated plots and animations
```

---

## Customization

Modify parameters in `zpinch_params.hh` or in `main.cpp`:

```cpp
params.Nr = 64; // radial points
params.Nz = 32; // axial points
params.dt = 1e-9; // time step
params.t_max = 1e-5; // total simulation time
params.n0 = 1e18; // plasma density
params.T0 = 20.0; // temperature
params.B0 = 0.05; // axial field
params.I0 = 3e5; // current
```

---

## Physics Validation

* Equilibrium force-balance verification
* Linear growth rates vs theory
* Nonlinear mode evolution
* Radial and axial profile checks
* Safety factor calculations

---

## License

MIT License © 2025 Juan Pablo Solís Ruiz

---

## Contact

* Email: [jp.sruiz18.tec@gmail.com](mailto:jp.sruiz18.tec@gmail.com)
* GitHub: [h4xter1612](https://github.com/h4xter1612)

