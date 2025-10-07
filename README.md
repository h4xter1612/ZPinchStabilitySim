Z-Pinch MHD Stability Simulator

[C++17] [Python] [MHD Plasma] [CMake]

------------------------------------------------------------------------

Description

This project implements a comprehensive Z-Pinch stability analysis tool
for studying magnetohydrodynamic (MHD) instabilities in cylindrical
plasma configurations.
The code solves the full nonlinear MHD equations with realistic boundary
conditions and provides both linear stability analysis and nonlinear
evolution capabilities.

Main features

-   Equilibrium solver for force-balanced Z-pinch configurations
-   Linear stability analysis for kink (m=1) and sausage (m=0) modes
-   Nonlinear MHD evolution with artificial viscosity for numerical
    stability
-   Multiple simulation modes: equilibrium, linear, nonlinear, and full
    analysis pipeline
-   Advanced diagnostics including growth rates, energy evolution, and
    instability detection
-   Comprehensive Python visualization for 2D animations and radial
    cross-section analysis
-   Automatic stop conditions for detecting explosive instability growth
-   Progress tracking with real-time growth factor display

------------------------------------------------------------------------

Physical Model

The simulator solves the ideal MHD equations in cylindrical coordinates:

Continuity Equation:
$$\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \vec{v}) = 0$$

Momentum Equation:
$$\rho \left( \frac{\partial \vec{v}}{\partial t} + \vec{v} \cdot \nabla \vec{v} \right) = -\nabla p + \vec{J} \times \vec{B}$$

Energy Equation:
$$\frac{\partial p}{\partial t} + \vec{v} \cdot \nabla p + \gamma p \nabla \cdot \vec{v} = 0$$

Magnetic Induction:
$$\frac{\partial \vec{B}}{\partial t} = \nabla \times (\vec{v} \times \vec{B})$$

Where:
- ρ: plasma density
- v⃗: velocity field
- p: plasma pressure
- B⃗: magnetic field
- J⃗: current density (∇ × B⃗/μ₀)
- γ: adiabatic index (5/3)

------------------------------------------------------------------------

Project Structure

    ZPinchStabilitySim/
    ├── include/
    │   ├── zpinch_params.hh
    │   ├── equilibrium_solver.hh
    │   ├── linear_stability.hh
    │   ├── nonlinear_evolution.hh
    │   └── diagnostics.hh
    ├── src/
    │   ├── main.cpp
    │   ├── equilibrium_solver.cpp
    │   ├── linear_stability.cpp
    │   ├── nonlinear_evolution.cpp
    │   └── diagnostics.cpp
    ├── pyscripts/
    │   ├── animation_2d.py
    │   ├── radial_analysis.py
    │   ├── instability_growth.py
    │   ├── equilibrium_2d.py
    │   ├── mode_structure.py
    │   └── stability_diagram.py
    ├── CMakeLists.txt
    └── LICENSE

------------------------------------------------------------------------

Features

Simulation Capabilities

-   Equilibrium calculation with self-consistent force balance
-   Linear stability analysis for kink and sausage modes
-   Nonlinear MHD evolution with Euler time integration
-   Multiple instability modes: kink (m=1), sausage (m=0), and mixed
    modes
-   Automatic growth detection with configurable thresholds
-   CFL-based time stepping for numerical stability
-   Artificial viscosity for damping numerical oscillations

Analysis Tools

-   2D field visualization (v_(r), B_(r), B_(θ), pressure)
-   Radial cross-section analysis at fixed axial positions
-   Instability growth tracking with exponential fit analysis
-   Energy evolution (kinetic, magnetic, internal, total)
-   Mode structure visualization from linear analysis
-   Stability diagrams for parameter space exploration

------------------------------------------------------------------------

Build Instructions

Prerequisites

-   C++17 compiler (GCC 7+, Clang 5+, MSVC 2019+)
-   CMake 3.12+
-   Python 3.6+ with NumPy, Matplotlib, Pandas

Building the Project

    mkdir build && cd build
    cmake ..
    cmake --build . --config Release

------------------------------------------------------------------------

Usage

Available Modes

    ./zpinch_sim equilibrium
    ./zpinch_sim linear
    ./zpinch_sim nonlinear
    ./zpinch_sim full

Command Line Options

-   --no-kink → Disable kink mode
-   --sausage → Enable sausage mode
-   --help → Show help message

------------------------------------------------------------------------

Simulation Output

Data Files (in data/)

-   equilibrium_profiles.csv
-   stability_diagram.csv
-   time_history.csv
-   mode_structures.csv
-   equilibrium_diagnostics.csv
-   full_analysis_report.txt

Snapshots (in data/snap/)

-   output_state_*.csv — 2D field snapshots for visualization

------------------------------------------------------------------------

Visualization

    python pyscripts/animation_2d.py
    python pyscripts/radial_analysis.py
    python pyscripts/instability_growth.py
    python pyscripts/equilibrium_2d.py
    python pyscripts/mode_structure.py
    python pyscripts/stability_diagram.py

------------------------------------------------------------------------

Example Results

-   zpinch_evolution.gif — 2D plasma evolution animation
-   radial_evolution.gif — Radial profile evolution
-   time_evolution_plots.png — Diagnostics
-   zpinch_static_comparison.png — 2D snapshots

------------------------------------------------------------------------

License

MIT License © 2025 Juan Pablo Solís Ruiz

------------------------------------------------------------------------

References

-   Freidberg, J. P. Ideal Magnetohydrodynamics (1987)
-   Bateman, G. MHD Instabilities (1978)
-   Shafranov, V. D. Reviews of Plasma Physics (1966)
-   Ryu & Jones. Numerical Magnetohydrodynamics (1995)

------------------------------------------------------------------------

Contact

Author: Juan Pablo Solís Ruiz
Email: jp.sruiz18.tec@gmail.com
GitHub: h4xter1612

------------------------------------------------------------------------

Future Work

-   Higher-order time integration (RK4, Adams-Bashforth)
-   Adaptive mesh refinement
-   Hall MHD and resistive effects
-   OpenMP/MPI parallelization
-   GPU acceleration (CUDA)
-   Experimental validation

