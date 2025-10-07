#include "zpinch_params.hh"
#include "equilibrium_solver.hh"
#include "linear_stability.hh"
#include "nonlinear_evolution.hh"
#include "diagnostics.hh"
#include <iostream>
#include <memory>
#include <string>
#include <direct.h> // For directory creation on Windows

// =============================================================================
// FUNCTION DECLARATIONS
// =============================================================================
void runEquilibriumAnalysis();
void runLinearStabilityAnalysis(bool enable_kink, bool enable_sausage);
void runNonlinearEvolution(bool enable_kink, bool enable_sausage);
void runFullAnalysis(bool enable_kink, bool enable_sausage);
void runStableMode();
void printHelp();
void createDataDirectory();

// =============================================================================
// PARAMETER SETTINGS FOR DIFFERENT MODES
// =============================================================================
ZPinchParameters getParametersForMode(bool enable_kink, bool enable_sausage, bool is_full_analysis = false) {
    ZPinchParameters params;
    
    // Common simulation parameters
    params.Nr = 64;
    params.Nz = 32;
    params.dt = 1e-9;
    params.t_max = is_full_analysis ? 5e-5 : 1e-5; // NEW NEW
    
    // Different time settings for different modes
    if (is_full_analysis) {
        params.t_max = 5e-5;  // Longer for full analysis
    } else {
        params.t_max = 1e-5;  // Shorter for individual modes
    }
    
    // Mode-specific physics parameters
    if (enable_sausage && !enable_kink) {
        // Sausage-only mode: parameters that favor sausage but are stable
        params.n0 = 1e18;        // Lower density for stability
        params.T0 = 30.0;        // Moderate temperature
        params.B0 = 0.05;        // Higher axial field for stability
        params.I0 = 2e4;         // Lower current
    } else if (!enable_kink && !enable_sausage) {
        // Stable mode: high stability parameters
        params.n0 = 1e16;        // Much lower density 1E17
        params.T0 = 1.0;         // Lower temperature  5
        params.B0 = 0.5;        // Higher axial field (stabilizing) 0.15
        params.I0 = 500;         // Lower current (stabilizing) 1e4
        params.a = 0.1;          // Radio pequeño NEW NEW
        params.L = 1.0;          // Longitud moderada NEW NEW
    } else if (enable_kink && !enable_sausage) {
        // Kink-only mode: classic unstable parameters
        params.n0 = 5e18;        // Higher density
        params.T0 = 20.0;        // Moderate temperature
        params.B0 = 0.05;        // Moderate axial field
        params.I0 = 3e5;         // Higher current (kink-unstable)
    } else {
        // Mixed mode: balanced parameters
        params.n0 = 2e18;        // Moderate density
        params.T0 = 25.0;        // Moderate temperature
        params.B0 = 0.04;        // Moderate axial field
        params.I0 = 1e5;         // Moderate current
    }
    
    return params;
}

// =============================================================================
// PERTURBATION SETTINGS - CORREGIDO
// =============================================================================
void applyPerturbations(std::shared_ptr<NonlinearEvolution> evolution, 
                       bool enable_kink, bool enable_sausage, 
                       double kz_realistic, double fastest_kink_k = 0.0) {
    // Use provided kink k or default to realistic k
    double kink_k = (fastest_kink_k > 0) ? fastest_kink_k : kz_realistic;
    
    std::cout << "Adding perturbations:" << std::endl;
    
    if (enable_kink) {
        double kink_amplitude = 0.000001;  // Small amplitude for stability 0.001
        evolution->addKinkPerturbation(kink_amplitude, kink_k, 1);
        std::cout << "  - Kink perturbation: amplitude=" << kink_amplitude 
                  << ", k=" << kink_k << std::endl;
    }
    
    if (enable_sausage) {
        double sausage_amplitude = 0.0000005; // Very small amplitude 0.0005
        evolution->addSausagePerturbation(sausage_amplitude, kz_realistic);
        std::cout << "  - Sausage perturbation: amplitude=" << sausage_amplitude << std::endl;
    }
    
    if (!enable_kink && !enable_sausage) {
        // STABLE MODE: Add VERY small noise or none at all
        double random_amplitude = 0.000001; // Extremely small - almost negligible
        evolution->addRandomPerturbation(random_amplitude);
        std::cout << "  - Minimal random noise: amplitude=" << random_amplitude 
                  << " (stable mode)" << std::endl;
    }
}

// =============================================================================
// MAIN FUNCTION
// =============================================================================
int main(int argc, char* argv[]) {
    std::cout << "==========================================" << std::endl;
    std::cout << "    Z-Pinch Stability Simulator" << std::endl;
    std::cout << "    MHD Stability Analysis Tool" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << std::endl;
    
    // Create data directory
    createDataDirectory();
    
    // Parse command line arguments
    std::string mode = "help";
    bool enable_kink = true;      // Default: kink enabled
    bool enable_sausage = false;  // Default: sausage disabled
    
    if (argc > 1) {
        mode = argv[1];
    }
    
    // Parse additional flags
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--no-kink") {
            enable_kink = false;
        } else if (arg == "--sausage") {
            enable_sausage = true;
        } else if (arg == "--help") {
            printHelp();
            return 0;
        }
    }
    
    // Special case: if no modes enabled, use stable mode for nonlinear
    if (mode == "nonlinear" && !enable_kink && !enable_sausage) {
        runStableMode();
        return 0;
    }
    
    // Run selected mode
    if (mode == "equilibrium") {
        runEquilibriumAnalysis();
    } else if (mode == "linear") {
        runLinearStabilityAnalysis(enable_kink, enable_sausage);
    } else if (mode == "nonlinear") {
        runNonlinearEvolution(enable_kink, enable_sausage);
    } else if (mode == "full") {
        runFullAnalysis(enable_kink, enable_sausage);
    } else if (mode == "help") {
        printHelp();
    } else {
        std::cout << "Unknown mode: " << mode << std::endl;
        std::cout << "Use 'help' for available options." << std::endl;
        return 1;
    }
    
    return 0;
}

// =============================================================================
// DIRECTORY CREATION
// =============================================================================
void createDataDirectory() {
    // Create data directory if it doesn't exist
    _mkdir("data");
    _mkdir("data/snap");
}

// =============================================================================
// EQUILIBRIUM ANALYSIS MODE
// =============================================================================
void runEquilibriumAnalysis() {
    std::cout << "=== RUNNING EQUILIBRIUM ANALYSIS ===" << std::endl;
    
    // Create parameters
    ZPinchParameters params = getParametersForMode(true, false);
    
    // Create and run equilibrium solver
    auto equilibrium = std::make_shared<EquilibriumSolver>(params);
    equilibrium->solveForceBalance();
    equilibrium->printEquilibriumSummary();
    equilibrium->saveProfiles("data/equilibrium_profiles.csv");
    
    // Create stability analyzer WITH equilibrium
    auto stability = std::make_shared<LinearStabilityAnalyzer>(params, equilibrium);
    
    // Create evolution and diagnostics
    auto evolution = std::make_shared<NonlinearEvolution>(params, equilibrium);
    auto diagnostics = std::make_unique<Diagnostics>(params, equilibrium, stability, evolution);
    
    // Run equilibrium diagnostics
    diagnostics->analyzeEquilibrium();
    diagnostics->exportEquilibriumData("data/equilibrium_diagnostics.csv");
    diagnostics->generateReport("data/equilibrium_report.txt");
    
    std::cout << "Equilibrium analysis completed successfully!" << std::endl;
}

// =============================================================================
// LINEAR STABILITY ANALYSIS MODE
// =============================================================================
void runLinearStabilityAnalysis(bool enable_kink, bool enable_sausage) {
    std::cout << "=== RUNNING LINEAR STABILITY ANALYSIS ===" << std::endl;
    std::cout << "Kink modes: " << (enable_kink ? "ENABLED" : "DISABLED") << std::endl;
    std::cout << "Sausage modes: " << (enable_sausage ? "ENABLED" : "DISABLED") << std::endl;
    
    // Create parameters
    ZPinchParameters params = getParametersForMode(enable_kink, enable_sausage);
    
    // Create equilibrium
    auto equilibrium = std::make_shared<EquilibriumSolver>(params);
    equilibrium->solveForceBalance();
    
    // Create stability analyzer WITH equilibrium
    auto stability = std::make_shared<LinearStabilityAnalyzer>(params, equilibrium);
    
    // Create diagnostics
    auto evolution = std::make_shared<NonlinearEvolution>(params, equilibrium);
    auto diagnostics = std::make_unique<Diagnostics>(params, equilibrium, stability, evolution);
    
    // Run stability analysis based on enabled modes
    if (enable_kink || enable_sausage) {
        diagnostics->analyzeLinearStability();
        diagnostics->generateStabilityDiagram("data/stability_diagram.csv");
        
        stability->exportEigenfunctions("data/mode_structures.csv");
    } else {
        std::cout << "No stability modes enabled. Skipping linear analysis." << std::endl;
    }
    
    std::cout << "Linear stability analysis completed!" << std::endl;
}

// =============================================================================
// NONLINEAR EVOLUTION MODE
// =============================================================================
void runNonlinearEvolution(bool enable_kink, bool enable_sausage) {
    std::cout << "=== RUNNING NONLINEAR EVOLUTION ===" << std::endl;
    std::cout << "Kink modes: " << (enable_kink ? "ENABLED" : "DISABLED") << std::endl;
    std::cout << "Sausage modes: " << (enable_sausage ? "ENABLED" : "DISABLED") << std::endl;
    
    // Get consistent parameters
    ZPinchParameters params = getParametersForMode(enable_kink, enable_sausage);
    
    auto equilibrium = std::make_shared<EquilibriumSolver>(params);
    equilibrium->solveForceBalance();
    
    auto stability = std::make_shared<LinearStabilityAnalyzer>(params, equilibrium);
    auto evolution = std::make_shared<NonlinearEvolution>(params, equilibrium);
    
    // Physical parameters info
    double vA = params.alfven_speed();
    std::cout << "Physical parameters:" << std::endl;
    std::cout << "  - Alfven speed: " << vA << " m/s" << std::endl;
    std::cout << "  - Plasma beta: " << params.beta() << std::endl;
    
    // Calculate safety factor for stability assessment
    double B_theta_surface = (4e-7 * M_PI * params.I0) / (2.0 * M_PI * params.a);
    double safety_factor = (params.a * params.B0) / (params.R0 * B_theta_surface);
    std::cout << "  - Safety factor: " << safety_factor << std::endl;
    
    double kz_realistic = 2.0 * M_PI / params.L;
    
    // Apply consistent perturbations
    applyPerturbations(evolution, enable_kink, enable_sausage, kz_realistic);
    
    auto diagnostics = std::make_unique<Diagnostics>(params, equilibrium, stability, evolution);
    
    std::cout << "Starting evolution..." << std::endl;
    evolution->evolve(params.t_max, params.t_max / 50.0);
    
    std::cout << "Nonlinear evolution completed!" << std::endl;
}

// =============================================================================
// STABLE MODE
// =============================================================================
void runStableMode() {
    std::cout << "=== RUNNING STABLE MODE ===" << std::endl;
    std::cout << "No instabilities enabled - testing stable configuration" << std::endl;
    
    // Get ultra-stable parameters
    ZPinchParameters params = getParametersForMode(false, false);
    
    // Additional stability enhancements
    params.n0 = 5e16;        // Very low density
    params.T0 = 2.0;         // Very low temperature  
    params.B0 = 0.3;         // Very high axial field
    params.I0 = 1e3;         // Very low current
    
    auto equilibrium = std::make_shared<EquilibriumSolver>(params);
    equilibrium->solveForceBalance();
    
    auto stability = std::make_shared<LinearStabilityAnalyzer>(params, equilibrium);
    auto evolution = std::make_shared<NonlinearEvolution>(params, equilibrium);

    // ✅ DISABLE STOP CONDITIONS FOR STABLE MODE
    evolution->disableStopConditions();
    
    std::cout << "Physical parameters:" << std::endl;
    std::cout << "  - Alfven speed: " << params.alfven_speed() << " m/s" << std::endl;
    std::cout << "  - Plasma beta: " << params.beta() << std::endl;
    
    // Calculate safety factor
    double B_theta_surface = (4e-7 * M_PI * params.I0) / (2.0 * M_PI * params.a);
    double safety_factor = (params.a * params.B0) / (params.R0 * B_theta_surface);
    std::cout << "  - Safety factor: " << safety_factor << " (should be > 2 for stability)" << std::endl;
    
    double kz_realistic = 2.0 * M_PI / params.L;
    
    // Apply minimal perturbations for stable mode
    applyPerturbations(evolution, false, false, kz_realistic);
    
    auto diagnostics = std::make_unique<Diagnostics>(params, equilibrium, stability, evolution);
    
    std::cout << "Starting evolution..." << std::endl;
    evolution->evolve(params.t_max, params.t_max / 50.0);
    
    std::cout << "Stable mode completed!" << std::endl;
}

// =============================================================================
// FULL ANALYSIS MODE
// =============================================================================
void runFullAnalysis(bool enable_kink, bool enable_sausage) {
    std::cout << "=== RUNNING FULL ANALYSIS PIPELINE ===" << std::endl;
    std::cout << "Kink modes: " << (enable_kink ? "ENABLED" : "DISABLED") << std::endl;
    std::cout << "Sausage modes: " << (enable_sausage ? "ENABLED" : "DISABLED") << std::endl;
    
    // Use consistent parameters with other modes
    ZPinchParameters params = getParametersForMode(enable_kink, enable_sausage, true);
    
    // Check stability criteria
    checkStabilityCriteria(params);
    
    // Step 1: Equilibrium
    std::cout << "Step 1/3: Equilibrium calculation..." << std::endl;
    auto equilibrium = std::make_shared<EquilibriumSolver>(params);
    equilibrium->solveForceBalance();
    equilibrium->printEquilibriumSummary();
    equilibrium->saveProfiles("data/equilibrium_profiles.csv");
    
    // Step 2: Linear stability (SOLO si algún modo está habilitado)
    std::cout << "Step 2/3: Linear stability analysis..." << std::endl;
    auto stability = std::make_shared<LinearStabilityAnalyzer>(params, equilibrium);

    // Find most unstable kink mode if enabled
    double fastest_kink_k = 2.0 * M_PI / params.L; // Default value
    double max_growth = 0.0;
    
    if (enable_kink) {
        auto kink_results = stability->kinkStabilityScan();
        
        if (!kink_results.empty()) {
            for (const auto& result : kink_results) {
                if (result.growth_rate > max_growth && result.growth_rate < 1e7) {
                    max_growth = result.growth_rate;
                    fastest_kink_k = result.k;
                }
            }
            std::cout << "Most unstable kink mode: k=" << fastest_kink_k 
                      << ", growth_rate=" << max_growth << " s^-1" << std::endl;
        } else {
            std::cout << "Using default k: " << fastest_kink_k << " (no unstable modes found)" << std::endl;
        }
    }
    
    // Step 3: Nonlinear evolution
    std::cout << "Step 3/3: Nonlinear evolution..." << std::endl;
    auto evolution = std::make_shared<NonlinearEvolution>(params, equilibrium);
    
    // Physical parameters info
    double vA = params.alfven_speed();
    std::cout << "Physical parameters:" << std::endl;
    std::cout << "  - Alfven speed: " << vA << " m/s" << std::endl;
    std::cout << "  - Plasma beta: " << params.beta() << std::endl;
    
    double kz_realistic = 2.0 * M_PI / params.L;
    
    // Apply consistent perturbations
    applyPerturbations(evolution, enable_kink, enable_sausage, kz_realistic, fastest_kink_k);
    
    // Create diagnostics
    auto diagnostics = std::make_unique<Diagnostics>(params, equilibrium, stability, evolution);
    
    // Run comprehensive analysis
    std::cout << "=== EQUILIBRIUM ANALYSIS ===" << std::endl;
    diagnostics->analyzeEquilibrium();
    
    // LINEAR ANALYSIS SOLO SI ALGÚN MODO ESTÁ HABILITADO
    if (enable_kink || enable_sausage) {
        std::cout << "=== LINEAR STABILITY ANALYSIS ===" << std::endl;
        diagnostics->analyzeLinearStability();
        diagnostics->generateStabilityDiagram("data/stability_diagram.csv");
        stability->exportEigenfunctions("data/mode_structures.csv");
    } else {
        std::cout << "=== SKIPPING LINEAR STABILITY ANALYSIS ===" << std::endl;
        std::cout << "No instability modes enabled for linear analysis." << std::endl;
        evolution->disableStopConditions();
    }
    
    // Run evolution
    std::cout << "Running nonlinear evolution..." << std::endl;
    std::cout << "Simulation time: " << params.t_max << " s" << std::endl;
    
    if (enable_kink && max_growth > 0) {
        std::cout << "Expected growth time: " << 1.0/max_growth << " s" << std::endl;
    }
    
    evolution->evolve(params.t_max, params.t_max / 100.0);
    
    // Final analysis
    std::cout << "=== NONLINEAR ANALYSIS ===" << std::endl;
    diagnostics->analyzeNonlinearGrowth();
    diagnostics->generateReport("data/full_analysis_report.txt");
    
    std::cout << "Full analysis pipeline completed successfully!" << std::endl;
    std::cout << "Generated data for Python scripts:" << std::endl;
    if (enable_kink || enable_sausage) {
        std::cout << "  - data/stability_diagram.csv (for stability_diagram.py)" << std::endl;
        std::cout << "  - data/mode_structures.csv (for mode_structure.py)" << std::endl;
    }
    std::cout << "  - data/equilibrium_profiles.csv (for equilibrium_2d.py)" << std::endl;
    std::cout << "  - data/time_history.csv (for instability_growth.py)" << std::endl;
    std::cout << "  - data/snap/output_state_*.csv (for animation_2d.py)" << std::endl;
    std::cout << "  - data/full_analysis_report.txt (comprehensive report)" << std::endl;
}
// =============================================================================
// HELP FUNCTION
// =============================================================================
void printHelp() {
    std::cout << "Z-Pinch Stability Simulator - Usage Guide" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << std::endl;
    std::cout << "Available modes:" << std::endl;
    std::cout << "  equilibrium  - Run equilibrium analysis only" << std::endl;
    std::cout << "  linear       - Run linear stability analysis" << std::endl;
    std::cout << "  nonlinear    - Run nonlinear evolution" << std::endl;
    std::cout << "  full         - Run complete analysis pipeline" << std::endl;
    std::cout << "  help         - Show this help message (default)" << std::endl;
    std::cout << std::endl;
    std::cout << "Mode options:" << std::endl;
    std::cout << "  --no-kink    - Disable kink mode analysis" << std::endl;
    std::cout << "  --sausage    - Enable sausage mode analysis" << std::endl;
    std::cout << "  --help       - Show this help message" << std::endl;
    std::cout << std::endl;
    std::cout << "Examples:" << std::endl;
    std::cout << "  ./zpinch_sim equilibrium" << std::endl;
    std::cout << "  ./zpinch_sim linear --sausage" << std::endl;
    std::cout << "  ./zpinch_sim nonlinear --no-kink --sausage" << std::endl;
    std::cout << "  ./zpinch_sim nonlinear --no-kink          (stable mode)" << std::endl;
    std::cout << "  ./zpinch_sim full --sausage" << std::endl;
    std::cout << std::endl;
    std::cout << "Default behavior:" << std::endl;
    std::cout << "  - Kink modes: ENABLED" << std::endl;
    std::cout << "  - Sausage modes: DISABLED" << std::endl;
    std::cout << std::endl;
    std::cout << "Output files in data/:" << std::endl;
    std::cout << "  - equilibrium_profiles.csv" << std::endl;
    std::cout << "  - stability_diagram.csv (only if modes enabled)" << std::endl;
    std::cout << "  - time_history.csv" << std::endl;
    std::cout << "  - mode_structures.csv (only if modes enabled)" << std::endl;
    std::cout << "  - Various diagnostic reports" << std::endl;
    std::cout << std::endl;
    std::cout << "Snapshots in data/snap/:" << std::endl;
    std::cout << "  - output_state_*.csv (used for Python visualization)" << std::endl;
    std::cout << std::endl;
    std::cout << "Python Analysis Scripts:" << std::endl;
    std::cout << "  animation_2d.py        - 2D animation of plasma evolution" << std::endl;
    std::cout << "  radial_analysis.py     - Radial cross-section analysis (NEW!)" << std::endl;
    std::cout << "  instability_growth.py  - Growth rate analysis" << std::endl;
    std::cout << "  equilibrium_2d.py      - Equilibrium visualization" << std::endl;
    std::cout << std::endl;
    std::cout << "Radial Analysis Features (radial_analysis.py):" << std::endl;
    std::cout << "  - Radial profile evolution animation" << std::endl;
    std::cout << "  - Cross-section analysis at fixed axial position" << std::endl;
    std::cout << "  - Instability growth at different radii" << std::endl;
    std::cout << "  - Pressure and magnetic field profile evolution" << std::endl;
    std::cout << "  - Static comparison of key time points" << std::endl;
    std::cout << std::endl;
    std::cout << "Generated Visualization Files:" << std::endl;
    std::cout << "  - zpinch_evolution.gif (2D animation)" << std::endl;
    std::cout << "  - radial_evolution.gif (Radial profiles animation)" << std::endl;
    std::cout << "  - radial_profile_comparison.png" << std::endl;
    std::cout << "  - radial_growth_analysis.png" << std::endl;
    std::cout << "  - pressure_magnetic_profiles.png" << std::endl;
    std::cout << "  - time_evolution_plots.png" << std::endl;
    std::cout << "  - zpinch_static_comparison.png" << std::endl;
    std::cout << std::endl;
    std::cout << "Quick Start:" << std::endl;
    std::cout << "  1. Run simulation: ./zpinch_sim nonlinear" << std::endl;
    std::cout << "  2. Analyze results: python radial_analysis.py" << std::endl;
    std::cout << "  3. View animations and plots in current directory" << std::endl;
    std::cout << std::endl;
    std::cout << "Advanced Analysis:" << std::endl;
    std::cout << "  - Use 'full' mode for complete stability analysis" << std::endl;
    std::cout << "  - Combine kink and sausage modes for mixed instability studies" << std::endl;
    std::cout << "  - Use stable mode (--no-kink) for stability verification" << std::endl;
    std::cout << std::endl;
}
