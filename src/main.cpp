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
void printHelp();
void createDataDirectory();

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
    ZPinchParameters params;
    
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
    ZPinchParameters params;
    
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
        std::cout << "No stability modes enabled. Skipping analysis." << std::endl;
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
    
    // Realistic parameters
    ZPinchParameters params;
    params.Nr = 64;  
    params.Nz = 32;
    params.dt = 1e-9;
    params.t_max = 1e-5;
    
    auto equilibrium = std::make_shared<EquilibriumSolver>(params);
    equilibrium->solveForceBalance();
    
    auto stability = std::make_shared<LinearStabilityAnalyzer>(params, equilibrium);
    auto evolution = std::make_shared<NonlinearEvolution>(params, equilibrium);
    
    // Realistic perturbations
    double vA = params.alfven_speed();
    std::cout << "Alfven speed: " << vA << " m/s" << std::endl;
    
    double kz_realistic = 2.0 * M_PI / params.L;
    
    // Add perturbations based on enabled modes
    if (enable_kink) {
        evolution->addKinkPerturbation(0.01, kz_realistic, 1);
        std::cout << "Added kink perturbation" << std::endl;
    }
    
    if (enable_sausage) {
        evolution->addRandomPerturbation(0.001);
        std::cout << "Added sausage perturbation" << std::endl;
    }
    
    // Add random perturbation if no specific modes enabled
    if (!enable_kink && !enable_sausage) {
        evolution->addRandomPerturbation(0.001);
        std::cout << "Added random perturbation (no specific modes enabled)" << std::endl;
    }
    
    auto diagnostics = std::make_unique<Diagnostics>(params, equilibrium, stability, evolution);
    
    std::cout << "Starting evolution..." << std::endl;
    evolution->evolve(params.t_max, params.t_max / 50.0);
    
    std::cout << "Nonlinear evolution completed!" << std::endl;
}

// =============================================================================
// FULL ANALYSIS MODE
// =============================================================================
void runFullAnalysis(bool enable_kink, bool enable_sausage) {
    std::cout << "=== RUNNING FULL ANALYSIS PIPELINE ===" << std::endl;
    std::cout << "Kink modes: " << (enable_kink ? "ENABLED" : "DISABLED") << std::endl;
    std::cout << "Sausage modes: " << (enable_sausage ? "ENABLED" : "DISABLED") << std::endl;
    
    // Realistic parameters
    ZPinchParameters params;
    params.Nr = 64;
    params.Nz = 32; 
    params.dt = 1e-9;
    params.t_max = 5e-5;
    
    // Check stability criteria
    checkStabilityCriteria(params);
    
    // Step 1: Equilibrium
    std::cout << "Step 1/3: Equilibrium calculation..." << std::endl;
    auto equilibrium = std::make_shared<EquilibriumSolver>(params);
    equilibrium->solveForceBalance();
    equilibrium->printEquilibriumSummary();
    equilibrium->saveProfiles("data/equilibrium_profiles.csv"); // Para equilibrium_2d.py
    
    // Step 2: Linear stability
    std::cout << "Step 2/3: Linear stability analysis..." << std::endl;
    auto stability = std::make_shared<LinearStabilityAnalyzer>(params, equilibrium);

    // Find most unstable mode (only if kink is enabled)
    double fastest_kink_k = 10.0;
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
        } else {
            fastest_kink_k = 2.0 * M_PI / params.L;
            std::cout << "Using default k: " << fastest_kink_k << " (no results from scan)" << std::endl;
        }
        
        std::cout << "Most unstable kink mode: k=" << fastest_kink_k 
                  << ", growth_rate=" << max_growth << " s^-1" << std::endl;
    }
    
    // Step 3: Nonlinear evolution
    std::cout << "Step 3/3: Nonlinear evolution..." << std::endl;
    auto evolution = std::make_shared<NonlinearEvolution>(params, equilibrium);

    // Realistic perturbations (dimensionless amplitudes)
    double kink_amplitude = 0.005;  // 0.5% 
    double random_amplitude = 0.001; // 0.1%
    
    std::cout << "Adding perturbations:" << std::endl;
    
    // Add perturbations based on enabled modes
    if (enable_kink) {
        evolution->addKinkPerturbation(kink_amplitude, fastest_kink_k, 1);
        std::cout << "  - Kink perturbation: amplitude=" << kink_amplitude << std::endl;
    }
    
    if (enable_sausage) {
        evolution->addRandomPerturbation(random_amplitude);
        std::cout << "  - Sausage perturbation: amplitude=" << random_amplitude << std::endl;
    }
    
    // Add random perturbation if no specific modes enabled
    if (!enable_kink && !enable_sausage) {
        evolution->addRandomPerturbation(random_amplitude);
        std::cout << "  - Random perturbation: amplitude=" << random_amplitude << " (no specific modes enabled)" << std::endl;
    }
    
    // Create diagnostics
    auto diagnostics = std::make_unique<Diagnostics>(params, equilibrium, stability, evolution);
    
    // Run analysis
    std::cout << "=== EQUILIBRIUM ANALYSIS ===" << std::endl;
    diagnostics->analyzeEquilibrium();
    
    if (enable_kink || enable_sausage) {
        std::cout << "=== LINEAR STABILITY ANALYSIS ===" << std::endl;
        diagnostics->analyzeLinearStability();
        diagnostics->generateStabilityDiagram("data/stability_diagram.csv"); // Para stability_diagram.py
        
        stability->exportEigenfunctions("data/mode_structures.csv");
    }
    
    // Run evolution
    std::cout << "Running nonlinear evolution..." << std::endl;
    std::cout << "Simulation time: " << params.t_max << " s" << std::endl;
    
    if (enable_kink && max_growth > 0) {
        std::cout << "Expected growth time: " << 1.0/max_growth << " s" << std::endl;
    }
    
    evolution->evolve(params.t_max, params.t_max / 100.0);
    
    std::cout << "=== NONLINEAR ANALYSIS ===" << std::endl;
    diagnostics->analyzeNonlinearGrowth();
    // diagnostics->exportEvolutionData("data/evolution_data.csv"); // Para instability_growth.py
    diagnostics->generateReport("data/full_analysis_report.txt");
    
    std::cout << "Full analysis pipeline completed successfully!" << std::endl;
    std::cout << "Generated data for ALL Python scripts:" << std::endl;
    std::cout << "  - data/equilibrium_profiles.csv (for equilibrium_2d.py)" << std::endl;
    std::cout << "  - data/stability_diagram.csv (for stability_diagram.py)" << std::endl;
    std::cout << "  - data/time_history.csv (for instability_growth.py)" << std::endl;
    std::cout << "  - data/snap/output_state_*.csv (for animation_2d.py)" << std::endl;
    std::cout << "  - data/mode_structures.csv (for mode_structure.py)" << std::endl;
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
    std::cout << "  ./zpinch_sim full --sausage" << std::endl;
    std::cout << std::endl;
    std::cout << "Default behavior:" << std::endl;
    std::cout << "  - Kink modes: ENABLED" << std::endl;
    std::cout << "  - Sausage modes: DISABLED" << std::endl;
    std::cout << std::endl;
    std::cout << "Output files in data/:" << std::endl;
    std::cout << "  - equilibrium_profiles.csv" << std::endl;
    std::cout << "  - stability_diagram.csv" << std::endl;
    std::cout << "  - time_history.csv" << std::endl;
    std::cout << "  - Various diagnostic reports" << std::endl;
    std::cout << std::endl;
    std::cout << "Snapshots in data/snap/:" << std::endl;
    std::cout << "  - output_state_*.csv" << std::endl;
    std::cout << std::endl;
}
