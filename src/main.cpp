#include "zpinch_params.hh"
#include "equilibrium_solver.hh"
#include "linear_stability.hh"
#include "nonlinear_evolution.hh"
#include "diagnostics.hh"
#include <iostream>
#include <memory>
#include <string>
#include <direct.h> // Para crear directorios en Windows

// =============================================================================
// FUNCTION DECLARATIONS
// =============================================================================
void runEquilibriumAnalysis();
void runLinearStabilityAnalysis();
void runNonlinearEvolution();
void runFullAnalysis();
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
    if (argc > 1) {
        mode = argv[1];
    }
    
    // Run selected mode
    if (mode == "equilibrium") {
        runEquilibriumAnalysis();
    } else if (mode == "linear") {
        runLinearStabilityAnalysis();
    } else if (mode == "nonlinear") {
        runNonlinearEvolution();
    } else if (mode == "full") {
        runFullAnalysis();
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
// EQUILIBRIUM ANALYSIS MODE - CORREGIDO
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
// LINEAR STABILITY ANALYSIS MODE - CORREGIDO
// =============================================================================
void runLinearStabilityAnalysis() {
    std::cout << "=== RUNNING LINEAR STABILITY ANALYSIS ===" << std::endl;
    
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
    
    // Run stability analysis
    diagnostics->analyzeLinearStability();
    diagnostics->generateStabilityDiagram("data/stability_diagram.csv");
    
    std::cout << "Linear stability analysis completed successfully!" << std::endl;
}

// =============================================================================
// NONLINEAR EVOLUTION MODE - CORREGIDO
// =============================================================================
void runNonlinearEvolution() {
    std::cout << "=== RUNNING NONLINEAR EVOLUTION ===" << std::endl;
    
    // Create parameters
    ZPinchParameters params;
    params.dt = 1e-10;  // Time step más pequeño para estabilidad
    params.t_max = 1e-6; // Tiempo más corto para pruebas
    
    // Create equilibrium
    auto equilibrium = std::make_shared<EquilibriumSolver>(params);
    equilibrium->solveForceBalance();
    
    // Create stability analyzer WITH equilibrium
    auto stability = std::make_shared<LinearStabilityAnalyzer>(params, equilibrium);
    
    // Create and run nonlinear evolution
    auto evolution = std::make_shared<NonlinearEvolution>(params, equilibrium);
    
    // Add perturbations más pequeñas para estabilidad
    evolution->addKinkPerturbation(1e-4, 2.0 * M_PI / params.L, 1); // m=1 kink (amplitud más pequeña)
    evolution->addRandomPerturbation(1e-5); // Random noise más pequeño
    
    // Create diagnostics
    auto diagnostics = std::make_unique<Diagnostics>(params, equilibrium, stability, evolution);
    
    std::cout << "Starting nonlinear evolution..." << std::endl;
    
    // Run evolution with monitoring
    evolution->evolve(params.t_max, params.t_max / 50.0); // Más outputs para debugging
    
    // Analyze results
    diagnostics->analyzeNonlinearGrowth();
    diagnostics->exportEvolutionData("data/evolution_data.csv");
    diagnostics->generateReport("data/nonlinear_report.txt");
    
    std::cout << "Nonlinear evolution completed successfully!" << std::endl;
}

// =============================================================================
// FULL ANALYSIS MODE - CORREGIDO
// =============================================================================
void runFullAnalysis() {
    std::cout << "=== RUNNING FULL ANALYSIS PIPELINE ===" << std::endl;
    
    // Create parameters
    ZPinchParameters params;
    
    // Step 1: Equilibrium
    std::cout << "Step 1/3: Equilibrium calculation..." << std::endl;
    auto equilibrium = std::make_shared<EquilibriumSolver>(params);
    equilibrium->solveForceBalance();
    
    // Step 2: Linear stability
    std::cout << "Step 2/3: Linear stability analysis..." << std::endl;
    auto stability = std::make_shared<LinearStabilityAnalyzer>(params, equilibrium);
    
    // Analyze linear stability to find most unstable mode
    auto kink_results = stability->kinkStabilityScan();
    double fastest_kink_k = 2.0 * M_PI / params.L; // Default to longest wavelength
    double max_growth = 0.0;
    
    for (const auto& result : kink_results) {
        if (result.growth_rate > max_growth) {
            max_growth = result.growth_rate;
            fastest_kink_k = result.k;
        }
    }
    
    // Step 3: Nonlinear evolution
    std::cout << "Step 3/3: Nonlinear evolution..." << std::endl;
    auto evolution = std::make_shared<NonlinearEvolution>(params, equilibrium);
    
    // Add perturbations based on linear analysis
    if (max_growth > 0.0) {
        evolution->addKinkPerturbation(1e-4, fastest_kink_k, 1);
        std::cout << "Added kink perturbation with k = " << fastest_kink_k << std::endl;
    } else {
        evolution->addKinkPerturbation(1e-4, 2.0 * M_PI / params.L, 1);
        std::cout << "Added default kink perturbation" << std::endl;
    }
    
    evolution->addRandomPerturbation(1e-5);
    
    // Create diagnostics
    auto diagnostics = std::make_unique<Diagnostics>(params, equilibrium, stability, evolution);
    
    // Run comprehensive analysis
    diagnostics->analyzeEquilibrium();
    diagnostics->analyzeLinearStability();
    
    // Run short nonlinear evolution for verification
    std::cout << "Running short nonlinear evolution for verification..." << std::endl;
    evolution->evolve(1e-6, 1e-6 / 20.0); // Even shorter for stability
    
    diagnostics->analyzeNonlinearGrowth();
    diagnostics->generateReport("data/full_analysis_report.txt");
    
    std::cout << "Full analysis pipeline completed successfully!" << std::endl;
    std::cout << "Check the generated reports in data/ directory for detailed results." << std::endl;
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
    std::cout << "Examples:" << std::endl;
    std::cout << "  ./zpinch_sim equilibrium" << std::endl;
    std::cout << "  ./zpinch_sim linear" << std::endl;
    std::cout << "  ./zpinch_sim full" << std::endl;
    std::cout << std::endl;
    std::cout << "Output files in data/:" << std::endl;
    std::cout << "  - equilibrium_profiles.csv" << std::endl;
    std::cout << "  - stability_diagram.csv" << std::endl;
    std::cout << "  - evolution_data.csv" << std::endl;
    std::cout << "  - Various diagnostic reports" << std::endl;
    std::cout << std::endl;
    std::cout << "Snapshots in data/snap/:" << std::endl;
    std::cout << "  - output_state_*.csv" << std::endl;
    std::cout << std::endl;
}
