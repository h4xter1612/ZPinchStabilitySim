#ifndef DIAGNOSTICS_HH
#define DIAGNOSTICS_HH

#include "zpinch_params.hh"
#include "equilibrium_solver.hh"
#include "linear_stability.hh"
#include "nonlinear_evolution.hh"
#include <vector>
#include <string>
#include <memory>

class Diagnostics {
private:
    ZPinchParameters params_;
    std::shared_ptr<EquilibriumSolver> equilibrium_;
    std::shared_ptr<LinearStabilityAnalyzer> stability_;
    std::shared_ptr<NonlinearEvolution> evolution_;
    
    // Diagnostic data storage
    std::vector<double> time_data_;
    std::vector<std::vector<double>> energy_data_;
    std::vector<std::vector<double>> mode_amplitude_data_;
    std::vector<double> growth_rate_data_;
    
public:
    // =========================================================================
    // CONSTRUCTOR AND DESTRUCTOR
    // =========================================================================
    Diagnostics(const ZPinchParameters& params,
                std::shared_ptr<EquilibriumSolver> equilibrium,
                std::shared_ptr<LinearStabilityAnalyzer> stability,
                std::shared_ptr<NonlinearEvolution> evolution);
    ~Diagnostics() = default;
    
    // =========================================================================
    // EQUILIBRIUM DIAGNOSTICS
    // =========================================================================
    
    // Comprehensive equilibrium analysis
    void analyzeEquilibrium();
    
    // Calculate global equilibrium parameters
    struct EquilibriumSummary {
        double total_current;
        double internal_inductance;
        double poloidal_beta;
        double toroidal_beta;
        double total_beta;
        double safety_factor_edge;
        double kink_stability_margin;
        double sausage_stability_margin;
    };
    
    EquilibriumSummary computeEquilibriumSummary() const;
    
    // =========================================================================
    // LINEAR STABILITY DIAGNOSTICS
    // =========================================================================
    
    // Perform complete linear stability analysis
    void analyzeLinearStability();
    
    // Generate stability diagrams
    void generateStabilityDiagram(const std::string& filename) const;
    
    // Calculate stability boundaries
    struct StabilityBoundaries {
        double kink_limit_current;
        double sausage_limit_beta;
        double kruskal_shafranov_length;
        double troyon_limit;
    };
    
    StabilityBoundaries computeStabilityBoundaries() const;
    
    // =========================================================================
    // NONLINEAR EVOLUTION DIAGNOSTICS
    // =========================================================================
    
    // Monitor nonlinear evolution
    void monitorEvolution(double output_interval);
    
    // Analyze growth rates from nonlinear simulation
    void analyzeNonlinearGrowth();
    
    // Extract mode amplitudes from nonlinear data
    std::vector<double> extractModeAmplitude(int m, double kz) const;
    
    // =========================================================================
    // ENERGY DIAGNOSTICS
    // =========================================================================
    
    // Track energy components over time
    struct EnergyComponents {
        double kinetic;
        double magnetic;
        double internal;
        double total;
    };
    
    EnergyComponents computeInstantaneousEnergy() const;
    std::vector<EnergyComponents> computeEnergyHistory() const;
    
    // Energy transfer analysis
    void analyzeEnergyTransfer();
    
    // =========================================================================
    // MODE ANALYSIS
    // =========================================================================
    
    // Fourier analysis of perturbations
    std::vector<double> computeFourierSpectrum(int variable_index) const;
    
    // Dominant mode identification
    struct DominantMode {
        int m;          // Azimuthal mode number
        double kz;      // Axial wavenumber
        double amplitude;
        double growth_rate;
    };
    
    DominantMode identifyDominantMode() const;
    
    // =========================================================================
    // CONSERVATION CHECKS
    // =========================================================================
    
    // Check conservation laws
    struct ConservationMetrics {
        double mass_conservation;
        double energy_conservation;
        double momentum_conservation;
        double magnetic_flux_conservation;
    };
    
    ConservationMetrics checkConservation() const;
    
    // =========================================================================
    // VISUALIZATION DATA EXPORT
    // =========================================================================
    
    // Export data for Python visualization
    void exportEquilibriumData(const std::string& filename) const;
    void exportStabilityData(const std::string& filename) const;
    void exportEvolutionData(const std::string& filename) const;
    void exportModeStructureData(const std::string& filename, int m, double kz) const;
    
    // =========================================================================
    // REPORT GENERATION
    // =========================================================================
    
    // Generate comprehensive diagnostic report
    void generateReport(const std::string& filename) const;
    
    // =========================================================================
    // REAL-TIME MONITORING
    // =========================================================================
    
    // Print status during simulation
    void printStatus(double current_time) const;
    
    // Check for simulation problems
    bool checkSimulationHealth() const;
    
private:
    // =========================================================================
    // INTERNAL METHODS
    // =========================================================================
    
    // Helper functions for specific calculations
    double computeMagneticShear() const;
    double computeMercierCriterion() const;
    double computeBallooningStability() const;
    
    // Data analysis utilities
    double exponentialFit(const std::vector<double>& x, const std::vector<double>& y) const;
    std::vector<double> movingAverage(const std::vector<double>& data, int window) const;
    
    // File I/O helpers
    void writeCSV(const std::string& filename, 
                  const std::vector<std::string>& headers,
                  const std::vector<std::vector<double>>& data) const;
};

#endif
