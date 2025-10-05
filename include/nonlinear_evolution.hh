#ifndef NONLINEAR_EVOLUTION_HH
#define NONLINEAR_EVOLUTION_HH

#include "zpinch_params.hh"
#include "equilibrium_solver.hh"
#include "linear_stability.hh"
#include <vector>
#include <memory>
#include <string>
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <chrono>

class NonlinearEvolution {
private:
    ZPinchParameters params_;
    std::shared_ptr<EquilibriumSolver> equilibrium_;
    
    // State variables for evolution
    std::vector<double> r_grid_;
    std::vector<double> z_grid_;
    
    // Primary variables: [ρ, v_r, v_z, v_θ, p, B_r, B_z, B_θ]
    std::vector<std::vector<double>> state_curr_;
    std::vector<std::vector<double>> state_next_;
    std::vector<std::vector<double>> state_temp_;
    
    // Time tracking
    double time_;
    int step_count_;
    
    // Diagnostic arrays
    std::vector<double> time_history_;
    std::vector<double> energy_history_;
    std::vector<double> growth_rate_history_;
    std::vector<double> kink_amplitude_history_;

    // Numerical parameters
    double cfl_number_ = 0.3;
    bool enable_divB_cleaning_ = true;
    double divB_cleaning_coeff_ = 0.1;
    double artificial_viscosity_ = 1e-6;
    
public:
    // Constructor and destructor
    NonlinearEvolution(const ZPinchParameters& params, 
                      std::shared_ptr<EquilibriumSolver> equilibrium);
    ~NonlinearEvolution() = default;
    
    // Initialization methods
    void initializeFromEquilibrium();
    void addKinkPerturbation(double amplitude, double kz, int m);
    void addSausagePerturbation(double amplitude, double kz);
    void addRandomPerturbation(double amplitude);

    // Time evolution methods
    void evolve(double t_max, double output_interval);
    bool stepEuler(double dt);
    bool stepRK2(double dt);
    
    // Diagnostic methods
    double computeKineticEnergy() const;
    double computeMagneticEnergy() const;
    double computeInternalEnergy() const;
    double computeTotalEnergy() const;
    double computeGrowthRate() const;
    double computeKinkAmplitude() const;
    bool checkStateFinite() const;

    // Historical energy methods
    double computeKineticEnergyAtStep(size_t step) const;
    double computeMagneticEnergyAtStep(size_t step) const;
    double computeInternalEnergyAtStep(size_t step) const;

    // Output methods
    void saveState(const std::string& filename) const;
    void saveTimeHistory(const std::string& filename) const;
    
    // Getters
    const std::vector<std::vector<double>>& getState() const { return state_curr_; }
    const std::vector<double>& getTimeHistory() const { return time_history_; }
    const std::vector<double>& getEnergyHistory() const { return energy_history_; }
    const std::vector<double>& getGrowthRateHistory() const { return growth_rate_history_; }
    const std::vector<double>& getKinkAmplitudeHistory() const { return kink_amplitude_history_; }
    
private:
    // Internal methods
    void initializeGrid();
    double computeCFLTimeStep() const;
    void storeDiagnostics();
    void applyBoundaryConditions();
    std::vector<std::vector<double>> computeRHS(const std::vector<std::vector<double>>& state);
    void debugGrowthCalculation() const;  // Added debug method
    
    // State variable indexing
    enum StateVariable {
        RHO = 0, VR = 1, VZ = 2, VTHETA = 3, 
        P = 4, BR = 5, BZ = 6, BTHETA = 7
    };
    
    // Physics constants
    const double mu0_ = 4.0e-7 * M_PI;
    const double gamma_ = 5.0 / 3.0;
};

#endif
