#ifndef NONLINEAR_EVOLUTION_HH
#define NONLINEAR_EVOLUTION_HH

#include "zpinch_params.hh"
#include "equilibrium_solver.hh"
#include "linear_stability.hh"
#include <vector>
#include <memory>
#include <functional>

class NonlinearEvolution {
private:
    ZPinchParameters params_;
    std::shared_ptr<EquilibriumSolver> equilibrium_;
    
    // State variables for evolution
    std::vector<double> r_grid_;
    std::vector<double> z_grid_;
    
    // Primary variables: [ρ, v_r, v_z, v_θ, p, B_r, B_z, B_θ]
    std::vector<std::vector<double>> state_curr_;  // Current state
    std::vector<std::vector<double>> state_next_;  // Next state
    std::vector<std::vector<double>> state_temp_;  // Temporary state for RK steps
    
    // Perturbation fields
    std::vector<std::vector<double>> perturbation_;
    
    // Time tracking
    double time_;
    int step_count_;
    
    // Diagnostic arrays
    std::vector<double> time_history_;
    std::vector<double> energy_history_;
    std::vector<double> growth_rate_history_;
    
public:
    // =========================================================================
    // CONSTRUCTOR AND DESTRUCTOR
    // =========================================================================
    NonlinearEvolution(const ZPinchParameters& params, 
                      std::shared_ptr<EquilibriumSolver> equilibrium);
    ~NonlinearEvolution() = default;
    
    // =========================================================================
    // INITIALIZATION METHODS
    // =========================================================================
    
    // Initialize from equilibrium solution
    void initializeFromEquilibrium();
    
    // Add perturbation for specific instability mode
    void addKinkPerturbation(double amplitude, double kz, int m);
    void addSausagePerturbation(double amplitude, double kz);
    void addRandomPerturbation(double amplitude);
    
    // =========================================================================
    // TIME EVOLUTION METHODS
    // =========================================================================
    
    // Main evolution driver
    void evolve(double t_max, double output_interval);
    
    // Single time step using split-step method
    void stepSplitStep(double dt);
    
    // Alternative time integrators (for comparison)
    void stepRK4(double dt);
    void stepPredictorCorrector(double dt);
    
    // =========================================================================
    // PHYSICS METHODS - SPLIT-STEP COMPONENTS
    // =========================================================================
    
    // Hyperbolic step (ideal MHD terms)
    void stepHyperbolic(double dt);
    
    // Source terms step (pressure gradient, Lorentz force)
    void stepSourceTerms(double dt);
    
    // Divergence cleaning step (maintain ∇·B = 0)
    void stepDivergenceCleaning(double dt);
    
    // Boundary conditions application
    void applyBoundaryConditions();
    
    // =========================================================================
    // RIGHT-HAND SIDE COMPUTATIONS
    // =========================================================================
    
    // Compute RHS for MHD equations
    void computeRHS(const std::vector<std::vector<double>>& state,
                   std::vector<std::vector<double>>& rhs);
    
    // Individual term computations
    void computeAdvectionTerms(const std::vector<std::vector<double>>& state,
                              std::vector<std::vector<double>>& advection);
    void computePressureTerms(const std::vector<std::vector<double>>& state,
                             std::vector<std::vector<double>>& pressure);
    void computeLorentzTerms(const std::vector<std::vector<double>>& state,
                            std::vector<std::vector<double>>& lorentz);
    
    // =========================================================================
    // DIAGNOSTIC METHODS
    // =========================================================================
    
    // Energy calculations
    double computeKineticEnergy() const;
    double computeMagneticEnergy() const;
    double computeInternalEnergy() const;
    double computeTotalEnergy() const;
    
    // Instability growth rates
    double computeGrowthRate() const;
    
    // Conservation checks
    double checkMassConservation() const;
    double checkEnergyConservation() const;
    double checkDivBCondition() const;
    
    // =========================================================================
    // OUTPUT AND VISUALIZATION
    // =========================================================================
    
    // Save current state to file
    void saveState(const std::string& filename) const;
    
    // Save time history
    void saveTimeHistory(const std::string& filename) const;
    
    // Get state for visualization
    const std::vector<std::vector<double>>& getState() const { return state_curr_; }
    const std::vector<double>& getTimeHistory() const { return time_history_; }
    const std::vector<double>& getEnergyHistory() const { return energy_history_; }
    
    // =========================================================================
    // PARAMETER CONTROL
    // =========================================================================
    void setCFLNumber(double cfl) { cfl_number_ = cfl; }
    void setDivBCleaning(bool enable) { enable_divB_cleaning_ = enable; }

    void addArtificialDissipation(double dt);
    
private:
    // =========================================================================
    // INTERNAL METHODS
    // =========================================================================
    
    // Grid initialization
    void initializeGrid();
    
    // State variable indexing
    enum StateVariable {
        RHO = 0,    // Density
        VR = 1,     // Radial velocity
        VZ = 2,     // Axial velocity
        VTHETA = 3, // Azimuthal velocity
        P = 4,      // Pressure
        BR = 5,     // Radial magnetic field
        BZ = 6,     // Axial magnetic field
        BTHETA = 7  // Azimuthal magnetic field
    };
    
    // Numerical parameters
    double cfl_number_ = 0.5;
    bool enable_divB_cleaning_ = true;
    double divB_cleaning_coeff_ = 1.0;
    
    // Helper functions
    double computeSoundSpeed(double pressure, double density) const;
    double computeAlfvenSpeed(double B, double density) const;
    double computeCFLTimeStep() const;
    
    // Numerical derivatives
    double derivativeR(const std::vector<double>& f, int i, int j) const;
    double derivativeZ(const std::vector<double>& f, int i, int j) const;
    
    // Interpolation functions
    double interpolateToFaceR(const std::vector<double>& f, int i, int j) const;
    double interpolateToFaceZ(const std::vector<double>& f, int i, int j) const;
};

#endif
