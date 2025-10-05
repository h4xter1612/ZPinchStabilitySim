#ifndef EQUILIBRIUM_SOLVER_HH
#define EQUILIBRIUM_SOLVER_HH

#include "zpinch_params.hh"
#include <vector>
#include <memory>

class EquilibriumSolver {
private:
    ZPinchParameters params_;
    
    // Grid and equilibrium profiles
    std::vector<double> r_grid_;
    std::vector<double> p_eq_;          // Equilibrium pressure
    std::vector<double> Bz_eq_;         // Equilibrium axial field
    std::vector<double> Btheta_eq_;     // Equilibrium azimuthal field  
    std::vector<double> jz_eq_;         // Equilibrium current density
    
    // Diagnostic quantities
    std::vector<double> q_profile_;     // Safety factor profile
    std::vector<double> beta_profile_;  // Beta parameter profile
    
public:
    // =========================================================================
    // CONSTRUCTOR AND DESTRUCTOR
    // =========================================================================
    EquilibriumSolver(const ZPinchParameters& params);
    ~EquilibriumSolver() = default;
    
    // =========================================================================
    // EQUILIBRIUM SOLVING METHODS
    // =========================================================================
    
    // Solves the Grad-Shafranov equation for Z-pinch equilibrium
    void solveGradShafranov();
    
    // Solves force balance equation directly
    void solveForceBalance();
    
    // =========================================================================
    // EQUILIBRIUM PROFILES
    // =========================================================================
    
    // Bennett equilibrium profile (analytical)
    void computeBennettEquilibrium();
    
    // Constant current density profile
    void computeConstantCurrentEquilibrium();
    
    // Parabolic pressure profile
    void computeParabolicEquilibrium();
    
    // =========================================================================
    // PROFILE GETTERS
    // =========================================================================
    const std::vector<double>& getRadialGrid() const { return r_grid_; }
    const std::vector<double>& getPressureProfile() const { return p_eq_; }
    const std::vector<double>& getBzProfile() const { return Bz_eq_; }
    const std::vector<double>& getBthetaProfile() const { return Btheta_eq_; }
    const std::vector<double>& getCurrentProfile() const { return jz_eq_; }
    const std::vector<double>& getSafetyFactorProfile() const { return q_profile_; }
    const std::vector<double>& getBetaProfile() const { return beta_profile_; }
    
    // =========================================================================
    // DIAGNOSTIC METHODS
    // =========================================================================
    
    // Calculate safety factor q(r)
    void computeSafetyFactor();
    
    // Calculate beta parameter β(r)
    void computeBetaProfile();
    
    // Calculate total plasma current
    double computeTotalCurrent() const;
    
    // Calculate plasma internal inductance
    double computeInternalInductance() const;
    
    // Check force balance residual
    double checkForceBalance() const;
    
    // =========================================================================
    // OUTPUT METHODS
    // =========================================================================
    void printEquilibriumSummary() const;
    void saveProfiles(const std::string& filename) const;
    
private:
    // =========================================================================
    // INTERNAL METHODS
    // =========================================================================
    
    // Initialize radial grid
    void initializeGrid();
    
    // Compute B_theta from current density using Ampere's law
    void computeBthetaFromCurrent();
    
    // Compute current density from B_theta
    void computeCurrentFromBtheta();
    
    // Compute pressure from force balance
    void computePressureFromForceBalance();
    
    // Boundary conditions
    void applyBoundaryConditions();
};

#endif
