#ifndef LINEAR_STABILITY_HH
#define LINEAR_STABILITY_HH

#include "zpinch_params.hh"
#include "equilibrium_solver.hh"
#include <vector>
#include <complex>
#include <memory>

class LinearStabilityAnalyzer {
private:
    ZPinchParameters params_;
    std::shared_ptr<EquilibriumSolver> equilibrium_;
    
    // Equilibrium fields
    std::vector<double> r_grid_;
    std::vector<double> p_eq_;      // Equilibrium pressure
    std::vector<double> Bz_eq_;     // Equilibrium axial field  
    std::vector<double> Btheta_eq_; // Equilibrium azimuthal field
    std::vector<double> jz_eq_;     // Axial current density
    
public:
    // =========================================================================
    // CONSTRUCTOR
    // =========================================================================
    LinearStabilityAnalyzer(const ZPinchParameters& params, std::shared_ptr<EquilibriumSolver> equilibrium);
    
    // =========================================================================
    // EQUILIBRIUM CALCULATION
    // =========================================================================
    
    // Calculates Z-pinch equilibrium (force balance j×B = ∇p)
    void computeEquilibrium();
    
    // Bennett pressure profile
    double bennettPressureProfile(double r) const;
    
    // Azimuthal field profile
    double bthetaProfile(double r) const;
    
    // =========================================================================
    // LINEAR STABILITY ANALYSIS
    // =========================================================================
    
    // Analyzes kink mode (m=1)
    struct KinkModeResult {
        double k;                    // Axial wavenumber
        std::complex<double> omega;  // Complex frequency
        double growth_rate;          // Growth rate (Im(omega))
        double frequency;            // Real frequency (Re(omega))
        std::vector<double> eigenfunction; // Radial eigenfunction
    };
    
    KinkModeResult analyzeKinkMode(double k);
    
    // Analyzes sausage mode (m=0)  
    struct SausageModeResult {
        double k;
        std::complex<double> omega;
        double growth_rate;
        double frequency;
        std::vector<double> eigenfunction;
    };
    
    SausageModeResult analyzeSausageMode(double k);
    
    // k-scan to find unstable modes
    std::vector<KinkModeResult> kinkStabilityScan();
    std::vector<SausageModeResult> sausageStabilityScan();
    
    // =========================================================================
    // HELPER METHODS
    // =========================================================================
    
    // Kruskal-Shafranov stability criterion
    double kruskalShafranovLimit() const;
    
    // Sausage mode stability criterion
    double sausageStabilityCriterion() const;
    
    // Get equilibrium grids and fields
    const std::vector<double>& getRadialGrid() const { return r_grid_; }
    const std::vector<double>& getPressureProfile() const { return p_eq_; }
    const std::vector<double>& getBthetaProfile() const { return Btheta_eq_; }
    
private:
    // =========================================================================
    // INTERNAL METHODS
    // =========================================================================
    
    // Initializes radial grid
    void initializeGrid();
    
    // Solves eigenvalue equation for given mode
    std::complex<double> solveEigenvalueEquation(int m, double k);
    
    // Computes eigenfunction for a mode
    std::vector<double> computeEigenfunction(int m, double k, 
                                           std::complex<double> omega);
};

#endif
