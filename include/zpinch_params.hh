#ifndef ZPINCH_PARAMS_HH
#define ZPINCH_PARAMS_HH

#include <vector>
#include <string>
#include <iostream>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct ZPinchParameters {
    // =========================================================================
    // PLASMA PHYSICAL PARAMETERS
    // =========================================================================
    
    // Geometric configuration
    double R0;           // Z-pinch radius [m]
    double L;            // System length [m]
    double a;            // Plasma radius [m] (a ≤ R0)
    
    // Plasma parameters
    double n0;           // Peak electron density [m^-3]
    double T0;           // Temperature [eV]
    double B0;           // Axial magnetic field [T]
    double I0;           // Total current [A]
    
    // Species parameters
    double mi_over_me;   // Ion/electron mass ratio (default: 1836 for proton)
    double Z;            // Average ion charge
    
    // =========================================================================
    // SIMULATION PARAMETERS
    // =========================================================================
    
    // Spatial discretization
    int Nr;              // Points in radial direction
    int Nz;              // Points in axial direction  
    int Ntheta;          // Points in azimuthal direction (for kink modes)
    
    // Temporal discretization
    double dt;           // Time step [s]
    double t_max;        // Maximum simulation time [s]
    
    // Stability parameters
    double k_min;        // Minimum wavenumber for analysis
    double k_max;        // Maximum wavenumber
    int n_modes;         // Number of modes to analyze
    
    // =========================================================================
    // CONSTRUCTORS AND METHODS
    // =========================================================================
    
    // Default constructor - parameters for INSTABILITY
    ZPinchParameters() :
        R0(0.1), L(1.0), a(0.08),
        n0(5e18),        // LOWER density for greater instability
        T0(20.0),        // LOWER temperature  
        B0(0.05),        // LOWER axial field
        I0(3e5),         // HIGHER current (300 kA)
        mi_over_me(1836.0), 
        Z(1.0),
        Nr(64), 
        Nz(32), 
        Ntheta(16),
        dt(1e-10), 
        t_max(1e-5),     // Longer time to see growth
        k_min(0.1), 
        k_max(50.0), 
        n_modes(10)
    {}
    
    // Method to calculate derived parameters
    double plasma_frequency() const {
        return sqrt(n0 * Z * Z * 1.6e-19 * 1.6e-19 / (8.854e-12 * 9.1e-31));
    }
    
    double alfven_speed() const {
        double rho = n0 * (1.67e-27 * mi_over_me + 9.1e-31); // mass density
        return B0 / sqrt(4e-7 * M_PI * rho);
    }
    
    double beta() const {
        double p = 2.0 * n0 * T0 * 1.6e-19; // pressure (Pa)
        double p_mag = B0 * B0 / (2.0 * 4e-7 * M_PI);
        return p / p_mag;
    }
    
    // Parameter validation
    bool validate() const {
        if (a > R0) return false;
        if (dt <= 0.0) return false;
        if (Nr <= 0 || Nz <= 0 || Ntheta <= 0) return false;
        return true;
    }
};

// =============================================================================
// STABILITY CHECK FUNCTION
// =============================================================================
inline void checkStabilityCriteria(const ZPinchParameters& params) {
    const double mu0 = 4e-7 * M_PI;
    double Btheta_surface = (mu0 * params.I0) / (2.0 * M_PI * params.a);
    double q_edge = (params.a * params.B0) / (params.R0 * Btheta_surface);
    
    double I_KS = (2.0 * M_PI * params.R0 * params.B0) / mu0;
    double kink_margin = I_KS / params.I0;
    
    std::cout << "=== STABILITY CRITERIA ===" << std::endl;
    std::cout << "Safety factor q: " << q_edge << std::endl;
    std::cout << "Kruskal-Shafranov current: " << I_KS << " A" << std::endl;
    std::cout << "Kink margin: " << kink_margin << std::endl;
    std::cout << "Plasma beta: " << params.beta() << std::endl;
    
    if (q_edge < 1.0) {
        std::cout << "→ UNSTABLE to kink modes (q < 1)" << std::endl;
    } else {
        std::cout << "→ STABLE to kink modes (q > 1)" << std::endl;
    }
    
    if (kink_margin < 1.0) {
        std::cout << "→ UNSTABLE by Kruskal-Shafranov (I > I_KS)" << std::endl;
    } else {
        std::cout << "→ STABLE by Kruskal-Shafranov (I < I_KS)" << std::endl;
    }
    
    // Additional criterion for sausage mode
    double critical_beta = 0.1; // Typical value
    if (params.beta() > critical_beta) {
        std::cout << "→ POTENTIALLY UNSTABLE to sausage modes (high beta)" << std::endl;
    } else {
        std::cout << "→ STABLE to sausage modes (low beta)" << std::endl;
    }
}

#endif
