#ifndef ZPINCH_PARAMS_HH
#define ZPINCH_PARAMS_HH

#include <vector>
#include <string>

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
    
    // Default constructor - typical laboratory parameters
    ZPinchParameters() :
        R0(0.1), L(1.0), a(0.05),
        n0(1e20),  // Mayor densidad
        T0(100.0), // Temperatura más realista (100 eV)
        B0(2.0),   // Campo magnético más fuerte
        I0(5e4),   // Corriente positiva y realista (50 kA)
        mi_over_me(1836.0), Z(1.0),
        Nr(64),    // Menos puntos para mayor estabilidad
        Nz(32), 
        Ntheta(16),
        dt(1e-10), // Paso temporal más pequeño
        t_max(1e-6),
        k_min(0.1), k_max(50.0), n_modes(10)
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

#endif
