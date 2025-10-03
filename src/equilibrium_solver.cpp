#include "equilibrium_solver.hh"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

// =============================================================================
// CONSTRUCTOR
// =============================================================================
EquilibriumSolver::EquilibriumSolver(const ZPinchParameters& params)
    : params_(params) {
    initializeGrid();
    
    // Initialize profiles
    p_eq_.resize(params_.Nr, 0.0);
    Bz_eq_.resize(params_.Nr, params_.B0);
    Btheta_eq_.resize(params_.Nr, 0.0);
    jz_eq_.resize(params_.Nr, 0.0);
    q_profile_.resize(params_.Nr, 0.0);
    beta_profile_.resize(params_.Nr, 0.0);
}

// =============================================================================
// GRID INITIALIZATION
// =============================================================================
void EquilibriumSolver::initializeGrid() {
    r_grid_.resize(params_.Nr);
    double dr = params_.a / (params_.Nr - 1);
    
    for (int i = 0; i < params_.Nr; ++i) {
        r_grid_[i] = i * dr;
    }
}
// =============================================================================
// FIELD PROFILE METHODS
// =============================================================================
double EquilibriumSolver::bthetaProfile(double r) const {
    // Campo azimutal: B_θ(r) = (μ₀ I₀ r) / (2π a²) para r < a
    // Asegurar dirección positiva para corriente positiva
    if (r <= params_.a) {
        return std::abs(4e-7 * params_.I0 * r) / (2.0 * params_.a * params_.a);
    } else {
        return std::abs(4e-7 * params_.I0) / (2.0 * r);
    }
}

// =============================================================================
// EQUILIBRIUM SOLVING METHODS
// =============================================================================
void EquilibriumSolver::solveGradShafranov() {
    std::cout << "Solving Grad-Shafranov equation for Z-pinch..." << std::endl;
    
    // For Z-pinch with symmetry, Grad-Shafranov simplifies to:
    // d²ψ/dr² + (1/r) dψ/dr = -μ₀ r j_φ(ψ)
    // where ψ is the flux function
    
    // We'll use a simple iterative approach
    computeBennettEquilibrium();  // Start with Bennett profile
    
    // Iterate to improve solution
    int max_iterations = 100;
    double tolerance = 1e-6;
    
    for (int iter = 0; iter < max_iterations; ++iter) {
        double max_residual = checkForceBalance();
        
        if (max_residual < tolerance) {
            std::cout << "Grad-Shafranov converged after " << iter << " iterations" << std::endl;
            break;
        }
        
        // Simple relaxation update
        computePressureFromForceBalance();
        computeBthetaFromCurrent();
        
        if (iter == max_iterations - 1) {
            std::cout << "Warning: Grad-Shafranov did not converge fully" << std::endl;
        }
    }
    
    computeSafetyFactor();
    computeBetaProfile();
}

void EquilibriumSolver::solveForceBalance() {
    std::cout << "Solving force balance equation..." << std::endl;
    
    // Force balance: dp/dr + (B_θ/μ₀r) d/dr(rB_θ) + (B_z/μ₀) dB_z/dr = 0
    
    computeBennettEquilibrium();  // Use Bennett as initial guess
    
    // Iterative solution of force balance
    int max_iterations = 50;
    
    for (int iter = 0; iter < max_iterations; ++iter) {
        // Update pressure from force balance
        computePressureFromForceBalance();
        
        // Update magnetic fields consistently
        computeBthetaFromCurrent();
        computeCurrentFromBtheta();
        
        applyBoundaryConditions();
    }
    
    computeSafetyFactor();
    computeBetaProfile();
}

// =============================================================================
// EQUILIBRIUM PROFILES
// =============================================================================
void EquilibriumSolver::computeBennettEquilibrium() {
    // Bennett equilibrium: p(r) = p0 / (1 + r²/a²)²
    // B_θ(r) = (μ₀ I₀ r) / (2π a²) for r < a
    
    double p0 = params_.n0 * params_.T0 * 1.6e-19; // Central pressure
    
    for (int i = 0; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        double r_a = r / params_.a;
        
        // Pressure profile
        p_eq_[i] = p0 / std::pow(1.0 + r_a * r_a, 2.0);
        
        // Azimuthal field profile
        if (r <= params_.a) {
            Btheta_eq_[i] = (4e-7 * M_PI * params_.I0 * r) / (2.0 * M_PI * params_.a * params_.a);
        } else {
            Btheta_eq_[i] = (4e-7 * M_PI * params_.I0) / (2.0 * M_PI * r);
        }
        
        // Axial field (constant for now)
        Bz_eq_[i] = params_.B0;
        
        // Current density from Ampere's law
        if (i == 0) {
            jz_eq_[i] = params_.I0 / (M_PI * params_.a * params_.a);
        } else {
            double r_plus = r_grid_[std::min(i+1, params_.Nr-1)];
            double r_minus = r_grid_[std::max(i-1, 0)];
            double dr = r_plus - r_minus;
            
            double d_rBtheta_dr = (r_plus * Btheta_eq_[std::min(i+1, params_.Nr-1)] - 
                                 r_minus * Btheta_eq_[std::max(i-1, 0)]) / dr;
            
            jz_eq_[i] = d_rBtheta_dr / (4e-7 * M_PI * r);
        }
    }
    
    applyBoundaryConditions();
}

void EquilibriumSolver::computeConstantCurrentEquilibrium() {
    // Constant current density: j_z = I₀ / (π a²) for r < a
    
    double j0 = params_.I0 / (M_PI * params_.a * params_.a);
    
    for (int i = 0; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        
        if (r <= params_.a) {
            jz_eq_[i] = j0;
            Btheta_eq_[i] = (4e-7 * M_PI * j0 * r) / 2.0;
        } else {
            jz_eq_[i] = 0.0;
            Btheta_eq_[i] = (4e-7 * M_PI * params_.I0) / (2.0 * M_PI * r);
        }
        
        // Compute pressure from force balance
        Bz_eq_[i] = params_.B0;
    }
    
    computePressureFromForceBalance();
    applyBoundaryConditions();
}

void EquilibriumSolver::computeParabolicEquilibrium() {
    // Parabolic profiles: p(r) = p0 (1 - r²/a²), j_z(r) = j0 (1 - r²/a²)
    
    double p0 = params_.n0 * params_.T0 * 1.6e-19;
    double j0 = 2.0 * params_.I0 / (M_PI * params_.a * params_.a); // Factor for parabolic profile
    
    for (int i = 0; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        double r_a = r / params_.a;
        
        if (r <= params_.a) {
            p_eq_[i] = p0 * (1.0 - r_a * r_a);
            jz_eq_[i] = j0 * (1.0 - r_a * r_a);
            Btheta_eq_[i] = (4e-7 * M_PI * j0 * r) * (1.0 - r_a * r_a / 2.0) / 2.0;
        } else {
            p_eq_[i] = 0.0;
            jz_eq_[i] = 0.0;
            Btheta_eq_[i] = (4e-7 * M_PI * params_.I0) / (2.0 * M_PI * r);
        }
        
        Bz_eq_[i] = params_.B0;
    }
    
    applyBoundaryConditions();
}

// =============================================================================
// FIELD COMPUTATION METHODS
// =============================================================================
void EquilibriumSolver::computeBthetaFromCurrent() {
    // B_θ(r) = (μ₀/r) ∫₀ʳ j_z(r') r' dr'
    
    for (int i = 0; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        if (r == 0.0) {
            Btheta_eq_[i] = 0.0;
            continue;
        }
        
        // Numerical integration using trapezoidal rule
        double integral = 0.0;
        for (int j = 1; j <= i; ++j) {
            double r_prev = r_grid_[j-1];
            double r_curr = r_grid_[j];
            double jz_prev = jz_eq_[j-1];
            double jz_curr = jz_eq_[j];
            
            integral += 0.5 * (jz_prev * r_prev + jz_curr * r_curr) * (r_curr - r_prev);
        }
        
        Btheta_eq_[i] = (4e-7 * M_PI * integral) / r;
    }
}

void EquilibriumSolver::computeCurrentFromBtheta() {
    // j_z(r) = (1/μ₀) (1/r) d/dr (r B_θ)
    
    for (int i = 0; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        
        if (i == 0) {
            // At r=0, use symmetric difference
            if (params_.Nr > 1) {
                double r1 = r_grid_[1];
                double f1 = r1 * Btheta_eq_[1];
                jz_eq_[0] = (2.0 * f1 / (r1 * r1)) / (4e-7 * M_PI);
            } else {
                jz_eq_[0] = 0.0;
            }
        } else {
            double dr = r_grid_[i] - r_grid_[i-1];
            double d_rBtheta_dr = (r * Btheta_eq_[i] - r_grid_[i-1] * Btheta_eq_[i-1]) / dr;
            jz_eq_[i] = d_rBtheta_dr / (4e-7 * M_PI * r);
        }
    }
}

void EquilibriumSolver::computePressureFromForceBalance() {
    // Force balance: dp/dr = - (B_θ/μ₀r) d/dr(rB_θ) - (B_z/μ₀) dB_z/dr
    
    // Assume B_z is constant for now, so dB_z/dr = 0
    p_eq_[0] = params_.n0 * params_.T0 * 1.6e-19; // Central pressure
    
    for (int i = 1; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        double dr = r - r_grid_[i-1];
        
        // Calculate derivative term
        double rBtheta_curr = r * Btheta_eq_[i];
        double rBtheta_prev = r_grid_[i-1] * Btheta_eq_[i-1];
        double d_rBtheta_dr = (rBtheta_curr - rBtheta_prev) / dr;
        
        // Pressure gradient
        double dp_dr = - (Btheta_eq_[i] / (4e-7 * M_PI * r)) * d_rBtheta_dr;
        
        // Update pressure
        p_eq_[i] = p_eq_[i-1] + dp_dr * dr;
        
        // Ensure non-negative pressure
        if (p_eq_[i] < 0.0) p_eq_[i] = 0.0;
    }
}

// =============================================================================
// DIAGNOSTIC METHODS
// =============================================================================
void EquilibriumSolver::computeSafetyFactor() {
    // Safety factor q(r) = (r B_z) / (R₀ B_θ) for cylindrical geometry
    
    for (int i = 0; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        if (r > 0.0 && std::abs(Btheta_eq_[i]) > 1e-12) {
            q_profile_[i] = (r * Bz_eq_[i]) / (params_.R0 * Btheta_eq_[i]);
        } else {
            q_profile_[i] = 1e12; // Large value at axis
        }
    }
}

void EquilibriumSolver::computeBetaProfile() {
    // Beta parameter: β(r) = 2μ₀ p(r) / B²(r)
    
    for (int i = 0; i < params_.Nr; ++i) {
        double B_sq = Bz_eq_[i] * Bz_eq_[i] + Btheta_eq_[i] * Btheta_eq_[i];
        if (B_sq > 1e-12) {
            beta_profile_[i] = 2.0 * 4e-7 * M_PI * p_eq_[i] / B_sq;
        } else {
            beta_profile_[i] = 0.0;
        }
    }
}

double EquilibriumSolver::computeTotalCurrent() const {
    // Total current: I = 2π ∫ j_z(r) r dr
    
    double total_current = 0.0;
    for (int i = 1; i < params_.Nr; ++i) {
        double r_prev = r_grid_[i-1];
        double r_curr = r_grid_[i];
        double jz_prev = jz_eq_[i-1];
        double jz_curr = jz_eq_[i];
        
        total_current += 0.5 * (jz_prev * r_prev + jz_curr * r_curr) * (r_curr - r_prev);
    }
    
    return 2.0 * M_PI * total_current;
}

double EquilibriumSolver::computeInternalInductance() const {
    // Internal inductance per unit length: l_i = (2/I²) ∫ (B_θ²/2μ₀) 2πr dr
    
    double magnetic_energy = 0.0;
    for (int i = 1; i < params_.Nr; ++i) {
        double r_prev = r_grid_[i-1];
        double r_curr = r_grid_[i];
        double Btheta_prev = Btheta_eq_[i-1];
        double Btheta_curr = Btheta_eq_[i];
        
        double energy_density_prev = Btheta_prev * Btheta_prev / (2.0 * 4e-7 * M_PI);
        double energy_density_curr = Btheta_curr * Btheta_curr / (2.0 * 4e-7 * M_PI);
        
        magnetic_energy += 0.5 * (energy_density_prev * r_prev + energy_density_curr * r_curr) * 
                          (r_curr - r_prev);
    }
    
    magnetic_energy *= 2.0 * M_PI;
    double I_total = computeTotalCurrent();
    
    return (2.0 * magnetic_energy) / (I_total * I_total);
}

double EquilibriumSolver::checkForceBalance() const {
    // Check how well force balance is satisfied
    // Returns maximum relative error
    
    double max_error = 0.0;
    
    for (int i = 1; i < params_.Nr - 1; ++i) {
        double r = r_grid_[i];
        double dr = r_grid_[i+1] - r_grid_[i];
        
        // Pressure gradient
        double dp_dr = (p_eq_[i+1] - p_eq_[i-1]) / (2.0 * dr);
        
        // Magnetic force term
        double rBtheta_plus = r_grid_[i+1] * Btheta_eq_[i+1];
        double rBtheta_minus = r_grid_[i-1] * Btheta_eq_[i-1];
        double d_rBtheta_dr = (rBtheta_plus - rBtheta_minus) / (2.0 * dr);
        
        double magnetic_force = (Btheta_eq_[i] / (4e-7 * M_PI * r)) * d_rBtheta_dr;
        
        // Total force (should be zero)
        double total_force = dp_dr + magnetic_force;
        double relative_error = std::abs(total_force) / (std::abs(dp_dr) + 1e-12);
        
        if (relative_error > max_error) {
            max_error = relative_error;
        }
    }
    
    return max_error;
}

// =============================================================================
// BOUNDARY CONDITIONS AND UTILITIES
// =============================================================================
void EquilibriumSolver::applyBoundaryConditions() {
    // Axis boundary conditions (r = 0)
    if (params_.Nr > 0) {
        // Pressure and B_z are finite at axis
        // B_θ = 0 at axis
        Btheta_eq_[0] = 0.0;
        
        // Current density finite at axis
        if (params_.Nr > 1) {
            jz_eq_[0] = jz_eq_[1]; // Approximate
        }
    }
    
    // Edge boundary (r = a)
    // Pressure goes to zero at plasma edge
    for (int i = 0; i < params_.Nr; ++i) {
        if (r_grid_[i] > params_.a) {
            p_eq_[i] = 0.0;
            jz_eq_[i] = 0.0;
        }
    }
}

double EquilibriumSolver::derivative(const std::vector<double>& f, int i) const {
    // Central difference for interior points
    if (i > 0 && i < params_.Nr - 1) {
        return (f[i+1] - f[i-1]) / (r_grid_[i+1] - r_grid_[i-1]);
    }
    // Forward/backward difference for boundaries
    else if (i == 0 && params_.Nr > 1) {
        return (f[1] - f[0]) / (r_grid_[1] - r_grid_[0]);
    }
    else if (i == params_.Nr - 1 && params_.Nr > 1) {
        return (f[params_.Nr-1] - f[params_.Nr-2]) / (r_grid_[params_.Nr-1] - r_grid_[params_.Nr-2]);
    }
    return 0.0;
}

double EquilibriumSolver::secondDerivative(const std::vector<double>& f, int i) const {
    if (i > 0 && i < params_.Nr - 1) {
        double h_plus = r_grid_[i+1] - r_grid_[i];
        double h_minus = r_grid_[i] - r_grid_[i-1];
        return 2.0 * (f[i+1] - f[i] - (f[i] - f[i-1]) * h_plus / h_minus) / 
               (h_plus * (h_plus + h_minus));
    }
    return 0.0;
}

// =============================================================================
// OUTPUT METHODS
// =============================================================================
void EquilibriumSolver::printEquilibriumSummary() const {
    std::cout << "=== Z-PINCH EQUILIBRIUM SUMMARY ===" << std::endl;
    std::cout << "Plasma radius: " << params_.a << " m" << std::endl;
    std::cout << "Total current: " << computeTotalCurrent() << " A" << std::endl;
    std::cout << "Central pressure: " << p_eq_[0] << " Pa" << std::endl;
    std::cout << "Central beta: " << beta_profile_[0] << std::endl;
    std::cout << "Safety factor at edge: " << q_profile_.back() << std::endl;
    std::cout << "Force balance error: " << checkForceBalance() << std::endl;
    std::cout << "Internal inductance: " << computeInternalInductance() << std::endl;
    std::cout << "===================================" << std::endl;
}

void EquilibriumSolver::saveProfiles(const std::string& filename) const {
    std::ofstream file(filename);
    file << "r,p,Bz,Btheta,jz,q,beta\n";
    file << std::scientific << std::setprecision(6);
    
    for (int i = 0; i < params_.Nr; ++i) {
        file << r_grid_[i] << ","
             << p_eq_[i] << ","
             << Bz_eq_[i] << ","
             << Btheta_eq_[i] << ","
             << jz_eq_[i] << ","
             << q_profile_[i] << ","
             << beta_profile_[i] << "\n";
    }
    
    file.close();
    std::cout << "Equilibrium profiles saved to: " << filename << std::endl;
}
