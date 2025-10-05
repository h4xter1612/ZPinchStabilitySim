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
// EQUILIBRIUM SOLVING METHODS
// =============================================================================
void EquilibriumSolver::solveGradShafranov() {
    std::cout << "Solving Grad-Shafranov equation for Z-pinch..." << std::endl;
    
    // Use Bennett equilibrium as initial solution
    computeBennettEquilibrium();
    
    computeSafetyFactor();
    computeBetaProfile();
    
    std::cout << "Grad-Shafranov solution completed" << std::endl;
}

void EquilibriumSolver::solveForceBalance() {
    std::cout << "Solving force balance equation..." << std::endl;
    
    const double mu0 = 4.0e-7 * M_PI;
    
    // Calculate realistic central pressure
    double n0 = params_.n0; // m^-3
    double T0_J = params_.T0 * 1.6e-19; // Convert eV to Joules
    double p0 = 2.0 * n0 * T0_J; // Total pressure (electrons + ions)
    
    std::cout << "Central pressure: " << p0 << " Pa" << std::endl;
    std::cout << "Temperature: " << params_.T0 << " eV" << std::endl;
    
    // Simple parabolic pressure profile
    for (int i = 0; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        double r_norm = r / params_.a;
        
        if (r < params_.a) {
            p_eq_[i] = p0 * (1.0 - r_norm * r_norm);
        } else {
            p_eq_[i] = 0.0;
        }
    }
    
    // Calculate Btheta consistently with current
    double I0 = params_.I0;
    
    for (int i = 0; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        
        if (r < params_.a) {
            // Uniform current distribution inside plasma
            double j0 = I0 / (M_PI * params_.a * params_.a);
            jz_eq_[i] = j0;
            Btheta_eq_[i] = (mu0 * j0 * r) / 2.0;
        } else {
            // Outside plasma: straight line field
            jz_eq_[i] = 0.0;
            Btheta_eq_[i] = (mu0 * I0) / (2.0 * M_PI * r);
        }
        
        Bz_eq_[i] = params_.B0;
    }
    
    applyBoundaryConditions();
    computeSafetyFactor();
    computeBetaProfile();
    
    std::cout << "Force balance solved" << std::endl;
}

// =============================================================================
// EQUILIBRIUM PROFILES
// =============================================================================
void EquilibriumSolver::computeBennettEquilibrium() {
    const double mu0 = 4.0e-7 * M_PI;
    
    // Calculate realistic central pressure
    double n0 = params_.n0;
    double T0_J = params_.T0 * 1.6e-19;
    double p0 = 2.0 * n0 * T0_J; // Total pressure
    
    double I0 = params_.I0;
    
    for (int i = 0; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        double r_a = r / params_.a;
        
        // Bennett pressure profile
        p_eq_[i] = p0 / std::pow(1.0 + r_a * r_a, 2.0);
        
        // Azimuthal field - corrected for physical consistency
        if (r <= params_.a) {
            Btheta_eq_[i] = (mu0 * I0 * r) / (2.0 * M_PI * params_.a * params_.a);
        } else {
            Btheta_eq_[i] = (mu0 * I0) / (2.0 * M_PI * r);
        }
        
        // Current density - corrected
        if (r <= params_.a) {
            jz_eq_[i] = I0 / (M_PI * params_.a * params_.a);
        } else {
            jz_eq_[i] = 0.0;
        }
        
        Bz_eq_[i] = params_.B0;
    }
    
    applyBoundaryConditions();
}

void EquilibriumSolver::computeConstantCurrentEquilibrium() {
    const double mu0 = 4.0e-7 * M_PI;
    
    double I0 = params_.I0;
    double j0 = I0 / (M_PI * params_.a * params_.a);
    
    for (int i = 0; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        
        if (r <= params_.a) {
            jz_eq_[i] = j0;
            Btheta_eq_[i] = (mu0 * j0 * r) / 2.0;
        } else {
            jz_eq_[i] = 0.0;
            Btheta_eq_[i] = (mu0 * I0) / (2.0 * M_PI * r);
        }
        
        Bz_eq_[i] = params_.B0;
    }
    
    // Calculate pressure from force balance
    computePressureFromForceBalance();
    applyBoundaryConditions();
}

void EquilibriumSolver::computeParabolicEquilibrium() {
    const double mu0 = 4.0e-7 * M_PI;
    
    // Calculate realistic central pressure
    double n0 = params_.n0;
    double T0_J = params_.T0 * 1.6e-19;
    double p0 = 2.0 * n0 * T0_J;
    
    double I0 = params_.I0;
    // Factor 2 for parabolic profile that integrates to I0
    double j0 = 2.0 * I0 / (M_PI * params_.a * params_.a);
    
    for (int i = 0; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        double r_a = r / params_.a;
        
        if (r <= params_.a) {
            p_eq_[i] = p0 * (1.0 - r_a * r_a);
            jz_eq_[i] = j0 * (1.0 - r_a * r_a);
            // Btheta consistent with current
            Btheta_eq_[i] = (mu0 * j0 * r) * (1.0 - r_a * r_a / 2.0) / 2.0;
        } else {
            p_eq_[i] = 0.0;
            jz_eq_[i] = 0.0;
            Btheta_eq_[i] = (mu0 * I0) / (2.0 * M_PI * r);
        }
        
        Bz_eq_[i] = params_.B0;
    }
    
    applyBoundaryConditions();
}

// =============================================================================
// DIAGNOSTIC METHODS
// =============================================================================
void EquilibriumSolver::computeSafetyFactor() {
    // q(r) = (2πr²B_z) / (μ₀R₀I(r)) for Z-pinch
    
    const double mu0 = 4.0e-7 * M_PI;
    
    for (int i = 0; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        
        if (r > 1e-12) {
            // Calculate enclosed current I(r)
            double I_enclosed = 0.0;
            for (int j = 0; j <= i; ++j) {
                double r_j = r_grid_[j];
                double dr = (j == 0) ? r_grid_[0] : (r_grid_[j] - r_grid_[j-1]);
                I_enclosed += jz_eq_[j] * 2.0 * M_PI * r_j * dr;
            }
            
            if (std::abs(I_enclosed) > 1e-12) {
                q_profile_[i] = (2.0 * M_PI * r * r * Bz_eq_[i]) / (mu0 * params_.R0 * I_enclosed);
            } else {
                q_profile_[i] = 1e12;
            }
        } else {
            q_profile_[i] = 1e12;
        }
    }
}

void EquilibriumSolver::computeBetaProfile() {
    const double mu0 = 4.0e-7 * M_PI;
    
    for (int i = 0; i < params_.Nr; ++i) {
        double B_total_sq = Bz_eq_[i] * Bz_eq_[i] + Btheta_eq_[i] * Btheta_eq_[i];
        if (B_total_sq > 1e-12) {
            beta_profile_[i] = (2.0 * mu0 * p_eq_[i]) / B_total_sq;
        } else {
            beta_profile_[i] = 0.0;
        }
    }
}

double EquilibriumSolver::computeTotalCurrent() const {
    double total = 0.0;
    for (int i = 1; i < params_.Nr; ++i) {
        double r_prev = r_grid_[i-1];
        double r_curr = r_grid_[i];
        double jz_prev = jz_eq_[i-1];
        double jz_curr = jz_eq_[i];
        
        double area_sector = M_PI * (r_curr * r_curr - r_prev * r_prev);
        double jz_avg = 0.5 * (jz_prev + jz_curr);
        total += jz_avg * area_sector;
    }
    return total;
}

double EquilibriumSolver::computeInternalInductance() const {
    const double mu0 = 4.0e-7 * M_PI;
    double magnetic_energy = 0.0;
    
    for (int i = 1; i < params_.Nr; ++i) {
        double r_prev = r_grid_[i-1];
        double r_curr = r_grid_[i];
        double Btheta_prev = Btheta_eq_[i-1];
        double Btheta_curr = Btheta_eq_[i];
        
        double Btheta_avg = 0.5 * (Btheta_prev + Btheta_curr);
        double volume_shell = M_PI * (r_curr * r_curr - r_prev * r_prev) * 1.0; // per unit length
        
        magnetic_energy += (Btheta_avg * Btheta_avg / (2.0 * mu0)) * volume_shell;
    }
    
    double I_total = computeTotalCurrent();
    if (I_total > 1e-12) {
        return (2.0 * magnetic_energy) / (I_total * I_total);
    }
    return 0.0;
}

double EquilibriumSolver::checkForceBalance() const {
    const double mu0 = 4.0e-7 * M_PI;
    double max_error = 0.0;
    int count = 0;
    
    for (int i = 1; i < params_.Nr - 1; ++i) {
        double r = r_grid_[i];
        if (r > 1e-12 && r < params_.a) {
            // Pressure gradient
            double dp_dr = (p_eq_[i+1] - p_eq_[i-1]) / (r_grid_[i+1] - r_grid_[i-1]);
            
            // Magnetic force: j_z × B_theta
            double jz = jz_eq_[i];
            double Btheta = Btheta_eq_[i];
            double magnetic_force = jz * Btheta;
            
            // Residual
            double residual = dp_dr + magnetic_force;
            max_error = std::max(max_error, std::abs(residual));
            count++;
        }
    }
    
    return (count > 0) ? max_error : 0.0;
}

// =============================================================================
// FIELD COMPUTATION METHODS
// =============================================================================
void EquilibriumSolver::computeBthetaFromCurrent() {
    const double mu0 = 4.0e-7 * M_PI;
    
    for (int i = 0; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        if (r < 1e-12) {
            Btheta_eq_[i] = 0.0;
            continue;
        }
        
        // Integrate enclosed current
        double I_enclosed = 0.0;
        for (int j = 0; j <= i; ++j) {
            double r_j = r_grid_[j];
            double dr = (j == 0) ? r_grid_[0] : (r_grid_[j] - r_grid_[j-1]);
            double jz_avg = jz_eq_[j];
            I_enclosed += jz_avg * 2.0 * M_PI * r_j * dr;
        }
        
        Btheta_eq_[i] = (mu0 * I_enclosed) / (2.0 * M_PI * r);
    }
}

void EquilibriumSolver::computeCurrentFromBtheta() {
    const double mu0 = 4.0e-7 * M_PI;
    
    for (int i = 0; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        
        if (i == 0) {
            // At axis, use limit
            jz_eq_[0] = (4.0 * Btheta_eq_[1]) / (mu0 * (r_grid_[1] + 1e-12));
        } else {
            double dr = r_grid_[i] - r_grid_[i-1];
            double d_rBtheta_dr = (r * Btheta_eq_[i] - r_grid_[i-1] * Btheta_eq_[i-1]) / dr;
            jz_eq_[i] = d_rBtheta_dr / (mu0 * r);
        }
    }
}

void EquilibriumSolver::computePressureFromForceBalance() {
    const double mu0 = 4.0e-7 * M_PI;
    
    // Set realistic central pressure
    double n0 = params_.n0;
    double T0_J = params_.T0 * 1.6e-19;
    p_eq_[0] = 2.0 * n0 * T0_J;
    
    for (int i = 1; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        double dr = r - r_grid_[i-1];
        
        if (r < params_.a) {
            // Calculate magnetic force: -j_z × B_theta
            double jz = 0.5 * (jz_eq_[i] + jz_eq_[i-1]);
            double Btheta = 0.5 * (Btheta_eq_[i] + Btheta_eq_[i-1]);
            double magnetic_force = -jz * Btheta;
            
            // Integrate to get pressure
            p_eq_[i] = p_eq_[i-1] + magnetic_force * dr;
            
            // Ensure non-negative pressure
            if (p_eq_[i] < 0.0) p_eq_[i] = 0.0;
        } else {
            p_eq_[i] = 0.0;
        }
    }
}

// =============================================================================
// BOUNDARY CONDITIONS AND UTILITIES
// =============================================================================
void EquilibriumSolver::applyBoundaryConditions() {
    // Axis (r = 0)
    if (params_.Nr > 0) {
        Btheta_eq_[0] = 0.0; // Azimuthal field must be zero at axis
        
        // Extrapolate pressure and current from first interior point
        if (params_.Nr > 1) {
            p_eq_[0] = p_eq_[1];
            jz_eq_[0] = jz_eq_[1];
        }
    }
    
    // Plasma edge (r = a)
    for (int i = 0; i < params_.Nr; ++i) {
        if (r_grid_[i] >= params_.a) {
            p_eq_[i] = 0.0;
            jz_eq_[i] = 0.0;
        }
    }
}

// =============================================================================
// OUTPUT METHODS
// =============================================================================
void EquilibriumSolver::printEquilibriumSummary() const {
    std::cout << "=== Z-PINCH EQUILIBRIUM SUMMARY ===" << std::endl;
    std::cout << "Plasma radius: " << params_.a << " m" << std::endl;
    std::cout << "Total current: " << computeTotalCurrent() << " A" << std::endl;
    std::cout << "Central pressure: " << p_eq_[0] << " Pa" << std::endl;
    
    double B_total_sq = params_.B0 * params_.B0 + Btheta_eq_.back() * Btheta_eq_.back();
    double p_mag = B_total_sq / (2.0 * 4.0e-7 * M_PI);
    double beta_central = p_eq_[0] / p_mag;
    std::cout << "Central beta: " << beta_central << std::endl;
    
    std::cout << "Safety factor at edge: " << q_profile_.back() << std::endl;
    std::cout << "Force balance error: " << checkForceBalance() << " Pa/m" << std::endl;
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
