#include "nonlinear_evolution.hh"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <random>

// =============================================================================
// CONSTRUCTOR
// =============================================================================
NonlinearEvolution::NonlinearEvolution(const ZPinchParameters& params,
                                     std::shared_ptr<EquilibriumSolver> equilibrium)
    : params_(params), equilibrium_(equilibrium), time_(0.0), step_count_(0) {
    
    initializeGrid();
    initializeFromEquilibrium();
}

// =============================================================================
// GRID INITIALIZATION
// =============================================================================
void NonlinearEvolution::initializeGrid() {
    // Radial grid
    r_grid_.resize(params_.Nr);
    double dr = params_.a / (params_.Nr - 1);
    for (int i = 0; i < params_.Nr; ++i) {
        r_grid_[i] = i * dr;
    }
    
    // Axial grid (periodic boundary conditions)
    z_grid_.resize(params_.Nz);
    double dz = params_.L / params_.Nz;
    for (int j = 0; j < params_.Nz; ++j) {
        z_grid_[j] = j * dz;
    }
    
    // Initialize state arrays (8 variables × Nr × Nz)
    int n_vars = 8; // ρ, v_r, v_z, v_θ, p, B_r, B_z, B_θ
    state_curr_.resize(n_vars);
    state_next_.resize(n_vars);
    state_temp_.resize(n_vars);
    perturbation_.resize(n_vars);
    
    for (int var = 0; var < n_vars; ++var) {
        state_curr_[var].resize(params_.Nr * params_.Nz, 0.0);
        state_next_[var].resize(params_.Nr * params_.Nz, 0.0);
        state_temp_[var].resize(params_.Nr * params_.Nz, 0.0);
        perturbation_[var].resize(params_.Nr * params_.Nz, 0.0);
    }
}

// =============================================================================
// INITIALIZATION FROM EQUILIBRIUM
// =============================================================================
void NonlinearEvolution::initializeFromEquilibrium() {
    // Get equilibrium profiles
    const auto& r_eq = equilibrium_->getRadialGrid(); // Non used variable
    const auto& p_eq = equilibrium_->getPressureProfile();
    const auto& Bz_eq = equilibrium_->getBzProfile();
    const auto& Btheta_eq = equilibrium_->getBthetaProfile();
    
    // Initialize state with equilibrium values
    for (int i = 0; i < params_.Nr; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;
            
            // Density from pressure and temperature (ideal gas law)
            double T0 = params_.T0 * 1.6e-19; // Convert eV to J
            state_curr_[RHO][idx] = p_eq[i] / T0;
            
            // Velocities (zero in equilibrium)
            state_curr_[VR][idx] = 0.0;
            state_curr_[VZ][idx] = 0.0;
            state_curr_[VTHETA][idx] = 0.0;
            
            // Pressure
            state_curr_[P][idx] = p_eq[i];
            
            // Magnetic fields
            state_curr_[BR][idx] = 0.0; // No radial field in equilibrium
            state_curr_[BZ][idx] = Bz_eq[i];
            state_curr_[BTHETA][idx] = Btheta_eq[i];
        }
    }
    
    std::cout << "Initialized from equilibrium solution" << std::endl;
}

// =============================================================================
// PERTURBATION METHODS
// =============================================================================
void NonlinearEvolution::addKinkPerturbation(double amplitude, double kz, int m) {
    // Kink mode (m=1) perturbation
    // Perturbs v_r, v_θ, and B_r, B_θ
    
    for (int i = 0; i < params_.Nr; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;
            double r = r_grid_[i];
            double z = z_grid_[j];
            
            if (r < params_.a) { // Only perturb inside plasma
                // Velocity perturbation
                perturbation_[VR][idx] = amplitude * std::sin(m * M_PI * r / params_.a) * 
                                        std::cos(kz * z);
                perturbation_[VTHETA][idx] = amplitude * std::cos(m * M_PI * r / params_.a) * 
                                           std::sin(kz * z);
                
                // Magnetic field perturbation
                perturbation_[BR][idx] = 0.1 * amplitude * state_curr_[BTHETA][idx] * 
                                       std::sin(m * M_PI * r / params_.a) * std::sin(kz * z);
                perturbation_[BTHETA][idx] = 0.1 * amplitude * state_curr_[BTHETA][idx] * 
                                          std::cos(m * M_PI * r / params_.a) * std::cos(kz * z);
            }
        }
    }
    
    // Apply perturbation to current state
    for (int var = 0; var < 8; ++var) {
        for (int idx = 0; idx < params_.Nr * params_.Nz; ++idx) {
            state_curr_[var][idx] += perturbation_[var][idx];
        }
    }
    
    std::cout << "Added kink mode perturbation: amplitude=" << amplitude 
              << ", kz=" << kz << ", m=" << m << std::endl;
}
// =============================================================================
// ENERGY COMPUTATIONS - IMPLEMENTACIONES FALTANTES
// =============================================================================
double NonlinearEvolution::computeInternalEnergy() const {
    double ie = 0.0;
    double gamma = 5.0 / 3.0; // adiabatic index
    
    for (int i = 0; i < params_.Nr; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;
            double r = r_grid_[i];
            double dr = (i < params_.Nr - 1) ? (r_grid_[i+1] - r_grid_[i]) : 0.0;
            double dz = params_.L / params_.Nz;
            
            ie += (state_curr_[P][idx] / (gamma - 1.0)) * r * dr * dz;
        }
    }
    return ie;
}

double NonlinearEvolution::computeTotalEnergy() const {
    return computeKineticEnergy() + computeMagneticEnergy() + computeInternalEnergy();
}

double NonlinearEvolution::computeSoundSpeed(double pressure, double density) const {
    if (density <= 0.0) return 0.0;
    double gamma = 5.0 / 3.0; // adiabatic index for ideal gas
    return std::sqrt(gamma * pressure / density);
}

double NonlinearEvolution::computeAlfvenSpeed(double B, double density) const {
    if (density <= 0.0) return 0.0;
    const double mu0 = 4.0e-7 * M_PI;
    return B / std::sqrt(mu0 * density);
}

void NonlinearEvolution::addSausagePerturbation(double amplitude, double kz) {
    // Sausage mode (m=0) perturbation
    // Perturbs v_r, p, and B_θ
    
    for (int i = 0; i < params_.Nr; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;
            double r = r_grid_[i];
            double z = z_grid_[j];
            
            if (r < params_.a) {
                // Radial velocity perturbation (compressive)
                perturbation_[VR][idx] = amplitude * (r / params_.a) * 
                                        std::sin(kz * z);
                
                // Pressure perturbation
                perturbation_[P][idx] = 0.5 * amplitude * state_curr_[P][idx] * 
                                      (r / params_.a) * std::cos(kz * z);
                
                // Magnetic field perturbation
                perturbation_[BTHETA][idx] = 0.1 * amplitude * state_curr_[BTHETA][idx] * 
                                          (r / params_.a) * std::sin(kz * z);
            }
        }
    }
    
    // Apply perturbation
    for (int var = 0; var < 8; ++var) {
        for (int idx = 0; idx < params_.Nr * params_.Nz; ++idx) {
            state_curr_[var][idx] += perturbation_[var][idx];
        }
    }
    
    std::cout << "Added sausage mode perturbation: amplitude=" << amplitude 
              << ", kz=" << kz << std::endl;
}

void NonlinearEvolution::addRandomPerturbation(double amplitude) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-amplitude, amplitude);
    
    for (int i = 0; i < params_.Nr; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;
            
            if (r_grid_[i] < params_.a) {
                // Add random perturbations to velocities and magnetic fields
                state_curr_[VR][idx] += dis(gen) * 0.01 * params_.alfven_speed();
                state_curr_[VZ][idx] += dis(gen) * 0.01 * params_.alfven_speed();
                state_curr_[VTHETA][idx] += dis(gen) * 0.01 * params_.alfven_speed();
                
                state_curr_[BR][idx] += dis(gen) * 0.01 * params_.B0;
                state_curr_[BTHETA][idx] += dis(gen) * 0.01 * params_.B0;
            }
        }
    }
    
    std::cout << "Added random perturbation with amplitude: " << amplitude << std::endl;
}

// =============================================================================
// TIME EVOLUTION - MAIN DRIVER
// =============================================================================
void NonlinearEvolution::evolve(double t_max, double output_interval) {
    std::cout << "Starting nonlinear evolution to t = " << t_max << std::endl;
    std::cout << "Output interval: " << output_interval << std::endl;
    
    double next_output_time = output_interval;
    int output_count = 0;
    
    while (time_ < t_max) {
        // Compute adaptive time step
        double dt = computeCFLTimeStep();
        if (time_ + dt > t_max) {
            dt = t_max - time_;
        }
        
        // Perform time step
        stepSplitStep(dt);
        
        // Update time and step counter
        time_ += dt;
        step_count_++;
        
        // Store diagnostics
        time_history_.push_back(time_);
        energy_history_.push_back(computeTotalEnergy());
        
        // Output if needed
        if (time_ >= next_output_time || time_ >= t_max) {
            std::string filename = "data/snap/output_state_" + std::to_string(output_count) + ".csv";
            saveState(filename);
            
            std::cout << "Step " << step_count_ << ", t = " << time_ 
                      << ", dt = " << dt << ", Energy = " << computeTotalEnergy() << std::endl;
            
            next_output_time += output_interval;
            output_count++;
        }
        
        // Check for instability growth
        if (step_count_ % 100 == 0) {
            double growth_rate = computeGrowthRate();
            growth_rate_history_.push_back(growth_rate);
            std::cout << "Growth rate estimate: " << growth_rate << " s^-1" << std::endl;
        }
    }
    
    std::cout << "Evolution completed. Total steps: " << step_count_ << std::endl;
}

// =============================================================================
// SPLIT-STEP TIME INTEGRATION
// =============================================================================
void NonlinearEvolution::stepSplitStep(double dt) {
    // Strang splitting for MHD equations:
    // Step 1: Half step for hyperbolic terms (advection)
    // Step 2: Full step for source terms (pressure + Lorentz force)
    // Step 3: Half step for hyperbolic terms (advection)
    
    // Step 1: Half hyperbolic step
    stepHyperbolic(dt / 2.0);
    
    // Step 2: Full source terms step
    stepSourceTerms(dt);
    
    // Step 3: Half hyperbolic step
    stepHyperbolic(dt / 2.0);
    
    // Step 4: Divergence cleaning (if enabled)
    if (enable_divB_cleaning_) {
        stepDivergenceCleaning(dt);
    }
    
    // Apply boundary conditions
    applyBoundaryConditions();
}

void NonlinearEvolution::stepHyperbolic(double dt) {
    // Hyperbolic step: advection terms using conservative update
    std::vector<std::vector<double>> flux_contrib(8, 
        std::vector<double>(params_.Nr * params_.Nz, 0.0));
    
    computeAdvectionTerms(state_curr_, flux_contrib);
    
    // Update state using flux
    for (int var = 0; var < 8; ++var) {
        for (int idx = 0; idx < params_.Nr * params_.Nz; ++idx) {
            state_curr_[var][idx] += dt * flux_contrib[var][idx];
        }
    }
}

void NonlinearEvolution::stepSourceTerms(double dt) {
    // Source terms: pressure gradient + Lorentz force
    std::vector<std::vector<double>> source_contrib(8, 
        std::vector<double>(params_.Nr * params_.Nz, 0.0));
    
    // Compute pressure terms
    computePressureTerms(state_curr_, source_contrib);
    
    // Compute Lorentz force terms
    computeLorentzTerms(state_curr_, source_contrib);
    
    // Update state using source terms
    for (int var = 0; var < 8; ++var) {
        for (int idx = 0; idx < params_.Nr * params_.Nz; ++idx) {
            state_curr_[var][idx] += dt * source_contrib[var][idx];
        }
    }
}

void NonlinearEvolution::stepDivergenceCleaning(double dt) {
    // Simple hyperbolic divergence cleaning for ∇·B = 0
    // Uses a diffusive approach to damp numerical divergence
    
    for (int i = 1; i < params_.Nr - 1; ++i) {
        for (int j = 1; j < params_.Nz - 1; ++j) {
            int idx = i * params_.Nz + j;
            
            // Compute divergence of B
            double divB = derivativeR(state_curr_[BR], i, j) + 
                         derivativeZ(state_curr_[BZ], i, j) +
                         state_curr_[BR][idx] / r_grid_[i]; // Cylindrical coordinates
            
            // Apply correction to magnetic field
            double correction = divB_cleaning_coeff_ * divB * dt;
            state_curr_[BR][idx] -= correction * derivativeR(state_curr_[BR], i, j);
            state_curr_[BZ][idx] -= correction * derivativeZ(state_curr_[BZ], i, j);
        }
    }
}

// =============================================================================
// RIGHT-HAND SIDE COMPUTATIONS
// =============================================================================
void NonlinearEvolution::computeAdvectionTerms(const std::vector<std::vector<double>>& state,
                                              std::vector<std::vector<double>>& advection) {
    // Compute advection terms using conservative finite differences
    // This is a simplified implementation - real code would use more sophisticated methods
    
    for (int i = 1; i < params_.Nr - 1; ++i) {
        for (int j = 1; j < params_.Nz - 1; ++j) {
            int idx = i * params_.Nz + j;
            double r = r_grid_[i];
            double dr = r_grid_[i+1] - r_grid_[i];
            double dz = z_grid_[j+1] - z_grid_[j];
            
            // Mass conservation
            double flux_r = state[RHO][idx] * state[VR][idx];
            double flux_z = state[RHO][idx] * state[VZ][idx];
            
            advection[RHO][idx] = -(1.0/r) * (r * flux_r - (r - dr) * 
                             state[RHO][(i-1)*params_.Nz + j] * state[VR][(i-1)*params_.Nz + j]) / dr
                             - (flux_z - state[RHO][i*params_.Nz + j-1] * 
                             state[VZ][i*params_.Nz + j-1]) / dz;
            
            // Similar implementations for other variables...
            // Note: This is a simplified placeholder - full MHD advection is complex
        }
    }
}

void NonlinearEvolution::computePressureTerms(const std::vector<std::vector<double>>& state,
                                             std::vector<std::vector<double>>& pressure) {
    // Pressure gradient terms in momentum equations
    
    for (int i = 1; i < params_.Nr - 1; ++i) {
        for (int j = 1; j < params_.Nz - 1; ++j) {
            int idx = i * params_.Nz + j;
            
            // Radial momentum: -∂p/∂r
            pressure[VR][idx] = -derivativeR(state[P], i, j);
            
            // Axial momentum: -∂p/∂z  
            pressure[VZ][idx] = -derivativeZ(state[P], i, j);
            
            // Azimuthal momentum: no direct pressure term in cylindrical coords
            pressure[VTHETA][idx] = 0.0;
        }
    }
}

void NonlinearEvolution::computeLorentzTerms(const std::vector<std::vector<double>>& state,
                                            std::vector<std::vector<double>>& lorentz) {
    // Lorentz force: j × B = (∇ × B) × B / μ₀
    
    const double mu0 = 4.0e-7 * M_PI;
    
    for (int i = 1; i < params_.Nr - 1; ++i) {
        for (int j = 1; j < params_.Nz - 1; ++j) {
            int idx = i * params_.Nz + j;
            double r = r_grid_[i];
            
            // Compute current density j = ∇ × B / μ₀
            double j_r, j_z, j_theta;
            
            // j_r = -∂B_θ/∂z (in cylindrical coordinates)
            j_r = -derivativeZ(state[BTHETA], i, j);
            
            // j_z = (1/r) ∂(r B_θ)/∂r
            j_z = (1.0/r) * derivativeR(state[BTHETA], i, j) + state[BTHETA][idx] / r;
            
            // j_θ = ∂B_r/∂z - ∂B_z/∂r
            j_theta = derivativeZ(state[BR], i, j) - derivativeR(state[BZ], i, j);
            
            // Lorentz force: (j × B)
            lorentz[VR][idx] = (j_theta * state[BZ][idx] - j_z * state[BTHETA][idx]) / mu0;
            lorentz[VZ][idx] = (j_r * state[BTHETA][idx] - j_theta * state[BR][idx]) / mu0;
            lorentz[VTHETA][idx] = (j_z * state[BR][idx] - j_r * state[BZ][idx]) / mu0;
        }
    }
}

// =============================================================================
// DIAGNOSTIC METHODS
// =============================================================================
double NonlinearEvolution::computeKineticEnergy() const {
    double ke = 0.0;
    for (int i = 0; i < params_.Nr; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;
            double r = r_grid_[i];
            double dr = (i < params_.Nr - 1) ? (r_grid_[i+1] - r_grid_[i]) : 0.0;
            double dz = params_.L / params_.Nz;
            
            double v2 = state_curr_[VR][idx] * state_curr_[VR][idx] +
                       state_curr_[VZ][idx] * state_curr_[VZ][idx] +
                       state_curr_[VTHETA][idx] * state_curr_[VTHETA][idx];
            
            ke += 0.5 * state_curr_[RHO][idx] * v2 * r * dr * dz;
        }
    }
    return ke;
}

double NonlinearEvolution::computeMagneticEnergy() const {
    double me = 0.0;
    const double mu0 = 4.0e-7 * M_PI;
    
    for (int i = 0; i < params_.Nr; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;
            double r = r_grid_[i];
            double dr = (i < params_.Nr - 1) ? (r_grid_[i+1] - r_grid_[i]) : 0.0;
            double dz = params_.L / params_.Nz;
            
            double B2 = state_curr_[BR][idx] * state_curr_[BR][idx] +
                       state_curr_[BZ][idx] * state_curr_[BZ][idx] +
                       state_curr_[BTHETA][idx] * state_curr_[BTHETA][idx];
            
            me += (B2 / (2.0 * mu0)) * r * dr * dz;
        }
    }
    return me;
}

double NonlinearEvolution::computeGrowthRate() const {
    // Estimate growth rate from recent energy history
    if (energy_history_.size() < 10) return 0.0;
    
    int n_points = std::min(20, static_cast<int>(energy_history_.size()));
    double sum_t = 0.0, sum_logE = 0.0, sum_t2 = 0.0, sum_t_logE = 0.0;
    
    for (int i = energy_history_.size() - n_points; i < energy_history_.size(); ++i) {
        double t = time_history_[i];
        double logE = std::log(std::max(energy_history_[i], 1e-30));
        
        sum_t += t;
        sum_logE += logE;
        sum_t2 += t * t;
        sum_t_logE += t * logE;
    }
    
    double denom = n_points * sum_t2 - sum_t * sum_t;
    if (std::abs(denom) < 1e-12) return 0.0;
    
    return (n_points * sum_t_logE - sum_t * sum_logE) / denom;
}

// =============================================================================
// NUMERICAL METHODS
// =============================================================================
// En nonlinear_evolution.cpp - función computeCFLTimeStep
double NonlinearEvolution::computeCFLTimeStep() const {
    double min_dt = 1e-12; // Más conservador
    
    for (int i = 0; i < params_.Nr; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;
            
            double cs = computeSoundSpeed(state_curr_[P][idx], state_curr_[RHO][idx]);
            double vA = computeAlfvenSpeed(sqrt(state_curr_[BR][idx]*state_curr_[BR][idx] + 
                                              state_curr_[BZ][idx]*state_curr_[BZ][idx] + 
                                              state_curr_[BTHETA][idx]*state_curr_[BTHETA][idx]), 
                                         state_curr_[RHO][idx]);
            
            double v_max = std::max(cs, vA);
            v_max = std::max(v_max, std::abs(state_curr_[VR][idx]));
            v_max = std::max(v_max, std::abs(state_curr_[VZ][idx]));
            v_max = std::max(v_max, std::abs(state_curr_[VTHETA][idx]));
            
            double dr = (i < params_.Nr - 1) ? (r_grid_[i+1] - r_grid_[i]) : 1e-6;
            double dz = params_.L / params_.Nz;
            
            // CFL más conservador para MHD
            double dt_local = 0.1 * std::min(dr, dz) / (v_max + 1e-12); // 0.1 en lugar de 0.5
            min_dt = std::min(min_dt, dt_local);
        }
    }
    
    return std::max(min_dt, 1e-12);
}

double NonlinearEvolution::derivativeR(const std::vector<double>& f, int i, int j) const {
    // Central difference in radial direction
    if (i > 0 && i < params_.Nr - 1) {
        int idx_prev = (i-1) * params_.Nz + j;
        int idx_next = (i+1) * params_.Nz + j;
        double dr = r_grid_[i+1] - r_grid_[i];
        return (f[idx_next] - f[idx_prev]) / (2.0 * dr);
    }
    return 0.0;
}

double NonlinearEvolution::derivativeZ(const std::vector<double>& f, int i, int j) const {
    // Central difference in axial direction
    if (j > 0 && j < params_.Nz - 1) {
        int idx_prev = i * params_.Nz + (j-1);
        int idx_next = i * params_.Nz + (j+1);
        double dz = z_grid_[j+1] - z_grid_[j];
        return (f[idx_next] - f[idx_prev]) / (2.0 * dz);
    }
    return 0.0;
}

// =============================================================================
// BOUNDARY CONDITIONS
// =============================================================================
void NonlinearEvolution::applyBoundaryConditions() {
    // Radial boundary conditions
    for (int j = 0; j < params_.Nz; ++j) {
        // Axis (r=0) - symmetry conditions
        state_curr_[VR][0*params_.Nz + j] = 0.0; // No flow through axis
        state_curr_[VTHETA][0*params_.Nz + j] = 0.0; // No swirl at axis
        state_curr_[BR][0*params_.Nz + j] = 0.0; // No radial field at axis
        
        // Outer boundary (r=R0) - conducting wall
        int outer_idx = (params_.Nr-1) * params_.Nz + j;
        state_curr_[VR][outer_idx] = 0.0; // No flow through wall
        state_curr_[BR][outer_idx] = 0.0; // No radial field penetration
    }
    
    // Axial boundary conditions (periodic)
    for (int i = 0; i < params_.Nr; ++i) {
        state_curr_[VZ][i*params_.Nz + 0] = state_curr_[VZ][i*params_.Nz + (params_.Nz-1)];
        state_curr_[BZ][i*params_.Nz + 0] = state_curr_[BZ][i*params_.Nz + (params_.Nz-1)];
    }
}

// =============================================================================
// OUTPUT METHODS
// =============================================================================
void NonlinearEvolution::saveState(const std::string& filename) const {
    std::ofstream file(filename);
    file << "r,z,rho,vr,vz,vtheta,p,Br,Bz,Btheta\n";
    file << std::scientific << std::setprecision(6);
    
    for (int i = 0; i < params_.Nr; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;
            file << r_grid_[i] << "," << z_grid_[j] << ","
                 << state_curr_[RHO][idx] << ","
                 << state_curr_[VR][idx] << ","
                 << state_curr_[VZ][idx] << ","
                 << state_curr_[VTHETA][idx] << ","
                 << state_curr_[P][idx] << ","
                 << state_curr_[BR][idx] << ","
                 << state_curr_[BZ][idx] << ","
                 << state_curr_[BTHETA][idx] << "\n";
        }
    }
    
    file.close();
}

void NonlinearEvolution::saveTimeHistory(const std::string& filename) const {
    std::ofstream file(filename);
    file << "time,energy,growth_rate\n";
    file << std::scientific << std::setprecision(6);
    
    for (size_t i = 0; i < time_history_.size(); ++i) {
        file << time_history_[i] << "," << energy_history_[i];
        if (i < growth_rate_history_.size()) {
            file << "," << growth_rate_history_[i];
        } else {
            file << ",0.0";
        }
        file << "\n";
    }
    
    file.close();
}

// Agregar al final de nonlinear_evolution.cpp
void NonlinearEvolution::addArtificialDissipation(double dt) {
    // Coeficientes de disipación artificial
    double nu_viscosity = 1e-3;
    double eta_resistivity = 1e-3;
    
    for (int var = 0; var < 8; ++var) {
        for (int i = 1; i < params_.Nr - 1; ++i) {
            for (int j = 1; j < params_.Nz - 1; ++j) {
                int idx = i * params_.Nz + j;
                
                // Disipación tipo Laplaciano
                double laplacian = 0.0;
                if (var >= VR && var <= VTHETA) {
                    // Disipación para velocidades
                    laplacian = (state_curr_[var][(i+1)*params_.Nz + j] - 2.0 * state_curr_[var][idx] + 
                                state_curr_[var][(i-1)*params_.Nz + j]) / pow(r_grid_[1] - r_grid_[0], 2) +
                                (state_curr_[var][i*params_.Nz + (j+1)] - 2.0 * state_curr_[var][idx] + 
                                state_curr_[var][i*params_.Nz + (j-1)]) / pow(params_.L/params_.Nz, 2);
                    state_curr_[var][idx] += dt * nu_viscosity * laplacian;
                }
                else if (var >= BR && var <= BTHETA) {
                    // Disipación para campos magnéticos (resistividad artificial)
                    laplacian = (state_curr_[var][(i+1)*params_.Nz + j] - 2.0 * state_curr_[var][idx] + 
                                state_curr_[var][(i-1)*params_.Nz + j]) / pow(r_grid_[1] - r_grid_[0], 2) +
                                (state_curr_[var][i*params_.Nz + (j+1)] - 2.0 * state_curr_[var][idx] + 
                                state_curr_[var][i*params_.Nz + (j-1)]) / pow(params_.L/params_.Nz, 2);
                    state_curr_[var][idx] += dt * eta_resistivity * laplacian;
                }
            }
        }
    }
}
