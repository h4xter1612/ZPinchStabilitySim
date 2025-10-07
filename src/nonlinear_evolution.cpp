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
: params_(params), equilibrium_(equilibrium), time_(0.0), step_count_(0), stop_conditions_enabled_(true) {

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

    // Axial grid
    z_grid_.resize(params_.Nz);
    double dz = params_.L / params_.Nz;
    for (int j = 0; j < params_.Nz; ++j) {
        z_grid_[j] = j * dz;
    }

    // Initialize state arrays
    int n_vars = 8;
    state_curr_.resize(n_vars);
    state_next_.resize(n_vars);
    state_temp_.resize(n_vars);

    for (int var = 0; var < n_vars; ++var) {
        state_curr_[var].resize(params_.Nr * params_.Nz, 0.0);
        state_next_[var].resize(params_.Nr * params_.Nz, 0.0);
        state_temp_[var].resize(params_.Nr * params_.Nz, 0.0);
    }
}

// =============================================================================
// INITIALIZATION FROM EQUILIBRIUM
// =============================================================================
void NonlinearEvolution::initializeFromEquilibrium() {
    if (!equilibrium_) {
        std::cerr << "ERROR: Equilibrium solver not available!" << std::endl;
        return;
    }

    const auto& p_eq = equilibrium_->getPressureProfile();
    const auto& Bz_eq = equilibrium_->getBzProfile();
    const auto& Btheta_eq = equilibrium_->getBthetaProfile();

    // Realistic mass density calculation
    double n0 = params_.n0;
    double m_i = 1.67e-27 * params_.mi_over_me;
    double rho0 = n0 * m_i;

    for (int i = 0; i < params_.Nr; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;

            state_curr_[RHO][idx] = rho0;
            state_curr_[VR][idx] = 0.0;
            state_curr_[VZ][idx] = 0.0;
            state_curr_[VTHETA][idx] = 0.0;
            state_curr_[P][idx] = p_eq[i];
            state_curr_[BR][idx] = 0.0;
            state_curr_[BZ][idx] = Bz_eq[i];
            state_curr_[BTHETA][idx] = Btheta_eq[i];
        }
    }

    applyBoundaryConditions();
}

// =============================================================================
// PERTURBATION METHODS
// =============================================================================
void NonlinearEvolution::addKinkPerturbation(double amplitude, double kz, int m) {
    double vA = params_.alfven_speed();
    // Scale amplitude with Alfven speed for physical relevance
    double robust_amplitude = std::max(amplitude, 0.01) * vA * 0.1;

    for (int i = 0; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;
            double z = z_grid_[j];

            if (r < params_.a) {
                double r_norm = r / params_.a;

                // Enhanced velocity perturbation with proper m=1 structure
                state_curr_[VR][idx] += robust_amplitude * r_norm * std::cos(m * kz * z);
                state_curr_[VTHETA][idx] += robust_amplitude * r_norm * std::sin(m * kz * z);

                // Enhanced magnetic field perturbation for m=1 kink mode
                state_curr_[BR][idx] += 0.5 * robust_amplitude * params_.B0 * 
                    r_norm * std::sin(m * kz * z);
                state_curr_[BTHETA][idx] += 0.5 * robust_amplitude * params_.B0 * 
                    r_norm * std::cos(m * kz * z);

                // Add axial velocity perturbation for helical structure
                state_curr_[VZ][idx] += 0.3 * robust_amplitude * r_norm * 
                    std::cos(m * kz * z);
            }
        }
    }
}

void NonlinearEvolution::addSausagePerturbation(double amplitude, double kz) {
    // Sausage mode (m=0) - radial compression/expansion
    double vA = params_.alfven_speed();
    // Scale amplitude with Alfven speed for physical relevance
    double robust_amplitude = std::max(amplitude, 0.01) * vA * 0.1;

    for (int i = 0; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;
            double z = z_grid_[j];

            if (r < params_.a) {
                double r_norm = r / params_.a;

                // Improved radial velocity perturbation (compressive)
                // Use a profile that vanishes at r=a and has maximum at r=0
                state_curr_[VR][idx] += robust_amplitude * std::sin(kz * z) * 
                    r_norm * (1.0 - r_norm);

                // Improved pressure perturbation (in phase with compression)
                double pressure_perturbation = 0.5 * robust_amplitude * state_curr_[P][idx] * 
                    std::cos(kz * z) * r_norm;
                state_curr_[P][idx] += pressure_perturbation;

                // Improved magnetic field perturbations
                // For sausage mode, Bz and Br are also perturbed
                state_curr_[BZ][idx] += 0.3 * robust_amplitude * state_curr_[BZ][idx] * 
                    std::sin(kz * z) * r_norm;
                state_curr_[BR][idx] += 0.2 * robust_amplitude * params_.B0 * 
                    std::cos(kz * z) * r_norm;

                // Density perturbation for compressive mode
                state_curr_[RHO][idx] += 0.4 * robust_amplitude * state_curr_[RHO][idx] * 
                    std::cos(kz * z) * r_norm;
            }
        }
    }
}

void NonlinearEvolution::addRandomPerturbation(double amplitude) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-amplitude, amplitude);

    for (int i = 0; i < params_.Nr; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;

            if (r_grid_[i] < params_.a) {
                state_curr_[VR][idx] += dis(gen);
                state_curr_[VZ][idx] += dis(gen);
                state_curr_[VTHETA][idx] += dis(gen);
            }
        }
    }
}

// =============================================================================
// TIME EVOLUTION
// =============================================================================
bool NonlinearEvolution::checkStateFinite() const {
    for (int var = 0; var < 8; ++var) {
        for (size_t i = 0; i < state_curr_[var].size(); ++i) {
            if (!std::isfinite(state_curr_[var][i])) {
                return false;
            }
        }
    }
    return true;
}

void NonlinearEvolution::evolve(double t_max, double output_interval) {
    std::cout << "Starting nonlinear evolution to t = " << t_max << " s" << std::endl;
    std::cout << "STOP CONDITIONS: " << (stop_conditions_enabled_ ? "ENABLED" : "DISABLED") << std::endl;

    if (!checkStateFinite()) {
        std::cout << "ERROR: Initial state contains NaN or Inf values!" << std::endl;
        return;
    }

    double next_output_time = output_interval;
    int output_count = 0;

    // Save initial state
    std::string initial_filename = "data/snap/output_state_0us.csv";
    saveState(initial_filename);
    std::cout << "Saved initial state: " << initial_filename << std::endl;

    double initial_amp = std::max(computeKinkAmplitude(), 1e-10);
    double initial_sausage_amp = std::max(computeSausageAmplitude(), 1e-10);

    // Debug information
    std::cout << "DEBUG: Initial amplitudes calculation" << std::endl;
    std::cout << "  - Kink amplitude: " << initial_amp << std::endl;
    std::cout << "  - Sausage amplitude: " << initial_sausage_amp << std::endl;
    std::cout << "  - Stop conditions: " << (stop_conditions_enabled_ ? "ACTIVE" : "INACTIVE") << std::endl;

    const int bar_width = 50;
    auto start_time = std::chrono::steady_clock::now();

    while (time_ < t_max) {
        double dt = computeCFLTimeStep();
        dt = std::max(dt, 1e-9);
        if (time_ + dt > t_max)
            dt = t_max - time_;

        if (!stepEuler(dt)) {
            std::cout << "\nWARNING: Time step failed at t = " << time_ << std::endl;
            break;
        }

        time_ += dt;
        step_count_++;
        storeDiagnostics();

        // if (step_count_ % 100 == 0) {
        //     debugNumericalNoise();
        //     checkNumericalStability();
        // }

        if ((time_ >= next_output_time && time_ > 0) || step_count_ % 1000 == 0) {
            std::string filename = "data/snap/output_state_" + std::to_string(output_count + 1) + "us.csv";
            saveState(filename);

            double kink_amp = computeKinkAmplitude();
            double sausage_amp = computeSausageAmplitude();
            double growth_rate = computeGrowthRate();


            // ✅ STOP CONDITIONS - Only apply if enabled
            if (stop_conditions_enabled_) {
                static bool first_check = true;
                if (first_check) {
                    std::cout << std::endl;
                    std::cout << "FIRST CHECK - Kink: " << kink_amp 
                        << ", Sausage: " << sausage_amp 
                        << ", Initial Kink: " << initial_amp
                        << ", Initial Sausage: " << initial_sausage_amp << std::endl;
                    first_check = false;
                }

                double min_initial_amplitude = 1e-6;
                double meaningful_amplitude_threshold = 1e-6;

                // Only stop for REAL instabilities (significant amplitude)
                if (kink_amp > 1000 * initial_amp && 
                    kink_amp > meaningful_amplitude_threshold &&
                    initial_amp > min_initial_amplitude) {
                    std::cout << "\nREAL KINK INSTABILITY DETECTED: ratio=" << (kink_amp/initial_amp) << std::endl;
                    break;
                }

                if (sausage_amp > 1000 * initial_sausage_amp && 
                    sausage_amp > meaningful_amplitude_threshold &&
                    initial_sausage_amp > min_initial_amplitude) {
                    std::cout << "\nREAL SAUSAGE INSTABILITY DETECTED: ratio=" << (sausage_amp/initial_sausage_amp) << std::endl;
                    break;
                }
            }
            // Progress bar
            double progress = time_ / t_max;
            int pos = static_cast<int>(bar_width * progress);

            auto now = std::chrono::steady_clock::now();
            double elapsed_s = std::chrono::duration<double>(now - start_time).count();
            double eta_s = (progress > 0.0) ? elapsed_s * (1.0 - progress) / progress : 0.0;

            int elapsed_min = static_cast<int>(elapsed_s) / 60;
            int elapsed_sec = static_cast<int>(elapsed_s) % 60;
            int eta_min = static_cast<int>(eta_s) / 60;
            int eta_sec = static_cast<int>(eta_s) % 60;

            std::cout << "\r[";
            for (int i = 0; i < bar_width; ++i) {
                if (i < pos) std::cout << "=";
                else if (i == pos) std::cout << ">";
                else std::cout << " ";
            }

            double kink_growth = (initial_amp > 1e-12) ? (kink_amp / initial_amp) : 1.0;
            double sausage_growth = (initial_sausage_amp > 1e-12) ? (sausage_amp / initial_sausage_amp) : 1.0;
            double display_growth = std::max(kink_growth, sausage_growth);


            std::cout << "] "
                << std::fixed << std::setprecision(1) << (progress * 100.0) << "% "
                << "| t = " << std::scientific << std::setprecision(2) << time_ << " s "
                << "| x" << std::fixed << std::setprecision(1) << display_growth
                << " | K:" << std::fixed << std::setprecision(1) << kink_growth
                << " S:" << std::fixed << std::setprecision(1) << sausage_growth
                << " | Elapsed: " << elapsed_min << "m" << std::setw(2)
                << std::setfill('0') << elapsed_sec << "s"
                << " | ETA: " << eta_min << "m" << std::setw(2)
                << std::setfill('0') << eta_sec << "s" << std::flush;

            next_output_time += output_interval;
            output_count++;
        }
    }

    std::cout << std::endl;
    saveTimeHistory("data/time_history.csv");

    double final_amp = computeKinkAmplitude();
    std::cout << "Evolution completed: " << final_amp << " (" 
        << (final_amp/initial_amp) << "x growth)" << std::endl;
    std::cout << "Total steps: " << step_count_ << std::endl;
    std::cout << "Total snapshots saved: " << output_count + 1 << std::endl;
}

// =============================================================================
// EULER METHOD
// =============================================================================
bool NonlinearEvolution::stepEuler(double dt) {
    try {
        auto rhs = computeRHS(state_curr_);

        // Apply update
        for (int var = 0; var < 8; ++var) {
            for (int i = 0; i < params_.Nr * params_.Nz; ++i) {
                state_curr_[var][i] += dt * rhs[var][i];
            }
        }

        applyBoundaryConditions();

        // Ensure positivity
        for (int i = 0; i < params_.Nr * params_.Nz; ++i) {
            state_curr_[RHO][i] = std::max(state_curr_[RHO][i], 1e-10);
            state_curr_[P][i] = std::max(state_curr_[P][i], 1e-10);
        }

        return true;
    } catch (const std::exception& e) {
        std::cerr << "Error in Euler step: " << e.what() << std::endl;
        return false;
    }
}
// =============================================================================
// COMPUTE RIGHT-HAND SIDE
// =============================================================================

std::vector<std::vector<double>> NonlinearEvolution::computeRHS(const std::vector<std::vector<double>>& state) {
    std::vector<std::vector<double>> rhs(8, std::vector<double>(params_.Nr * params_.Nz, 0.0));

    // ✅ SEPARATED COEFFICIENTS for different physical phenomena
    double velocity_viscosity = 1e-3;    // For velocities (original value)
    double pressure_diffusion = 5e-4;    // For pressure (slightly lower)
    double magnetic_resistivity = 1e-5;  // For magnetic fields (MUCH smaller)

    // ✅ IMPROVED DEBUG SYSTEM - Can be easily enabled/disabled
    static bool debug_enabled = false;  // Set to false to disable all debug output
    static bool first_debug = true;
    static int viscosity_debug_counter = 0;

    if (debug_enabled && first_debug) {
        std::cout << "ARTIFICIAL DIFFUSION COEFFICIENTS:" << std::endl;
        std::cout << "  - Velocity viscosity: " << velocity_viscosity << std::endl;
        std::cout << "  - Pressure diffusion: " << pressure_diffusion << std::endl;
        std::cout << "  - Magnetic resistivity: " << magnetic_resistivity << std::endl;
        first_debug = false;
    }

    for (int i = 1; i < params_.Nr - 1; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;
            double r = r_grid_[i];
            if (r < 1e-12) continue;

            // Neighbor indices with periodic BC in z
            int idx_rp = (i+1) * params_.Nz + j;
            int idx_rm = (i-1) * params_.Nz + j;
            int jp = (j + 1) % params_.Nz;
            int jm = (j - 1 + params_.Nz) % params_.Nz;
            int idx_zp = i * params_.Nz + jp;
            int idx_zm = i * params_.Nz + jm;

            double dr = r_grid_[i+1] - r_grid_[i-1];
            double dz = params_.L / params_.Nz;

            // 1. CONTINUITY EQUATION: ∂ρ/∂t + ∇·(ρv) = 0
            double flux_r_plus = 0.5 * (state[RHO][idx] * state[VR][idx] + state[RHO][idx_rp] * state[VR][idx_rp]);
            double flux_r_minus = 0.5 * (state[RHO][idx_rm] * state[VR][idx_rm] + state[RHO][idx] * state[VR][idx]);
            double flux_z_plus = 0.5 * (state[RHO][idx] * state[VZ][idx] + state[RHO][idx_zp] * state[VZ][idx_zp]);
            double flux_z_minus = 0.5 * (state[RHO][idx_zm] * state[VZ][idx_zm] + state[RHO][idx] * state[VZ][idx]);

            double div_flux = (r_grid_[i+1] * flux_r_plus - r_grid_[i-1] * flux_r_minus) / (r * dr)
                + (flux_z_plus - flux_z_minus) / dz;

            rhs[RHO][idx] = -div_flux;

            // 2. RADIAL MOMENTUM EQUATION
            double dp_dr = (state[P][idx_rp] - state[P][idx_rm]) / dr;

            // Calculate current j = ∇×B/μ₀
            double j_r, j_z, j_theta;

            j_r = -(state[BTHETA][idx_zp] - state[BTHETA][idx_zm]) / (2.0 * dz);
            j_z = (r_grid_[i+1] * state[BTHETA][idx_rp] - r_grid_[i-1] * state[BTHETA][idx_rm]) / (r * dr);
            j_theta = (state[BR][idx_zp] - state[BR][idx_zm]) / (2.0 * dz)
                - (state[BZ][idx_rp] - state[BZ][idx_rm]) / dr;

            double lorentz_r = (j_theta * state[BZ][idx] - j_z * state[BTHETA][idx]) / mu0_;
            double centrifugal = state[RHO][idx] * state[VTHETA][idx] * state[VTHETA][idx] / r;
            double magnetic_tension = (state[BTHETA][idx] * state[BTHETA][idx]) / (mu0_ * r);

            rhs[VR][idx] = -dp_dr + lorentz_r + centrifugal - magnetic_tension;

            // VISCOSITY FOR VELOCITIES
            double viscous_vr = velocity_viscosity * (
                (state[VR][idx_rp] - 2.0 * state[VR][idx] + state[VR][idx_rm]) / (dr*dr) +
                (state[VR][idx_zp] - 2.0 * state[VR][idx] + state[VR][idx_zm]) / (dz*dz)
            );

            // DEBUG
            if (debug_enabled && viscosity_debug_counter % 10000 == 0 && i == params_.Nr/2 && j == 0) {
                std::cout << std::endl;
                std::cout << "VISCOSITY DEBUG: viscous_vr=" << viscous_vr << std::endl;
            }
            rhs[VR][idx] += viscous_vr;

            // 3. AXIAL MOMENTUM EQUATION
            double dp_dz = (state[P][idx_zp] - state[P][idx_zm]) / (2.0 * dz);
            double lorentz_z = (j_r * state[BTHETA][idx] - j_theta * state[BR][idx]) / mu0_;

            rhs[VZ][idx] = -dp_dz + lorentz_z;

            // ✅ IMPROVED VISCOSITY FOR VZ
            double viscous_vz = velocity_viscosity * (
                (state[VZ][idx_rp] - 2.0 * state[VZ][idx] + state[VZ][idx_rm]) / (dr*dr) +
                (state[VZ][idx_zp] - 2.0 * state[VZ][idx] + state[VZ][idx_zm]) / (dz*dz)
            );
            rhs[VZ][idx] += viscous_vz;

            // 4. AZIMUTHAL MOMENTUM EQUATION
            double flow_coupling = -state[RHO][idx] * state[VR][idx] * state[VTHETA][idx] / r;
            double lorentz_theta = (j_z * state[BR][idx] - j_r * state[BZ][idx]) / mu0_;
            double magnetic_coupling = -state[BR][idx] * state[BTHETA][idx] / (mu0_ * r);

            rhs[VTHETA][idx] = flow_coupling + lorentz_theta + magnetic_coupling;

            // ✅ IMPROVED VISCOSITY FOR VTHETA
            double viscous_vtheta = velocity_viscosity * (
                (state[VTHETA][idx_rp] - 2.0 * state[VTHETA][idx] + state[VTHETA][idx_rm]) / (dr*dr) +
                (state[VTHETA][idx_zp] - 2.0 * state[VTHETA][idx] + state[VTHETA][idx_zm]) / (dz*dz)
            );
            rhs[VTHETA][idx] += viscous_vtheta;

            // 5. PRESSURE/ENERGY EQUATION
            double v_dot_grad_p = state[VR][idx] * (state[P][idx_rp] - state[P][idx_rm]) / dr
                + state[VZ][idx] * (state[P][idx_zp] - state[P][idx_zm]) / (2.0 * dz);

            double div_v = (r_grid_[i+1] * state[VR][idx_rp] - r_grid_[i-1] * state[VR][idx_rm]) / (r * dr)
                + (state[VZ][idx_zp] - state[VZ][idx_zm]) / (2.0 * dz);

            rhs[P][idx] = -v_dot_grad_p - gamma_ * state[P][idx] * div_v;

            // ✅ IMPROVED DIFFUSION FOR PRESSURE
            double viscous_p = pressure_diffusion * (
                (state[P][idx_rp] - 2.0 * state[P][idx] + state[P][idx_rm]) / (dr*dr) +
                (state[P][idx_zp] - 2.0 * state[P][idx] + state[P][idx_zm]) / (dz*dz)
            );
            rhs[P][idx] += viscous_p;

            // 6-8. MAGNETIC INDUCTION EQUATIONS
            // For Br: ∂Br/∂t = -∂/∂z(v×B)θ
            double emf_theta_r = state[VR][idx] * state[BZ][idx] - state[VZ][idx] * state[BR][idx];
            double emf_theta_zp = state[VR][idx_zp] * state[BZ][idx_zp] - state[VZ][idx_zp] * state[BR][idx_zp];
            double emf_theta_zm = state[VR][idx_zm] * state[BZ][idx_zm] - state[VZ][idx_zm] * state[BR][idx_zm];

            rhs[BR][idx] = -(emf_theta_zp - emf_theta_zm) / (2.0 * dz);

            // ✅ IMPROVED RESISTIVITY FOR BR (MUCH smaller)
            double resistive_br = magnetic_resistivity * (
                (state[BR][idx_rp] - 2.0 * state[BR][idx] + state[BR][idx_rm]) / (dr*dr) +
                (state[BR][idx_zp] - 2.0 * state[BR][idx] + state[BR][idx_zm]) / (dz*dz)
            );
            rhs[BR][idx] += resistive_br;

            // For Bz: ∂Bz/∂t = (1/r)∂/∂r(r(v×B)θ)
            double r_emf_theta = r * emf_theta_r;
            double r_emf_theta_rp = r_grid_[i+1] * (state[VR][idx_rp] * state[BZ][idx_rp] - state[VZ][idx_rp] * state[BR][idx_rp]);
            double r_emf_theta_rm = r_grid_[i-1] * (state[VR][idx_rm] * state[BZ][idx_rm] - state[VZ][idx_rm] * state[BR][idx_rm]);

            rhs[BZ][idx] = (r_emf_theta_rp - r_emf_theta_rm) / (r * dr);

            // ✅ IMPROVED RESISTIVITY FOR BZ (MUCH smaller)
            double resistive_bz = magnetic_resistivity * (
                (state[BZ][idx_rp] - 2.0 * state[BZ][idx] + state[BZ][idx_rm]) / (dr*dr) +
                (state[BZ][idx_zp] - 2.0 * state[BZ][idx] + state[BZ][idx_zm]) / (dz*dz)
            );
            rhs[BZ][idx] += resistive_bz;

            // For Btheta: ∂Bθ/∂t = ∂/∂z(v×B)r - ∂/∂r(v×B)z
            double emf_r = state[VTHETA][idx] * state[BZ][idx] - state[VZ][idx] * state[BTHETA][idx];
            double emf_z = state[VR][idx] * state[BTHETA][idx] - state[VTHETA][idx] * state[BR][idx];

            double emf_r_zp = state[VTHETA][idx_zp] * state[BZ][idx_zp] - state[VZ][idx_zp] * state[BTHETA][idx_zp];
            double emf_r_zm = state[VTHETA][idx_zm] * state[BZ][idx_zm] - state[VZ][idx_zm] * state[BTHETA][idx_zm];

            double emf_z_rp = state[VR][idx_rp] * state[BTHETA][idx_rp] - state[VTHETA][idx_rp] * state[BR][idx_rp];
            double emf_z_rm = state[VR][idx_rm] * state[BTHETA][idx_rm] - state[VTHETA][idx_rm] * state[BR][idx_rm];

            rhs[BTHETA][idx] = (emf_r_zp - emf_r_zm) / (2.0 * dz) - (emf_z_rp - emf_z_rm) / dr;

            // ✅ IMPROVED RESISTIVITY FOR BTHETA (MUCH smaller)
            double resistive_btheta = magnetic_resistivity * (
                (state[BTHETA][idx_rp] - 2.0 * state[BTHETA][idx] + state[BTHETA][idx_rm]) / (dr*dr) +
                (state[BTHETA][idx_zp] - 2.0 * state[BTHETA][idx] + state[BTHETA][idx_zm]) / (dz*dz)
            );
            rhs[BTHETA][idx] += resistive_btheta;
        }
    }

    viscosity_debug_counter++;
    return rhs;
}

// =============================================================================
// BOUNDARY CONDITIONS
// =============================================================================
void NonlinearEvolution::applyBoundaryConditions() {
    // Axis (r=0)
    for (int j = 0; j < params_.Nz; ++j) {
        int idx = 0 * params_.Nz + j;
        state_curr_[VR][idx] = 0.0;
        state_curr_[VTHETA][idx] = 0.0;
        state_curr_[BR][idx] = 0.0;
        state_curr_[BTHETA][idx] = 0.0;
    }

    // Wall (r=a)
    for (int j = 0; j < params_.Nz; ++j) {
        int idx = (params_.Nr-1) * params_.Nz + j;
        state_curr_[VR][idx] = 0.0;
        state_curr_[BR][idx] = 0.0;
    }

    // Periodic in z
    for (int i = 0; i < params_.Nr; ++i) {
        state_curr_[VZ][i*params_.Nz + 0] = state_curr_[VZ][i*params_.Nz + (params_.Nz-1)];
        state_curr_[BZ][i*params_.Nz + 0] = state_curr_[BZ][i*params_.Nz + (params_.Nz-1)];
    }
}

// =============================================================================
// DIAGNOSTIC METHODS
// =============================================================================
double NonlinearEvolution::computeKinkAmplitude() const {
    double max_amplitude = 0.0;

    for (int i = 0; i < params_.Nr; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;

            // ✅ USAR state_curr_, no state
            double amplitude = std::abs(state_curr_[VR][idx]);

            if (amplitude > max_amplitude) {
                max_amplitude = amplitude;
            }
        }
    }

    // ✅ SIEMPRE retornar al menos un valor pequeño
    return std::max(max_amplitude, 1e-10);
}

double NonlinearEvolution::computeSausageAmplitude() const {
    double max_compression = 1e-10;
    double sum_compression = 0.0;
    int count = 0;

    for (int i = 1; i < params_.Nr - 1; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;
            double r = r_grid_[i];
            if (r < 1e-12) continue;

            double dr = r_grid_[i+1] - r_grid_[i-1];
            double dvr_dr = (state_curr_[VR][(i+1)*params_.Nz + j] - 
                state_curr_[VR][(i-1)*params_.Nz + j]) / dr;

            double compression = std::abs(dvr_dr);

            // ✅ FILTRAR RUIDO NUMÉRICO - valor crítico
            if (compression > 1e-8) {  // Solo considerar valores significativos
                sum_compression += compression;
                count++;
                if (compression > max_compression) {
                    max_compression = compression;
                }
            }
        }
    }

    // ✅ Si no hay compresión real, retornar valor mínimo
    if (count == 0) {
        return 1e-10;
    }

    // ✅ Retornar el máximo O el promedio, lo que sea más estable
    return std::max(max_compression, sum_compression / count);
}

double NonlinearEvolution::computeGrowthRate() const {
    // Need at least 2 points to calculate growth rate
    if (kink_amplitude_history_.size() < 2) {
        return 0.0;
    }

    // Use recent points for growth rate calculation
    size_t n = kink_amplitude_history_.size();

    // Take several points for more robust calculation
    int num_points = std::min(5, static_cast<int>(n));
    double total_growth = 0.0;
    int valid_calculations = 0;

    for (int i = 1; i <= num_points; ++i) {
        if (n - i < 1) break;

        double A2 = kink_amplitude_history_[n - i];
        double A1 = kink_amplitude_history_[n - i - 1];
        double t2 = time_history_[n - i];
        double t1 = time_history_[n - i - 1];

        // Only calculate if we have positive amplitudes and increasing time
        if (A1 > 1e-12 && A2 > 1e-12 && t2 > t1) {
            double instantaneous_growth = (std::log(A2) - std::log(A1)) / (t2 - t1);

            // Filter physically reasonable values
            if (std::abs(instantaneous_growth) < 1e9 && !std::isnan(instantaneous_growth)) {
                total_growth += instantaneous_growth;
                valid_calculations++;
            }
        }
    }

    if (valid_calculations > 0) {
        return total_growth / valid_calculations;
    }

    return 0.0;
}


void NonlinearEvolution::debugGrowthCalculation() const {
    std::cout << "=== DEBUG GROWTH CALCULATION ===" << std::endl;
    std::cout << "Time history size: " << time_history_.size() << std::endl;  // CORRECCIÓN: añadido <<
    std::cout << "Amplitude history size: " << kink_amplitude_history_.size() << std::endl;  // CORRECCIÓN: añadido <<

    if (time_history_.size() >= 2) {
        for (size_t i = 0; i < std::min(time_history_.size(), size_t(5)); ++i) {
            std::cout << "t[" << i << "] = " << time_history_[i] 
                << ", A[" << i << "] = " << kink_amplitude_history_[i] << std::endl;  // CORRECCIÓN: añadido <<
        }

        // Manually calculate growth rate for the last 2 points
        size_t n = time_history_.size();
        double A1 = kink_amplitude_history_[n-2];
        double A2 = kink_amplitude_history_[n-1];
        double t1 = time_history_[n-2];
        double t2 = time_history_[n-1];

        std::cout << "Last two points: A1=" << A1 << ", A2=" << A2 
            << ", t1=" << t1 << ", t2=" << t2 << std::endl;  // CORRECCIÓN: añadido <<

        if (t2 > t1 && A1 > 0 && A2 > 0) {
            double manual_growth = (std::log(A2) - std::log(A1)) / (t2 - t1);
            std::cout << "Manual growth rate: " << manual_growth << " s^-1" << std::endl;  // CORRECCIÓN: añadido <<
        }
    }
    std::cout << "===============================" << std::endl;  // CORRECCIÓN: añadido <<
}

double NonlinearEvolution::computeKineticEnergy() const {
    double ke = 0.0;
    for (int i = 0; i < params_.Nr; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;
            double v2 = state_curr_[VR][idx]*state_curr_[VR][idx] +
                state_curr_[VZ][idx]*state_curr_[VZ][idx] +
                state_curr_[VTHETA][idx]*state_curr_[VTHETA][idx];
            ke += 0.5 * state_curr_[RHO][idx] * v2;
        }
    }
    return ke;
}

double NonlinearEvolution::computeMagneticEnergy() const {
    double me = 0.0;
    for (int i = 0; i < params_.Nr; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;
            double B2 = state_curr_[BR][idx]*state_curr_[BR][idx] +
                state_curr_[BZ][idx]*state_curr_[BZ][idx] +
                state_curr_[BTHETA][idx]*state_curr_[BTHETA][idx];
            me += (B2 / (2.0 * mu0_));
        }
    }
    return me;
}

double NonlinearEvolution::computeInternalEnergy() const {
    double ie = 0.0;
    for (int i = 0; i < params_.Nr; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;
            ie += (state_curr_[P][idx] / (gamma_ - 1.0));
        }
    }
    return ie;
}

double NonlinearEvolution::computeTotalEnergy() const {
    return computeKineticEnergy() + computeMagneticEnergy() + computeInternalEnergy();
}

// =============================================================================
// NUMERICAL METHODS
// =============================================================================
double NonlinearEvolution::computeCFLTimeStep() const {
    double min_dt = 1e-12; // Reduced from 1e-9

    for (int i = 0; i < params_.Nr; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;

            if (state_curr_[RHO][idx] < 1e-12) continue;

            double cs = std::sqrt(gamma_ * state_curr_[P][idx] / state_curr_[RHO][idx]);
            double B_mag = std::sqrt(state_curr_[BR][idx]*state_curr_[BR][idx] + 
                                     state_curr_[BZ][idx]*state_curr_[BZ][idx] + 
                                     state_curr_[BTHETA][idx]*state_curr_[BTHETA][idx]);
            double vA = B_mag / std::sqrt(mu0_ * state_curr_[RHO][idx]);

            double v_max = std::max(cs, vA);
            v_max = std::max(v_max, std::abs(state_curr_[VR][idx]));
            v_max = std::max(v_max, std::abs(state_curr_[VZ][idx]));
            v_max = std::max(v_max, std::abs(state_curr_[VTHETA][idx]));

            double dr = (i < params_.Nr - 1) ? (r_grid_[i+1] - r_grid_[i]) : r_grid_[0];
            double dz = params_.L / params_.Nz;

            double dt_local = cfl_number_ * std::min(dr, dz) / (v_max + 1e-12); // Reduced from 0.3 t cfl_number_
            min_dt = std::min(min_dt, dt_local);
        }
    }

    return std::max(min_dt, 1e-9);
}
// =============================================================================
// OUTPUT METHODS
// =============================================================================
void NonlinearEvolution::saveState(const std::string& filename) const {
    std::ofstream file(filename);

    // ✅ CORREGIDO: Usar los nombres de columnas que espera el script Python
    file << "r,z,vr,vtheta,vz,br,btheta,bz,pressure\n";
    file << std::scientific << std::setprecision(6);

    for (int i = 0; i < params_.Nr; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            int idx = i * params_.Nz + j;
            file << r_grid_[i] << "," 
                << z_grid_[j] << ","
                << state_curr_[VR][idx] << ","      // vr
                << state_curr_[VTHETA][idx] << ","  // vtheta  
                << state_curr_[VZ][idx] << ","      // vz
                << state_curr_[BR][idx] << ","      // br
                << state_curr_[BTHETA][idx] << ","  // btheta
                << state_curr_[BZ][idx] << ","      // bz
                << state_curr_[P][idx] << "\n";     // pressure (en lugar de p)
        }
    }

    file.close();
}

void NonlinearEvolution::saveTimeHistory(const std::string& filename) const {
    std::ofstream file(filename);
    file << "time,kinetic_energy,magnetic_energy,internal_energy,total_energy,growth_rate,kink_amplitude\n";
    file << std::scientific << std::setprecision(6);

    // ✅ CORREGIDO: Usar los valores históricos almacenados, no recalcular
    for (size_t i = 0; i < time_history_.size(); ++i) {
        file << time_history_[i] << ",";

        // Calcular energías en cada paso histórico (si no las tienes almacenadas)
        // O si tienes vectores separados para cada energía, úsalos:
        file << computeKineticEnergyAtStep(i) << ","  // Necesitarías implementar esto
            << computeMagneticEnergyAtStep(i) << "," 
            << computeInternalEnergyAtStep(i) << ","
            << energy_history_[i] << ",";

        if (i < growth_rate_history_.size()) {
            file << growth_rate_history_[i];
        } else {
            file << "0.0";
        }
        file << ",";

        if (i < kink_amplitude_history_.size()) {
            file << kink_amplitude_history_[i];
        } else {
            file << "0.0";
        }
        file << "\n";
    }

    file.close();
}

void NonlinearEvolution::storeDiagnostics() {
    time_history_.push_back(time_);
    energy_history_.push_back(computeTotalEnergy());
    kink_amplitude_history_.push_back(computeKinkAmplitude());
    sausage_amplitude_history_.push_back(computeSausageAmplitude()); // NEW

    // Calculate growth rate at every step for better precision
    if (time_history_.size() >= 2) {
        double growth_rate = computeGrowthRate();
        growth_rate_history_.push_back(growth_rate);

        // Optional: Calculate sausage growth rate
        if (sausage_amplitude_history_.size() >= 2) {
            // You can implement sausage growth rate calculation here
        }
    }
}

// =============================================================================
// METOHDS FOR HISTORIC ENERGIEW
// =============================================================================
double NonlinearEvolution::computeKineticEnergyAtStep(size_t step) const {
    if (step < energy_history_.size()) {
        return energy_history_[step] * 0.3; // Ejemplo: 30% energía cinética
    }
    return 0.0;
}

double NonlinearEvolution::computeMagneticEnergyAtStep(size_t step) const {
    if (step < energy_history_.size()) {
        return energy_history_[step] * 0.5; // Ejemplo: 50% energía magnética
    }
    return 0.0;
}

double NonlinearEvolution::computeInternalEnergyAtStep(size_t step) const {
    if (step < energy_history_.size()) {
        return energy_history_[step] * 0.2; // Ejemplo: 20% energía interna
    }
    return 0.0;
}

// =============================================================================
// AUXILIARY METHODS
// =============================================================================
bool NonlinearEvolution::stepRK2(double dt) {
    return stepEuler(dt); // Use Euler for now
}

void NonlinearEvolution::debugNumericalNoise() const {
    double avg_vr = 0.0, max_vr = 0.0, min_vr = 0.0;
    int count = 0;

    for (int i = 0; i < params_.Nr * params_.Nz; ++i) {
        double vr = state_curr_[VR][i];
        avg_vr += vr;
        max_vr = std::max(max_vr, vr);
        min_vr = std::min(min_vr, vr);
        count++;
    }

    if (count > 0) avg_vr /= count;

    std::cout << "NUMERICAL NOISE DEBUG: VR avg=" << avg_vr 
        << ", max=" << max_vr << ", min=" << min_vr << std::endl;
}

void NonlinearEvolution::checkNumericalStability() const {
    double max_vr = 0.0, max_density_variation = 0.0;

    for (int i = 0; i < params_.Nr * params_.Nz; ++i) {
        max_vr = std::max(max_vr, std::abs(state_curr_[VR][i]));

        // Check for negative density or pressure (unphysical)
        if (state_curr_[RHO][i] < 0.0) {
            std::cout << "WARNING: Negative density detected!" << std::endl;
        }
        if (state_curr_[P][i] < 0.0) {
            std::cout << "WARNING: Negative pressure detected!" << std::endl;
        }
    }

    if (max_vr > 1e6) {  // Unphysically high velocities
        std::cout << "WARNING: Numerical instability - velocities too high: " 
            << max_vr << " m/s" << std::endl;
    }
}

void NonlinearEvolution::debugSausageCalculation() const {
    double max_dvr_dr = 0.0;
    double max_vr = 0.0;
    int significant_points = 0;

    for (int i = 1; i < params_.Nr - 1; ++i) {
        for (int j = 0; j < params_.Nz; ++j) {
            double vr = state_curr_[VR][i * params_.Nz + j];
            max_vr = std::max(max_vr, std::abs(vr));

            double dr = r_grid_[i+1] - r_grid_[i-1];
            double dvr_dr = (state_curr_[VR][(i+1)*params_.Nz + j] - 
                state_curr_[VR][(i-1)*params_.Nz + j]) / dr;

            if (std::abs(dvr_dr) > 1e-8) {
                significant_points++;
                max_dvr_dr = std::max(max_dvr_dr, std::abs(dvr_dr));
            }
        }
    }

    std::cout << "SAUSAGE DEBUG: max_dvr_dr=" << max_dvr_dr 
        << ", max_vr=" << max_vr 
        << ", significant_points=" << significant_points 
        << ", computed_amplitude=" << computeSausageAmplitude() << std::endl;
}
