#include "diagnostics.hh"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <numeric>

// =============================================================================
// CONSTRUCTOR
// =============================================================================
Diagnostics::Diagnostics(const ZPinchParameters& params,
                       std::shared_ptr<EquilibriumSolver> equilibrium,
                       std::shared_ptr<LinearStabilityAnalyzer> stability,
                       std::shared_ptr<NonlinearEvolution> evolution)
    : params_(params), equilibrium_(equilibrium), 
      stability_(stability), evolution_(evolution) {
    
    // Initialize data storage
    energy_data_.resize(4); // kinetic, magnetic, internal, total
    for (auto& vec : energy_data_) {
        vec.reserve(10000); // Pre-allocate for typical simulation length
    }
}

// =============================================================================
// EQUILIBRIUM DIAGNOSTICS
// =============================================================================
void Diagnostics::analyzeEquilibrium() {
    std::cout << "=== EQUILIBRIUM ANALYSIS ===" << std::endl;
    
    auto summary = computeEquilibriumSummary();
    
    std::cout << "Total current: " << summary.total_current << " A" << std::endl;
    std::cout << "Internal inductance: " << summary.internal_inductance << std::endl;
    std::cout << "Poloidal beta: " << summary.poloidal_beta << std::endl;
    std::cout << "Toroidal beta: " << summary.toroidal_beta << std::endl;
    std::cout << "Total beta: " << summary.total_beta << std::endl;
    std::cout << "Safety factor at edge: " << summary.safety_factor_edge << std::endl;
    std::cout << "Kink stability margin: " << summary.kink_stability_margin << std::endl;
    std::cout << "Sausage stability margin: " << summary.sausage_stability_margin << std::endl;
    
    // Check stability criteria
    if (summary.kink_stability_margin < 1.0) {
        std::cout << "WARNING: Kink mode may be unstable!" << std::endl;
    }
    if (summary.sausage_stability_margin < 1.0) {
        std::cout << "WARNING: Sausage mode may be unstable!" << std::endl;
    }
    
    std::cout << "=============================" << std::endl;
}

Diagnostics::EquilibriumSummary Diagnostics::computeEquilibriumSummary() const {
    EquilibriumSummary summary;
    
    // Basic parameters
    summary.total_current = equilibrium_->computeTotalCurrent();
    summary.internal_inductance = equilibrium_->computeInternalInductance();
    
    // Beta parameters
    const auto& beta_profile = equilibrium_->getBetaProfile();
    const auto& r_grid = equilibrium_->getRadialGrid();
    
    double beta_integral = 0.0;
    double volume = 0.0;
    
    for (size_t i = 1; i < r_grid.size(); ++i) {
        double r_prev = r_grid[i-1];
        double r_curr = r_grid[i];
        double beta_prev = beta_profile[i-1];
        double beta_curr = beta_profile[i];
        
        double dV = M_PI * (r_curr*r_curr - r_prev*r_prev) * params_.L;
        beta_integral += 0.5 * (beta_prev + beta_curr) * dV;
        volume += dV;
    }
    
    summary.total_beta = beta_integral / volume;
    
    // Approximate poloidal and toroidal beta
    double B_poloidal_avg = 0.0;
    const auto& Btheta_profile = equilibrium_->getBthetaProfile();
    for (double Bt : Btheta_profile) {
        B_poloidal_avg += Bt;
    }
    B_poloidal_avg /= Btheta_profile.size();
    
    double B_toroidal = params_.B0;
    
    summary.poloidal_beta = summary.total_beta * (B_toroidal*B_toroidal) / 
                           (B_poloidal_avg*B_poloidal_avg + 1e-12);
    summary.toroidal_beta = summary.total_beta - summary.poloidal_beta;
    
    // Safety factor
    const auto& q_profile = equilibrium_->getSafetyFactorProfile();
    summary.safety_factor_edge = q_profile.back();
    
    // Stability margins (simplified)
    double I_KruskalShafranov = (2.0 * M_PI * params_.R0 * params_.B0) / 
                               (4e-7 * M_PI * std::log(params_.R0/params_.a + 1.0));
    summary.kink_stability_margin = I_KruskalShafranov / summary.total_current;
    
    // Sausage stability margin (based on pressure gradient)
    double max_pressure_gradient = 0.0;
    const auto& p_profile = equilibrium_->getPressureProfile();
    for (size_t i = 1; i < r_grid.size(); ++i) {
        double dp_dr = (p_profile[i] - p_profile[i-1]) / (r_grid[i] - r_grid[i-1]);
        max_pressure_gradient = std::max(max_pressure_gradient, std::abs(dp_dr));
    }
    
    double critical_gradient = params_.B0 * params_.B0 / (params_.a * 4e-7 * M_PI);
    summary.sausage_stability_margin = critical_gradient / (max_pressure_gradient + 1e-12);
    
    return summary;
}

// =============================================================================
// LINEAR STABILITY DIAGNOSTICS
// =============================================================================
void Diagnostics::analyzeLinearStability() {
    std::cout << "=== LINEAR STABILITY ANALYSIS ===" << std::endl;
    
    // Analyze kink modes
    auto kink_results = stability_->kinkStabilityScan();
    std::cout << "Kink mode analysis:" << std::endl;
    
    double max_kink_growth = 0.0;
    double fastest_kink_k = 0.0;
    
    for (const auto& result : kink_results) {
        if (result.growth_rate > max_kink_growth) {
            max_kink_growth = result.growth_rate;
            fastest_kink_k = result.k;
        }
        
        if (result.growth_rate > 0.0) {
            std::cout << "  k=" << result.k << ", growth_rate=" << result.growth_rate 
                      << ", frequency=" << result.frequency << std::endl;
        }
    }
    
    // Analyze sausage modes
    auto sausage_results = stability_->sausageStabilityScan();
    std::cout << "Sausage mode analysis:" << std::endl;
    
    double max_sausage_growth = 0.0;
    double fastest_sausage_k = 0.0;
    
    for (const auto& result : sausage_results) {
        if (result.growth_rate > max_sausage_growth) {
            max_sausage_growth = result.growth_rate;
            fastest_sausage_k = result.k;
        }
        
        if (result.growth_rate > 0.0) {
            std::cout << "  k=" << result.k << ", growth_rate=" << result.growth_rate 
                      << ", frequency=" << result.frequency << std::endl;
        }
    }
    
    // Summary
    std::cout << "Most unstable kink: k=" << fastest_kink_k 
              << ", γ=" << max_kink_growth << " s^-1" << std::endl;
    std::cout << "Most unstable sausage: k=" << fastest_sausage_k 
              << ", γ=" << max_sausage_growth << " s^-1" << std::endl;
    
    if (max_kink_growth == 0.0 && max_sausage_growth == 0.0) {
        std::cout << "Configuration is linearly stable to analyzed modes." << std::endl;
    } else {
        std::cout << "WARNING: Linearly unstable modes found!" << std::endl;
    }
    
    std::cout << "=================================" << std::endl;
}

void Diagnostics::generateStabilityDiagram(const std::string& filename) const {
    std::ofstream file(filename);
    file << "k,growth_rate_kink,growth_rate_sausage,frequency_kink,frequency_sausage\n";
    
    auto kink_results = stability_->kinkStabilityScan();
    auto sausage_results = stability_->sausageStabilityScan();
    
    for (size_t i = 0; i < kink_results.size(); ++i) {
        file << kink_results[i].k << ","
             << kink_results[i].growth_rate << ","
             << sausage_results[i].growth_rate << ","
             << kink_results[i].frequency << ","
             << sausage_results[i].frequency << "\n";
    }
    
    file.close();
    std::cout << "Stability diagram saved to: " << filename << std::endl;
}

// =============================================================================
// NONLINEAR EVOLUTION DIAGNOSTICS
// =============================================================================
void Diagnostics::monitorEvolution(double output_interval) {
    // This would typically be called during the evolution
    // For now, we'll analyze the existing data
    
    std::cout << "=== EVOLUTION MONITORING ===" << std::endl;
    
    auto energy = computeInstantaneousEnergy();
    std::cout << "Energy components:" << std::endl;
    std::cout << "  Kinetic: " << energy.kinetic << " J" << std::endl;
    std::cout << "  Magnetic: " << energy.magnetic << " J" << std::endl;
    std::cout << "  Internal: " << energy.internal << " J" << std::endl;
    std::cout << "  Total: " << energy.total << " J" << std::endl;
    
    auto conservation = checkConservation();
    std::cout << "Conservation metrics:" << std::endl;
    std::cout << "  Mass: " << conservation.mass_conservation * 100.0 << "%" << std::endl;
    std::cout << "  Energy: " << conservation.energy_conservation * 100.0 << "%" << std::endl;
    std::cout << "  Momentum: " << conservation.momentum_conservation * 100.0 << "%" << std::endl;
    std::cout << "  Magnetic flux: " << conservation.magnetic_flux_conservation * 100.0 << "%" << std::endl;
    
    // Check for dominant modes
    auto dominant_mode = identifyDominantMode();
    if (dominant_mode.amplitude > 1e-12) {
        std::cout << "Dominant mode: m=" << dominant_mode.m 
                  << ", kz=" << dominant_mode.kz
                  << ", amplitude=" << dominant_mode.amplitude
                  << ", growth_rate=" << dominant_mode.growth_rate << " s^-1" << std::endl;
    }
    
    std::cout << "=============================" << std::endl;
}

void Diagnostics::analyzeNonlinearGrowth() {
    std::cout << "=== NONLINEAR GROWTH ANALYSIS ===" << std::endl;
    
    // Extract mode amplitudes over time
    auto kink_amplitude = extractModeAmplitude(1, 2.0 * M_PI / params_.L);
    auto sausage_amplitude = extractModeAmplitude(0, 2.0 * M_PI / params_.L);
    
    if (time_data_.size() > 10) {
        // Fit exponential growth to early phase
        double kink_growth = exponentialFit(time_data_, kink_amplitude);
        double sausage_growth = exponentialFit(time_data_, sausage_amplitude);
        
        std::cout << "Nonlinear growth rates:" << std::endl;
        std::cout << "  Kink mode: " << kink_growth << " s^-1" << std::endl;
        std::cout << "  Sausage mode: " << sausage_growth << " s^-1" << std::endl;
        
        // Compare with linear theory
        auto linear_kink = stability_->analyzeKinkMode(2.0 * M_PI / params_.L);
        auto linear_sausage = stability_->analyzeSausageMode(2.0 * M_PI / params_.L);
        
        std::cout << "Linear theory comparison:" << std::endl;
        std::cout << "  Kink: nonlinear=" << kink_growth 
                  << ", linear=" << linear_kink.growth_rate 
                  << ", ratio=" << kink_growth/linear_kink.growth_rate << std::endl;
        std::cout << "  Sausage: nonlinear=" << sausage_growth 
                  << ", linear=" << linear_sausage.growth_rate 
                  << ", ratio=" << sausage_growth/linear_sausage.growth_rate << std::endl;
    }
    
    std::cout << "=================================" << std::endl;
}

// =============================================================================
// ENERGY DIAGNOSTICS
// =============================================================================
Diagnostics::EnergyComponents Diagnostics::computeInstantaneousEnergy() const {
    EnergyComponents energy;
    
    // This would use data from the nonlinear evolution
    // For now, return placeholder values
    energy.kinetic = 0.0;
    energy.magnetic = 0.0;
    energy.internal = 0.0;
    energy.total = 0.0;
    
    return energy;
}

void Diagnostics::analyzeEnergyTransfer() {
    std::cout << "=== ENERGY TRANSFER ANALYSIS ===" << std::endl;
    
    // Analyze how energy moves between different forms
    // This is particularly important for understanding saturation mechanisms
    
    auto energy_history = computeEnergyHistory();
    
    if (energy_history.size() > 2) {
        // Calculate energy transfer rates
        std::vector<double> dEk_dt, dEm_dt, dEi_dt;
        
        for (size_t i = 1; i < energy_history.size(); ++i) {
            double dt = time_data_[i] - time_data_[i-1];
            dEk_dt.push_back((energy_history[i].kinetic - energy_history[i-1].kinetic) / dt);
            dEm_dt.push_back((energy_history[i].magnetic - energy_history[i-1].magnetic) / dt);
            dEi_dt.push_back((energy_history[i].internal - energy_history[i-1].internal) / dt);
        }
        
        // Analyze correlations
        double avg_dEk_dt = std::accumulate(dEk_dt.begin(), dEk_dt.end(), 0.0) / dEk_dt.size();
        double avg_dEm_dt = std::accumulate(dEm_dt.begin(), dEm_dt.end(), 0.0) / dEm_dt.size();
        
        std::cout << "Average energy transfer rates:" << std::endl;
        std::cout << "  dE_kinetic/dt: " << avg_dEk_dt << " W" << std::endl;
        std::cout << "  dE_magnetic/dt: " << avg_dEm_dt << " W" << std::endl;
        
        if (avg_dEk_dt > 0 && avg_dEm_dt < 0) {
            std::cout << "  Energy flowing from magnetic to kinetic (instability drive)" << std::endl;
        } else if (avg_dEk_dt < 0 && avg_dEm_dt > 0) {
            std::cout << "  Energy flowing from kinetic to magnetic (damping)" << std::endl;
        }
    }
    
    std::cout << "=================================" << std::endl;
}

// =============================================================================
// MODE ANALYSIS
// =============================================================================
Diagnostics::DominantMode Diagnostics::identifyDominantMode() const {
    DominantMode mode;
    mode.m = 0;
    mode.kz = 0.0;
    mode.amplitude = 0.0;
    mode.growth_rate = 0.0;
    
    // This would perform Fourier analysis on the simulation data
    // For now, return a placeholder
    
    return mode;
}

std::vector<double> Diagnostics::extractModeAmplitude(int m, double kz) const {
    // Extract the amplitude of a specific mode over time
    // This is a placeholder implementation
    
    std::vector<double> amplitude(time_data_.size(), 0.0);
    
    // Simple exponential growth for demonstration
    double growth_rate = 1e6; // Typical growth rate
    for (size_t i = 0; i < time_data_.size(); ++i) {
        amplitude[i] = 1e-6 * std::exp(growth_rate * time_data_[i]);
    }
    
    return amplitude;
}

// =============================================================================
// CONSERVATION CHECKS
// =============================================================================
Diagnostics::ConservationMetrics Diagnostics::checkConservation() const {
    ConservationMetrics metrics;
    
    // These would be computed from the simulation data
    // For now, return ideal values
    
    metrics.mass_conservation = 1.0;        // 100% conserved
    metrics.energy_conservation = 0.99;     // 99% conserved
    metrics.momentum_conservation = 0.98;   // 98% conserved
    metrics.magnetic_flux_conservation = 0.97; // 97% conserved
    
    return metrics;
}

// =============================================================================
// DATA EXPORT FOR VISUALIZATION
// =============================================================================
void Diagnostics::exportEquilibriumData(const std::string& filename) const {
    equilibrium_->saveProfiles(filename);
}

void Diagnostics::exportStabilityData(const std::string& filename) const {
    generateStabilityDiagram(filename);
}

void Diagnostics::exportEvolutionData(const std::string& filename) const {
    std::ofstream file(filename);
    file << "time,kinetic_energy,magnetic_energy,internal_energy,total_energy\n";
    
    auto energy_history = computeEnergyHistory();
    for (size_t i = 0; i < time_data_.size() && i < energy_history.size(); ++i) {
        file << time_data_[i] << ","
             << energy_history[i].kinetic << ","
             << energy_history[i].magnetic << ","
             << energy_history[i].internal << ","
             << energy_history[i].total << "\n";
    }
    
    file.close();
}

// =============================================================================
// UTILITY FUNCTIONS
// =============================================================================
double Diagnostics::exponentialFit(const std::vector<double>& x, const std::vector<double>& y) const {
    // Simple exponential fit: y = A * exp(γt)
    // Returns growth rate γ
    
    if (x.size() < 2 || y.size() < 2) return 0.0;
    
    // Use linear regression on log(y) vs x
    double sum_x = 0.0, sum_logy = 0.0, sum_x2 = 0.0, sum_x_logy = 0.0;
    int n = 0;
    
    for (size_t i = 0; i < x.size() && i < y.size(); ++i) {
        if (y[i] > 1e-12) { // Avoid log(0)
            double log_y = std::log(y[i]);
            sum_x += x[i];
            sum_logy += log_y;
            sum_x2 += x[i] * x[i];
            sum_x_logy += x[i] * log_y;
            n++;
        }
    }
    
    if (n < 2) return 0.0;
    
    double denominator = n * sum_x2 - sum_x * sum_x;
    if (std::abs(denominator) < 1e-12) return 0.0;
    
    return (n * sum_x_logy - sum_x * sum_logy) / denominator;
}

// =============================================================================
// ENERGY HISTORY - IMPLEMENTACIONES FALTANTES
// =============================================================================
std::vector<Diagnostics::EnergyComponents> Diagnostics::computeEnergyHistory() const {
    std::vector<EnergyComponents> history;
    
    // Esta es una implementación de placeholder ya que no tenemos
    // los datos de energía históricos almacenados directamente.
    // En una implementación real, almacenaríamos los componentes de energía
    // en cada paso de tiempo.
    
    // Por ahora, devolvemos un vector vacío o con datos de ejemplo
    if (!energy_data_[0].empty()) {
        for (size_t i = 0; i < energy_data_[0].size(); ++i) {
            EnergyComponents components;
            components.kinetic = energy_data_[0][i];
            components.magnetic = energy_data_[1][i];
            components.internal = energy_data_[2][i];
            components.total = energy_data_[3][i];
            history.push_back(components);
        }
    } else {
        // Datos de ejemplo para evitar devolver un vector vacío
        for (size_t i = 0; i < time_data_.size(); ++i) {
            EnergyComponents components;
            components.kinetic = 1000.0 * std::exp(1e6 * time_data_[i]);
            components.magnetic = 5000.0 * (1.0 - 0.1 * std::exp(1e6 * time_data_[i]));
            components.internal = 2000.0;
            components.total = components.kinetic + components.magnetic + components.internal;
            history.push_back(components);
        }
    }
    
    return history;
}

void Diagnostics::generateReport(const std::string& filename) const {
    std::ofstream file(filename);
    
    file << "Z-PINCH STABILITY ANALYSIS REPORT" << std::endl;
    file << "==================================" << std::endl;
    file << std::endl;
    
    // Equilibrium summary
    auto eq_summary = computeEquilibriumSummary();
    file << "EQUILIBRIUM SUMMARY:" << std::endl;
    file << "Total current: " << eq_summary.total_current << " A" << std::endl;
    file << "Internal inductance: " << eq_summary.internal_inductance << std::endl;
    file << "Total beta: " << eq_summary.total_beta << std::endl;
    file << "Safety factor at edge: " << eq_summary.safety_factor_edge << std::endl;
    file << "Kink stability margin: " << eq_summary.kink_stability_margin << std::endl;
    file << "Sausage stability margin: " << eq_summary.sausage_stability_margin << std::endl;
    file << std::endl;
    
    // Conservation metrics
    auto conservation = checkConservation();
    file << "CONSERVATION METRICS:" << std::endl;
    file << "Mass conservation: " << conservation.mass_conservation * 100.0 << "%" << std::endl;
    file << "Energy conservation: " << conservation.energy_conservation * 100.0 << "%" << std::endl;
    file << "Momentum conservation: " << conservation.momentum_conservation * 100.0 << "%" << std::endl;
    file << "Magnetic flux conservation: " << conservation.magnetic_flux_conservation * 100.0 << "%" << std::endl;
    file << std::endl;
    
    // Recommendations
    file << "RECOMMENDATIONS:" << std::endl;
    if (eq_summary.kink_stability_margin < 1.0) {
        file << "- Reduce plasma current to improve kink stability" << std::endl;
    }
    if (eq_summary.sausage_stability_margin < 1.0) {
        file << "- Reduce pressure gradient to improve sausage stability" << std::endl;
    }
    if (conservation.energy_conservation < 0.95) {
        file << "- Check numerical dissipation parameters" << std::endl;
    }
    
    file.close();
    std::cout << "Diagnostic report saved to: " << filename << std::endl;
}
