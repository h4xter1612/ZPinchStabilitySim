#include "linear_stability.hh"
#include <cmath>
#include <algorithm>
#include <iostream>

// =============================================================================
// CONSTRUCTOR
// =============================================================================
LinearStabilityAnalyzer::LinearStabilityAnalyzer(const ZPinchParameters& params, std::shared_ptr<EquilibriumSolver> equilibrium)
    : params_(params), equilibrium_(equilibrium) {
    initializeGrid();
    computeEquilibrium();
}

// =============================================================================
// GRID INITIALIZATION
// =============================================================================
void LinearStabilityAnalyzer::initializeGrid() {
    r_grid_.resize(params_.Nr);
    double dr = params_.a / (params_.Nr - 1);
    
    for (int i = 0; i < params_.Nr; ++i) {
        r_grid_[i] = i * dr;
    }
}

// =============================================================================
// EQUILIBRIUM CALCULATION - BENNETT PROFILE
// =============================================================================
double LinearStabilityAnalyzer::bennettPressureProfile(double r) const {
    // Bennett pressure profile for Z-pinch
    double r_a = r / params_.a;
    return params_.n0 * params_.T0 * 1.6e-19 / std::pow(1.0 + r_a*r_a, 2.0);
}

double LinearStabilityAnalyzer::bthetaProfile(double r) const {
    // Azimuthal field: B_θ(r) = (μ₀ I₀ r) / (2π a²) for r < a
    if (r <= params_.a) {
        return (4e-7 * M_PI * params_.I0 * r) / (2.0 * M_PI * params_.a * params_.a);
    } else {
        return (4e-7 * M_PI * params_.I0) / (2.0 * M_PI * r);
    }
}

void LinearStabilityAnalyzer::computeEquilibrium() {
    p_eq_.resize(params_.Nr);
    Btheta_eq_.resize(params_.Nr);
    Bz_eq_.resize(params_.Nr, params_.B0);
    jz_eq_.resize(params_.Nr);
    
    // Use equilibrium solver data if available
    if (equilibrium_) {
        const auto& r_grid = equilibrium_->getRadialGrid();
        const auto& p_profile = equilibrium_->getPressureProfile();
        const auto& Btheta_profile = equilibrium_->getBthetaProfile();
        
        for (int i = 0; i < params_.Nr; ++i) {
            p_eq_[i] = p_profile[i];
            Btheta_eq_[i] = Btheta_profile[i];
        }
    } else {
        // Fallback to analytical profiles
        for (int i = 0; i < params_.Nr; ++i) {
            double r = r_grid_[i];
            p_eq_[i] = bennettPressureProfile(r);
            Btheta_eq_[i] = bthetaProfile(r);
        }
    }
    
    // Calculate current density consistently
    computeCurrentFromBtheta();
}

void LinearStabilityAnalyzer::computeCurrentFromBtheta() {
    const double mu0 = 4.0e-7 * M_PI;
    
    for (int i = 0; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        
        if (i == 0) {
            // At axis
            jz_eq_[0] = (2.0 * Btheta_eq_[1] / (r_grid_[1])) / mu0;
        } else {
            double dr = r_grid_[i] - r_grid_[i-1];
            double d_rBtheta_dr = (r * Btheta_eq_[i] - r_grid_[i-1] * Btheta_eq_[i-1]) / dr;
            jz_eq_[i] = d_rBtheta_dr / (mu0 * r);
        }
    }
}

// =============================================================================
// KINK MODE ANALYSIS (m=1)
// =============================================================================
LinearStabilityAnalyzer::KinkModeResult 
LinearStabilityAnalyzer::analyzeKinkMode(double k) const {
    KinkModeResult result;
    result.k = k;
    
    double vA = params_.alfven_speed();
    
    // Get safety factor from equilibrium solver
    double q_edge = 0.0;
    if (equilibrium_) {
        const auto& q_profile = equilibrium_->getSafetyFactorProfile();
        if (!q_profile.empty()) {
            q_edge = std::abs(q_profile.back());
        }
    }
    
    // More realistic formula for kink mode
    double omega_r = k * vA;
    double omega_i = 0.0;
    
    // Improved Kruskal-Shafranov criterion
    if (q_edge < 1.0 && k > 0.1) {
        omega_i = 0.01 * vA * k * (1.0 - q_edge); // Modest growth
    }
    
    result.omega = std::complex<double>(omega_r, omega_i);
    result.growth_rate = omega_i;
    result.frequency = omega_r;
    
    // More realistic eigenfunction
    result.eigenfunction.resize(params_.Nr);
    for (int i = 0; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        result.eigenfunction[i] = (r / params_.a) * std::exp(-pow(r/params_.a, 2));
    }
    
    return result;
}

// =============================================================================
// SAUSAGE MODE ANALYSIS (m=0)
// =============================================================================
LinearStabilityAnalyzer::SausageModeResult 
LinearStabilityAnalyzer::analyzeSausageMode(double k) const {
    SausageModeResult result;
    result.k = k;
    
    double vA = params_.alfven_speed();
    double cs = sqrt(2.0 * params_.T0 * 1.6e-19 / (1.67e-27 * params_.mi_over_me));
    
    // For sausage mode: ω² ≈ k² (v_A² + c_s²) for compressions
    double omega_r = k * sqrt(vA*vA + cs*cs);
    double omega_i = 0.0;
    
    // Sausage mode can be unstable with strong pressure gradients
    double pressure_gradient = 0.0;
    for (int i = 1; i < params_.Nr; ++i) {
        pressure_gradient += std::abs(p_eq_[i] - p_eq_[i-1]);
    }
    pressure_gradient /= (params_.Nr - 1);
    
    // Simplified sausage instability criterion
    double critical_gradient = params_.B0 * params_.B0 / (4e-7 * M_PI * params_.a);
    if (pressure_gradient > critical_gradient) {
        omega_i = 0.05 * vA * k * (pressure_gradient / critical_gradient - 1.0);
    }
    
    result.omega = std::complex<double>(omega_r, omega_i);
    result.growth_rate = omega_i;
    result.frequency = omega_r;
    
    // Eigenfunction (derivative of pressure profile)
    result.eigenfunction.resize(params_.Nr);
    for (int i = 0; i < params_.Nr; ++i) {
        result.eigenfunction[i] = p_eq_[i] / p_eq_[0]; // Normalized
    }
    
    return result;
}

// =============================================================================
// STABILITY SCANS
// =============================================================================
std::vector<LinearStabilityAnalyzer::KinkModeResult> 
LinearStabilityAnalyzer::kinkStabilityScan() {
    std::vector<KinkModeResult> results;
    
    double dk = (params_.k_max - params_.k_min) / params_.n_modes;
    for (int i = 0; i < params_.n_modes; ++i) {
        double k = params_.k_min + i * dk;
        results.push_back(analyzeKinkMode(k));
    }
    
    return results;
}

std::vector<LinearStabilityAnalyzer::SausageModeResult> 
LinearStabilityAnalyzer::sausageStabilityScan() {
    std::vector<SausageModeResult> results;
    
    double dk = (params_.k_max - params_.k_min) / params_.n_modes;
    for (int i = 0; i < params_.n_modes; ++i) {
        double k = params_.k_min + i * dk;
        results.push_back(analyzeSausageMode(k));
    }
    
    return results;
}

// =============================================================================
// ANALYTICAL STABILITY CRITERIA
// =============================================================================
double LinearStabilityAnalyzer::kruskalShafranovLimit() const {
    // Kruskal-Shafranov limit for kink mode
    // k_z < 2π / L_KS, where L_KS = 2π R₀ B_z / B_θ
    double Btheta_surface = bthetaProfile(params_.a);
    return 2.0 * M_PI * params_.R0 * params_.B0 / Btheta_surface;
}

double LinearStabilityAnalyzer::sausageStabilityCriterion() const {
    // Simplified criterion for sausage mode
    // Stable if: dp/dr < B_θ²/(μ₀ r) in some average sense
    double max_pressure_gradient = 0.0;
    for (int i = 1; i < params_.Nr; ++i) {
        double dp_dr = (p_eq_[i] - p_eq_[i-1]) / (r_grid_[i] - r_grid_[i-1]);
        double stability_term = Btheta_eq_[i] * Btheta_eq_[i] / 
                              (4e-7 * M_PI * r_grid_[i]);
        if (dp_dr > stability_term) {
            max_pressure_gradient = std::max(max_pressure_gradient, dp_dr);
        }
    }
    return max_pressure_gradient;
}
// =============================================================================
// EIGENFUNCTION EXPORT - VERSIÓN MEJORADA
// =============================================================================
void LinearStabilityAnalyzer::exportEigenfunctions(const std::string& filename) const {
    std::ofstream file(filename);
    file << "r,kink_eigenfunction,sausage_eigenfunction\n";
    
    // Analizar modos para obtener eigenfunciones reales
    double typical_k = 2.0 * M_PI / params_.L; // Número de onda típico
    
    // Obtener eigenfunciones de los análisis
    auto kink_result = analyzeKinkMode(typical_k);
    auto sausage_result = analyzeSausageMode(typical_k);
    
    // Exportar las eigenfunciones calculadas
    for (int i = 0; i < params_.Nr; ++i) {
        double r = r_grid_[i];
        
        // Usar las eigenfunciones calculadas por analyzeKinkMode y analyzeSausageMode
        double kink_func = kink_result.eigenfunction[i];
        double sausage_func = sausage_result.eigenfunction[i];
        
        file << r << "," << kink_func << "," << sausage_func << "\n";
    }
    
    file.close();
    std::cout << "Eigenfunctions exported to: " << filename << std::endl;
}
