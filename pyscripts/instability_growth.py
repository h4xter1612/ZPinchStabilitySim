import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os

def exponential_growth(t, A, gamma):
    return A * np.exp(gamma * t)

def analyze_instability_growth():
    # ✅ USAR time_history.csv EN LUGAR DE evolution_data.csv
    data_file = 'data/time_history.csv'
    
    if not os.path.exists(data_file):
        print(f"ERROR: File {data_file} not found!")
        print("Please run the simulation first to generate time_history.csv")
        return
    
    try:
        df = pd.read_csv(data_file)
        
        # Verificar que el DataFrame no está vacío
        if len(df) == 0:
            print("ERROR: CSV file is empty!")
            return
            
        print(f"✅ Successfully loaded {len(df)} data points from time_history.csv")
        print(f"📊 Columns available: {list(df.columns)}")
        
    except Exception as e:
        print(f"Error loading CSV file: {e}")
        return
    
    # ✅ USAR kink_amplitude PARA EL ANÁLISIS (más significativo que kinetic_energy)
    time = df['time'].values
    amplitude = df['kink_amplitude'].values
    
    # Análisis de crecimiento exponencial
    try:
        # Encontrar fase de crecimiento (excluir valores muy pequeños)
        mask = amplitude > 1e-10 * np.max(amplitude)
        if np.sum(mask) > 10:
            popt, pcov = curve_fit(exponential_growth, time[mask], amplitude[mask], 
                                 p0=[amplitude[mask][0], 1e5])
            growth_rate = popt[1]
            print(f"✅ Exponential growth rate: {growth_rate:.2e} s⁻¹")
            
            # Calcular R² para evaluar calidad del ajuste
            predictions = exponential_growth(time[mask], *popt)
            ss_res = np.sum((amplitude[mask] - predictions) ** 2)
            ss_tot = np.sum((amplitude[mask] - np.mean(amplitude[mask])) ** 2)
            r_squared = 1 - (ss_res / ss_tot)
            print(f"📊 Exponential fit R²: {r_squared:.4f}")
            
        else:
            growth_rate = 0
            print("⚠️ No significant growth detected")
    except Exception as e:
        growth_rate = 0
        print(f"❌ Could not fit exponential growth: {e}")
    
    # Graficar
    plt.figure(figsize=(14, 10))
    
    # 1. Kink amplitude con ajuste exponencial
    plt.subplot(2, 3, 1)
    plt.semilogy(time, amplitude, 'b-', linewidth=2, label='Kink Amplitude')
    if growth_rate > 0:
        t_fit = np.linspace(time[mask][0], time[mask][-1], 100)
        plt.semilogy(t_fit, exponential_growth(t_fit, *popt), 'r--', 
                    label=f'Exponential fit (γ={growth_rate:.2e} s⁻¹)')
    plt.xlabel('Time (s)')
    plt.ylabel('Kink Amplitude')
    plt.title('Kink Mode Growth')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 2. Growth rate history
    plt.subplot(2, 3, 2)
    if 'growth_rate' in df.columns:
        plt.plot(time, df['growth_rate'], 'g-', linewidth=2)
        plt.xlabel('Time (s)')
        plt.ylabel('Growth Rate (s⁻¹)')
        plt.title('Growth Rate Evolution')
        plt.grid(True, alpha=0.3)
    
    # 3. Energías en escala lineal
    plt.subplot(2, 3, 3)
    if 'kinetic_energy' in df.columns:
        plt.plot(time, df['kinetic_energy'], label='Kinetic', linewidth=2)
    if 'magnetic_energy' in df.columns:
        plt.plot(time, df['magnetic_energy'], label='Magnetic', linewidth=2)
    if 'internal_energy' in df.columns:
        plt.plot(time, df['internal_energy'], label='Internal', linewidth=2)
    if 'total_energy' in df.columns:
        plt.plot(time, df['total_energy'], label='Total', linewidth=2, color='black')
    plt.xlabel('Time (s)')
    plt.ylabel('Energy (J)')
    plt.title('Energy Evolution')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 4. Conservación de energía
    plt.subplot(2, 3, 4)
    if 'total_energy' in df.columns and len(df) > 0:
        energy_conservation = df['total_energy'] / df['total_energy'].iloc[0]
        plt.plot(time, energy_conservation)
        plt.xlabel('Time (s)')
        plt.ylabel('Total Energy / Initial Energy')
        energy_error = (energy_conservation.max() - 1.0) * 100
        plt.title(f'Energy Conservation\nMax error: {energy_error:.4f}%')
        plt.grid(True, alpha=0.3)
    
    # 5. Comparación growth rates
    plt.subplot(2, 3, 5)
    if 'growth_rate' in df.columns and growth_rate > 0:
        plt.plot(time, df['growth_rate'], 'g-', label='Instantaneous', linewidth=2)
        plt.axhline(y=growth_rate, color='r', linestyle='--', 
                   label=f'Exponential fit: {growth_rate:.2e} s⁻¹')
        plt.xlabel('Time (s)')
        plt.ylabel('Growth Rate (s⁻¹)')
        plt.title('Growth Rate Comparison')
        plt.legend()
        plt.grid(True, alpha=0.3)
    
    # 6. Factor de crecimiento
    plt.subplot(2, 3, 6)
    initial_amp = amplitude[0]
    growth_factor = amplitude / initial_amp
    plt.semilogy(time, growth_factor, 'purple', linewidth=2)
    plt.xlabel('Time (s)')
    plt.ylabel('Growth Factor (A/A₀)')
    final_growth = growth_factor[-1]
    plt.title(f'Amplitude Growth\nTotal: {final_growth:.2e}x')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('instability_growth_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Análisis adicional
    print("\n" + "="*50)
    print("COMPREHENSIVE ANALYSIS")
    print("="*50)
    
    if len(df) > 0:
        initial_amp = amplitude[0]
        final_amp = amplitude[-1]
        total_growth = final_amp / initial_amp
        
        print(f"📈 Kink amplitude: {initial_amp:.2e} → {final_amp:.2e}")
        print(f"🔺 Total growth factor: {total_growth:.2e}x")
        
        if 'growth_rate' in df.columns:
            valid_growth = df[df['growth_rate'] > 100]  # Filtrar valores significativos
            if len(valid_growth) > 0:
                avg_growth = valid_growth['growth_rate'].mean()
                max_growth = valid_growth['growth_rate'].max()
                print(f"🌡️  Average growth rate: {avg_growth:.2e} s⁻¹")
                print(f"📊 Maximum growth rate: {max_growth:.2e} s⁻¹")
                print(f"📐 Exponential fit growth rate: {growth_rate:.2e} s⁻¹")
        
        if 'total_energy' in df.columns:
            energy_variation = (df['total_energy'].max() - df['total_energy'].min()) / df['total_energy'].iloc[0]
            print(f"⚡ Energy conservation: {energy_variation*100:.4f}% variation")

if __name__ == "__main__":
    analyze_instability_growth()
