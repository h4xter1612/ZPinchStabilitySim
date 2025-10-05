import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_stability_diagram():
    # Cargar datos de estabilidad
    df = pd.read_csv('data/stability_diagram.csv')
    
    plt.figure(figsize=(12, 5))
    
    # Tasas de crecimiento
    plt.subplot(1, 2, 1)
    plt.plot(df['k'], df['growth_rate_kink'], 'bo-', label='Kink mode', linewidth=2, markersize=4)
    plt.plot(df['k'], df['growth_rate_sausage'], 'rs-', label='Sausage mode', linewidth=2, markersize=4)
    plt.xlabel('Wavenumber k (m⁻¹)')
    plt.ylabel('Growth Rate (s⁻¹)')
    plt.title('Linear Growth Rates')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Frecuencias
    plt.subplot(1, 2, 2)
    plt.plot(df['k'], df['frequency_kink'], 'bo-', label='Kink mode', linewidth=2, markersize=4)
    plt.plot(df['k'], df['frequency_sausage'], 'rs-', label='Sausage mode', linewidth=2, markersize=4)
    plt.xlabel('Wavenumber k (m⁻¹)')
    plt.ylabel('Frequency (rad/s)')
    plt.title('Oscillation Frequencies')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('stability_diagram.png', dpi=300)
    plt.show()
    
    # Encontrar modos más inestables
    max_kink_idx = df['growth_rate_kink'].idxmax()
    max_sausage_idx = df['growth_rate_sausage'].idxmax()
    
    print("Most unstable modes:")
    print(f"Kink: k = {df.loc[max_kink_idx, 'k']:.3f} m⁻¹, γ = {df.loc[max_kink_idx, 'growth_rate_kink']:.2e} s⁻¹")
    print(f"Sausage: k = {df.loc[max_sausage_idx, 'k']:.3f} m⁻¹, γ = {df.loc[max_sausage_idx, 'growth_rate_sausage']:.2e} s⁻¹")

if __name__ == "__main__":
    plot_stability_diagram()
