import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import re
from matplotlib.animation import PillowWriter
from matplotlib import cm
from tqdm import tqdm
from scipy.stats import linregress

# -------------------------------
# ANÁLISIS COMPLETO DE TIME_HISTORY
# -------------------------------
def analyze_time_history():
    """Comprehensive analysis of time_history.csv"""
    time_history_file = 'data/time_history.csv'
    
    if not os.path.exists(time_history_file):
        print(f"Time history file {time_history_file} not found!")
        return None
    
    df = pd.read_csv(time_history_file)
    print("=" * 60)
    print("ANÁLISIS COMPLETO DE TIME_HISTORY.CSV")
    print("=" * 60)
    
    # Información básica
    print(f"📊 Datos: {len(df)} puntos temporales")
    print(f"⏰ Tiempo total de simulación: {df['time'].iloc[-1]:.2e} s")
    print(f"📈 Amplitud inicial: {df['kink_amplitude'].iloc[0]:.2e}")
    print(f"🚀 Amplitud final: {df['kink_amplitude'].iloc[-1]:.2e}")
    print(f"🔺 Crecimiento total: {df['kink_amplitude'].iloc[-1]/df['kink_amplitude'].iloc[0]:.2f}x")
    
    # Análisis de crecimiento
    mask = df['kink_amplitude'] > df['kink_amplitude'].iloc[0] * 1.1  # Solo donde crece significativamente
    if mask.sum() > 10:
        growth_fit = linregress(df['time'][mask], np.log(df['kink_amplitude'][mask]))
        print(f"📐 Growth rate del ajuste exponencial: {growth_fit.slope:.2e} s⁻¹")
        print(f"✅ Calidad del ajuste (R²): {growth_fit.rvalue**2:.4f}")
    
    # Análisis energético
    energy_std = df['total_energy'].std() / df['total_energy'].mean()
    print(f"⚡ Estabilidad energía total: {energy_std:.2e} (ideal < 0.01)")
    
    # Growth rate estadísticas
    valid_growth = df[df['growth_rate'] > 100]  # Filtrar valores significativos
    if len(valid_growth) > 0:
        print(f"🌡️  Growth rate promedio: {valid_growth['growth_rate'].mean():.2e} s⁻¹")
        print(f"📊 Growth rate máximo: {valid_growth['growth_rate'].max():.2e} s⁻¹")
    
    return df

# -------------------------------
# ANIMACIÓN 2D OPTIMIZADA
# -------------------------------
def create_evolution_animation():
    snap_dir = 'data/snap'
    snap_files = sorted(
        [f for f in os.listdir(snap_dir) if f.startswith('output_state_') and f.endswith('.csv')],
        key=lambda x: int(re.findall(r'\d+', x)[0])
    )
    
    first_snap = pd.read_csv(os.path.join(snap_dir, snap_files[0]))
    r_unique = np.unique(first_snap['r'])
    z_unique = np.unique(first_snap['z'])
    Nr, Nz = len(r_unique), len(z_unique)
    R, Z = np.meshgrid(r_unique, z_unique, indexing='ij')

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    titles = ['Radial Velocity (m/s)', 'Radial Magnetic Field (T)', 'Pressure (Pa)', 'Azimuthal Magnetic Field (T)']
    cmaps = ['RdBu_r', 'RdBu_r', 'viridis', 'plasma']

    df = first_snap
    vr_2d = df['vr'].values.reshape(Nr, Nz)
    br_2d = df['br'].values.reshape(Nr, Nz)
    pressure_2d = df['pressure'].values.reshape(Nr, Nz)
    btheta_2d = df['btheta'].values.reshape(Nr, Nz)
    variables = [vr_2d, br_2d, pressure_2d, btheta_2d]

    # Inicialización de contornos
    contours = []
    for i, ax in enumerate(axes.flat):
        cont = ax.contourf(Z, R, variables[i], levels=50, cmap=cmaps[i])
        ax.set_title(titles[i])
        ax.set_xlabel('Z (m)')
        if i % 2 == 0:
            ax.set_ylabel('Radius (m)')
        fig.colorbar(cont, ax=ax)
        contours.append(cont)

    writer = PillowWriter(fps=5)
    print("Creating animation...")

    with writer.saving(fig, "zpinch_evolution.gif", dpi=120):
        for i in tqdm(range(len(snap_files)), desc="Saving frames"):
            df = pd.read_csv(os.path.join(snap_dir, snap_files[i]))
            vr_2d = df['vr'].values.reshape(Nr, Nz)
            br_2d = df['br'].values.reshape(Nr, Nz)
            pressure_2d = df['pressure'].values.reshape(Nr, Nz)
            btheta_2d = df['btheta'].values.reshape(Nr, Nz)
            variables = [vr_2d, br_2d, pressure_2d, btheta_2d]

            for j, ax in enumerate(axes.flat):
                ax.clear()
                contours[j] = ax.contourf(Z, R, variables[j], levels=50, cmap=cmaps[j])
                ax.set_title(f'{titles[j]}\nTime: {i*1e-6:.2e} s')
                ax.set_xlabel('Z (m)')
                if j % 2 == 0:
                    ax.set_ylabel('Radius (m)')

            writer.grab_frame()

    print("Animation saved as 'zpinch_evolution.gif'")

# -------------------------------
# GRÁFICOS ESTÁTICOS
# -------------------------------
def create_static_plots():
    snap_dir = 'data/snap'
    snap_files = sorted([f for f in os.listdir(snap_dir) if f.startswith('output_state_') and f.endswith('.csv')],
                        key=lambda x: int(re.findall(r'\d+', x)[0]))
    
    first_snap = pd.read_csv(os.path.join(snap_dir, snap_files[0]))
    r_unique = np.unique(first_snap['r'])
    z_unique = np.unique(first_snap['z'])
    Nr, Nz = len(r_unique), len(z_unique)
    R, Z = np.meshgrid(r_unique, z_unique, indexing='ij')

    indices = [0, len(snap_files)//2, -1]
    snapshot_names = ['Initial', 'Middle', 'Final']
    
    fig, axes = plt.subplots(3, 4, figsize=(18, 12))
    fig.suptitle('Z-Pinch Evolution - Key Snapshots', fontsize=16)
    
    titles = ['Radial Velocity (m/s)', 'Radial Magnetic Field (T)', 
              'Pressure (Pa)', 'Azimuthal Magnetic Field (T)']
    cmaps = ['RdBu_r', 'RdBu_r', 'viridis', 'plasma']

    for row, idx in enumerate(indices):
        if idx >= len(snap_files):
            continue
        df = pd.read_csv(os.path.join(snap_dir, snap_files[idx]))
        time = int(re.findall(r'\d+', snap_files[idx])[0]) * 1e-6

        vr_2d = df['vr'].values.reshape(Nr, Nz)
        br_2d = df['br'].values.reshape(Nr, Nz)
        pressure_2d = df['pressure'].values.reshape(Nr, Nz)
        btheta_2d = df['btheta'].values.reshape(Nr, Nz)
        variables = [vr_2d, br_2d, pressure_2d, btheta_2d]

        for col in range(4):
            ax = axes[row, col]
            cont = ax.contourf(Z, R, variables[col], levels=50, cmap=cmaps[col])
            ax.set_title(f'{titles[col]}')
            ax.set_xlabel('Z (m)')
            if col == 0:
                ax.set_ylabel(f'{snapshot_names[row]}\nR (m)\nt = {time:.2e} s')
            plt.colorbar(cont, ax=ax)
    
    plt.tight_layout()
    plt.savefig('zpinch_static_comparison.png', dpi=150, bbox_inches='tight')
    plt.show()
    print("Static comparison saved as 'zpinch_static_comparison.png'")

# -------------------------------
# PLOT TIME EVOLUTION
# -------------------------------
def plot_time_evolution():
    time_history_file = 'data/time_history.csv'
    
    if not os.path.exists(time_history_file):
        print(f"Time history file {time_history_file} not found!")
        return
    
    df = pd.read_csv(time_history_file)
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Z-Pinch Time Evolution Diagnostics', fontsize=16)
    
    ax1.plot(df['time'], df['kink_amplitude'], 'b-', linewidth=2, label='Kink Amplitude')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Kink Amplitude', color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.grid(True, alpha=0.3)
    ax1.set_yscale('log')
    ax1.legend(loc='upper left')
    
    ax1_twin = ax1.twinx()
    ax1_twin.plot(df['time'], df['growth_rate'], 'r--', linewidth=2, label='Growth Rate')
    ax1_twin.set_ylabel('Growth Rate (s⁻¹)', color='r')
    ax1_twin.tick_params(axis='y', labelcolor='r')
    ax1_twin.legend(loc='upper right')
    
    ax2.plot(df['time'], df['kinetic_energy'], label='Kinetic Energy')
    ax2.plot(df['time'], df['magnetic_energy'], label='Magnetic Energy')
    ax2.plot(df['time'], df['internal_energy'], label='Internal Energy')
    ax2.plot(df['time'], df['total_energy'], 'k-', linewidth=2, label='Total Energy')
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Energy (J)')
    ax2.set_yscale('log')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_title('Energy Evolution')
    
    ax3.plot(df['time'], df['growth_rate'], 'g-', linewidth=2)
    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel('Growth Rate (s⁻¹)')
    ax3.grid(True, alpha=0.3)
    ax3.set_title('Growth Rate Evolution')
    
    initial_amp = df['kink_amplitude'].iloc[0] if len(df) > 0 else 1
    growth_factor = df['kink_amplitude'] / initial_amp
    ax4.plot(df['time'], growth_factor, 'purple', linewidth=2)
    ax4.set_xlabel('Time (s)')
    ax4.set_ylabel('Growth Factor (A/A₀)')
    ax4.set_yscale('log')
    ax4.grid(True, alpha=0.3)
    ax4.set_title(f'Kink Amplitude Growth (Initial: {initial_amp:.2e})')
    
    plt.tight_layout()
    plt.savefig('time_evolution_plots.png', dpi=150, bbox_inches='tight')
    plt.show()
    print("Time evolution plots saved as 'time_evolution_plots.png'")

# -------------------------------
# MAIN
# -------------------------------
if __name__ == "__main__":
    print("=" * 60)
    print("Z-Pinch Stability Simulation - Complete Analysis")
    print("=" * 60)
    
    # Análisis detallado del time_history
    df = analyze_time_history()
    
    # Crear animación
    create_evolution_animation()
    
    # Crear gráficos estáticos
    # create_static_plots()
    
    # Plot time evolution
    # plot_time_evolution()
    
    print("All analyses completed!")

