import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import re
from matplotlib.animation import PillowWriter
from matplotlib import cm
from tqdm import tqdm
from scipy.interpolate import griddata

# -------------------------------
# RADIAL CROSS-SECTION ANALYSIS
# -------------------------------
def analyze_radial_evolution():
    """Analyze radial cross-sections at fixed z position"""
    snap_dir = 'data/snap'
    snap_files = sorted(
        [f for f in os.listdir(snap_dir) if f.startswith('output_state_') and f.endswith('.csv')],
        key=lambda x: int(re.findall(r'\d+', x)[0])
    )
    
    if not snap_files:
        print("No snapshot files found!")
        return
    
    # Read first snapshot to get grid structure
    first_snap = pd.read_csv(os.path.join(snap_dir, snap_files[0]))
    
    # Get unique radial and axial coordinates
    r_unique = np.unique(first_snap['r'])
    z_unique = np.unique(first_snap['z'])
    
    print("=" * 60)
    print("RADIAL CROSS-SECTION ANALYSIS")
    print("=" * 60)
    print(f"Radial grid points: {len(r_unique)}")
    print(f"Axial grid points: {len(z_unique)}")
    print(f"Radial range: {r_unique.min():.3f} to {r_unique.max():.3f} m")
    print(f"Axial range: {z_unique.min():.3f} to {z_unique.max():.3f} m")
    
    # Choose fixed z position for cross-section (middle of the plasma)
    z_cross_section = z_unique[len(z_unique) // 2]
    print(f"Cross-section at z = {z_cross_section:.3f} m")
    
    return r_unique, z_unique, z_cross_section, snap_files

# -------------------------------
# RADIAL PROFILE ANIMATION
# -------------------------------
def create_radial_profile_animation():
    """Create animation of radial profiles evolution"""
    r_unique, z_unique, z_cross_section, snap_files = analyze_radial_evolution()
    
    # Set up the figure
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle('Z-Pinch Radial Cross-Section Evolution', fontsize=16, y=0.95)
    
    # Storage for time evolution data
    time_points = []
    radial_profiles = {
        'vr': [], 'br': [], 'pressure': [], 
        'btheta': [], 'vtheta': [], 'bz': []
    }
    
    print("Extracting radial profiles...")
    
    # Extract radial profiles for each snapshot
    for snap_file in tqdm(snap_files, desc="Processing snapshots"):
        df = pd.read_csv(os.path.join('data/snap', snap_file))
        
        # Extract time from filename (microseconds)
        time_us = int(re.findall(r'\d+', snap_file)[0])
        time_points.append(time_us * 1e-6)  # Convert to seconds
        
        # Get data at fixed z position
        z_slice = df[np.abs(df['z'] - z_cross_section) < 1e-6]
        
        # Sort by radius and store profiles
        z_slice = z_slice.sort_values('r')
        radial_profiles['vr'].append(z_slice['vr'].values)
        radial_profiles['br'].append(z_slice['br'].values)
        radial_profiles['pressure'].append(z_slice['pressure'].values)
        radial_profiles['btheta'].append(z_slice['btheta'].values)
        radial_profiles['vtheta'].append(z_slice['vtheta'].values)
        radial_profiles['bz'].append(z_slice['bz'].values)
    
    # Convert to numpy arrays
    for key in radial_profiles:
        radial_profiles[key] = np.array(radial_profiles[key])
    
    # Create initial plots
    lines = []
    titles = [
        'Radial Velocity (m/s)', 'Radial Magnetic Field (T)', 'Pressure (Pa)',
        'Azimuthal Magnetic Field (T)', 'Azimuthal Velocity (m/s)', 'Axial Magnetic Field (T)'
    ]
    variables = ['vr', 'br', 'pressure', 'btheta', 'vtheta', 'bz']
    colors = ['blue', 'red', 'green', 'purple', 'orange', 'brown']
    
    for i, (ax, var, title, color) in enumerate(zip(axes.flat, variables, titles, colors)):
        line, = ax.plot(r_unique, radial_profiles[var][0], color=color, linewidth=2)
        lines.append(line)
        ax.set_xlabel('Radius (m)')
        ax.set_ylabel(title)
        ax.grid(True, alpha=0.3)
        ax.set_title(title)
        
        # Set appropriate scales
        if var in ['pressure']:
            ax.set_yscale('log')
    
    # Add time annotation
    time_text = fig.text(0.5, 0.01, f'Time: {time_points[0]:.2e} s', 
                        ha='center', fontsize=12, transform=fig.transFigure)
    
    # Create animation
    writer = PillowWriter(fps=10)
    
    print("Creating radial profile animation...")
    
    with writer.saving(fig, "radial_evolution.gif", dpi=120):
        for i in tqdm(range(len(time_points)), desc="Rendering frames"):
            # Update all lines
            for j, var in enumerate(variables):
                lines[j].set_ydata(radial_profiles[var][i])
                
                # Auto-scale y-axis for each frame
                data = radial_profiles[var][i]
                if len(data) > 0:
                    margin = 0.1 * (np.max(data) - np.min(data)) if np.max(data) != np.min(data) else 0.1
                    axes.flat[j].set_ylim(np.min(data) - margin, np.max(data) + margin)
            
            # Update time
            time_text.set_text(f'Time: {time_points[i]:.2e} s')
            
            writer.grab_frame()
    
    plt.close(fig)
    print("Radial evolution animation saved as 'radial_evolution.gif'")
    
    return radial_profiles, time_points, r_unique

# -------------------------------
# RADIAL PROFILE STATIC COMPARISON
# -------------------------------
def create_radial_static_comparison():
    """Create static comparison of radial profiles at key times"""
    r_unique, z_unique, z_cross_section, snap_files = analyze_radial_evolution()
    
    # Select key time points
    indices = [0, len(snap_files)//3, 2*len(snap_files)//3, -1]
    time_labels = ['Initial', 'Early', 'Mid', 'Final']
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle('Radial Profile Comparison at Key Times', fontsize=16)
    
    titles = [
        'Radial Velocity (m/s)', 'Radial Magnetic Field (T)', 'Pressure (Pa)',
        'Azimuthal Magnetic Field (T)', 'Azimuthal Velocity (m/s)', 'Axial Magnetic Field (T)'
    ]
    variables = ['vr', 'br', 'pressure', 'btheta', 'vtheta', 'bz']
    colors = ['blue', 'red', 'green', 'purple', 'orange', 'brown']
    linestyles = ['-', '--', '-.', ':']
    
    for idx, (time_idx, label) in enumerate(zip(indices, time_labels)):
        if time_idx >= len(snap_files):
            continue
            
        df = pd.read_csv(os.path.join('data/snap', snap_files[time_idx]))
        z_slice = df[np.abs(df['z'] - z_cross_section) < 1e-6].sort_values('r')
        time_us = int(re.findall(r'\d+', snap_files[time_idx])[0])
        time_sec = time_us * 1e-6
        
        for i, (ax, var, title) in enumerate(zip(axes.flat, variables, titles)):
            ax.plot(r_unique, z_slice[var].values, 
                   linestyle=linestyles[idx], color=colors[i], 
                   linewidth=2, label=f'{label} (t={time_sec:.1e}s)')
            ax.set_xlabel('Radius (m)')
            ax.set_ylabel(title)
            ax.grid(True, alpha=0.3)
            ax.set_title(title)
            ax.legend()
            
            if var in ['pressure']:
                ax.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig('radial_profile_comparison.png', dpi=150, bbox_inches='tight')
    plt.show()
    print("Radial profile comparison saved as 'radial_profile_comparison.png'")

# -------------------------------
# RADIAL INSTABILITY GROWTH ANALYSIS
# -------------------------------
def analyze_radial_instability_growth():
    """Analyze growth of instabilities at different radii"""
    r_unique, z_unique, z_cross_section, snap_files = analyze_radial_evolution()
    
    print("Analyzing radial instability growth...")
    
    # Storage for amplitude evolution at different radii
    radial_positions = [0.2, 0.5, 0.8]  # Normalized radial positions
    actual_radii = [r_unique[int(pos * len(r_unique))] for pos in radial_positions]
    
    time_points = []
    vr_evolution = {f'r_{r:.3f}': [] for r in actual_radii}
    br_evolution = {f'r_{r:.3f}': [] for r in actual_radii}
    
    for snap_file in tqdm(snap_files, desc="Analyzing growth"):
        df = pd.read_csv(os.path.join('data/snap', snap_file))
        z_slice = df[np.abs(df['z'] - z_cross_section) < 1e-6].sort_values('r')
        
        time_us = int(re.findall(r'\d+', snap_file)[0])
        time_points.append(time_us * 1e-6)
        
        for r_actual in actual_radii:
            # Find closest radial point
            r_idx = np.argmin(np.abs(z_slice['r'] - r_actual))
            key = f'r_{r_actual:.3f}'
            
            vr_evolution[key].append(z_slice['vr'].iloc[r_idx])
            br_evolution[key].append(z_slice['br'].iloc[r_idx])
    
    # Plot radial growth comparison
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    fig.suptitle('Instability Growth at Different Radii', fontsize=16)
    
    colors = ['red', 'blue', 'green']
    for (r_actual, color) in zip(actual_radii, colors):
        key = f'r_{r_actual:.3f}'
        
        # Radial velocity growth
        ax1.semilogy(time_points, np.abs(vr_evolution[key]), 
                    color=color, linewidth=2, label=f'r = {r_actual:.3f} m')
        
        # Magnetic field growth  
        ax2.semilogy(time_points, np.abs(br_evolution[key]),
                    color=color, linewidth=2, label=f'r = {r_actual:.3f} m')
    
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('|Radial Velocity| (m/s)')
    ax1.set_title('Radial Velocity Growth')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('|Radial Magnetic Field| (T)')
    ax2.set_title('Radial Magnetic Field Growth')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('radial_growth_analysis.png', dpi=150, bbox_inches='tight')
    plt.show()
    print("Radial growth analysis saved as 'radial_growth_analysis.png'")
    
    return vr_evolution, br_evolution, time_points, actual_radii

# -------------------------------
# PRESSURE AND MAGNETIC FIELD PROFILES
# -------------------------------
def plot_pressure_magnetic_profiles():
    """Plot detailed pressure and magnetic field radial profiles"""
    r_unique, z_unique, z_cross_section, snap_files = analyze_radial_evolution()
    
    # Select key snapshots
    indices = [0, len(snap_files)//2, -1]
    labels = ['Initial', 'Mid', 'Final']
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('Pressure and Magnetic Field Radial Profiles', fontsize=16)
    
    colors = ['blue', 'orange', 'red']
    
    for idx, (time_idx, label, color) in enumerate(zip(indices, labels, colors)):
        if time_idx >= len(snap_files):
            continue
            
        df = pd.read_csv(os.path.join('data/snap', snap_files[time_idx]))
        z_slice = df[np.abs(df['z'] - z_cross_section) < 1e-6].sort_values('r')
        time_us = int(re.findall(r'\d+', snap_files[time_idx])[0])
        
        # Pressure profile
        axes[0,0].plot(z_slice['r'], z_slice['pressure'], 
                      color=color, linewidth=2, label=f'{label} (t={time_us}μs)')
        
        # Total magnetic field
        b_total = np.sqrt(z_slice['br']**2 + z_slice['btheta']**2 + z_slice['bz']**2)
        axes[0,1].plot(z_slice['r'], b_total, 
                      color=color, linewidth=2, label=f'{label} (t={time_us}μs)')
        
        # Magnetic field components
        axes[1,0].plot(z_slice['r'], z_slice['btheta'], 
                      color=color, linewidth=2, label=f'Bθ {label}')
        axes[1,1].plot(z_slice['r'], z_slice['bz'], 
                      color=color, linewidth=2, label=f'Bz {label}')
    
    # Pressure plot
    axes[0,0].set_xlabel('Radius (m)')
    axes[0,0].set_ylabel('Pressure (Pa)')
    axes[0,0].set_title('Pressure Radial Profile')
    axes[0,0].legend()
    axes[0,0].grid(True, alpha=0.3)
    
    # Total B field
    axes[0,1].set_xlabel('Radius (m)')
    axes[0,1].set_ylabel('|B| (T)')
    axes[0,1].set_title('Total Magnetic Field')
    axes[0,1].legend()
    axes[0,1].grid(True, alpha=0.3)
    
    # Btheta
    axes[1,0].set_xlabel('Radius (m)')
    axes[1,0].set_ylabel('Bθ (T)')
    axes[1,0].set_title('Azimuthal Magnetic Field')
    axes[1,0].legend()
    axes[1,0].grid(True, alpha=0.3)
    
    # Bz
    axes[1,1].set_xlabel('Radius (m)')
    axes[1,1].set_ylabel('Bz (T)')
    axes[1,1].set_title('Axial Magnetic Field')
    axes[1,1].legend()
    axes[1,1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('pressure_magnetic_profiles.png', dpi=150, bbox_inches='tight')
    plt.show()
    print("Pressure and magnetic field profiles saved as 'pressure_magnetic_profiles.png'")

# -------------------------------
# MAIN EXECUTION
# -------------------------------
if __name__ == "__main__":
    print("=" * 60)
    print("Z-Pinch Radial Cross-Section Analysis")
    print("=" * 60)
    
    try:
        # 1. Create radial profile animation
        radial_profiles, time_points, r_unique = create_radial_profile_animation()
        
        # 2. Create static comparison
        create_radial_static_comparison()
        
        # 3. Analyze instability growth at different radii
        vr_evolution, br_evolution, time_points, actual_radii = analyze_radial_instability_growth()
        
        # 4. Plot detailed pressure and magnetic profiles
        plot_pressure_magnetic_profiles()
        
        print("\n" + "=" * 60)
        print("ANALYSIS COMPLETE!")
        print("Generated files:")
        print("  - radial_evolution.gif (Animation)")
        print("  - radial_profile_comparison.png (Static profiles)")
        print("  - radial_growth_analysis.png (Growth rates)")
        print("  - pressure_magnetic_profiles.png (Detailed profiles)")
        print("=" * 60)
        
    except Exception as e:
        print(f"Error during analysis: {e}")
        import traceback
        traceback.print_exc()
