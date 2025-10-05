import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_mode_structures():
    """Plot the radial structure of linear instability modes"""
    try:
        # Load mode structures data
        df = pd.read_csv('data/mode_structures.csv')
        
        # Create figure with subplots
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        
        # Plot kink mode eigenfunction
        ax1.plot(df['r'], df['kink_eigenfunction'], 'b-', linewidth=2, label='Kink mode (m=1)')
        ax1.set_xlabel('Radius (m)')
        ax1.set_ylabel('Eigenfunction Amplitude')
        ax1.set_title('Kink Mode Structure')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Plot sausage mode eigenfunction  
        ax2.plot(df['r'], df['sausage_eigenfunction'], 'r-', linewidth=2, label='Sausage mode (m=0)')
        ax2.set_xlabel('Radius (m)')
        ax2.set_ylabel('Eigenfunction Amplitude')
        ax2.set_title('Sausage Mode Structure')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig('mode_structures.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        print("Mode structure plots saved to 'mode_structures.png'")
        
    except FileNotFoundError:
        print("Mode structures file not found. Run linear stability analysis first.")
        print("Usage: ./zpinch_sim linear --sausage")

def analyze_mode_properties():
    """Analyze properties of the instability modes"""
    try:
        df = pd.read_csv('../data/mode_structures.csv')
        
        # Find maximum amplitudes and their locations
        kink_max_idx = df['kink_eigenfunction'].idxmax()
        sausage_max_idx = df['sausage_eigenfunction'].idxmax()
        
        kink_max_r = df.loc[kink_max_idx, 'r']
        kink_max_amp = df.loc[kink_max_idx, 'kink_eigenfunction']
        
        sausage_max_r = df.loc[sausage_max_idx, 'r']
        sausage_max_amp = df.loc[sausage_max_idx, 'sausage_eigenfunction']
        
        print("\n=== MODE PROPERTIES ===")
        print(f"Kink mode:")
        print(f"  - Maximum amplitude: {kink_max_amp:.4f} at r = {kink_max_r:.4f} m")
        print(f"  - Peak location: {kink_max_r/df['r'].max()*100:.1f}% of plasma radius")
        
        print(f"Sausage mode:")
        print(f"  - Maximum amplitude: {sausage_max_amp:.4f} at r = {sausage_max_r:.4f} m")
        print(f"  - Peak location: {sausage_max_r/df['r'].max()*100:.1f}% of plasma radius")
        
        # Calculate mode widths (FWHM approximation)
        kink_half_max = kink_max_amp / 2
        kink_crossings = df[df['kink_eigenfunction'] > kink_half_max]
        if len(kink_crossings) > 1:
            kink_width = kink_crossings['r'].max() - kink_crossings['r'].min()
            print(f"  - Approximate width: {kink_width:.4f} m")
        
    except FileNotFoundError:
        print("Mode structures file not found.")

if __name__ == "__main__":
    plot_mode_structures()
    analyze_mode_properties()
