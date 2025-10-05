import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_equilibrium_profiles():
    # Cargar datos de equilibrio
    df = pd.read_csv('data/equilibrium_profiles.csv')
    
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    
    # Presión
    axes[0,0].plot(df['r'], df['p'])
    axes[0,0].set_title('Pressure Profile')
    axes[0,0].set_ylabel('p (Pa)')
    axes[0,0].grid(True)
    
    # Campo magnético axial
    axes[0,1].plot(df['r'], df['Bz'])
    axes[0,1].set_title('Axial Magnetic Field')
    axes[0,1].set_ylabel('Bz (T)')
    axes[0,1].grid(True)
    
    # Campo magnético azimutal
    axes[0,2].plot(df['r'], df['Btheta'])
    axes[0,2].set_title('Azimuthal Magnetic Field')
    axes[0,2].set_ylabel('Bθ (T)')
    axes[0,2].grid(True)
    
    # Densidad de corriente
    axes[1,0].plot(df['r'], df['jz'])
    axes[1,0].set_title('Current Density')
    axes[1,0].set_ylabel('jz (A/m²)')
    axes[1,0].grid(True)
    
    # Factor de seguridad
    axes[1,1].plot(df['r'], df['q'])
    axes[1,1].set_title('Safety Factor')
    axes[1,1].set_ylabel('q')
    axes[1,1].grid(True)
    
    # Parámetro beta
    axes[1,2].plot(df['r'], df['beta'])
    axes[1,2].set_title('Beta Parameter')
    axes[1,2].set_ylabel('β')
    axes[1,2].grid(True)
    
    for ax in axes.flat:
        ax.set_xlabel('Radius (m)')
    
    plt.tight_layout()
    plt.savefig('equilibrium_profiles.png', dpi=300)
    plt.show()

if __name__ == "__main__":
    plot_equilibrium_profiles()
