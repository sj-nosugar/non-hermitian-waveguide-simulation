"""
main.py
Example script: Compute band structure for a square lattice of dielectric rods.
"""

import numpy as np
from lattice import SquareLattice
from structure import CircularRods
from pwe_solver import PWESolver
from plotter import plot_band_structure

def main():
    # ====== USER PARAMETERS ======
    a = 1.0                      # lattice constant (arbitrary units)
    eps_rod = 11.9               # dielectric constant of rods (e.g., silicon)
    eps_bg = 1.0                 # background (air)
    radius = 0.2 * a             # rod radius
    N_plane_waves = 5            # plane wave cutoff (gives ~121 waves)
    num_bands = 8                # number of bands to compute
    n_kpoints_per_seg = 30       # resolution along each segment
    # =============================
    
    print("Initializing lattice and structure...")
    lattice = SquareLattice(a=a)
    structure = CircularRods(eps_rod, eps_bg, radius, lattice_constant=a)
    
    print("Generating k-point path (Γ-X-M-Γ)...")
    kpoints, labels, label_indices = lattice.get_kpoints_path(n_points=n_kpoints_per_seg)
    
    print(f"Initializing PWE solver with { (2*N_plane_waves+1)**2 } plane waves...")
    solver = PWESolver(lattice, structure, N_plane_waves=N_plane_waves, polarization='TE')
    
    print("Computing band structure...")
    bands = solver.solve_path(kpoints, num_bands=num_bands)
    
    print("Plotting...")
    plot_band_structure(bands, labels, label_indices, save_path='band_structure.png')
    
    # Optional: print bandgap info
    if bands.shape[1] >= 2:
        max_band1 = np.max(bands[:, 0])
        min_band2 = np.min(bands[:, 1])
        if min_band2 > max_band1:
            print(f"Found bandgap: {max_band1:.4f} – {min_band2:.4f} (a/λ)")
        else:
            print("No complete bandgap for these parameters.")
    
    print("Done.")

if __name__ == "__main__":
    main()
