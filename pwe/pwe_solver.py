"""
Module: pwe_solver.py
Implements the plane wave expansion eigenvalue solver for TE modes.
"""

import numpy as np
from scipy.linalg import eigh
from typing import List, Tuple, Optional

from lattice import SquareLattice
from structure import CircularRods

class PWESolver:
    """
    Plane wave expansion solver for 2D photonic crystals.
    """
    def __init__(self, lattice: SquareLattice, structure: CircularRods,
                 N_plane_waves: int = 5, polarization: str = 'TE'):
        """
        Parameters:
        -----------
        lattice : SquareLattice
            Lattice object.
        structure : CircularRods
            Structure object (provides Fourier coefficients).
        N_plane_waves : int
            Cutoff for plane waves (m,n from -N to N).
        polarization : str
            'TE' or 'TM' (currently only TE implemented).
        """
        self.lattice = lattice
        self.structure = structure
        self.N = N_plane_waves
        self.polarization = polarization.upper()
        
        # Generate all plane wave vectors
        self.G_vectors, self.G_indices = lattice.generate_reciprocal_vectors(N_plane_waves)
        self.n_waves = len(self.G_vectors)
        
        # Pre‑compute all Fourier coefficients η_G
        self.eta = structure.compute_all_coeffs(self.G_vectors)
        
        # Create a dictionary for fast lookup of η(G - G')
        self.eta_dict = {}
        for i, Gi in enumerate(self.G_vectors):
            for j, Gj in enumerate(self.G_vectors):
                G_diff = Gi - Gj
                # Use tuple of rounded floats as key (approximate, but fine for small N)
                key = tuple(np.round(G_diff, decimals=8))
                if key not in self.eta_dict:
                    # Find index of G_diff in the list of G_vectors
                    # For numerical stability, we search for the closest vector
                    diff_norms = np.linalg.norm(self.G_vectors - G_diff, axis=1)
                    idx = np.argmin(diff_norms)
                    if diff_norms[idx] < 1e-6:
                        self.eta_dict[key] = self.eta[idx]
                    else:
                        # G_diff not in our set – its Fourier coefficient is zero
                        self.eta_dict[key] = 0.0
    
    def build_matrix(self, k: np.ndarray) -> np.ndarray:
        """
        Build the Hamiltonian matrix for a given k-point.
        
        For TE modes: H_{G,G'} = |k+G| * |k+G'| * η(G-G')
        
        Parameters:
        -----------
        k : ndarray of shape (2,)
            Wavevector in the first Brillouin zone.
            
        Returns:
        --------
        H : ndarray of shape (n_waves, n_waves), complex
        """
        H = np.zeros((self.n_waves, self.n_waves), dtype=complex)
        
        for i, Gi in enumerate(self.G_vectors):
            kGi = k + Gi
            norm_i = np.linalg.norm(kGi)
            for j, Gj in enumerate(self.G_vectors):
                G_diff = Gi - Gj
                key = tuple(np.round(G_diff, decimals=8))
                eta_ij = self.eta_dict.get(key, 0.0)
                H[i, j] = norm_i * np.linalg.norm(k + Gj) * eta_ij
        
        return H
    
    def solve_kpoint(self, k: np.ndarray, num_bands: int = 8) -> np.ndarray:
        """
        Solve for the eigenfrequencies at a single k-point.
        
        Returns:
        --------
        frequencies : ndarray of shape (num_bands,)
            Normalized frequencies a/λ = (ω a)/(2πc).
        """
        H = self.build_matrix(k)
        # Solve hermitian eigenvalue problem
        eigvals, _ = eigh(H)
        # Take the lowest `num_bands` eigenvalues (sorted)
        eigvals = np.sort(eigvals)[:num_bands]
        # Convert ω²/c² to a/λ
        # ω²/c² = eigenvalue, so ω/c = sqrt(eigenvalue)
        # a/λ = (ω a)/(2πc) = (a/(2π)) * (ω/c)
        a = self.lattice.a
        freq_norm = np.sqrt(np.real(eigvals)) * a / (2 * np.pi)
        return freq_norm
    
    def solve_path(self, kpoints: List[np.ndarray], num_bands: int = 8) -> np.ndarray:
        """
        Solve for bands along a path of k-points.
        
        Returns:
        --------
        bands : ndarray of shape (len(kpoints), num_bands)
        """
        bands = []
        for k in kpoints:
            bands.append(self.solve_kpoint(k, num_bands))
        return np.array(bands)
