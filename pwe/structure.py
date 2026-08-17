"""
Module: structure.py
Defines the dielectric structure and computes Fourier coefficients of 1/ε.
"""

import numpy as np
from scipy.special import j1
from typing import Union

class CircularRods:
    """
    2D photonic crystal consisting of circular rods in a square lattice.
    """
    def __init__(self, eps_rod: float, eps_background: float,
                 radius: float, lattice_constant: float = 1.0):
        """
        Parameters:
        -----------
        eps_rod : float
            Dielectric constant of the rod material.
        eps_background : float
            Dielectric constant of the background.
        radius : float
            Rod radius (in units of lattice constant if lattice_constant=1).
        lattice_constant : float
            Lattice constant (used to compute filling fraction).
        """
        self.eps_rod = eps_rod
        self.eps_bg = eps_background
        self.r = radius
        self.a = lattice_constant
        
        # Area of unit cell and rod
        self.area_cell = self.a**2
        self.area_rod = np.pi * self.r**2
        self.filling_frac = self.area_rod / self.area_cell
    
    def epsilon_inv_fourier(self, G: np.ndarray) -> complex:
        """
        Fourier coefficient η_G = (1/ε)_G for a square lattice of circular rods.
        
        Parameters:
        -----------
        G : ndarray of shape (2,)
            Reciprocal vector.
            
        Returns:
        --------
        eta : complex
            Fourier coefficient.
        """
        G_norm = np.linalg.norm(G)
        if G_norm < 1e-12:   # G = 0
            return (self.filling_frac / self.eps_rod +
                    (1 - self.filling_frac) / self.eps_bg)
        else:
            # Form factor for circular rod: 2 * J₁(|G|r) / (|G|r)
            form_factor = 2.0 * j1(G_norm * self.r) / (G_norm * self.r)
            return ((1.0 / self.eps_rod - 1.0 / self.eps_bg) *
                    self.filling_frac * form_factor)
    
    def compute_all_coeffs(self, G_vectors: np.ndarray) -> np.ndarray:
        """
        Compute Fourier coefficients for a list of reciprocal vectors.
        
        Parameters:
        -----------
        G_vectors : ndarray of shape (n_vectors, 2)
        
        Returns:
        --------
        coeffs : ndarray of shape (n_vectors,), complex
            η_G for each G.
        """
        n = len(G_vectors)
        coeffs = np.zeros(n, dtype=complex)
        for i, G in enumerate(G_vectors):
            coeffs[i] = self.epsilon_inv_fourier(G)
        return coeffs
