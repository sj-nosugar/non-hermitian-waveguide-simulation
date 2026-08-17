"""
Module: lattice.py
Defines the reciprocal lattice and Brillouin zone for a 2D square lattice.
"""

import numpy as np
from typing import List, Tuple

class SquareLattice:
    """
    Square lattice with lattice constant a.
    """
    def __init__(self, a: float = 1.0):
        """
        Parameters:
        -----------
        a : float
            Lattice constant (arbitrary units).
        """
        self.a = a
        
        # Real space lattice vectors
        self.a1 = np.array([a, 0.0])
        self.a2 = np.array([0.0, a])
        
        # Reciprocal lattice vectors
        self.b1 = (2 * np.pi / a) * np.array([1.0, 0.0])
        self.b2 = (2 * np.pi / a) * np.array([0.0, 1.0])
        
        # High-symmetry points in the first Brillouin zone (Cartesian coordinates)
        self.Gamma = np.array([0.0, 0.0])
        self.X = np.array([np.pi / a, 0.0])
        self.M = np.array([np.pi / a, np.pi / a])
    
    def get_kpoints_path(self, n_points: int = 50) -> Tuple[List[np.ndarray], List[str], List[int]]:
        """
        Generate a list of k-points along the path Gamma-X-M-Gamma.
        
        Parameters:
        -----------
        n_points : int
            Number of points per segment.
            
        Returns:
        --------
        kpoints : list of (2,) ndarray
            List of k-vectors.
        labels : list of str
            Labels for the high-symmetry points.
        label_indices : list of int
            Indices in the kpoints list corresponding to the labels.
        """
        # Segments
        seg1 = np.linspace(0, 1, n_points)   # Gamma -> X
        seg2 = np.linspace(0, 1, n_points)   # X -> M
        seg3 = np.linspace(0, 1, n_points)   # M -> Gamma
        
        kpoints = []
        # Gamma to X
        for t in seg1:
            k = (1 - t) * self.Gamma + t * self.X
            kpoints.append(k)
        # X to M
        for t in seg2:
            k = (1 - t) * self.X + t * self.M
            kpoints.append(k)
        # M to Gamma
        for t in seg3:
            k = (1 - t) * self.M + t * self.Gamma
            kpoints.append(k)
        
        labels = ['Γ', 'X', 'M', 'Γ']
        label_indices = [0, n_points - 1, 2 * n_points - 1, 3 * n_points - 1]
        
        return kpoints, labels, label_indices
    
    def generate_reciprocal_vectors(self, N: int) -> Tuple[np.ndarray, List[Tuple[int, int]]]:
        """
        Generate a list of reciprocal lattice vectors G = m*b1 + n*b2.
        
        Parameters:
        -----------
        N : int
            Cutoff: m,n run from -N to N.
            
        Returns:
        --------
        G_vectors : ndarray of shape (n_vectors, 2)
            List of reciprocal vectors.
        indices : list of (int,int)
            The integer pairs (m,n) for each vector.
        """
        G_vectors = []
        indices = []
        for m in range(-N, N + 1):
            for n in range(-N, N + 1):
                G = m * self.b1 + n * self.b2
                G_vectors.append(G)
                indices.append((m, n))
        return np.array(G_vectors), indices
