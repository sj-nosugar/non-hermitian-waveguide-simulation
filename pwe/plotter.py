"""
Module: plotter.py
Plotting utilities for band structures.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Optional, Tuple

def plot_band_structure(bands: np.ndarray,
                        labels: List[str],
                        label_indices: List[int],
                        gap_highlight: bool = True,
                        save_path: Optional[str] = None):
    """
    Plot the photonic band structure.
    
    Parameters:
    -----------
    bands : ndarray of shape (n_kpoints, n_bands)
    labels : list of str
        Labels for high-symmetry points.
    label_indices : list of int
        Indices where labels should be placed.
    gap_highlight : bool
        If True, highlight the bandgap between first and second band.
    save_path : str or None
        If provided, save figure to this path.
    """
    plt.figure(figsize=(8, 6))
    n_k = bands.shape[0]
    k_index = np.arange(n_k)
    
    # Plot each band
    for b in range(bands.shape[1]):
        plt.plot(k_index, bands[:, b], 'b-', linewidth=1, alpha=0.8)
    
    # Set x-ticks at high-symmetry points
    plt.xticks(label_indices, labels)
    plt.xlim(0, n_k - 1)
    plt.ylabel(r'Normalized frequency $a/\lambda$')
    plt.title('Photonic Band Structure')
    plt.grid(True, alpha=0.3)
    
    # Highlight bandgap if present
    if gap_highlight and bands.shape[1] >= 2:
        max_band1 = np.max(bands[:, 0])
        min_band2 = np.min(bands[:, 1])
        if min_band2 > max_band1:
            plt.axhspan(max_band1, min_band2, alpha=0.2, color='green',
                        label=f'Bandgap: {max_band1:.3f} – {min_band2:.3f}')
            plt.legend()
    
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
    # show only when an interactive display exists (headless-safe)
    import os
    if os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"):
        plt.show()
    else:
        print("[plotter] no display detected; figure saved to", save_path or "memory")  # noqa: E402 (headless environments skip via DISPLAY check)
