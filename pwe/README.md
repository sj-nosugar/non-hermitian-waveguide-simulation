# Plane-Wave Expansion (PWE): 2D Photonic-Crystal Band Structures

A self-contained Python implementation of the **plane-wave expansion** method for
computing the band structure of 2D photonic crystals (TE modes, square lattice of
circular rods, Γ–X–M–Γ path).

## Relationship to this repository

This repository covers **two complementary levels of photonic simulation**:

| Perspective | Method | Model | Here |
|-------------|--------|-------|------|
| Discrete lattice | tight-binding / SSH / cSSH | hopping amplitudes, on-site terms | `matlab/`, `python/` (4N+1 chain) |
| Continuous medium | plane-wave expansion (PWE) | periodic dielectric profile ε(r) | `pwe/` (this directory) |

The two are connected: a periodic dielectric structure (solved by PWE) can be
coarse-grained into a tight-binding model (SSH/cSSH) near the band edge — the
coupling constants extracted from COMSOL in `comsol/MODELS.md` play exactly that
role. Studying both builds the full chain from Maxwell's equations to lattice
Hamiltonians.

## Files

| File | Purpose |
|------|---------|
| `lattice.py` | Square lattice, reciprocal vectors, Γ–X–M–Γ k-path |
| `structure.py` | Circular-rod geometry, Fourier coefficients η_G (Bessel form factor) |
| `pwe_solver.py` | Matrix assembly and Hermitian eigensolve along the k-path |
| `plotter.py` | Band-structure plot with bandgap highlighting (headless-safe) |
| `main.py` | Example: silicon rods (ε=11.9) in air, r=0.2a |

## Run

```bash
python main.py        # -> band_structure.png
```

Default parameters give a **complete bandgap between band 1 and 2** at
Δ(a/λ) ≈ 0.10 (silicon rods, r = 0.2a), consistent with textbook results.

## Notes / limitations

- **TE polarization only** (TM needs a different matrix form — see section 2.2 of
  the original write-up).
- Plane-wave cutoff N=5 (121 waves) is fine for the lowest bands; increase for
  converged higher bands.
- This is a **learning implementation** of the standard method (Joannopoulos et
  al., *Photonic Crystals*, 2008); it is not part of the published review article.

## Extensions

- TM polarization, hexagonal/triangular lattices
- Bandgap engineering via parameter sweeps (radius, ε)
- ML-driven inverse design (predict bandgaps from geometry) — the original
  write-up flags this as a natural portfolio direction
