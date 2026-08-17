# Python port (NumPy/SciPy)

A clean-room Python reimplementation of the MATLAB toolbox, with **English figure
labels** throughout. It reproduces the Model-a pipeline and adds an interactive
tutorial.

## Requirements

```bash
pip install numpy scipy matplotlib jupyter
```

## What's here

| File | Purpose |
|------|---------|
| `ssh_lattice.py` | Core library: 4N+1 Hamiltonian (variants a/b/c), open-boundary SSH, spectrum, time evolution, IPR |
| `run_demo.py` | Reproduces the Model-a pipeline: spectrum, PT breaking, time evolution, heatmap, zero modes |
| `compare_variants.py` | Complex spectra + final-state localization for variants a/b/c |
| `tutorial.ipynb` | Interactive Jupyter tutorial (physics → code → figures) |
| `tutorial_defect_states.ipynb` | Deep dive: mid-gap defect state vs γ (Lang 2018), non-unitary evolution, skin effect vs gain/loss |

## Run

```bash
# self-check (Hamiltonian dimension, Hermitian limit, localization)
python ssh_lattice.py

# full Model-a demo -> figures_py/*.png
python run_demo.py

# variant comparison -> figures_py/variants_*.png
python compare_variants.py

# interactive tutorial
jupyter notebook tutorial.ipynb
```

## Variant summary

| Variant | PT phase pattern | b_N diagonal |
|---------|-----------------|--------------|
| `a` | same phase both halves (a: −iγ, b: +iγ) | 0 (lossless) |
| `b` | anti-phase (second half flipped) | 0 (lossless) |
| `c` | anti-phase | +iγ (gain) |

The two-domain geometry and its defect-state physics follow the cSSH model of
Lang et al., PRB 98, 094307 (2018) — see [`paper/CITATIONS.md`](../paper/CITATIONS.md).

## Numerical check

The self-check verifies: (1) Hamiltonian dimension 4N+1; (2) real spectrum at
γ=0 (Hermitian limit); (3) a state initialized on the central site b_N stays
localized (final probability ≈ 0.64, IPR ≈ 0.44 for N=10, γ=0, J1=0.5, J2=1.5),
consistent with the MATLAB results.
