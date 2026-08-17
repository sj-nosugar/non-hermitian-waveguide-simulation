# Non-Hermitian Waveguide Simulation

**Numerical toolbox for non-Hermitian lattice models on the 4N+1-site SSH chain** — eigenvalues in the complex plane, PT-symmetry breaking, time evolution, and inverse participation ratio (IPR) analysis. Written in MATLAB during my M.Sc. at South China Normal University.

> Built while I was writing my first-author review article *Electrical circuit simulation of non-Hermitian lattice models* (Acta Phys. Sin. 72, 200301 (2023)). The simulations here are **my own independent implementation** for understanding the physics behind the review — they do not reproduce the figures of the article itself.

## What's inside

| Module | Files | What it does |
|--------|-------|-------------|
| `matlab/ssh_4nplus1/` | `m_a/b/c_4Nplus1.m`, `OBC.m` | Core: 4N+1-site SSH Hamiltonian with gain/loss (`γ`) and staggered hopping (`J1, J2`); complex-spectrum analysis, time evolution via matrix exponential, IPR dynamics; open-boundary SSH baseline |
| `matlab/coupling/` | `kappa_cal.m`, `单链t+a,t-a循环变t.m` | Coupling-constant calibration and single-chain hopping modulation |
| `matlab/model-compare/` | `Model_cmp_dimensionless.m`, `model_comparision.m` | Dimensionless comparison of the three model variants (final-state distributions, evolution heatmaps) |
| `comsol/` | `MODELS.md` | COMSOL multiphysics models of 1/2/8/9-waveguide arrays (`.mph`, kept out of git) |
| `python/` | NumPy/SciPy port + Jupyter tutorial | Clean-room Python reimplementation with English labels |
| `paper/` | `CITATION.md` | Citation info of the companion review article |

## Physics in one paragraph

The optics–quantum analogy maps paraxial light propagation in waveguide arrays onto the Schrödinger equation, making photonic lattices a classical emulator of tight-binding quantum models. Non-Hermitian extensions — complex on-site potentials (PT symmetry) or non-reciprocal hoppings (skin effect) — break the usual bulk-boundary correspondence and produce phenomena absent in Hermitian systems. The SSH chain is the canonical starting point; the 4N+1-site construction used here adds gain/loss terms and studies how the complex spectrum, edge states, and localization (IPR) respond.

## Run it

Requires MATLAB (or GNU Octave, which runs the pure-numerical parts without modification).

```bash
# Spectrum, time evolution, IPR — model a
octave matlab/ssh_4nplus1/m_a_4Nplus1.m

# Open-boundary SSH baseline
octave matlab/ssh_4nplus1/OBC.m

# Compare model variants (dimensionless)
octave matlab/model-compare/Model_cmp_dimensionless.m
```

Each script self-contained (Hamiltonian construction and solvers are local functions inside the file).

**Python port** (NumPy/SciPy, English labels, Jupyter tutorial): see [`python/`](python/README.md).

## Sample results

| | |
|---|---|
| **Model a — final-state probability** (N=10, γ=0, J1=0.5, J2=1.5): probability localizes at the central lattice site (0.70), characteristic of the bound state in the 4N+1 geometry | ![Model a final state](figures/Model-a-final.png) |
| **Model comparison** — final-state distributions of variants a/b/c side by side | ![Model comparison](figures/Model-final-comparison.png) |

## Why this matters for my research direction

This project sits at the intersection of **wave physics + numerical modeling + machine learning**: (1) non-Hermitian/photonic systems are my M.Sc. core; (2) the toolbox demonstrates hands-on Hamiltonian-level simulation; (3) the natural next step is learning the dynamics *data-driven* — see [my Neural ODE project](https://github.com/sj-nosugar/neutral-ode-double-pendulum) (chaotic double pendulum, physics→data→NODE→MPC) and the planned extension: *Neural ODE for non-Hermitian waveguide dynamics*.

## License

MIT — see [LICENSE](LICENSE).
