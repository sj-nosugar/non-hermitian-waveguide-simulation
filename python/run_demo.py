"""Reproduce the Model-a analysis pipeline (MATLAB m_a_4Nplus1.m) in Python.

All figure labels are in English. Outputs PNGs into the current directory.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ssh_lattice import build_hamiltonian, spectrum, time_evolve, demo_defaults

OUT = "figures_py"
import os
os.makedirs(OUT, exist_ok=True)


def main():
    p = demo_defaults("a")
    N, gamma, J1, J2, t_max, dt = (p[k] for k in "N gamma J1 J2 t_max dt".split())
    n_sites = 4 * N + 1

    H = build_hamiltonian(N, gamma, J1, J2, "a")
    eigvals, eigvecs = spectrum(H)

    # --- 1. spectrum (sorted by real part) ---
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(np.arange(1, n_sites + 1), eigvals.real, "b.", ms=12)
    ax.set_xlabel("Eigenstate index")
    ax.set_ylabel(r"$\mathrm{Re}\,E$")
    ax.set_title(f"Eigenvalue spectrum (N={N}, sites={n_sites}, $\\gamma$={gamma})")
    ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(f"{OUT}/spectrum.png", dpi=200); plt.close(fig)

    # --- 2. complex-plane spectrum vs gamma (PT breaking) ---
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    for ax, g in zip(axes, [0.0, 1.0, 3.0]):
        Hg = build_hamiltonian(N, g, J1, J2, "a")
        ev, _ = spectrum(Hg)
        ax.plot(ev.real, ev.imag, "o", ms=5, mfc="none", mec="C0")
        ax.axhline(0, color="k", lw=0.6)
        ax.set_title(rf"$\gamma={g}$")
        ax.set_xlabel(r"$\mathrm{Re}\,E$"); ax.set_ylabel(r"$\mathrm{Im}\,E$")
        ax.grid(alpha=0.3)
    fig.suptitle("Complex spectrum: PT-symmetry breaking with gain/loss")
    fig.tight_layout(); fig.savefig(f"{OUT}/spectrum_complex.png", dpi=200); plt.close(fig)

    # --- 3. time evolution from b_N ---
    psi0 = np.zeros(n_sites, dtype=complex)
    bN = 2 * N                      # b_N site (0-based)
    psi0[bN] = 1.0
    t, psi_t, prob_bN, ipr = time_evolve(H, psi0, t_max, dt, bN)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(t, prob_bN, lw=2)
    ax.set_xlabel("Time $t$"); ax.set_ylabel(r"$|\psi(b_N)|^2$")
    ax.set_title(f"Probability on site b_N (N={N}, $\\gamma$={gamma}, J1={J1}, J2={J2})")
    ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(f"{OUT}/prob_evolution.png", dpi=200); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(t, ipr, lw=2, color="tab:red")
    ax.set_xlabel("Time $t$"); ax.set_ylabel("IPR")
    ax.set_title("Inverse participation ratio vs time")
    ax.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(f"{OUT}/ipr_evolution.png", dpi=200); plt.close(fig)

    # --- 4. density heatmap ---
    prob_density = np.abs(psi_t) ** 2
    fig, ax = plt.subplots(figsize=(9, 4))
    im = ax.imshow(prob_density.T, aspect="auto", origin="lower",
                   extent=[0, n_sites, 0, t_max], cmap="magma")
    ax.set_xlabel("Lattice site"); ax.set_ylabel("Time")
    ax.set_title("Probability density evolution")
    ax.axvline(bN + 0.5, color="w", ls="--", lw=1, alpha=0.6)
    plt.colorbar(im, ax=ax, label="probability")
    fig.tight_layout(); fig.savefig(f"{OUT}/density_heatmap.png", dpi=200); plt.close(fig)

    # --- 5. final-state distribution ---
    final_prob = prob_density[:, -1]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(np.arange(1, n_sites + 1), final_prob, color="steelblue")
    ax.plot(bN + 1, final_prob[bN], "o", color="tab:red", ms=10, mec="k")
    ax.annotate("b_N", (bN + 1, final_prob[bN]), textcoords="offset points",
                xytext=(0, 10), ha="center", color="tab:red")
    ax.set_xlabel("Lattice site"); ax.set_ylabel("Probability")
    ax.set_title(f"Final-state distribution at t={t_max}")
    ax.set_xlim(0, n_sites + 1)
    ax.grid(alpha=0.3, axis="y")
    fig.tight_layout(); fig.savefig(f"{OUT}/final_distribution.png", dpi=200); plt.close(fig)

    # --- 6. zero-energy modes at gamma=0 ---
    H0 = build_hamiltonian(N, 0.0, J1, J2, "a")
    ev0, evec0 = spectrum(H0)
    zero_idx = np.where(np.abs(ev0.real) < 1e-10)[0]
    if len(zero_idx):
        nz = len(zero_idx)
        rows = int(np.ceil(np.sqrt(nz))); cols = int(np.ceil(nz / rows))
        fig, axes = plt.subplots(rows, cols, figsize=(4 * cols, 3 * rows),
                                 squeeze=False)
        for ax, idx in zip(axes.ravel(), zero_idx):
            prob = np.abs(evec0[:, idx]) ** 2
            prob /= prob.sum()
            ax.bar(np.arange(1, n_sites + 1), prob, color="steelblue", width=1)
            ax.set_title(f"$E={ev0[idx]:.3f}$")
            ax.set_xlabel("site"); ax.set_ylabel("probability")
        for ax in axes.ravel()[len(zero_idx):]:
            ax.axis("off")
        fig.suptitle(rf"Zero-energy modes at $\gamma=0$ ({len(zero_idx)} found)")
        fig.tight_layout(); fig.savefig(f"{OUT}/zero_modes.png", dpi=200); plt.close(fig)

    print(f"final b_N probability: {final_prob[bN]:.3f} | final IPR: {ipr[-1]:.3f}")
    print(f"figures written to {OUT}/")


if __name__ == "__main__":
    main()
