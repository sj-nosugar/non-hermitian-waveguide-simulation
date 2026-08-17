"""Compare the three 4N+1 model variants (a / b / c) at finite gain-loss.

- 'a': same-phase PT parts (a: -iγ, b: +iγ throughout; b_N lossless)
- 'b': anti-phase PT parts (second half flipped)
- 'c': anti-phase PT parts with b_N keeping +iγ

All labels in English. Outputs PNGs into figures_py/.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ssh_lattice import build_hamiltonian, spectrum, time_evolve, demo_defaults

OUT = "figures_py"
import os
os.makedirs(OUT, exist_ok=True)

VARIANTS = ["a", "b", "c"]
COLORS = {"a": "C0", "b": "C1", "c": "C2"}


def main():
    # complex-plane spectra at gamma = 2.5 (non-Hermitian regime)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4), sharex=True, sharey=True)
    for ax, v in zip(axes, VARIANTS):
        p = demo_defaults("c" if v == "c" else "a")
        N = p["N"]
        H = build_hamiltonian(N, p["gamma"], p["J1"], p["J2"], v)
        ev, _ = spectrum(H)
        ax.plot(ev.real, ev.imag, "o", ms=4, mfc="none", mec=COLORS[v])
        ax.axhline(0, color="k", lw=0.6)
        ax.set_title(f"variant {v} (N={N}, $\\gamma$={p['gamma']})")
        ax.set_xlabel(r"$\mathrm{Re}\,E$")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel(r"$\mathrm{Im}\,E$")
    fig.suptitle("Complex spectra of the three 4N+1 variants")
    fig.tight_layout(); fig.savefig(f"{OUT}/variants_spectra.png", dpi=200)
    plt.close(fig)

    # time evolution from b_N: final-state localization per variant
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
    for ax, v in zip(axes, VARIANTS):
        p = demo_defaults("c" if v == "c" else "a")
        N = p["N"]
        H = build_hamiltonian(N, p["gamma"], p["J1"], p["J2"], v)
        n_sites = 4 * N + 1
        psi0 = np.zeros(n_sites, dtype=complex)
        psi0[2 * N] = 1.0
        t, psi_t, _, _ = time_evolve(H, psi0, p["t_max"], p["dt"], 2 * N)
        final_prob = np.abs(psi_t[:, -1]) ** 2
        ax.bar(np.arange(1, n_sites + 1), final_prob, color=COLORS[v], width=1)
        ax.set_title(f"variant {v} (N={N}, $\\gamma$={p['gamma']}, t={p['t_max']})")
        ax.set_xlabel("Lattice site")
        ax.grid(alpha=0.3, axis="y")
    axes[0].set_ylabel("Final probability")
    fig.suptitle("Final-state distribution after time evolution from b_N")
    fig.tight_layout(); fig.savefig(f"{OUT}/variants_final.png", dpi=200)
    plt.close(fig)
    print("variant comparison written to figures_py/variants_*.png")


if __name__ == "__main__":
    main()
