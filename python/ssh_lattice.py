"""Non-Hermitian SSH lattice toolbox (Python port of the MATLAB code).

4N+1-site SSH chain with gain/loss terms, built while writing
    Xu et al., "Electrical circuit simulation of non-Hermitian lattice
    models", Acta Phys. Sin. 72, 200301 (2023).

Model variants
--------------
'a' : same-phase PT parts (a_n: -iγ, b_n: +iγ throughout; b_N lossless)
'b' : anti-phase PT parts (second part flipped: a_n: +iγ, b_n: -iγ)
'c' : anti-phase PT parts, b_N keeps +iγ (no special lossless site)

Only NumPy/SciPy are used; every figure is generated with English labels.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import expm


def build_hamiltonian(N: int, gamma: float, J1: float, J2: float,
                      variant: str = "a") -> np.ndarray:
    """Hamiltonian of the 4N+1-site SSH chain.

    Lattice sites (0-based indices):
        site 0          : a_0  (gain +iγ)
        sites 2n-1, 2n  : a_n, b_n for n = 1..2N  (one-indexed convention)
    Total sites = 4N + 1.

    Parameters
    ----------
    N : number of unit cells per half
    gamma : gain/loss strength (imaginary diagonal part)
    J1, J2 : intra-/inter-cell hoppings
    variant : 'a' | 'b' | 'c' — see module docstring
    """
    n_sites = 4 * N + 1
    H = np.zeros((n_sites, n_sites), dtype=complex)

    # site a_0 (index 0)
    H[0, 0] = +1j * gamma

    # --- first half: n = 1..N, sites a_n=2n-1, b_n=2n ---
    for n in range(1, N + 1):
        a = 2 * n - 1          # a_n (0-based)
        b = 2 * n              # b_n
        # diagonal gain/loss
        H[a, a] = -1j * gamma  # a_n loss
        if variant == "c" or n < N:
            H[b, b] = +1j * gamma          # b_n gain (c: b_N also gains)
        else:
            H[b, b] = 0.0                  # b_N lossless (a, b only)
        # intra-cell hopping J1
        H[a, b] = H[b, a] = J1
        # inter-cell hopping J2
        if n < N:
            a_next = 2 * (n + 1) - 1
            H[b, a_next] = H[a_next, b] = J2

    # a_0 <-> a_1 hopping J2
    H[0, 1] = H[1, 0] = J2

    # --- second half: n = N+1..2N ---
    for n in range(N + 1, 2 * N + 1):
        a = 2 * n - 1
        b = 2 * n
        if variant == "a":
            # same phase as first half: a loss, b gain
            H[a, a] = -1j * gamma
            H[b, b] = +1j * gamma
        else:
            # anti phase: a gain, b loss
            H[a, a] = +1j * gamma
            H[b, b] = -1j * gamma
        # intra-cell hopping J2 (swapped vs first half)
        H[a, b] = H[b, a] = J2
        # inter-cell hopping J1
        if n < 2 * N:
            a_next = 2 * (n + 1) - 1
            H[b, a_next] = H[a_next, b] = J1

    # connect b_N (index 2N) to a_{N+1} (index 2N+1) with J1
    H[2 * N, 2 * N + 1] = H[2 * N + 1, 2 * N] = J1
    return H


def open_boundary_ssh(N: int, t: float, delta: float) -> np.ndarray:
    """Standard open-boundary SSH chain (2N sites), hopping t±delta."""
    H = np.zeros((2 * N, 2 * N), dtype=complex)
    for i in range(N):
        H[2 * i, 2 * i + 1] = H[2 * i + 1, 2 * i] = t + delta
    for i in range(N - 1):
        H[2 * i + 1, 2 * i + 2] = H[2 * i + 2, 2 * i + 1] = t - delta
    return H


def spectrum(H: np.ndarray):
    """Eigenvalues (sorted by real part) and eigenvectors."""
    eigvals, eigvecs = np.linalg.eig(H)
    order = np.argsort(eigvals.real)
    return eigvals[order], eigvecs[:, order]


def time_evolve(H: np.ndarray, psi0: np.ndarray, t_max: float, dt: float,
                target_index: int | None = None):
    """Time evolution via the matrix exponential U = exp(-iHΔt).

    Returns t, psi_t (sites × steps), prob_target, IPR(t).
    """
    H = np.asarray(H, dtype=complex)
    psi0 = np.asarray(psi0, dtype=complex)
    psi0 = psi0 / np.linalg.norm(psi0)
    U = expm(-1j * H * dt)
    t = np.arange(0, t_max + dt / 2, dt)
    n_steps = len(t)

    psi = psi0.copy()
    psi_t = np.zeros((H.shape[0], n_steps), dtype=complex)
    prob_target = np.zeros(n_steps)
    ipr = np.zeros(n_steps)
    for i in range(n_steps):
        psi_t[:, i] = psi
        prob = np.abs(psi) ** 2
        prob_target[i] = prob[target_index] if target_index is not None else 0.0
        ipr[i] = (prob @ prob) / (prob.sum() ** 2)
        psi = U @ psi
    return t, psi_t, prob_target, ipr


def ipr_of_state(psi: np.ndarray) -> float:
    """Inverse participation ratio of a normalized state."""
    prob = np.abs(psi) ** 2
    return float((prob @ prob) / (prob.sum() ** 2))


def demo_defaults(variant: str = "a") -> dict:
    """Default run parameters used in the original MATLAB scripts."""
    if variant == "c":
        return dict(N=20, gamma=2.5, J1=0.5, J2=1.5, t_max=100.0, dt=0.01)
    return dict(N=10, gamma=0.0, J1=0.5, J2=1.5, t_max=1.0, dt=0.01)


if __name__ == "__main__":
    # self-check: Hamiltonian dimension, Hermitian limit, b_N localization
    p = demo_defaults("a")
    N = p["N"]
    H = build_hamiltonian(N, p["gamma"], p["J1"], p["J2"], "a")
    assert H.shape == (4 * N + 1, 4 * N + 1), "wrong Hamiltonian dimension"
    # gamma=0 -> Hermitian -> all eigenvalues real
    eigvals, _ = spectrum(H)
    assert np.allclose(eigvals.imag, 0, atol=1e-12), "gamma=0 must give real spectrum"
    # time evolution from b_N: probability must stay localized near b_N
    psi0 = np.zeros(4 * N + 1, dtype=complex)
    psi0[2 * N] = 1.0
    t, psi_t, prob_bN, ipr = time_evolve(H, psi0, p["t_max"], p["dt"], 2 * N)
    final_prob = np.abs(psi_t[:, -1]) ** 2
    assert final_prob[2 * N] > 0.5, "b_N localization lost"
    assert 0.1 < ipr[-1] < 1.0, "IPR out of range"
    print(f"H {H.shape} | real spectrum: OK | b_N final prob: {final_prob[2 * N]:.3f} "
          f"| final IPR: {ipr[-1]:.3f}")
    print("self-check passed")
