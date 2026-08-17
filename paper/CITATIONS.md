# Companion Article

## Xu, C. et al. (2023) — Electrical circuit simulation of non-Hermitian lattice models

This repository was built while preparing the following first-author **review article**:

> C. Xu, Z. Xu, Z. Zhou, E. Cheng, L. Lang, "Electrical circuit simulation of non-Hermitian lattice models" (非厄米格点模型的经典电路模拟), *Acta Physica Sinica* **72**(20), 200301 (2023). DOI: [10.7498/aps.72.20230914](https://doi.org/10.7498/aps.72.20230914)

BibTeX:

```bibtex
@article{xu2023electrical,
  title   = {Electrical circuit simulation of non-Hermitian lattice models},
  author  = {Xu, Canhong and Xu, Zhicong and Zhou, Ziyu and Cheng, Enhong and Lang, Lijun},
  journal = {Acta Physica Sinica},
  volume  = {72},
  number  = {20},
  pages   = {200301},
  year    = {2023},
  doi     = {10.7498/aps.72.20230914}
}
```

## Relationship to this repository

- The article is a **review** of classical-circuit emulation of non-Hermitian lattice models (mathematical foundations, circuit-mapping theory, experimental progress on PT symmetry, skin effect, and non-Hermitian topology).
- **The simulation code in this repository is an independent implementation** created to build hands-on understanding of the physics surveyed in the article. It does **not** reproduce the figures of the article.
- The simulations (4N+1-site SSH chain with gain/loss) are the author's own numerical work and can be cited independently as learning/research materials.

## Model source — the cSSH model

The 4N+1-site chains in `matlab/` and `python/` are numerical realizations of the
**complex SSH (cSSH)** model — SSH lattice with alternating on-site gain/loss — studied in:

> L.-J. Lang, Y. Wang, H. Wang, Y. D. Chong, "Effects of non-Hermiticity on Su-Schrieffer-Heeger defect states", *Phys. Rev. B* **98**, 094307 (2018). DOI: [10.1103/PhysRevB.98.094307](https://doi.org/10.1103/PhysRevB.98.094307)

BibTeX:

```bibtex
@article{lang2018effects,
  title   = {Effects of non-Hermiticity on Su-Schrieffer-Heeger defect states},
  author  = {Lang, Li-Jun and Wang, You and Wang, Hailong and Chong, Y. D.},
  journal = {Physical Review B},
  volume  = {98},
  number  = {9},
  pages   = {094307},
  year    = {2018},
  doi     = {10.1103/PhysRevB.98.094307}
}
```

**Relationship:** Lang et al. (2018) analyze the mid-gap defect state localized at
the interface between two cSSH domains and its fate under strong gain/loss
(disappearance into the complex continuum; spontaneous ST-symmetry breaking at an
exceptional point producing a pair of defect states). The 4N+1-site geometry in this
repository implements that two-domain configuration; variants `a`/`b`/`c` differ in
the gain/loss phase pattern (same-phase PT / anti-phase PT / anti-phase with $b_N$
gaining). The author's review (Xu et al. 2023, above) situates the cSSH model within
the broader experimental landscape of classical-circuit emulation of non-Hermitian
lattices.
