# COMSOL Models

The `.mph` model files are large binary files and are **not committed to git**. They are listed here for reproducibility — contact the repository owner if you need a specific model.

## Waveguide array models (COMSOL)

| Model | Description |
|-------|-------------|
| `cssh-1wg-beta.mph` | Single waveguide, propagation-constant calibration |
| `cssh-1wg-model.mph` | Single waveguide baseline |
| `cssh-2-coupling-constant.mph` | Two-waveguide coupling-constant extraction (feeds `kappa_cal.m`) |
| `cssh-8wg-4-figc.mph` / `cssh-8wg-45-figc.mph` / `cssh-8wg-5-figc.mph` | 8-waveguide CSSH arrays (model-c figures) |
| `cssh-9wg-5-figa.mph` / `cssh-9wg-5-figb.mph` / `cssh-9wg-5-figa-longtime.mph` | 9-waveguide arrays incl. long-time evolution |
| `microstrip_line_crosstalk_test1.mph` | Microstrip crosstalk test |

## Workflow

1. COMSOL extracts the coupling constants between adjacent waveguides (overlap/coupled-mode analysis) → `kappa_cal.m`
2. MATLAB builds the effective tight-binding Hamiltonian (4N+1 SSH) using those couplings
3. COMSOL full-wave models (8/9-waveguide arrays) validate the effective-model predictions
