# Sparse Precision Matrix Estimation Benchmark (glasso / QUIC / dpglasso / hubglasso)

This repo is an *industry-facing* simulation benchmark for sparse precision matrix estimators under multiple data-generating regimes:
- **Structured sparsity** (banded)
- **Dense** stress test
- **Unstructured sparsity** (random sparse; varying sparsity levels)
- **Hub-dominated networks** (hubglasso)

The emphasis is on **inference-facing evaluation** (structure recovery, estimation error), **tuning stability**, and **computational behavior**.

## Quickstart

### 1) Restore R environment
```r
install.packages("renv")
renv::restore()
```

### 2) Run simulations (writes outputs to `results/`)
```bash
Rscript scripts/run_case12.R
Rscript scripts/run_case3.R
Rscript scripts/run_case4_hub.R
```

### 3) Generate figures (writes to `figures/`)
```bash
Rscript scripts/make_figures.R
```

## Repository layout

- `R/` — reusable functions (generators, methods wrappers, metrics, tuning)
- `scripts/` — executable entrypoints (one script per experiment family)
- `results/` — **generated** outputs (do not commit large files; see `artifacts/`)
- `figures/` — **generated** figures
- `artifacts/` — small curated results/figures safe to version-control
- `archive/` — historical drafts / old scripts (not part of main pipeline)

## Key findings (fill in once reproduced)

- **Structured sparsity**: methods recover structure well; speed differs substantially.
- **Unstructured sparsity**: all methods struggle to recover zeros; tuning can become unstable.
- **Hub networks**: hubglasso can detect hubs in small p but degrades with dimensionality and is sensitive to tuning.

## Reproducibility notes

- All scripts set seeds and write metadata (session info, run config).
- Prefer storing outputs as `.rds` (structured lists) instead of `save.image()`.

## License
MIT (or your choice).
