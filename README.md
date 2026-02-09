# Precision Matrix Benchmarking

This repository provides a reproducible, industry-facing benchmarking pipeline
for precision matrix estimation under multiple structural regimes. The project
is designed as a work sample for **data science, biostatistics, and applied
statistical research roles**, with an emphasis on **simulation design,
statistical evaluation, and engineering trade-offs**, rather than methodological
novelty.

The codebase demonstrates how to design, implement, debug, and deliver a
statistical simulation pipeline that is robust, interpretable, and directly
usable in **industry settings such as pharma, biotech, and applied research**.


## Overview

Estimating sparse precision (inverse covariance) matrices is a core task in
modern data science and statistical modeling, with applications including:

- Graphical models and network inference
- High-dimensional covariance estimation
- Exploratory structure discovery
- Downstream modeling and uncertainty analysis

This repository benchmarks several established estimation strategies across
controlled simulation settings. The focus is **not** on proposing new methods,
but on:

- Principled simulation design
- Careful evaluation logic
- Reproducibility and robustness
- Computational and numerical practicality


## Key Design Principles

1. **Reproducibility over interactivity**  
   All experiments are executed via non-interactive batch entrypoints
   (`Rscript scripts/run_*.R`). No manual console steps are required.

2. **Simulation-driven evaluation**  
   Ground-truth precision matrices are explicitly constructed, enabling
   controlled assessment of recovery accuracy, sparsity, and estimation error.

3. **Engineering realism**  
   Numerical stability, convergence behavior, and implementation-level issues
   are treated as first-class concerns, reflecting real applied data science
   workflows.

4. **Minimal but extensible configuration**  
   A quick configuration (small tuning grids, moderate dimensions) is used by
   default for fast reproducibility. Full-scale runs can be enabled with minor
   parameter changes.


## Methods Benchmarked

Across all simulation cases, the following three estimators are evaluated:

- **glasso**  
  Standard graphical lasso with an ℓ₁ penalty, implemented via the
  `glasso` R package.

- **dpglasso (custom implementation)**  
  A self-implemented dpglasso estimator based on coordinate-wise updates of the
  precision matrix. This implementation follows the dpglasso formulation and is
  written explicitly in this repository due to a known issue in the
  `dpglasso` package.

- **QUIC**  
  Sparse inverse covariance estimation using the QUIC algorithm, implemented via
  the `QUIC` R package.

**Implementation note:**  
The dpglasso estimator is implemented directly in this repository. Additional
numerical safeguards (e.g. symmetry enforcement, ridge stabilization, and
finite-value checks) are included to ensure stable execution across simulation
settings.


## Simulation Cases

**Case 1: Banded precision matrix**
- Strong local conditional dependence
- Canonical sparse structure
- Computationally favorable regime

**Case 2: Dense precision matrix**
- Near-complete conditional dependence
- Worst-case computational regime
- Highlights scalability limits

**Case 3: Random sparse precision matrices**
- Controlled sparsity levels
- Stochastic graph structure
- Evaluates robustness across random designs

**Case 4: Hub-based network structure**
- Few high-degree hub nodes
- Structured heterogeneity
- Motivated by network and systems biology settings


## Evaluation Metrics

For each estimator and configuration, the following metrics are computed:

- BIC-style objective
- True Positive Rate (TPR)
- False Positive Rate (FPR)
- Sparsity ratio
- Frobenius norm error
- Wall-clock runtime

Metrics are computed against known ground truth, emphasizing interpretability
and diagnostic value rather than single-number performance.


## Repository Structure

```
scripts/
  run_*.R        Runnable entrypoints (batch execution)
  impl_*.R       Implementation and simulation logic

figures/
  plots/         Generated figures from quick configurations

results/
  (ignored)      Reproducible intermediate outputs
```

The separation between entrypoints and implementation mirrors production data
science workflows, where pipelines must be executable, inspectable, and
maintainable.


## How to Run

From the repository root:

```
Rscript scripts/run_case12.R
Rscript scripts/run_case3.R
Rscript scripts/run_case4_hub.R
Rscript scripts/make_figures.R
```

All scripts are non-interactive and can be executed in sequence. Generated
figures are saved under `figures/plots/`.


## Relation to the Final Report

This repository is an engineering-oriented companion to the associated final
report. While the report emphasizes methodological motivation and empirical
findings, this codebase focuses on:

- Executable simulation design
- Reproducible benchmarking
- Clear separation of concerns
- Practical handling of numerical and implementation-level issues

Together, the report and this repository illustrate both statistical reasoning
and applied data science execution.


## Acknowledgement

This final project was completed jointly with **Mingshuo Liu**, who primarily
contributed to literature review and experimental design.  
The GitHub repository, code implementation, and reproducible benchmarking
pipeline were developed, refactored, and uploaded by **Zijie Tian**.

