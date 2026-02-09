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

This repository benchmarks several widely used approaches across controlled
simulation settings. The focus is **not** on proposing new methods, but on:

- Principled simulation design
- Careful evaluation logic
- Reproducibility and robustness
- Computational practicality


## Key Design Principles

1. **Reproducibility over interactivity**  
   All experiments are executed via non-interactive batch entrypoints
   (`Rscript scripts/run_*.R`). No manual console steps are required.

2. **Simulation-driven evaluation**  
   Ground-truth precision matrices are explicitly constructed, enabling
   controlled assessment of recovery accuracy, sparsity, and error metrics.

3. **Engineering realism**  
   Computational cost, numerical stability, and package-level issues are treated
   as first-class concerns, reflecting real industry workflows.

4. **Minimal but extensible configuration**  
   A quick configuration (small lambda grids, moderate dimensions) is used by
   default for fast reproducibility. Full-scale runs can be enabled with minor
   parameter changes.


## Methods Benchmarked

Across all cases, the following estimators are evaluated:

- **glasso** (graphical lasso)
- **QUIC**-based precision estimation
- **dpglasso-style** estimation via a QUIC backend

**Note:**  
On some platforms, the `dpglasso` package exhibits a known symmetry-check bug.
To ensure robustness and reproducibility, dpglasso-style estimation is
implemented via a QUIC backend in this repository. This reflects practical
engineering considerations rather than theoretical preference.


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
- Practical handling of numerical and package-level issues

Together, the report and this repository illustrate both statistical reasoning
and applied data science execution.


## Note

This final project was completed jointly with Mingshuo Liu.  
The GitHub repository and code organization were curated, refactored,
and uploaded by Zijie Tian.
