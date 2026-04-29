# Environments

## Overview

This project uses **two types of environments** depending on the analysis step:

### 1. Main analysis environment

* **Conda**: provides R and system-level dependencies
* **renv**: manages R package versions

This environment is used for:

* Figure generation (Panels 1E, 1F, 1G)
* All R-based analysis scripts

---

### 2. Contamination scan environment

* Defined in `contam_scan_env.yaml`
* **Conda-only environment** (no renv)

This environment is used for:

* Preprocessing / QC steps
* Contamination detection

---

## 1. Main Analysis Environment

### Step 1 — Create conda environment

```bash
conda env create -f envs/environment.yml
conda activate pherez-farah-2026
```

This installs:

* R
* system libraries required by R packages (e.g. Seurat dependencies)

---

### Step 2 — Restore R packages with renv

Start R in the **repository root**, then run:

```r
install.packages("renv")  # if not already installed
renv::restore()
```

This installs all R packages at the exact versions used in the analysis.

---

### Usage

Before running any analysis script:

```bash
conda activate pherez-farah-2026
```

Then in R:

```r
renv::activate()
```

---

### Lockfiles

* `renv.lock` → exact R package versions (source of truth for R)
* `environment.yml` → base system + R version

---

## 2. Contamination Scan Environment

Environment file:

```
envs/contam_scan_env.yaml
```

### Setup

```bash
conda env create -f envs/contam_scan_env.yaml
conda activate contam-scan
```

---

### Usage

Use this environment **only** for contamination detection / QC steps.

It is **not required** for figure generation unless explicitly stated in a script.

---

## Environment Summary

| Task                      | Environment              |
| ------------------------- | ------------------------ |
| Figure generation (1E–1G) | main (conda + renv)      |
| R-based analysis          | main (conda + renv)      |
| QC / contamination scan   | contam-scan (conda only) |

---

## Notes

* Always activate the correct environment before running scripts
* Do not mix environments
* Scripts should indicate which environment they require
* renv is only used in the main analysis environment
* The R version in `environment.yml` should be compatible with `renv.lock`
