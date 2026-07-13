# Context-specific cytokine regulators as therapeutic targets

Finding genes that control disease-relevant cytokines **only in stimulated
T cells, not at rest**, and ranking them by how they could be drugged.

Context-specific regulators are attractive drug targets: acting only in the
activated state, they can dampen (or boost) a cytokine response without
disturbing resting immune function.

This codebase has two stages:

```
cytokine_regulators/   discover context-specific regulators (128-gene list)
druggability/          score those genes for protein- vs RNA-level druggability
```

## Data

Built on the genome-scale CD4+ T-cell Perturb-seq screen of
Zhu, Dann et al. 2025 (Marson & Pritchard labs), distributed via the CZI
Virtual Cells Platform:

```
s3://genome-scale-tcell-perturb-seq/marson2025_data/
```

The analysis uses the authors' **precomputed differential-expression effect
sizes** (`GWCD4i.DE_stats.h5ad`, DESeq2 fit per perturbation x culture
condition), not raw counts. The 16.8 GB file is streamed from S3 and sliced to
the cytokine genes on the fly, so nothing large is written to disk
(`cytokine_regulators/01_extract_cytokine_DE.py`).

## Pipelines

### 1. Cytokine-regulator discovery (`cytokine_regulators/`)

Tests each perturbation for a **perturbation x condition interaction** on 30
cytokines (Rest vs 8-hour stimulation), calls context-specific hits, scores
their selectivity, and applies a reproducibility gate to freeze a **128-gene
regulator list**. See `cytokine_regulators/README.md`.

### 2. Druggability scoring (`druggability/`)

Scores each of the 128 genes on two independent axes, protein-level (small
molecule / antibody) and RNA-level (siRNA / ASO), places them in a
protein x RNA druggability matrix, and nominates a therapeutic modality per
gene. See `druggability/README.md`.

## Running

Each stage is a set of numbered scripts run in order from its own directory,
e.g.:

```bash
cd cytokine_regulators && for s in 0*.py; do python "$s"; done
cd ../druggability     && for s in 0*.py; do python "$s"; done
```

Dependencies: `numpy`, `pandas`, `pyarrow`, `scipy`, `statsmodels`, `s3fs`,
`h5py`, `matplotlib`, `adjustText`.

## Citation

Please cite the source dataset:

> Zhu R., Dann E. et al. (2025) Genome-scale perturb-seq in primary human CD4+
> T cells maps context-specific regulators of T cell programs and human immune
> traits. *bioRxiv.*

## License

Code in this repository is released under the MIT License (see `LICENSE`). The
source Perturb-seq dataset is the property of its authors and is governed by
its own license and terms of use.
