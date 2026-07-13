# Cytokine-regulator discovery

Discovers genes whose knockdown changes a specific cytokine **only in
stimulated T cells, not at rest**, from the genome-scale CD4+ T-cell
Perturb-seq screen.

## Approach

A context-specific regulator is defined by an interaction, not by comparing
two states separately: we test whether a knockdown's effect on a cytokine
**differs between resting and stimulated cells**, then keep only effects that
are present under stimulation, absent at rest, focused on individual cytokines,
and reproducible.

## Steps

| Script | Does |
|---|---|
| `01_extract_cytokine_DE.py` | Stream precomputed DESeq2 effect sizes from S3, slice to the cytokine module. Writes `cytokine_DE_long.parquet` + `DE_obs_metadata.parquet`. |
| `02_interaction_test.py` | Per gene x cytokine, Wald test on the difference of Rest vs Stim8hr log fold-changes (perturbation x condition interaction), Benjamini-Hochberg FDR. |
| `03_hits_and_selectivity.py` | Call context-specific hits and score each gene's selectivity (focused regulators ranked above broad master regulators). |
| `04_robustness_gate.py` | Keep only genes with confirmed on-target knockdown and cross-donor reproducibility -> **frozen 128-gene list**. |
| `05_figures.py` | Selectivity-vs-effect landscape and per-cytokine panels. |

## Key definitions

**Interaction test** (per-condition DESeq2 fits treated as independent):

```
delta = log_fc(Stim8hr) - log_fc(Rest)
z     = delta / sqrt(lfcSE_rest^2 + lfcSE_stim^2)
p     = two-sided Wald;  FDR = Benjamini-Hochberg across all tests
```

**Context-specific hit** (per gene x cytokine): `adj_p_stim < 0.05` **and**
`adj_p_rest >= 0.05` **and** `fdr_interaction < 0.05` **and** cytokine
expressed under stimulation (`baseMean_stim > 5`) **and** the effect
strengthens in stimulation. Direction is recorded: a *suppressor* lowers the
cytokine on knockdown (anti-inflammatory lead), an *inducer* raises it
(immuno-activating lead).

**Selectivity index** (per gene): `0.7 * focus + 0.3 * (1 - footprint)`, where
focus is the share of the largest single-cytokine effect and footprint is the
genome-wide count of genes the knockdown perturbs. Rewards clean,
single-cytokine regulators over broad master regulators.

**Robustness gate**: `ontarget_significant` **and**
`donor_correlation_all_mean > 0.05`.

## Cytokine module

The "Cytokine" category of the dataset's `immune_effector_genes.csv`, keeping
the 30 cytokines that are testable in this dataset (IL17F and IL6 are dropped
for trace expression, following the source study). The resolved list is in
`data/cytokine_module.csv`.

## Result

312 context-specific hits across 217 genes before gating; **202 robust hits
across 128 genes** after (148 inducer, 54 suppressor); 22 of 30 cytokines have
at least one robust regulator. This 128-gene list is the input to the
druggability stage.

## Inputs

- Streamed from S3: `GWCD4i.DE_stats.h5ad` (step 01).
- `immune_effector_genes.csv` from the source dataset's analysis repo (step
  01); the derived module is cached in `data/cytokine_module.csv`.
