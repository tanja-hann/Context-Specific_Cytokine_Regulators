# Druggability scoring

Scores the 128 context-specific cytokine regulators on two independent axes,
**protein-level** and **RNA-level** druggability, and nominates a therapeutic
modality for each gene.

The premise: a gene that is hard to drug as a protein (no small-molecule
pocket, not on the cell surface) may still be an excellent RNA-level target
(siRNA/ASO knockdown), and vice versa. Scoring both axes routes every gene to
its best-supported modality rather than discarding the "undruggable" ones.

## Steps

| Script | Does |
|---|---|
| `01_protein_score.py` | Protein-level score: small-molecule route (pocket x foldability) vs antibody/biologic route (surface localization x confidence). Writes `protein_scores.tsv`. |
| `02_rna_score.py` | RNA-level score: siRNA/ASO target-site availability x transcript stability x (for siRNA) cytoplasmic localization. Writes `rna_scores.tsv`. |
| `03_matrix_and_nomination.py` | Merge both axes, assign quadrants, nominate a modality, and draw the two-axis matrix. Writes `druggability_summary_table.csv` + figure. |

## Axis definitions

**Protein score** = max(small-molecule, antibody) sub-score, 0–1:

- small molecule = max(OpenTargets pocket, fpocket pocket druggability) x foldability
- antibody/biologic = min(OpenTargets surface localization x confidence, 0.8)

**RNA score** = max(siRNA, ASO) sub-score, 0–1:

- siRNA = site_score x stability x cytoplasmic localization
- ASO   = site_score x stability

`site_score` is the number of accessible, specificity-filtered target sites
found by scanning each transcript, capped and scaled to 0–1.

Both absolute scores are then renormalized across the 128-gene cohort
(z-standardized, centered on 0.5) so the two matrix axes share a common,
cohort-relative scale.

## Quadrants and modality nomination

Each relative axis is split at **0.30**:

| Quadrant | Rule | n |
|---|---|---|
| Both viable | protein >= 0.3 and RNA >= 0.3 | 86 |
| RNA-preferred | protein < 0.3 and RNA >= 0.3 | 20 |
| Protein-preferred | protein >= 0.3 and RNA < 0.3 | 16 |
| Neither viable | both < 0.3 | 6 |

The higher axis wins and its sub-route is nominated. Modality calls across the
128 genes: **ASO 59, small molecule 46, siRNA 14, antibody/biologic 9**.

RNA-preferred genes, those hard to drug as proteins but strong RNA targets,
are the leads for a nucleic-acid strategy (siRNA / ASO), which is the core
translational output of this analysis.

## Feature provenance

The scripts consume two precomputed per-gene feature tables in `data/`;
generating them is upstream of this codebase:

- `protein_axis_features.csv`: OpenTargets tractability, plus structure-based
  pocket druggability (fpocket on AlphaFold models) and foldability/surface
  features derived from AlphaFold structures.
- `rna_axis_features.csv`: per-transcript siRNA/ASO target-site counts (from
  accessibility- and specificity-filtered transcript scans), predicted
  transcript stability, and cytoplasmic localization (lncATLAS-derived CN-RCI).

## Known limitation

A validated master transcription factor (GATA3) scores low on both relative
axes (protein 0.12, RNA 0.09) and lands in the "Neither viable" quadrant:
neither axis clears the 0.30 viability threshold. The nomination rule still
reports a nominal best route by argmax (here small molecule), so GATA3 is
counted in the modality tally above, but "Neither viable" is the call that
matters. This reflects an information gap in the three scored mechanisms
(small-molecule pocket, surface accessibility, RNA-targeting sites), not a
biological verdict; mechanisms outside these axes (e.g. targeted protein
degradation) are not scored here.
