"""
Step 2 (RNA axis) — RNA-level druggability score.

Scores each of the 128 regulators for druggability by an RNA-directed modality
(siRNA or ASO), from a precomputed transcript feature table (see README for
feature provenance). Each transcript contributes three factors:

  site_score   fraction of a per-transcript cap (10) of high-quality target
               sites found by scanning the transcript for accessible,
               specificity-filtered windows (siRNA and ASO scored separately)
  stability    predicted transcript stability (favours a durable target)
  localization  cytoplasmic fraction (siRNA/RISC acts in the cytoplasm; scaled
               so that a cytoplasmic proportion >= 0.75 saturates to 1.0)

Modality sub-scores:
    siRNA_score = site_score_siRNA * stability * localization
    ASO_score   = site_score_ASO   * stability
    RNA_score   = max(siRNA_score, ASO_score)
    RNA_modality = argmax(siRNA_score, ASO_score)

A relative score renormalizes RNA_score across the 128-gene cohort exactly as
on the protein axis, so both matrix axes share a cohort-relative [0,1] scale.

Input : rna_axis_features.csv   (per-transcript site counts, stability, localization)
Output: rna_scores.tsv
"""
import numpy as np
import pandas as pd

FEATURES = "data/rna_axis_features.csv"
SITE_CAP = 10
LOC_SATURATION = 0.75
K_SIGMA = 5.0


def main():
    feat = pd.read_csv(FEATURES)

    # Per-transcript site scores: fraction of the cap of usable target sites
    feat["site_score_siRNA"] = (feat["n_sirna_sites"].clip(upper=SITE_CAP) / SITE_CAP).round(4)
    feat["site_score_ASO"] = (feat["n_aso_sites"].clip(upper=SITE_CAP) / SITE_CAP).round(4)

    # Cytoplasmic localization scaled so >= 0.75 cytoplasmic proportion saturates
    loc = (feat["cyto_proportion"] / LOC_SATURATION).clip(upper=1.0)

    feat["siRNA_score"] = (feat["site_score_siRNA"] * feat["stability_score"] * loc).round(4)
    feat["ASO_score"] = (feat["site_score_ASO"] * feat["stability_score"]).round(4)
    feat["RNA_score"] = feat[["siRNA_score", "ASO_score"]].max(axis=1).round(4)
    feat["RNA_modality"] = np.where(feat["siRNA_score"] >= feat["ASO_score"], "siRNA", "ASO")
    feat["localization"] = loc.round(4)

    raw = feat["RNA_score"].values
    feat["RNA_score_relative"] = np.clip(
        0.5 + (raw - raw.mean()) / raw.std(ddof=0) / K_SIGMA, 0, 1).round(4)

    cols = ["gene", "transcript_id", "site_score_siRNA", "site_score_ASO",
            "stability_score", "localization", "siRNA_score", "ASO_score",
            "RNA_score", "RNA_score_relative", "RNA_modality"]
    out = feat[cols].sort_values("RNA_score_relative", ascending=False).reset_index(drop=True)
    for c in ["site_score_siRNA", "site_score_ASO", "stability_score", "localization",
              "siRNA_score", "ASO_score", "RNA_score", "RNA_score_relative"]:
        out[c] = out[c].round(3)

    out.to_csv("rna_scores.tsv", sep="\t", index=False)
    print(f"genes: {len(out)} | RNA_modality: {out['RNA_modality'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
