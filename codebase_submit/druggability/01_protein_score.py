"""
Step 1 (protein axis) — Protein-level druggability score.

Scores each of the 128 regulators for druggability by a protein-directed
modality (small molecule or antibody/biologic), from a precomputed feature
table (see README for feature provenance).

Small-molecule sub-score combines a pocket term with structural foldability:
    SM = max(OpenTargets pocket, fpocket druggability) * foldability
Antibody/biologic sub-score combines surface localization with confidence,
capped at 0.8 (antibodies cannot exceed a practical ceiling here):
    AB = min(OpenTargets localization * surface_confidence, 0.8)
Absolute protein score is the better of the two routes:
    protein_score_abs = max(SM, AB)

A relative score renormalizes the absolute score across the 128-gene cohort
(z-standardized, centered on 0.5; +/-2.5 sigma maps to the [0,1] edges) so the
two axes of the final matrix are on a common, cohort-relative scale.

Input : protein_axis_features.csv
Output: protein_scores.tsv
"""
import numpy as np
import pandas as pd

FEATURES = "data/protein_axis_features.csv"
AB_CAP = 0.8
K_SIGMA = 5.0   # z-score / K, so +/-2.5 sigma spans [0, 1] after the 0.5 offset


def main():
    ax = pd.read_csv(FEATURES)

    fp_raw = ax["fpocket_drug_raw"].fillna(0.0).values
    ax["SM_score"] = (np.maximum(ax["SM_ot"].values, fp_raw) * ax["foldability"].values).round(3)
    ax["AB_eff"] = np.minimum(ax["AB_score"].values * ax["surface_confidence"].values, AB_CAP).round(3)

    ax["protein_score_abs"] = np.maximum(ax["SM_score"], ax["AB_eff"]).round(3)
    route = ax[["SM_score", "AB_eff"]].idxmax(axis=1).map(
        {"SM_score": "small molecule", "AB_eff": "antibody/biologic"})
    ax["best_route"] = np.where(ax["protein_score_abs"] > 0, route, "none")

    x = ax["protein_score_abs"].values
    ax["protein_score"] = np.clip(0.5 + (x - x.mean()) / x.std() / K_SIGMA, 0, 1).round(3)

    out = pd.DataFrame({
        "gene": ax["gene"], "uniprot": ax["acc"],
        "protein_score": ax["protein_score"].round(3),
        "protein_score_abs": ax["protein_score_abs"].round(3),
        "best_route": ax["best_route"],
        "sm_score": ax["SM_score"].round(3), "ab_score": ax["AB_eff"].round(3),
        "sm_ot_pocket": ax["SM_ot"].round(3), "fpocket_drug": ax["fpocket_drug_raw"].round(3),
        "disorder_frac": ax["disorder_frac"].round(3), "foldability": ax["foldability"].round(3),
        "struct_available": np.where(ax["af"], "yes", "no"),
        "ab_ot_loc": ax["AB_score"].round(3),
        "surface_confidence": ax["surface_confidence"].round(3),
    }).sort_values("protein_score", ascending=False).reset_index(drop=True)

    out.to_csv("protein_scores.tsv", sep="\t", index=False)
    print(f"genes: {len(out)} | best_route: {out['best_route'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
