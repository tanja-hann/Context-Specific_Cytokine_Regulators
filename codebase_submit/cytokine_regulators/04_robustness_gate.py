"""
Step 4 — Robustness gate -> frozen 128-gene regulator list.

Selectivity alone can reward genes with a single, barely-significant,
non-reproducible cytokine hit. A gene-level robustness gate removes these
artifacts, keeping only regulators whose knockdown is both real and
reproducible (using QC metrics the dataset authors precomputed per
perturbation in the stimulated condition):

    KEEP gene IF  ontarget_significant == True        (guide knocked down its
                                                        intended target)
             AND  donor_correlation_all_mean > 0.05   (effect reproducible
                                                        across donors)

The surviving 128 genes are the frozen regulator list used by all downstream
analysis (selectivity landscape figure, druggability scoring).

Inputs : curated_hits_selectivity_v2.csv, DE_obs_metadata.parquet
Output : curated_hits_robust_selectivity.csv
"""
import pandas as pd

HITS = "curated_hits_selectivity_v2.csv"
OBS_META = "DE_obs_metadata.parquet"
DONOR_CORR_MIN = 0.05


def main():
    final = pd.read_csv(HITS)
    obs_meta = pd.read_parquet(OBS_META)
    stim = obs_meta[obs_meta["culture_condition"] == "Stim8hr"].set_index(
        "target_contrast_gene_name")

    q = stim.reindex(final["gene"].unique())
    gate = (q["ontarget_significant"] == True) & (q["donor_correlation_all_mean"] > DONOR_CORR_MIN)
    robust_genes = set(q.index[gate])

    final["ontarget_sig"] = final["gene"].map(stim["ontarget_significant"])
    final["donor_corr"] = final["gene"].map(stim["donor_correlation_all_mean"]).round(3)
    final["robust"] = final["gene"].isin(robust_genes)

    gated = final[final["robust"]].copy().sort_values("selectivity_index", ascending=False)
    gated.to_csv("curated_hits_robust_selectivity.csv", index=False)
    print(f"robust hits: {len(gated)} | genes: {gated['gene'].nunique()} "
          f"| {gated['direction'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
