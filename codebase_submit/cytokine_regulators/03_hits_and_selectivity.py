"""
Step 3 — Context-specific hit calling and selectivity scoring.

A perturbation x cytokine pair is a context-specific hit when the knockdown:
  - significantly changes the cytokine in stimulated cells   (adj_p_stim  < 0.05)
  - does NOT change it at rest                               (adj_p_rest >= 0.05)
  - has a significant interaction                            (fdr_interaction < 0.05)
  - acts on a cytokine actually expressed in stimulation     (baseMean_stim > 5)
  - strengthens the effect in stimulation      (sign(delta) == sign(lfc_stim))

Direction is recorded, not filtered: a suppressor lowers the cytokine on
knockdown (anti-inflammatory lead); an inducer raises it (immuno-activating lead).

Each hit gene is then scored for selectivity, to reward focused regulators of
individual cytokines over broad master regulators that reshape global T-cell
state:

    focus_ratio = max|log2FC| / sum|log2FC|   over the gene's significant cytokines
    footprint   = n_downstream                (genome-wide count, precomputed by
                                               the dataset authors; not recomputed)
    s_focus     = minmax(focus_ratio)
    s_footprint = 1 - minmax(log1p(n_downstream))
    selectivity_index = 0.7 * s_focus + 0.3 * s_footprint

Min-max normalization is relative to the hit-gene set. By design this ranks
focused single-cytokine regulators above broad lineage transcription factors.

Inputs : interaction_stim8hr_30cyto.parquet, DE_obs_metadata.parquet
Outputs: curated_hits_selectivity_v2.csv  (one row per hit, with selectivity)
         selectivity_scores_v2.csv         (per-gene selectivity components)
"""
import numpy as np
import pandas as pd

INTERACTION = "interaction_stim8hr_30cyto.parquet"
OBS_META = "DE_obs_metadata.parquet"
AP, FDR, BM_MIN = 0.05, 0.05, 5.0


def minmax(s, invert=False):
    s = s.astype(float)
    r = (s - s.min()) / (s.max() - s.min() + 1e-9)
    return 1 - r if invert else r


def score_selectivity(w, footprint):
    """Per-gene focus and footprint terms over the hit-gene set."""
    hits = call_hits(w)
    rows = []
    for g in sorted(hits["perturbation"].unique()):
        gc = w[w["perturbation"] == g]
        sig = gc.loc[gc["ap_stim"] < AP, "lfc_stim"].abs()
        focus = (sig.max() / sig.sum()) if sig.sum() > 0 else np.nan
        rows.append({"gene": g,
                     "n_sig_cytokines": int((gc["ap_stim"] < AP).sum()),
                     "focus_ratio": round(focus, 4),
                     "n_downstream": footprint["n_downstream"].get(g, np.nan)})
    sel = pd.DataFrame(rows).set_index("gene")
    sel["s_focus"] = minmax(sel["focus_ratio"])
    sel["s_footprint"] = minmax(np.log1p(sel["n_downstream"]), invert=True)
    sel["selectivity_index"] = (0.7 * sel["s_focus"] + 0.3 * sel["s_footprint"]).round(4)
    return sel


def call_hits(w):
    base = ((w["ap_stim"] < AP) & (w["ap_rest"] >= AP)
            & (w["fdr_interaction"] < FDR) & (w["bm_stim"] > BM_MIN))
    hits = w[base].copy()
    return hits[np.sign(hits["delta"]) == np.sign(hits["lfc_stim"])].copy()


def main():
    w = pd.read_parquet(INTERACTION)
    obs_meta = pd.read_parquet(OBS_META)
    footprint = obs_meta[obs_meta["culture_condition"] == "Stim8hr"].set_index(
        "target_contrast_gene_name")

    w["direction"] = np.where(w["lfc_stim"] < 0, "suppressor", "inducer")
    hits = call_hits(w)
    sel = score_selectivity(w, footprint)

    final = hits[["perturbation", "cytokine", "direction", "lfc_rest", "lfc_stim",
                  "ap_rest", "ap_stim", "delta", "fdr_interaction", "bm_stim"]].merge(
        sel[["n_sig_cytokines", "focus_ratio", "n_downstream", "selectivity_index"]],
        left_on="perturbation", right_index=True, how="left")
    final = final.rename(columns={"perturbation": "gene"}).sort_values(
        ["direction", "selectivity_index"], ascending=[True, False]).reset_index(drop=True)
    for c in ["lfc_rest", "lfc_stim", "delta"]:
        final[c] = final[c].round(3)

    final.to_csv("curated_hits_selectivity_v2.csv", index=False)
    sel.reset_index().to_csv("selectivity_scores_v2.csv", index=False)
    print(f"hits: {len(final)} | genes: {final['gene'].nunique()} "
          f"| {final['direction'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
