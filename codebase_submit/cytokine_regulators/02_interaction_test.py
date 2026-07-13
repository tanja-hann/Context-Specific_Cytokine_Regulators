"""
Step 2 — Perturbation x condition interaction test (context-specificity).

For each perturbation x cytokine, tests whether the knockdown's effect on the
cytokine differs between resting and stimulated (Stim8hr) T cells. Because the
DE object provides per-condition log fold-changes and standard errors (DESeq2
fit separately per condition), the interaction is evaluated as a Wald test on
the difference of two independently-estimated log fold-changes:

    delta    = log_fc(Stim8hr) - log_fc(Rest)
    se_delta = sqrt(lfcSE_rest^2 + lfcSE_stim^2)
    z        = delta / se_delta
    p        = 2 * Phi(-|z|)                      (two-sided)
    fdr      = Benjamini-Hochberg(p)              (across all tests)

Rest vs Stim8hr is used because effector-cytokine transcription peaks early
after TCR activation and this condition carries the strongest cytokine effects.

Two cytokines with trace expression (IL17F, IL6) are excluded following the
dataset authors' own cytokine analysis, leaving a 30-cytokine module.

Input : cytokine_DE_long.parquet
Output: interaction_stim8hr_30cyto.parquet
"""
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

DE_LONG = "cytokine_DE_long.parquet"
STIM_CONDITION = "Stim8hr"
EXCLUDE_CYTOKINES = ["IL17F", "IL6"]


def wide_contrast(df, stim=STIM_CONDITION):
    """Pivot Rest and stimulated effect sizes side by side per gene x cytokine."""
    keep = df[df["culture_condition"].isin(["Rest", stim])]
    idx = ["perturbation", "cytokine"]
    piv = lambda v: keep.pivot_table(index=idx, columns="culture_condition", values=v)
    lfc, se, ap, bm = piv("log_fc"), piv("lfcSE"), piv("adj_p_value"), piv("baseMean")
    out = pd.DataFrame(index=lfc.index)
    out["lfc_rest"], out["lfc_stim"] = lfc["Rest"], lfc[stim]
    out["se_rest"], out["se_stim"] = se["Rest"], se[stim]
    out["ap_rest"], out["ap_stim"] = ap["Rest"], ap[stim]
    out["bm_rest"], out["bm_stim"] = bm["Rest"], bm[stim]
    return out.reset_index()


def main():
    long = pd.read_parquet(DE_LONG)
    long = long[~long["cytokine"].isin(EXCLUDE_CYTOKINES)].copy()
    long = long[np.isfinite(long["log_fc"]) & np.isfinite(long["lfcSE"])
                & (long["lfcSE"] > 0)].copy()
    print(f"cytokines: {long['cytokine'].nunique()}")

    w = wide_contrast(long)
    # Drop self-pairs (a gene tested against its own transcript) and incomplete rows
    w = w[w["perturbation"] != w["cytokine"]].dropna(
        subset=["lfc_rest", "lfc_stim", "se_rest", "se_stim"]).copy()

    w["delta"] = w["lfc_stim"] - w["lfc_rest"]
    w["se_delta"] = np.sqrt(w["se_rest"] ** 2 + w["se_stim"] ** 2)
    w["z_interaction"] = w["delta"] / w["se_delta"]
    w["p_interaction"] = 2 * stats.norm.sf(np.abs(w["z_interaction"]))
    w["z_rest"] = w["lfc_rest"] / w["se_rest"]
    w["z_stim"] = w["lfc_stim"] / w["se_stim"]
    w["fdr_interaction"] = multipletests(w["p_interaction"], method="fdr_bh")[1]

    w.to_parquet("interaction_stim8hr_30cyto.parquet", index=False)
    print(f"tests: {len(w)} | interaction FDR<0.05: {(w['fdr_interaction'] < 0.05).sum()}")


if __name__ == "__main__":
    main()
