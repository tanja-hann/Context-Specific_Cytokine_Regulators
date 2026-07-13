"""
Step 5 — Figures for the context-specific regulator landscape.

Produces two figures from the frozen 128-gene robust list:

  fig_focused_regulator_landscape.png
      One point per robust gene. x = selectivity index, y = log2 fold-change on
      the target cytokine in Stim8hr. Colour encodes knockdown direction
      (red = suppressor / lowers cytokine, blue = inducer / raises cytokine).
      Focused, high-selectivity regulators are direct-labelled.

  fig_30cytokine_panels.png
      One panel per cytokine (30 panels), same axes. Grey points are the
      selectivity-scored genes (the hit-candidate set from
      selectivity_scores_v2.csv, i.e. genes that qualified as a context-specific
      hit for some cytokine and so have a selectivity index); coloured points
      are the robust hits. Shows at a glance which cytokines have focused
      regulators and in which direction.

Inputs : curated_hits_robust_selectivity.csv, selectivity_scores_v2.csv,
         interaction_stim8hr_30cyto.parquet
Outputs: figures/fig_focused_regulator_landscape.png
         figures/fig_30cytokine_panels.png
"""
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

ROBUST = "curated_hits_robust_selectivity.csv"
SELECTIVITY = "selectivity_scores_v2.csv"
INTERACTION = "interaction_stim8hr_30cyto.parquet"
SUPPRESSOR_COLOR, INDUCER_COLOR = "#c1272d", "#2166ac"

mpl.rcParams.update({"font.family": "Arial", "axes.titleweight": "normal",
                     "axes.labelweight": "normal"})


def landscape(robust):
    """Single-panel selectivity vs effect-size landscape, one point per gene."""
    g = robust.sort_values("selectivity_index", ascending=False).drop_duplicates("gene")
    fig, ax = plt.subplots(figsize=(8.6, 6.4))
    color = np.where(g["lfc_stim"] < 0, SUPPRESSOR_COLOR, INDUCER_COLOR)
    ax.axhline(0, color="0.4", lw=1.0, zorder=1)
    ax.scatter(g["selectivity_index"], g["lfc_stim"], s=30, c=color, alpha=0.8,
               edgecolors="white", linewidths=0.5, zorder=3)

    # Manual label offsets for the crowded high-selectivity cluster
    offsets = {"CCDC134": (4, -9), "NFKB1": (4, 4), "TXNDC12": (4, -9),
               "SLC7A6": (4, 5), "H2AZ1": (4, 3), "ANKRD11": (4, 3),
               "WRAP53": (4, 3), "AVEN": (4, 3)}
    for _, r in g[g["selectivity_index"] > 0.55].iterrows():
        dx, dy = offsets.get(r["gene"], (4, 3))
        ax.annotate(f"{r['gene']}\u2192{r['cytokine']}",
                    (r["selectivity_index"], r["lfc_stim"]),
                    xytext=(dx, dy), textcoords="offset points",
                    fontsize=6.3, style="italic", zorder=5)

    ax.set_xlabel("Selectivity index  (0.7\u00b7focus + 0.3\u00b7footprint)")
    ax.set_ylabel("log\u2082FC on target cytokine in Stim8hr\n"
                  "(KD effect: \u2212 lowers cytokine, + raises)")
    ax.set_title("Focused, robust context-specific cytokine regulators",
                 loc="left", fontsize=10, weight="bold", pad=16)
    ax.text(0.0, 1.02, "n = 128 genes; cross-donor reproducible + confirmed "
            "on-target knockdown", transform=ax.transAxes, fontsize=7, color="0.35")
    ax.margins(0.04)
    ax.set_xlim(right=0.97)
    fig.tight_layout()
    fig.savefig("figures/fig_focused_regulator_landscape.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def cytokine_panels(robust, sel, inter):
    """Grid of per-cytokine panels: grey = all scored genes, colour = robust hits."""
    lf = (inter[inter["perturbation"].isin(sel["gene"])]
          [["perturbation", "cytokine", "lfc_stim"]]
          .rename(columns={"perturbation": "gene"})
          .merge(sel, on="gene", how="left"))
    all_cyto = sorted(inter["cytokine"].unique())
    hitset = set(zip(robust["gene"], robust["cytokine"]))

    fs = 24
    fig, axes = plt.subplots(6, 5, figsize=(24, 27), sharex=True, sharey=True)
    for ax, cyt in zip(axes.ravel(), all_cyto):
        d = lf[lf["cytokine"] == cyt]
        is_hit = np.array([(g, cyt) in hitset for g in d["gene"]])
        hit, non = d[is_hit], d[~is_hit]
        ax.axhline(0, color="0.75", lw=1.2, zorder=1)
        ax.scatter(non["selectivity_index"], non["lfc_stim"], s=27, c="0.8",
                   edgecolors="none", zorder=2)
        col = np.where(hit["lfc_stim"] < 0, SUPPRESSOR_COLOR, INDUCER_COLOR)
        ax.scatter(hit["selectivity_index"], hit["lfc_stim"], s=48, c=col,
                   edgecolors="white", linewidths=0.8, zorder=4)
        ax.set_title(f"{cyt}  (n={len(hit)})", fontsize=fs)
        ax.set_xlim(0, 1)
        ax.set_ylim(-9, 6)
        ax.set_xlabel("Selectivity index", fontsize=fs)
        ax.set_ylabel("log2FC (Stim8hr)", fontsize=fs)
        ax.tick_params(labelsize=fs, length=6, width=1.2)
        for s in ax.spines.values():
            s.set_linewidth(1.2)
    fig.tight_layout()
    fig.savefig("figures/fig_30cytokine_panels.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    robust = pd.read_csv(ROBUST)
    sel = pd.read_csv(SELECTIVITY)[["gene", "selectivity_index"]]
    inter = pd.read_parquet(INTERACTION)
    landscape(robust)
    cytokine_panels(robust, sel, inter)
    print("saved figures/fig_focused_regulator_landscape.png, figures/fig_30cytokine_panels.png")


if __name__ == "__main__":
    main()
