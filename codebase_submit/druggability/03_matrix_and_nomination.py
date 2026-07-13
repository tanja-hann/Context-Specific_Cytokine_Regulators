"""
Step 3 — Two-axis druggability matrix and modality nomination.

Merges the protein and RNA scores for the 128 regulators, places each gene in
the protein x RNA druggability plane (relative axes), assigns it to a quadrant,
and nominates a therapeutic modality.

Quadrants split each relative axis at 0.30 (below 0.30 on an axis = that
modality class is not a viable route):
    Both viable        protein >= 0.3 and RNA >= 0.3
    Protein-preferred  protein >= 0.3 and RNA <  0.3
    RNA-preferred      protein <  0.3 and RNA >= 0.3
    Neither viable     protein <  0.3 and RNA <  0.3

Modality nomination — the higher relative axis wins, and the winning axis's
sub-route is taken (protein -> small molecule / antibody-biologic; RNA -> siRNA / ASO).

Outputs the full per-gene table and the matrix figure.

Inputs : protein_scores.tsv, rna_scores.tsv
Outputs: druggability_summary_table.csv
         figures/druggability_quadrants_128.png
"""
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from adjustText import adjust_text

PROTEIN = "protein_scores.tsv"
RNA = "rna_scores.tsv"
THRESHOLD = 0.30
MODALITY_COLOR = {"Small molecule": "#1F5FA8", "Antibody/biologic": "#3EA9BF",
                  "ASO": "#E87DB4", "siRNA": "#B22222"}
NONE_COLOR = "#FFFFFF"
PROTEIN_ROUTE = {"small molecule": "Small molecule",
                 "antibody/biologic": "Antibody/biologic", "none": "none"}
FS, LW = 9, 1.0


def build_table(prot, rna):
    m = prot[["gene", "uniprot", "protein_score", "protein_score_abs", "best_route"]].merge(
        rna[["gene", "transcript_id", "siRNA_score", "ASO_score", "RNA_score",
             "RNA_score_relative", "RNA_modality"]], on="gene", how="inner")
    assert len(m) == 128, f"expected 128 genes, got {len(m)}"
    m = m.rename(columns={"protein_score": "P", "RNA_score_relative": "R"})

    m["nominated_modality"] = m.apply(
        lambda r: PROTEIN_ROUTE.get(r.best_route, "Small molecule") if r.P >= r.R
        else r.RNA_modality, axis=1)
    m["quadrant"] = m.apply(
        lambda r: ("Both viable" if r.P >= THRESHOLD and r.R >= THRESHOLD
                   else "Protein-preferred" if r.P >= THRESHOLD
                   else "RNA-preferred" if r.R >= THRESHOLD
                   else "Neither viable"), axis=1)
    m["point_hex"] = [NONE_COLOR if (q == "Neither viable" or md == "none")
                      else MODALITY_COLOR[md]
                      for q, md in zip(m.quadrant, m.nominated_modality)]
    m["winning_axis"] = np.where(m.P >= m.R, "protein", "RNA")
    m["axis_margin"] = (m.P - m.R).round(3)
    return m


def figure(m):
    mpl.rcParams.update({
        "font.family": "Arial", "font.size": FS, "axes.linewidth": LW,
        "axes.edgecolor": "black", "axes.labelcolor": "black", "text.color": "black",
        "xtick.color": "black", "ytick.color": "black", "legend.frameon": False,
    })
    fig, ax = plt.subplots(figsize=(8.2, 7.4))
    ax.axvline(THRESHOLD, color="black", lw=LW, zorder=1)
    ax.axhline(THRESHOLD, color="black", lw=LW, zorder=1)
    ax.scatter(m.P, m.R, s=46, c=m.point_hex, marker="o", edgecolors="black",
               linewidths=LW, zorder=3)

    for x, y, ha, va, txt in [(0.015, 0.955, "left", "top", "RNA-PREFERRED"),
                              (0.985, 0.955, "right", "top", "BOTH VIABLE"),
                              (0.0, -0.105, "left", "top", "NEITHER VIABLE"),
                              (1.0, -0.105, "right", "top", "PROTEIN-PREFERRED")]:
        ax.text(x, y, txt, transform=ax.transAxes, ha=ha, va=va, fontweight="bold")

    # Label all RNA-preferred genes plus the top 3 on each axis
    lab = set(m[m.quadrant == "RNA-preferred"].gene)
    lab |= set(m.nlargest(3, "P").gene) | set(m.nlargest(3, "R").gene)
    sub = m[m.gene.isin(lab)]
    texts = [ax.text(r.P, r.R, r.gene, fontsize=FS, fontstyle="italic") for _, r in sub.iterrows()]
    adjust_text(texts, x=list(sub.P), y=list(sub.R), ax=ax, expand=(1.15, 1.4),
                arrowprops=dict(arrowstyle="-", color="black", lw=LW))

    ax.set_xlabel("Protein-level druggability  (relative)")
    ax.set_ylabel("RNA-level druggability  (relative)")
    ax.set_title("Two-axis druggability of 128 context-specific cytokine regulators", loc="left")
    ax.set_xlim(-0.03, 1.03)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xticks(np.arange(0, 1.01, 0.2))
    ax.set_yticks(np.arange(0, 1.01, 0.2))

    handles = [Line2D([], [], marker="o", ls="", mfc=c, mec="black", mew=LW, ms=8, label=k)
               for k, c in MODALITY_COLOR.items()]
    handles.append(Line2D([], [], marker="o", ls="", mfc=NONE_COLOR, mec="black",
                          mew=LW, ms=8, label="Neither viable"))
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0),
              fontsize=FS, title="Top modality", title_fontsize=FS, handletextpad=0.4)
    fig.savefig("figures/druggability_quadrants_128.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main():
    prot = pd.read_csv(PROTEIN, sep="\t")
    rna = pd.read_csv(RNA, sep="\t")
    m = build_table(prot, rna)

    labeled = set(m[m.quadrant == "RNA-preferred"].gene) | \
        set(m.nlargest(3, "P").gene) | set(m.nlargest(3, "R").gene)
    order = {"RNA-preferred": 0, "Both viable": 1, "Protein-preferred": 2, "Neither viable": 3}
    m["_q"] = m.quadrant.map(order)
    m["winning_score"] = np.where(m.P >= m.R, m.P, m.R)
    tab = m.sort_values(["_q", "winning_score"], ascending=[True, False]).drop(columns="_q")

    summary = tab[["gene", "uniprot", "transcript_id", "protein_score_abs", "P",
                   "RNA_score", "R", "best_route", "siRNA_score", "ASO_score",
                   "RNA_modality", "nominated_modality", "winning_axis", "axis_margin",
                   "quadrant", "point_hex"]].rename(
        columns={"P": "protein_score_relative", "R": "RNA_score_relative"})
    summary["labeled_in_figure"] = summary.gene.isin(labeled)
    summary.to_csv("druggability_summary_table.csv", index=False)

    figure(m)
    print(f"quadrants: {m['quadrant'].value_counts().to_dict()}")
    print(f"nominated: {m['nominated_modality'].value_counts().to_dict()}")


if __name__ == "__main__":
    main()
