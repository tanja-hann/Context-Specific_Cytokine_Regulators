# Context-Specific Cytokine Regulators

An interactive target browser for **128 context-specific cytokine regulators** in
human CD4+ T cells — genes that control disease-relevant cytokines specifically
when T cells are activated, and not at rest. These are candidate drug targets
whose interference could calm an over-active immune response while leaving the
resting immune system largely untouched.

The browser presents, for each target: a **selectivity** score, **two-axis
druggability** (protein-level vs RNA-level), a **nominated therapeutic modality**
(small molecule, antibody/biologic, ASO, or siRNA), the **AlphaFold structure**
with its predicted druggable pockets highlighted, real **ASO/siRNA target sites**
mapped onto the transcript, the **cytokines it significantly regulates**, and
links out to UniProt, AlphaFold, Open Targets, and other resources.

## The science, briefly

Immune diseases such as rheumatoid arthritis, psoriasis, and inflammatory bowel
disease are driven by T cells releasing too much of certain cytokines. Building
on a genome-scale **Perturb-seq** screen in primary human CD4+ T cells (Zhu,
Dann et al. 2025), this project asks which genes regulate a disease-relevant
cytokine **only in the stimulated state** — using a statistical
*interaction* test between resting and stimulated conditions — and then, for the
128 regulators that pass, which could realistically be turned into a drug and by
what modality. The full analysis and write-up are on the site's **Background**
and **Data support** pages.

## Viewing the site

It's a fully static website (plain HTML/CSS/JS + JSON) — no backend, no build
step needed to view it. It is served here via GitHub Pages.

To run it locally, from the repository root:

```bash
python3 -m http.server 8000
# then open http://localhost:8000
```

Serve over HTTP (not by double-clicking the files) — the pages fetch JSON. The
3D viewer (Mol*), fonts, and all AlphaFold structures are bundled in the repo,
so it works without an internet connection.

## What's in the site

| Page | Contents |
|---|---|
| **Home** (`index.html`) | Interactive protein-vs-RNA druggability matrix + a sortable, filterable table of all 128 targets |
| **Gene pages** (`gene.html`) | Per-target scores, 3D AlphaFold structure with druggable pockets, mRNA map with ASO/siRNA sites, landscape position, and cytokine effects |
| **Data support** (`data-support.html`) | The 30-cytokine effect-vs-selectivity showcase |
| **Background** (`background.html`) | Research write-up |

## Data provenance

Druggability scores are cohort-relative values (protein and RNA axes each
renormalized across the 128-gene set). Pocket-lining residues are computed from
fpocket alpha-spheres on the AlphaFold models. RNA sites are qualifying ASO/siRNA
target windows (accessibility + efficacy + specificity gates). Cytokine effects
are the robust per-gene × cytokine log₂ fold-changes (8 h stimulation) at
interaction FDR < 0.05.

The `build/` folder and `build_data.py` regenerate the site's `data/` files from
the source tables; the front-end reads only `data/`, so refreshing results never
requires touching the site code.

## Citation

Zhu R., Dann E. et al. (2025) *Genome-scale perturb-seq in primary human CD4+ T
cells maps context-specific regulators of T cell programs and human immune
traits.* bioRxiv.

## Credits

Developed by Tanja Hann, with Claude Science.
