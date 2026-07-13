# T-cell Cytokine Regulator Target Browser

A lightweight, fully static website for browsing **128 context-specific cytokine
regulators** identified in genome-scale CD4+ T-cell Perturb-seq (Zhu, Dann et
al. 2025). Each target carries selectivity and two-axis (protein / RNA)
druggability scores, a nominated therapeutic modality, an AlphaFold structure
with its fpocket-detected druggable pockets highlighted, real ASO/siRNA target
sites on the transcript, external database links, and the cytokines it
significantly regulates.

**No backend. No build server.** Everything is plain HTML/CSS/JS + JSON files,
so it hosts directly on GitHub Pages.

The landing page opens with an **interactive two-axis druggability matrix**
(protein-level vs RNA-level relative druggability, colored by nominated
modality); click any point — or any table row — to open that gene's page.

---

## Quick start (local preview)

```bash
cd site
python -m http.server 8000
# open http://localhost:8000
```

> You must serve over HTTP (not `file://`) — the pages `fetch()` JSON.

**Runs offline.** Everything this site loads is vendored locally — Mol* (3D
viewer), Roboto / Roboto Mono fonts, and the AlphaFold structures — so the demo
runs with no internet. (Vendored assets: `assets/vendor/molstar/`,
`assets/fonts/`, `data/structures/`.)

Note: the bundled `molstar.js` library contains inert external URLs (license
comments and its built-in default data-provider endpoints, e.g. RCSB/PDBe).
These are **not** contacted by this site — structures are loaded from local
`data/structures/*.pdb` paths, so no remote provider is invoked. If you later
add a feature that fetches structures by PDB ID from a remote source, that
specific action would require network access.

---

## How the data flows

```
build/config.json         ← score columns, flags, links, detail fields (EDIT HERE)
build/source_128.csv      ← one row per gene: symbol, uniprot, scores, flags, modality, details, transcript
build/structures/<SYM>.pdb← AlphaFold model; occupancy col = 1.0 for pocket-lining residues
build/pockets/<SYM>.json  ← top-3 fpocket pockets (druggability, volume, proximal residues)
build/rna_real/<SYM>.json ← real ASO/siRNA target sites (position, window, accessibility, modality)
build/cds_cache/<TX>.json ← transcript-relative CDS bounds (Ensembl, cached)
build/perturb/<SYM>.json  ← cytokines this gene significantly regulates (log2FC, direction)
build/landscape.json      ← shared effect-vs-selectivity scatter (all 128 genes)
        │
        ▼   python build_data.py --source build/source_128.csv
data/index.json           → landing-page table + interactive matrix
data/genes/<SYM>.json     → per-gene detail page
data/structures/<SYM>.pdb → AlphaFold model (copied from build/, pocket residues in occupancy col)
data/landscape.json       → landscape scatter shared by all gene pages
```

The frontend reads only the `data/` files. **Swapping in new results never
touches the site code** — you regenerate `data/` and push. All inputs are local;
`build_data.py` needs the network only to refresh UniProt display metadata
(protein name / length / function), which is cached under `build/cache/` — pass
`--no-net` to build entirely offline from the cache.

### Data provenance

Scores are the **final, cohort-relative** values from the druggability analysis
(`druggability_summary_table.csv`): `protein_score_relative` and
`RNA_score_relative` (each renormalized across the 128-gene cohort, centered
0.5), plus the per-gene selectivity index. Pocket-lining residues are computed
from fpocket alpha-spheres (≤4.5 Å to any sphere) on the AlphaFold model and
encoded in the PDB occupancy column so Mol*'s built-in `occupancy` color theme
highlights them. RNA sites are the qualifying ASO/siRNA target windows
(accessibility + efficacy + specificity gates). Cytokine effects are the robust
per-gene×cytokine log2FC (Stim8hr) at interaction FDR < 0.05.

---

## Deploy to GitHub Pages

```bash
git init && git add . && git commit -m "target browser"
git branch -M main
git remote add origin git@github.com:<you>/<repo>.git
git push -u origin main
```
Then in the repo: **Settings → Pages → Source: `main` / root** (or move the
`site/` contents to `/docs` and point Pages at `/docs`). Live in ~1 minute.

---

## Files

| Path | Purpose |
|---|---|
| `index.html` / `assets/js/table.js` | Sortable, filterable score table with color bars |
| `index.html` / `assets/js/table.js` | Interactive druggability matrix + sortable, filterable score table |
| `gene.html` / `assets/js/gene.js` | Per-gene page: scores, Mol* 3D viewer with pockets, mRNA map, landscape, cytokine effects |
| `assets/js/theme.js` | Light/dark theme toggle (light default) + footer |
| `assets/css/style.css` | Warm Anthropic-inspired theme (light + dark) |
| `build_data.py` | Data pipeline (source → `data/`) |
| `build/config.json` | Column / flag / link / detail-field definitions |

Current data: the **final 128-gene target set** (128 context-specific cytokine
regulators, cohort-relative druggability scores, real modality calls,
structures, RNA sites, and cytokine effects).
