#!/usr/bin/env python3
"""
Build the static-site data layer for the T-cell cytokine-regulator target browser.

    python build_data.py --source build/source_128.csv

Reads:
  build/config.json          score/flag/link/detail definitions (edit here to add columns)
  <source CSV>               one row per gene: symbol, uniprot, <score keys>, <flag keys>,
                             modality, <detail fields>, transcript_id
  build/structures/<SYM>.pdb AlphaFold model (protein atoms extracted from the fpocket PDB)
  build/pockets/<SYM>.json   {"pockets":[{"rank","druggability","volume","residues":[...]}]}
  build/rna_real/<SYM>.json  {"transcript","mrna_len","sites":[{"pos","start","end","modality","score","seq"}], ...}
  build/cds_cache/<TX>.json  {"transcript","mrna_len","cds_start","cds_end"}  (transcript-relative CDS)
  build/perturb/<SYM>.json   {"cytokines":[{"name","effect","dir","fdr"}], "placeholder":false}
  build/landscape.json       {"points":[{"gene","sel","lfc","cyt","dir"}]}  (shared by all gene pages)

Writes:
  data/index.json            table rows + column/link config for the landing page
  data/genes/<SYM>.json       full per-gene record for the detail page
  data/structures/<SYM>.pdb   local AlphaFold model (copied from build/structures)
  data/landscape.json         shared landscape scatter data

All inputs are local — no network required. UniProt display metadata (protein name,
length, function) is fetched once and cached to build/cache/uniprot_<ACC>.json; a missing
fetch degrades gracefully (falls back to cache, then to the bare symbol).
"""
import argparse, csv, json, os, sys, shutil, urllib.request, urllib.error

HERE = os.path.dirname(os.path.abspath(__file__))
BUILD = os.path.join(HERE, "build")
DATA = os.path.join(HERE, "data")
STRUCT = os.path.join(DATA, "structures")
GENES = os.path.join(DATA, "genes")
PERTURB = os.path.join(BUILD, "perturb")
SRC_STRUCT = os.path.join(BUILD, "structures")
POCKETS = os.path.join(BUILD, "pockets")
RNA_REAL = os.path.join(BUILD, "rna_real")
CDS_CACHE = os.path.join(BUILD, "cds_cache")
UP_CACHE = os.path.join(BUILD, "cache")


def log(*a):
    print("[build]", *a, file=sys.stderr)


def load_local(path):
    return json.load(open(path)) if os.path.exists(path) else None


def fetch_json(url, timeout=30):
    try:
        with urllib.request.urlopen(url, timeout=timeout) as r:
            return json.load(r)
    except Exception as e:
        log("  fetch failed:", url, "->", e)
        return None


def uniprot_record(acc, use_net=True):
    """Return display fields from UniProt; cached to build/cache/uniprot_<ACC>.json.
    Falls back to cache when offline, then to all-None."""
    os.makedirs(UP_CACHE, exist_ok=True)
    cache = os.path.join(UP_CACHE, f"uniprot_{acc}.json")
    out = {"protein_name": None, "length": None, "function": None, "ensembl_gene": None}
    up = fetch_json(f"https://rest.uniprot.org/uniprotkb/{acc}.json") if use_net else None
    if not up:
        cached = load_local(cache)
        if cached:
            log("  uniprot from cache")
            return cached
        return out
    try:
        out["protein_name"] = up["proteinDescription"]["recommendedName"]["fullName"]["value"]
    except Exception:
        pass
    out["length"] = up.get("sequence", {}).get("length")
    for c in up.get("comments", []):
        if c.get("commentType") == "FUNCTION" and c.get("texts"):
            out["function"] = c["texts"][0].get("value")
            break
    for xr in up.get("uniProtKBCrossReferences", []):
        if xr.get("database") == "Ensembl":
            for p in xr.get("properties", []):
                if p.get("key") == "GeneId":
                    out["ensembl_gene"] = p.get("value", "").split(".")[0]
    json.dump(out, open(cache, "w"))
    return out


def rna_record(symbol, transcript_id):
    """Real RNA record: transcript length + real ASO/siRNA sites from build/rna_real,
    merged with transcript-relative CDS coords from build/cds_cache. All local."""
    rec = load_local(os.path.join(RNA_REAL, f"{symbol}.json"))
    if rec is None:
        return None
    tx = rec.get("transcript") or transcript_id
    cds = load_local(os.path.join(CDS_CACHE, f"{tx}.json")) or {}
    rec["cds_start"] = cds.get("cds_start")
    rec["cds_end"] = cds.get("cds_end")
    # prefer the Ensembl mRNA length when available (authoritative full transcript)
    if cds.get("mrna_len"):
        rec["mrna_len"] = cds["mrna_len"]
    rec.setdefault("sites_placeholder", False)
    return rec


def copy_structure(symbol):
    """Ensure data/structures/<SYM>.pdb exists. Prefer build/structures as source;
    if that (large, dev-only) folder isn't shipped, the already-built data/ copy is
    the source of truth. Returns relative path or None."""
    src = os.path.join(SRC_STRUCT, f"{symbol}.pdb")
    dst = os.path.join(STRUCT, f"{symbol}.pdb")
    rel = f"data/structures/{symbol}.pdb"
    if os.path.exists(src):
        shutil.copyfile(src, dst)
        return rel
    if os.path.exists(dst):        # data/ is self-sustaining without the build/ duplicate
        return rel
    return None


def pocket_record(symbol):
    """Top-ranked druggable pockets with proximal residue lists (local)."""
    rec = load_local(os.path.join(POCKETS, f"{symbol}.json"))
    return rec["pockets"] if rec else []


def build_links(cfg, acc, symbol, ensembl_gene):
    links = []
    for L in cfg["external_links"]:
        try:
            url = L["url"].format(acc=acc or "", symbol=symbol,
                                  ensembl_gene=ensembl_gene or "")
        except Exception:
            continue
        # skip links whose required field is missing
        if "{ensembl_gene}" in L["url"] and not ensembl_gene:
            continue
        if "{acc}" in L["url"] and not acc:
            continue
        links.append({"key": L["key"], "label": L["label"], "url": url})
    return links


def render_writeup():
    """Convert build/content/research_writeup.md → data/writeup.json {title, html}.
    Minimal, dependency-free: handles #/##/### headings, **bold**, *italic*, and
    paragraphs (the constructs this document uses). No external markdown lib —
    keeps the shipped package small and the build offline."""
    import re, html as _html
    path = os.path.join(BUILD, "content", "research_writeup.md")
    if not os.path.exists(path):
        return None
    lines = open(path).read().splitlines()
    title, blocks, para = None, [], []

    def inline(t):
        t = _html.escape(t)
        t = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", t)
        t = re.sub(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)", r"<em>\1</em>", t)
        return t

    def flush():
        if para:
            blocks.append("<p>" + inline(" ".join(para)) + "</p>")
            para.clear()

    for ln in lines:
        s = ln.strip()
        if not s:
            flush(); continue
        m = re.match(r"^(#{1,3})\s+(.*)$", s)
        if m:
            flush()
            level = len(m.group(1)); txt = m.group(2)
            if level == 1 and title is None:
                title = txt
            else:
                blocks.append(f"<h{level}>{inline(txt)}</h{level}>")
        else:
            para.append(s)
    flush()
    rec = {"title": title or "Research write-up", "html": "\n".join(blocks)}
    with open(os.path.join(DATA, "writeup.json"), "w") as f:
        json.dump(rec, f)
    return rec


def load_perturb(symbol):
    """Real per-gene Perturb-seq effect payload: the cytokines this gene significantly
    regulates (robust hits only), each with log2FC in Stim8hr and direction."""
    rec = load_local(os.path.join(PERTURB, f"{symbol}.json"))
    return rec if rec else {"cytokines": [], "placeholder": False}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=os.path.join(BUILD, "source_128.csv"))
    ap.add_argument("--no-net", action="store_true", help="skip UniProt metadata fetch (uses cache)")
    args = ap.parse_args()

    with open(os.path.join(BUILD, "config.json")) as f:
        cfg = json.load(f)
    os.makedirs(STRUCT, exist_ok=True)
    os.makedirs(GENES, exist_ok=True)

    with open(args.source) as f:
        rows = list(csv.DictReader(f))
    log(f"{len(rows)} genes from {os.path.basename(args.source)}")

    score_keys = [s["key"] for s in cfg["scores"]]
    flag_keys = [fl["key"] for fl in cfg["flags"]]
    detail_keys = [d["key"] for d in cfg.get("detail_fields", [])]
    index_rows = []
    n_struct = n_rna = n_pockets = n_perturb = 0

    for row in rows:
        symbol = row["symbol"].strip()
        acc = row.get("uniprot", "").strip() or None
        tx_id = row.get("transcript_id", "").strip() or None
        up = uniprot_record(acc, use_net=not args.no_net) if acc else \
            {"protein_name": None, "length": None, "function": None, "ensembl_gene": None}

        struct_rel = copy_structure(symbol)
        pockets = pocket_record(symbol)
        rna = rna_record(symbol, tx_id)
        perturb = load_perturb(symbol)
        scores = {k: _num(row.get(k)) for k in score_keys}
        flags = {k: _bool(row.get(k)) for k in flag_keys}
        details = {k: row.get(k, "").strip() for k in detail_keys}
        links = build_links(cfg, acc, symbol, up["ensembl_gene"])

        n_struct += bool(struct_rel); n_pockets += bool(pockets)
        n_rna += bool(rna and rna.get("sites")); n_perturb += bool(perturb.get("cytokines"))

        record = {
            "symbol": symbol,
            "uniprot": acc,
            "protein_name": up["protein_name"] or symbol,
            "length": up["length"],
            "function": up["function"],
            "ensembl_gene": up["ensembl_gene"],
            "scores": scores,
            "flags": flags,
            "details": details,
            "modality": row.get("modality", "Undetermined").strip() or "Undetermined",
            "structure": struct_rel,
            "pockets": pockets,
            "links": links,
            "rna": rna,
            "perturb": perturb,
        }
        with open(os.path.join(GENES, f"{symbol}.json"), "w") as f:
            json.dump(record, f)

        index_rows.append({
            "symbol": symbol,
            "protein_name": record["protein_name"],
            "modality": record["modality"],
            "scores": scores,
            "flags": flags,
            "has_structure": bool(struct_rel),
        })

    # shared landscape scatter (all 128 genes) — copied into data/ for the gene pages
    land = load_local(os.path.join(BUILD, "landscape.json"))
    if land:
        with open(os.path.join(DATA, "landscape.json"), "w") as f:
            json.dump(land, f)

    writeup = render_writeup()

    # per-cytokine showcase panels (effect vs selectivity, one panel per cytokine)
    panels = load_local(os.path.join(BUILD, "cytokine_panels.json"))
    if panels:
        with open(os.path.join(DATA, "cytokine_panels.json"), "w") as f:
            json.dump(panels, f)

    index = {
        "site": cfg["site"],
        "scores": cfg["scores"],
        "flags": cfg["flags"],
        "detail_fields": cfg.get("detail_fields", []),
        "modalities": cfg["modalities"],
        "genes": index_rows,
        "n_genes": len(index_rows),
    }
    with open(os.path.join(DATA, "index.json"), "w") as f:
        json.dump(index, f)
    log(f"wrote index.json ({len(index_rows)} genes) + per-gene records")
    log(f"coverage: structures {n_struct}/{len(rows)} · pockets {n_pockets} · "
        f"RNA-sites {n_rna} · cytokine-effects {n_perturb} · landscape {'yes' if land else 'no'} · "
        f"writeup {'yes' if writeup else 'no'} · cytokine-panels {'yes' if panels else 'no'}")


def _num(v):
    try:
        return round(float(v), 4)
    except (TypeError, ValueError):
        return None


def _bool(v):
    return str(v).strip() in ("1", "true", "True", "yes", "Y")


if __name__ == "__main__":
    main()
