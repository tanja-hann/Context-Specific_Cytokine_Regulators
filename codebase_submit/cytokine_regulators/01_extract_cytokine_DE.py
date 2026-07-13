"""
Step 1 — Extract cytokine differential-expression effect sizes.

Streams the genome-scale CD4+ T-cell Perturb-seq DE statistics
(Zhu, Dann et al. 2025; DESeq2 fit per perturbation x culture condition)
directly from its public S3 bucket, keeping only the cytokine-module columns.
The full 16.8 GB matrix is never written to disk: the file is read in
row-blocks and sliced to the cytokine genes on the fly.

Two products are written:
  cytokine_DE_long.parquet  long-format effect sizes (log_fc, lfcSE,
                            adj_p_value, baseMean, zscore) for every
                            perturbation x cytokine x culture condition
  DE_obs_metadata.parquet   per-perturbation QC/summary metadata precomputed
                            genome-wide by the dataset authors, including
                            n_downstream (genes changed across all 10,282
                            measured genes), on-target knockdown, and
                            cross-donor reproducibility

The cytokine module is the "Cytokine" category of immune_effector_genes.csv
(from the dataset's companion analysis repo).

Inputs : immune_effector_genes.csv
Outputs: cytokine_DE_long.parquet, DE_obs_metadata.parquet
"""
import h5py
import s3fs
import numpy as np
import pandas as pd

IMMUNE_GENES_CSV = "immune_effector_genes.csv"
S3_KEY = "genome-scale-tcell-perturb-seq/marson2025_data/GWCD4i.DE_stats.h5ad"
EFFECT_LAYERS = ["log_fc", "lfcSE", "adj_p_value", "baseMean"]
ROW_BLOCK = 3000

# Per-perturbation metadata columns to carry through from .obs
OBS_CATEGORICAL = ["target_contrast_gene_name", "culture_condition",
                   "target_contrast", "n_total_genes_category"]
OBS_NUMERIC = ["n_cells_target", "ontarget_effect_size", "ontarget_significant",
               "target_baseMean", "n_guides", "n_downstream",
               "donor_correlation_all_mean", "donor_correlation_hits_mean",
               "low_target_gex", "neighboring_gene_KD", "distal_offtarget_flag",
               "single_guide_estimate", "n_total_de_genes"]


def read_categorical(group):
    cats = group["categories"][:].astype(str)
    codes = group["codes"][:]
    return pd.Categorical.from_codes(codes, cats)


def load_obs(h):
    obs = pd.DataFrame(index=h["obs"]["index"][:].astype(str))
    for col in OBS_CATEGORICAL:
        obs[col] = read_categorical(h["obs"][col]).astype(str)
    for col in OBS_NUMERIC:
        obs[col] = h["obs"][col][:]
    return obs


def main():
    immune = pd.read_csv(IMMUNE_GENES_CSV)
    cytokines = immune.loc[immune["Category"] == "Cytokine", "gene_name"].tolist()

    fs = s3fs.S3FileSystem(anon=True)
    with fs.open(S3_KEY, "rb", block_size=2 ** 20) as f:
        h = h5py.File(f, "r")

        obs = load_obs(h)
        var = pd.DataFrame(
            {"gene_ids": h["var"]["gene_ids"][:].astype(str),
             "gene_name": h["var"]["gene_name"][:].astype(str)},
            index=h["var"]["_index"][:].astype(str))

        # Column indices of the cytokines present in the measured gene set
        present = [g for g in cytokines if g in set(var["gene_name"])]
        name_to_idx = {n: i for i, n in enumerate(var["gene_name"].values)}
        cyto_idx = np.array([name_to_idx[g] for g in present])

        # Stream each effect-size layer in row-blocks, keeping cytokine columns
        n_rows = h["layers"]["log_fc"].shape[0]
        layer_data = {L: np.empty((n_rows, len(cyto_idx)), dtype=np.float64)
                      for L in EFFECT_LAYERS}
        for L in EFFECT_LAYERS:
            ds = h["layers"][L]
            for r0 in range(0, n_rows, ROW_BLOCK):
                r1 = min(r0 + ROW_BLOCK, n_rows)
                layer_data[L][r0:r1, :] = ds[r0:r1, :][:, cyto_idx]

    # Reshape (perturbation x cytokine) grid to long format
    n_cyto = len(present)
    rep = np.repeat(np.arange(n_rows), n_cyto)
    tile = np.tile(np.arange(n_cyto), n_rows)
    long = pd.DataFrame({
        "perturbation": obs["target_contrast_gene_name"].values[rep],
        "culture_condition": obs["culture_condition"].values[rep],
        "target_contrast": obs["target_contrast"].values[rep],
        "cytokine": np.array(present)[tile],
        "log_fc": layer_data["log_fc"].reshape(-1),
        "lfcSE": layer_data["lfcSE"].reshape(-1),
        "adj_p_value": layer_data["adj_p_value"].reshape(-1),
        "baseMean": layer_data["baseMean"].reshape(-1),
        "n_cells_target": obs["n_cells_target"].values[rep],
        "ontarget_effect_size": obs["ontarget_effect_size"].values[rep],
        "ontarget_significant": obs["ontarget_significant"].values[rep],
        "donor_corr_hits_mean": obs["donor_correlation_hits_mean"].values[rep],
        "donor_corr_all_mean": obs["donor_correlation_all_mean"].values[rep],
        "low_target_gex": obs["low_target_gex"].values[rep],
        "neighboring_gene_KD": obs["neighboring_gene_KD"].values[rep],
        "n_guides": obs["n_guides"].values[rep],
    })
    long["zscore"] = long["log_fc"] / long["lfcSE"]
    long.to_parquet("cytokine_DE_long.parquet", index=False)

    obs.reset_index().rename(columns={"index": "obs_name"}).to_parquet(
        "DE_obs_metadata.parquet", index=False)

    print(f"cytokines present: {len(present)}")
    print(f"cytokine_DE_long.parquet rows: {len(long)}")
    print(f"DE_obs_metadata.parquet rows: {len(obs)}")


if __name__ == "__main__":
    main()
