import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
import numpy as np
import scanpy as sc
from sklearn.decomposition import NMF
from scipy.optimize import linear_sum_assignment
from scipy.stats import pearsonr
from cnmf import cNMF

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Add phase2/analyses to path to import pseudobulk_convert
sys.path.append(str(Path(__file__).resolve().parent.parent.parent / "analyses"))
try:
    from pseudobulk_convert import load_and_prep_data, VAR_MODAL_DICT, MODAL_TYPES_DICT
except ImportError as e:
    logger.error(
        f"Could not import pseudobulk_convert. Make sure it exists in phase2/analyses. Error: {e}"
    )
    sys.exit(1)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Test and compare NMF against cNMF outputs."
    )
    parser.add_argument("--project", type=str, default="aging_phase2")
    parser.add_argument(
        "--work-dir",
        type=str,
        default="/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2",
    )
    parser.add_argument("--modality", type=str, default="rna", choices=["rna", "atac"])
    parser.add_argument(
        "--cell-type", type=str, required=True, help="Specific cell type to process."
    )
    parser.add_argument("--num-highvar-genes", type=int, default=2000)
    parser.add_argument("--use-aaf", action="store_true")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--latent-dir-name", type=str, default="cnmf")
    parser.add_argument("--density-threshold", type=float, default=0.1)
    parser.add_argument("--out-dir", type=str, default="nmf_comparison")
    parser.add_argument(
        "--random-state", type=int, default=42, help="Random state for sklearn NMF."
    )
    return parser.parse_args()


def load_df_from_npz(filename):
    with np.load(filename, allow_pickle=True) as f:
        obj = pd.DataFrame(**f)
    return obj


def auto_select_k(cnmf_dir, run_name):
    stats_file = Path(cnmf_dir) / run_name / f"{run_name}.k_selection_stats.df.npz"
    if not stats_file.exists():
        raise FileNotFoundError(f"k_selection_stats file not found: {stats_file}")

    stats = load_df_from_npz(stats_file)
    k_vals = stats["k"].values
    stability = stats["silhouette"].values
    error = stats["prediction_error"].values

    S_norm = (stability - np.min(stability)) / (
        np.max(stability) - np.min(stability) + 1e-9
    )
    E_norm = (error - np.min(error)) / (np.max(error) - np.min(error) + 1e-9)
    dist = np.sqrt((1 - S_norm) ** 2 + (E_norm - 0) ** 2)
    optimal_idx = np.argmin(dist)
    optimal_k = int(k_vals[optimal_idx])
    return optimal_k, stats_file


def compute_correlations_and_match(nmf_df, cnmf_df):
    # Align rows
    common_idx = nmf_df.index.intersection(cnmf_df.index)
    nmf_df_aligned = nmf_df.loc[common_idx]
    cnmf_df_aligned = cnmf_df.loc[common_idx]

    K_nmf = nmf_df.shape[1]
    K_cnmf = cnmf_df.shape[1]

    if K_nmf != K_cnmf:
        logger.warning(f"K mismatch: NMF={K_nmf}, cNMF={K_cnmf}")

    K = min(K_nmf, K_cnmf)
    corr_matrix = np.zeros((K_nmf, K_cnmf))

    for i in range(K_nmf):
        for j in range(K_cnmf):
            col_i = nmf_df_aligned.iloc[:, i]
            col_j = cnmf_df_aligned.iloc[:, j]
            corr, _ = pearsonr(col_i, col_j)
            corr_matrix[i, j] = corr

    # Match using absolute correlation to handle sign flips
    row_ind, col_ind = linear_sum_assignment(-np.abs(corr_matrix))

    matching = []
    for i in range(len(row_ind)):
        nmf_factor = nmf_df.columns[row_ind[i]]
        cnmf_factor = cnmf_df.columns[col_ind[i]]
        r_val = corr_matrix[row_ind[i], col_ind[i]]
        matching.append(
            {
                "NMF_Factor": nmf_factor,
                "cNMF_Factor": cnmf_factor,
                "Pearson_R": r_val,
                "Abs_Pearson_R": abs(r_val),
            }
        )

    match_df = pd.DataFrame(matching)
    corr_df = pd.DataFrame(corr_matrix, index=nmf_df.columns, columns=cnmf_df.columns)

    return corr_df, match_df


def main():
    args = parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)

    work_dir = Path(args.work_dir)
    quants_dir = work_dir / "quants"
    results_dir = work_dir / "results"
    info_dir = work_dir / "sample_info"

    safe_ct = args.cell_type.replace(" ", "_").replace("/", "-")
    run_name = f"{args.project}_{safe_ct}_{args.modality}"

    cnmf_dir = results_dir / "latents" / args.latent_dir_name
    out_dir = results_dir / "latents" / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Auto-select K
    logger.info(f"Auto-selecting K for {run_name} from cNMF outputs...")
    try:
        selected_k, stats_file = auto_select_k(cnmf_dir, run_name)
        logger.info(
            f"Auto-selected K={selected_k} based on distance to ideal point in {stats_file}"
        )
    except Exception as e:
        logger.error(f"Failed to auto-select K: {e}")
        sys.exit(1)

    # 2. Load and prep data
    logger.info("Loading and preprocessing data...")
    raw_file = quants_dir / f"{args.project}.raw.multivi_prep.h5ad"
    annot_file = quants_dir / f"{args.project}.multivi.annotated.h5ad"

    adata_raw, _ = load_and_prep_data(raw_file, annot_file, args.debug)

    modal_full = None
    for k, v in VAR_MODAL_DICT.items():
        if v == args.modality:
            modal_full = k
            break

    modality_list = MODAL_TYPES_DICT.get(args.modality)
    adata_modal = adata_raw[
        adata_raw.obs.modality.isin(modality_list),
        adata_raw.var.modality == modal_full,
    ].copy()

    adata_ct = adata_modal[adata_modal.obs["cell_label"] == args.cell_type].copy()

    # Valid donors filter
    covars_file = (
        info_dir / f"{args.project}.{safe_ct}.{args.modality}.final_covariates.csv"
    )
    if covars_file.exists():
        logger.info("Filtering by valid donors...")
        covars_df = pd.read_csv(covars_file, index_col=0)
        valid_donors = set(covars_df.index.astype(str))
        mask = adata_ct.obs["sample_id"].astype(str).isin(valid_donors)
        adata_ct = adata_ct[mask].copy()

    # Features filter
    reg_file = results_dir / f"{args.project}.all_celltypes.{args.modality}.wls.age.csv"
    if reg_file.exists():
        logger.info("Filtering by regression features...")
        reg_df = pd.read_csv(reg_file)
        ct_reg_df = reg_df[reg_df["tissue"] == safe_ct]

        valid_features = set(ct_reg_df["feature"].unique())

        if args.use_aaf:
            fdr_reg_file = (
                results_dir
                / f"{args.project}.{args.modality}.all_celltypes.wls_fdr_filtered.age.csv"
            )
            if fdr_reg_file.exists():
                fdr_reg_df = pd.read_csv(fdr_reg_file)
                union_features = set(fdr_reg_df["feature"].unique())
                valid_features = valid_features.intersection(union_features)

        intersecting_features = valid_features.intersection(adata_ct.var_names)
        adata_ct = adata_ct[:, list(intersecting_features)].copy()

    logger.info("Applying standard normalizations and HVG selection...")
    sc.pp.filter_cells(adata_ct, min_genes=200)
    sc.pp.filter_cells(adata_ct, min_counts=200)
    sc.pp.filter_genes(adata_ct, min_cells=3)

    temp_adata = adata_ct.copy()
    sc.pp.normalize_total(temp_adata, target_sum=10000)
    sc.pp.log1p(temp_adata)
    sc.pp.highly_variable_genes(
        temp_adata, n_top_genes=args.num_highvar_genes, flavor="seurat"
    )

    hvg_mask = temp_adata.var["highly_variable"].values
    if hasattr(adata_ct.X, "toarray"):
        hvg_counts = np.array(adata_ct.X[:, hvg_mask].sum(axis=1)).flatten()
    else:
        hvg_counts = np.array(adata_ct.X[:, hvg_mask].sum(axis=1)).flatten()

    cells_to_keep = hvg_counts > 0
    adata_ct = adata_ct[cells_to_keep, :].copy()

    # We will use TPM-like normalized counts for NMF (without log) as is standard
    sc.pp.normalize_total(adata_ct, target_sum=10000)
    custom_hvgs = temp_adata.var_names[hvg_mask]
    adata_hvg = adata_ct[:, custom_hvgs].copy()

    # 3. Run Sklearn NMF
    logger.info(f"Running Sklearn NMF with K={selected_k} and init='nndsvda'...")
    # nndsvda does not use random_state
    nmf = NMF(
        n_components=selected_k,
        init="nndsvda",
        max_iter=500,
    )

    X = adata_hvg.X
    if hasattr(X, "toarray"):
        X = X.toarray()

    W = nmf.fit_transform(X)
    H = nmf.components_
    
    logger.info(f"Sklearn NMF completed after {nmf.n_iter_} iterations.")

    nmf_usage = pd.DataFrame(
        W,
        index=adata_hvg.obs_names,
        columns=[f"NMF_{i + 1}" for i in range(selected_k)],
    )
    nmf_spectra = pd.DataFrame(
        H.T,
        index=adata_hvg.var_names,
        columns=[f"NMF_{i + 1}" for i in range(selected_k)],
    )

    # 4. Load cNMF
    logger.info("Loading cNMF results...")
    cnmf_obj = cNMF(output_dir=str(cnmf_dir), name=run_name)
    try:
        cnmf_usage, cnmf_spectra, cnmf_spectra_tpm, top_genes = cnmf_obj.load_results(
            K=selected_k, density_threshold=args.density_threshold
        )
    except Exception as e:
        logger.error(f"Failed to load cNMF results: {e}")
        sys.exit(1)

    # Append 'cNMF_' prefix for clarity
    cnmf_usage.columns = [f"cNMF_{c}" for c in cnmf_usage.columns]
    cnmf_spectra_tpm.columns = [f"cNMF_{c}" for c in cnmf_spectra_tpm.columns]
    cnmf_spectra_df = cnmf_spectra_tpm

    # 5. Compare Usage (W)
    logger.info("Computing Usage (W) correlations and matching factors...")
    usage_corr, usage_match = compute_correlations_and_match(nmf_usage, cnmf_usage)

    # 6. Compare Spectra (H)
    logger.info("Computing Spectra (H) correlations and matching factors...")
    spectra_corr, spectra_match = compute_correlations_and_match(
        nmf_spectra, cnmf_spectra_df
    )

    # 7. Export Results
    logger.info(f"Saving results to {out_dir}...")

    prefix = f"{run_name}_K{selected_k}"
    usage_match.to_csv(out_dir / f"{prefix}_usage_matching.csv", index=False)
    spectra_match.to_csv(out_dir / f"{prefix}_spectra_matching.csv", index=False)

    usage_corr.to_csv(out_dir / f"{prefix}_usage_correlation.csv")
    spectra_corr.to_csv(out_dir / f"{prefix}_spectra_correlation.csv")

    logger.info("Done.")


if __name__ == "__main__":
    main()
