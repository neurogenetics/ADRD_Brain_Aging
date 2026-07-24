import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
from cnmf import cNMF
import scanpy as sc


# Import functions from pseudobulk_convert to reuse loading logic
from pseudobulk_convert import load_and_prep_data, VAR_MODAL_DICT, MODAL_TYPES_DICT

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run cNMF on raw single-cell data per cell type."
    )
    parser.add_argument("--project", type=str, default=DEFAULT_PROJECT)
    parser.add_argument("--work-dir", type=str, default=DEFAULT_WRK_DIR)
    parser.add_argument("--modality", type=str, default="rna", choices=["rna", "atac"])
    parser.add_argument(
        "--components",
        type=int,
        nargs="+",
        default=None,
        help="List of components (K) to evaluate. If not provided, generates a range from 2 to (donors // 2).",
    )
    parser.add_argument("--n-iter", type=int, default=100)
    parser.add_argument("--num-highvar-genes", type=int, default=2000)
    parser.add_argument(
        "--covariates",
        type=str,
        nargs="+",
        default=None,
        help="List of technical covariates in obs to regress out using cNMF Preprocess (Harmony).",
    )
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--debug", action="store_true")
    parser.add_argument(
        "--cell-type", type=str, required=True, help="Specific cell type to process."
    )
    parser.add_argument(
        "--cnmf-dir-name",
        type=str,
        default="cnmf",
        help="Name of the cnmf output directory within the latents path.",
    )
    parser.add_argument(
        "--use-aaf",
        action="store_true",
        help="Use the union of age-associated features across all cell types that were tested in the target cell type.",
    )
    parser.add_argument(
        "--exclude-ids",
        type=str,
        default="",
        help="Comma separated list of sample IDs to exclude.",
    )
    return parser.parse_args()


def process_cell_type(
    adata_ct, cell_type, args, cnmf_dir, tmp_dir, results_dir, info_dir
):
    logger.info(f"--- Processing Cell Type: {cell_type} ---")
    if adata_ct.n_obs < 50:
        logger.warning(f"Skipping {cell_type}: only {adata_ct.n_obs} cells.")
        return

    safe_ct = cell_type.replace(" ", "_").replace("/", "-")
    run_name = f"{args.project}_{safe_ct}_{args.modality}"
    out_base = str(tmp_dir / run_name)

    genes_file = None
    tpm_fn = None

    # Filter cells based on valid donors from covariates
    covars_file = (
        info_dir / f"{args.project}.{safe_ct}.{args.modality}.final_covariates.csv"
    )
    if not covars_file.exists():
        logger.error(f"Covariates file not found: {covars_file}")
        sys.exit(1)

    logger.info(f"Loading valid donors from covariates file: {covars_file}...")
    try:
        # Load the index (first column) which should contain the donor/sample IDs
        covars_df = pd.read_csv(covars_file, index_col=0)
        valid_donors = set(covars_df.index.astype(str))

        # Check against 'sample_id' in obs
        if "sample_id" not in adata_ct.obs.columns:
            logger.error("'sample_id' column not found in AnnData obs.")
            sys.exit(1)

        # Filter the AnnData object
        mask = adata_ct.obs["sample_id"].astype(str).isin(valid_donors)
        logger.info(
            f"Filtering cells by valid donors. Keeping {mask.sum()} of {len(mask)} cells."
        )
        adata_ct = adata_ct[mask].copy()

        unique_donors = adata_ct.obs["sample_id"].nunique()
        logger.info(f"Analysis will proceed with {unique_donors} unique donors.")

        components_to_use = args.components
        if components_to_use is None:
            max_k = max(2, int(unique_donors / 2))
            components_to_use = list(range(2, max_k + 1))
            logger.info(
                f"Dynamically generated components range to test: {components_to_use}"
            )

        if adata_ct.n_obs < 50:
            logger.warning(
                f"Skipping {cell_type}: only {adata_ct.n_obs} cells remaining after donor filter."
            )
            return

    except Exception as e:
        logger.error(f"Failed to filter cells by valid donors: {e}")
        sys.exit(1)

    # Filter features based on age association results (wls)
    reg_file = results_dir / f"{args.project}.all_celltypes.{args.modality}.wls.age.csv"
    if not reg_file.exists():
        logger.error(f"Regression results file not found: {reg_file}")
        sys.exit(1)

    logger.info(f"Loading regression results from {reg_file} to filter features...")
    try:
        reg_df = pd.read_csv(reg_file)

        if "tissue" not in reg_df.columns or "feature" not in reg_df.columns:
            logger.error(
                "Regression results must contain 'tissue' and 'feature' columns."
            )
            sys.exit(1)

        # Filter for the specific cell type (spaces replaced with underscores)
        safe_ct_tissue = cell_type.replace(" ", "_")
        ct_reg_df = reg_df[reg_df["tissue"] == safe_ct_tissue]

        if len(ct_reg_df) == 0:
            logger.error(
                f"No regression results found for cell type '{safe_ct_tissue}' in {reg_file}."
            )
            sys.exit(1)

        valid_features = set(ct_reg_df["feature"].unique())
        logger.info(
            f"Found {len(valid_features)} tested features for {safe_ct_tissue}."
        )

        if args.use_aaf:
            fdr_reg_file = (
                results_dir
                / f"{args.project}.{args.modality}.all_celltypes.wls_fdr_filtered.age.csv"
            )
            if not fdr_reg_file.exists():
                logger.error(f"FDR Regression results file not found: {fdr_reg_file}")
                sys.exit(1)

            logger.info(
                f"Loading FDR regression results from {fdr_reg_file} to construct union of age features..."
            )
            fdr_reg_df = pd.read_csv(fdr_reg_file)
            if "feature" not in fdr_reg_df.columns:
                logger.error("FDR Regression results must contain 'feature' column.")
                sys.exit(1)

            union_features = set(fdr_reg_df["feature"].unique())
            logger.info(
                f"Found {len(union_features)} unique age-associated features across all cell types."
            )

            valid_features = valid_features.intersection(union_features)
            logger.info(
                f"Filtered to {len(valid_features)} tested features that are in the age-associated union set."
            )

        # Intersect with the var_names in adata
        intersecting_features = valid_features.intersection(adata_ct.var_names)
        logger.info(
            f"Filtering AnnData to {len(intersecting_features)} overlapping features."
        )

        if len(intersecting_features) == 0:
            logger.error(
                "No overlapping features found between regression results and AnnData."
            )
            sys.exit(1)

        adata_ct = adata_ct[:, list(intersecting_features)].copy()

    except Exception as e:
        logger.error(f"Failed to filter features based on regression results: {e}")
        sys.exit(1)

    if args.covariates:
        # Use cNMF Preprocess module to integrate known technical covariates via Harmony
        from cnmf.preprocess import Preprocess
        import cnmf.preprocess
        import numpy as np

        # --- Monkeypatch cNMF's buggy moe_correct_ridge ---
        def patched_moe_correct_ridge(
            Z_orig, Z_cos, Z_corr, R, W, K, Phi_Rk, Phi_moe, lamb
        ):
            Z_corr = Z_orig.copy()
            if Phi_moe.shape[0] == Z_orig.shape[1]:
                Phi_moe = Phi_moe.T
            if R.shape[0] == Z_orig.shape[1]:
                R = R.T
            lamb_mat = np.diag(lamb) if len(np.shape(lamb)) == 1 else lamb
            for i in range(K):
                Phi_Rk = np.multiply(Phi_moe, R[i, :])
                x = np.dot(Phi_Rk, Phi_moe.T) + lamb_mat
                W = np.dot(np.dot(np.linalg.inv(x), Phi_Rk), Z_orig.T)
                W[0, :] = 0
                Z_corr -= np.dot(W.T, Phi_Rk)
            Z_cos = Z_corr / np.linalg.norm(Z_corr, ord=2, axis=0)
            return Z_cos, Z_corr, W, Phi_Rk

        cnmf.preprocess.moe_correct_ridge = patched_moe_correct_ridge

        # --- Monkeypatch harmony_correct_X to handle transposed outputs & negative strides ---
        def patched_harmony_correct_X(
            self, X, obs, pca, harmony_vars, theta=1, max_iter_harmony=20
        ):
            import harmonypy

            harmony_res = harmonypy.run_harmony(
                pca.copy(),
                obs,
                harmony_vars,
                max_iter_harmony=max_iter_harmony,
                theta=theta,
            )

            X_pca_harmony = harmony_res.Z_corr
            if X_pca_harmony.shape[1] == pca.shape[0]:
                X_pca_harmony = (
                    X_pca_harmony.T.copy()
                )  # .copy() avoids negative stride issues

            _, X_corr, _, _ = cnmf.preprocess.moe_correct_ridge(
                X.T.copy(),
                None,
                None,
                harmony_res.R.copy(),
                None,
                harmony_res.K,
                None,
                harmony_res.Phi_moe.copy(),
                harmony_res.lamb,
            )

            X_corr = np.array(X_corr.T).copy()
            X_corr[X_corr < 0] = 0
            return (X_corr, X_pca_harmony)

        Preprocess.harmony_correct_X = patched_harmony_correct_X
        # ----------------------------------------------------

        logger.info(
            f"[{cell_type}] Running cNMF Preprocess with covariates: {args.covariates}"
        )
        pp = Preprocess()
        try:
            pp.preprocess_for_cnmf(
                adata_ct,
                harmony_vars=args.covariates,
                n_top_rna_genes=args.num_highvar_genes,
                save_output_base=out_base,
            )
            # The Preprocess module automatically saves these files
            counts_fn = f"{out_base}.Corrected.HVG.Varnorm.h5ad"
            genes_file = f"{out_base}.Corrected.HVGs.txt"
            tpm_fn = f"{out_base}.TP10K.h5ad"
            logger.info(
                f"[{cell_type}] Preprocess completed. Normalized counts saved to {counts_fn}"
            )
        except Exception as e:
            logger.error(
                f"[{cell_type}] cNMF Preprocess failed (potentially due to package compatibility): {e}"
            )
            return
    else:
        import numpy as np
        # Standard approach without preprocessing covariates

        sc.pp.filter_cells(adata_ct, min_genes=200)
        sc.pp.filter_cells(adata_ct, min_counts=200)
        sc.pp.filter_genes(adata_ct, min_cells=3)

        # cNMF crashes if any cell has 0 counts in the highly variable genes (HVGs).
        # In highly sparse ATAC data, this happens often. To prevent it, we pre-calculate
        # the HVGs ourselves, filter cells that have 0 counts across just those HVGs,
        # and then pass our HVG list explicitly to cNMF via `genes_file` so it bypasses
        # its internal HVG variance selection and uses our EXACT safe list.

        logger.info(
            f"[{cell_type}] Pre-calculating HVGs to guarantee no zero-count cells..."
        )
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
        n_dropped = np.sum(~cells_to_keep)

        if n_dropped > 0:
            logger.info(
                f"[{cell_type}] Dropped {n_dropped} cells with zero counts across the {args.num_highvar_genes} HVGs."
            )
            adata_ct = adata_ct[cells_to_keep, :].copy()

            if adata_ct.n_obs < 50:
                logger.error(f"[{cell_type}] Too few cells remaining. Aborting.")
                return

        # Save our custom HVGs to a file and update the `genes_file` variable!
        custom_hvgs = temp_adata.var_names[hvg_mask].tolist()
        genes_file = f"{out_base}_custom_hvgs.txt"
        with open(genes_file, "w") as f:
            f.write("\n".join(custom_hvgs))

        counts_fn = f"{out_base}.counts.h5ad"
        import anndata as ad

        clean_adata = ad.AnnData(
            X=adata_ct.X, obs=adata_ct.obs.copy(), var=adata_ct.var.copy()
        )
        clean_adata.write_h5ad(counts_fn)
        logger.info(
            f"[{cell_type}] Saved raw counts to {counts_fn} and custom HVGs to {genes_file}"
        )

    # Initialize cNMF
    cnmf_obj = cNMF(output_dir=str(cnmf_dir), name=run_name)

    # Prepare
    logger.info(f"[{cell_type}] Running cNMF prepare...")
    try:
        cnmf_obj.prepare(
            counts_fn=str(counts_fn),
            components=components_to_use,
            n_iter=args.n_iter,
            num_highvar_genes=args.num_highvar_genes,
            genes_file=genes_file,
            tpm_fn=tpm_fn,
        )
    except Exception as e:
        logger.error(f"[{cell_type}] cNMF prepare failed: {e}")
        return

    # Factorize
    logger.info(f"[{cell_type}] Running cNMF factorize with {args.workers} workers...")
    try:
        cnmf_obj.factorize_multi_process(total_workers=args.workers)
    except Exception as e:
        logger.error(f"[{cell_type}] cNMF factorize failed: {e}")
        return

    # Combine
    logger.info(f"[{cell_type}] Running cNMF combine...")
    try:
        cnmf_obj.combine()
    except Exception as e:
        logger.error(f"[{cell_type}] cNMF combine failed: {e}")
        return

    # Plot K selection
    logger.info(f"[{cell_type}] Generating K selection plot...")
    try:
        cnmf_obj.k_selection_plot()
    except Exception as e:
        logger.error(f"[{cell_type}] cNMF K selection plot failed: {e}")

    # Run consensus for all K
    for k in components_to_use:
        logger.info(f"[{cell_type}] Running consensus for K={k}...")
        try:
            cnmf_obj.consensus(k=k, density_threshold=0.1, show_clustering=True)
        except Exception as e:
            logger.error(f"[{cell_type}] cNMF consensus for K={k} failed: {e}")

    logger.info(f"[{cell_type}] cNMF pipeline complete.")


def main():
    args = parse_args()
    debug = args.debug

    work_dir = Path(args.work_dir)
    quants_dir = work_dir / "quants"
    results_dir = work_dir / "results"
    info_dir = work_dir / "sample_info"
    logs_dir = work_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    cnmf_dir = results_dir / "latents" / args.cnmf_dir_name
    cnmf_dir.mkdir(parents=True, exist_ok=True)

    tmp_dir = work_dir / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    safe_ct = args.cell_type.replace(" ", "_").replace("/", "-")
    log_filename = logs_dir / f"{args.project}_{args.modality}_{safe_ct}_run_cnmf.log"
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        level=logging.DEBUG if debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler(sys.stdout)],
    )

    logger.info(f"Command line: {' '.join(sys.argv)}")

    raw_file = quants_dir / f"{args.project}.raw.multivi_prep.h5ad"
    annot_file = quants_dir / f"{args.project}.multivi.annotated.h5ad"

    # Load and Sync Data
    adata_raw, _ = load_and_prep_data(raw_file, annot_file, debug)

    if args.exclude_ids:
        exclude_ids = [x.strip() for x in args.exclude_ids.split(",")]
        logger.info(f"Excluding sample IDs: {exclude_ids}")
        if "sample_id" in adata_raw.obs:
            adata_raw = adata_raw[~adata_raw.obs["sample_id"].isin(exclude_ids)].copy()
        else:
            logger.warning("sample_id not found in adata_raw.obs, skipping exclusion.")

    # Find modality name
    modal_full = None
    for k, v in VAR_MODAL_DICT.items():
        if v == args.modality:
            modal_full = k
            break

    if not modal_full:
        logger.error(f"Modality {args.modality} not found in VAR_MODAL_DICT.")
        sys.exit(1)

    logger.info(f"Subsetting anndata by modality: {args.modality}")
    modality_list = MODAL_TYPES_DICT.get(args.modality)
    adata_modal = adata_raw[
        adata_raw.obs.modality.isin(modality_list),
        adata_raw.var.modality == modal_full,
    ].copy()

    unique_cell_types = adata_modal.obs["cell_label"].unique()
    logger.info(f"Available cell types: {unique_cell_types}")

    if args.cell_type not in unique_cell_types:
        logger.error(f"Cell type '{args.cell_type}' not found in the dataset.")
        sys.exit(1)

    adata_ct = adata_modal[adata_modal.obs["cell_label"] == args.cell_type].copy()

    import shutil

    safe_ct = args.cell_type.replace(" ", "_").replace("/", "-")
    run_name = f"{args.project}_{safe_ct}_{args.modality}"
    ct_tmp_dir = tmp_dir / run_name
    ct_tmp_dir.mkdir(parents=True, exist_ok=True)

    try:
        process_cell_type(
            adata_ct, args.cell_type, args, cnmf_dir, ct_tmp_dir, results_dir, info_dir
        )
    finally:
        if ct_tmp_dir.exists():
            shutil.rmtree(ct_tmp_dir)
            logger.info(f"Cleaned up temporary directory: {ct_tmp_dir}")

    logger.info(f"Cell type {args.cell_type} processed.")


if __name__ == "__main__":
    main()
