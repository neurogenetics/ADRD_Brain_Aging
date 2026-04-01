import argparse
import logging
import sys
import multiprocessing
from pathlib import Path

import scanpy as sc
from cnmf import cNMF

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
        "--components", type=int, nargs="+", default=[10, 15, 20, 25, 30]
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
    parser.add_argument("--workers", type=int, default=multiprocessing.cpu_count())
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()


def process_cell_type(adata_ct, cell_type, args, cnmf_dir, tmp_dir):
    logger.info(f"--- Processing Cell Type: {cell_type} ---")
    if adata_ct.n_obs < 50:
        logger.warning(f"Skipping {cell_type}: only {adata_ct.n_obs} cells.")
        return

    # cNMF expects dense count arrays or efficiently structured sparse arrays without all zero columns
    # sc.pp.filter_genes(adata_ct, min_cells=1)

    safe_ct = cell_type.replace(" ", "_").replace("/", "-")
    run_name = f"{args.project}_{safe_ct}_{args.modality}"
    out_base = str(tmp_dir / run_name)

    genes_file = None
    tpm_fn = None

    if args.covariates:
        # Use cNMF Preprocess module to integrate known technical covariates via Harmony
        from cnmf.preprocess import Preprocess
        import cnmf.preprocess
        import numpy as np

        # --- Monkeypatch cNMF's buggy moe_correct_ridge ---
        def patched_moe_correct_ridge(Z_orig, Z_cos, Z_corr, R, W, K, Phi_Rk, Phi_moe, lamb):
            Z_corr = Z_orig.copy()
            if Phi_moe.shape[0] == Z_orig.shape[1]:  
                Phi_moe = Phi_moe.T
            if R.shape[0] == Z_orig.shape[1]:        
                R = R.T
            lamb_mat = np.diag(lamb) if len(np.shape(lamb)) == 1 else lamb
            for i in range(K):
                Phi_Rk = np.multiply(Phi_moe, R[i,:])
                x = np.dot(Phi_Rk, Phi_moe.T) + lamb_mat
                W = np.dot(np.dot(np.linalg.inv(x), Phi_Rk), Z_orig.T)
                W[0,:] = 0 
                Z_corr -= np.dot(W.T, Phi_Rk)
            Z_cos = Z_corr / np.linalg.norm(Z_corr, ord=2, axis=0)
            return Z_cos, Z_corr, W, Phi_Rk

        cnmf.preprocess.moe_correct_ridge = patched_moe_correct_ridge

        # --- Monkeypatch harmony_correct_X to handle transposed outputs & negative strides ---
        def patched_harmony_correct_X(self, X, obs, pca, harmony_vars, theta=1, max_iter_harmony=20):
            import harmonypy
            harmony_res = harmonypy.run_harmony(pca.copy(), obs, harmony_vars, max_iter_harmony=max_iter_harmony, theta=theta)

            X_pca_harmony = harmony_res.Z_corr
            if X_pca_harmony.shape[1] == pca.shape[0]:
                X_pca_harmony = X_pca_harmony.T.copy() # .copy() avoids negative stride issues

            _, X_corr, _, _ = cnmf.preprocess.moe_correct_ridge(
                X.T.copy(), None, None, harmony_res.R.copy(), None, harmony_res.K, None, harmony_res.Phi_moe.copy(), harmony_res.lamb
            )

            X_corr = np.array(X_corr.T).copy()
            X_corr[X_corr<0] = 0
            return (X_corr, X_pca_harmony)

        Preprocess.harmony_correct_X = patched_harmony_correct_X
        # ----------------------------------------------------

        logger.info(f"[{cell_type}] Running cNMF Preprocess with covariates: {args.covariates}")
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
        # Standard approach without preprocessing covariates
        counts_fn = f"{out_base}.counts.h5ad"
        adata_ct.write_h5ad(counts_fn)
        logger.info(
            f"[{cell_type}] Saved raw counts to {counts_fn} ({adata_ct.shape[0]} cells, {adata_ct.shape[1]} features)"
        )

    # Initialize cNMF
    cnmf_obj = cNMF(output_dir=str(cnmf_dir), name=run_name)

    # Prepare
    logger.info(f"[{cell_type}] Running cNMF prepare...")
    try:
        cnmf_obj.prepare(
            counts_fn=str(counts_fn),
            components=args.components,
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
    for k in args.components:
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
    logs_dir = work_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    cnmf_dir = results_dir / "cnmf"
    cnmf_dir.mkdir(parents=True, exist_ok=True)

    tmp_dir = work_dir / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    log_filename = logs_dir / f"{args.project}_{args.modality}_run_cnmf.log"
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
    logger.info(
        f"Found {len(unique_cell_types)} cell types to process: {unique_cell_types}"
    )

    for ct in unique_cell_types:
        adata_ct = adata_modal[adata_modal.obs["cell_label"] == ct].copy()
        process_cell_type(adata_ct, ct, args, cnmf_dir, tmp_dir)

    logger.info("All cell types processed.")


if __name__ == "__main__":
    main()

