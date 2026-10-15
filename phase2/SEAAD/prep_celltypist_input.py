import os
import sys
import logging
import argparse
import warnings
import scanpy as sc
import mudata as md
import hdf5plugin

# Silence warnings to keep console clean
warnings.filterwarnings("ignore", category=FutureWarning, module="mudata")
warnings.filterwarnings("ignore", category=UserWarning, module="scanpy")

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(
        description="Prepare a combined MuData object or raw AnnData file for use with CellTypist cell-type prediction models."
    )
    parser.add_argument(
        "--input-file",
        type=str,
        required=True,
        help="Path to the input combined MuData (.h5mu) or AnnData (.h5ad) file",
    )
    parser.add_argument(
        "--output-file",
        type=str,
        required=True,
        help="Path where the prepped CellTypist-ready AnnData (.h5ad) file will be saved",
    )
    parser.add_argument(
        "--target-sum",
        type=float,
        default=1e4,
        help="Target sum for normalization (default: 10000 / 1e4)",
    )
    parser.add_argument(
        "--detect-hv-features",
        action="store_true",
        help="Whether to detect highly variable features (without subsetting by default unless --filter-hv-features is also set)",
    )
    parser.add_argument(
        "--filter-hv-features",
        action="store_true",
        help="Whether to subset/filter the AnnData to highly variable features",
    )
    parser.add_argument(
        "--top-features-percent",
        type=float,
        default=0.15,
        help="Percentage of top highly variable genes to identify (default: 0.15)",
    )
    parser.add_argument(
        "--batch-key",
        type=str,
        default=None,
        help="Optional batch key to use for highly variable genes detection",
    )
    parser.add_argument(
        "--min-genes",
        type=int,
        default=None,
        help="Optional minimum genes per cell filter",
    )
    parser.add_argument(
        "--max-mito-percent",
        type=float,
        default=None,
        help="Optional maximum mitochondrial percentage filter",
    )
    return parser.parse_args()


def peek_anndata(adata, message=""):
    logger.info("%s", message)
    logger.info("  Shape: %s", adata.shape)
    logger.info("  Obs columns: %s", list(adata.obs.columns))
    logger.info("  Var columns: %s", list(adata.var.columns))


def main():
    args = parse_args()

    # Ensure output parent directory exists
    out_dir = os.path.dirname(os.path.abspath(args.output_file))
    os.makedirs(out_dir, exist_ok=True)

    # Configure file logging in addition to terminal
    log_file_path = os.path.join(out_dir, "prep_celltypist_input.log")
    file_handler = logging.FileHandler(log_file_path, mode="a")
    file_handler.setFormatter(
        logging.Formatter(
            "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )
    )
    logging.getLogger().addHandler(file_handler)

    logger.info("=== Starting CellTypist Input Prep ===")
    logger.info("Reading input file: %s", args.input_file)

    try:
        if args.input_file.endswith(".h5mu"):
            logger.info("Loading as MuData object...")
            mdata = md.read(args.input_file)
            if "rna" not in mdata.mod:
                logger.error("Modality 'rna' not found in MuData object.")
                sys.exit(1)
            adata = mdata.mod["rna"].copy()
            # Copy global obs to rna obs
            for col in mdata.obs.columns:
                if col not in adata.obs.columns:
                    adata.obs[col] = mdata.obs[col]
        else:
            logger.info("Loading as AnnData object...")
            adata = sc.read_h5ad(args.input_file)
    except Exception as e:
        logger.exception("Failed to load input file: %s", str(e))
        sys.exit(1)

    peek_anndata(adata, "Loaded data structure:")

    # Subset to just gene expression features if it contains multimodal features (and wasn't parsed from h5mu)
    if "modality" in adata.var.columns:
        logger.info("Subsetting to Gene Expression features based on 'modality' column...")
        adata = adata[:, adata.var.modality == "Gene Expression"].copy()
        peek_anndata(adata, "Data after subsetting to GEX features:")

    # Calculate QC metrics
    logger.info("Calculating QC metrics...")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")
    
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )

    # Optional QC filtering if specified
    if args.max_mito_percent is not None:
        logger.info("Filtering cells with >%.2f%% mitochondrial counts...", args.max_mito_percent)
        adata = adata[adata.obs.pct_counts_mt < args.max_mito_percent, :].copy()

    if args.min_genes is not None:
        logger.info("Filtering cells with < %d genes...", args.min_genes)
        sc.pp.filter_cells(adata, min_genes=args.min_genes)

    # Identify highly variable genes
    if args.detect_hv_features:
        n_top_genes = int(adata.n_vars * args.top_features_percent)
        logger.info(
            "Detecting top %d highly variable features (subset=%s)...",
            n_top_genes,
            args.filter_hv_features,
        )
        try:
            sc.pp.highly_variable_genes(
                adata,
                n_top_genes=n_top_genes,
                batch_key=args.batch_key if args.batch_key in adata.obs.columns else None,
                flavor="seurat_v3",
                subset=args.filter_hv_features,
            )
            logger.info("Shape after HVG step: %s", adata.shape)
        except Exception as e:
            logger.warning("Failed to run highly_variable_genes (likely due to non-integer counts). Error: %s", str(e))

    # Normalization (CellTypist requires 10K normalization and log1p)
    logger.info("Normalizing to total counts = %s...", args.target_sum)
    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    
    logger.info("Applying log1p transformation...")
    sc.pp.log1p(adata)

    logger.info("Freezing state in .raw layer...")
    adata.raw = adata

    # Remove log1p from uns to prevent older anndata versions in Singularity/Docker containers 
    # (such as CellTypist's container) from failing to read the file due to the encoding_type='null' mismatch on uns['log1p']['base']
    logger.info("Removing log1p metadata from .uns to ensure compatibility with older anndata readers in containers...")
    adata.uns.pop("log1p", None)

    peek_anndata(adata, "Final prepped AnnData structure:")

    logger.info("Saving prepped AnnData to: %s", args.output_file)
    try:
        # Sanitize obs to prevent h5py write exceptions for object/mixed columns
        for col in adata.obs.columns:
            if adata.obs[col].dtype == "object":
                adata.obs[col] = adata.obs[col].fillna("").astype(str)
                
        adata.write(args.output_file)
        logger.info("=== CellTypist Input Prep Completed Successfully ===")
    except Exception as e:
        logger.exception("Failed to write output file: %s", str(e))
        sys.exit(1)


if __name__ == "__main__":
    main()
