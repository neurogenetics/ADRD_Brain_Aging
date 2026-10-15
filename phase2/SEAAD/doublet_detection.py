import argparse
import logging
import os
import sys
import scanpy as sc

# Configure logging to stdout initially
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run Scrublet doublet detection on Scanpy AnnData."
    )
    parser.add_argument(
        "--input-file",
        type=str,
        required=True,
        help="Path to input 10x h5 file.",
    )
    parser.add_argument(
        "--output-file",
        type=str,
        required=True,
        help="Path to save the output AnnData file (h5ad).",
    )
    parser.add_argument(
        "--full-doublet-obs-file",
        type=str,
        required=True,
        help="Path to save the full predictions CSV.",
    )
    parser.add_argument(
        "--max-mito-percent",
        type=float,
        default=10.0,
        help="Maximum mitochondrial percentage (default: 10.0).",
    )
    parser.add_argument(
        "--min-genes-count",
        type=int,
        default=200,
        help="Minimum genes count per cell (default: 200).",
    )
    parser.add_argument(
        "--min-cells-count",
        type=int,
        default=3,
        help="Minimum cells count per gene (default: 3).",
    )
    parser.add_argument(
        "--max-doublet-rate",
        type=float,
        default=0.08,
        help="Expected doublet rate (default: 0.08).",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Random state for reproducibility (default: 42).",
    )
    parser.add_argument(
        "--gex-only",
        action="store_true",
        default=True,
        help="Only load gene expression data from h5 (default: True).",
    )
    parser.add_argument(
        "--no-gex-only",
        action="store_false",
        dest="gex_only",
        help="Do not restrict h5 load to gene expression only.",
    )
    parser.add_argument(
        "--plot-file",
        type=str,
        default="",
        help="Path to save the doublet score distribution plot (e.g., png/pdf).",
    )
    parser.add_argument(
        "--log-file",
        type=str,
        default="",
        help="Path to save execution log output in addition to stdout.",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Convert paths to absolute paths immediately to prevent permission/relative path errors
    args.input_file = os.path.abspath(args.input_file)
    args.output_file = os.path.abspath(args.output_file)
    args.full_doublet_obs_file = os.path.abspath(args.full_doublet_obs_file)
    if args.plot_file:
        args.plot_file = os.path.abspath(args.plot_file)
    if args.log_file:
        args.log_file = os.path.abspath(args.log_file)

    # Configure additional log file handler if specified
    if args.log_file:
        log_dir = os.path.dirname(args.log_file)
        if log_dir:
            os.makedirs(log_dir, exist_ok=True)
        file_handler = logging.FileHandler(args.log_file)
        file_handler.setFormatter(
            logging.Formatter(
                "%(asctime)s [%(levelname)s] %(message)s",
                datefmt="%Y-%m-%d %H:%M:%S",
            )
        )
        file_handler.setLevel(logging.INFO)
        logging.getLogger().addHandler(file_handler)
        logger.info(f"Logging configured to write to file: {args.log_file}")

    # Configure matplotlib to run headlessly using Agg backend
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # 1. Input validation
    if not os.path.exists(args.input_file):
        logger.error(f"Input file does not exist: {args.input_file}")
        raise FileNotFoundError(f"Input file not found: {args.input_file}")

    # Ensure output directories exist for all specified target files
    for file_path in [args.output_file, args.full_doublet_obs_file, args.plot_file]:
        if file_path:
            dir_name = os.path.dirname(file_path)
            if dir_name:
                os.makedirs(dir_name, exist_ok=True)

    logger.info(f"Loading 10x h5 data from {args.input_file} (gex_only={args.gex_only})...")
    adata = sc.read_10x_h5(args.input_file, gex_only=args.gex_only)

    logger.info(f"Loaded AnnData with shape: {adata.shape}")

    # Ensure variable names are unique
    adata.var_names_make_unique()

    # Annotate mitochondrial genes to assess cell viability
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    # Filter out empty droplets and dead cells
    logger.info(f"Filtering cells with < {args.min_genes_count} genes...")
    sc.pp.filter_cells(adata, min_genes=args.min_genes_count)

    logger.info(f"Filtering genes with < {args.min_cells_count} cells...")
    sc.pp.filter_genes(adata, min_cells=args.min_cells_count)

    logger.info(f"Filtering cells with >= {args.max_mito_percent}% mitochondrial counts...")
    adata = adata[adata.obs.pct_counts_mt < args.max_mito_percent, :]

    logger.info(f"AnnData shape after cell/gene QC filtering: {adata.shape}")

    # Filter universally unexpressed genes
    # (Scrublet needs the rest of the broad gene space to build its KNN graph)
    logger.info("Filtering universally unexpressed genes (min_cells=3) for Scrublet KNN graph...")
    sc.pp.filter_genes(adata, min_cells=3)

    # Save the raw integer counts to a dedicated layer
    adata.layers["counts"] = adata.X.copy()

    # Run Scrublet directly through Scanpy
    logger.info(
        f"Running Scrublet (expected_doublet_rate={args.max_doublet_rate}, "
        f"random_state={args.random_state})..."
    )
    sc.pp.scrublet(
        adata,
        expected_doublet_rate=args.max_doublet_rate,
        random_state=args.random_state,
    )

    # Report doublet detection results
    num_cells = adata.n_obs
    if num_cells > 0:
        num_doublets = (adata.obs["predicted_doublet"] == True).sum()
        pct_doublets = (num_doublets / num_cells) * 100
        logger.info(
            f"Doublet calling results: identified {num_doublets} doublets out of "
            f"{num_cells} total cells ({pct_doublets:.2f}%)."
        )
    else:
        logger.warning("No cells remaining to run doublet calling calculations.")

    # Save score distribution plot if file path is provided
    if args.plot_file:
        logger.info("Generating and saving Scrublet score distribution plot...")
        sc.pl.scrublet_score_distribution(adata, show=False)
        plt.savefig(args.plot_file, bbox_inches="tight")
        logger.info(f"Saved score distribution plot to {args.plot_file}")

    # Save the full predictions
    logger.info(f"Saving full doublet predictions observation metadata to {args.full_doublet_obs_file}...")
    adata.obs.to_csv(args.full_doublet_obs_file)

    # Filter the doublets out of the dataset
    logger.info("Filtering out predicted doublets...")
    adata = adata[adata.obs["predicted_doublet"] == False, :].copy()
    logger.info(f"AnnData shape after doublet filtering: {adata.shape}")

    # Save the non-doublet adata
    logger.info(f"Saving filtered non-doublet AnnData to {args.output_file}...")
    adata.write(args.output_file)
    logger.info("Doublet detection workflow successfully completed.")


if __name__ == "__main__":
    main()
