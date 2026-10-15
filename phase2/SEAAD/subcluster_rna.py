import os
import sys
import logging
import argparse
import warnings
import scanpy as sc
import mudata as md
import scvi
import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from sklearn.metrics import silhouette_score
import hdf5plugin

# Silence warnings to keep logs clean
from pandas.errors import PerformanceWarning
warnings.filterwarnings("ignore", category=FutureWarning, module="pandas")
warnings.filterwarnings("ignore", category=FutureWarning, module="mudata")
warnings.filterwarnings("ignore", category=UserWarning, module="scanpy")
warnings.filterwarnings("ignore", category=PerformanceWarning)

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
        description="Perform broad cell-type subclustering on GEX (RNA) data using scvi-tools."
    )
    # File inputs/outputs
    parser.add_argument(
        "--raw-input",
        type=str,
        required=True,
        help="Path to the raw/unprocessed RNA AnnData (.h5ad) or MuData (.h5mu) file",
    )
    parser.add_argument(
        "--processed-h5ad",
        type=str,
        required=True,
        help="Path to the processed/clustered GEX AnnData (.h5ad) file containing cell labels",
    )
    parser.add_argument(
        "--resolution-dir",
        type=str,
        required=True,
        help="Directory where all subclustering output files (h5ad, plot, logs, markers) will be saved",
    )
    parser.add_argument(
        "--name-prefix",
        type=str,
        required=True,
        help="Prefix name for the output files (e.g. 'ExN' or 'InN') to prevent files overwriting each other",
    )
    parser.add_argument(
        "--model-dir",
        type=str,
        default=None,
        help="Directory to save the newly trained scVI subclustering model (optional)",
    )

    # Cohort subsetting
    parser.add_argument(
        "--subset-col",
        type=str,
        required=True,
        help="Obs column in processed AnnData to filter cells by (e.g., 'leiden_scvi' or 'celltype_major')",
    )
    parser.add_argument(
        "--subset-values",
        type=str,
        required=True,
        help="Comma-separated list of values in subset-col to keep for subclustering (e.g. '2,3,4' or 'Excitatory Neurons')",
    )
    parser.add_argument(
        "--annotation-col",
        type=str,
        default="celltype_major",
        help="Obs column in processed AnnData to copy over as the 'initial_label' metadata (default: 'celltype_major')",
    )

    # Preprocessing
    parser.add_argument(
        "--min-cell-percent",
        type=float,
        default=0.005,
        help="Minimum percentage of cells a gene must be expressed in to be kept (default: 0.005)",
    )
    parser.add_argument(
        "--max-mito-percent",
        type=float,
        default=10.0,
        help="Maximum percentage of mitochondrial counts per cell (default: 10.0)",
    )
    parser.add_argument(
        "--detect-hv-features",
        action="store_true",
        help="Whether to detect highly variable features for scVI",
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
        help="Percentage of top highly variable genes to keep if --detect-hv-features is used (default: 0.15)",
    )

    # scVI parameters
    parser.add_argument(
        "--batch-key",
        type=str,
        default="sample_id",
        help="Obs column to use as batch_key in scVI (default: 'sample_id')",
    )
    parser.add_argument(
        "--categorical-covariates",
        type=str,
        default="donor_id,sex,race",
        help="Comma-separated list of categorical covariates for scVI (default: 'donor_id,sex,race')",
    )
    parser.add_argument(
        "--continuous-covariates",
        type=str,
        default="pct_counts_mt,pct_counts_ribo",
        help="Comma-separated list of continuous covariates for scVI (default: 'pct_counts_mt,pct_counts_ribo')",
    )
    parser.add_argument(
        "--n-latent",
        type=int,
        default=12,
        help="Number of latent dimensions for scVI (default: 12)",
    )
    parser.add_argument(
        "--n-layers",
        type=int,
        default=2,
        help="Number of layers for scVI (default: 2)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=10000,
        help="Batch size for scVI training (default: 10000)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for scvi and clustering (default: 42)",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Set seeds
    scvi.settings.seed = args.seed

    # Create directory
    os.makedirs(args.resolution_dir, exist_ok=True)

    # Define prefix-aware output file paths inside resolution_dir
    output_h5ad_path = os.path.join(args.resolution_dir, f"{args.name_prefix}_subclustered.h5ad")
    output_plot_path = os.path.join(args.resolution_dir, f"{args.name_prefix}_subclustering_sweep.png")
    log_file_path = os.path.join(args.resolution_dir, f"{args.name_prefix}_subcluster_rna.log")

    # Configure file logging
    file_handler = logging.FileHandler(log_file_path, mode="a")
    file_handler.setFormatter(
        logging.Formatter(
            "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )
    )
    logging.getLogger().addHandler(file_handler)

    logger.info("=== Starting GEX Broadtype Subclustering ===")

    # 1. Load data
    logger.info("Loading raw input: %s", args.raw_input)
    try:
        if args.raw_input.endswith(".h5mu"):
            logger.info("Loading raw as MuData object...")
            mdata_raw = md.read(args.raw_input)
            if "rna" not in mdata_raw.mod:
                logger.error("Modality 'rna' not found in raw MuData object.")
                sys.exit(1)
            adata_raw = mdata_raw.mod["rna"].copy()
            # Copy global obs to rna obs
            for col in mdata_raw.obs.columns:
                if col not in adata_raw.obs.columns:
                    adata_raw.obs[col] = mdata_raw.obs[col]
        else:
            logger.info("Loading raw as AnnData object...")
            adata_raw = sc.read_h5ad(args.raw_input)
    except Exception as e:
        logger.exception("Failed to load raw input: %s", str(e))
        sys.exit(1)

    logger.info("Loading processed AnnData: %s", args.processed_h5ad)
    try:
        adata_proc = sc.read_h5ad(args.processed_h5ad)
    except Exception as e:
        logger.exception("Failed to load processed AnnData: %s", str(e))
        sys.exit(1)

    # 2. Subset cells
    subset_vals = [x.strip() for x in args.subset_values.split(",")]
    logger.info(
        "Filtering processed cells where %s is in %s", args.subset_col, subset_vals
    )

    if args.subset_col not in adata_proc.obs.columns:
        logger.error(
            "Subset column '%s' not found in processed AnnData obs columns.",
            args.subset_col,
        )
        sys.exit(1)

    # Cast categories or objects to string for matching
    matching_barcodes = adata_proc.obs[
        adata_proc.obs[args.subset_col].astype(str).isin(subset_vals)
    ].index
    logger.info(
        "Found %d matching cell barcodes in processed dataset.", len(matching_barcodes)
    )

    if len(matching_barcodes) == 0:
        logger.error("No barcodes matched the filter. Aborting.")
        sys.exit(1)

    # Subset the raw object
    adata_gex = adata_raw[adata_raw.obs.index.isin(matching_barcodes)].copy()
    logger.info("Subset raw AnnData shape: %s", adata_gex.shape)

    # Bring over 'initial_label' from the processed AnnData
    if args.annotation_col in adata_proc.obs.columns:
        logger.info(
            "Mapping initial annotations from processed obs['%s'] as obs['initial_label']...",
            args.annotation_col,
        )
        adata_gex.obs["initial_label"] = adata_gex.obs.index.map(
            adata_proc.obs[args.annotation_col]
        )
    else:
        logger.warning(
            "Annotation column '%s' not found. Setting initial_label to 'Unknown'.",
            args.annotation_col,
        )
        adata_gex.obs["initial_label"] = "Unknown"

    # Merge overall obs metadata from processed to ensure other study covariates exist
    for col in adata_proc.obs.columns:
        if col not in adata_gex.obs.columns:
            adata_gex.obs[col] = adata_proc.obs[col]

    # 3. Quality Control
    logger.info("Calculating QC metrics...")
    adata_gex.var["mt"] = adata_gex.var_names.str.startswith("MT-")
    adata_gex.var["ribo"] = adata_gex.var_names.str.startswith(("RPS", "RPL"))
    adata_gex.var["hb"] = adata_gex.var_names.str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(
        adata_gex, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )

    logger.info(
        "Filtering cells with >%.2f%% mitochondrial counts...", args.max_mito_percent
    )
    adata_gex = adata_gex[adata_gex.obs.pct_counts_mt < args.max_mito_percent, :].copy()

    logger.info("Filtering cells with <200 genes...")
    sc.pp.filter_cells(adata_gex, min_genes=200)

    min_cells = int(adata_gex.shape[0] * args.min_cell_percent)
    logger.info(
        "Filtering genes expressed in <%.2f%% of cells (min_cells=%d)...",
        args.min_cell_percent * 100,
        min_cells,
    )
    sc.pp.filter_genes(adata_gex, min_cells=min_cells)

    logger.info("Shape after QC and filtering: %s", adata_gex.shape)

    # 4. Highly Variable Features
    if args.detect_hv_features:
        n_top_genes = int(adata_gex.n_vars * args.top_features_percent)
        logger.info("Detecting top %d highly variable genes...", n_top_genes)
        try:
            sc.pp.highly_variable_genes(
                adata_gex,
                n_top_genes=n_top_genes,
                batch_key=args.batch_key if args.batch_key in adata_gex.obs else None,
                flavor="seurat_v3",
                subset=True,
            )
            logger.info("Shape after highly variable genes step: %s", adata_gex.shape)

            # Save the variable features for downstream use
            hvg_file = os.path.join(args.resolution_dir, f"{args.name_prefix}_variable_features.txt")
            logger.info("Saving highly variable features index list to: %s", hvg_file)
            variable_genes = adata_gex.var[
                adata_gex.var["highly_variable"]
            ].index.to_frame()
            variable_genes.to_csv(hvg_file, index=False, header=False)
        except Exception as e:
            logger.warning(
                "Failed to calculate highly variable genes. Error: %s", str(e)
            )

    # 5. scVI Setup
    cat_covs = (
        [x.strip() for x in args.categorical_covariates.split(",") if x.strip()]
        if args.categorical_covariates
        else None
    )
    cont_covs = (
        [x.strip() for x in args.continuous_covariates.split(",") if x.strip()]
        if args.continuous_covariates
        else None
    )

    # Validate batch_key and covariates
    batch_key = args.batch_key if args.batch_key in adata_gex.obs.columns else None
    valid_cat_covs = [c for c in (cat_covs or []) if c in adata_gex.obs.columns]
    valid_cont_covs = [c for c in (cont_covs or []) if c in adata_gex.obs.columns]

    # Convert categoricals
    for c in valid_cat_covs:
        adata_gex.obs[c] = adata_gex.obs[c].astype("category")

    logger.info("Setting up scVI model...")
    scvi.model.SCVI.setup_anndata(
        adata_gex,
        batch_key=batch_key,
        categorical_covariate_keys=valid_cat_covs if valid_cat_covs else None,
        continuous_covariate_keys=valid_cont_covs if valid_cont_covs else None,
    )

    # 6. Train scVI Model
    logger.info("Initializing scVI subclustering model...")
    model = scvi.model.SCVI(adata_gex, n_latent=args.n_latent, n_layers=args.n_layers)

    logger.info("Training model with batch_size=%d...", args.batch_size)
    model.train(batch_size=args.batch_size)

    if args.model_dir:
        logger.info("Saving trained scVI subclustering model to: %s", args.model_dir)
        model.save(args.model_dir, overwrite=True)

    # 7. Latent representation and layers
    logger.info("Extracting scVI latent representation...")
    adata_gex.obsm["scvi_latent"] = model.get_latent_representation()

    logger.info("Extracting scVI normalized expression...")
    adata_gex.layers["scvi_normalized"] = model.get_normalized_expression(
        library_size=10e4
    )

    logger.info("Computing log-transformed normalized layer...")
    adata_gex.layers["scvi_normalized_log1p"] = np.log1p(
        adata_gex.layers["scvi_normalized"]
    )

    logger.info("Computing neighborhood graph and UMAP...")
    sc.pp.neighbors(adata_gex, use_rep="scvi_latent")
    sc.tl.umap(adata_gex)

    # 8. Resolution Sweep and Silhouette Analysis
    resolutions_to_try = np.arange(0.1, 1.0, 0.1)
    logger.info(
        "Running Leiden clustering sweep over resolutions: %s",
        list(np.round(resolutions_to_try, 2)),
    )

    clust_assign_by_res = None
    mean_scores = {}
    largest_score = -1.0
    best_res = 0.1
    new_leiden_key = "leiden_scvi"

    # Setup igraph if available
    leiden_kwargs = {}
    try:
        import igraph
        import leidenalg

        leiden_kwargs = {"flavor": "igraph", "n_iterations": 2}
    except ImportError:
        logger.warning("igraph or leidenalg not found. Using default scanpy settings.")

    for res in resolutions_to_try:
        res = round(res, 2)
        logger.info("  Testing resolution: %.1f", res)

        # Run Leiden
        sc.tl.leiden(
            adata_gex, key_added=new_leiden_key, resolution=res, **leiden_kwargs
        )

        # Silhouette
        score = silhouette_score(
            adata_gex.obsm["scvi_latent"], adata_gex.obs[new_leiden_key]
        )
        num_clusters = adata_gex.obs[new_leiden_key].nunique()
        logger.info("    -> Silhouette: %.3f, Clusters: %d", score, num_clusters)

        # Donor counts check (safely check if sample_id is present)
        mean_sample_per_cluster = 0.0
        less_than_half = 0
        if "sample_id" in adata_gex.obs.columns:
            df_grouped = adata_gex.obs.groupby([new_leiden_key])[
                "sample_id"
            ].value_counts()
            df_grouped = df_grouped[df_grouped >= 30].to_frame().reset_index()
            df_grouped = df_grouped.groupby(new_leiden_key)["sample_id"].nunique()
            mean_sample_per_cluster = df_grouped.mean()
            less_than_half = df_grouped[
                df_grouped < adata_gex.obs.sample_id.nunique() / 3
            ].shape[0]

        mean_cell_per_cluster = adata_gex.obs[new_leiden_key].value_counts().mean()
        mean_scores[res] = [
            score,
            num_clusters,
            mean_sample_per_cluster,
            mean_cell_per_cluster,
            less_than_half,
        ]

        # Retain cluster assignments
        res_col_name = f"leiden_{res}"
        if clust_assign_by_res is None:
            clust_assign_by_res = (
                adata_gex.obs[[new_leiden_key]]
                .copy()
                .rename(columns={new_leiden_key: res_col_name})
            )
        else:
            clust_assign_by_res = pd.concat(
                [
                    clust_assign_by_res,
                    adata_gex.obs[[new_leiden_key]]
                    .copy()
                    .rename(columns={new_leiden_key: res_col_name}),
                ],
                axis="columns",
            )

        # Save normalized expression averages per resolution
        logger.info(
            "    -> Writing average expression and marker CSV files for resolution %.1f...",
            res,
        )
        try:
            avgexp = (
                sc.get.obs_df(
                    adata_gex, keys=list(adata_gex.var_names), layer="scvi_normalized"
                )
                .groupby(adata_gex.obs[new_leiden_key])
                .mean()
            )

            res_avg_exp_file = os.path.join(args.resolution_dir, f"{args.name_prefix}_avgexp_res{res}.csv")
            avgexp.to_csv(res_avg_exp_file)

            # Marker calculation
            sc.tl.rank_genes_groups(
                adata_gex,
                groupby=new_leiden_key,
                method="wilcoxon",
                pts=True,
                layer="scvi_normalized_log1p",
            )
            markers_df = sc.get.rank_genes_groups_df(adata_gex, group=None)
            res_markers_file = os.path.join(
                args.resolution_dir, f"{args.name_prefix}_markers_res{res}.csv"
            )
            markers_df.to_csv(res_markers_file, index=False)
        except Exception as e:
            logger.warning(
                "    -> Failed to write markers/averages for resolution %.1f: %s",
                res,
                str(e),
            )

        if score > largest_score:
            largest_score = score
            best_res = res

    # Save the aggregated resolutions metadata CSV
    res_obs_file = os.path.join(args.resolution_dir, f"{args.name_prefix}_resolution_assignments_obs.csv")
    logger.info("Saving multi-resolution assignments and metadata to: %s", res_obs_file)
    try:
        res_obs_info = pd.concat([clust_assign_by_res, adata_gex.obs], axis="columns")
        res_obs_info.to_csv(res_obs_file)
    except Exception as e:
        logger.warning(
            "Failed to save final merged obs assignments CSV. Error: %s", str(e)
        )

    # Run final clustering with the best resolution
    best_res = round(best_res, 2)
    logger.info(
        "Best resolution selected by Silhouette Score: %.2f (Score: %.3f)",
        best_res,
        largest_score,
    )
    sc.tl.leiden(
        adata_gex, key_added=new_leiden_key, resolution=best_res, **leiden_kwargs
    )

    # 9. Visualization (Multi-panel Plot)
    logger.info("Generating final multi-panel visualization plot...")
    try:
        scores_df = pd.DataFrame(index=mean_scores.keys(), data=mean_scores.values())
        scores_df.columns = [
            "score",
            "num_clusters",
            "mean_samples",
            "mean_cells",
            "less_than_half",
        ]

        fig, axes = plt.subplots(1, 3, figsize=(24, 7), dpi=300)

        # Panel 1: UMAP colored by final Leiden subclusters
        sc.pl.umap(
            adata_gex,
            color=new_leiden_key,
            ax=axes[0],
            show=False,
            legend_loc="on data",
            legend_fontsize=8,
            legend_fontweight="bold",
        )
        axes[0].set_title(
            f"Subclusters (Res = {best_res})", fontsize=14, fontweight="bold"
        )

        # Panel 2: UMAP colored by initial cell labels
        sc.pl.umap(
            adata_gex, color="initial_label", ax=axes[1], show=False, legend_fontsize=8
        )
        axes[1].set_title("Initial Cell Labels", fontsize=14, fontweight="bold")

        # Panel 3: Silhouette score plot
        axes[2].plot(
            scores_df.index,
            scores_df["score"],
            marker="o",
            color="darkred",
            linewidth=2,
        )
        axes[2].set_xlabel("Resolution", fontsize=11)
        axes[2].set_ylabel("Average Silhouette Score", fontsize=11)
        axes[2].set_title("Silhouette Sweep Analysis", fontsize=14, fontweight="bold")
        axes[2].grid(True, linestyle="--", alpha=0.5)

        plt.tight_layout()
        plt.savefig(output_plot_path, bbox_inches="tight")
        plt.close()
        logger.info("UMAP multi-panel plot saved to: %s", output_plot_path)
    except Exception as e:
        logger.exception("Failed to write visualization plots: %s", str(e))

    # 10. Save the final AnnData
    logger.info("Saving final subclustered AnnData to: %s", output_h5ad_path)
    try:
        # Sanitize obs to prevent h5py write exceptions for object/mixed columns
        for col in adata_gex.obs.columns:
            if adata_gex.obs[col].dtype == "object":
                adata_gex.obs[col] = adata_gex.obs[col].fillna("").astype(str)

        adata_gex.write(output_h5ad_path)
        logger.info("=== GEX Subclustering Completed Successfully ===")
    except Exception as e:
        logger.exception("Failed to write subclustered AnnData file: %s", str(e))
        sys.exit(1)


if __name__ == "__main__":
    main()
