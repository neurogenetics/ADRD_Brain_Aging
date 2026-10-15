import os
import sys
import logging
import argparse
import warnings
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import hdf5plugin

# Silence warnings
warnings.filterwarnings("ignore", category=FutureWarning, module="pandas")
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
        description="Integrate CellTypist majority voting predictions with cluster assignments and annotate clusters."
    )
    parser.add_argument(
        "--input-h5ad",
        type=str,
        required=True,
        help="Path to the clustered RNA AnnData (.h5ad) file",
    )
    parser.add_argument(
        "--predictions-csv",
        type=str,
        required=True,
        help="Path to the CellTypist predicted labels CSV file",
    )
    parser.add_argument(
        "--output-h5ad",
        type=str,
        required=True,
        help="Path where the finalized annotated RNA AnnData file will be saved",
    )
    parser.add_argument(
        "--output-plot",
        type=str,
        required=True,
        help="Path where the final UMAP visualization plot (e.g. .png) will be saved",
    )
    parser.add_argument(
        "--leiden-key",
        type=str,
        default="leiden_scvi",
        help="AnnData .obs column containing the cluster IDs (default: 'leiden_scvi')",
    )
    parser.add_argument(
        "--predictions-column",
        type=str,
        default="majority_voting",
        help="Column in predictions CSV representing the predicted labels (default: 'majority_voting')",
    )
    parser.add_argument(
        "--palette",
        type=str,
        default="tab20b",
        help="Color palette or colormap for plotting. Can be a standard name (like 'tab20b', 'tab20c', 'zeileis', 'godtlib') or a comma-separated list of colors/hex-codes (default: 'tab20b')",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Ensure output parent directories exist
    out_dir = os.path.dirname(os.path.abspath(args.output_h5ad))
    os.makedirs(out_dir, exist_ok=True)
    plot_dir = os.path.dirname(os.path.abspath(args.output_plot))
    os.makedirs(plot_dir, exist_ok=True)

    # Configure file logging
    log_file_path = os.path.join(out_dir, "annotate_clusters.log")
    file_handler = logging.FileHandler(log_file_path, mode="a")
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    )
    logging.getLogger().addHandler(file_handler)

    logger.info("=== Starting Cluster Annotation ===")
    logger.info("Reading input AnnData: %s", args.input_h5ad)
    try:
        adata = sc.read_h5ad(args.input_h5ad)
    except Exception as e:
        logger.exception("Failed to read input AnnData: %s", str(e))
        sys.exit(1)

    logger.info("Reading CellTypist predictions: %s", args.predictions_csv)
    try:
        pred_df = pd.read_csv(args.predictions_csv, index_col=0)
    except Exception as e:
        logger.exception("Failed to read predictions CSV: %s", str(e))
        sys.exit(1)

    if args.predictions_column not in pred_df.columns:
        logger.error("Predictions column '%s' not found in CSV. Available columns: %s", args.predictions_column, list(pred_df.columns))
        sys.exit(1)

    if args.leiden_key not in adata.obs.columns:
        logger.error("Leiden key '%s' not found in AnnData obs. Available columns: %s", args.leiden_key, list(adata.obs.columns))
        sys.exit(1)

    logger.info("Integrating predicted cell-type labels into AnnData...")
    adata.obs["predicted_celltype"] = adata.obs.index.map(pred_df[args.predictions_column])
    
    # Handle cells that might have been filtered or missing in predictions
    missing_cnt = adata.obs["predicted_celltype"].isna().sum()
    if missing_cnt > 0:
        logger.warning("%d cells in AnnData are missing predictions. Labeling as 'Unknown'.", missing_cnt)
        adata.obs["predicted_celltype"] = adata.obs["predicted_celltype"].fillna("Unknown")

    logger.info("Mapping majority CellTypist labels to Leiden clusters...")
    
    # Dictionaries to hold the mappings
    major_labels_dict = {}
    detailed_labels_dict = {}

    for cluster in sorted(adata.obs[args.leiden_key].unique()):
        # Subset to cells in this cluster
        cluster_obs = adata.obs[adata.obs[args.leiden_key] == cluster]
        
        # Calculate cell type value counts
        cell_type_counts = cluster_obs["predicted_celltype"].value_counts()
        
        if cell_type_counts.empty:
            logger.warning("Cluster %s has no cell type counts.", cluster)
            major_labels_dict[cluster] = "Unknown"
            detailed_labels_dict[cluster] = f"Unknown-(0.00%)-{cluster}"
            continue

        # Top majority cell type
        majority_celltype = cell_type_counts.index[0]
        majority_count = cell_type_counts.iloc[0]
        total_cluster_cells = cell_type_counts.sum()
        percentage = (majority_count / total_cluster_cells) * 100

        logger.info("  Cluster %s: Majority CellType = %s (%.2f%% of %d cells)", 
                    cluster, majority_celltype, percentage, total_cluster_cells)

        major_labels_dict[cluster] = majority_celltype
        detailed_labels_dict[cluster] = f"{majority_celltype}-(I={percentage:.1f}%)-C{cluster}"

    # Map the cluster-level annotations back to all cells in .obs
    adata.obs["celltype_major"] = adata.obs[args.leiden_key].map(major_labels_dict).astype("category")
    adata.obs["celltype_detailed"] = adata.obs[args.leiden_key].map(detailed_labels_dict).astype("category")

    logger.info("Generating UMAP visualization...")
    try:
        # Determine the palette to use (supports comma-separated custom list of colors)
        palette_to_use = [c.strip() for c in args.palette.split(",")] if "," in args.palette else args.palette
        logger.info("Using color palette: %s", palette_to_use)

        # We will create a 2-panel plot comparing original Leiden clusters with annotated major cell types
        fig, axes = plt.subplots(1, 2, figsize=(18, 8), dpi=300)
        
        # 1. Plot the Leiden clusters first (this populates f"{args.leiden_key}_colors" in adata.uns)
        sc.pl.umap(
            adata, 
            color=args.leiden_key, 
            ax=axes[0], 
            show=False, 
            legend_loc="on data", 
            legend_fontsize=8, 
            legend_fontweight="bold",
            palette=palette_to_use
        )
        axes[0].set_title("Leiden Clusters", fontsize=14, fontweight="bold")
        
        # Synchronize colors from leiden_scvi to celltype_major so colors don't change between panels
        leiden_colors_key = f"{args.leiden_key}_colors"
        if leiden_colors_key in adata.uns:
            leiden_colors = adata.uns[leiden_colors_key]
            leiden_categories = list(adata.obs[args.leiden_key].cat.categories)
            celltype_categories = list(adata.obs["celltype_major"].cat.categories)
            
            celltype_colors = []
            for ct in celltype_categories:
                # Find all clusters that have this celltype as their majority
                matching_clusters = [cl for cl, maj_ct in major_labels_dict.items() if maj_ct == ct]
                if matching_clusters:
                    # Take the first cluster
                    rep_cluster = sorted(matching_clusters)[0]
                    try:
                        rep_idx = leiden_categories.index(str(rep_cluster))
                        celltype_colors.append(leiden_colors[rep_idx])
                    except ValueError:
                        celltype_colors.append("#808080")
                else:
                    celltype_colors.append("#808080")
            
            adata.uns["celltype_major_colors"] = celltype_colors
            logger.info("Synchronized celltype_major colors with Leiden cluster colors.")

        # 2. Plot the cell-type annotations on the right panel with synchronized colors and legend on data
        sc.pl.umap(
            adata, 
            color="celltype_major", 
            ax=axes[1], 
            show=False,
            legend_loc="on data",
            legend_fontsize=7,
            legend_fontweight="bold"
        )
        axes[1].set_title("Assigned Cell-types (Majority Voting)", fontsize=14, fontweight="bold")
        
        plt.tight_layout()
        plt.savefig(args.output_plot, bbox_inches="tight")
        plt.close()
        logger.info("Successfully saved UMAP plot to: %s", args.output_plot)
    except Exception as e:
        logger.exception("Failed to generate/save visualization plot: %s", str(e))

    logger.info("Saving annotated AnnData to: %s", args.output_h5ad)
    try:
        # Sanitize obs strings for safe HDF5 serialization
        for col in ["predicted_celltype", "celltype_major", "celltype_detailed"]:
            if col in adata.obs.columns:
                # Convert to string first to safely unpack any Categoricals, then replace "nan"
                adata.obs[col] = adata.obs[col].astype(str).replace("nan", "")
                
        adata.write(args.output_h5ad)
        logger.info("=== Cluster Annotation Completed Successfully ===")
    except Exception as e:
        logger.exception("Failed to write annotated AnnData file: %s", str(e))
        sys.exit(1)


if __name__ == "__main__":
    main()
