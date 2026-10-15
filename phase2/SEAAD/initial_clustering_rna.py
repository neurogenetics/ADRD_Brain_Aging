import os
import sys
import logging
import argparse
import warnings
import mudata as md
import scanpy as sc
import scvi
import torch
import numpy as np
from sklearn.metrics import silhouette_score
import hdf5plugin

# Silence warnings to keep console and logs clean
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
        description="Perform initial clustering on the GEX (RNA) portion of a combined MuData object using scvi-tools."
    )
    parser.add_argument(
        "--input-file",
        type=str,
        required=True,
        help="Path to the input combined MuData (.h5mu) file",
    )
    parser.add_argument(
        "--output-file",
        type=str,
        required=True,
        help="Path to save the resulting RNA-only AnnData (.h5ad) file",
    )
    parser.add_argument(
        "--model-dir",
        type=str,
        default=None,
        help="Directory to save the trained scVI model (optional)",
    )
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
        help="Whether to detect and subset to highly variable features",
    )
    parser.add_argument(
        "--top-features-percent",
        type=float,
        default=0.15,
        help="Percentage of top highly variable genes to keep if --detect-hv-features is used (default: 0.15)",
    )
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

    # Configure logging to a file in the output directory
    out_dir = os.path.dirname(os.path.abspath(args.output_file))
    os.makedirs(out_dir, exist_ok=True)
    log_file_path = os.path.join(out_dir, "initial_clustering_rna.log")
    file_handler = logging.FileHandler(log_file_path, mode="a")
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    )
    logging.getLogger().addHandler(file_handler)

    logger.info("=== Starting RNA Initial Clustering ===")
    logger.info("Reading input MuData: %s", args.input_file)
    
    try:
        mdata = md.read(args.input_file)
    except Exception as e:
        logger.exception("Failed to read input file: %s", str(e))
        sys.exit(1)

    if "rna" not in mdata.mod:
        logger.error("Modality 'rna' not found in MuData object. Available modalities: %s", list(mdata.mod.keys()))
        sys.exit(1)

    logger.info("Extracting RNA modality...")
    adata_gex = mdata.mod["rna"].copy()

    # Copy global .obs to RNA .obs to ensure covariates (like donor_id, sex, race) are present
    logger.info("Syncing metadata from global MuData obs to RNA obs...")
    for col in mdata.obs.columns:
        if col not in adata_gex.obs.columns:
            adata_gex.obs[col] = mdata.obs[col]

    logger.info("Initial RNA shape: %s", adata_gex.shape)

    # 1. Quality Control
    logger.info("Calculating QC metrics...")
    adata_gex.var["mt"] = adata_gex.var_names.str.startswith("MT-")
    adata_gex.var["ribo"] = adata_gex.var_names.str.startswith(("RPS", "RPL"))
    adata_gex.var["hb"] = adata_gex.var_names.str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(
        adata_gex, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )

    logger.info("Filtering cells with >%.2f%% mitochondrial counts...", args.max_mito_percent)
    adata_gex = adata_gex[adata_gex.obs.pct_counts_mt < args.max_mito_percent, :].copy()

    logger.info("Filtering cells with <200 genes...")
    sc.pp.filter_cells(adata_gex, min_genes=200)

    min_cells = int(adata_gex.shape[0] * args.min_cell_percent)
    logger.info("Filtering genes expressed in <%.2f%% of cells (min_cells=%d)...", args.min_cell_percent * 100, min_cells)
    sc.pp.filter_genes(adata_gex, min_cells=min_cells)
    
    logger.info("Shape after QC and filtering: %s", adata_gex.shape)

    # 2. Highly Variable Features
    if args.detect_hv_features:
        n_top_genes = int(adata_gex.n_vars * args.top_features_percent)
        logger.info("Detecting top %d highly variable genes...", n_top_genes)
        
        # seurat_v3 flavor requires counts
        try:
            sc.pp.highly_variable_genes(
                adata_gex,
                n_top_genes=n_top_genes,
                batch_key=args.batch_key if args.batch_key in adata_gex.obs else None,
                flavor="seurat_v3",
                subset=True
            )
            logger.info("Shape after subsetting HVGs: %s", adata_gex.shape)
        except Exception as e:
            logger.warning("Failed to calculate highly variable genes (possibly due to non-integer counts). Skipping HVG step. Error: %s", str(e))

    # 3. scVI Setup
    cat_covs = [x.strip() for x in args.categorical_covariates.split(",") if x.strip()] if args.categorical_covariates else None
    cont_covs = [x.strip() for x in args.continuous_covariates.split(",") if x.strip()] if args.continuous_covariates else None

    # Validate covariates exist
    if args.batch_key and args.batch_key not in adata_gex.obs.columns:
        logger.warning("Batch key '%s' not found in obs. Proceeding without batch key.", args.batch_key)
        batch_key = None
    else:
        batch_key = args.batch_key

    valid_cat_covs = []
    if cat_covs:
        for c in cat_covs:
            if c in adata_gex.obs.columns:
                valid_cat_covs.append(c)
                # Ensure they are categorical
                adata_gex.obs[c] = adata_gex.obs[c].astype("category")
            else:
                logger.warning("Categorical covariate '%s' not found. Skipping.", c)

    valid_cont_covs = []
    if cont_covs:
        for c in cont_covs:
            if c in adata_gex.obs.columns:
                valid_cont_covs.append(c)
            else:
                logger.warning("Continuous covariate '%s' not found. Skipping.", c)

    logger.info("Setting up scVI AnnData...")
    logger.info("Batch key: %s", batch_key)
    logger.info("Categorical covariates: %s", valid_cat_covs)
    logger.info("Continuous covariates: %s", valid_cont_covs)

    scvi.model.SCVI.setup_anndata(
        adata_gex,
        batch_key=batch_key,
        categorical_covariate_keys=valid_cat_covs if valid_cat_covs else None,
        continuous_covariate_keys=valid_cont_covs if valid_cont_covs else None,
    )

    # 4. Train Model
    logger.info("Initializing scVI model...")
    model = scvi.model.SCVI(adata_gex, n_latent=args.n_latent, n_layers=args.n_layers)
    
    logger.info("Training scVI model (batch_size=%d)...", args.batch_size)
    model.train(batch_size=args.batch_size)

    if args.model_dir:
        logger.info("Saving trained scVI model to: %s", args.model_dir)
        os.makedirs(args.model_dir, exist_ok=True)
        model.save(args.model_dir, overwrite=True)

    # 5. Extract Latent Representation
    logger.info("Extracting latent representation...")
    adata_gex.obsm["scvi_latent"] = model.get_latent_representation()

    # 6. Neighborhood graph and UMAP
    logger.info("Computing neighborhood graph...")
    sc.pp.neighbors(adata_gex, use_rep="scvi_latent")
    
    logger.info("Computing UMAP...")
    sc.tl.umap(adata_gex)

    # 7. Clustering (Auto-search for best resolution based on silhouette score)
    logger.info("Auto-searching for best Leiden resolution (0.3 to 0.7)...")
    resolutions_to_try = np.arange(0.3, 0.8, 0.1)
    
    best_res = 0.3
    largest_score = -1.0
    leiden_key = "leiden_scvi"

    # Determine if igraph flavor is available
    leiden_kwargs = {}
    try:
        import igraph
        import leidenalg
        # Modern scanpy flavor defaults are changing, but let's safely try igraph flavor if requested in notebook
        leiden_kwargs = {"flavor": "igraph", "n_iterations": 2}
    except ImportError:
        logger.warning("igraph or leidenalg not found. Using default scanpy leiden settings.")

    for res in resolutions_to_try:
        res = round(res, 2)
        logger.info("  Testing resolution: %.2f", res)
        try:
            sc.tl.leiden(adata_gex, key_added=leiden_key, resolution=res, **leiden_kwargs)
            score = silhouette_score(adata_gex.obsm["scvi_latent"], adata_gex.obs[leiden_key])
            num_clusters = adata_gex.obs[leiden_key].nunique()
            logger.info("    -> Silhouette score: %.3f (Clusters: %d)", score, num_clusters)
            
            if score > largest_score:
                largest_score = score
                best_res = res
        except Exception as e:
            logger.warning("    -> Failed at resolution %.2f: %s", res, str(e))

    logger.info("Best resolution found: %.2f (Score: %.3f)", best_res, largest_score)
    
    # Run final clustering with the best resolution
    logger.info("Running final Leiden clustering at resolution %.2f...", best_res)
    sc.tl.leiden(adata_gex, key_added=leiden_key, resolution=best_res, **leiden_kwargs)

    # 8. Save Output
    logger.info("Saving processed RNA AnnData to: %s", args.output_file)
    try:
        # Sanitize obs to prevent h5py write errors with mixed types (like concat_mudata does)
        for col in adata_gex.obs.columns:
            if adata_gex.obs[col].dtype == "object":
                adata_gex.obs[col] = adata_gex.obs[col].fillna("").astype(str)
                
        adata_gex.write(args.output_file)
        logger.info("=== Run Completed Successfully ===")
    except Exception as e:
        logger.exception("Failed to write output file: %s", str(e))
        sys.exit(1)

if __name__ == "__main__":
    main()
