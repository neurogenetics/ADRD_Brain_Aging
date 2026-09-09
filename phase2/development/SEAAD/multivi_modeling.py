#!/usr/bin/env python3
"""
Train a MultiVI model on a 10X Multiome MuData object, extract joint latent representations,
compute a joint UMAP, and generate UMAP visualization plots for specified obs metadata columns.
"""

import os
import sys
import logging
import argparse
import warnings
import pandas as pd
import mudata as md
import scanpy as sc
import scvi

# Silence warnings to keep logs clean
warnings.filterwarnings("ignore", category=FutureWarning, module="mudata")
warnings.filterwarnings("ignore", category=UserWarning, module="scvi")

# Setup default logger to stdout
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Train a MultiVI model on a Multiome MuData object and generate joint UMAPs."
    )
    parser.add_argument(
        "-i", "--input-file",
        type=str,
        required=True,
        help="Path to the input Multiome MuData (.h5mu) file.",
    )
    parser.add_argument(
        "-o", "--output-file",
        type=str,
        required=True,
        help="Path to save the modeled MuData (.h5mu) file.",
    )
    parser.add_argument(
        "-c", "--plot-cols",
        type=str,
        nargs="+",
        required=True,
        help="List of one or more obs column names to visualize on the joint UMAP.",
    )
    parser.add_argument(
        "-b", "--batch-key",
        type=str,
        default=None,
        help="Obs column to use as batch_key in MultiVI setup (optional).",
    )
    parser.add_argument(
        "--categorical-covariates",
        type=str,
        nargs="+",
        default=None,
        help="List of categorical covariate column names for MultiVI (optional).",
    )
    parser.add_argument(
        "--continuous-covariates",
        type=str,
        nargs="+",
        default=None,
        help="List of continuous covariate column names for MultiVI (optional).",
    )
    parser.add_argument(
        "--detect-hv-features",
        action="store_true",
        help="Whether to detect and subset to highly variable features before training.",
    )
    parser.add_argument(
        "--top-genes",
        type=int,
        default=4000,
        help="Number of top highly variable genes to keep if --detect-hv-features is used (default: 4000).",
    )
    parser.add_argument(
        "--top-peaks",
        type=int,
        default=20000,
        help="Number of top highly variable peaks to keep if --detect-hv-features is used (default: 20000).",
    )
    parser.add_argument(
        "-e", "--max-epochs",
        type=int,
        default=250,
        help="Maximum training epochs for MultiVI (default: 250).",
    )
    parser.add_argument(
        "-l", "--n-latent",
        type=int,
        default=30,
        help="Number of latent dimensions for MultiVI (default: 30).",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=10000,
        help="Minibatch size to use during training (default: 10000).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for scvi and training reproducibility (default: 42).",
    )
    return parser.parse_args()


def sanitize_dataframe(df):
    """
    Sanitizes a pandas DataFrame to ensure HDF5-compatible dtypes,
    preventing h5py write exceptions for object/mixed columns.
    """
    if df is None or df.empty:
        return df

    for col in df.columns:
        if df[col].dtype == "object":
            non_null = df[col].dropna()
            if not non_null.empty and all(
                isinstance(val, bool) for col_val in non_null for val in [col_val]
            ):
                logger.info(
                    "  Sanitizing boolean column: %s -> filling NaNs with False and casting to bool",
                    col,
                )
                df[col] = df[col].fillna(False).astype(bool)
            else:
                logger.info("  Sanitizing object column: %s -> casting to string", col)
                df[col] = df[col].fillna("").astype(str)

    return df


def main():
    args = parse_args()

    # Set training reproducibility seeds
    scvi.settings.seed = args.seed

    # Ensure parent output directory exists
    out_dir = os.path.dirname(os.path.abspath(args.output_file))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Derive log file path from output file
    out_base, _ = os.path.splitext(os.path.abspath(args.output_file))
    log_file_path = f"{out_base}.log"

    # Set up FileHandler for logging to the derived log file
    file_handler = logging.FileHandler(log_file_path, mode="w")
    file_handler.setFormatter(
        logging.Formatter(
            "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )
    )
    logging.getLogger().addHandler(file_handler)

    logger.info("=== Starting MultiVI Modeling Pipeline ===")
    logger.info("Logging dynamically routed to file: %s", log_file_path)

    logger.info("Loading MuData object from: %s", args.input_file)
    try:
        mdata = md.read(args.input_file)
    except Exception as e:
        logger.error("Failed to read input MuData file: %s", str(e))
        sys.exit(1)

    # Validate modalities
    if "rna" not in mdata.mod or "atac" not in mdata.mod:
        logger.error(
            "MultiVI requires both 'rna' and 'atac' modalities in the input MuData. Modalities found: %s",
            list(mdata.mod.keys()),
        )
        sys.exit(1)

    # Validate plot columns exist in global .obs
    missing_plot_cols = [col for col in args.plot_cols if col not in mdata.obs.columns]
    if missing_plot_cols:
        logger.error(
            "The following requested plot columns were not found in MuData .obs: %s",
            missing_plot_cols,
        )
        sys.exit(1)

    # Validate covariate keys if provided
    all_covs = []
    if args.batch_key:
        all_covs.append(args.batch_key)
    if args.categorical_covariates:
        all_covs.extend(args.categorical_covariates)
    if args.continuous_covariates:
        all_covs.extend(args.continuous_covariates)

    missing_covs = [
        cov for col in all_covs for cov in [col] if cov not in mdata.obs.columns
    ]
    if missing_covs:
        logger.error(
            "The following requested covariate columns were not found in MuData .obs: %s",
            missing_covs,
        )
        sys.exit(1)

    # Detect and subset to highly variable features if requested
    if args.detect_hv_features:
        logger.info("Detecting and subsetting to highly variable features...")
        try:
            # Subsetting genes (RNA modality)
            rna_mod = mdata.mod["rna"]
            logger.info("Current RNA features: %d", rna_mod.n_vars)
            hvg_flavor = "seurat"
            hvg_layer = None
            if "counts" in rna_mod.layers:
                hvg_flavor = "seurat_v3"
                hvg_layer = "counts"
                logger.info("  -> Found 'counts' layer in RNA modality, using seurat_v3 flavor.")
            
            try:
                sc.pp.highly_variable_genes(
                    rna_mod,
                    n_top_genes=min(args.top_genes, rna_mod.n_vars),
                    flavor=hvg_flavor,
                    layer=hvg_layer,
                    batch_key=args.batch_key if (args.batch_key and args.batch_key in rna_mod.obs.columns) else None,
                    subset=True,
                )
            except Exception as e:
                logger.warning("Scanpy highly_variable_genes failed with flavor '%s': %s. Falling back to robust 'seurat' flavor.", hvg_flavor, str(e))
                sc.pp.highly_variable_genes(
                    rna_mod,
                    n_top_genes=min(args.top_genes, rna_mod.n_vars),
                    flavor="seurat",
                    batch_key=args.batch_key if (args.batch_key and args.batch_key in rna_mod.obs.columns) else None,
                    subset=True,
                )
            logger.info("  -> Subsetted RNA modality to %d highly variable genes.", rna_mod.n_vars)

            # Subsetting peaks (ATAC modality)
            atac_mod = mdata.mod["atac"]
            logger.info("Current ATAC features: %d", atac_mod.n_vars)
            haf_flavor = "seurat"
            haf_layer = None
            if "counts" in atac_mod.layers:
                haf_flavor = "seurat_v3"
                haf_layer = "counts"
                logger.info("  -> Found 'counts' layer in ATAC modality, using seurat_v3 flavor.")
            
            try:
                sc.pp.highly_variable_genes(
                    atac_mod,
                    n_top_genes=min(args.top_peaks, atac_mod.n_vars),
                    flavor=haf_flavor,
                    layer=haf_layer,
                    batch_key=args.batch_key if (args.batch_key and args.batch_key in atac_mod.obs.columns) else None,
                    subset=True,
                )
            except Exception as e:
                logger.warning("Scanpy highly_variable_genes failed with flavor '%s': %s. Falling back to robust 'seurat' flavor.", haf_flavor, str(e))
                sc.pp.highly_variable_genes(
                    atac_mod,
                    n_top_genes=min(args.top_peaks, atac_mod.n_vars),
                    flavor="seurat",
                    batch_key=args.batch_key if (args.batch_key and args.batch_key in atac_mod.obs.columns) else None,
                    subset=True,
                )
            logger.info("  -> Subsetted ATAC modality to %d highly variable peaks.", atac_mod.n_vars)

            # Synchronize parent MuData object shape with modalities
            mdata.update()
            logger.info("MuData object successfully updated. New shape: %s", mdata.shape)

        except Exception as e:
            logger.error("Failed to detect highly variable features: %s", str(e))
            sys.exit(1)

    # Ensure required input layer structure is compatible
    # Check if "counts" layer exists in each modality to use for setup_mudata
    rna_layer_key = "counts" if "counts" in mdata.mod["rna"].layers else None
    atac_layer_key = "counts" if "counts" in mdata.mod["atac"].layers else None

    if rna_layer_key:
        logger.info("  -> Found 'counts' layer in RNA modality, will use for MultiVI model setup.")
    if atac_layer_key:
        logger.info("  -> Found 'counts' layer in ATAC modality, will use for MultiVI model setup.")

    # Set up MuData for MultiVI modeling
    logger.info("Setting up MuData object for MultiVI...")
    try:
        scvi.model.MULTIVI.setup_mudata(
            mdata,
            modalities={"rna_layer": "rna", "atac_layer": "atac"},
            rna_layer=rna_layer_key,
            atac_layer=atac_layer_key,
            batch_key=args.batch_key,
            categorical_covariate_keys=args.categorical_covariates,
            continuous_covariate_keys=args.continuous_covariates,
        )
    except Exception as e:
        logger.error("Failed to run MULTIVI.setup_mudata: %s", str(e))
        sys.exit(1)

    # Initialize MultiVI model
    logger.info("Initializing MULTIVI model with n_latent=%d...", args.n_latent)
    try:
        model = scvi.model.MULTIVI(mdata, n_latent=args.n_latent)
    except Exception as e:
        logger.error("Failed to initialize MultiVI model: %s", str(e))
        sys.exit(1)

    # Train model
    logger.info("Training MultiVI model for up to %d epochs (batch_size=%d)...", args.max_epochs, args.batch_size)
    try:
        model.train(max_epochs=args.max_epochs, batch_size=args.batch_size)
    except Exception as e:
        logger.error("Error during MultiVI model training: %s", str(e))
        sys.exit(1)

    # Retrieve joint latent representation and store in obsm
    logger.info("Extracting MultiVI joint latent representation...")
    try:
        mdata.obsm["X_multivi"] = model.get_latent_representation()
    except Exception as e:
        logger.error("Failed to retrieve latent representation: %s", str(e))
        sys.exit(1)

    # Compute joint UMAP using Scanpy
    logger.info("Computing nearest neighbors on MultiVI latent representation...")
    try:
        sc.pp.neighbors(mdata, use_rep="X_multivi")
        logger.info("Computing UMAP coordinates...")
        sc.tl.umap(mdata)
    except Exception as e:
        logger.error("Failed to compute UMAP embedding: %s", str(e))
        sys.exit(1)

    # Sanitize and write out MuData object
    logger.info("Sanitizing and saving modeled MuData object to: %s", args.output_file)
    try:
        mdata.obs = sanitize_dataframe(mdata.obs)
        for mod_name, mod_obj in mdata.mod.items():
            mod_obj.obs = sanitize_dataframe(mod_obj.obs)
        mdata.write(args.output_file)
        logger.info("Successfully saved modeled MuData object.")
    except Exception as e:
        logger.error("Failed to save MuData object: %s", str(e))
        sys.exit(1)

    # Generate and save UMAP visualizations
    logger.info("Generating UMAP plots for specified observation columns...")
    import matplotlib

    matplotlib.use(
        "Agg"
    )  # Set headless backend for saving figures without window display
    import matplotlib.pyplot as plt

    for col in args.plot_cols:
        plot_file_path = f"{out_base}_umap_{col}.png"
        logger.info("Generating plot for '%s' -> Saving to: %s", col, plot_file_path)
        try:
            # Dynamically configure arguments based on column data type
            plot_kwargs = {
                "color": col,
                "show": False,
            }
            is_numeric = pd.api.types.is_numeric_dtype(mdata.obs[col])
            if not is_numeric:
                plot_kwargs["legend_loc"] = "on data"
                plot_kwargs["palette"] = "tab20"
                logger.info("  -> Using 'on data' legend placement and 'tab20' color palette.")

            fig, ax = plt.subplots(figsize=(8, 8))
            plot_kwargs["ax"] = ax  # Associate with the newly created axes
            sc.pl.umap(mdata, **plot_kwargs)
            plt.title(f"MultiVI Joint UMAP - Colored by {col}")
            plt.savefig(plot_file_path, bbox_inches="tight", dpi=150)
            plt.close(fig)
            logger.info("Successfully saved plot for '%s'", col)
        except Exception as e:
            logger.error("Failed to generate UMAP plot for '%s': %s", col, str(e))

    logger.info("=== MultiVI Modeling Pipeline Finished Successfully ===")


if __name__ == "__main__":
    main()
