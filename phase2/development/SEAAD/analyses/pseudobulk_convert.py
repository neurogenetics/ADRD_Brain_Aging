#!/usr/bin/env python3
"""
Convert single-cell profiles in a MuData object to pseudobulk profiles for disease analysis.
"""

import sys
import logging
import argparse
import warnings
from pathlib import Path
import pandas as pd
import numpy as np
import scanpy as sc
import mudata as md
from tabulate import tabulate
from anndata import AnnData

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)

# Thresholds
MIN_CELLS_RNA = 10
MIN_CELLS_ATAC = 20
MIN_CELLS_PROP = 0.30
TARGET_SUM_NORM = 1e6


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert single-cell MuData to pseudobulk profiles for disease analysis."
    )
    parser.add_argument(
        "-i",
        "--input-file",
        type=str,
        required=True,
        help="Path to the input MuData (.h5mu) file.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=str,
        default=None,
        help="Path to save the output parquet files. Defaults to the 'quants' folder inside --work-dir.",
    )
    parser.add_argument(
        "--project",
        type=str,
        default="seaad_ec_multiome",
        help="Project name used for file prefixes.",
    )
    parser.add_argument(
        "--work-dir",
        type=str,
        default="/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/seaad",
        help="Base working directory.",
    )
    parser.add_argument(
        "--cell-type-col",
        type=str,
        default="broad_cell_type",
        help="The obs column containing cell type labels (default: 'broad_cell_type').",
    )
    parser.add_argument(
        "--sample-col",
        type=str,
        default="sample_id",
        help="The obs column containing sample/donor IDs (default: 'sample_id').",
    )
    parser.add_argument(
        "--aggregate-type",
        type=str,
        default="sum",
        choices=["mean", "sum"],
        help="Type of pseudobulk aggregation to use (default: 'sum').",
    )
    parser.add_argument("--debug", action="store_true", help="Enable debug output.")
    parser.add_argument(
        "--exclude-ids",
        type=str,
        default="",
        help="Comma separated list of sample IDs to exclude.",
    )
    return parser.parse_args()


def peek_anndata(adata: AnnData, message: str = None, verbose: bool = False):
    if message:
        logger.info(message)
    print(adata)
    if verbose:
        print(tabulate(adata.obs.head(), headers="keys", tablefmt="psql"))
        print(tabulate(adata.var.head(), headers="keys", tablefmt="psql"))


def process_modality(
    ct_data: AnnData,
    non_pb_obs: pd.DataFrame,
    ct: str,
    modal_short: str,
    output_dir: Path,
    project_name: str,
    cell_type_col: str,
    sample_col: str,
    verbose: bool = False,
):
    """Filters, normalizes, and saves pseudobulk data for a specific modality."""
    if ct_data.n_vars == 0:
        logger.warning(f"No features found for {ct} - {modal_short}")
        return

    # Calculate donor-level cell counts to identify low-count samples
    adata_cell_info = non_pb_obs.loc[non_pb_obs[cell_type_col] == ct]
    donor_counts = adata_cell_info.groupby(sample_col, observed=True).size()

    min_cells = MIN_CELLS_RNA if modal_short == "rna" else MIN_CELLS_ATAC
    low_count_samples = list(donor_counts[donor_counts < min_cells].index.values)

    # Mask donors with low cell counts
    # The index in pseudobulk is usually "{sample}_{celltype}"
    ct_low_count_ids = [f"{x}_{ct}" for x in low_count_samples]
    donor_mask = ct_data.obs_names.isin(ct_low_count_ids)

    if verbose:
        logger.info(
            "low count masked donors: %s", ct_data.obs[donor_mask][sample_col].tolist()
        )

    # For low count samples, change their values to missing
    ct_data.X[donor_mask, :] = np.nan

    # Filter features
    min_cells_threshold = int(ct_data.n_obs * MIN_CELLS_PROP)
    pre_filter = ct_data.n_vars
    sc.pp.filter_genes(ct_data, min_cells=min_cells_threshold)
    post_filter = ct_data.n_vars

    logger.info(
        "[%s - %s] Filtered %d features. Masked %d low-count samples.",
        ct,
        modal_short,
        pre_filter - post_filter,
        len(low_count_samples),
    )

    # Convert to DataFrame
    df_modal = ct_data.to_df()
    # Clean index names: remove the cell type suffix added by aggregation if present
    df_modal.index = df_modal.index.str.removesuffix(f"_{ct}")

    # Save
    out_file = (
        output_dir / f"{project_name}.{ct.replace(' ', '_')}.{modal_short}.parquet"
    )
    df_modal.to_parquet(out_file)
    logger.info("Saved %s (Shape: %s)", out_file, df_modal.shape)


def main():
    args = parse_args()
    debug = args.debug
    aggr_type = args.aggregate_type

    # Setup directories
    work_dir = Path(args.work_dir)
    quants_dir = Path(args.output_dir) if args.output_dir else (work_dir / "quants")
    quants_dir.mkdir(parents=True, exist_ok=True)

    # Derive log file from output directory / project name
    log_dir = work_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file_path = log_dir / f"{args.project}_disease_pseudobulk_convert.log"

    # Set up FileHandler for logging to the derived log file
    file_handler = logging.FileHandler(log_file_path, mode="w")
    file_handler.setFormatter(
        logging.Formatter(
            "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )
    )
    logging.getLogger().addHandler(file_handler)

    logger.info("=== Starting Disease Pseudobulk Conversion ===")
    logger.info("Logging dynamically routed to file: %s", log_file_path)

    logger.info("Loading MuData object from: %s", args.input_file)
    try:
        mdata = md.read(args.input_file)
    except Exception as e:
        logger.error("Failed to read input MuData file: %s", str(e))
        sys.exit(1)

    # Exclude sample IDs if specified
    if args.exclude_ids:
        exclude_ids = [x.strip() for x in args.exclude_ids.split(",")]
        logger.info("Excluding sample IDs: %s", exclude_ids)
        # Apply filter to the main obs
        if args.sample_col in mdata.obs.columns:
            mdata = mdata[~mdata.obs[args.sample_col].isin(exclude_ids)].copy()
        else:
            logger.warning(
                "sample_col '%s' not found in mdata.obs. Skipping exclusion.",
                args.sample_col,
            )

    # Process each modality
    for modality in ["rna", "atac"]:
        if modality not in mdata.mod:
            logger.warning(
                "Modality '%s' not found in MuData object. Skipping.", modality
            )
            continue

        logger.info("Processing modality: %s", modality)
        adata_modal = mdata.mod[modality].copy()

        # Ensure cell type and sample id are present in modality obs
        for col in [args.cell_type_col, args.sample_col]:
            if col in mdata.obs.columns:
                adata_modal.obs[col] = mdata.obs[col].reindex(adata_modal.obs.index)

        # Validate that required columns exist
        if args.cell_type_col not in adata_modal.obs.columns:
            logger.error(
                "Cell type column '%s' not found in %s modality obs.",
                args.cell_type_col,
                modality,
            )
            sys.exit(1)
        if args.sample_col not in adata_modal.obs.columns:
            logger.error(
                "Sample column '%s' not found in %s modality obs.",
                args.sample_col,
                modality,
            )
            sys.exit(1)

        # Apply exclusion list at modality level just in case
        if args.exclude_ids:
            adata_modal = adata_modal[
                ~adata_modal.obs[args.sample_col].isin(exclude_ids)
            ].copy()

        if debug:
            peek_anndata(adata_modal, f"Subsetted AnnData {modality}", debug)

        # Perform aggregation
        if aggr_type == "mean":
            logger.info("Normalizing modality %s before mean aggregation...", modality)
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message="Some cells have zero counts",
                    category=UserWarning,
                )
                sc.pp.normalize_total(
                    adata_modal,
                    target_sum=TARGET_SUM_NORM,
                    exclude_highly_expressed=True,
                )
            sc.pp.log1p(adata_modal)

            logger.info(
                "Mean aggregating data by %s and %s...",
                args.sample_col,
                args.cell_type_col,
            )
            pb_adata = sc.get.aggregate(
                adata_modal, by=[args.sample_col, args.cell_type_col], func="mean"
            )
            pb_adata.X = pb_adata.layers["mean"].copy()

        elif aggr_type == "sum":
            logger.info(
                "Sum aggregating data by %s and %s...",
                args.sample_col,
                args.cell_type_col,
            )
            pb_adata = sc.get.aggregate(
                adata_modal, by=[args.sample_col, args.cell_type_col], func="sum"
            )
            pb_adata.X = pb_adata.layers["sum"].copy()

            logger.info("Normalizing modality %s after sum aggregation...", modality)
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message="Some cells have zero counts",
                    category=UserWarning,
                )
                sc.pp.normalize_total(
                    pb_adata, target_sum=TARGET_SUM_NORM, exclude_highly_expressed=True
                )
            sc.pp.log1p(pb_adata)

        if debug:
            peek_anndata(pb_adata, f"Transformed Pseudobulk {modality}", debug)

        # Process each cell type
        unique_cell_types = pb_adata.obs[args.cell_type_col].dropna().unique()
        logger.info(
            "Processing %d cell types for %s: %s",
            len(unique_cell_types),
            modality,
            list(unique_cell_types),
        )

        for ct in unique_cell_types:
            logger.info("--- Processing Cell Type: %s ---", ct)
            ct_data = pb_adata[pb_adata.obs[args.cell_type_col] == ct].copy()

            process_modality(
                ct_data=ct_data,
                non_pb_obs=adata_modal.obs,
                ct=ct,
                modal_short=modality,
                output_dir=quants_dir,
                project_name=args.project,
                cell_type_col=args.cell_type_col,
                sample_col=args.sample_col,
                verbose=debug,
            )

    logger.info("=== Disease Pseudobulk Conversion Finished Successfully ===")


if __name__ == "__main__":
    main()
