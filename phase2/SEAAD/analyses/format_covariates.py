#!/usr/bin/env python3
"""
Format covariate tables for use with data prep and regression analysis for disease study.
"""

import os
import sys
import logging
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import mudata as md
from tabulate import tabulate

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
        description="Format covariates from unified MuData object for disease study."
    )
    parser.add_argument(
        "-i",
        "--input-file",
        type=str,
        required=True,
        help="Path to the input MuData (.h5mu) file.",
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
    parser.add_argument("--debug", action="store_true", help="Enable debug output.")
    parser.add_argument(
        "--exclude-ids",
        type=str,
        default="",
        help="Comma separated list of sample IDs to exclude.",
    )
    return parser.parse_args()


def peek_dataframe(df: pd.DataFrame, message: str = None, verbose: bool = False):
    if message:
        logger.info(message)
    logger.info("DataFrame shape: %s", df.shape)
    if verbose:
        print(tabulate(df.head(), headers="keys", tablefmt="psql"))


def main():
    args = parse_args()
    debug = args.debug

    # Setup directories
    work_dir = Path(args.work_dir)
    info_dir = work_dir / "sample_info"
    info_dir.mkdir(parents=True, exist_ok=True)

    # Setup dynamic logging file
    log_dir = work_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file_path = log_dir / f"{args.project}_disease_format_covariates.log"

    # Set up FileHandler for logging to the derived log file
    file_handler = logging.FileHandler(log_file_path, mode="w")
    file_handler.setFormatter(
        logging.Formatter(
            "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )
    )
    logging.getLogger().addHandler(file_handler)

    logger.info("=== Starting Disease Covariates Formatting ===")
    logger.info("Logging dynamically routed to file: %s", log_file_path)

    logger.info("Loading MuData object from: %s", args.input_file)
    try:
        mdata = md.read(args.input_file)
    except Exception as e:
        logger.error("Failed to read input MuData file: %s", str(e))
        sys.exit(1)

    # Ensure required obs columns are present in global mdata.obs by pulling from modalities if needed
    required_cols = [
        "diagnosis",
        args.sample_col,
        args.cell_type_col,
        "sex",
        "race",
        "ageDeath",
        "PMI",
        "pH",
        "brainWeight",
    ]
    for col in required_cols:
        if col not in mdata.obs.columns:
            for mod_name, mod_obj in mdata.mod.items():
                if col in mod_obj.obs.columns:
                    mdata.obs[col] = mod_obj.obs[col].reindex(mdata.obs.index)
                    logger.info(
                        "  -> Pulled missing column '%s' from modality '%s' into global .obs",
                        col,
                        mod_name,
                    )
                    break

    # Exclude specified sample IDs
    if args.exclude_ids:
        exclude_ids = [x.strip() for x in args.exclude_ids.split(",")]
        logger.info("Excluding sample IDs: %s", exclude_ids)
        if args.sample_col in mdata.obs.columns:
            mdata = mdata[~mdata.obs[args.sample_col].isin(exclude_ids)].copy()
        else:
            logger.warning(
                "sample_col '%s' not found in mdata.obs, skipping exclusion.",
                args.sample_col,
            )

    # Feature Engineer 'dx' column from 'diagnosis'
    if "diagnosis" not in mdata.obs.columns:
        logger.error("Required 'diagnosis' column not found in MuData global .obs.")
        sys.exit(1)

    logger.info("Engineering 'dx' binary variable from 'diagnosis'...")
    mdata.obs["dx"] = np.nan
    mdata.obs.loc[mdata.obs["diagnosis"] == "control,no cognitive impairment", "dx"] = 0
    mdata.obs.loc[
        mdata.obs["diagnosis"].isin(
            [
                "Alzheimer disease,dementia",
                "Alzheimer disease,dementia,Lewy body disease",
            ]
        ),
        "dx",
    ] = 1
    # Convert 'dx' to nullable integer Int64 to prevent it from saving as float (0.0 / 1.0)
    mdata.obs["dx"] = mdata.obs["dx"].astype("Int64")

    # Extract sample covariates
    sample_col = args.sample_col
    keep_terms = [
        sample_col,
        "dx",
        "sex",
        "race",
        "ageDeath",
        "PMI",
        "pH",
        "brainWeight",
    ]

    # Filter terms actually present in the mudata obs
    present_terms = [col for col in keep_terms if col in mdata.obs.columns]
    missing_terms = [col for col in keep_terms if col not in mdata.obs.columns]
    if missing_terms:
        logger.warning(
            "The following requested covariates were not found in MuData obs: %s",
            missing_terms,
        )

    logger.info("Extracting and deduplicating covariates at sample level...")
    covars_df = (
        mdata.obs[present_terms]
        .drop_duplicates(subset=[sample_col])
        .reset_index(drop=True)
    )
    covars_df = covars_df.set_index(sample_col)
    peek_dataframe(covars_df, "Sample covariates extracted", debug)

    # Impute missing values for numeric covariates using column means to prevent sample dropping
    numeric_cols = ["ageDeath", "PMI", "pH", "brainWeight"]
    for col in numeric_cols:
        if col in covars_df.columns:
            # Handle special privacy-censored '90+' values for age
            col_str = covars_df[col].astype(str).str.strip()
            
            # Identify '90+' entries
            mask = col_str == "90+"
            if mask.any():
                # Generate unique random ages in range [91.0, 96.0]
                random_ages = np.random.uniform(91.0, 96.0, size=mask.sum()).round(3)
                col_str.loc[mask] = random_ages.astype(str)
                logger.info("Replaced %d privacy-censored '90+' values in '%s' with random ages in range [91, 96]", mask.sum(), col)
            
            # Force conversion to numeric (coercing non-numeric strings to NaNs)
            covars_df[col] = pd.to_numeric(col_str, errors="coerce")
            
            col_mean = covars_df[col].mean()
            if pd.notna(col_mean):
                missing_cnt = covars_df[col].isna().sum()
                if missing_cnt > 0:
                    logger.info(
                        "Imputing %d missing values in '%s' with mean %.3f",
                        missing_cnt,
                        col,
                        col_mean,
                    )
                    covars_df[col] = covars_df[col].fillna(col_mean).round(3)

    if debug:
        print(covars_df.describe())
        print(covars_df["dx"].value_counts(dropna=False))

    # Calculate modality-specific cell counts
    for modality in ["rna", "atac"]:
        if modality not in mdata.mod:
            logger.warning(
                "Modality '%s' not found in MuData object. Skipping count formatting.",
                modality,
            )
            continue

        logger.info("Formatting covariates file for modality: %s", modality)
        adata_modal = mdata.mod[modality].copy()

        # Ensure cell type and sample id are present in modality obs
        for col in [args.cell_type_col, args.sample_col]:
            if col in mdata.obs.columns:
                adata_modal.obs[col] = mdata.obs[col].reindex(adata_modal.obs.index)

        # Apply exclusion list at modality level just in case
        if args.exclude_ids:
            adata_modal = adata_modal[
                ~adata_modal.obs[args.sample_col].isin(exclude_ids)
            ].copy()

        # Calculate cell counts per sample per cell type
        donor_counts = adata_modal.obs.groupby(
            [args.cell_type_col, args.sample_col], observed=True
        ).size()
        donor_counts_df = donor_counts.unstack(level=0).fillna(0)
        donor_counts_df.columns.name = None
        donor_counts_df.index.name = None
        # Replace spaces/special chars and append suffix
        donor_counts_df.columns = [
            f"{x.replace(' ', '_').replace('/', '-')}_counts"
            for x in donor_counts_df.columns
        ]

        # Merge pivoted counts with covariates
        out_df = covars_df.merge(
            donor_counts_df, how="inner", left_index=True, right_index=True
        )
        peek_dataframe(out_df, f"Merged covariates for {modality}", debug)

        # Save covariates file
        out_file = info_dir / f"{args.project}.covariates.{modality}.csv"
        out_df.to_csv(out_file)
        logger.info("Saved %s (Shape: %s)", out_file, out_df.shape)

    logger.info("=== Disease Covariates Formatting Finished Successfully ===")


if __name__ == "__main__":
    main()
