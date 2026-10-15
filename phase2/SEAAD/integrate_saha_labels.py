#!/usr/bin/env python3
"""
Integrate externally generated cell-type labels into a MuData object's .obs dataframe.
"""

import os
import sys
import logging
import argparse
import warnings
import pandas as pd
import mudata as md

# Silence mudata FutureWarnings regarding pull_on_update to keep console and logs clean
warnings.filterwarnings("ignore", category=FutureWarning, module="mudata")

# Configure logging to stdout
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Integrate external cell-type labels into the .obs dataframe of a MuData object."
    )
    parser.add_argument(
        "-i", "--input-file",
        type=str,
        required=True,
        help="Path to the input MuData (.h5mu) file.",
    )
    parser.add_argument(
        "-o", "--output-file",
        type=str,
        required=True,
        help="Path to save the modified MuData (.h5mu) file.",
    )
    parser.add_argument(
        "-l", "--labels-file",
        type=str,
        required=True,
        help="Path to the cell labels CSV file.",
    )
    parser.add_argument(
        "-b", "--barcode-col",
        type=str,
        required=True,
        help="Column name specifying cell barcodes in the labels CSV file.",
    )
    parser.add_argument(
        "-c", "--label-cols",
        type=str,
        nargs="+",
        required=True,
        help="List of one or more column names from the cell labels CSV to be added to the MuData .obs.",
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
        # Check if dtype is object, which can contain mixed types/NaNs
        if df[col].dtype == "object":
            # Check if column is boolean-like (contains True/False/NaN)
            non_null = df[col].dropna()
            if not non_null.empty and all(isinstance(val, bool) for col_val in non_null for val in [col_val]):
                logger.info("  Sanitizing boolean column: %s -> filling NaNs with False and casting to bool", col)
                df[col] = df[col].fillna(False).astype(bool)
            else:
                # Cast all other object columns to string to prevent generic h5py serialization failures
                logger.info("  Sanitizing object column: %s -> casting to string", col)
                df[col] = df[col].fillna("").astype(str)

    return df


def main():
    args = parse_args()

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

    logger.info("Logging to terminal and file: %s", log_file_path)

    logger.info("Loading MuData object from: %s", args.input_file)
    try:
        mdata = md.read(args.input_file)
    except Exception as e:
        logger.error("Failed to read input MuData file: %s", str(e))
        sys.exit(1)

    logger.info("Loading cell labels CSV from: %s", args.labels_file)
    try:
        labels_df = pd.read_csv(args.labels_file)
    except Exception as e:
        logger.error("Failed to read labels CSV file: %s", str(e))
        sys.exit(1)

    # Validate barcode column
    if args.barcode_col not in labels_df.columns:
        logger.error(
            "Barcode column '%s' not found in labels CSV columns: %s",
            args.barcode_col,
            list(labels_df.columns),
        )
        sys.exit(1)

    # Validate label columns
    missing_cols = [col for col in args.label_cols if col not in labels_df.columns]
    if missing_cols:
        logger.error(
            "The following requested label columns were not found in labels CSV: %s",
            missing_cols,
        )
        sys.exit(1)

    # Process and align labels
    logger.info("Aligning and integrating labels...")

    # Drop duplicates in the barcode column to ensure unique mapping
    initial_len = len(labels_df)
    labels_df = labels_df.drop_duplicates(subset=[args.barcode_col])
    final_len = len(labels_df)
    if initial_len != final_len:
        logger.warning(
            "Dropped %d duplicate cell barcodes from the labels file.",
            initial_len - final_len,
        )

    # Calculate barcode overlap statistics
    mdata_barcodes = set(mdata.obs.index)
    csv_barcodes = set(labels_df[args.barcode_col])
    overlap_barcodes = mdata_barcodes.intersection(csv_barcodes)

    num_mdata_cells = len(mdata.obs)
    num_csv_labels = len(labels_df)
    num_overlap = len(overlap_barcodes)

    pct_mdata_updated = (num_overlap / num_mdata_cells * 100) if num_mdata_cells > 0 else 0.0
    pct_csv_added = (num_overlap / num_csv_labels * 100) if num_csv_labels > 0 else 0.0

    logger.info("Alignment Statistics:")
    logger.info("  - Total cells in input MuData: %d", num_mdata_cells)
    logger.info("  - Total unique cell labels in CSV: %d", num_csv_labels)
    logger.info(
        "  - Cells in MuData updated with a label: %d (%.2f%%)",
        num_overlap,
        pct_mdata_updated,
    )
    logger.info(
        "  - Labels from CSV matched and added to MuData: %d (%.2f%%)",
        num_overlap,
        pct_csv_added,
    )

    # Reindex labels matching mdata.obs.index
    labels_subset = labels_df.set_index(args.barcode_col)[args.label_cols]

    # Add columns to mdata.obs
    for col in args.label_cols:
        series = labels_subset[col].reindex(mdata.obs.index)

        # Pre-convert object series with string values to category for efficiency if possible
        if series.dtype == "object":
            series = series.fillna("").astype(str)
            if series.nunique() < 0.5 * len(series):
                series = series.astype("category")

        mdata.obs[col] = series
        logger.info(
            "Added column '%s' to global MuData .obs. Non-null/non-empty values: %d/%d",
            col,
            mdata.obs[col].notna().sum(),
            len(mdata.obs),
        )

        # Log detailed summary stats for the newly integrated column
        if pd.api.types.is_numeric_dtype(series):
            stats = series.describe()
            logger.info("Summary statistics for numeric column '%s':\n%s", col, stats.to_string())
        else:
            counts = series.value_counts(dropna=False)
            logger.info("Value counts for categorical column '%s':\n%s", col, counts.to_string())

        # Also add to individual modalities if their indices match
        for mod_name, mod_obj in mdata.mod.items():
            mod_series = labels_subset[col].reindex(mod_obj.obs.index)
            if mod_series.dtype == "object":
                mod_series = mod_series.fillna("").astype(str)
                if mod_series.nunique() < 0.5 * len(mod_series):
                    mod_series = mod_series.astype("category")
            mod_obj.obs[col] = mod_series
            logger.info(
                "  -> Added to modality '%s' .obs. Non-null/non-empty values: %d/%d",
                mod_name,
                mod_obj.obs[col].notna().sum(),
                len(mod_obj.obs),
            )

    # Sanitize dataframes for compatibility
    logger.info("Sanitizing metadata dataframes for HDF5 compatibility...")
    mdata.obs = sanitize_dataframe(mdata.obs)
    for mod_name, mod_obj in mdata.mod.items():
        mod_obj.obs = sanitize_dataframe(mod_obj.obs)

    # Write modified MuData out
    logger.info("Saving modified MuData object to: %s", args.output_file)
    try:
        mdata.write(args.output_file)
        logger.info("Successfully saved modified MuData object.")
    except Exception as e:
        logger.error("Failed to write output MuData file: %s", str(e))
        sys.exit(1)


if __name__ == "__main__":
    main()
