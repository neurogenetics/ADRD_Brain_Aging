import os
import sys
import logging
import argparse
import warnings
import pandas as pd
import mudata as md
import scanpy as sc
import hdf5plugin

# Silence mudata FutureWarnings regarding pull_on_update to keep console and logs clean
warnings.filterwarnings("ignore", category=FutureWarning, module="mudata")

# Configure logging to stdout initially
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stdout,
)
logger = logging.getLogger(__name__)

# Pre-defined metadata and data directories as specified in the original file
suffix_map_file = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/seaad/SEA-AD_file_manifest_metadata.csv"
donor_info_file = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/seaad/SEA-AD_individual_metadata_harmonized.csv"
donor_info_columns = [
    "individualID",
    "sex",
    "race",
    "ageDeath",
    "PMI",
    "pH",
    "brainWeight",
    "diagnosis",
]
data_input_dir = (
    "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/seaad/doublet_det"
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine RNA and ATAC modalities into a single MuData object per sample, integrating donor information."
    )
    parser.add_argument(
        "--data-input-dir",
        type=str,
        default=data_input_dir,
        help=f"Directory containing filtered rna and atac h5ad files (default: {data_input_dir})",
    )
    parser.add_argument(
        "--suffix-map-file",
        type=str,
        default=suffix_map_file,
        help=f"Path to file manifest metadata CSV (default: {suffix_map_file})",
    )
    parser.add_argument(
        "--donor-info-file",
        type=str,
        default=donor_info_file,
        help=f"Path to individual metadata CSV (default: {donor_info_file})",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,
        help="Directory to save combined .h5mu files (default: same as data-input-dir)",
    )
    parser.add_argument(
        "--suffix",
        type=str,
        default=None,
        help="Process only this specific suffix ID instead of all matching files",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    in_dir = args.data_input_dir
    out_dir = args.output_dir if args.output_dir else in_dir

    os.makedirs(out_dir, exist_ok=True)

    # Configure logging to a file in the output directory in addition to terminal
    log_file_path = os.path.join(out_dir, "combine_modalities.log")
    file_handler = logging.FileHandler(log_file_path, mode="a")
    file_handler.setFormatter(
        logging.Formatter(
            "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )
    )
    logging.getLogger().addHandler(file_handler)

    logger.info("Logging to terminal and file: %s", log_file_path)

    logger.info("Reading suffix map file: %s", args.suffix_map_file)
    if not os.path.exists(args.suffix_map_file):
        logger.error("Suffix map file not found!")
        sys.exit(1)

    manifest_df = pd.read_csv(args.suffix_map_file).dropna(subset=["suffix"])

    # Clean the suffix representation to be string of integer (no trailing .0)
    def clean_suffix(x):
        s = str(x).strip()
        if s.endswith(".0"):
            return s[:-2]
        return s

    manifest_df["suffix_clean"] = manifest_df["suffix"].apply(clean_suffix)

    # Create mapping dictionary from suffix ID to donor ID (donor_id)
    suffix_to_donor = (
        manifest_df.drop_duplicates(subset=["suffix_clean"])
        .set_index("suffix_clean")["donor_id"]
        .to_dict()
    )

    logger.info("Reading donor info file: %s", args.donor_info_file)
    if not os.path.exists(args.donor_info_file):
        logger.error("Donor info file not found!")
        sys.exit(1)

    donor_df = pd.read_csv(args.donor_info_file)
    donor_df = donor_df.set_index("individualID", drop=False)

    # Identify suffixes to process
    if args.suffix:
        suffixes = [args.suffix]
        logger.info("Processing single specified suffix ID: %s", args.suffix)
    else:
        # Find all suffix IDs based on *_rna_filtered.h5ad files in in_dir
        if not os.path.exists(in_dir):
            logger.error("Data input directory does not exist: %s", in_dir)
            sys.exit(1)

        all_files = os.listdir(in_dir)
        suffixes = sorted(
            list(
                set(
                    [
                        f.split("_")[0]
                        for f in all_files
                        if f.endswith("_rna_filtered.h5ad")
                    ]
                )
            )
        )
        logger.info("Found %d suffixes with RNA files in %s", len(suffixes), in_dir)

    processed_count = 0
    failed_count = 0

    for suffix_id in suffixes:
        rna_file = f"{suffix_id}_rna_filtered.h5ad"
        atac_file = f"{suffix_id}_atac_filtered.h5ad"

        rna_path = os.path.join(in_dir, rna_file)
        atac_path = os.path.join(in_dir, atac_file)

        if not os.path.exists(rna_path):
            logger.warning(
                "RNA file not found for suffix %s: %s. Skipping.", suffix_id, rna_path
            )
            failed_count += 1
            continue

        if not os.path.exists(atac_path):
            logger.warning(
                "ATAC file not found for suffix %s: %s. Skipping.", suffix_id, atac_path
            )
            failed_count += 1
            continue

        logger.info("--- Processing Sample: %s ---", suffix_id)

        # Check if suffix ID has mapping
        donor_id = suffix_to_donor.get(suffix_id)
        if not donor_id:
            logger.warning(
                "No donor mapping found in suffix map for suffix ID %s. Skipping.",
                suffix_id,
            )
            failed_count += 1
            continue

        if donor_id not in donor_df.index:
            logger.warning(
                "Donor ID %s not found in donor metadata for suffix ID %s. Skipping.",
                donor_id,
                suffix_id,
            )
            failed_count += 1
            continue

        # Get donor record
        donor_row = donor_df.loc[donor_id]

        try:
            logger.info("Reading RNA data: %s", rna_path)
            adata_rna = sc.read_h5ad(rna_path)

            logger.info("Reading ATAC data: %s", atac_path)
            adata_atac = sc.read_h5ad(atac_path)

            logger.info(
                "RNA shape: %s, ATAC shape: %s", adata_rna.shape, adata_atac.shape
            )

            # 1. Intersect barcodes to ensure paired nuclei
            common_barcodes = adata_rna.obs_names.intersection(adata_atac.obs_names)

            num_common = len(common_barcodes)
            num_rna = adata_rna.shape[0]
            num_atac = adata_atac.shape[0]

            if num_rna > 0 and num_atac > 0:
                pct_rna = (num_common / num_rna) * 100
                pct_atac = (num_common / num_atac) * 100
                logger.info(
                    "Number of common barcodes: %d (%.2f%% of RNA cells, %.2f%% of ATAC cells)",
                    num_common,
                    pct_rna,
                    pct_atac,
                )
            else:
                logger.info("Number of common barcodes: %d", num_common)

            if num_common == 0:
                logger.warning(
                    "No common barcodes found between RNA and ATAC for suffix %s. Skipping.",
                    suffix_id,
                )
                failed_count += 1
                continue

            adata_rna_sub = adata_rna[common_barcodes].copy()
            adata_atac_sub = adata_atac[common_barcodes].copy()

            # 2. Combine into a MuData object first
            mdata = md.MuData({"rna": adata_rna_sub, "atac": adata_atac_sub})

            # 3 & 4. Populate the global MuData .obs directly from donor metadata
            for col in donor_info_columns:
                if col in donor_df.columns:
                    mdata.obs[col] = donor_row[col]

            mdata.obs["donor_id"] = donor_id
            mdata.obs["sample_id"] = suffix_id

            # 5. Save to disk
            output_file_name = f"{suffix_id}_combined.h5mu"
            output_file_path = os.path.join(out_dir, output_file_name)
            logger.info("Saving combined MuData to: %s", output_file_path)
            mdata.write(output_file_path)

            logger.info("Successfully processed suffix %s", suffix_id)
            processed_count += 1

        except Exception as e:
            logger.exception("Error processing suffix %s: %s", suffix_id, str(e))
            failed_count += 1

    logger.info("=== Run Completed ===")
    logger.info("Successfully processed: %d samples", processed_count)
    if failed_count > 0:
        logger.warning("Failed or skipped: %d samples", failed_count)


if __name__ == "__main__":
    main()
