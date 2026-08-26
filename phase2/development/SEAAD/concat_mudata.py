import os
import sys
import logging
import argparse
import warnings
import glob
import mudata as md
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

default_input_dir = (
    "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/public/seaad/doublet_det"
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Concatenate multiple per-sample MuData (.h5mu) files into a single unified MuData object."
    )
    parser.add_argument(
        "--input-dir",
        type=str,
        default=default_input_dir,
        help=f"Directory containing per-sample h5mu files to combine (default: {default_input_dir})",
    )
    parser.add_argument(
        "--file-pattern",
        type=str,
        default="*_combined.h5mu",
        help="Glob pattern to match individual h5mu files (default: '*_combined.h5mu')",
    )
    parser.add_argument(
        "--output-file",
        type=str,
        required=True,
        help="Path where the final combined h5mu file will be saved",
    )
    parser.add_argument(
        "--join",
        type=str,
        choices=["inner", "outer"],
        default="outer",
        help="How to handle variables (var) not present in all objects (default: 'outer')",
    )
    parser.add_argument(
        "--index-unique",
        type=str,
        default="-",
        help="Character to join suffix to ensure unique observation (barcode) names across samples (default: '-')",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    
    # Ensure parent output directory exists
    out_dir = os.path.dirname(os.path.abspath(args.output_file))
    os.makedirs(out_dir, exist_ok=True)
    
    # Configure logging to a file in the output directory in addition to terminal
    log_file_path = os.path.join(out_dir, "concat_mudata.log")
    file_handler = logging.FileHandler(log_file_path, mode="a")
    file_handler.setFormatter(logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    ))
    logging.getLogger().addHandler(file_handler)
    
    logger.info("Logging to terminal and file: %s", log_file_path)
    
    # Search for files matching the pattern in the input directory
    search_path = os.path.join(args.input_dir, args.file_pattern)
    logger.info("Searching for h5mu files using pattern: %s", search_path)
    
    h5mu_files = sorted(glob.glob(search_path))
    num_files = len(h5mu_files)
    
    if num_files == 0:
        logger.error("No matching h5mu files found!")
        sys.exit(1)
        
    logger.info("Found %d h5mu files to concatenate.", num_files)
    
    mdata_list = []
    keys = []
    
    # Iteratively load all MuData objects
    for i, file_path in enumerate(h5mu_files):
        logger.info("[%d/%d] Reading MuData: %s", i + 1, num_files, file_path)
        # Extract suffix ID from filename to use as the key for cell barcode uniqueness
        suffix_id = os.path.basename(file_path).split("_")[0]
        keys.append(suffix_id)
        
        try:
            mdata = md.read(file_path)
            logger.info("  Shapes - obs: %d, rna: %s, atac: %s", 
                        mdata.shape[0], 
                        mdata.mod["rna"].shape if "rna" in mdata.mod else "N/A", 
                        mdata.mod["atac"].shape if "atac" in mdata.mod else "N/A")
            mdata_list.append(mdata)
        except Exception as e:
            logger.exception("Failed to read %s: %s. Aborting.", file_path, str(e))
            sys.exit(1)
            
    # Handle single file edge case
    if len(mdata_list) == 1:
        logger.info("Only 1 MuData object found. No concatenation required. Saving directly to: %s", args.output_file)
        try:
            mdata_list[0].write(args.output_file)
            logger.info("Successfully saved MuData!")
        except Exception as e:
            logger.exception("Error saving single MuData file: %s", str(e))
            sys.exit(1)
        logger.info("=== Run Completed Successfully ===")
        return

    # Perform concatenation for multiple files
    logger.info("Concatenating %d MuData objects (join='%s', index_unique='%s')...", 
                len(mdata_list), args.join, args.index_unique)
                
    try:
        mdata_concat = md.concat(
            mdata_list,
            join=args.join,
            keys=keys,
            index_unique=args.index_unique
        )
        
        logger.info("Concatenation complete. Combined shape: %s", str(mdata_concat.shape))
        if "rna" in mdata_concat.mod:
            logger.info("Combined RNA shape: %s", str(mdata_concat.mod["rna"].shape))
        if "atac" in mdata_concat.mod:
            logger.info("Combined ATAC shape: %s", str(mdata_concat.mod["atac"].shape))
            
        # Write output to disk
        logger.info("Saving combined MuData to: %s", args.output_file)
        mdata_concat.write(args.output_file)
        logger.info("Successfully saved combined MuData!")
        
    except Exception as e:
        logger.exception("Error during concatenation or write: %s", str(e))
        sys.exit(1)
        
    logger.info("=== Concatenation Completed Successfully ===")


if __name__ == "__main__":
    main()
