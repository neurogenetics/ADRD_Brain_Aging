import argparse
import logging
import os
import sys
import snapatac2 as snap

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
        description="Process chromatin accessibility peaks using SnapATAC2 and a reference BED file."
    )
    parser.add_argument(
        "--fragment-file",
        type=str,
        required=True,
        help="Path to 10x fragments tsv.gz file.",
    )
    parser.add_argument(
        "--cpeaks-bed",
        type=str,
        required=True,
        help="Path to reference peaks BED file (e.g., cPeaks_hg38.bed).",
    )
    parser.add_argument(
        "--output-file",
        type=str,
        required=True,
        help="Path to save the final cell-by-peak AnnData file (h5ad).",
    )
    parser.add_argument(
        "--genome",
        type=str,
        default="hg38",
        help="Reference genome name (e.g., 'hg38', 'mm10') or custom annotation file path (default: 'hg38').",
    )
    parser.add_argument(
        "--min-counts",
        type=int,
        default=1000,
        help="Minimum fragment counts per cell (default: 1000).",
    )
    parser.add_argument(
        "--min-tsse",
        type=float,
        default=5.0,
        help="Minimum TSS enrichment (TSSe) score (default: 5.0).",
    )
    parser.add_argument(
        "--max-mito",
        type=float,
        default=10.0,
        help="Maximum mitochondrial contamination fraction (default: 0.05).",
    )
    parser.add_argument(
        "--n-features",
        type=int,
        default=100000,
        help="Number of features to select for doublet calling/downstream (default: 100000).",
    )
    parser.add_argument(
        "--doublet-threshold",
        type=float,
        default=0.5,
        help="Probability threshold for doublet filtering (default: 0.5).",
    )
    parser.add_argument(
        "--sorted-by-barcode",
        action="store_true",
        default=False,
        help="Set to True if input fragments are already sorted by barcode (default: False).",
    )
    parser.add_argument(
        "--temp-dir",
        type=str,
        default="",
        help="Path to temporary directory for SnapATAC2 sorting/caching (defaults to the output directory).",
    )
    parser.add_argument(
        "--log-file",
        type=str,
        default="",
        help="Path to save execution log output in addition to stdout.",
    )
    return parser.parse_args()


def get_genome_robust(genome_str):
    """
    Robustly resolves reference genome annotations in SnapATAC2.
    """
    if os.path.exists(genome_str):
        logger.info(f"Using custom genome annotation file path: {genome_str}")
        return genome_str

    genome_name = genome_str.lower()
    if genome_name in ["hg38", "grch38"]:
        for attr in ["hg38", "HG38", "GRCh38"]:
            if hasattr(snap.genome, attr):
                resolved = getattr(snap.genome, attr)
                logger.info(f"Resolved reference genome hg38 via snap.genome.{attr}")
                return resolved
    elif genome_name in ["hg19", "grch37"]:
        for attr in ["hg19", "HG19", "GRCh37"]:
            if hasattr(snap.genome, attr):
                resolved = getattr(snap.genome, attr)
                logger.info(f"Resolved reference genome hg19 via snap.genome.{attr}")
                return resolved
    elif genome_name == "mm10":
        for attr in ["mm10", "MM10"]:
            if hasattr(snap.genome, attr):
                resolved = getattr(snap.genome, attr)
                logger.info(f"Resolved reference genome mm10 via snap.genome.{attr}")
                return resolved

    logger.warning(
        f"Could not find matching pre-defined genome for '{genome_str}'. Using raw value."
    )
    return genome_str


def import_fragments_robust(
    fragment_file, chrom_sizes, file, sorted_by_barcode, chrM, tempdir
):
    """
    Robustly imports fragment data using import_fragments or legacy import_data.
    """
    if hasattr(snap.pp, "import_fragments"):
        logger.info(f"Using snap.pp.import_fragments with tempdir: {tempdir}")
        return snap.pp.import_fragments(
            fragment_file=fragment_file,
            chrom_sizes=chrom_sizes,
            file=file,
            sorted_by_barcode=sorted_by_barcode,
            chrM=chrM,
            tempdir=tempdir,
        )
    elif hasattr(snap.pp, "import_data"):
        logger.info("Using legacy snap.pp.import_data...")
        try:
            return snap.pp.import_data(
                fragment_file=fragment_file,
                chrom_sizes=chrom_sizes,
                file=file,
                sorted_by_barcode=sorted_by_barcode,
                chrM=chrM,
                tempdir=tempdir,
            )
        except TypeError:
            return snap.pp.import_data(
                fragment_file=fragment_file,
                chrom_sizes=chrom_sizes,
                file=file,
                sorted_by_barcode=sorted_by_barcode,
                chrM=chrM,
            )
    else:
        raise AttributeError(
            "Neither import_fragments nor import_data was found in snapatac2.pp."
        )


def subset_robust(adata, mask):
    """
    Robustly subsets backed or in-memory AnnData objects.
    """
    if hasattr(adata, "subset"):
        logger.info("Using adata.subset for on-disk backed filtering...")
        adata.subset(obs_indices=mask, inplace=True)
        return adata
    else:
        logger.info("Using standard slicing copy for in-memory filtering...")
        return adata[mask, :].copy()


def make_peak_matrix_robust(adata, peak_file, output_file):
    """
    Robustly creates a cell-by-peak matrix using either peak_file or legacy use_rep.
    """
    logger.info(f"Creating peak matrix using reference BED: {peak_file}")
    try:
        return snap.pp.make_peak_matrix(
            adata,
            peak_file=peak_file,
            file=output_file,
        )
    except TypeError:
        logger.warning(
            "peak_file parameter failed, attempting legacy use_rep fallback..."
        )
        return snap.pp.make_peak_matrix(
            adata,
            use_rep=peak_file,
            file=output_file,
        )


def main():
    args = parse_args()

    # Convert paths to absolute paths immediately to prevent permission/relative path errors
    args.fragment_file = os.path.abspath(args.fragment_file)
    args.cpeaks_bed = os.path.abspath(args.cpeaks_bed)
    args.output_file = os.path.abspath(args.output_file)
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

    # Validate input file paths
    if not os.path.exists(args.fragment_file):
        logger.error(f"Input fragment file does not exist: {args.fragment_file}")
        raise FileNotFoundError(f"Fragment file not found: {args.fragment_file}")

    if not os.path.exists(args.cpeaks_bed):
        logger.error(f"Input BED file does not exist: {args.cpeaks_bed}")
        raise FileNotFoundError(f"BED file not found: {args.cpeaks_bed}")

    # Ensure final output directory exists
    out_dir = os.path.dirname(args.output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Resolve and ensure temporary directory for SnapATAC2 sorting/caching exists
    if args.temp_dir:
        temp_dir = os.path.abspath(args.temp_dir)
    else:
        temp_dir = out_dir if out_dir else os.getcwd()
    os.makedirs(temp_dir, exist_ok=True)
    logger.info(f"Using temporary directory for SnapATAC2 sorting: {temp_dir}")

    # Temporary intermediate h5ad path to avoid corrupting/mixing raw and peak matrices
    temp_h5ad = args.output_file + ".temp_fragments.h5ad"

    try:
        # Resolve genome object
        genome_obj = get_genome_robust(args.genome)

        # 1. Import fragment data
        logger.info(
            f"Importing fragments from {args.fragment_file} to intermediate disk file {temp_h5ad}..."
        )
        adata = import_fragments_robust(
            fragment_file=args.fragment_file,
            chrom_sizes=genome_obj,
            file=temp_h5ad,
            sorted_by_barcode=args.sorted_by_barcode,
            chrM=["chrM", "M"],
            tempdir=temp_dir,
        )
        logger.info(f"Successfully imported fragment AnnData: {adata}")

        # 2. Compute TSS enrichment score
        logger.info("Computing TSS enrichment (TSSe) scores...")
        snap.metrics.tsse(adata, genome_obj)

        # 3. Filter cells based on TSSe and unique fragment counts
        logger.info(
            f"Filtering cells (min_counts={args.min_counts}, min_tsse={args.min_tsse})..."
        )
        snap.pp.filter_cells(adata, min_counts=args.min_counts, min_tsse=args.min_tsse)

        # 4. Filter cells based on mitochondrial contamination fraction
        if "frac_mito" in adata.obs:
            logger.info(
                f"Filtering cells by mitochondrial contamination (max_mito={args.max_mito})..."
            )
            mito_mask = adata.obs["frac_mito"] < args.max_mito
            adata = subset_robust(adata, mito_mask)
        else:
            logger.warning(
                "frac_mito column not found in adata.obs; skipping mitochondrial filtering."
            )

        # 5. Create cell-by-cPeak matrix, saving it directly to final output file
        logger.info(f"Creating peak matrix at: {args.output_file}")
        peak_adata = make_peak_matrix_robust(
            adata=adata,
            peak_file=args.cpeaks_bed,
            output_file=args.output_file,
        )
        logger.info(f"Created peak matrix AnnData: {peak_adata}")

        # 6. Feature selection on peak matrix
        logger.info(f"Selecting top {args.n_features} features on peak matrix...")
        snap.pp.select_features(peak_adata, n_features=args.n_features)

        # 7. Doublet calling
        logger.info("Running Scrublet doublet calling on peak matrix...")
        snap.pp.scrublet(peak_adata, features="selected")

        # 8. Filter doublets
        logger.info(
            f"Filtering predicted doublets with probability threshold {args.doublet_threshold}..."
        )
        snap.pp.filter_doublets(
            peak_adata, probability_threshold=args.doublet_threshold
        )
        # Save the raw counts to a layer for MultiVI downstream
        peak_adata.layers["counts"] = peak_adata.X[:]

        # Close final h5ad file handles to write successfully
        logger.info("Closing file handles and flushing to disk...")
        if hasattr(peak_adata, "close"):
            peak_adata.close()
        if hasattr(adata, "close"):
            adata.close()

        logger.info("Workflow successfully completed!")

    finally:
        # Clean up temporary fragment AnnData file
        if os.path.exists(temp_h5ad):
            try:
                os.remove(temp_h5ad)
                logger.info(f"Cleaned up temporary file: {temp_h5ad}")
            except Exception as e:
                logger.warning(f"Could not clean up temporary file {temp_h5ad}: {e}")


if __name__ == "__main__":
    main()
