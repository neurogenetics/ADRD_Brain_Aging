import argparse
import logging
import sys
import warnings
from pathlib import Path

import scanpy as sc
from anndata import AnnData
from tabulate import tabulate
from pandas import DataFrame
import numpy as np

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Constants
DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"
VAR_MODAL_DICT = {"Gene Expression": "rna", "Peaks": "atac"}
MODAL_TYPES_DICT = {
    "rna": ["expression", "paired"],
    "atac": ["accessibility", "paired"],
}

# Thresholds
MIN_CELLS_RNA = 10
MIN_CELLS_ATAC = 20
MIN_CELLS_PROP = 0.30
TARGET_SUM_NORM = 1e6


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert single-cell AnnData to pseudobulk profiles."
    )
    parser.add_argument(
        "--project",
        type=str,
        default=DEFAULT_PROJECT,
        help="Project name used for file prefixes.",
    )
    parser.add_argument(
        "--work-dir",
        type=str,
        default=DEFAULT_WRK_DIR,
        help="Base working directory.",
    )
    parser.add_argument("--debug", action="store_true", help="Enable debug output.")
    return parser.parse_args()


def peek_anndata(adata: AnnData, message: str = None, verbose: bool = False):
    if message:
        logger.info(message)
    print(adata)
    if verbose:
        print(tabulate(adata.obs.head(), headers="keys", tablefmt="psql"))
        print(tabulate(adata.var.head(), headers="keys", tablefmt="psql"))


def load_and_prep_data(
    raw_path: Path, annot_path: Path, debug: bool = False
) -> tuple[AnnData, AnnData]:
    """Loads raw and annotated AnnData files and synchronizes them."""
    logger.info(f"Loading annotated AnnData: {annot_path}")
    if not annot_path.exists():
        logger.error(f"File not found: {annot_path}")
        sys.exit(1)
    adata_annot = sc.read(annot_path, cache=True)

    logger.info(f"Loading raw AnnData: {raw_path}")
    if not raw_path.exists():
        logger.error(f"File not found: {raw_path}")
        sys.exit(1)
    adata_raw = sc.read(raw_path)

    if debug:
        peek_anndata(adata_annot, "Loaded annotated AnnData", debug)
        peek_anndata(adata_raw, "Loaded raw AnnData", debug)

    # Subset and align raw data to annotated data
    logger.info("Aligning raw AnnData to annotated AnnData cells")
    # Using strict label transfer based on index intersection/ordering
    adata_raw = adata_raw[adata_annot.obs_names, :].copy()

    # Transfer necessary columns
    logger.info("Transferring metadata (cell_label, label_prob, sample_id)")
    adata_raw.obs["cell_label"] = adata_annot.obs["cell_label"]
    adata_raw.obs["label_prob"] = adata_annot.obs["label_prob"]

    # Drop unnecessary columns
    drop_columns = ["phase1_cluster", "phase1_celltype"]
    adata_raw.obs = adata_raw.obs.drop(columns=drop_columns, errors="ignore")

    if debug:
        peek_anndata(adata_raw, "Synced Raw AnnData", debug)
        print(adata_raw.obs.modality.value_counts())

    return adata_raw, adata_annot


def process_modality(
    ct_data: AnnData,
    non_pb_obs: DataFrame,
    ct: str,
    modal_short: str,
    output_dir: Path,
    project_name: str,
    verbose: bool = False,
):
    """Filters, normalizes, and saves pseudobulk data for a specific modality."""
    if ct_data.n_vars == 0:
        logger.warning(f"No features found for {ct} - {modal_short}")
        return

    # Calculate donor-level cell counts to identify low-count samples
    adata_cell_info = non_pb_obs.loc[non_pb_obs.cell_label == ct]
    donor_counts = adata_cell_info.groupby("sample_id", observed=True).size()

    min_cells = MIN_CELLS_RNA if modal_short == "rna" else MIN_CELLS_ATAC
    low_count_samples = list(donor_counts[donor_counts < min_cells].index.values)

    # Mask donors with low cell counts
    # The index in pseudobulk is usually "{sample}_{celltype}"
    ct_low_count_ids = [f"{x}_{ct}" for x in low_count_samples]
    donor_mask = ct_data.obs_names.isin(ct_low_count_ids)

    if verbose:
        logger.info(f"low count masked donors: {ct_data.obs[donor_mask]['sample_id']}")

    # Reset X to raw sums
    ct_data.X = ct_data.layers["sum"].copy()
    ct_data.X[donor_mask, :] = np.nan

    # Filter features
    min_cells_threshold = int(ct_data.n_obs * MIN_CELLS_PROP)
    pre_filter = ct_data.n_vars
    sc.pp.filter_genes(ct_data, min_cells=min_cells_threshold)
    post_filter = ct_data.n_vars

    logger.info(
        f"[{ct} - {modal_short}] Filtered {pre_filter - post_filter} features. "
        f"Masked {len(low_count_samples)} low-count samples."
    )

    # Normalize
    # RNA: CPM + Log1p | ATAC: RPM + Log1p
    # Suppress warning about zero-count cells (intentionally masked)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message="Some cells have zero counts", category=UserWarning
        )
        sc.pp.normalize_total(ct_data, target_sum=TARGET_SUM_NORM)

    sc.pp.log1p(ct_data)

    if verbose:
        logger.debug(ct_data, f"transformed pseudobulk {modal_short}", verbose)

    # Convert to DataFrame
    df_modal = ct_data.to_df()
    # Clean index names: remove the cell type suffix added by aggregation if present
    df_modal.index = df_modal.index.str.removesuffix(f"_{ct}")

    # Save
    out_file = (
        output_dir / f"{project_name}.{ct.replace(' ', '_')}.{modal_short}.parquet"
    )
    df_modal.to_parquet(out_file)
    logger.info(f"Saved {out_file} (Shape: {df_modal.shape})")


def main():
    args = parse_args()
    debug = args.debug

    # Setup directories
    work_dir = Path(args.work_dir)
    quants_dir = work_dir / "quants"
    logs_dir = work_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    log_filename = logs_dir / f"{args.project}_pseudobulk_convert.log"

    # Reconfigure logging to include file output
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        level=logging.DEBUG if debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_filename),
            logging.StreamHandler(sys.stdout)
        ]
    )

    logger.info(f"Command line: {' '.join(sys.argv)}")
    logger.info(f"Logging configured. Writing to {log_filename}")

    raw_file = quants_dir / f"{args.project}.raw.multivi_prep.h5ad"
    annot_file = quants_dir / f"{args.project}.multivi.annotated.h5ad"

    # Load and Sync Data
    adata_raw, _ = load_and_prep_data(raw_file, annot_file, debug)

    for modal_full, modal_short in VAR_MODAL_DICT.items():
        # subset for the modality
        logger.info(f"Subsetting anndata by modality: {modal_short}")
        modality_list = MODAL_TYPES_DICT.get(modal_short)
        adata_modal = adata_raw[
            adata_raw.obs.modality.isin(modality_list),
            adata_raw.var.modality == modal_full,
        ].copy()
        logger.debug(
            peek_anndata(adata_modal, f"subsetted anndata {modal_short}", debug)
        )
        # Pseudobulk Aggregation
        logger.info("Aggregating data by sample_id and cell_label...")
        pb_adata = sc.get.aggregate(
            adata_modal, by=["sample_id", "cell_label"], func="sum"
        )
        logger.debug(peek_anndata(pb_adata, "Pseudobulked AnnData", debug))

        # 3. Process each Cell Type
        unique_cell_types = pb_adata.obs["cell_label"].unique()
        logger.info(
            f"Processing {len(unique_cell_types)} cell types: {unique_cell_types}"
        )

        for ct in unique_cell_types:
            logger.info(f"--- Processing Cell Type: {ct} ---")
            ct_data = pb_adata[pb_adata.obs["cell_label"] == ct].copy()

            process_modality(
                ct_data=ct_data,
                non_pb_obs=adata_modal.obs,
                ct=ct,
                modal_short=modal_short,
                output_dir=quants_dir,
                project_name=args.project,
                verbose=debug,
            )


if __name__ == "__main__":
    main()
