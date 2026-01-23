from anndata import AnnData
from tabulate import tabulate
import logging
import argparse
from pathlib import Path
from scanpy import read_h5ad
from pandas import DataFrame

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Constants
DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"
VAR_MODAL_DICT = {"Gene Expression": "rna", "Peaks": "atac"}


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


def peek_dataframe(df: DataFrame, message: str = None, verbose: bool = False):
    if message:
        logger.info(message)
    print(f"{df.shape=}")
    if verbose:
        print(tabulate(df.head(), headers="keys", tablefmt="psql"))


def main():
    args = parse_args()
    debug = args.debug

    # Setup directories
    work_dir = Path(args.work_dir)
    quants_dir = work_dir / "quants"
    info_dir = work_dir / "sample_info"

    anndata_file = quants_dir / f"{args.project}.multivi.annotated.h5ad"

    # the final annotated anndata object should have all neccesary covariate information needed
    # load the annotated anndata
    adata = read_h5ad(anndata_file)
    peek_anndata(adata, "loaded annotated anndata from {anndata_file}", debug)

    # format sample covariates: sex, ancestry, age, (gex_pool or atac_pool), pmi, ph, smoker, bmi
    keep_terms = [
        "sample_id",
        "sex",
        "ancestry",
        "age",
        "gex_pool",
        "atac_pool",
        "pmi",
        "ph",
        "smoker",
        "bmi",
    ]
    covars_df = adata.obs[keep_terms].drop_duplicates().reset_index(drop=True)
    covars_df = covars_df.set_index("sample_id")
    peek_dataframe(covars_df, "copied covariates", debug)

    if debug:
        print(tabulate(covars_df.info(), headers="keys", tablefmt="psql"))
        print(covars_df.smoker.value_counts())
        print(covars_df.bmi.describe())

    # fill any missing covariate terms
    covars_df.loc[covars_df.smoker.isna(), "smoker"] = covars_df.smoker.mean().round(1)
    covars_df.loc[covars_df.bmi.isna(), "bmi"] = covars_df.bmi.mean().round(1)
    peek_dataframe(covars_df, "filled missing covariate values", debug)

    if debug:
        print(tabulate(covars_df.info(), headers="keys", tablefmt="psql"))
        print(covars_df.smoker.value_counts())
        print(covars_df.bmi.describe())

    # need to set the pool term based on modality being analyzed
    for modality in ["rna", "atac"]:
        if modality == "rna":
            covars_df["pool"] = covars_df.gex_pool
        elif modality == "atac":
            covars_df["pool"] = covars_df.atac_pool
        out_df = covars_df.drop(columns=["gex_pool", "atac_pool"])
        peek_dataframe(out_df, f"set pool covariate for {modality}", debug)

        # save the final covariate file
        out_file = info_dir / f"{args.project}.covariates.{modality}.csv"
        out_df.to_csv(out_file)
        logger.info(f"Saved {out_file} (Shape: {out_df.shape})")


if __name__ == "__main__":
    main()
