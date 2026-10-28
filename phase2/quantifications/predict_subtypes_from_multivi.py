#!/usr/bin/env python3
"""
Predict granular subtype labels for ATAC (and unlabeled) cells within a specific
coarse cell type using MultiVI joint latent representations and AutoGluon.

Based on phase2/quantifications/predict_celltypes_from_multivi.ipynb.
"""

import argparse
import logging
import sys
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from autogluon.tabular import TabularDataset, TabularPredictor
from matplotlib.pyplot import rc_context
from pandas import DataFrame, Series, concat, read_csv
from seaborn import barplot
from sklearn.model_selection import train_test_split
from tabulate import tabulate

warnings.filterwarnings("ignore")

logger = logging.getLogger(__name__)

# Constants
DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"
DEFAULT_CELL_TYPE = "ExN CUX2"
DEFAULT_CELL_LABEL_COL = "cell_label"
DEFAULT_LABELS_SUBPATH = (
    "quants/resolution_selection/"
    "ExN_res_obs_annotated__ALLDATASETS_refresh_res_0.3_subtype.csv"
)
DEFAULT_BARCODE_COL = "X"
DEFAULT_LABEL_COL = "final_anno_subclust"


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Subset an annotated MultiVI AnnData object to a coarse cell type, "
            "integrate granular RNA subtype labels from a CSV, train an AutoGluon "
            "classifier on MultiVI latents, and predict granular subtypes for ATAC cells."
        )
    )
    parser.add_argument(
        "--project",
        type=str,
        default=DEFAULT_PROJECT,
        help=f"Project name used for default file prefixes (default: {DEFAULT_PROJECT}).",
    )
    parser.add_argument(
        "--work-dir",
        type=str,
        default=DEFAULT_WRK_DIR,
        help=f"Base working directory (default: {DEFAULT_WRK_DIR}).",
    )
    parser.add_argument(
        "--input-file",
        type=str,
        default=None,
        help=(
            "Path to input annotated MultiVI .h5ad file. "
            "Defaults to {work_dir}/quants/{project}.multivi.annotated.h5ad."
        ),
    )
    parser.add_argument(
        "--output-file",
        type=str,
        default=None,
        help=(
            "Path to output subsetted and granularly annotated .h5ad file. "
            "Defaults to {work_dir}/quants/{project}.{cell_type}.multivi.annotated.h5ad."
        ),
    )
    parser.add_argument(
        "--cell-type",
        type=str,
        default=DEFAULT_CELL_TYPE,
        help=f"Coarse cell type to subset to (default: '{DEFAULT_CELL_TYPE}').",
    )
    parser.add_argument(
        "--cell-label-col",
        type=str,
        default=DEFAULT_CELL_LABEL_COL,
        help=f"Column in adata.obs holding cell-type labels (default: '{DEFAULT_CELL_LABEL_COL}').",
    )
    parser.add_argument(
        "--labels-file",
        type=str,
        default=None,
        help=(
            "Path to CSV file with granular RNA cell subtype labels. "
            f"Defaults to {{work_dir}}/{DEFAULT_LABELS_SUBPATH}."
        ),
    )
    parser.add_argument(
        "--barcode-col",
        type=str,
        default=DEFAULT_BARCODE_COL,
        help=f"Column in --labels-file containing cell barcodes (default: '{DEFAULT_BARCODE_COL}').",
    )
    parser.add_argument(
        "--label-col",
        type=str,
        default=DEFAULT_LABEL_COL,
        help=f"Column in --labels-file containing granular subtype labels (default: '{DEFAULT_LABEL_COL}').",
    )
    parser.add_argument(
        "--model-dir",
        type=str,
        default=None,
        help=(
            "Directory path to save AutoGluon model artifacts. "
            "Defaults to {work_dir}/models/{project}_{cell_type}_trained_cellpred."
        ),
    )
    parser.add_argument(
        "--preset",
        type=str,
        default="good",
        help="AutoGluon model quality preset (default: 'good').",
    )
    parser.add_argument(
        "--eval-metric",
        type=str,
        default="mcc",
        help="AutoGluon evaluation metric (default: 'mcc').",
    )
    parser.add_argument(
        "--split-mode",
        type=str,
        default="random",
        choices=["random", "modality"],
        help=(
            "Train/test split strategy: 'random' for stratified random split, "
            "or 'modality' for expression=train and paired=test (default: 'random')."
        ),
    )
    parser.add_argument(
        "--test-size",
        type=float,
        default=0.3,
        help="Test set fraction when --split-mode is 'random' (default: 0.3).",
    )
    parser.add_argument(
        "--num-gpus",
        type=int,
        default=1,
        help="Number of GPUs to allocate to AutoGluon fit (default: 1).",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42).",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable verbose debug logging and table output.",
    )
    return parser.parse_args()


def peek_anndata(adata: AnnData, message: str = None, verbose: bool = False):
    if message:
        logger.info(message)
    logger.info(f"AnnData: {adata}")
    if verbose:
        print(tabulate(adata.obs.head(), headers="keys", tablefmt="psql"))
        print(tabulate(adata.var.head(), headers="keys", tablefmt="psql"))


def peek_dataframe(df: DataFrame, message: str = None, verbose: bool = False):
    if message:
        logger.info(message)
    logger.info(f"DataFrame shape: {df.shape}")
    if verbose:
        if len(df.columns) <= 15:
            print(tabulate(df.head(), headers="keys", tablefmt="psql"))
        else:
            print(tabulate(df.iloc[:5, :10], headers="keys", tablefmt="psql"))


def subset_by_cell_type(
    adata: AnnData, cell_type: str, cell_label_col: str, debug: bool = False
) -> AnnData:
    if cell_label_col not in adata.obs.columns:
        raise KeyError(
            f"Column '{cell_label_col}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    available_labels = adata.obs[cell_label_col].astype(str).unique().tolist()
    target_label = cell_type
    if target_label not in available_labels:
        alt_space = cell_type.replace("_", " ")
        alt_underscore = cell_type.replace(" ", "_")
        if alt_space in available_labels:
            target_label = alt_space
        elif alt_underscore in available_labels:
            target_label = alt_underscore
        else:
            raise ValueError(
                f"Requested cell type '{cell_type}' not found in adata.obs['{cell_label_col}']. "
                f"Available cell types: {sorted(available_labels)}"
            )

    logger.info(
        f"Subsetting AnnData to {cell_label_col} == '{target_label}' "
        f"(from {adata.n_obs} total cells)..."
    )
    adata_sub = adata[adata.obs[cell_label_col] == target_label].copy()
    peek_anndata(
        adata_sub,
        f"Subsetted AnnData to '{target_label}' ({adata_sub.n_obs} cells)",
        debug,
    )
    if "modality" in adata_sub.obs.columns:
        logger.info(
            "Modality breakdown in subset:\n%s",
            adata_sub.obs["modality"].value_counts().to_string(),
        )
    return adata_sub


def integrate_granular_labels(
    adata: AnnData,
    labels_file: Path,
    barcode_col: str,
    label_col: str,
    cell_label_col: str,
    debug: bool = False,
) -> AnnData:
    if not labels_file.exists():
        raise FileNotFoundError(f"Granular labels CSV file not found: {labels_file}")

    logger.info(f"Loading granular RNA labels from {labels_file}...")
    labels_df = read_csv(labels_file, low_memory=False)
    peek_dataframe(labels_df, "Loaded granular labels CSV", debug)

    if barcode_col not in labels_df.columns:
        raise KeyError(
            f"Barcode column '{barcode_col}' not found in {labels_file}. "
            f"Available columns: {list(labels_df.columns)}"
        )
    if label_col not in labels_df.columns:
        raise KeyError(
            f"Label column '{label_col}' not found in {labels_file}. "
            f"Available columns: {list(labels_df.columns)}"
        )

    # Deduplicate barcodes
    initial_len = len(labels_df)
    labels_df = labels_df.drop_duplicates(subset=[barcode_col], keep="first")
    if len(labels_df) < initial_len:
        logger.warning(
            f"Dropped {initial_len - len(labels_df)} duplicate barcodes from {labels_file}."
        )

    label_series = labels_df.set_index(barcode_col)[label_col]

    # Align barcodes: check direct match vs stripping MultiVI modality suffix
    obs_idx = pd.Index(adata.obs.index.astype(str))
    stripped_idx = obs_idx.str.replace(
        r"_(paired|expression|accessibility)$", "", regex=True
    )

    direct_matches = obs_idx.isin(label_series.index).sum()
    stripped_matches = stripped_idx.isin(label_series.index).sum()

    if stripped_matches >= direct_matches and stripped_matches > 0:
        logger.info(
            f"Matching barcodes after stripping modality suffixes "
            f"({stripped_matches}/{len(obs_idx)} cells matched)."
        )
        mapped_labels = label_series.reindex(stripped_idx)
        mapped_labels.index = adata.obs.index
    elif direct_matches > 0:
        logger.info(
            f"Matching barcodes directly against adata.obs.index "
            f"({direct_matches}/{len(obs_idx)} cells matched)."
        )
        mapped_labels = label_series.reindex(obs_idx)
        mapped_labels.index = adata.obs.index
    else:
        raise ValueError(
            f"No cell barcodes in adata.obs matched column '{barcode_col}' in {labels_file}."
        )

    # Assign mapped labels and set unmatched cells to 'Unknown'
    mapped_str = mapped_labels.fillna("Unknown").astype(str)
    mapped_str = mapped_str.replace({"": "Unknown", "nan": "Unknown", "NA": "Unknown"})

    adata.obs[cell_label_col] = pd.Categorical(mapped_str)
    adata.obs["label_prob"] = np.nan

    logger.info(
        "Granular label distribution after CSV integration:\n%s",
        adata.obs[cell_label_col].value_counts(dropna=False).to_string(),
    )
    if "modality" in adata.obs.columns:
        cross_tab = pd.crosstab(
            adata.obs[cell_label_col], adata.obs["modality"], dropna=False
        )
        logger.info("Granular labels by modality:\n%s", cross_tab.to_string())

    return adata


def eval_classifier_model(
    model: TabularPredictor,
    this_data: TabularDataset,
    target: str,
    dataset_name: str = "eval",
    verbose: bool = False,
) -> tuple[Series, DataFrame]:
    x_pred = model.predict(this_data)
    eval_results = model.evaluate_predictions(
        this_data[target], x_pred, detailed_report=True
    )
    acc = eval_results.get("accuracy")
    bal_acc = eval_results.get("balanced_accuracy")
    mcc = eval_results.get("mcc")
    logger.info(f"[{dataset_name}] {target=} accuracy={acc:.5f}")
    logger.info(f"[{dataset_name}] {target=} balanced_accuracy={bal_acc:.5f}")
    logger.info(
        f"[{dataset_name}] {target=} Matthews Correlation Coefficient: {mcc:.5f}"
    )

    ret_df = (
        DataFrame(eval_results.get("classification_report"))
        .transpose()
        .sort_values("f1-score", ascending=False)
    )
    ret_df["cell_type"] = ret_df.index.values
    ret_df = ret_df.loc[
        ~ret_df.cell_type.isin(["accuracy", "macro avg", "weighted avg"])
    ].reset_index(drop=True)

    if verbose:
        logger.info(
            f"[{dataset_name}] True label counts:\n{this_data[target].value_counts().to_string()}"
        )
        logger.info(
            f"[{dataset_name}] Predicted label counts:\n{x_pred.value_counts().to_string()}"
        )
        print(tabulate(ret_df, headers="keys", tablefmt="psql"))

    return x_pred, ret_df


def save_f1_barplot(scores_df: DataFrame, out_file: Path):
    with rc_context({"figure.figsize": (8, 8), "figure.dpi": 100}):
        plt.style.use("seaborn-v0_8-talk")
        fig, ax = plt.subplots(figsize=(8, 8), dpi=100)
        barplot(
            data=scores_df,
            x="cell_type",
            y="f1-score",
            hue="dataset",
            palette="colorblind",
            ax=ax,
        )
        ax.legend(title="dataset", bbox_to_anchor=(1.05, 1), loc="upper left")
        ax.tick_params(axis="x", rotation=90)
        ax.set_xlabel("Cell Types")
        ax.set_ylabel("F1 Score")
        fig.tight_layout()
        fig.savefig(out_file, bbox_inches="tight")
        plt.close(fig)
    logger.info(f"Saved F1 score barplot to {out_file}")


def save_counts_barplot(combined_counts: DataFrame, cell_label_col: str, out_file: Path):
    with rc_context({"figure.figsize": (8, 8), "figure.dpi": 100}):
        plt.style.use("seaborn-v0_8-talk")
        fig, ax = plt.subplots(figsize=(8, 8), dpi=100)
        barplot(
            data=combined_counts,
            x=cell_label_col,
            y="Cell_counts",
            hue="Dataset",
            palette="colorblind",
            ax=ax,
        )
        ax.legend(title="Dataset", bbox_to_anchor=(1.05, 1), loc="upper left")
        ax.tick_params(axis="x", rotation=90)
        ax.set_xlabel("Cell Types")
        ax.set_ylabel("Cell Counts")
        fig.tight_layout()
        fig.savefig(out_file, bbox_inches="tight")
        plt.close(fig)
    logger.info(f"Saved cell counts barplot to {out_file}")


def save_umap_plots(
    adata: AnnData, figs_dir: Path, prefix: str, cell_label_col: str
):
    if "X_umap" not in adata.obsm:
        logger.warning("X_umap not found in adata.obsm; skipping UMAP plots.")
        return

    with rc_context({"figure.figsize": (8, 8), "figure.dpi": 100}):
        plt.style.use("seaborn-v0_8-talk")
        for col in ["sample_id", "modality", cell_label_col]:
            if col in adata.obs.columns:
                fig = sc.pl.umap(adata, color=[col], show=False, return_fig=True)
                out_path = figs_dir / f"{prefix}.umap.{col}.png"
                fig.savefig(out_path, bbox_inches="tight")
                plt.close(fig)
                logger.info(f"Saved UMAP plot to {out_path}")

        if cell_label_col in adata.obs.columns:
            fig = sc.pl.umap(
                adata,
                color=[cell_label_col],
                frameon=False,
                legend_loc="on data",
                show=False,
                return_fig=True,
            )
            out_path = figs_dir / f"{prefix}.umap.{cell_label_col}_ondata.png"
            fig.savefig(out_path, bbox_inches="tight")
            plt.close(fig)
            logger.info(f"Saved on-data UMAP plot to {out_path}")


def main():
    args = parse_args()
    debug = args.debug

    work_dir = Path(args.work_dir)
    quants_dir = work_dir / "quants"
    models_dir = work_dir / "models"
    figs_dir = work_dir / "figures"
    logs_dir = work_dir / "logs"

    for d in [quants_dir, models_dir, figs_dir, logs_dir]:
        d.mkdir(parents=True, exist_ok=True)

    cell_type_clean = args.cell_type.replace(" ", "_")
    prefix = f"{args.project}.{cell_type_clean}"

    # Configure logging to file and stdout
    log_file = logs_dir / f"{prefix}_predict_subtypes.log"
    logging.basicConfig(
        level=logging.DEBUG if debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_file, mode="w"),
            logging.StreamHandler(sys.stdout),
        ],
        force=True,
    )
    logger.info(f"Command line: {' '.join(sys.argv)}")
    logger.info(f"Logging to {log_file}")

    # Resolve file paths
    in_h5ad_file = (
        Path(args.input_file)
        if args.input_file
        else quants_dir / f"{args.project}.multivi.annotated.h5ad"
    )
    out_h5ad_file = (
        Path(args.output_file)
        if args.output_file
        else quants_dir / f"{args.project}.{cell_type_clean}.multivi.annotated.h5ad"
    )
    labels_file = (
        Path(args.labels_file)
        if args.labels_file
        else work_dir / DEFAULT_LABELS_SUBPATH
    )
    trained_model_path = (
        Path(args.model_dir)
        if args.model_dir
        else models_dir / f"{args.project}_{cell_type_clean}_trained_cellpred"
    )

    logger.info(f"Input h5ad file: {in_h5ad_file}")
    logger.info(f"Output h5ad file: {out_h5ad_file}")
    logger.info(f"Granular labels CSV: {labels_file}")
    logger.info(f"AutoGluon model path: {trained_model_path}")

    if not in_h5ad_file.exists():
        raise FileNotFoundError(f"Input AnnData file not found: {in_h5ad_file}")

    # 1. Load the MultiVI AnnData and subset to the target coarse cell type
    logger.info(f"Loading MultiVI AnnData from {in_h5ad_file}...")
    adata = sc.read_h5ad(in_h5ad_file)
    peek_anndata(adata, "Loaded input MultiVI AnnData", debug)

    adata = subset_by_cell_type(adata, args.cell_type, args.cell_label_col, debug)

    if "MultiVI_latent" not in adata.obsm:
        raise KeyError(
            f"'MultiVI_latent' not found in adata.obsm. Available obsm keys: {list(adata.obsm.keys())}"
        )

    # 2. Integrate granular RNA subtype labels from CSV
    adata = integrate_granular_labels(
        adata=adata,
        labels_file=labels_file,
        barcode_col=args.barcode_col,
        label_col=args.label_col,
        cell_label_col=args.cell_label_col,
        debug=debug,
    )

    # 3. Split into train/test (known RNA labels) and inference (Unknown / ATAC) datasets
    target_var = args.cell_label_col
    adata_train_test = adata[adata.obs[target_var] != "Unknown"]
    adata_infer = adata[adata.obs[target_var] == "Unknown"]

    if adata_train_test.n_obs == 0:
        raise ValueError("No labeled cells (cell_label != 'Unknown') available for training.")
    if adata_infer.n_obs == 0:
        logger.warning("No 'Unknown' cells found for inference.")

    train_test_data = TabularDataset(adata_train_test.obsm["MultiVI_latent"])
    train_test_data[target_var] = adata_train_test.obs[target_var].astype(str).values
    train_test_data.index = adata_train_test.obs.index.values

    inference_data = TabularDataset(adata_infer.obsm["MultiVI_latent"])
    inference_data[target_var] = adata_infer.obs[target_var].astype(str).values
    inference_data.index = adata_infer.obs.index.values

    random_split = args.split_mode == "random"
    if random_split:
        X_train, X_test, _, _ = train_test_split(
            train_test_data.drop(columns=[target_var]),
            train_test_data[target_var],
            stratify=train_test_data[target_var],
            test_size=args.test_size,
            random_state=args.seed,
        )
        train_data = train_test_data.loc[X_train.index]
        test_data = train_test_data.loc[X_test.index]
    else:
        gex_temp = adata[
            (adata.obs.modality == "expression") & (adata.obs[target_var] != "Unknown")
        ].obs
        arc_temp = adata[
            (adata.obs.modality == "paired") & (adata.obs[target_var] != "Unknown")
        ].obs
        train_data = train_test_data.loc[gex_temp.index]
        test_data = train_test_data.loc[arc_temp.index]

    # Ensure no unlabeled cells in train or test
    train_data = train_data.loc[train_data[target_var] != "Unknown"]
    test_data = test_data.loc[test_data[target_var] != "Unknown"]

    peek_dataframe(train_data, "Training dataframe", debug)
    peek_dataframe(test_data, "Testing dataframe", debug)
    logger.info(f"Training set subtype counts:\n{train_data[target_var].value_counts().to_string()}")
    logger.info(f"Testing set subtype counts:\n{test_data[target_var].value_counts().to_string()}")
    logger.info(f"Inference set size: {len(inference_data)}")

    # 4. Initialize and fit AutoGluon TabularPredictor
    predictor = TabularPredictor(
        label=target_var,
        path=str(trained_model_path),
        verbosity=2,
        log_to_file=True,
        eval_metric=args.eval_metric,
    )
    logger.info(
        f"Fitting AutoGluon TabularPredictor (presets='{args.preset}', num_gpus={args.num_gpus})..."
    )
    predictor.fit(train_data, presets=args.preset, num_gpus=args.num_gpus)

    try:
        predictor.plot_ensemble_model()
    except Exception as e:
        logger.debug(f"Skipping plot_ensemble_model: {e}")

    # 5. Evaluate model on train and test sets
    logger.info("Evaluating predictions on training set...")
    train_predictions, train_scores = eval_classifier_model(
        predictor, train_data, target_var, dataset_name="train", verbose=debug
    )

    logger.info("Evaluating predictions on test set...")
    test_predictions, test_scores = eval_classifier_model(
        predictor, test_data, target_var, dataset_name="test", verbose=debug
    )

    logger.info(f"Best model: {predictor.model_best}")
    if debug:
        print(tabulate(predictor.leaderboard(), headers="keys", tablefmt="psql"))

    if random_split:
        train_scores["dataset"] = "train"
        test_scores["dataset"] = "test"
    else:
        train_scores["dataset"] = "expression"
        test_scores["dataset"] = "paired"

    scores_df = concat([test_scores, train_scores], ignore_index=True)
    save_f1_barplot(scores_df, figs_dir / f"{prefix}.f1_scores_train_test.png")

    # 6. Infer subtype labels and probabilities for Unknown (ATAC) cells
    logger.info(f"Predicting granular subtype labels for {len(inference_data)} unknown cells...")
    unknown_labels = predictor.predict(inference_data)
    unknown_label_probs = predictor.predict_proba(inference_data).max(axis=1)

    logger.info(
        f"Predicted inference subtype counts:\n{unknown_labels.value_counts().to_string()}"
    )

    # Plot cell counts by dataset split
    if random_split:
        col_names = ["Train", "Test", "Inference"]
    else:
        col_names = ["Expression", "Paired", "Accessibility"]

    combined_counts = concat(
        [
            train_data[target_var].value_counts(),
            test_data[target_var].value_counts(),
            unknown_labels.value_counts(),
        ],
        axis="columns",
    )
    combined_counts.columns = col_names
    combined_counts.index.name = target_var
    combined_counts = combined_counts.reset_index()
    combined_counts = combined_counts.melt(
        id_vars=[target_var],
        value_vars=col_names,
        var_name="Dataset",
        value_name="Cell_counts",
    ).sort_values(target_var)
    save_counts_barplot(
        combined_counts, target_var, figs_dir / f"{prefix}.cell_counts_split.png"
    )

    # 7. Integrate predicted labels and probabilities back into the AnnData object
    adata.obs[target_var] = adata.obs[target_var].astype(str)
    adata.obs.loc[unknown_labels.index, target_var] = unknown_labels.astype(str)
    adata.obs[target_var] = pd.Categorical(adata.obs[target_var])
    adata.obs[target_var] = adata.obs[target_var].cat.remove_unused_categories()

    adata.obs["label_prob"] = np.nan
    adata.obs.loc[unknown_labels.index, "label_prob"] = unknown_label_probs.values
    adata.obs["label_prob"] = adata.obs["label_prob"].astype(float)

    peek_anndata(adata, "AnnData with predicted granular subtype labels included", debug)
    logger.info(
        f"Final subtype counts:\n{adata.obs[target_var].value_counts().to_string()}"
    )
    logger.info(
        "Mean predicted label probability by subtype:\n%s",
        adata.obs.groupby(target_var, observed=True)["label_prob"].mean().to_string(),
    )

    # Save UMAPs and modality-level counts/evaluations when using random split
    save_umap_plots(adata, figs_dir, prefix, target_var)

    if random_split and "modality" in adata.obs.columns:
        mod_counts = concat(
            [
                adata.obs[adata.obs.modality == "expression"][target_var].value_counts(),
                adata.obs[adata.obs.modality == "paired"][target_var].value_counts(),
                adata.obs[adata.obs.modality == "accessibility"][target_var].value_counts(),
            ],
            axis="columns",
        )
        mod_counts.columns = ["Expression", "Paired", "Accessibility"]
        mod_counts.index.name = target_var
        mod_counts = mod_counts.reset_index().melt(
            id_vars=[target_var],
            value_vars=["Expression", "Paired", "Accessibility"],
            var_name="Dataset",
            value_name="Cell_counts",
        ).sort_values(target_var)
        save_counts_barplot(
            mod_counts, target_var, figs_dir / f"{prefix}.cell_counts_modality.png"
        )

        # Evaluate GEX and ARC separately as in the reference notebook
        gex_idx = adata.obs[adata.obs.modality == "expression"].index.intersection(
            train_test_data.index
        )
        arc_idx = adata.obs[adata.obs.modality == "paired"].index.intersection(
            train_test_data.index
        )
        if len(gex_idx) > 0 and len(arc_idx) > 0:
            gex_data = train_test_data.loc[gex_idx]
            arc_data = train_test_data.loc[arc_idx]
            logger.info("Evaluating predictions on GEX (expression) set...")
            _, gex_scores = eval_classifier_model(
                predictor, gex_data, target_var, dataset_name="expression", verbose=debug
            )
            logger.info("Evaluating predictions on ARC (paired) set...")
            _, arc_scores = eval_classifier_model(
                predictor, arc_data, target_var, dataset_name="paired", verbose=debug
            )
            gex_scores["dataset"] = "expression"
            arc_scores["dataset"] = "paired"
            assay_scores_df = concat([gex_scores, arc_scores], ignore_index=True)
            save_f1_barplot(
                assay_scores_df, figs_dir / f"{prefix}.f1_scores_modality.png"
            )

    # 8. Save the annotated subset AnnData object
    out_h5ad_file.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Saving annotated subset AnnData to {out_h5ad_file}...")
    adata.write(out_h5ad_file)
    logger.info(f"Successfully saved {out_h5ad_file}")


if __name__ == "__main__":
    main()
