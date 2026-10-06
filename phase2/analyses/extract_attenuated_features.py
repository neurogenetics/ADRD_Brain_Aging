#!/usr/bin/env python3
"""
Extract and inspect attenuated features from cis-conditioned regression analysis.

Identifies endogenous features (e.g. genes) whose association with a target variable
(e.g., age or diagnosis/dx) is attenuated (i.e. loses significance, p > alpha)
when conditioned on a proximal cis-regulatory exogenous feature (e.g. ATAC peak).
"""

import argparse
import logging
from pathlib import Path
import sys

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract and inspect attenuated features from cis-conditioned regression analysis."
    )
    parser.add_argument(
        "--project",
        type=str,
        default=DEFAULT_PROJECT,
        help=f"Project prefix name (default: '{DEFAULT_PROJECT}').",
    )
    parser.add_argument(
        "--work-dir",
        type=str,
        default=DEFAULT_WRK_DIR,
        help=f"Base working directory containing results/ and figures/ (default: '{DEFAULT_WRK_DIR}').",
    )
    parser.add_argument(
        "--target-variable",
        type=str,
        default="age",
        help="Primary target variable column (e.g. 'age', 'dx') (default: 'age').",
    )
    parser.add_argument(
        "--endo-modality",
        type=str,
        default="rna",
        help="Endogenous modality (default: 'rna').",
    )
    parser.add_argument(
        "--exog-modality",
        type=str,
        default="atac",
        help="Exogenous modality (default: 'atac').",
    )
    parser.add_argument(
        "--regression-type",
        type=str,
        default="wls",
        help="Regression model type used (default: 'wls').",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Alpha significance threshold for exposure; exposure_pval > alpha indicates loss of significance (default: 0.05).",
    )
    parser.add_argument(
        "--output-suffix",
        type=str,
        default=None,
        help="Optional suffix appended to filename runs (e.g. 'dxage').",
    )
    parser.add_argument(
        "--cell-types",
        type=str,
        nargs="+",
        default=None,
        help="Optional list of specific cell types to filter for (e.g. ExN_CUX2 ExN_SEMA3E).",
    )
    parser.add_argument(
        "--cell-type-map",
        type=str,
        default=None,
        help="Comma-separated mapping of cell-type names (e.g. 'Source1:Target1,Source2:Target2').",
    )
    parser.add_argument(
        "--conditioned-file",
        type=str,
        default=None,
        help="Explicit path to conditioned results CSV (overrides auto-resolution).",
    )
    parser.add_argument(
        "--endo-baseline-file",
        type=str,
        default=None,
        help="Explicit path to baseline unconditioned endogenous FDR-filtered CSV (overrides auto-resolution).",
    )
    parser.add_argument(
        "--cis-file",
        type=str,
        default=None,
        help="Explicit path to cis-correlation results CSV to incorporate cis correlation statistics.",
    )
    parser.add_argument(
        "--output-file",
        "--out-file",
        type=str,
        default=None,
        help="Optional path to write the extracted attenuated features CSV. If omitted, prints table to stdout.",
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Auto-save the extracted attenuated table into <work-dir>/results/ if --output-file is not specified.",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Enable debug logging output.",
    )
    return parser.parse_args()


def parse_cell_type_map(map_str: str) -> dict:
    if not map_str:
        return {}
    mapping = {}
    for item in map_str.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" in item:
            k, v = item.split(":", 1)
            mapping[k.strip()] = v.strip()
    return mapping


def find_file(candidates: list) -> Path:
    for cand in candidates:
        if cand is not None and Path(cand).exists():
            return Path(cand)
    return None


def main():
    args = parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler()],
        force=True,
    )
    logger.info(f"Command line: {' '.join(sys.argv)}")

    work_dir = Path(args.work_dir)
    results_dir = work_dir / "results"

    # 1. Resolve Conditioned Regression File
    if args.conditioned_file:
        cond_file = Path(args.conditioned_file)
    else:
        cond_candidates = []
        if args.output_suffix:
            cond_candidates.append(
                results_dir
                / f"{args.project}.{args.endo_modality}-{args.exog_modality}.all_celltypes.{args.regression_type}.conditioned.{args.target_variable}.{args.output_suffix}.csv"
            )
        cond_candidates.append(
            results_dir
            / f"{args.project}.{args.endo_modality}-{args.exog_modality}.all_celltypes.{args.regression_type}.conditioned.{args.target_variable}.csv"
        )
        cond_candidates.append(
            results_dir
            / f"{args.project}.{args.endo_modality}-{args.exog_modality}.all_celltypes.{args.regression_type}.conditioned.csv"
        )
        cond_file = find_file(cond_candidates)

    if not cond_file or not cond_file.exists():
        logger.error(
            f"Conditioned results file could not be found. Checked candidates:\n"
            + "\n".join(str(c) for c in (cond_candidates if not args.conditioned_file else [args.conditioned_file]))
        )
        sys.exit(1)

    logger.info(f"Loading conditioned regression results from: {cond_file}")
    cond_df = pd.read_csv(cond_file)

    # 2. Resolve Baseline Endogenous FDR File
    if args.endo_baseline_file:
        endo_file = Path(args.endo_baseline_file)
    else:
        endo_candidates = [
            results_dir
            / f"{args.project}.{args.endo_modality}.all_celltypes.{args.regression_type}_fdr_filtered.{args.target_variable}.csv",
            results_dir
            / f"{args.project}.{args.endo_modality}.all_celltypes.{args.regression_type}_fdr_filtered.csv",
        ]
        endo_file = find_file(endo_candidates)

    endo_df = None
    if endo_file and endo_file.exists():
        logger.info(f"Loading baseline endogenous FDR results from: {endo_file}")
        endo_df = pd.read_csv(endo_file)
    else:
        logger.warning(
            "Baseline endogenous results file not found; proceeding without baseline comparison."
        )

    # 3. Resolve Cis Correlation File (Optional)
    if args.cis_file:
        cis_file = Path(args.cis_file)
    else:
        cis_candidates = [
            results_dir
            / f"{args.project}.{args.endo_modality}-{args.exog_modality}.all_celltypes.{args.regression_type}.{args.target_variable}.cis.csv",
            results_dir
            / f"{args.project}.{args.endo_modality}-{args.exog_modality}.all_celltypes.{args.regression_type}.cis.csv",
        ]
        cis_file = find_file(cis_candidates)

    cis_df = None
    if cis_file and cis_file.exists():
        logger.info(f"Loading cis correlation results from: {cis_file}")
        cis_df = pd.read_csv(cis_file)

    # Apply cell type map if provided
    cell_type_map = parse_cell_type_map(args.cell_type_map)
    if cell_type_map:
        cond_df["tissue"] = cond_df["tissue"].replace(cell_type_map)
        if endo_df is not None and "tissue" in endo_df.columns:
            endo_df["tissue"] = endo_df["tissue"].replace(cell_type_map)
        if cis_df is not None and "tissue" in cis_df.columns:
            cis_df["tissue"] = cis_df["tissue"].replace(cell_type_map)

    # Filter for attenuated features (loss of significance: exposure_pval > alpha)
    attenuated_mask = cond_df["exposure_pval"] > args.alpha
    attenuated_df = cond_df[attenuated_mask].copy()

    # Optional filter by specified cell types
    if args.cell_types:
        attenuated_df = attenuated_df[attenuated_df["tissue"].isin(args.cell_types)].copy()

    if attenuated_df.empty:
        logger.info(
            f"No attenuated features found with exposure_pval > {args.alpha}"
            + (f" in cell types {args.cell_types}" if args.cell_types else "")
            + "."
        )
        return

    logger.info(
        f"Found {len(attenuated_df)} attenuated pair(s) across "
        f"{attenuated_df['tissue'].nunique()} cell type(s) and "
        f"{attenuated_df['endo_feature'].nunique()} unique gene feature(s)."
    )

    # Rename conditioned columns for clarity
    rename_dict = {
        "exposure_coef": "cond_coef",
        "exposure_stderr": "cond_stderr",
        "exposure_tval": "cond_tval",
        "exposure_pval": "cond_pval",
        "exposure_fdr": "cond_fdr",
    }
    attenuated_df = attenuated_df.rename(columns=rename_dict)

    # Merge baseline endogenous metrics if available
    if endo_df is not None:
        endo_cols_to_merge = ["tissue", "feature"]
        for c in ["coef", "stderr", "p-value", "fdr_bh", "log2fc", "percentchange"]:
            if c in endo_df.columns:
                endo_cols_to_merge.append(c)

        merged = attenuated_df.merge(
            endo_df[endo_cols_to_merge],
            left_on=["tissue", "endo_feature"],
            right_on=["tissue", "feature"],
            how="left",
        )
        if "feature" in merged.columns:
            merged = merged.drop(columns=["feature"])

        merged = merged.rename(
            columns={
                "coef": "baseline_coef",
                "stderr": "baseline_stderr",
                "p-value": "baseline_pval",
                "fdr_bh": "baseline_fdr",
                "log2fc": "baseline_log2fc",
                "percentchange": "baseline_pctchange",
            }
        )

        # Calculate attenuation effect reduction metrics
        if "baseline_coef" in merged.columns and "cond_coef" in merged.columns:
            merged["delta_coef"] = merged["cond_coef"] - merged["baseline_coef"]
            # Percentage reduction in magnitude of effect size: (1 - |cond| / |base|) * 100
            valid_base = merged["baseline_coef"].replace(0, np.nan).abs()
            merged["pct_effect_reduction"] = (
                (1.0 - (merged["cond_coef"].abs() / valid_base)) * 100.0
            ).round(2)
    else:
        merged = attenuated_df

    # Merge cis-correlation metrics if available
    if cis_df is not None:
        cis_cols = ["tissue", "endo_feature", "exog_feature"]
        for c in ["coef", "p-value", "bh_fdr"]:
            if c in cis_df.columns:
                cis_cols.append(c)
        cis_subset = cis_df[cis_cols].drop_duplicates(subset=["tissue", "endo_feature", "exog_feature"])
        merged = merged.merge(
            cis_subset,
            on=["tissue", "endo_feature", "exog_feature"],
            how="left",
            suffixes=("", "_cis"),
        )
        merged = merged.rename(
            columns={
                "coef": "cis_coef",
                "p-value": "cis_pval",
                "bh_fdr": "cis_fdr",
            }
        )

    # Sort results
    sort_cols = [c for c in ["tissue", "cond_pval", "endo_feature"] if c in merged.columns]
    merged = merged.sort_values(by=sort_cols).reset_index(drop=True)

    # Organize display columns
    key_order = [
        "tissue",
        "endo_feature",
        "exog_feature",
        "baseline_coef",
        "baseline_pval",
        "baseline_fdr",
        "cond_coef",
        "cond_pval",
        "cond_fdr",
        "pct_effect_reduction",
        "cis_coef",
        "cis_pval",
    ]
    ordered_cols = [c for c in key_order if c in merged.columns]
    remaining_cols = [c for c in merged.columns if c not in ordered_cols]
    final_df = merged[ordered_cols + remaining_cols]

    # Print summary to stdout
    print("\n" + "=" * 90)
    print(f"ATTENUATED FEATURES SUMMARY (Target: {args.target_variable}, Alpha threshold: {args.alpha})")
    print("=" * 90)

    # Pretty-print format for float columns
    display_df = final_df[ordered_cols].copy()
    for col in display_df.select_dtypes(include=[np.number]).columns:
        if "pval" in col or "fdr" in col:
            display_df[col] = display_df[col].apply(lambda x: f"{x:.4e}" if pd.notnull(x) else "")
        elif "reduction" in col:
            display_df[col] = display_df[col].apply(lambda x: f"{x:.2f}%" if pd.notnull(x) else "")
        else:
            display_df[col] = display_df[col].apply(lambda x: f"{x:.4f}" if pd.notnull(x) else "")

    print(display_df.to_string(index=False))
    print("=" * 90 + "\n")

    # Output file handling
    out_path = None
    if args.output_file:
        out_path = Path(args.output_file)
    elif args.save:
        suffix_part = f".{args.output_suffix}" if args.output_suffix else ""
        out_path = (
            results_dir
            / f"{args.project}.{args.endo_modality}-{args.exog_modality}.{args.regression_type}.conditioned.{args.target_variable}{suffix_part}.attenuated.csv"
        )

    if out_path:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        final_df.to_csv(out_path, index=False)
        logger.info(f"Saved attenuated features to: {out_path}")


if __name__ == "__main__":
    main()
