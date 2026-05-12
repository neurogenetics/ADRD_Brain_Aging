import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pybedtools

logger = logging.getLogger(__name__)

DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run feature enrichment permutation testing for age-associated features."
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
    parser.add_argument(
        "--modality",
        type=str,
        required=True,
        choices=["rna", "atac"],
        help="Data modality (rna or atac).",
    )
    parser.add_argument(
        "--name",
        type=str,
        default="",
        help="Optional custom name to append to the output and log files (e.g. 'centromere').",
    )
    parser.add_argument(
        "--metrics-file",
        type=str,
        help="Path to a pre-computed CSV containing feature metrics.",
    )
    parser.add_argument(
        "--annotation-csv",
        type=str,
        help="Path to a UCSC Table Browser exported CSV containing genomic annotations.",
    )
    parser.add_argument(
        "--distance",
        action="store_true",
        help="Calculate enrichment based on average distance to nearest annotation rather than overlap.",
    )
    parser.add_argument(
        "--permutations",
        type=int,
        default=1000,
        help="Number of permutations to run (default: 1000).",
    )
    parser.add_argument(
        "--fdr-threshold",
        type=float,
        default=0.05,
        help="FDR threshold for defining significant features (default: 0.05).",
    )
    return parser.parse_args()


def empirical_p_value(observed, permuted):
    """Calculates two-tailed empirical p-value."""
    permuted = np.array(permuted)
    # Proportion of permuted values greater/less than the observed value
    p_greater = np.sum(permuted >= observed) / len(permuted)
    p_less = np.sum(permuted <= observed) / len(permuted)
    
    # Two-tailed p-value
    p_val = min(p_greater, p_less) * 2
    
    # Ensure p-value bounds
    return min(max(p_val, 1 / (len(permuted) + 1)), 1.0)


def run_metric_permutations(
    background_df: pd.DataFrame,
    sig_features: list,
    metric_cols: list,
    n_perms: int,
):
    """Permutation test for numeric metrics."""
    results = {}
    n_sig = len(sig_features)
    
    if n_sig == 0:
        return results

    # Get observed means
    obs_means = background_df.loc[sig_features, metric_cols].mean()

    # Pre-allocate array for permutations
    perm_means = {col: np.zeros(n_perms) for col in metric_cols}

    feature_indices = background_df.index.values

    for i in range(n_perms):
        # Sample without replacement
        sampled_features = np.random.choice(feature_indices, size=n_sig, replace=False)
        sampled_means = background_df.loc[sampled_features, metric_cols].mean()
        
        for col in metric_cols:
             perm_means[col][i] = sampled_means[col]

    for col in metric_cols:
        obs = obs_means[col]
        perms = perm_means[col]
        exp = np.mean(perms)
        fc = obs / exp if exp > 0 else np.nan
        pval = empirical_p_value(obs, perms)
        
        results[col] = {
            "observed": obs,
            "expected": exp,
            "fold_change": fc,
            "p_value": pval,
        }
        
    return results


def run_categorical_permutations(
    background_df: pd.DataFrame,
    sig_features: list,
    cat_cols: list,
    n_perms: int,
):
    """Permutation test for categorical metrics (proportions)."""
    results = {}
    n_sig = len(sig_features)
    
    if n_sig == 0:
        return results

    feature_indices = background_df.index.values

    for col in cat_cols:
        # Get all unique categories for this column in the background
        categories = background_df[col].dropna().unique()
        
        # Observed proportions in significant features
        obs_counts = background_df.loc[sig_features, col].value_counts()
        obs_props = {cat: obs_counts.get(cat, 0) / n_sig for cat in categories}
        
        # Pre-allocate arrays for permutations
        perm_props = {cat: np.zeros(n_perms) for cat in categories}
        
        for i in range(n_perms):
            sampled_features = np.random.choice(feature_indices, size=n_sig, replace=False)
            sampled_counts = background_df.loc[sampled_features, col].value_counts()
            for cat in categories:
                 perm_props[cat][i] = sampled_counts.get(cat, 0) / n_sig
                 
        for cat in categories:
             obs = obs_props[cat]
             perms = perm_props[cat]
             exp = np.mean(perms)
             fc = obs / exp if exp > 0 else np.nan
             pval = empirical_p_value(obs, perms)
             
             test_name = f"{col}::{cat}"
             results[test_name] = {
                 "observed": obs,
                 "expected": exp,
                 "fold_change": fc,
                 "p_value": pval,
             }

    return results


def prep_bedtool(df: pd.DataFrame, name_col=None) -> pybedtools.BedTool:
    """Converts a DataFrame with chrom, chromStart, chromEnd to a pybedtools.BedTool."""
    # Ensure standard ordering for BedTool
    cols = ["chrom", "chromStart", "chromEnd"]
    if name_col and name_col in df.columns:
         cols.append(name_col)
    
    # Filter, clean, and sort
    bed_df = df[cols].copy()
    
    # Drop any rows with missing coordinates
    bed_df.dropna(subset=["chrom", "chromStart", "chromEnd"], inplace=True)
    
    # Ensure coordinates are integers (required by bedtools)
    bed_df["chromStart"] = bed_df["chromStart"].astype(int)
    bed_df["chromEnd"] = bed_df["chromEnd"].astype(int)
    
    # Convert to strings and create a BedTool from the string representation
    bed_str = bed_df.to_csv(sep="\t", index=False, header=False)
    return pybedtools.BedTool(bed_str, from_string=True)


def run_overlap_permutations(
    background_df: pd.DataFrame,
    sig_features: list,
    annotation_bed: pybedtools.BedTool,
    n_perms: int,
):
    """Permutation test for genomic overlaps using pybedtools."""
    # Create background BedTool
    # We need the feature ID as the name to map intersections back
    # Assuming the index is the feature ID
    bg_bed_df = background_df.copy()
    bg_bed_df["feature_id"] = bg_bed_df.index
    
    # Map common coordinate columns if they exist
    col_mapping = {
        "chr": "chrom",
        "start": "chromStart",
        "end": "chromEnd"
    }
    bg_bed_df.rename(columns=col_mapping, inplace=True)
    
    # Check for coordinate columns
    coord_cols = ["chrom", "chromStart", "chromEnd"]
    if not all(col in bg_bed_df.columns for col in coord_cols):
        logger.warning(f"Background features missing coordinate columns {coord_cols}. Found: {bg_bed_df.columns.tolist()}")
        return None

    # Pre-calculate overlaps for the entire background to speed up permutations
    bg_bed = prep_bedtool(bg_bed_df, name_col="feature_id")
    
    # wa=True keeps the original A feature, u=True means return A once if it overlaps B
    intersect_res = bg_bed.intersect(annotation_bed, wa=True, u=True)
    
    # Extract the IDs of features that overlapped
    overlapped_feature_ids = set([interval.name for interval in intersect_res])
    logger.info(f"Total overlapping features in background: {len(overlapped_feature_ids)} / {len(background_df)}")

    # We now have a boolean vector of overlaps for the background
    background_df["is_overlapped"] = background_df.index.isin(overlapped_feature_ids)
    
    n_sig = len(sig_features)
    if n_sig == 0:
        return None

    # Observed overlap count
    obs_overlap = background_df.loc[sig_features, "is_overlapped"].sum()
    obs_rate = obs_overlap / n_sig

    perm_rates = np.zeros(n_perms)
    feature_indices = background_df.index.values

    for i in range(n_perms):
        sampled_features = np.random.choice(feature_indices, size=n_sig, replace=False)
        sampled_overlap = background_df.loc[sampled_features, "is_overlapped"].sum()
        perm_rates[i] = sampled_overlap / n_sig

    exp_rate = np.mean(perm_rates)
    fc = obs_rate / exp_rate if exp_rate > 0 else np.nan
    pval = empirical_p_value(obs_rate, perm_rates)

    return {
        "observed_rate": obs_rate,
        "expected_rate": exp_rate,
        "fold_change": fc,
        "p_value": pval,
        "observed_count": obs_overlap,
    }


def run_distance_permutations(
    background_df: pd.DataFrame,
    sig_features: list,
    annotation_bed: pybedtools.BedTool,
    n_perms: int,
):
    """Permutation test for average distance to nearest genomic annotation."""
    bg_bed_df = background_df.copy()
    bg_bed_df["feature_id"] = bg_bed_df.index
    
    col_mapping = {"chr": "chrom", "start": "chromStart", "end": "chromEnd"}
    bg_bed_df.rename(columns=col_mapping, inplace=True)
    
    coord_cols = ["chrom", "chromStart", "chromEnd"]
    if not all(col in bg_bed_df.columns for col in coord_cols):
        logger.warning(f"Background features missing coordinate columns {coord_cols}.")
        return None

    # Sort beds before finding closest to ensure proper execution
    bg_bed = prep_bedtool(bg_bed_df, name_col="feature_id").sort()
    ann_bed = annotation_bed.sort()
    
    # closest -d adds distance as the last column. t="first" handles multiple ties.
    closest_res = bg_bed.closest(ann_bed, d=True, t="first")
    
    feature_distances = {}
    for interval in closest_res:
        try:
            dist = float(interval[-1])
            # Ignore cases where closest wasn't found (-1)
            if dist >= 0:
                feature_distances[interval.name] = dist
        except (ValueError, IndexError):
            pass
            
    # Attach distances to background df
    background_df["annotation_distance"] = background_df.index.map(feature_distances)
    
    # Drop features with no valid distance (e.g., no annotations on chromosome)
    valid_bg = background_df.dropna(subset=["annotation_distance"]).copy()
    valid_sig = [f for f in sig_features if f in valid_bg.index]
    
    if len(valid_sig) == 0:
        return None
        
    # We can just reuse run_metric_permutations for this continuous metric!
    res = run_metric_permutations(valid_bg, valid_sig, ["annotation_distance"], n_perms)
    if "annotation_distance" in res:
        return res["annotation_distance"]
    return None


def main():
    args = parse_args()
    work_dir = Path(args.work_dir)
    results_dir = work_dir / "results"
    quants_dir = work_dir / "quants"
    logs_dir = work_dir / "logs"
    
    name_suffix = f".{args.name}" if args.name else ""
    
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_filename = logs_dir / f"{args.project}.{args.modality}{name_suffix}.feature_enrichment.log"
    
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_filename),
            logging.StreamHandler(sys.stdout)
        ],
        force=True,
    )
    
    logger.info(f"Command line: {' '.join(sys.argv)}")
    logger.info(f"Logging configured. Writing to {log_filename}")

    # 1. Load Data
    wls_file = results_dir / f"{args.project}.all_celltypes.{args.modality}.wls.age.csv"
    if not wls_file.exists():
        logger.error(f"WLS results file not found: {wls_file}")
        sys.exit(1)

    logger.info(f"Loading WLS results from {wls_file}")
    wls_df = pd.read_csv(wls_file)
    
    # 2. Load Feature Coordinates
    features_file = quants_dir / f"{args.project}.features.csv"
    if not features_file.exists():
        logger.error(f"Features file not found: {features_file}")
        sys.exit(1)
        
    logger.info(f"Loading feature coordinates from {features_file}")
    features_df = pd.read_csv(features_file, index_col=0) # Assuming first column is the feature ID

    # Merge WLS results with features to get coordinates
    # WLS results usually have 'feature' or 'gene' column. Let's assume 'feature'
    feature_col = "feature" if "feature" in wls_df.columns else wls_df.columns[0]
    logger.info(f"Assuming '{feature_col}' is the feature identifier in WLS results.")
    
    # 3. Load Metrics (Optional)
    metrics_cols_num = []
    metrics_cols_cat = []
    if args.metrics_file:
        logger.info(f"Loading metrics from {args.metrics_file}")
        metrics_df = pd.read_csv(args.metrics_file, index_col=0)
        
        # Identify numeric vs categorical columns for permutation testing
        metrics_cols_num = metrics_df.select_dtypes(include=np.number).columns.tolist()
        metrics_cols_cat = metrics_df.select_dtypes(exclude=np.number).columns.tolist()
        
        logger.info(f"Found {len(metrics_cols_num)} numeric metrics: {metrics_cols_num}")
        logger.info(f"Found {len(metrics_cols_cat)} categorical metrics: {metrics_cols_cat}")
        
        # Merge metrics into features_df
        features_df = features_df.join(metrics_df, how="left")

    # 4. Load Annotations (Optional)
    annotation_bed = None
    annotation_name = "overlap"
    if args.annotation_csv:
        logger.info(f"Loading annotations from {args.annotation_csv}")
        annot_df = pd.read_csv(args.annotation_csv)
        annotation_name = Path(args.annotation_csv).stem
        
        req_cols = ["chrom", "chromStart", "chromEnd"]
        if not all(col in annot_df.columns for col in req_cols):
             logger.error(f"Annotation CSV must contain columns: {req_cols}")
             sys.exit(1)
             
        annotation_bed = prep_bedtool(annot_df)
        logger.info(f"Loaded {len(annot_df)} annotation regions.")

    if not metrics_cols_num and not metrics_cols_cat and annotation_bed is None:
         logger.warning("Neither metrics-file nor annotation-csv provided. Nothing to test.")
         sys.exit(0)

    # 5. Permutation Testing per Cell Type
    all_results = []
    
    # Group WLS results by cell type ('tissue' column)
    if "tissue" not in wls_df.columns:
        logger.error("WLS results must contain a 'tissue' column for cell types.")
        sys.exit(1)
        
    for cell_type, ct_group in wls_df.groupby("tissue"):
        logger.info(f"Processing cell type: {cell_type}")
        
        # Background pool: all tested features in this cell type
        bg_features = ct_group[feature_col].unique()
        
        # Significant pool
        sig_features = ct_group[ct_group["fdr_bh"] < args.fdr_threshold][feature_col].unique()
        
        logger.info(f"  Background features: {len(bg_features)}")
        logger.info(f"  Significant features: {len(sig_features)}")
        
        if len(sig_features) == 0:
            logger.warning(f"  No significant features found for {cell_type}. Skipping.")
            continue
            
        # Subset features_df to only the background features
        # Ensure we only keep features that exist in features_df
        valid_bg = [f for f in bg_features if f in features_df.index]
        valid_sig = [f for f in sig_features if f in features_df.index]
        
        if len(valid_bg) == 0:
            logger.warning(f"  None of the tested features exist in the features file. Skipping.")
            continue
            
        ct_features_df = features_df.loc[valid_bg].copy()
        
        # Test Numeric Metrics
        if metrics_cols_num:
             logger.info("  Running numeric metric permutations...")
             metric_res = run_metric_permutations(
                 ct_features_df, valid_sig, metrics_cols_num, args.permutations
             )
             for metric, res in metric_res.items():
                 all_results.append({
                     "cell_type": cell_type,
                     "modality": args.modality,
                     "test_type": "numeric_metric",
                     "test_name": metric,
                     **res
                 })
                 
        # Test Categorical Metrics
        if metrics_cols_cat:
             logger.info("  Running categorical metric permutations...")
             cat_res = run_categorical_permutations(
                 ct_features_df, valid_sig, metrics_cols_cat, args.permutations
             )
             for metric, res in cat_res.items():
                 all_results.append({
                     "cell_type": cell_type,
                     "modality": args.modality,
                     "test_type": "categorical_metric",
                     "test_name": metric,
                     **res
                 })
                 
        # Test Overlaps or Distance
        if annotation_bed is not None:
             if args.distance:
                 logger.info(f"  Running distance permutations for {annotation_name}...")
                 dist_res = run_distance_permutations(
                     ct_features_df, valid_sig, annotation_bed, args.permutations
                 )
                 if dist_res:
                     all_results.append({
                         "cell_type": cell_type,
                         "modality": args.modality,
                         "test_type": "distance",
                         "test_name": annotation_name,
                         **dist_res
                     })
             else:
                 logger.info(f"  Running overlap permutations for {annotation_name}...")
                 overlap_res = run_overlap_permutations(
                     ct_features_df, valid_sig, annotation_bed, args.permutations
                 )
                 if overlap_res:
                     all_results.append({
                         "cell_type": cell_type,
                         "modality": args.modality,
                         "test_type": "overlap",
                         "test_name": annotation_name,
                         **overlap_res
                     })

    # 6. Save Results
    if all_results:
        out_df = pd.DataFrame(all_results)
        out_file = results_dir / f"{args.project}.{args.modality}{name_suffix}.feature_enrichment.csv"
        out_df.to_csv(out_file, index=False)
        logger.info(f"Saved enrichment results to {out_file}")
    else:
        logger.info("No results generated.")


if __name__ == "__main__":
    main()
