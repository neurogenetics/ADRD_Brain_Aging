import argparse
import logging
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from upsetplot import from_contents, plot as upset_plot

logger = logging.getLogger(__name__)

DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"

def parse_args():
    parser = argparse.ArgumentParser(description="Summarize and integrate QML pipeline results across cell types and modalities.")
    parser.add_argument("--project", type=str, default=DEFAULT_PROJECT)
    parser.add_argument("--work-dir", type=str, default=DEFAULT_WRK_DIR)
    parser.add_argument("--cnmf-dir-name", type=str, default="cnmf", help="Name of the cnmf output directory within the latents path.")
    parser.add_argument("--max-dist", type=int, default=1_000_000, help="Max distance to map ATAC peaks to RNA genes.")
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()

def jaccard_similarity(list1, list2):
    s1 = set(list1)
    s2 = set(list2)
    if not s1 and not s2:
        return 0.0
    return len(s1.intersection(s2)) / len(s1.union(s2))

def plot_jaccard_heatmap(features_dict, out_file, title):
    """Generates a pairwise Jaccard similarity heatmap for a dictionary of feature lists."""
    labels = list(features_dict.keys())
    n = len(labels)
    matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            matrix[i, j] = jaccard_similarity(features_dict[labels[i]], features_dict[labels[j]])
            
    df = pd.DataFrame(matrix, index=labels, columns=labels)
    
    plt.figure(figsize=(10, 8))
    sns.heatmap(df, annot=True, cmap="YlGnBu", vmin=0, vmax=1, fmt=".2f")
    plt.title(title)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    plt.close()
    logger.info(f"Saved Jaccard heatmap to {out_file}")

def generate_upset_plot(features_dict, out_file, title):
    """Generates an UpSet plot for a dictionary of feature lists."""
    if len(features_dict) < 2:
        logger.warning(f"Not enough groups to generate UpSet plot for {title}")
        return
        
    # upsetplot expects a dictionary where keys are group names and values are lists of items
    upset_data = from_contents(features_dict)
    
    plt.figure(figsize=(12, 8))
    upset_plot(upset_data, sort_by='cardinality')
    plt.title(title)
    plt.savefig(out_file, dpi=300)
    plt.close()
    logger.info(f"Saved UpSet plot to {out_file}")

def find_multiomic_pairs(rna_dict, atac_dict, features_df, max_dist, out_file):
    """
    Maps ATAC peaks to proximal RNA genes using the features_df coordinates.
    Finds pairs where BOTH the peak and the proximal gene are in the QML key features for the same cell type.
    """
    results = []
    
    # Common cell types between RNA and ATAC
    common_cts = set(rna_dict.keys()).intersection(set(atac_dict.keys()))
    
    if not common_cts:
        logger.warning("No common cell types found between RNA and ATAC results. Skipping multi-omic pair generation.")
        return
        
    for ct in common_cts:
        rna_features = set(rna_dict[ct])
        atac_features = set(atac_dict[ct])
        
        # Filter features_df to just the features relevant to this cell type
        ct_rna_df = features_df.loc[features_df.index.intersection(rna_features)].copy()
        ct_atac_df = features_df.loc[features_df.index.intersection(atac_features)].copy()
        
        for chrom in ct_rna_df['chr'].unique():
            chrom_rna = ct_rna_df[ct_rna_df['chr'] == chrom]
            chrom_atac = ct_atac_df[ct_atac_df['chr'] == chrom]
            
            if chrom_atac.empty:
                continue
                
            for rna_id, rna_row in chrom_rna.iterrows():
                start_boundary = rna_row["start"] - max_dist
                end_boundary = rna_row["end"] + max_dist
                
                # Find ATAC peaks within the window
                proximal_atac = chrom_atac[
                    (chrom_atac["start"] >= start_boundary) & 
                    (chrom_atac["end"] <= end_boundary)
                ]
                
                for atac_id, _ in proximal_atac.iterrows():
                    results.append({
                        "cell_type": ct,
                        "rna_gene": rna_id,
                        "atac_peak": atac_id,
                        "chrom": chrom,
                        "distance": min(abs(rna_row["start"] - proximal_atac.loc[atac_id, "end"]), 
                                        abs(rna_row["end"] - proximal_atac.loc[atac_id, "start"]))
                    })
                    
    if results:
        df = pd.DataFrame(results)
        df = df.sort_values(by=["cell_type", "distance"])
        df.to_csv(out_file, index=False)
        logger.info(f"Found {len(df)} multi-omic QML pairs. Saved to {out_file}")
    else:
        logger.info("No multi-omic QML pairs found within the specified distance.")

def main():
    args = parse_args()
    
    work_dir = Path(args.work_dir)
    results_dir = work_dir / "results"
    quants_dir = work_dir / "quants"
    figs_dir = work_dir / "figures"
    logs_dir = work_dir / "logs"
    
    figs_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    
    log_filename = logs_dir / f"{args.project}_qml_summary.log"
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler()],
        force=True,
    )
    
    logger.info(f"Command line: {' '.join(sys.argv)}")
    
    qml_dir = results_dir / "latents" / args.cnmf_dir_name / "qml_results"
    
    if not qml_dir.exists():
        logger.error(f"QML results directory not found: {qml_dir}. Please run run_qml_pipeline.py first.")
        sys.exit(1)
        
    # Load all key feature files
    # Dictionary structure: dict[modality][walk_type][cell_type] = list of features
    extracted_features = {
        "rna": {"quantum": {}, "classical": {}},
        "atac": {"quantum": {}, "classical": {}}
    }
    
    pattern = f"{args.project}_*_key_features.csv"
    file_paths = list(qml_dir.glob(pattern))
    
    if not file_paths:
        logger.error(f"No key feature files found matching pattern {pattern} in {qml_dir}")
        sys.exit(1)
        
    for fp in file_paths:
        try:
            filename = fp.name
            # Expected format: project_cell_type_modality_walktype_key_features.csv
            parts = filename.replace("_key_features.csv", "").split("_")
            # This parsing is tricky because cell_type can have underscores.
            # We know project is parts[0:2] usually (aging_phase2)
            # Walk type is the last part
            walk_type = parts[-1]
            # Modality is the second to last part
            modality = parts[-2]
            
            # The rest is the cell type (excluding project prefix)
            project_prefix_len = len(args.project.split("_"))
            cell_type = "_".join(parts[project_prefix_len:-2])
            
            df = pd.read_csv(fp)
            features = df["feature"].tolist()
            
            if modality in extracted_features and walk_type in extracted_features[modality]:
                extracted_features[modality][walk_type][cell_type] = features
                logger.debug(f"Loaded {len(features)} {walk_type} features for {cell_type} ({modality})")
                
        except Exception as e:
            logger.error(f"Error parsing file {fp}: {e}")
            
    # --- 1. Cross-Cell-Type Jaccard Heatmaps ---
    for modality in ["rna", "atac"]:
        for walk_type in ["quantum", "classical"]:
            data_dict = extracted_features[modality][walk_type]
            if len(data_dict) > 1:
                out_name = figs_dir / f"{args.project}_qml_{modality}_{walk_type}_jaccard.png"
                title = f"{walk_type.capitalize()} QML Features Jaccard Similarity ({modality.upper()})"
                plot_jaccard_heatmap(data_dict, out_name, title)
                
    # --- 2. UpSet Plots for Conserved Core Features ---
    for modality in ["rna", "atac"]:
        for walk_type in ["quantum", "classical"]:
            data_dict = extracted_features[modality][walk_type]
            if len(data_dict) > 1:
                out_name = figs_dir / f"{args.project}_qml_{modality}_{walk_type}_upset.png"
                title = f"{walk_type.capitalize()} QML Conserved Features ({modality.upper()})"
                generate_upset_plot(data_dict, out_name, title)
                
    # --- 3. Multi-Omic (RNA-ATAC) Cis-Regulatory Bridging ---
    logger.info("Generating multi-omic cross-modality pairs...")
    features_file = quants_dir / f"{args.project}.features.csv"
    
    if features_file.exists():
        features_df = pd.read_csv(features_file, index_col=0)
        
        for walk_type in ["quantum", "classical"]:
            rna_dict = extracted_features["rna"][walk_type]
            atac_dict = extracted_features["atac"][walk_type]
            
            if rna_dict and atac_dict:
                out_file = qml_dir / f"{args.project}_multiomic_pairs_{walk_type}.csv"
                find_multiomic_pairs(rna_dict, atac_dict, features_df, args.max_dist, out_file)
    else:
        logger.warning(f"features.csv not found at {features_file}. Skipping multi-omic integration.")
        
    logger.info("QML Summarization complete.")

if __name__ == "__main__":
    main()
