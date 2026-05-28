import argparse
import logging
from pathlib import Path
import sys

import pandas as pd
import numpy as np
import torch
import networkx as nx

# Add src to path so we can import qml_pipeline
sys.path.append(str(Path(__file__).resolve().parent.parent / "src"))

from qml_pipeline import compute_cosine_similarity

logger = logging.getLogger(__name__)

DEFAULT_PROJECT = "aging_phase2"
DEFAULT_WRK_DIR = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"

def parse_args():
    parser = argparse.ArgumentParser(description="Export divergent QML subnetworks as GEXF for Gephi/Cytoscape.")
    parser.add_argument("--project", type=str, default=DEFAULT_PROJECT)
    parser.add_argument("--work-dir", type=str, default=DEFAULT_WRK_DIR)
    parser.add_argument("--cell-type", type=str, required=True, help="Target cell type.")
    parser.add_argument("--modality", type=str, default="rna", choices=["rna", "atac"])
    parser.add_argument("--cnmf-dir-name", type=str, default="cnmf")
    parser.add_argument("--threshold", type=float, default=0.85, help="Minimum cosine similarity to retain an edge.")
    parser.add_argument("--debug", action="store_true")
    return parser.parse_args()

def main():
    args = parse_args()
    
    work_dir = Path(args.work_dir)
    results_dir = work_dir / "results"
    figs_dir = work_dir / "figures"
    logs_dir = work_dir / "logs"
    
    figs_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    safe_ct = args.cell_type.replace(" ", "_").replace("/", "-")
    
    log_filename = logs_dir / f"{args.project}_{args.modality}_{safe_ct}_export_topology.log"
    logging.basicConfig(
        level=logging.DEBUG if args.debug else logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.FileHandler(log_filename), logging.StreamHandler()],
        force=True,
    )
    
    logger.info(f"Command line: {' '.join(sys.argv)}")
    
    cnmf_dir = results_dir / "latents" / args.cnmf_dir_name
    qml_dir = cnmf_dir / "qml_results"
    
    # 1. Load Key Features
    q_file = qml_dir / f"{args.project}_{safe_ct}_{args.modality}_quantum_key_features.csv"
    c_file = qml_dir / f"{args.project}_{safe_ct}_{args.modality}_classical_key_features.csv"
    
    if not q_file.exists() or not c_file.exists():
        logger.error("Key feature files not found. Run run_qml_pipeline.py first.")
        sys.exit(1)
        
    q_df = pd.read_csv(q_file)
    c_df = pd.read_csv(c_file)
    
    q_features = set(q_df["feature"].tolist())
    c_features = set(c_df["feature"].tolist())
    
    target_nodes = list(q_features.union(c_features))
    logger.info(f"Found {len(target_nodes)} unique nodes across both walks (Quantum: {len(q_features)}, Classical: {len(c_features)})")
    
    if len(target_nodes) == 0:
        logger.warning("No features to export.")
        return
        
    # Create dicts for easy probability lookup
    q_probs = dict(zip(q_df["feature"], q_df["probability"]))
    c_probs = dict(zip(c_df["feature"], c_df["probability"]))
    
    # 2. Determine K and Load H Matrix
    regression_type = "pb_wls"
    reg_file = cnmf_dir / "wls_pb_results" / f"{args.project}_{safe_ct}_{args.modality}_{regression_type}.csv"
    
    if not reg_file.exists():
        logger.error(f"Regression results file not found: {reg_file}")
        sys.exit(1)
        
    reg_df = pd.read_csv(reg_file)
    selected_k = int(reg_df["k"].iloc[0])
    
    run_name = f"{args.project}_{safe_ct}_{args.modality}"
    run_dir = cnmf_dir / run_name
    score_files = list(run_dir.glob(f"*spectra_score.k_{selected_k}.dt_*.txt"))
    
    if not score_files:
        logger.error(f"No spectra score file found for K={selected_k}")
        sys.exit(1)
        
    logger.info(f"Loading cNMF spectra scores from: {score_files[0]}")
    h_df = pd.read_csv(score_files[0], sep='\t', index_col=0)
    all_feature_names = h_df.columns.tolist()
    
    # 3. Compute Similarity Matrix
    logger.info("Computing Cosine Similarity Matrix...")
    h_tensor = torch.tensor(h_df.values, dtype=torch.float32)
    # Using the same function as the pipeline to ensure exact parity
    similarity_matrix = compute_cosine_similarity(h_tensor)
    
    # 4. Extract Induced Subgraph
    logger.info("Extracting Induced Subgraph for target nodes...")
    sim_df = pd.DataFrame(similarity_matrix.numpy(), index=all_feature_names, columns=all_feature_names)
    
    # Ensure targets are actually in the dataframe
    valid_targets = [n for n in target_nodes if n in sim_df.columns]
    subgraph_df = sim_df.loc[valid_targets, valid_targets]
    
    # 5. Thresholding
    logger.info(f"Applying threshold >= {args.threshold}...")
    # Set diagonals and low values to 0 to remove them from networkx edgelist
    np.fill_diagonal(subgraph_df.values, 0)
    subgraph_df[subgraph_df < args.threshold] = 0
    
    # 6. Build NetworkX Graph
    logger.info("Building NetworkX graph...")
    G = nx.from_pandas_adjacency(subgraph_df)
    logger.info(f"Graph constructed with {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
    
    # Remove isolated nodes (nodes with 0 edges after thresholding)
    isolated = list(nx.isolates(G))
    G.remove_nodes_from(isolated)
    logger.info(f"Removed {len(isolated)} isolated nodes. Final size: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges.")
    
    # 7. Add Node Attributes
    for node in G.nodes():
        if node in q_features and node in c_features:
            category = "Shared"
        elif node in q_features:
            category = "Quantum Specific"
        else:
            category = "Classical Specific"
            
        G.nodes[node]["category"] = category
        G.nodes[node]["quantum_prob"] = float(q_probs.get(node, 0.0))
        G.nodes[node]["classical_prob"] = float(c_probs.get(node, 0.0))
        
    # 8. Export
    out_file = figs_dir / f"{args.project}_{safe_ct}_{args.modality}_topology.gexf"
    nx.write_gexf(G, out_file, version="1.2draft") # The warning in Gephi is usually about 1.2 vs 1.3, but networkx supports 1.2draft
    # Wait, networkx write_gexf supports version="1.2draft" (default) or "1.1draft".
    # Let me check the documentation. Actually, modern networkx supports '1.2draft'
    # Gephi prefers .graphml for modern formats without deprecation warnings.
    # Let's switch to graphml which is widely supported and never throws this warning in Gephi.
    out_file_gml = figs_dir / f"{args.project}_{safe_ct}_{args.modality}_topology.graphml"
    nx.write_graphml(G, out_file_gml)
    logger.info(f"Successfully exported topology to {out_file_gml}")

if __name__ == "__main__":
    main()
