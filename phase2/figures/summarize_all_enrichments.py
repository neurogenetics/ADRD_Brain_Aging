import pandas as pd
import glob
import os
import numpy as np
import warnings

warnings.filterwarnings('ignore')

files = glob.glob("/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/results/aging_phase2.*.*.feature_enrichment.csv")

# Categorize files
categories = {
    "atac": {
        "Metrics & Specificity": ["peak_metrics", "cell_specific"],
        "Genomic Geography": ["centromere", "telomere"],
        "Evolutionary Elements": ["har", "haqer"],
        "Transposable Elements": ["TEs"],
        "Epigenetic Clocks": ["clocksites", "clocksites_Horvath2013", "clocksites_Shireby2020", "clocksites_Tong2024_BrainClock", "clocksites_Tong2024_Glia-In", "clocksites_Tong2024_Glia-Sin", "clocksites_Tong2024_Neu-In", "clocksites_Tong2024_Neu-Sin"],
        "ENCODE cCREs": ["Encode4_CA", "Encode4_CA-CTCF", "Encode4_CA-H3K4me34", "Encode4_CA-TF", "Encode4_Distal_enhancer", "Encode4_Promoter", "Encode4_Proximal_enhancer", "Encode4_TF"]
    },
    "rna": {
        "Metrics & Specificity": ["gene_metrics", "cell_specific"],
        "Genomic Geography": ["centromere", "telomere"],
        "Evolutionary Elements": ["har", "haqer"],
        "Aging Hallmarks": ["gene_aging_hallmarks", "gene_aging_hallmarks_level1", "gene_aging_hallmarks_level2", "gene_aging_hallmarks_level3", "gene_aging_hallmarks_level4", "gene_aging_hallmarks_level5"]
    }
}

def analyze_category(modality, cat_name, file_tags):
    dfs = []
    for tag in file_tags:
        file_path = f"/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/results/aging_phase2.{modality}.{tag}.feature_enrichment.csv"
        if os.path.exists(file_path):
            df = pd.read_csv(file_path)
            # Add file source tag to help trace back
            df['source'] = tag
            dfs.append(df)
            
    if not dfs:
        return
        
    combined = pd.concat(dfs, ignore_index=True)
    
    # Filter to significant
    sig = combined[combined['p_value'] < 0.05]
    if len(sig) == 0:
        return
        
    print(f"\n### {cat_name} ({modality.upper()})")
    print(f"**Total tests:** {len(combined)} | **Significant (p < 0.05):** {len(sig)}")
    
    # Analyze by test_name (which represents the specific enrichment tested)
    metric_summary = []
    
    for (t_name, t_type, source), group in sig.groupby(['test_name', 'test_type', 'source']):
        avg_fc = group['fold_change'].mean()
        min_p = group['p_value'].min()
        sig_cells = group['cell_type'].unique()
        
        # Sort cell types by most significant to print a few
        top_cells = group.sort_values('p_value').head(3)['cell_type'].tolist()
        cell_str = ", ".join(top_cells)
        if len(sig_cells) > 3:
            cell_str += f" (+{len(sig_cells)-3} others)"
            
        metric_summary.append({
            'source': source,
            'name': t_name,
            'type': t_type,
            'sig_count': len(sig_cells),
            'avg_fc': avg_fc,
            'min_p': min_p,
            'top_cells': cell_str
        })
        
    res_df = pd.DataFrame(metric_summary).sort_values(by=['sig_count', 'avg_fc'], ascending=[False, False])
    
    # Print Top 15 robust enrichments in this category
    for _, row in res_df.head(15).iterrows():
        fc_str = f"Over-represented (FC={row['avg_fc']:.2f})" if row['avg_fc'] >= 1 else f"Under-represented (FC={row['avg_fc']:.2f})"
        print(f"- **{row['name']}** *[{row['source']}]* -> Significant in **{row['sig_count']} cell types** (e.g., {row['top_cells']}). {fc_str}, Min p={row['min_p']:.4f}")

    # Print extreme outliers in Fold Change (abs log2 FC) to catch rare but intense effects
    print("\n*Notable extreme fold-changes (p < 0.01):*")
    extreme_sig = sig[sig['p_value'] < 0.01].copy()
    extreme_sig['abs_l2fc'] = np.abs(np.log2(extreme_sig['fold_change'].replace(0, np.nan)))
    extreme_sig = extreme_sig.sort_values('abs_l2fc', ascending=False).head(5)
    
    for _, row in extreme_sig.iterrows():
         if pd.notna(row['abs_l2fc']):
            dir_str = "Over-represented" if row['fold_change'] > 1 else "Under-represented"
            print(f"  - **{row['cell_type']}** | {row['test_name']} *[{row['source']}]*: FC={row['fold_change']:.2f} (p={row['p_value']:.4f})")
            
for mod in ["rna", "atac"]:
    print(f"\n# {'='*10} Modality: {mod.upper()} {'='*10}")
    for cat, tags in categories[mod].items():
        analyze_category(mod, cat, tags)

