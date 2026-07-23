import pandas as pd

atac_file = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/results/aging_phase2.atac.peak_metrics.feature_enrichment.csv"
rna_file = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2/results/aging_phase2.rna.gene_metrics.feature_enrichment.csv"

atac_df = pd.read_csv(atac_file)
rna_df = pd.read_csv(rna_file)

def summarize_df(df, name):
    sig = df[df['p_value'] < 0.05]
    print(f"\n# {name} Summary")
    print(f"Total tests: {len(df)}")
    print(f"Significant (p < 0.05): {len(sig)}")
    
    print("\n## Top Enriched Metrics (by number of significant cell types)")
    metric_counts = sig.groupby(['test_type', 'test_name']).size().sort_values(ascending=False).head(15)
    for (t_type, t_name), count in metric_counts.items():
        # avg FC for this metric
        avg_fc = sig[(sig['test_type']==t_type) & (sig['test_name']==t_name)]['fold_change'].mean()
        print(f"- {t_name} ({t_type}): Significant in {count} cell types. Mean FC (when sig): {avg_fc:.2f}")

    print("\n## Cell Types with Most Significant Enrichments")
    ct_counts = sig['cell_type'].value_counts().head(10)
    for ct, count in ct_counts.items():
        print(f"- {ct}: {count} significant enrichments")

    print("\n## Notable specific enrichments (p < 0.005, extreme FC)")
    # Filter highly significant
    high_sig = df[(df['p_value'] < 0.01)].copy()
    
    # Sort by absolute log2 fold change to get extreme over/under representation
    import numpy as np
    high_sig['abs_log2_fc'] = np.abs(np.log2(high_sig['fold_change'].replace(0, np.nan)))
    high_sig = high_sig.sort_values('abs_log2_fc', ascending=False).head(20)
    
    for _, row in high_sig.iterrows():
        print(f"- {row['cell_type']}: {row['test_name']} (FC={row['fold_change']:.2f}, p={row['p_value']:.4f})")

summarize_df(rna_df, "RNA (Gene Metrics)")
summarize_df(atac_df, "ATAC (Peak Metrics)")
