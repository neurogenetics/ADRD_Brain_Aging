# Latent factor analysis of aging results
1. Analysis
   - Latents
       - Per broad and specific cell-types generate latent factors of the age associated gene expression and ATAC peak peak features using [ICA](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.FastICA.html), [NMF](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.NMF.html), and [PCA](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html). Test each of these latent factors for association with age. latent_factors_age_analysis.ipynb and pm_run_latent_age_analysis.ipynb
       - Combine, aggregate, and summarize the results across cell-types. post_latent_factor_analyses.ipynb
3. Figures
    - Age associated features
        - Create a graph of age associated features and partition the graph with Leiden clustering based on [Leiden ModularityVertexPartition](https://leidenalg.readthedocs.io/en/stable/reference.html#modularityvertexpartition). association_graph.ipynb
        - Run GSEA analysis per age associated feature graph partitions. gsea_partition_features.ipynb
        - Generate of Sankey diagram of the cell-type age associated features and their graph paritions. partitioned_features_sankey.ipynb
    - Latents
        - Run GSEA analysis per age associated latent factor for all broad and specific cell-types. gsea_celltype_latents_age_factors.ipynb
        - Generate Sankey diagram per broad cell-types and their specific cell-types latent factors and the GSEA pathways detected. celltype_latents_sankey.ipynb and pm_run_celltype_latents_sankey.ipynb
        - Create a graph of age associated latent factors and their feature loadings and partition the graph with Leiden clustering based on [Leiden ModularityVertexPartition](https://leidenalg.readthedocs.io/en/stable/reference.html#modularityvertexpartition). latents_graph.ipynb
        - Find the high degree features in the graph of age associated latent factors. latent_graphs_high_degree_features.ipynb
        - Run GSEA Enrichr analysis for each of the partitions in the graph of age associated latent factors. gsea_partitioned_latents_age_factors.ipynb
        - Generate of Sankey diagram of the age associated features and paritions of the age associated latent factors along with there GSEA pathways. partitioned_latents_sankey.ipynb