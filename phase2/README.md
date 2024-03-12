# ADRD_Brain_Aging Phase 2
Single-nuclei atlas of aging in the entorhinal cortex using RNA and ATAC 10X multiome assay

## Quantification, QC, and Clustering steps
1. NISC ran Cellranger (gex v7.1.0, -arc v2.0.0, -atac v2.0.0) thru count
2. Cellranger (-arc v2.0.1, -atac v2.1.0) aggr was run to aggregate the count outputs across the samples; cellranger-arc_aggr.ipynb, cellranger-atac_aggr.ipynb, create_consensus_peaks_bed.ipynb
3. Demultiplex pooled cell data based on genotypes using demuxlet; 
    - prepare pool genotypes; prep_genotypes.ipynb
    - run Demuxlet WDL vis GCP Life Sciences; run_demuxlet_wdl_job.ipynb
    - create genotype demultiplexed AnnData files
4. Detect ambient RNA in GEX pools and non-pooled ARC samples, using CellBender
    - run Cellbender per sample, run_cellbender_wdl_job.ipynb
    - take a look at the ambient RNA results, Cellcellbender_results.ipynb
5. Create combined demultiplexed AnnData files for each modality    
    - ATAC use aggregated pool data; create_aggr_anndata_with_demuxlet_identified_donors.ipynb
    - GEX
        - for each GEX pool create demultiplexed AnnData files; create_anndata_with_demuxlet_identified_donors.ipynb and pm_run_create_demultiplexed_gex_anndatas.ipynb
        - migrate the phase1 GEX pools for phase1 pools 4 and 5; migrate_phase1_gex_pools.ipynb
        - concatenated the GEX AnnData objects into single AnnData files using Scanpy; combine_demultiplexed_pool_anndatas.ipynb
7. Convert ARC data to AnnData object and populate the 'obs' info; convert_arc_data.ipynb
8. Check for additional non-genotype doublet cells using Scrublet; scrublet.ipynb
9. scvi-tools MultiVi was used to generate latent variables across ARC, GEX, and ATAC and cluster the data; MultiVI_analysis.ipynb
10. Preliminary automated cell-type labeling using Phase1 labels and CellAssign predictions based on scTypes and Bakken et al marker sets; pm_run_cellassign.ipynb, scvi_cellassign.ipynb
     - compare cell-types assigments between label sets and with Leiden clusters; compare_celltype_predictions.ipynb
11. Populate the multiVI clustering and CellAssign to an anndata object that still retains the full features data instead of just the high variance features used for clustering; populate_full_anndata.ipynb
12. Annotate the curated cell-type assignments back into the MultiVI anndata object; annotate_curated_cluster.ipynb
