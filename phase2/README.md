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
    - take a look at the ambient RNA results, cellbender_results.ipynb
5. Create combined demultiplexed AnnData files for each modality    
    - ATAC use aggregated pool data; create_aggr_anndata_with_demuxlet_identified_donors.ipynb
    - GEX
        - migrate the phase1 GEX pools for phase1 pools 4 and 5; migrate_phase1_gex_pools.ipynb
        - for each GEX pool create demultiplexed AnnData files; create_anndata_with_demuxlet_identified_donors.ipynb and pm_run_create_demultiplexed_gex_anndatas.ipynb 
        - concatenated the GEX AnnData objects into single AnnData files using Scanpy; combine_demultiplexed_pool_anndatas.ipynb
7. Convert ARC data to AnnData object and populate the 'obs' info; convert_arc_data.ipynb
8. Check for additional non-genotype doublet cells using Scrublet; scrublet.ipynb
9. Use Celltypist with human brain models to predict cell-type labels as a starting point; prep_celltypist_input.ipynb
10. Perform initial clustering and annotation of cell types using only the mRNA (GEX and ARC cells) features; initial_clustering_rna.ipynb
11. Subcluster by broad cell type classes, exictatory, inhibitory, and non-neuronal to refine and manually curate appropriate cell type labels
    - For each broad cell type class perform clustering for a range of resolutions; broadtype_subcluster_rna.ipynb
    - Use SAHA for cluster and resolution refinement and to curate final cell-type labels; ?
12. Perform joint analysis of RNA and ATAC using MultiVI so final RNA annotation labels can be transferred to ATAC cells; multivi_joint_analysis.ipynb
13. Transfer RNA cell-type annotation labels to the ATAC cells; predict_celltypes_from_multivi.ipynb
## Analysis
1. Identify gene expression and chromatin accessibility features associated with age per broad cell-type and cluster specfic cell-type
    - Convert the single-cell data to pseudobulk (mean) values for each broad and cluster specific cell-type for both GEX and ATAC data; pseudobulk_convert.py
    - Format covariate tables for use with data prep and regression analysis; format_covariates.py
    - Prepare the pseudobulk data for analysis by analysing and removing non-target variable variance; prep_pb_data.py
    - Regression analysis between quantified features (expression and accessibility) and age; pseudobulk_regression_analysis.ipynb. Where the possible regression methods include GLM, GLM with Tweedie distribution, and RLM.
    - Post-proceesing of the regression analysis across cell-types to apply B&H FDR and identify the statistally significant linear correlations between feature quantification and age; post_pseudobulk_regression.ipynb
    - Filter outlier effects from age regression from the GLM Tweedie results based on the RLM results; filter_regression_type_differences.ipynb
2. Identify chromatin accessibility features correlated with cis proximal age associated gene expression features
    - Regression analysis between quantified age associated gene expression features and cis proximal chromatin accessibility features; cis_correlation.ipynb. Where the possible regression methods include GLM, GLM with Tweedie distribution, and RLM.
    - Post-proceesing of the regression analysis across cell-types to apply B&H FDR and identify the statistally significant linear correlations between age associated gene expression features and cis proximal chromatin accessibility features; post_cis_correlation.ipynb
    - Filter outlier effects from age regression from the GLM Tweedie results based on the RLM results; filter_regression_type_differences.ipynb
3. Conditioned age regression analysis. Rerun age regression analysis for the age associated GEX features conditioned on <i>cis</i> proximal correlated ATAC features that are also age associated. cis_conditioned_regression_analysis.ipynb.
    - Summarize the conditioned age regression differences between cell-types. post_cis_conditioned_regression.ipynb
## Visualization
1. View a cell-type specific feature ~ age result show the model summary and some scatter plots; specific_age_result.ipynb
2. View a cell-type specific cis proximal GEX ~ ATAC result show the model summary and some scatter plots; specific_cis_result.ipynb
3. View to proportion of features tested in a cell-type that are age associated; proportion_age_effects.ipynb
