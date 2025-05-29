# ADRD_Brain_Aging Phase 1
<img src="figures/P1_Figure1_Draft_MEM.png" width="320px" />

- Single-nuclei RNA atlas between young and aged subjects in four brain regions
- Where 'young' subjects were between 20-30 years of age at time of death and 'aged' were 60-85 years of age at time of death
- Brain regions included for each subject are: Enorhinal cortex (EC), Middle Temporal Gyrus (MTG), Subventricular Zone (SVZ), and Putamen (PUT)

## Demultiplexing
    - use subject genotypes to deconvolute cells from pools and assign each cell to a subject using Demuxlet to eastablish identity based on genotype
    1. quantifications/prep_hbcc_genos.ipynb
    2. quantifications/prep_pool_sample_info.ipynb
    3. quantifications/run_demuxlet_wdl_job.ipynb
    4. quantifications/create_anndata_with_demuxlet_identified_donors.ipynb via pm_nb_runners/create_demuxed_anndata_runner.ipynb  
    5. quantifications/combine_demultiplexed_pool_anndatas.ipynb
    
## Clustering and cell-type identification
    1. quantifications/pegasus_analysis.ipynb
    2. quantifications/subset_recluster.ipynb
    3. quantifications/manually_add_broad_celltypes.ipynb
    
## Differential expression analysis by age group
    1. analyses/frmt_glmmtmb_diffexp.ipynb via pm_nb_runners/glmmtmb_diffexp_runner.ipynb
        a. analyses/frmt_broad_types_by_region_glmmtmb_prep.ipynb
    2. analyses/manually_run_glmmTMB_GCP.md and glmmTMB.R
    3. analyses/post_glmmtmb_diffexp.ipynb
    
## Check ambient RNA and non-genotype doublets relative to doublets and ambigious reads from demux and 'uncertain' cell-type clusters
    1. quantifications/run_cellbender_wdl_job.ipynb
    2. quantifications/cellbender_results.ipynb
    3. quantifications/scrublet.ipynb
