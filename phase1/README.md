# ADRD_Brain_Aging Phase 1 
- single-nuclei atlas between young and old subjects in four brain regions
- Where 'young' subjects were between 20-30 years of age at time of death and 'old' were 60-80 years of age at time of death
- brain regions included for each subject are: putamen, entorhinal cortex, middle temporal gyrus, and subventricular zone

## Analyses
### Demultiplexing
    - use subject genotypes to deconvolute cells from pools and assign each cell to a subject using Demuxlet to eastablish identity based on genotype
    1. prep_hbcc_genos.ipynb
    2. prep_pool_sample_info.ipynb
    3. run_demuxlet_wdl_job.ipynb
    4. create_anndata_with_demuxlet_identified_donors.ipynb
    5. combine_demultiplexed_pool_anndatas.ipynb
### Clustering and cell-type identification
    1. pegasus_analysis.ipynb
    2. subset_recluster.ipynb
    N. other clustering and processing testing not used for analysis; scanpy_cluster_analysis.ipynb and scvi_id_zi_features.ipynb
### Differential expression analysis by age group
    1. frmt_glmmtmb_diffexp.ipynb
    2. manually_run_glmmTMB_GCP.md and glmmTMB.R
    3. post_glmmtmb_diffexp.ipynb
    N. other analyses types run for testing but not used for analysis; glmm_diffexp.ipynb, glmmtmb_diffexp.ipynb, and glm_pb_diffexp.ipynb with corrsponding post processing notebooks and comparisons. post_glm_pb_diffexp.ipynb, post_glmm_diffexp.ipynb, post_glmm_zi_diffexp.ipynb, and compare_de_methods.ipynb