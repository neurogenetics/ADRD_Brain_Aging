# ADRD_Brain_Aging Phase 2
Single-nuclei atlas of aging in the entorhinal cortex using RNA and ATAC 10X multiome assay

## Processing steps
1. NISC ran Cellranger (gex v7.1.0, -arc v2.0.0, -atac v2.0.0) thru count
2. Cellranger (-arc v2.0.1, -atac v2.1.0) aggr was run to aggregate the count outputs across the samples; cellranger-arc_aggr.ipynb, cellranger-atac_aggr.ipynb, create_consensus_peaks_bed.ipynb
3. Demultiplex pooled cell data based on genotypes using demuxlet; prep_genotypes.ipynb


3. scvi-tools MultiVi was used to generate latent variables and format the data into AnnData objects; MultiVI_analysis.ipynb
4. Detect doublet cells using Scrublet; scrublet.ipynb
5. Preliminary cell-type assignment was done using MACA and CNS markers from PangloaDB and Bakken motor cortex. These are not fully appropriate for the brain regions of this experiment but work as a good starting point for cell assignments prior to further curation. label_cells_maca.ipynb via Papermill runner generating and running a notebook per marker set, label_cells_runner.ipynb
6. QC, clustering, and visualization was done with SCANPY, the clustering is performed with multiple Leiden resolutions; scanpy_runner.ipynb and scanpy_processing.ipynb
7. Clustering and currated cell assignment done with SCANPY; scanpy_tuning.ipynb