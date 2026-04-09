1. Summarization and Visualization of Correlated Pairs (cis_correlation_summary.py)
   * Objective: To quantify the proportion of age-associated features that are also cis-correlated and to visually summarize these proportions across
     distinct cell types.
   * Data Integration and Filtering:
       * File Aggregation: The script ingests the cell-type aggregated age-association results for both modalities (e.g., RNA_fdr_filtered.age.csv,
         ATAC_fdr_filtered.age.csv) and the raw, unfiltered cis-correlation results.
       * Age-Association Restraint: The raw cis-correlation results are strictly filtered to include only those pairs where both the endogenous and
         exogenous features are established age-associated features.
   * Re-evaluation of Significance:
       * Dynamic FDR Adjustment: Because the correlation space is restricted to only age-associated features, the BH FDR is dynamically re-calculated
         solely on this reduced subset of p-values.
       * Significance Thresholding: Correlated pairs are deemed significant based on this updated, context-specific FDR threshold (e.g., FDR < 0.05).
   * Quantification of Correlation Rates:
       * For each cell type, the pipeline calculates the absolute number of significant features in both modalities.
       * It determines the number of significant cis-correlated pairs.
       * It calculates the percentage of age-associated features in each modality that participate in at least one significant cis-correlation
         ([Correlated Features] / [Total Age-Associated Features] * 100).
   * Reporting and Plotting:
       * Summary Logging: A detailed statistical table containing these counts and percentages is generated and logged for descriptive reporting.
       * Visualization: A grouped bar plot is generated (using Seaborn) illustrating the percentage of cis-correlated features per cell type,
         stratified by modality. The legend is positioned externally to prevent occlusion of the data.

2. Cell-Type Age Effect Similarity Visualization (celltype_age_effect_similarity.py)
   * Objective: To visually assess and quantify the similarities between different cell types based on their age-associated feature effect sizes within a specific modality.
   * Feature Selection:
       * The script identifies a universal set of significant age-associated features by reading the FDR-filtered results (`_fdr_filtered.age.csv`). A feature is included if it is significant in at least one cell type.
   * Data Aggregation and Matrix Construction:
       * Full regression results (`age.csv`) for the selected features are loaded across all cell types.
       * A feature-by-cell-type matrix is constructed using a user-specified effect size metric (e.g., `coef`, `z`, `log2fc`).
   * Similarity Computation and Visualization:
       * A Spearman rank correlation matrix is computed to determine the pairwise similarity between cell types.
       * A clustered heatmap (using Seaborn's `clustermap`) is generated to visualize these correlations. The dendrograms are hidden for a cleaner aesthetic, and the color bar is positioned externally.
   * Output:
       * The resulting heatmap is saved in both PNG and SVG formats within the figures directory.

3. Summarization and Visualization of Cis Conditioned Regression (cis_conditioned_regression_summary.py)
   * Objective: To summarize and visualize the results of the cis conditioned regression analysis, quantifying the extent to which age-associated genes are mediated by cis-correlated age-associated ATAC peaks.
   * Data Integration and Mediation Analysis:
       * File Ingestion: Loads the FDR-filtered age-associated features for the endogenous modality (e.g., RNA) and the regression results of the conditioned pairs.
       * Mediation Determination: For each endogenous feature, the script checks if conditioning on a cis-correlated ATAC peak causes a loss of significance for the age effect (conditioned exposure p-value > alpha threshold, default 0.05). This loss of significance indicates that the peak mediates the age effect.
   * Quantification of Mediation Rates:
       * Calculates the number and percentage of age-associated genes that are mediated by at least one cis-correlated ATAC peak per cell type.
       * Computes the mean distances between multiple mediating peaks for genes that are mediated by more than one peak.
   * Reporting and Visualization:
       * Summary Data: Outputs a tabular CSV summary of the mediation counts and proportions per cell type.
       * Bar Plot: Generates a bar plot showing the percentage of age-associated genes mediated by a cis-correlated ATAC peak, ordered by cell type.
       * Distance Boxen Plot: Visualizes the distribution of mean distances (in Kb) between mediating cis proximal peaks across cell types.
       * Peak Count Boxen Plot: Displays the distribution of the number of mediating cis proximal ATAC peaks per gene across cell types.
