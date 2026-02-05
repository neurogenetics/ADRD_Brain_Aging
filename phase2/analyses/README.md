# Analysis
1. Prepare pseudobulk converted data for age regression analysis, done per cell-type and modality: prep_pb_data.py
  Core Workflow:
    1. Data Loading: Loads quantified data (RNA/ATAC) and sample covariates (age, sex, ancestry, technical metrics) for a specific cell type.
    2. Initial Variance Partition: Models gene expression variance against known covariates (e.g., age, sex, BMI, pool) to understand initial drivers of variation.
    3. Latent Variable Discovery:
        * Identifies high-variance autosomal features.
        * Imputes missing values (KNN, Iterative, or Simple).
        * Regresses out "age" effects to focus on non-age-related variance.
        * Performs PCA on residuals and automatically selects the optimal number of components (using knee-point detection on R² and RMSE curves) to capture unobserved technical/biological noise.
    4. Final Variance Partition: Re-evaluates variance partitioning including the newly generated PCA components to verify they capture the residual variance.
    5. Covariate Correction: Regresses out all non-target covariates (known technical factors + latent PCA components), leaving only age-related variance in the final residuals.
    6. Output: Saves the "clean" residuals and updated covariate files for downstream regression analysis.
  Key Technical Details:
    * Parallelization: Uses concurrent.futures.ProcessPoolExecutor for computationally intensive variance partitioning.
    * Modeling: Relies on statsmodels for GLM/correlation checks and sklearn for PCA, imputation, and linear regression.
    * Visualization: Generates heatmaps (covariate vs. PCA correlations) and boxen plots (variance explained distributions) to quality check the process.

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

