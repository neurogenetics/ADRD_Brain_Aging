# Analysis
1. Prepare pseudobulk converted data for age regression analysis, done per cell-type and modality: prep_pb_data.py
  Core Workflow:
    1. Data Loading: Loads quantified data (RNA/ATAC) and sample covariates (age, sex, ancestry, technical metrics) for a specific cell type.
    2. Latent Variable Discovery:
        * Identifies high-variance autosomal features.
        * Imputes missing values (KNN, Iterative, or Simple).
        * Regresses out "age" and other known covariates (sex, ancestry, technical metrics) to focus on unknown variance.
        * Performs PCA on residuals and automatically selects the optimal number of components (using knee-point detection on R² and RMSE curves) to capture unobserved technical/biological noise.
    3. Covariate Correction: Regresses out all non-target covariates (known technical factors + latent PCA components), leaving only age-related variance in the final residuals.
    4. Scaling: Scales the residuals (MinMax) to a 0-1 range to standardize the data distribution.
    5. Output: Saves the "clean" residuals and updated covariate files for downstream regression analysis.
  Key Technical Details:
    * Modeling: Relies on sklearn for PCA, imputation, and linear regression.

2. Analyze variance partition (Optional), done per cell-type and modality: run_variance_partition.py
  Core Workflow:
    1. Data Loading: Loads covariates and quantified data.
    2. Initial Variance Partition: Models gene expression variance against known covariates (e.g., age, sex, BMI, pool) to understand initial drivers of variation.
    3. Final Variance Partition: If final covariates from prep_pb_data.py exist, re-evaluates variance partitioning including the newly generated PCA components.
  Key Technical Details:
    * Parallelization: Uses concurrent.futures.ProcessPoolExecutor for computationally intensive variance partitioning.
    * Modeling: Relies on statsmodels for GLM/correlation checks.
    * Visualization: Generates heatmaps (covariate vs. PCA correlations) and boxen plots (variance explained distributions) to quality check the process.

1. Identify gene expression and chromatin accessibility features associated with age per broad cell-type and cluster specfic cell-type
    - Convert the single-cell data to pseudobulk (mean) values for each broad and cluster specific cell-type for both GEX and ATAC data; pseudobulk_convert.py
    - Format covariate tables for use with data prep and regression analysis; format_covariates.py
    - Prepare the pseudobulk data for analysis by removing non-target variable variance; prep_pb_data.py
    - Analyze sources of variance in the data (Optional); run_variance_partition.py
    - Regression analysis between quantified features (expression and accessibility) and age; pseudobulk_regression_analysis.ipynb. Where the possible regression methods include GLM, GLM with Tweedie distribution, and RLM.
    - Post-proceesing of the regression analysis across cell-types to apply B&H FDR and identify the statistally significant linear correlations between feature quantification and age; post_pseudobulk_regression.ipynb
    - Filter outlier effects from age regression from the GLM Tweedie results based on the RLM results; filter_regression_type_differences.ipynb
1. Identification of Cis-Correlated Features (cis_correlation.py)
   * Objective: To identify significant correlations between endogenous modalities (e.g., RNA expression) and proximal exogenous modalities (e.g.,
     ATAC-seq peaks) within a defined genomic window.
   * Data Preparation:
       * Feature Filtration: Feature Filtration: The analysis begins by loading significant endogenous features that have been previously identified as being associated with age (FDR < 0.05) and all detected exogenous features.
       * Proximity Mapping: For each cell type, endogenous features (e.g., genes) are mapped to proximal exogenous features (e.g., peaks) that reside
         within a specified maximum genomic distance (default: 1 Mb). This establishes the candidate cis-regulatory pairs to be tested.
       * Covariate Integration: Biological and technical covariates (e.g., age, BMI, PMI, pH, sex, ancestry, smoking status, and sequencing pool) are
         harmonized across both modalities. Cell-type-specific sample counts are incorporated.
   * Statistical Modeling:
       * Regression Framework: For each identified cis-proximal pair within each cell type, a regression model (e.g., Weighted Least Squares, WLS) is
         applied to evaluate the relationship between the endogenous and exogenous quantifications.
       * Model Formula: The regression model controls for the unified set of covariates. The formula takes the form: Endogenous_Feature ~
         Exogenous_Feature + Covariates.
       * Parallelization: Regressions are executed in parallel across cell types to optimize computational efficiency.
   * Significance Thresholding:
       * FDR Calculation: The resulting p-values for the exogenous terms across all tested pairs are adjusted for multiple comparisons using the
         Benjamini-Hochberg (BH) procedure.
       * Output Generation: The pipeline outputs both the complete unfiltered regression results and a subset filtered at a predefined false discovery
         rate (e.g., FDR ≤ 0.05).
5. Cis-Conditioned Age Regression (cis_conditioned_regression.py)
   * Objective: To evaluate whether the age effect on an endogenous feature (e.g., gene expression) persists when conditioning on a proximal, cis-correlated exogenous feature (e.g., ATAC peak) where both are independently age-associated.
   * Data Selection:
       * Feature Restraint: Loads previously computed cis-correlations and restricts the analysis solely to pairs where both the endogenous and exogenous features are significantly associated with age (FDR < 0.05).
       * Re-evaluates FDR strictly on this subset of dual-age-associated pairs to establish the baseline of significantly correlated pairs to test.
   * Statistical Modeling:
       * Outcome Model: Implements a Weighted Least Squares (WLS) regression model mirroring the outcome phase of a mediation analysis.
       * Formula: `Endogenous ~ Age + Exogenous + [Covariates]`. The model controls for covariates (e.g., sex, technical metrics) from both modalities and is weighted by the endogenous feature's sample counts.
   * Execution & Outputs:
       * The analysis is parallelized across cell types.
       * A new FDR correction is applied to the age exposure term.
       * Outputs full results and an FDR-filtered list (significantly conditioned pairs) into the `results` directory.
       * Logs the total pairs tested, nominally significant pairs, and the final FDR-significant count.

6. Latent Factor Generation and Age Association Pipeline
   * Objective: To discover coordinated networks of features (latent factors) within single-cell data using Consensus Non-negative Matrix Factorization (cNMF), and to subsequently evaluate these factors for significant associations with chronological age using linear mixed-effects modeling.
   * Latent Factor Generation (cnmf_latent_generation.py):
       * Data Preparation: Single-cell data is subset by cell type and filtered to retain only cells from valid donors and features used for age-associated regression analysis. Highly variable genes (HVGs) are pre-calculated to ensure stability during factorization.
       * Covariate Integration: If specified, the cNMF Preprocess module is utilized to regress out known technical covariates via Harmony prior to factorization.
       * Factorization (cNMF): The pipeline iterates across a defined range of components (K), or dynamically determines the range based on the number of unique donors. It performs multi-processed non-negative matrix factorization, combines the replicates, and executes consensus clustering to define stable latent factors (spectra) and cellular usage scores.
   * Mixed-Effects Regression (cnmf_latent_regressions.py):
       * Optimal K Selection: Automatically identifies the optimal number of components (K) by minimizing the distance to the ideal trade-off point between clustering stability (silhouette score) and prediction error, unless a specific K is manually provided.
       * Statistical Modeling: For each latent factor defined at the selected K, a Linear Mixed Effects Model (MixedLM) is fitted to test the association between cellular usage scores and age. 
       * Formula Framework: The model takes the form "latent_factor ~ age + [covariates]", incorporating the donor sample ID as a random intercept to account for donor-level baseline variations. Collinear covariates are dynamically identified and dropped to ensure model convergence.
   * Post-Processing and Feature Extraction (post_cnmf_latent_regressions.py):
       * Significance Thresholding: Regression results are aggregated across all tested cell types. P-values for the age coefficient are adjusted for multiple comparisons using the Benjamini-Hochberg (BH) procedure to generate a False Discovery Rate (FDR). Factors with an FDR <= 0.05 are deemed significantly associated with age.
       * Dynamic Feature Selection: For each significant latent factor, the corresponding cNMF feature spectra scores are extracted. A knee-point detection algorithm is applied to the descending curve of positive feature scores to empirically determine the inflection point (elbow). Features with scores exceeding this data-driven threshold are selected as the "top features" driving the latent network.
       * Visualization: Automatically generates comprehensive scatter plots detailing the spectra scores for all features within significant factors. Selected top features are highlighted, and the top five highest-scoring features are explicitly annotated utilizing collision-avoidance algorithms to prevent label overlap. Outputs are saved as both long-format CSVs and high-resolution figures.
