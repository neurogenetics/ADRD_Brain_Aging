"""
Python implementation of the R package variancePartition
Repository: https://github.com/GabrielHoffman/variancePartition

This module provides tools to quantify and interpret multiple sources of
biological and technical variation in gene expression experiments using
Linear Mixed Models (LMMs). It mimics the core functionality of the
fitExtractVarPartModel function from the original R package.
"""

import re
import warnings
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from patsy import dmatrices
from joblib import Parallel, delayed

# Optional: suppress statsmodels convergence warnings for clean output
warnings.filterwarnings("ignore", category=sm.tools.sm_exceptions.ConvergenceWarning)


def parse_r_formula(formula_str):
    """
    Parses an R-style mixed model formula into fixed and random components.
    e.g., '~ Age + (1|Batch) + Sex' -> '~ Age + Sex', ['Batch']
    """
    random_pattern = re.compile(r"\(\s*1\s*\|\s*([^)]+)\s*\)")
    random_effects = [x.strip() for x in random_pattern.findall(formula_str)]

    fixed_str = random_pattern.sub("", formula_str)

    if "~" in fixed_str:
        lhs, rhs = fixed_str.split("~", 1)
        rhs_parts = [p.strip() for p in rhs.split("+") if p.strip()]
        if not rhs_parts:
            rhs_parts = ["1"]
        fixed_formula = (
            f"{lhs.strip()} ~ {' + '.join(rhs_parts)}"
            if lhs.strip()
            else f"~ {' + '.join(rhs_parts)}"
        )
    else:
        rhs_parts = [p.strip() for p in fixed_str.split("+") if p.strip()]
        if not rhs_parts:
            rhs_parts = ["1"]
        fixed_formula = f"~ {' + '.join(rhs_parts)}"

    return fixed_formula, random_effects


def fit_single_gene(gene_name, y, metadata, fixed_formula, random_effects, maxiter=200):
    """
    Fits the LMM for a single gene and calculates variance fractions.
    """
    df = metadata.copy()
    df["Expression"] = y

    # Drop missing values
    df = df.dropna(
        subset=["Expression"]
        + [
            col
            for col in metadata.columns
            if col in fixed_formula or col in random_effects
        ]
    )

    if len(df) < 5:
        return gene_name, None

    df["_group"] = 1

    if "~" in fixed_formula:
        _, rhs = fixed_formula.split("~", 1)
        model_formula = f"Expression ~ {rhs}"
    else:
        model_formula = f"Expression {fixed_formula}"

    vc_formula = {re_name: f"0 + C({re_name})" for re_name in random_effects}

    try:
        y_design, X_design = dmatrices(model_formula, df, return_type="dataframe")
        design_info = X_design.design_info

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if vc_formula:
                model = smf.mixedlm(
                    model_formula,
                    df,
                    groups="_group",
                    vc_formula=vc_formula,
                    re_formula="0",
                )
                fit = model.fit(reml=False, maxiter=maxiter)
            else:
                model = smf.ols(model_formula, df)
                fit = model.fit()

        variances = {}

        # 1. Random Effects
        if vc_formula and hasattr(fit, "vcomp"):
            try:
                vc_names = fit.model.exog_vc.names
            except AttributeError:
                vc_names = random_effects
            
            # Use zip to map the NumPy array to the component names safely
            vcomp_dict = dict(zip(vc_names, fit.vcomp))
            for re_name in random_effects:
                variances[re_name] = vcomp_dict.get(re_name, 0.0)

        # 2. Residuals
        variances["Residuals"] = fit.scale

        # 3. Fixed Effects
        params = fit.params
        for term_name in design_info.term_names:
            if term_name == "Intercept":
                continue

            slice_idx = design_info.term_name_slices[term_name]
            col_names = design_info.column_names[slice_idx]

            beta_term = params[col_names]
            X_term = X_design[col_names]

            fitted_term = X_term.dot(beta_term.values)
            variances[term_name] = np.var(fitted_term, ddof=1)

        total_var = sum(variances.values())
        if total_var > 0:
            fractions = {k: v / total_var for k, v in variances.items()}
        else:
            fractions = {k: np.nan for k in variances}

        return gene_name, fractions

    except Exception as e:
        return gene_name, None


def fitExtractVarPartModel(expr_obj, formula, data, n_jobs=-1, maxiter=200):
    """
    Fits Linear Mixed Model to each feature (gene) and extracts variance fractions.

    Parameters:
    -----------
    expr_obj : pd.DataFrame
        Gene expression matrix (features x samples).
    formula : str
        R-style mixed model formula. e.g. '~ Age + (1|Batch) + Sex'
    data : pd.DataFrame
        Sample metadata (samples x covariates). Rows must align with expr_obj columns.
    n_jobs : int
        Number of parallel jobs to run. -1 means use all processors.

    Returns:
    --------
    pd.DataFrame
        Fractions of variance explained by each component for each feature.
    """
    fixed_formula, random_effects = parse_r_formula(formula)

    print(f"Parsed Formula:")
    print(f"  Fixed effects:  {fixed_formula}")
    print(f"  Random effects: {random_effects}")

    common_samples = expr_obj.columns.intersection(data.index)
    if len(common_samples) < len(expr_obj.columns):
        print(
            f"Warning: Dropped {len(expr_obj.columns) - len(common_samples)} samples due to missing metadata."
        )

    expr_data = expr_obj[common_samples]
    metadata = data.loc[common_samples]

    print(
        f"Fitting models for {len(expr_data)} features across {len(common_samples)} samples..."
    )

    try:
        from tqdm.auto import tqdm

        iterable = tqdm(expr_data.index, desc="Genes")
    except ImportError:
        iterable = expr_data.index

    results = Parallel(n_jobs=n_jobs)(
        delayed(fit_single_gene)(
            gene,
            expr_data.loc[gene].values,
            metadata,
            fixed_formula,
            random_effects,
            maxiter,
        )
        for gene in iterable
    )

    records = []
    for gene, res in results:
        if res is not None:
            res["Feature"] = gene
            records.append(res)
        else:
            nan_record = {re: np.nan for re in random_effects}
            nan_record["Residuals"] = np.nan
            nan_record["Feature"] = gene
            records.append(nan_record)

    df_res = pd.DataFrame(records).set_index("Feature")

    cols = [c for c in df_res.columns if c != "Residuals"] + ["Residuals"]
    return df_res[cols]


def plotVarPart(var_part_df, title="Variance Partition", filename="variance_partition.png"):
    """
    Generates a violin plot of variance fractions similar to the R package
    and saves it to a file. Requires seaborn and matplotlib.
    """
    try:
        import seaborn as sns
        import matplotlib.pyplot as plt
    except ImportError:
        print("seaborn and matplotlib are required for plotting.")
        return

    plot_df = var_part_df.reset_index().melt(
        id_vars=["Feature"], var_name="Component", value_name="Variance Fraction"
    )

    plt.figure(figsize=(10, 6))
    sns.violinplot(
        x="Component", y="Variance Fraction", data=plot_df, inner="quartile", cut=0
    )
    plt.xticks(rotation=45, ha="right")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Plot saved successfully to {filename}")


if __name__ == "__main__":
    # --- Example Usage ---
    np.random.seed(42)
    N_SAMPLES = 100
    N_GENES = 20

    # 1. Create Mock Metadata
    metadata = pd.DataFrame(
        {
            "Age": np.random.normal(50, 10, N_SAMPLES),
            "Disease": np.random.choice(["Control", "Case"], N_SAMPLES),
            "Batch": np.random.choice(["B1", "B2", "B3"], N_SAMPLES),
            "Tissue": np.random.choice(["Brain", "Blood"], N_SAMPLES),
        },
        index=[f"S{i}" for i in range(N_SAMPLES)],
    )

    # 2. Create Mock Expression Data
    # Genes (rows) x Samples (columns)
    expr = pd.DataFrame(
        np.random.randn(N_GENES, N_SAMPLES),
        columns=metadata.index,
        index=[f"Gene{i}" for i in range(N_GENES)],
    )

    # Add some structured variance
    expr.loc["Gene0", metadata["Batch"] == "B1"] += 2.0  # Batch effect
    expr.loc["Gene1"] += metadata["Age"].values * 0.1  # Age effect

    # 3. Define the mixed model formula
    # Categorical variables driven by grouping should be random effects: (1|Variable)
    # Continuous or categorical variables of direct biological interest are fixed effects
    formula = "~ Age + Disease + (1|Batch) + (1|Tissue)"

    # 4. Run fitExtractVarPartModel
    var_part = fitExtractVarPartModel(expr, formula, metadata, n_jobs=4)

    print("Variance Partitioning Results (First 5 Genes):")
    print(var_part.head())
    
    # 5. Save plot
    plotVarPart(var_part, filename="example_variance_partition.png")
