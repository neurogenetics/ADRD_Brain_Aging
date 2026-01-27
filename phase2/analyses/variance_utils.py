from pandas import DataFrame, Series
from numpy import where, cumsum, arange
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.decomposition import PCA, FastICA, NMF
from kneed import KneeLocator
from umap import UMAP
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
import statsmodels.formula.api as smf
import logging
import warnings

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


# small function to perform the and minmax scale on a pandas dataframe
def scale_dataframe(this_df: DataFrame):
    scaledX = MinMaxScaler().fit_transform(this_df)
    scaled_df = DataFrame(data=scaledX, columns=this_df.columns, index=this_df.index)
    return scaled_df


# function for dimensionality reduction analyses
def iterate_model_component_counts(
    max_count: int, data_df: DataFrame, model_type: str = ["PCA", "NMF", "ICA"]
) -> (list, list):
    r2_rets = []
    rmse_rets = []
    for comp_num in arange(1, max_count + 1):
        _, _, r2, rmse = generate_selected_model(comp_num, data_df, model_type)
        r2_rets.append(r2)
        rmse_rets.append(rmse)
    return r2_rets, rmse_rets


def generate_selected_model(
    n_comps: int, data_df: DataFrame, model_type: str = ["PCA", "NMF", "ICA"]
) -> (object, DataFrame, float, float):
    if model_type == "PCA":
        model = PCA(n_components=n_comps, random_state=42)
    if model_type == "NMF":
        model = NMF(n_components=n_comps, init="random", random_state=42, max_iter=500)
    if model_type == "ICA":
        model = FastICA(n_components=n_comps, random_state=42)
    components = model.fit_transform(data_df)
    recon_input = model.inverse_transform(components)
    r2 = r2_score(y_true=data_df, y_pred=recon_input)
    rmse = mean_squared_error(data_df, recon_input, squared=False)
    print(
        f"{model_type} with {n_comps} components accuracy is {r2:.4f}, RMSE is {rmse:.4f}"
    )
    ret_df = DataFrame(data=components, index=data_df.index).round(4)
    ret_df = ret_df.add_prefix(f"{model_type}_")
    return model, ret_df, r2, rmse


def component_from_max_curve(scores, label: str = ["R2", "RMSE", "EVR"]) -> int:
    if label == "R2" or label == "EVR":
        data_curve = "concave"
        data_direction = "increasing"
    if label == "RMSE":
        data_curve = "convex"
        data_direction = "decreasing"
    knee = KneeLocator(
        arange(1, len(scores) + 1),
        scores,
        S=1.0,
        curve=data_curve,
        direction=data_direction,
    )
    print(f"best curve at knee {knee.knee}")
    num_comp = int(knee.knee)
    exp_value = scores[num_comp - 1]
    print(f"best number of components is {num_comp} at {label} of {exp_value}")
    knee.plot_knee()
    plt.show()
    knee.plot_knee_normalized()
    plt.show()
    return num_comp


# small function to generate umap from pandas dataframe,
# for all features (columns) \
# and return back as dataframe with source index intact
def generate_umap_covs_df(this_df, other_covs_df=None, rnd_digits=3, merge_input=False):
    # run UMAP on the data frame features
    umap_results = UMAP(random_state=42).fit_transform(this_df)
    umap_df = DataFrame(
        umap_results, columns=["x_umap", "y_umap"], index=this_df.index
    ).round(rnd_digits)
    if merge_input:
        umap_df = umap_df.merge(this_df, left_index=True, right_index=True)
    if other_covs_df is not None:
        umap_df = umap_df.merge(
            other_covs_df, how="left", left_index=True, right_index=True
        )
    print(f"The dimensions of the umap df and the traits are {umap_df.shape}")
    return umap_df


def get_high_variance_features(
    df: DataFrame, top_percent: float = 0.15, min_presence_percent: float = 0.8
) -> list[str]:
    """
    Identifies the top n% features with the highest variance,
    filtering out features with high missingness.
    
    Args:
        df: Input DataFrame (samples x features)
        top_percent: Percentage of top features to keep (0.0 to 1.0)
        min_presence_percent: Minimum fraction of non-null values required (0.0 to 1.0)
    """
    if top_percent <= 0 or top_percent > 1:
        logger.warning(f"Invalid percentage {top_percent}, using all features.")
        return df.columns.tolist()

    # 1. Filter based on missingness
    # count() returns non-NA count
    presence_fraction = df.count() / len(df)
    valid_features = presence_fraction[presence_fraction >= min_presence_percent].index
    
    n_dropped = len(df.columns) - len(valid_features)
    if n_dropped > 0:
        logger.info(f"Dropped {n_dropped} features with < {min_presence_percent*100:.0f}% presence.")
        
    if len(valid_features) == 0:
        logger.warning("No features met the presence criteria.")
        return []

    # 2. Calculate variance for valid features
    # df.var() skips NaNs by default, which is what we want
    variances = df[valid_features].var()

    # Sort descending
    sorted_vars = variances.sort_values(ascending=False)

    # Determine cutoff
    n_features = len(valid_features)
    n_keep = int(n_features * top_percent)
    n_keep = max(1, n_keep)  # Ensure at least one

    top_features = sorted_vars.head(n_keep).index.tolist()
    logger.info(
        f"Selected {len(top_features)} high-variance features (Top {top_percent*100:.1f}% of {n_features} valid features)"
    )

    return top_features


def perform_variance_partition(
    data_df: DataFrame,
    feature_name: str,
    fixed_effects: list[str],
    random_effects: list[str],
):
    """
    Mimics variancePartition by modeling continuous variables as fixed effects
    and categorical variables as random effects using statsmodels MixedLM.
    """
    try:
        # Create model columns list
        model_cols = fixed_effects + random_effects + [feature_name]

        # Create a temporary dataframe for modeling
        fit_df = data_df[model_cols].dropna()
        # Check if enough data remains
        if fit_df.empty:
            return None

        fit_df["dummy_group"] = 1

        # formula for fixed effects
        # Patsy handles basic formula parsing
        fixed_formula = f"Q('{feature_name}') ~ " + " + ".join(fixed_effects)

        # formula for random effects (categorical)
        vc_formula = {effect: f"0 + C({effect})" for effect in random_effects}

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # Fit the Linear Mixed Model
            model = smf.mixedlm(
                fixed_formula,
                data=fit_df,
                groups="dummy_group",
                re_formula="0",  # No random intercept for the dummy group itself
                vc_formula=vc_formula,
            )
            result = model.fit()
            # logger.info(result.summary()) # Silenced for batch processing

        # --- CALCULATE VARIANCE FRACTIONS ---

        # 1. Random Effects Variance (from vc_formula)
        try:
            vc_names = result.model.exog_vc.names
        except AttributeError:
            vc_names = [f"Var_Comp_{i}" for i in range(len(result.vcomp))]

        random_vars = dict(zip(vc_names, result.vcomp))

        # 2. Fixed Effects Variance
        # Var_explained = beta^2 * Var(X)
        fixed_vars = {}
        for term in fixed_effects:
            # Statsmodels params keys might differ slightly (e.g., Q('feature'))
            # But the predictors are in fixed_effects
            if term in result.params:
                beta = result.params[term]
                var_x = fit_df[term].var()
                fixed_vars[term] = (beta**2) * var_x

        # 3. Residual Variance
        residual_var = result.scale

        # 4. Combine and Normalize
        all_variances = {**random_vars, **fixed_vars, "Residual": residual_var}
        all_variances_series = Series(all_variances)

        total_var = all_variances_series.sum()
        fractions = all_variances_series / total_var

        return (feature_name, fractions)

    except Exception as e:
        # logger.warning(f"Failed to process {feature_name}: {e}")
        return None
