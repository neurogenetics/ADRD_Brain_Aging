from pandas import DataFrame
from numpy import where, cumsum, arange
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.decomposition import PCA, FastICA, NMF
from kneed import KneeLocator
from umap import UMAP
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler


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
