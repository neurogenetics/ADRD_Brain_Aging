import scanpy as sc
from anndata import AnnData
from pandas import DataFrame
from tabulate import tabulate

# variables and constants
project = "aging_phase2"
DEBUG = True
TESTING = False
TEST_FEATURE_SIZE = 1000
var_modal_dict = {"Gene Expression": "rna", "Peaks": "atac"}
# directories
wrk_dir = "/mnt/labshare/raph/datasets/adrd_neuro/brain_aging/phase2"
quants_dir = f"{wrk_dir}/quants"
results_dir = f"{wrk_dir}/results"
figures_dir = f"{wrk_dir}/figures"
sc.settings.figdir = f"{figures_dir}/"

# in files
raw_anndata_file = f"{quants_dir}/{project}.raw.multivi_prep.h5ad"
anndata_file = f"{quants_dir}/{project}.multivi.annotated.h5ad"

# out files


# helper functions
def peek_anndata(adata: AnnData, message: str = None, verbose: bool = False):
    if message is not None and len(message) > 0:
        print(message)
    print(adata)
    if verbose:
        print(tabulate(adata.obs.head(), headers="keys", tablefmt="psql"))
        print(tabulate(adata.var.head(), headers="keys", tablefmt="psql"))


def peek_dataframe(df: DataFrame, message: str = None, verbose: bool = False):
    if message is not None and len(message) > 0:
        print(message)
    print(f"{df.shape=}")
    if verbose:
        print(tabulate(df.head(), headers="keys", tablefmt="psql"))


# load the annatated anndate file
print("loading the annotated multivi anndata file")
adata_annot = sc.read(anndata_file, cache=True)
if DEBUG:
    peek_anndata(adata_annot, f"loaded {anndata_file}", DEBUG)
# load the full raw anndata file
print("loading the full raw multivi anndata file")
adata_raw = sc.read(raw_anndata_file)
if DEBUG:
    peek_anndata(adata_raw, f"loaded {raw_anndata_file}", DEBUG)

print(adata_raw.obs.modality.value_counts())
print(adata_raw.obs.Study.value_counts())
print(adata_raw.var.modality.value_counts())

# 1. Create the aggregate AnnData
pb_adata = sc.get.aggregate(adata_annot, by=["sample_id", "cell_label"], func="sum")
peek_anndata(pb_adata, "pseudobulked anndata", DEBUG)

# 2. Loop through cell types for your regression analysis
unique_cell_types = pb_adata.obs["cell_label"].unique()
if DEBUG:
    print(unique_cell_types)

for ct in unique_cell_types:
    print(f"{ct=}")
    # Extract just the samples for this cell type
    ct_data = pb_adata[pb_adata.obs["cell_label"] == ct].copy()

    # Loop through modalities (rna/atac)
    for modal_full, modal_short in var_modal_dict.items():
        # Subset AnnData by modality using the var table
        ct_modal_data = ct_data[:, ct_data.var["modality"] == modal_full].copy()
        # find the cell counts per donor on the subset
        adata_cell_info = adata_annot[
            adata_annot.obs.cell_label == ct, adata_annot.var.modality == modal_full
        ]
        donor_counts = adata_cell_info.obs.groupby("sample_id").size()

        if ct_modal_data.n_vars > 0:
            # Determine min cells threshold based on modality
            min_cells = 10 if modal_short == "rna" else 20
            low_count_samples = list(
                donor_counts[donor_counts < min_cells].index.values
            )
            ct_low_count_samples = [f"{x}_{ct}" for x in low_count_samples]
            donor_mask = ct_modal_data.obs_names.isin(ct_low_count_samples)

            # Initialize X with sums
            ct_modal_data.X = ct_modal_data.layers["sum"].copy()

            # mask donors with low cell counts
            ct_modal_data.X[donor_mask, :] = 0

            # now filter features that aren't well detected
            min_cells_threshold = int(ct_modal_data.n_obs * 0.30)
            pre_filter_var_count = ct_modal_data.n_vars
            sc.pp.filter_genes(ct_modal_data, min_cells=min_cells_threshold)
            post_filter_var_counts = ct_modal_data.n_vars
            print(
                (
                    f"## for {ct} {len(low_count_samples)} samples had a low cell cell count and "
                    f"{pre_filter_var_count - post_filter_var_counts} features were filtered"
                )
            )

            # For RNA: CPM + Log1p
            # For ATAC: RPM + Log1p (standard for linear regression input)
            sc.pp.normalize_total(ct_modal_data, target_sum=1e6)
            sc.pp.log1p(ct_modal_data)

            # Convert to DataFrame for statsmodels/regression
            df_modal = ct_modal_data.to_df()

            # correct the scanpy agrregated index names since we are split by cell-type already
            df_modal.index = df_modal.index.str.replace(f"_{ct}", "")

            # Run regression...
            print(
                f"Pseudobulk converted {ct} {modal_short} with shape {df_modal.shape}"
            )

            # Save the converted data to file
            out_file = (
                f"{quants_dir}/{project}.{ct.replace(' ', '_')}.{modal_short}.parquet"
            )
            df_modal.to_parquet(out_file)
            print(f"Saved {ct} {modal_short} to {out_file}")
