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

    # The aggregation put the data in layers['sum'], but normalize_total expects X
    ct_data.X = ct_data.layers["sum"].copy()

    # Now you have a standard AnnData for this cell type
    # You can compute CPM and Log2 here easily:
    sc.pp.normalize_total(ct_data, target_sum=1e6)  # CPM
    sc.pp.log1p(ct_data)  # Log2(CPM+1)

    # Loop through modalities (rna/atac)
    for modal_full, modal_short in var_modal_dict.items():
        # Subset AnnData by modality using the var table
        ct_modal_data = ct_data[:, ct_data.var["modality"] == modal_full]

        if ct_modal_data.n_vars > 0:
            # Convert to DataFrame for statsmodels/regression
            df_reg = ct_modal_data.to_df()

            # Run regression...
            print(f"Ready to regress {ct} {modal_short} with shape {df_reg.shape}")
            # peek_dataframe(df_reg, f"{ct} {modal_short} dataframes", DEBUG)

            # Save the converted data to file
            out_file = f"{quants_dir}/{project}.{ct.replace(' ', '_')}.{modal_short}.parquet"
            df_reg.to_parquet(out_file)
            print(f"Saved {ct} {modal_short} to {out_file}")
