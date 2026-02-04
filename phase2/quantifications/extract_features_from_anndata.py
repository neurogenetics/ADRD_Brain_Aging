from scanpy import read_h5ad
from tabulate import tabulate

debug = True

adata = read_h5ad(
    "../../datasets/adrd_neuro/brain_aging/phase2/quants/aging_phase2.raw.multivi_prep.h5ad"
)

print(adata)

if debug:
    print(tabulate(adata.var.sample(20), headers="keys", tablefmt="psql"))

adata.var.index.name = "gene"

adata.var.to_csv(
    "../../datasets/adrd_neuro/brain_aging/phase2/quants/aging_phase2.features.csv",
    index=True,
    header=True,
)
