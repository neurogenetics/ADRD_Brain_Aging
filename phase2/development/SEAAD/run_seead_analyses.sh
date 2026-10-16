# perform latent based regression analyses
# generate latent features per-celltype using cNMF
CELLTYPES="Astro BEC_venous ExN_CUX2 ExN_LAMP5 ExN_RELN ExN_RORB ExN_SEMA3E ExN_THEMIS InN_LAMP5 InN_PAX6 InN_PVALB InN_SST InN_VIP Micro OPC Oligo Peri"
MODALITIES="rna atac"
for CELLTYPE in ${CELLTYPES[@]}; do
  for MODALITY in ${MODALITIES[@]}; do
    echo "$CELLTYPE" "$MODALITY"
    --input-file "$DATADIR"/public/seaad/seaad_ec_multiome_labeled.h5mu \
      uv run phase2/development/SEAAD/analyses/cnmf_latent_generation.py \
      --cell-type "$CELLTYPE" \
      --cell-type-col anno_coarse \
      --sample-col donor_id \
      --modality "$MODALITY" \
      --exclude-ids H21.33.018,H21.33.001,H20.33.036 \
      --workers 48
  done
done

# run regression analsysis for AD using the generated latents
for CELLTYPE in ${CELLTYPES[@]}; do
  # review cNMF stability figures and run latent based analysis using the determined K
  for MODALITY in ${MODALITIES[@]}; do
    echo "$CELLTYPE" "$MODALITY"
    uv run phase2/development/SEAAD/analyses/cnmf_latent_regressions.py \
      --modality "$MODALITY" \
      --cell-type "$CELLTYPE" \
      --cell-type-col anno_coarse \
      --sample-col donor_id \
      --k auto \
      --covariates PCA_0 PCA_1 PCA_2 PCA_3
  done
done

# combine the cNMF latent regression output and compute FDRs
uv run phase2/development/SEAAD/analyses/post_cnmf_latent_regressions.py --modality rna
uv run phase2/development/SEAAD/analyses/post_cnmf_latent_regressions.py --modality atac

# compare the age association latent factor with others across cell-types and modalities
uv run phase2/development/SEAAD/analyses/compare_target_latent_factors.py --cell-type-col anno_coarse --sample-col donor_id
