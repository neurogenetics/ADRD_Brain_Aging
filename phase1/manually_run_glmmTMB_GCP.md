### instruction for installing R and needed packages for running glmmTMB.R in GCP debian instance

using n1-highcpu-96 machine type

## so what was working two days ago doesn't seem to work anymore for setting up the R environment, on new instances I get a seg fault when running glmmTMB. The instance I set up on Oct 6th are still running fine, so I removed the auto-update settings for those instances. Tried from both python and R 4.0 notebook instances, also get the same thing from lngnode7 setup.

# for adding R to the python instance
add 'deb http://cloud.r-project.org/bin/linux/debian buster-cran40/' to /etc/apt/sources.list
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys FCAE2A0E115C3D8A
sudo apt update
sudo apt list --upgradable
sudo apt-get install r-base r-base-dev libssl-dev libxml2 libxml2-dev

# once R is installed or if starting from an R 4.0 instance
R --version
update.packages(ask = FALSE, checkBuilt = TRUE, Ncpus = 6)
install.packages("data.table", Ncpus = 6)
install.packages("dplyr", Ncpus = 6)
install.packages("glmmTMB", Ncpus = 6)
install.packages("broom.mixed", Ncpus = 6)

mkdir results
gsutil -mq cp -r gs://nihnialng-aging-brain/analysis/phase1_de/* .

Then can run regions and cell-types, for example
CELLTYPES="Oligodendrocyte-1 Oligodendrocyte-2 Astrocyte OPC Radial_Glia Microglia Endothelial Astrocyte-GFAP-Hi"

CELLTYPES="ExN_THEMIS ExN_CUX2_LAMP5 ExN_FEZF2 ExN_RORB_THEMIS ExN_CUX2_ADARB2 ExN_RORB ExN_LAMP5"

CELLTYPES="SPN_D1 SPN_D1-2 SPN_D2 SPN_D2-2 InN_ADARB2_VIP InN_LHX6_PVALB InN_ADARB2_LAMP5 InN_LHX6_SST"

for CELLTYPE in ${CELLTYPES[@]}
do
nohup R --no-echo --no-restore --file=/home/jupyter/glmmTMB.R --args /home/jupyter/quants/${CELLTYPE}_glmmtmb_in_df_temp.csv /home/jupyter/results/aging.${CELLTYPE}_glmmtmb_results_temp.csv &
done