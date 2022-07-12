#!Rscript

# options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

library(glmmTMB)
library(data.table)
library(dplyr)
library(broom.mixed)

in_file <- as.character(args[1])
out_file <-  as.character(args[2])

# this_df = fread('/labshare/raph/datasets/adrd_neuro/aging/demux/putamen_df.csv')
this_df = fread(in_file)
core_cnt = parallel::detectCores() #/2)
if (core_cnt > 48) {
    core_cnt = 48
}
ctrl <- glmmTMB::glmmTMBControl(parallel=core_cnt)

compute_feature_mixmodel <- function(feature, df, ctrl) {
    if (grepl('-', feature)) {
      ret_value = NULL
    }
    else {
      this_formula = as.formula(sprintf('%1s ~ old + pool_name + (1 | Sample_id)', feature))
      this_result <- glmmTMB(this_formula, data=df, family=glmmTMB::tweedie, 
      ziformula= ~0, control = ctrl)
      ret_value = tidy(this_result)
      ret_value$feature = feature
    }
    return(ret_value)
}

features = colnames(this_df)[2:(dim(this_df)[2]-9)]

Sys.time()
fitmixed <- lapply(features, compute_feature_mixmodel, this_df, ctrl) %>% bind_rows()
Sys.time()

write.table(fitmixed,out_file, quote=F, row.names=F, sep = ",")
