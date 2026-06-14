library(tidyverse)
library(readr)

setwd("/data/ADRD/brain_aging/")

SEA_AD_file_manifest_metadata <- read_csv("./exploration/data/Phase2/SEAAD_EC/SEA-AD_file_manifest_metadata.csv")
fastqs=SEA_AD_file_manifest_metadata%>%filter(tissue=="lateral entorhinal cortex"& fileFormat=="fastq")
glimpse(fastqs)
length(unique(fastqs$donor_id))
length(unique(fastqs$specimen_id))

assay_meta <- read.csv('./exploration/data/Phase2/SEAAD_EC/SEA-AD_assay_multiome_metadata.csv')
assay_meta<-assay_meta%>%
  filter(specimenID%in%fastqs$specimen_id)#only LEC - but then only 34??? should be 37 people, 40 specimens? - I guess only 34 survived their cutoffs, even though they had these extra samples....

# donor_meta <- read.csv('./exploration/data/Phase2/SEAAD_EC/SEA-AD_individual_metadata.csv')#this has 38 columns
donor_meta <- read.csv('./exploration/data/Phase2/SEAAD_EC/SEA-AD_individual_metadata_harmonized.csv') #this added 4 columns?
donor_meta<-donor_meta%>%
  filter(individualID%in%fastqs$donor_id)#so this still says 37... i guess 3 will fail some QC somewhere if they didn't consider it in their sample table??

sample_meta <- read.csv('./exploration/data/Phase2/SEAAD_EC/SEA-AD_biospecimen_metadata.csv') 
sample_meta<- sample_meta%>%
  filter(individualID%in% donor_meta$individualID & specimenID%in%fastqs$specimen_id)



##############################
##############################
##############################

fastq_map <- fastqs %>%
  mutate(
    base = str_remove(file_name, "_S\\d+_L\\d{3}_[A-Z0-9]+_001\\.fastq\\.gz$"),
    read_type = str_match(file_name, "_(I1|I2|R1|R2|R3)_001\\.fastq\\.gz$")[, 2],
    library_type = case_when(
      read_type == "R3" ~ "ATAC",
      read_type %in% c("R1", "R2") ~ "possibly_GEX_or_ATAC",
      read_type %in% c("I1", "I2") ~ "index",
      TRUE ~ "unknown"
    )
  )

fastq_summary <- fastq_map %>%
  group_by(specimen_id, donor_id, assay, base, library_type) %>%
  summarise(
    files = paste(sort(file_name), collapse = "\n"),
    reads = paste(sort(unique(read_type)), collapse = ", "),
    n_files = n(),
    .groups = "drop"
  ) %>%
  arrange(specimen_id, base, library_type)

fastq_summary

library_check <- fastq_map %>%
  group_by(specimen_id, base) %>%
  summarise(
    has_r3 = any(read_type == "R3"),
    has_only_r1_r2 = all(read_type %in% c("R1", "R2", "I1", "I2")),
    reads_present = paste(sort(unique(read_type)), collapse = ", "),
    .groups = "drop"
  )

library_check

##############################
##############################
##############################
files <- list.files(
  "/data/ADRD/brain_aging/exploration/data/Phase2/SEAAD_EC/SEAAD_FASTQ",
  recursive = TRUE,
  full.names = FALSE
)

length(files)

head(files)

fastq_tbl <- tibble(file_name = files) %>%
  mutate(
    mtx = str_extract(file_name, "MX\\d+"),
    sequencingBatch_fastq = str_replace(mtx, "^MX", "MTX-")
  )

# count(fastq_tbl, sequencingBatch_fastq, sort = TRUE)
count(assay_meta, sequencingBatch, sort = TRUE)


##############################
##############################
##############################

fastq_tbl <- fastq_tbl %>%
  mutate(
    sample_prefix =
      str_remove(
        file_name,
        "_S\\d+_L\\d{3}_(I1|I2|R1|R2|R3)_001.fastq.gz"
      ),
    read =
      str_match(
        file_name,
        "_(I1|I2|R1|R2|R3)_001.fastq.gz"
      )[,2]
  )

fastq_tbl %>%
  count(sample_prefix, sort = TRUE)

##############################
##############################
##############################
#what types of files are there per sample
fastq_tbl %>%
  group_by(sample_prefix) %>%
  summarise(
    reads = paste(sort(unique(read)), collapse = ","),
    n_files = n()
  ) %>%
  arrange(sample_prefix)

#it looks like all the ATAC are coming first (NY-AT) then later will be the GEX (NY-MX)

##############################
##############################
##############################

