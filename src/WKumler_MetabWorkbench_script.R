library(tidyverse)

## Original code from WKumler

# Script to make files for Metabolomics Workbench upload
# Should be run from the *project* directory, not the metab_workbench_metadata folder


# QE data for GBT fate ----
## Study design table ----
all_ids <- c("pos", "neg", "cyano")
all_files <- sapply(paste0("mzXMLs_", all_ids), list.files, pattern="mzXML", full.names=TRUE)
all_file_df <- data.frame(filename=unlist(all_files))
study_table <- all_file_df %>%
  mutate(sample_id=str_replace(basename(filename), ".mzXML", "")) %>%
  mutate(column_type=ifelse(str_detect(filename, "^mzXMLs_cyano"), "CYANO", "HILIC")) %>%
  mutate(MS_instrument="QE HF Orbitrap") %>%
  mutate(MS_polarity=ifelse(str_detect(filename, "^mzXMLs_neg"), "neg", "pos")) %>%
  mutate(unique_samp_id=paste(sample_id, MS_polarity, column_type, sep = "_")) %>%
  mutate(sample_type=str_extract(filename, paste0("(?<=200612_|200605_).*?(?=_)"))) %>%
  mutate(sample_type=c(Smp="Sample", Blk="Blank", Poo="Pooled", Std="Standard")[sample_type]) %>%
  mutate(dilution_factor=ifelse(str_detect(filename, "Poo.*Half"), 0.5, 1)) %>%
  mutate(experiment_number=str_extract(filename, "(?<=Fate)\\d(?=M)")) %>%
  mutate(timepoint=str_extract(filename, "(?<=MT)\\d+|(?<=MT)long")) %>%
  mutate(timepoint=as.numeric(ifelse(timepoint=="long", 96, timepoint))) %>%
  mutate(has_IS=str_detect(filename, "-IS")) %>%
  mutate(triplicate_id=str_extract(filename, "[A-C](?=.mzXML)")) %>%
  select(`Raw file name`=filename, `Sample name`=unique_samp_id,
         `LC column type`=column_type, `MS instrument`=MS_instrument, 
         `MS polarity`=MS_polarity, 
         `Sample type`=sample_type, `Dilution factor`=dilution_factor, 
         `Experiment number`=experiment_number, `Timepoint`=timepoint, 
         `Has internal standards`=has_IS, `Triplicate ID`=triplicate_id)

tqs_files <- list.files("TQS_data_kinetics", pattern=".raw", full.names=TRUE, include.dirs = TRUE)
tqs_file_df <- data.frame(filename=unlist(tqs_files))
tqs_metadata <- read.csv("metab_workbench_metadata/raw_files/GBT Kintetics Gradients 3 Sample List TQS for R code use.csv")

tqs_study_table <- tqs_file_df %>%
  mutate(unique_samp_id=str_replace(basename(filename), ".raw", "")) %>%
  mutate(column_type="HILIC") %>%
  mutate(MS_instrument="Waters TQS triple-quadrupole") %>%
  mutate(MS_polarity="pos") %>%
  mutate(sample_type=str_extract(filename, "Smp|Std|Blk|Poo")) %>%
  mutate(dilution_factor=ifelse(str_detect(filename, "Poo.*Half"), 0.5, 1)) %>%
  mutate(sample_type=c(Smp="Sample", Blk="Blank", Poo="Pooled", Std="Standard")[sample_type]) %>%
  mutate(experiment_number=0) %>%
  mutate(timepoint=NA, has_IS=NA) %>%
  mutate(triplicate_id=str_extract(filename, "[A-C](\\d+)?(?=.raw)")) %>%
  right_join(tqs_metadata, by=c(unique_samp_id="Replicate.Name")) %>%
  select(`Raw file name`=filename, `Sample name`=unique_samp_id,
         `LC column type`=column_type, `MS instrument`=MS_instrument, 
         `MS polarity`=MS_polarity, 
         `Sample type`=sample_type, `Dilution factor`=dilution_factor, 
         `Experiment number`=experiment_number, `Timepoint`=timepoint, 
         `Has internal standards`=has_IS, `Triplicate ID`=triplicate_id)

study_table <- rbind(study_table, tqs_study_table)


write.table(study_table, file = "metab_workbench_metadata/created_files/study_design_table.txt", 
            row.names = FALSE, sep = "\t", quote = FALSE)


## Targeted metabolites ----
### Data ----
all_targ <- read.csv("metab_workbench_metadata/raw_files/nM.concentrations.allFractions.csv")
clean_targ <- all_targ %>%
  select(`Raw file name`=Replicate.Name, compound_name=Precursor.Ion.Name, 
         area=rawArea, nm=nM.in.sample) %>%
  mutate(`Raw file name`=paste0(`Raw file name`, ".mzXML")) %>%
  left_join(study_table %>% mutate(`Raw file name`=basename(`Raw file name`)))
pos_targ <- clean_targ %>%
  filter(`LC column type`=="HILIC" & `MS polarity`=="pos") %>%
  select(`Sample name`, `Metabolite name`=compound_name, nm) %>%
  pivot_wider(names_from = `Sample name`, values_from = nm)
write.table(pos_targ, "metab_workbench_metadata/created_files/DATASET_2_HILIC_POS_targ.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)
neg_targ <- clean_targ %>%
  filter(`LC column type`=="HILIC" & `MS polarity`=="neg") %>%
  select(`Sample name`, `Metabolite name`=compound_name, nm) %>%
  pivot_wider(names_from = `Sample name`, values_from = nm)
write.table(neg_targ, "metab_workbench_metadata/created_files/DATASET_3_HILIC_NEG_targ.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)
cyano_targ <- clean_targ %>%
  filter(`LC column type`=="CYANO") %>%
  select(`Sample name`, `Metabolite name`=compound_name, nm) %>%
  pivot_wider(names_from = `Sample name`, values_from = nm)
write.table(cyano_targ, "metab_workbench_metadata/created_files/DATASET_4_CYANO_POS_targ.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)
### Metadata ----
stan_metadata_url <- "https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/430992e10c0464a817c69d22e3905d73f2ea196f/Ingalls_Lab_Standards.csv"
stan_metadata <- read.csv(stan_metadata_url) %>% 
  select(compound_name=Compound_Name, pubchemcode=PubChem_Code, keggcode=KEGG_Code,
         mz=mz, rt=RT_minute, polarity=z, `LC column type`=Column) %>%
  mutate(polarity=ifelse(as.numeric(polarity)<0, "neg", "pos")) %>%
  mutate(`LC column type`=ifelse(`LC column type`=="RP", "CYANO", `LC column type`))
cmpd_metadata <- clean_targ %>%
  distinct(compound_name, `MS polarity`, `LC column type`) %>%
  mutate(compound_name_clean=str_replace(compound_name, "_13C-0( 15N-0)?", "")) %>%
  left_join(stan_metadata, by=c("compound_name_clean"="compound_name", 
                                `MS polarity`="polarity", "LC column type")) %>%
  select(`Metabolite name`=compound_name, `PubChem CID`=pubchemcode, 
         `KEGG entry`=keggcode, `m/z ratio`=mz, `Retention time in min`=rt, 
         `MS polarity`, `LC column type`)

pos_meta <- cmpd_metadata %>% filter(`LC column type`=="HILIC" & `MS polarity`=="pos")
write.table(pos_meta, "metab_workbench_metadata/created_files/DATASET_2_HILIC_POS_targ_METADATA.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)
neg_meta <- cmpd_metadata %>% filter(`LC column type`=="HILIC" & `MS polarity`=="neg")
write.table(neg_meta, "metab_workbench_metadata/created_files/DATASET_3_HILIC_NEG_targ_METADATA.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)
cyano_meta <- cmpd_metadata %>% filter(`LC column type`=="CYANO")
write.table(cyano_meta, "metab_workbench_metadata/created_files/DATASET_4_CYANO_POS_targ_METADATA.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)



## Untargeted metabolites ----
# Data must be in wide format
# First column: metabolite name, formatted as 'm/z underscore retention time', e.g. 645.5327_24.91
# First row: name of metabolite column, followed by all the sample names 
# exactly matching those that you submitted in the 'Study Design' section
# See https://www.metabolomicsworkbench.org/data/ds_examples.php
# Files are the raw peaklists exported from XCMS after the "peakpicking" step

### HILIC pos ----
pos_exp2 <- read.csv("output_folder_pos/raw_peaks.csv") %>%
  group_by(feature) %>%
  mutate(meanmz=round(mean(mz), 5), meanrt=round(mean(rt), 2)) %>%
  mutate(`Metabolite name`=paste(meanmz, meanrt, sep = "_")) %>%
  ungroup() %>%
  select(`Metabolite name`, filename, into) %>%
  left_join(study_table %>% mutate(filename=basename(`Raw file name`))) %>%
  select(`Metabolite name`, `Sample name`, into) %>%
  arrange(desc(into)) %>%
  group_by(`Metabolite name`, `Sample name`) %>%
  slice(1)
pos_exp1 <- read.csv("output_folder_pos_exp1/raw_peaks.csv") %>%
  group_by(feature) %>%
  mutate(meanmz=round(mean(mz), 5), meanrt=round(mean(rt), 2)) %>%
  mutate(`Metabolite name`=paste(meanmz, meanrt, sep = "_")) %>%
  ungroup() %>%
  select(`Metabolite name`, filename, into) %>%
  left_join(study_table %>% mutate(filename=basename(`Raw file name`))) %>%
  select(`Metabolite name`, `Sample name`, into) %>%
  arrange(desc(into)) %>%
  group_by(`Metabolite name`, `Sample name`) %>%
  slice(1)
all_pos <- pos_exp1 %>%
  rbind(pos_exp2) %>%
  pivot_wider(names_from = `Sample name`, values_from = into)
write.table(all_pos, "metab_workbench_metadata/created_files/DATASET_2_HILIC_POS_untarg.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)


### HILIC neg ----
neg_exp2 <- read.csv("output_folder_neg/raw_peaks.csv") %>%
  group_by(feature) %>%
  mutate(meanmz=round(mean(mz), 5), meanrt=round(mean(rt), 2)) %>%
  mutate(`Metabolite name`=paste(meanmz, meanrt, sep = "_")) %>%
  ungroup() %>%
  select(`Metabolite name`, filename, into) %>%
  left_join(study_table %>% mutate(filename=basename(`Raw file name`))) %>%
  select(`Metabolite name`, `Sample name`, into) %>%
  arrange(desc(into)) %>%
  group_by(`Metabolite name`, `Sample name`) %>%
  slice(1)
neg_exp1 <- read.csv("output_folder_neg_exp1/raw_peaks.csv") %>%
  group_by(feature) %>%
  mutate(meanmz=round(mean(mz), 5), meanrt=round(mean(rt), 2)) %>%
  mutate(`Metabolite name`=paste(meanmz, meanrt, sep = "_")) %>%
  ungroup() %>%
  select(`Metabolite name`, filename, into) %>%
  left_join(study_table %>% mutate(filename=basename(`Raw file name`))) %>%
  select(`Metabolite name`, `Sample name`, into) %>%
  arrange(desc(into)) %>%
  group_by(`Metabolite name`, `Sample name`) %>%
  slice(1)
all_neg <- all_pos <- neg_exp1 %>%
  rbind(neg_exp2) %>%
  pivot_wider(names_from = `Sample name`, values_from = into)
write.table(all_neg, "metab_workbench_metadata/created_files/DATASET_3_HILIC_NEG_untarg.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)


### RP ----
cyano_exp2 <- read.csv("output_folder_cyano/raw_peaks.csv") %>%
  group_by(feature) %>%
  mutate(meanmz=round(mean(mz), 5), meanrt=round(mean(rt), 2)) %>%
  mutate(`Metabolite name`=paste(meanmz, meanrt, sep = "_")) %>%
  ungroup() %>%
  select(`Metabolite name`, filename, into) %>%
  left_join(study_table %>% mutate(filename=basename(`Raw file name`))) %>%
  select(`Metabolite name`, `Sample name`, into) %>%
  arrange(desc(into)) %>%
  group_by(`Metabolite name`, `Sample name`) %>%
  slice(1)
cyano_exp1 <- read.csv("output_folder_cyano_exp1/raw_peaks.csv") %>%
  group_by(feature) %>%
  mutate(meanmz=round(mean(mz), 5), meanrt=round(mean(rt), 2)) %>%
  mutate(`Metabolite name`=paste(meanmz, meanrt, sep = "_")) %>%
  ungroup() %>%
  select(`Metabolite name`, filename, into) %>%
  left_join(study_table %>% mutate(filename=basename(`Raw file name`))) %>%
  select(`Metabolite name`, `Sample name`, into) %>%
  arrange(desc(into)) %>%
  group_by(`Metabolite name`, `Sample name`) %>%
  slice(1)
all_cyano <- cyano_exp1 %>%
  rbind(cyano_exp2) %>%
  group_by(`Sample name`, `Metabolite name`) %>%
  arrange(desc(into)) %>%
  slice(1) %>%
  pivot_wider(names_from = `Sample name`, values_from = into)
write.table(all_cyano, "metab_workbench_metadata/created_files/DATASET_4_CYANO_POS_untarg.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)


# TQS data for GBT kinetics ----
## Peak areas ----
peak_areas <- read.csv("metab_workbench_metadata/raw_files/QC_output_All_TQS_GBT_IS-GBT.csv")
TQS_study_table <- study_table %>%
  filter(`MS instrument`=="Waters TQS triple-quadrupole")
clean_TQS_areas <- peak_areas %>%
  full_join(TQS_study_table, by=c("Replicate.Name"="Sample name")) %>%
  select(`Metabolite name`=Compound.Name, `Sample name`=Replicate.Name, rawArea) %>%
  mutate(`Metabolite name`=str_replace(`Metabolite name`, "13C-15N Glycine Betaine", "13C5 15N1 Glycine Betaine")) %>%
  arrange(desc(rawArea)) %>%
  group_by(`Metabolite name`, `Sample name`) %>%
  slice(1) %>%
  pivot_wider(names_from = `Sample name`, values_from = rawArea) %>%
  head(-1)
write.table(clean_TQS_areas, "metab_workbench_metadata/created_files/DATASET_1_HILIC_POS_TQS.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)

## Compound metadata ----
gbt_metadata <- stan_metadata %>% 
  filter(compound_name=="Glycine betaine"&`LC column type`=="HILIC") %>%
  mutate(compound_name=str_to_title(compound_name))
tqs_meta <- clean_TQS_areas %>%
  select(`Metabolite name`) %>%
  mutate(`Metabolite name`=str_replace(`Metabolite name`, "13C-15N Glycine Betaine", "13C5 15N1 Glycine Betaine")) %>%
  left_join(gbt_metadata, by=c(`Metabolite name`="compound_name"))
write.table(tqs_meta, "metab_workbench_metadata/created_files/DATASET_1_HILIC_TQS_POS_targ_METADATA.txt",
            row.names = FALSE, sep = "\t", quote = FALSE)


# Raw data gzip folder creation (takes about a minute per file) ----
study_table %>%
  group_by(`LC column type`, `MS instrument`, `MS polarity`) %>%
  group_split() %>%
  sapply(function(zip_df){
    files_to_zip <- zip_df$`Raw file name`
    output_dirname <- paste0(
      "metab_workbench_metadata/created_files/",
      paste(unique(zip_df$`LC column type`), 
            gsub(" ", "_", unique(zip_df$`MS instrument`)), 
            toupper(unique(zip_df$`MS polarity`)),
      sep = "_")
    )
    zip(zipfile = output_dirname, files = files_to_zip)
  })
