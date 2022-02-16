## EXAMPLE OF AN OLD SCRIPT


# This script is intended to be applied to data that have already been processed through 
# the Ingalls Lab targeted pipeline:
#
# 1. Import of raw data csvs
# 2. Quality control
# 3. BMIS
# 4. If applicable, quantification
#
# Once the data are a single, processed, long dataframe, they can be put through this script,
# which will format them for upload to CMAP.

source("src/Functions.R")

library(tidyverse)
options(scipen = 999)

file.pattern <- "BMISd"

# Import all files --------------------------------------------------
filenames <- RemoveCsv(list.files(path = "data_processed/", pattern = file.pattern))

for (i in filenames) {
  filepath <- file.path("data_processed/", paste(i,".csv", sep = ""))
  assign(make.names(i), read.csv(filepath, stringsAsFactors = FALSE, check.names = TRUE) %>%
           select(-X) %>%
           rename(compound = Mass.Feature))
}

original.template <- read.csv("data_raw/SCOPE_B12Incubation_Example.csv", 
                        stringsAsFactors = FALSE)


# Combine averages --------------------------------------------------------
averages.combined <- BMISd_CmpdNamesFixed %>%
  separate(Replicate.Name, into = c("Date", "runtype", "Supergroup", "replicate"), remove = FALSE) %>%
  group_by(compound, Supergroup) %>%
  mutate(Adjusted.Area = mean(Adjusted.Area, na.rm = TRUE)) %>%
  select(compound, Replicate.Name, Supergroup, Adjusted.Area)
  
  
# Add treatment type ------------------------------------------------------
add.treatment <- averages.combined %>%
  mutate(Control.Status = ifelse(str_detect(Supergroup, "IT0"),
                                 "Incubation", ifelse(str_detect(Supergroup, "DSW"), "DeepSeaWater", 
                                                      ifelse(str_detect(Supergroup, "Control"), "Control", "Treatments")))) %>%
  mutate(treatment_type = ifelse(Control.Status == "Control", "Control",
                                   ifelse(Control.Status == "DeepSeaWater", "DeepSeaWater",
                                          ifelse(Control.Status == "Incubation", "Time0",
                                                 ifelse(str_detect(Supergroup, "DMBnoBT"), "DMBnoB12",
                                                        ifelse(str_detect(Supergroup, "WBT"), "B12",
                                                               ifelse(str_detect(Supergroup, "DMB"), "DMB", "noB12"))))))) %>%
  select(c(compound, treatment_type, Adjusted.Area)) %>%
  unique()

# Add eddy vorticity, lat, lon, size fraction ------------------------------------------------------
add.variables <- add.treatment %>%
  mutate(eddy_vorticity = ifelse(str_detect(Supergroup, "IL1"), "Cyclonic", "Anticyclonic"),
         lat = ifelse(str_detect(Supergroup, "IL1"), 25, 27),
         lon = ifelse(str_detect(Supergroup, "IL1"), 158.5, 157.0),
         depth = ifelse(treatment_type == "Time0", 15, 25),
         size_fraction_um = ifelse(str_detect(Supergroup, "5um"), 5, 0.2)) 


# Add time ----------------------------------------------------------------

time.df <- original.template %>%
  filter(compound == "2-Hydroxy-4-(methylthio)butyric acid") %>%
  select(time)

split.by.compound <- split(add.variables, f = add.variables$compound)
add.time <- lapply(split.by.compound, function(x){
  y <- x %>% cbind(time.df)
  return(y)
})

full.data <- bind_rows(add.time, .id = "compound")

final.data <- full.data %>%
  rename(value = Adjusted.Area) %>%
  ungroup() %>%
  select(time, lat, lon, depth, compound, treatment_type, eddy_vorticity, size_fraction_um, value) %>%
  as.data.frame()

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

final.data[is.nan(final.data)] <- NA

view(final.data)