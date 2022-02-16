library(readxl)
library(tidyverse)

## Creation of the Metabolomcis Workbench study table for 2017 MESOSCOPE B12 Incubations
# If it’s a case where all the samples are the same and you’re applying 2 
# chromatography methods(RP and HILIC)
# to the same samples then you can incorporate all the analyses into 1 study. 
# Is filter size an experimental condition? 
# If so, you can use filter size (1 or 2) as a factor in your study-design table. 
# Then when you get to the chromatography step specify 2 methods (RP and HILIC). 
# When you get to the MS step for HILIC specify 2 MS methods (pos and neg) and 1 for RP. 
# Then you’ll have 6 different analyses within the same study.

## Example of a previous study table created by the Ingalls Lab:
Example_StudyTable <- read_excel("data_extra/StudyDesign_HD_20200907_Laura.xlsx")

## Incubation #1, 0.2um filter: Sample Key from Shared Google Drive
SampKey_Original <- read.csv("data_extra/MS_CyanoAq_B12Inc_SampKey.csv")

SampKey <- SampKey_Original %>%
  rename(Sample_ID = 1)

## Import list of raw filenames, obtained using shell scripts from the local drive.
# Note that the below list contains both .RAW and .mzXML files (only need the .mzXML),
# as well as unnecessary .sld files.
RawNames <- read.csv("data_raw/MS_B12Incubation1_RP.txt", header = FALSE)
