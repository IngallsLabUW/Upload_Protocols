library(tidyverse)

## This script is an example of how columns were added for CMAP's required metadata.
# The information below (like the depths for the stations) was taken from the SCOPE cruise logs.

MESOSCOPE_EddyTransect <- read_csv("data_processed/MESOSCOPE_Ingalls_EddyTransect_Dec2021.csv")

## Info for lat, lon, time, depth
Metadata <- MESOSCOPE_EddyTransect %>%
  select(Replicate.Name) %>%
  unique()

Depth <- Metadata %>%
  mutate(depth = case_when(str_detect(Replicate.Name, "MS4C1DCM") ~ "100",
                           str_detect(Replicate.Name, "MS5C1DCM") ~ "100",
                           str_detect(Replicate.Name, "MS6C3DCM") ~ "119",
                           str_detect(Replicate.Name, "MS7C1DCM") ~ "133",
                           str_detect(Replicate.Name, "MS8C1DCM") ~ "122",
                           str_detect(Replicate.Name, "MS9C2DCM") ~ "110",
                           str_detect(Replicate.Name, "MS10C2DCM") ~ "105",
                           str_detect(Replicate.Name, "MS11C2DCM") ~ "110",
                           str_detect(Replicate.Name, "MS12C1DCM") ~ "103",
                           str_detect(Replicate.Name, "MS13C2DCM") ~ "112",
                           str_detect(Replicate.Name, "MS14C2DCM") ~ "105",
                           str_detect(Replicate.Name, "Poo|Std|Blk") ~ "0",
                           str_detect(Replicate.Name, "15m") ~ "15",
                           str_detect(Replicate.Name, "175m") ~ "175"))
Lat <- Metadata %>%
  mutate(lat = case_when(str_detect(Replicate.Name, "MS4") ~ "28",
                         str_detect(Replicate.Name, "MS5") ~ "27.25",
                         str_detect(Replicate.Name, "MS6") ~ "26.5",
                         str_detect(Replicate.Name, "MS7") ~ "26.25",
                         str_detect(Replicate.Name, "MS8") ~ "26",
                         str_detect(Replicate.Name, "MS9") ~ "25.75",
                         str_detect(Replicate.Name, "MS10") ~ "25.5",
                         str_detect(Replicate.Name, "MS11") ~ "25.25",
                         str_detect(Replicate.Name, "MS12") ~ "25",
                         str_detect(Replicate.Name, "MS13") ~ "24.5",
                         str_detect(Replicate.Name, "MS14") ~ "24",
                         str_detect(Replicate.Name, "Poo|Std|Blk") ~ "0"))
Lon <- Metadata %>%
  mutate(lon = case_when(str_detect(Replicate.Name, "MS4") ~ "157",
                         str_detect(Replicate.Name, "MS5") ~ "157.5",
                         str_detect(Replicate.Name, "MS6") ~ "157.75",
                         str_detect(Replicate.Name, "MS7") ~ "157.9",
                         str_detect(Replicate.Name, "MS8") ~ "158",
                         str_detect(Replicate.Name, "MS9") ~ "158.25",
                         str_detect(Replicate.Name, "MS10") ~ "158.3",
                         str_detect(Replicate.Name, "MS11") ~ "158.4",
                         str_detect(Replicate.Name, "MS12") ~ "158.5",
                         str_detect(Replicate.Name, "MS13") ~ "158.75",
                         str_detect(Replicate.Name, "MS14") ~ "159",
                         str_detect(Replicate.Name, "Poo|Std|Blk") ~ "0"))
Time <- Metadata %>%
  mutate(time = case_when(str_detect(Replicate.Name, "MS4C115m") ~ "2017-06-29T02:06:00",
                          str_detect(Replicate.Name, "MS4C1DCM") ~ "2017-06-29T02:11:00",
                          str_detect(Replicate.Name, "MS4C1175m") ~ "2017-06-29T02:56:00",
                          str_detect(Replicate.Name, "MS5C115m") ~ "2017-06-29T13:21:00",
                          str_detect(Replicate.Name, "MS5C1DCM") ~ "2017-06-29T13:23:00",
                          str_detect(Replicate.Name, "MS5C1175m") ~ "2017-06-29T14:09:00",
                          str_detect(Replicate.Name, "MS6C315m") ~ "2017-06-29T00:30:00",
                          str_detect(Replicate.Name, "MS6C3DCM") ~ "2017-06-29T00:35:00",
                          str_detect(Replicate.Name, "MS6C3175m") ~ "2017-06-29T01:15:00",
                          str_detect(Replicate.Name, "MS7C115m") ~ "2017-06-30T02:41:00",
                          str_detect(Replicate.Name, "MS7C1DCM") ~ "2017-06-30T02:46:00",
                          str_detect(Replicate.Name, "MS7C1175m") ~ "2017-06-30T03:27:00",
                          str_detect(Replicate.Name, "MS8C115m") ~ "2017-06-30T09:20:00",
                          str_detect(Replicate.Name, "MS8C1DCM") ~ "2017-06-30T09:25:00",
                          str_detect(Replicate.Name, "MS8C1175m") ~ "2017-06-30T10:05:00",
                          str_detect(Replicate.Name, "MS9C215m") ~ "2017-06-30T14:24:00",
                          str_detect(Replicate.Name, "MS9C2DCM") ~ "2017-06-30T14:31:00",
                          str_detect(Replicate.Name, "MS9C2175m") ~ "2017-06-30T15:15:00",
                          str_detect(Replicate.Name, "MS10C215m") ~ "2017-06-30T21:32:00",
                          str_detect(Replicate.Name, "MS10C2DCM") ~ "2017-06-30T21:34:00",
                          str_detect(Replicate.Name, "MS10C2175m") ~ "2017-06-30T22:13:00",
                          str_detect(Replicate.Name, "MS11C215m") ~ "2017-07-01T02:32:00",
                          str_detect(Replicate.Name, "MS11C2DCM") ~ "2017-07-01T02:35:00",
                          str_detect(Replicate.Name, "MS11C2175m") ~ "2017-07-01T03:15:00",
                          str_detect(Replicate.Name, "MS12C115m") ~ "2017-07-01T06:48:00",
                          str_detect(Replicate.Name, "MS12C1DCM") ~ "2017-07-01T06:52:00",
                          str_detect(Replicate.Name, "MS12C1175m") ~ "2017-07-01T07:34:00",
                          str_detect(Replicate.Name, "MS13C315m") ~ "2017-07-01T16:03:00",
                          str_detect(Replicate.Name, "MS13C2DCM") ~ "2017-07-01T15:13:00",
                          str_detect(Replicate.Name, "MS13C2175m") ~ "2017-07-01T15:13:00",
                          str_detect(Replicate.Name, "MS14C215m") ~ "2017-07-01T21:29:00",
                          str_detect(Replicate.Name, "MS14C2DCM") ~ "2017-07-01T21:32:00",
                          str_detect(Replicate.Name, "MS14C2175m") ~ "2017-07-01T22:11:00",
                          str_detect(Replicate.Name, "Poo|Std|Blk") ~ "2017-06-28T00:51:00"))

Metadata <- Metadata %>%
  left_join(Lat, by = "Replicate.Name") %>%
  left_join(Lon, by = "Replicate.Name") %>%
  left_join(Depth, by = "Replicate.Name") %>%
  left_join(Time, by = "Replicate.Name")

Full_Data <- MESOSCOPE_EddyTransect %>%
  left_join(Metadata, by = "Replicate.Name") 