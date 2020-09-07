## code to prepare `arctic`_2019 dataset goes here
library(readr)
library(dplyr)

arctic_raw_colnames <- read_csv("ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/north/daily/data/N_seaice_extent_daily_v3.0.csv",
                       col_names = TRUE) %>%
  colnames()

arctic <- read_csv("ftp://sidads.colorado.edu/DATASETS/NOAA/G02135/north/daily/data/N_seaice_extent_daily_v3.0.csv",
                   col_names = FALSE, skip = 2) %>%
  set_colnames(arctic_raw_colnames)

arctic_2019 <- arctic %>%
  filter(Year <= 2019)

usethis::use_data(arctic_2019, overwrite = TRUE)
