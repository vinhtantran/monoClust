## code to prepare `wind_sensit_2007` dataset goes here
library(readr)
library(dplyr)

wind <- read_csv("data-raw/sensit.csv",
                 col_types = cols(
                   LOC = col_factor(levels=NULL),
                   WDIR = col_factor(levels=NULL),
                   Date = col_datetime(format = "%m/%d/%Y %H:%M")
                 )) %>%
  select(LOC, Date, Sensit, WS, WDIR = WDIR....)

LOCATION <- "BR"
TIMERANGE1 <- c("2007-08-04", "2007-08-11")

wind_subset_2007 <- wind %>%
  filter(LOC == LOCATION,
         Date >= TIMERANGE1[1] & Date <= TIMERANGE1[2]) %>%
  select(Sensit, WS, WDIR) %>%
  drop_na()

wind_sensit_2007 <- wind_subset_2007 %>%
  mutate(has.sensit = as.numeric(Sensit > 0)) %>%
  select(has.sensit, WS, WDIR)

usethis::use_data(wind_sensit_2007, overwrite = TRUE)
