## code to prepare `wind_sensit_2008` dataset goes here
library(readr)
library(dplyr)

wind <- read_csv("data-raw/sensit.csv",
                 col_types = cols(
                   LOC = col_factor(levels = NULL),
                   WDIR = col_factor(levels = NULL),
                   Date = col_datetime(format = "%m/%d/%Y %H:%M")
                 )) %>%
  select(LOC, Date, Sensit, WS, WDIR = WDIR....)

LOCATION <- "BR"
TIMERANGE2 <- c("2008-07-07", "2008-07-14")

wind_subset_2008 <- wind %>%
  filter(LOC == LOCATION,
         Date >= TIMERANGE2[1] & Date <= TIMERANGE2[2]) %>%
  select(Sensit, WS, WDIR) %>%
  drop_na()

wind_sensit_2008 <- wind_subset_2008 %>%
  mutate(has.sensit = as.numeric(Sensit > 0)) %>%
  select(has.sensit, WS, WDIR)

usethis::use_data(wind_sensit_2008, overwrite = TRUE)
