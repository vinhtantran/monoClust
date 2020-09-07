#' NOAA's Arctic Sea Daily Ice Extend Data
#'
#' A data set containing the daily ice extent at Arctic Sea from 1978 to 2019,
#' collected by National Oceanic and Atmospheric Administration (NOAA)
#'
#' @format A data frame with 13391 rows and 6 variables:
#' \describe{
#'   \item{Year}{years of available data (1978--2019)}
#'   \item{Month}{month (01--12)}
#'   \item{Day}{day of the month indicated in Column Month (Fair, Good, Very Good, Premium, Ideal)}
#'   \item{Extent}{daily ice extent, to three decimal places}
#'   \item{Missing}{whether a day is missing (1) or not (0))}
#'   \item{Source}{data source in NOAA database}
#' }
#' @source <https://nsidc.org/data/G02135/versions/3>
"arctic"
