#' partition ET by TEA
#'
#' The Transpiration Estimation Algorithm (TEA) partitions ET by first
#' constraining the data to dry conditions where T/ET~1
#' (see \code{\link{tea_filter}}) and then learning water use efficiency (WUE)
#' relationship with predictors from the resulting training data
#' (see \code{link{tea_fit_wue}}).
#' The learned model is then used to predict WUE, T, and ET for each
#'
#' @param data data.frame with required columns: XX
#' @param control list with configuration options see \code{\link{tea_config}}
#'
#' @return see \code{\link{tea_predict}}
#' @export
partition_tea <- function(data,  control = tea_config()) {
  tea_checkvars(data)
  dff <- tea_filter(data, control)
  rf <- tea_fit_wue(dff, control)
  datap <- tea_predict_wue(rf, data, control)
}

#' configure hyperparameters of the TEA algorithm
#'
#' @param CSWIlimit filter for conservative surface
#'     wetness index to be lower than this limit
#' @param perc numeric vector of percentiles of the prediction distribution
#'     do average over
#' @param smax numeric scalar of maximum water storage
#'     (see \code{\link{compute_cswi}})
#' @param s0  numeric scalar of initial water storage
#'     (see \code{\link{compute_cswi}})
#' @param GPPlimit numeric scalar (mumol/m2/s): filter for half-hours
#'     with carbon fluxes larger than this half-hourly threshold
#' @param GPPdaylimit numeric scalar (mumol/m2/day): filter for days
#'     with carbon fluxes larger than this daily threshold
#' @param Tairlimit numeric scalar (degree Celsius): filter for half-hours
#'     with air temperatures larger than this threshold
#' @param Rglimit numeric scalar (W/m2): filter for half-hours
#'     with incoming radiation larger than this threshold
#'
#' @return list with the arguments
#' @export
tea_config <- function(
  CSWIlimit = -0.5
  ,perc = c(0.75)
  ,smax = 5.0
  ,GPPlimit = 0.05
  ,GPPdaylimit = 0.5
  ,Tairlimit = 5
  ,Rglimit = 0
) {
  list(
    CSWIlimit = CSWIlimit
    ,perc = perc
    ,smax = smax
    ,GPPlimit = GPPlimit
    ,Tairlimit = Tairlimit
    ,Rglimit = Rglimit
  )
}

tea_checkvars <- function(data) {
  required_num_vars <- c("ET","precip","Tair","Rg","GPP","RH","u")
  required_timestamp_vars <- c("timestamp")
  required_vars <- c(required_num_vars, required_timestamp_vars)
  imissing <- which(!(required_vars %in% colnames(data)))
  if (length(imissing)) stop(
    "Expected the following columns in data but were missing: "
    ,paste(required_vars[imissing], sep = ","))
  imissing <- c() #which(!(is.numeric(data))
  if (length(imissing)) stop(
    "Expected the following columns in data to be numeric but were not: "
    ,paste(required_num_vars[imissing], sep = ","))
  # check for half-hourly consistent time step
  if (!inherits(data$timestamp, "POSIXct")) stop(
    "timestamp needs to be of class POSIXct")
  dt <- diff(as.numeric(data$timestamp)) * 60
  ifailure <- which(dt != 30)
  if (length(ifailure)) stop(
    "need equidistant timestep of 30 minutes, but step at position ",
    ifailure[1]," was ", dt[ifailure[1]], "minutes.")
  # check for complete days (starting at 00:30 and ending at 00:00)
  # if (as.POSIXlt(data$timestep[1])$hour != 0 |
  #     as.POSIXlt(data$timestep[1])$min != 30) stop(
  #   "Expected first time step to be at 00:30, but was ", data$timestep[1])
  # nrec <- nrow(data)
  # if (as.POSIXlt(data$timestep[nrec])$hour != 0 |
  #     as.POSIXlt(data$timestep[nrec])$min != 0) stop(
  #       "Expected last time step to be at 00:00, but was ", data$timestep[nrec])
}

tea_filter <- function(data, control) {
  data$cswi <- compute_cswi(data, control$smax)
  data$GPPday <- compute_daily_GPP(data$GPP)
  dff <- data %>% filter(
    GPP > control$GPPlimit &
    GPPday > control$GPPdaylimit &
    Tair > control$Tairlimit &
    Rg > control$RgLimit &
    cswi < control$CSWIlimit
  )
}

#' Compute the conservative soil moisture index (CSWI)
#'
#' @param data data.frame with numeric columns \code{precip} and \code{ET} in mm
#' @param smax maximum water storage capacity in mm
#' @param s0 initial water storage
#'
#' @return numeric vector of CSWI in mm
#' @export
compute_cswi <- function(data, smax) {
  nrec <- nrow(data)
  s <- numeric(nrec)
  ds <- data$precip - data$ET
  s[1] <- smax #min(s0 + ds, smax)
  for (i in 2:nrec) {
    s[i] <- min(s[i-1] + ds[i], smax)
  }
  cswi <- pmax(s, pmin(data$precip, smax))
}

#' sum GPP by day of year
#'
#' @param GPP numeric vector to be summed
#' @param timestamp POSIXct vector of times
#' @param doy factor alternative to timestamp to specify the day of year
#'
#' @return numeric vector (length GPP) with corresponding daily sum of GPP
#' @export
compute_daily_GPP <- function(GPP, timestamp) {
  yr <- as.POSIXlt(timestamp)$year
  doy = as.POSIXlt(timestamp)$yday
  GPPday0 <- aggregate(GPP, list(doy = doy, yr = yr), sum)$x
  ndoyt <- table(doy, yr)
  ndoy <- as.integer(ndoyt)[1:length(GPPday0)]
  GPPday <- rep(GPPday0, times = ndoy )
}





