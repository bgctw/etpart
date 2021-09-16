#' preprocess for TEA partitioning
#'
#' Builds all the derived variables used in the partitioning such as CSWI,
#' DWCI, etc., as well as filters which remove night time, low air temp/GPP,
#' etc. periods.
#'
#' @param data data.frame of full days, must contain at least timestamp,
#'   ET, GPP, RH, Rg, Rg_pot, Tair, VPD, precip,
#' @param control
#'
#' @return ds with filters and additional indices
#' \itemize{
#'    \item \code{cswi}: see \code{\link{compute_cswi}},
#'    \item \code{C_ET}, \code{C_Rg}: diurnal centroids,
#'      see \code{\link{compute_diurnal_centroid}},
#'    \item \code{dcwi}: see \code{\link{compute_DWCI}},
#'    \item \code{Rg_pot_daily}: daily daily Rg_pot in MJ m-2 d-1
#'    \item \code{year}: the year of the time-stamp
#'       use mid-time-stamps. If using end-timestampes, NewYear is already next
#'    \item \code{GPPgrad}: daily smoothed GPP gradient in umol C m-2 s-1 d-1
#'      , which gives and indication of phenology
#'    \item \code{Rpotgrad} and \code{Rpotgrad_day}: gradients in \code{Rpot}
#'       and daily means of \code{Rpot}, which give an indication of time
#'       in the year, ie. season
#'    \item \code{tempFlag}: records with minimum temperature
#'       (see \code{\link{tea_config}})
#'    \item \code{GPPFlag}: records with minimum GPP and minimum daily GPP
#'       (see \code{\link{tea_config}})
#'    \item \code{seasonFlag}: combined \code{tempFlag} and \code{GPPFlag}
#' }
#' @export
tea_preprocess <- function(data,  control = tea_config()) {
  nrec = nrow(data)
  nday = nrec / control$units_per_day
  nrecday = control$units_per_day
  df <- data %>% mutate(
    cswi = compute_cswi(data, smax = control$smax)
    , iday = rep(1:nday, each = nrecday)
    , year = as.POSIXlt(timestamp)$year + 1900
    , C_ET = rep(compute_diurnal_centroid(ET, nrecday), each = nrecday)
    , C_Rg = rep(compute_diurnal_centroid(Rg, nrecday), each = nrecday)
    #, dcwi = compute_DWCI(data, nrecday)
    , GPPgrad = compute_GPPgrad(GPP, nrecday)
    , Rpotgrad = ETPart:::gradient_equi(Rg_pot)
    , Rpotgrad_day = rep(ETPart:::gradient_equi(
        colMeans(matrix(Rg_pot, nrow = nrecday))), each = nrecday)
    , DayNightFlag = Rg_pot > 0
    , posFlag = GPP > 0 & ET > 0
    , tempFlag = Tair > control$tempdaymin
    , GPPday = rep(colMeans(matrix(GPP, nrow = nrecday)), each = nrecday)
    , GPPFlag = GPPday > control$GPPdaymin & GPP > control$GPPmin
    , seasonFlag = tempFlag * GPPFlag
  )
  dfday <- df %>% group_by(iday) %>%
    summarise(
      Rg_pot_daily = sum(Rg_pot) *((3600*(24/nrecday))/1000000)
      ,.groups = "drop")
  df <- df %>% mutate(
    Rg_pot_daily = rep(dfday$Rg_pot_daily, each = nrecday)
  )
}

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
#' @param GPPlimit numeric scalar (mumol/m2/s): filter for half-hours
#'     with carbon fluxes larger than this half-hourly threshold
#' @param GPPdaylimit numeric scalar (mumol/m2/day): filter for days
#'     with carbon fluxes larger than this daily threshold
#' @param Tairlimit numeric scalar (degree Celsius): filter for half-hours
#'     with air temperatures larger than this threshold
#' @param Rglimit numeric scalar (W/m2): filter for half-hours
#'     with incoming radiation larger than this threshold
#' @param units_per_day number of records per day, 48 for half-hourly data
#' @param tempdaymin flag records of days having at list this minimum
#'     daily air temperature in degree Celsius
#' @param GPPdaymin flag records of days having at least this daily GPP in
#'     gC m-1 d-1
#' @param GPPmin flag days having at least this GPP in umol m-2 s-1
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
  ,units_per_day = 48
  ,tempdaymin = 5
  ,GPPdaymin=0.5
  ,GPPmin=0.05
) {
  list(
    CSWIlimit = CSWIlimit
    ,perc = perc
    ,smax = smax
    ,GPPlimit = GPPlimit
    ,Tairlimit = Tairlimit
    ,Rglimit = Rglimit
    ,units_per_day = units_per_day
    ,tempdaymin = tempdaymin
    ,GPPdaymin = GPPdaymin
    ,GPPmin = GPPmin
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
compute_cswi <- function(data, smax = 5) {
  nrec <- nrow(data)
  s <- numeric(nrec)
  # conservative: not knowing precip - refill upper soil storage
  precip_f <- data$precip
  precip_f[!is.finite(data$precip)] <- smax
  ET_f <- data$ET
  ET_f[!is.finite(data$ET)] <- -9999 # to comply to Nelson18
  ds <- precip_f - ET_f
  s[1] <- smax #min(s0 + ds, smax)
  for (i in 2:nrec) {
    s[i] <- min(s[i-1] + ds[i], smax)
  }
  cswi <- pmax(s, pmin(precip_f, smax))
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

#' Normalized diurnal centroid of latent energy (LE)
#'
#' Calculates the diurnal centroid of LE relative to the diurnal centroid of
#' incoming radiation (Rg).
#'
#' @param flux numeric vector of sub-daily flux (usuall LE)
#'   that must be condinuous and regular
#' @param Rg incoming radiation, can be any unit
#' @param units_per_day integer: frequency of the sub-daily measurements,
#'    48 for half hourly measurements
#'
#' @return The normalized diurnal centroid, in hours,
#'  at a daily frequency
#' @export
compute_norm_diurnal_centroid <- function(flux, Rg, units_per_day = 48){
  C_LE = compute_diurnal_centroid(flux, units_per_day=units_per_day)
  C_Rg = compute_diurnal_centroid(Rg, units_per_day=units_per_day)
  C_LE - C_Rg
}

#' Diurnal centroid of sub-daily fluxes
#'
#' Calculates the daily flux weighted time of a sub-daily flux.
#'
#' @param flux numeric vector of sub-daily flux that must be continuous
#'    and regular of full days
#' @param units_per_day integer: frequency of the sub-daily measurements,
#'    48 for half hourly measurements
#'
#' @return The diurnal centroid, in hours at a daily frequency
#' @export
compute_diurnal_centroid <- function(flux, units_per_day = 48) {
  nday = length(flux) / units_per_day
  fluxm <- matrix(flux, nrow = units_per_day)
  hour = (1:units_per_day)/units_per_day*24
  dci <- apply(fluxm, 2, function(x){
    sum(x*hour)/sum(x)
  })
  # calculate the total number of days
  # days,UPD=flux.reshape(-1,UnitsPerDay).shape
  # # create a 2D matrix providing a UPD time series for each day, used in the
  # matrix operations.
  # hours=np.tile(np.arange(UPD),days).reshape(days,UPD)
  # # calculate the diurnal centroid
  # C=np.sum(hours*flux.reshape(-1,48),axis=1)/np.sum(flux.reshape(-1,48),axis=1)
  # C=C*(24/UnitsPerDay)
  # return(C)
}


#' Diurnal water:carbon index (DWCI)
#'
#' DWCI measures the probability that the carbon and water are coupled
#' within a given day. Method takes the correlation between
#' evapotranspiration (LE) and gross primary productivity (GPP) and
#' calculates the correlation within each day. This correlation is then
#' compared to a distribution of correlations between artificial datasets
#' built from the signal of potential radiation and the uncertainty in
#' the LE and GPP.
#'
#' @param data data.frame of sub-daily timeseries with variables
#' \describe{
#' \item{Rg_pot}{Potential radiation}
#' \item{ET}{evapotranspiration or latent energy}
#' \item{GPP}{ross primary productivity}
#' \item{VPD}{vapor pressure deficit}
#' \item{NEE}{net ecosystem exchange}
#' \item{ET_sd}{estimation of the uncertainty of ET}
#' \item{GPP_sd}{eestimation of the uncertainty of GPP}
#' \item{NEE_fall}{ Modeled net ecosystem exchange i.e. no noise}
#' \item{ET_fall}{Modeled evapotranspiration or latent energy i.e. no noise}
#' }
#' @param units_per_day integer: frequency of the sub-daily measurements,
#'    48 for half hourly measurements'
#'
#' @return The diurnal water:carbon index (DWCI)
#' @export
compute_DWCI <- function(data, units_per_day = 48) {
  nrep = 100 # the number of artificial datasets to construct
  nday = nrow(data) / units_per_day
  #   # creates an empty 2D dataset to hold the artificial distributions
  #   StN=np.zeros([repeats,days])*np.nan
  #   corrDev=np.zeros([days,2,2])
  #
  dfg <- data %>%
    mutate(iday = rep(1:nday, each = units_per_day)) %>%
    group_by(iday) %>%
    mutate(
      #   # create the daily cycle by dividing Rg_pot by the daily mean
      #   daily_cycle=Rg_pot/Rg_pot.mean(axis=1)[:,None]
      daily_cycle = Rg_pot/mean(Rg_pot, na.rm = FALSE),
      # Isolate the error of the carbon and water fluxes.
      NEE_err = NEE_fall - NEE,
      ET_err = ET_fall - ET
    ) %>%
    nest()
  dfday <- dfg %>%
    mutate(df_corr_syn = map(data, correlations_from_artificial)) %>%
    unnest(cols = c(iday, df_corr_syn))
  dfday$dwci
}

#' Compute Artificial correlations from daily data
#'
#' Correlation between ET and GPP is obscured by noise. The effect size
#' depends on signal (flux) to noise (errors) ratio.
#' This function bootstaps the correlation of two perfectly correlated signals
#' with random noise added corresponding to errors in ET and NEE.
#' and then gives the rank
#'
#' @param dfs tibble with columns GPP, ET, GPP_sd, ET_sd, NEE_err, ET_err
#'
#' @return data.frame by day with columns mean_GPP, mean_ET,
#'   corr_err (pearson correlation coefficient between errors of ET and NEE),
#'   corr_syn (vector of perason correlation coefficients between GPP and ET)
#'   dwci (rank of observed correlation within artificial correlations in %)
correlations_from_artificial <- function(dfs, nboot = 100){
  mean_GPP = mean(dfs$GPP, na.rm = TRUE)
  mean_ET = mean(dfs$ET, na.rm = TRUE)
  daily_cylce = dfs$Rg_pot / mean(dfs$Rg_pot, na.rm = TRUE)
  corr_err = ifelse(
    !is.finite(mean_GPP) | !is.finite(mean_ET) |
      any(is.na(dfs$NEE_err)) | any(is.na(dfs$ET_err)) |
      any(is.na(dfs$GPP_sd)) | any(is.na(dfs$ET_sd))
    #,NA_real_,
    , 0.0 , # assume no correlation in errors
    ifelse(all(dfs$NEE_err == 0) | all(dfs$ET_err == 0),
           1.0, cor(-dfs$NEE_err, dfs$ET_err)))
  # generate artificial GPP and ET daily vectors
  GPP_syn <- ET_syn <- matrix(NA_real_, nrow = nrow(dfs), ncol = nboot)
  if (corr_err == 0) {
    for (irow in 1:nrow(dfs)) {
      GPP_syn[irow,] <- daily_cylce[irow] * mean_GPP +
        rnorm(nboot, dfs$GPP_sd[irow])
      ET_syn[irow,] <- daily_cylce[irow] * mean_ET +
        rnorm(nboot, dfs$ET_sd[irow])
    }
  } else {
    corm = matrix(c(1, corr_err, corr_err, 1), nrow = 2)
    for (irow in 1:nrow(dfs)) {
      sdm <- diag(c(dfs$GPP_sd[irow], dfs$ET_sd[irow]))
      covm <- sdm %*% corm %*% sdm
      binorm <- rbinorm(nboot, covm)
      GPP_syn[irow,] <- daily_cylce[irow] * mean_GPP + binorm[,1]
      ET_syn[irow,] <- daily_cylce[irow] * mean_ET + binorm[,2]
    }
  }
  # nboot artificial correlations - without VPD effect because R_pot is was used
  corr_syn = daily_corr(
    as.vector(ET_syn), as.vector(GPP_syn), rep(dfs$Rg_pot, nboot), nrow(dfs))
  # real correlation
  corr_obs <- daily_corr(dfs$ET, dfs$GPP*sqrt(dfs$VPD), dfs$Rg_pot, nrow(dfs))
  # rank
  dwci = sum(corr_syn < corr_obs)* 100/nboot
  nrow = tibble(
    mean_GPP = mean_GPP,
    mean_ET = mean_ET,
    corr_err = corr_err,
    dwci = dwci
    #,corr_art = corr_art
  )
}

rbinorm <- function(N, sigma) {
  # https://blog.revolutionanalytics.com/2016/08/simulating-form-the-bivariate-normal-distribution-in-r-1.html
  M <- t(chol(sigma))
  # M %*% t(M)
  Z <- matrix(rnorm(2*N),2,N) # 2 rows, N/2 columns
  bvn2 <- t(M %*% Z) #+ matrix(rep(mu,N), byrow=TRUE,ncol=2)
}

#' Daily correlation coefficient
#'
#' Calculates a daily correlation coefficient between two sub-daily timeseries
#' \code{x} and \code{y} considering only daytime, i.e. \code{Rg > 0}.
#' Vectors \code{x}, \code{y}, and \code{Rg} must have the same length.
#' Only complete cases, i.e. no NAs are considered.
#'
#' @param x
#' @param y
#' @param Rg_pot potential incoming radiation
#' @param units_per_day integer: frequency of the sub-daily measurements,
#'    48 for half hourly measurements
#' @return correlation coefficents at daily timescale
daily_corr <- function(x, y, Rg_pot, units_per_day = 48) {
  nday = length(x) / units_per_day
  dff <- data.frame(
    x = x, y = y, Rg_pot = Rg_pot,
    iday = rep(1:nday, each = units_per_day)
  )
  # corf <- function(x,y) {
  #   # pearson product-moment correlation coefficient without na checking
  #   mean( (x-mean(x))*(y-mean(y)) ) / (sd(x)*sd(y))
  # }
  ans <- dff %>%
    filter(complete.cases(.)) %>%
    filter(Rg_pot > 0) %>%
    group_by(iday) %>%
    summarise(
      r2 = cor(x,y)^2
    )
  ans$r2
  # x=x.reshape(-1,48)
  # y=y.reshape(-1,48)
  # Rg_pot=Rg_pot.reshape(-1,48)
  # mask=Rg_pot<=0
  # x=np.ma.MaskedArray(x,mask=mask)
  # y=np.ma.MaskedArray(y,mask=mask)
  # x=x/x.max(axis=1)[:,None]
  # y=y/y.max(axis=1)[:,None]
  # mx = x.mean(axis=1)
  # my = y.mean(axis=1)
  # xm, ym = x - mx[..., None], y - my[..., None]
  # r_num = np.ma.add.reduce(xm * ym, axis=1)
  # r_den = np.ma.sqrt(np.ma.sum(xm**2, axis=1) * np.ma.sum(ym**2, axis=1))
  # r = r_num / r_den
  # return(r**2)
}

#' Compute the seasonal gradiant from smoothed daily GPP.
#'
#' @param GPP numeric vector of gross primary productivity (umol m-2 s-1)
#' @param nStepsPerDay number of records per day
#'
#' @return numeric vector of gradients
#' @export
compute_GPPgrad <- function(GPP, units_per_day = 48, sigma = 20.0){
  gradGPP=GPP
  gradGPP[is.na(GPP)] <- 0
  gradGPP[(GPP < -9000)] <- 0
  #GPPgrad <- rep=np.repeat(np.gradient(gaussian_filter(gradGPP.reshape(-1,nStepsPerDay).mean(axis=1),sigma=[20])),nStepsPerDay)
  gppDay <- colMeans(matrix(GPP, nrow = units_per_day))
  gppDay_smooth <- gaussianSmooth(gppDay, sigma) # from mmand
  gppDay_grad <- ETPart:::gradient_equi(gppDay_smooth)
  GPPgrad <- rep(gppDay_grad, each = units_per_day)
  GPPgrad[(GPP < -9000)] <- NA
  GPPgrad[is.na(GPP)] <- NA
  GPPgrad[0]=0
  return(GPPgrad)
}

gradient_equi <- function(x){
  # https://numpy.org/doc/stable/reference/generated/numpy.gradient.html
  # second order central differences in the interior points and
  # first order at boundaries
  # for equidistant series
  # (x[i+1] - x[i-1])/2
  gr_inner <- (x[-(1:2)] - x[-(length(x)-(1:0))])/2
  # single gradients at the ends
  gr <- c(diff(x[1:2]), gr_inner, diff(x[length(x)-(1:0)]))
}





