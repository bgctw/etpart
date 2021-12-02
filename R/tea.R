#' preprocess for TEA partitioning
#'
#' Builds all the derived variables used in the partitioning such as CSWI,
#' DWCI, etc., as well as filters which remove night time, low air temp/GPP,
#' etc. periods.
#'
#' @param data data.frame of full days, must contain at least timestamp,
#'   ET, GPP, RH, Rg, Rg_pot, Tair, VPD, precip,
#' @param control see \code{\link{tea_config}}
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
#'    \item \code{Rgpotgrad} and \code{Rpotgrad_day}: gradients in \code{Rpot}
#'       and daily means of \code{Rpot}, which give an indication of time
#'       in the year, ie. season
#'    \item \code{tempFlag}: records with minimum temperature
#'       (see \code{\link{tea_config}})
#'    \item \code{GPPFlag}: records with minimum GPP and minimum daily GPP
#'       (see \code{\link{tea_config}})
#'    \item \code{seasonFlag}: combined \code{tempFlag} and \code{GPPFlag}
#'    \item \code{inst_WUE}: instant water use efficiency (GPP/ET)
#'        in g C per kg H2O
#' }
#' @export
tea_preprocess <- function(data,  control = tea_config()) {
  nrec = nrow(data)
  nday = nrec / control$nrecday
  nrecday = control$nrecday
  df <- data %>% mutate(
    CSWI = compute_cswi(data, smax = control$smax)
    , iday = rep(1:nday, each = nrecday)
    , year = as.POSIXlt(.data$timestamp)$year + 1900
    , C_ET = rep(compute_diurnal_centroid(.data$ET, nrecday), each = nrecday)
    , C_Rg = rep(compute_diurnal_centroid(.data$Rg, nrecday), each = nrecday)
    , C_Rg_ET = .data$C_ET - .data$C_Rg
    , dcwi =
    , GPPgrad = compute_GPPgrad(.data$GPP, nrecday)
    , Rgpotgrad = gradient_equi(.data$Rg_pot)
    , Rpotgrad_day = rep_daily_aggregate(
        .data$Rg_pot, compose(gradient_equi, colMeans), nrecday)
    , DayNightFlag = .data$Rg_pot > 0
    , posFlag = .data$GPP > 0 & .data$ET > 0
    , tempFlag = .data$Tair > control$tempdaymin
    , GPPday = rep_daily_aggregate(.data$GPP, colMeans, nrecday)
    , GPPFlag = .data$GPPday > control$GPPdaymin & .data$GPP > control$GPPmin
    , seasonFlag = .data$tempFlag * .data$GPPFlag
    , inst_WUE = ifelse(.data$ET <= 0, 0, .data$GPP/.data$ET) * (12*1800)/1000
  )
  df <- if (all(c("NEE_fall","ET_fall") %in% names(data))) {
      df %>% mutate(DWCI = rep(compute_DWCI(data, nrecday), each = nrecday))
  } else {
    warning(
      "Missing columns NEE_fall or ET_fall, therefore computing simplified ",
      "DWCI with neglecting daily correlation between NEE and ET errors.")
    df %>% mutate(
      DWCI = rep(compute_simplifiedDWCI(data, nrecday), each = nrecday))
  }
  if (!("quality_flag" %in% colnames(data))) df$quality_flag <- TRUE
  dfday <- df %>% group_by(.data$iday) %>%
    summarise(
      Rg_pot_daily = sum(.data$Rg_pot) *((3600*(24/nrecday))/1000000)
      ,.groups = "drop")
  df <- df %>% mutate(
    Rg_pot_daily = rep(dfday$Rg_pot_daily, each = nrecday)
  )
}

#' Compute daily aggregate and repeat it for each record
#'
#' x must have the same number of records for each day and contain
#' only complete days.
#'
#' @param x vector
#' @param FUN \code{function(x) -> numeric scalar} to aggregate over
#'  each column of a matrix with a day in each row, i.e. \code{colSums} or
#'  \code{applycols(sum)}
#'   daily numeric vector
#' @param nrecday number of records fore each day
#'
#' @return vector of the same length as x
rep_daily_aggregate <- function(x, FUN, nrecday) {
  rep(FUN(matrix(x, nrow = nrecday)), each = nrecday)
}

applycols <- function(FUN){
  function(x) apply(x, 2, FUN)
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
#' @return see \code{\link{tea_predict}}, \code{data} with predictions
#'   percentiles of WUE, E and T appended.
#' @export
partition_tea <- function(data, control = tea_config()) {
  tea_checkvars(data)
  data_train <- tea_filter(data, control)
  rf <- tea_fit_wue(data_train, control)
  datap <- tea_predict(rf, data, control)
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
#' @param nrecday number of records per day, 48 for half-hourly data
#' @param tempdaymin flag records of days having at list this minimum
#'     daily air temperature in degree Celsius
#' @param GPPdaymin flag records of days having at least this daily GPP in
#'     gC m-1 d-1
#' @param GPPmin flag days having at least this GPP in umol m-2 s-1
#' @param rfseed random generator seed used for random-forest fit
#' @param quantiles_wue quantiles for which precitions of WUE
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
  ,nrecday = 48
  ,tempdaymin = 5
  ,GPPdaymin=0.5
  ,GPPmin=0.05
  ,rfseed = 63233
  ,quantiles_wue = seq(0.5,1.0,length.out = 11)
) {
  list(
    CSWIlimit = CSWIlimit
    ,perc = perc
    ,smax = smax
    ,GPPlimit = GPPlimit
    ,GPPdaylimit = GPPdaylimit
    ,Tairlimit = Tairlimit
    ,Rglimit = Rglimit
    ,nrecday = nrecday
    ,tempdaymin = tempdaymin
    ,GPPdaymin = GPPdaymin
    ,GPPmin = GPPmin
    ,rfseed = rfseed
    ,quantiles_wue = quantiles_wue
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

#' Filter conditions for training the TEA WUE model.
#'
#' @param data data.frame with columns GPP, Tair, Rg, and cswi
#' @param control list with entries GPPlimit, GPPdaylimit, Tairlimit, Rglimit,
#'   and CWSIlimit. See \code{\link{tea_config}}.
#'
#' @return data.frame containing only rows matching the filter criteria.
#' @export
tea_filter <- function(data, control) {
  #data$cswi <- compute_cswi(data, control$smax)
  data$GPPday <- compute_daily_GPP(data$GPP, data$timestamp)
  dff <- data %>% filter(
    ,.data$GPP > control$GPPlimit &
      .data$GPPday > control$GPPdaylimit &
      .data$Tair > control$Tairlimit &
      .data$Rg > control$Rglimit &
      .data$CSWI < control$CSWIlimit
  )
}

#' Compute the conservative soil moisture index (CSWI)
#'
#' @param data data.frame with numeric columns \code{precip} and \code{ET} in mm
#' @param smax maximum water storage capacity in mm
#'
#' @return numeric vector of CSWI in mm
#' @export
compute_cswi <- function(data, smax = 5) {
  nrec <- nrow(data)
  # conservative: not knowing precip - refill upper soil storage
  precip_f <- data$precip
  precip_f[!is.finite(data$precip)] <- smax
  precip_f[data$precip < 0] <- smax
  ET_f <- data$ET
  ET_f[!is.finite(data$ET)] <- -9999 # to comply to Nelson18 ?
  ds <- precip_f - ET_f
  s <- numeric(nrec)
  s[1] <- smax #min(s0 + ds, smax)
  # implementatation according to Nelson18 paper does not give negative values?
  # for (i in 2:nrec) {
  #   s[i] <- min(s[i-1] + ds[i], smax)
  # }
  # cswi <- pmax(s, pmin(precip_f, smax)) # wrong: only if precip > 0
  # implementation according to TEA/CWSI.py
  for (i in 2:nrec) {
    # bound current water balance by smax
    stepval <- min(s[i-1] + ds[i], smax)
    # in case of a positive precip value, the current CSWI is the max between the previous
    # CSWI and either the value of the precip or the s0 depending on which is smaller
    s[i] <- if (precip_f[i] > 0) {
      max(stepval, min(precip_f[i-1], smax))  # i-1?
    } else {
      # if there is no precip, the CSWI is according to the stepVal,
      # causing simple water balance behaviour
      stepval
    }
  }
  s
}

#' sum GPP by day of year
#'
#' @param GPP numeric vector to be summed
#' @param timestamp POSIXct vector of times
#'
#' @return numeric vector (length GPP) with corresponding daily sum of GPP
#' @export
compute_daily_GPP <- function(GPP, timestamp) {
  yr <- as.POSIXlt(timestamp)$year
  doy = as.POSIXlt(timestamp)$yday
  GPPday0 <- aggregate(GPP, list(doy = doy, yr = yr), sum)$x
  ndoyt <- table(yr, doy)
  # skip 0-entries at start and end
  ndoy <- ndoyt[ndoyt > 0]
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
#' @param nrecday integer: frequency of the sub-daily measurements,
#'    48 for half hourly measurements
#'
#' @return The normalized diurnal centroid, in hours,
#'  at a daily frequency
#' @export
compute_norm_diurnal_centroid <- function(flux, Rg, nrecday = 48){
  C_LE = compute_diurnal_centroid(flux, nrecday=nrecday)
  C_Rg = compute_diurnal_centroid(Rg, nrecday=nrecday)
  C_LE - C_Rg
}

#' Diurnal centroid of sub-daily fluxes
#'
#' Calculates the daily flux weighted time of a sub-daily flux.
#'
#' @param flux numeric vector of sub-daily flux that must be continuous
#'    and regular of full days
#' @param nrecday integer: frequency of the sub-daily measurements,
#'    48 for half hourly measurements
#'
#' @return The diurnal centroid, in hours at a daily frequency
#' @export
compute_diurnal_centroid <- function(flux, nrecday = 48) {
  nday = length(flux) / nrecday
  fluxm <- matrix(flux, nrow = nrecday)
  hour = (1:nrecday)/nrecday*24
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

#' simplified Diurnal water:carbon index (DWCI)
#'
#' similar to DWCI it
#' measures the probability that the carbon and water are coupled
#' within a given day. In difference to the full DWCI the correlation
#' only the standard deviation of ET and GPP is used but not the
#' daily correlation between NEE and ET errors.
#' Hence, it can be used when noisefree, i.e. modelled, NEE_fall and ET_fall
#' are not available.
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
#' }
#' @param nrecday integer: frequency of the sub-daily measurements,
#'    48 for half hourly measurements'
#' @param na_value numeric scalar: to replace NA values in correlation
#'    Defaults to 0.0 - i.e. no probability of correlation in DWCI.
#'
#' @return numeric vector length(data)/nrecday: diurnal water:carbon index (DWCI)
#' @export
compute_simplifiedDWCI <- function(data, nrecday = 48, na_value = 0) {
  nboot = 100 # the number of artificial datasets to construct
  nday = nrow(data) / nrecday
  #   # creates an empty 2D dataset to hold the artificial distributions
  #   StN=np.zeros([repeats,days])*np.nan
  #   corrDev=np.zeros([days,2,2])
  #
  dfg <- data %>%
    mutate(iday = rep(1:nday, each = nrecday)) %>%
    group_by(.data$iday) %>%
    mutate(
      #   # create the daily cycle by dividing Rg_pot by the daily mean
      #   daily_cycle=Rg_pot/Rg_pot.mean(axis=1)[:,None]
      daily_cycle = .data$Rg_pot/mean(.data$Rg_pot, na.rm = FALSE),
      # note that dfg is grouped by iday, so mean is across records of a day
      GPP_mean = mean(.data$GPP),
      ET_mean = mean(.data$ET)
    )
  sd_GPP <- dfg %>%
    summarise(sd_GPP = sd(.data$GPP * dfg$GPP*sqrt(dfg$VPD))) %>% chuck("sd_GPP")
  # bootstrap daily correlation by nrep, each column is a sample of all days
  corr_syn <- do.call(cbind, map(1:nboot, function(irep){
    dfg <- dfg %>% mutate(
      #NEE_err = suppressWarnings(rnorm(nrecday, sd = .data$NEE_sd)),
      NEE_err = suppressWarnings(rnorm(nrecday, sd = .data$GPP_sd)),
      ET_err = suppressWarnings(rnorm(nrecday, sd = .data$ET_sd)),
      GPP_DayCycle = .data$daily_cycle * .data$GPP_mean + .data$NEE_err,
      ET_DayCycle = .data$daily_cycle * .data$ET_mean + .data$ET_err
    )
    dcorr = daily_corr(
      dfg$ET_DayCycle, dfg$GPP_DayCycle*sqrt(dfg$VPD),
      Rg_pot = dfg$Rg_pot, nrecday = nrecday)
  }))
  corr_obs <- daily_corr(dfg$ET, dfg$GPP*sqrt(dfg$VPD), dfg$Rg_pot, nrecday)
  # rank
  dwci <- rowSums(corr_syn < corr_obs)* 100/nboot # use vector recycling of corr_obs
  dwci[sd_GPP == 0] <- 0 # prob of coupling is zero fi sd_GPP == 0 instead of NA
  dwci[is.na(dwci)] <- na_value
  dwci
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
#' @param nrecday integer: frequency of the sub-daily measurements,
#'    48 for half hourly measurements'
#' @param na_value numeric scalar: to replace NA values in correlation
#'    Defaults to 0.0 - i.e. no probability of correlation in DWCI.
#'
#' @return numeric vector length(data)/nrecday: diurnal water:carbon index (DWCI)
#' @export
compute_DWCI <- function(data, nrecday = 48, na_value = 0.0) {
  nrep = 100 # the number of artificial datasets to construct
  nday = nrow(data) / nrecday
  #   # creates an empty 2D dataset to hold the artificial distributions
  #   StN=np.zeros([repeats,days])*np.nan
  #   corrDev=np.zeros([days,2,2])
  #
  dfg <- data %>%
    mutate(iday = rep(1:nday, each = nrecday)) %>%
    group_by(.data$iday) %>%
    mutate(
      #   # create the daily cycle by dividing Rg_pot by the daily mean
      #   daily_cycle=Rg_pot/Rg_pot.mean(axis=1)[:,None]
      daily_cycle = .data$Rg_pot/mean(.data$Rg_pot, na.rm = FALSE),
      # Isolate the error of the carbon and water fluxes.
      NEE_err = .data$NEE_fall - .data$NEE,
      ET_err = .data$ET_fall - .data$ET
    ) %>%
    nest()
  dfday <- dfg %>%
    mutate(df_corr_syn = map(data, correlations_from_artificial)) %>%
    unnest(cols = c(.data$iday, .data$df_corr_syn))
  dwci <- dfday$dwci
  dwci[data$sd_GPP == 0] <- 0 # prob of coupling is zero fi sd_GPP == 0 instead of NA
  dwci[is.na(dwci)] <- na_value
  dwci
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
#' @param nboot number of bootstrap samples used
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
  dwci = sum(corr_syn < corr_obs)* 100/nboot # use vector recycling of corr_obs
  ans = tibble(
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
#' @param x numeric vector to correlate
#' @param y numeric vector to correlate
#' @param Rg_pot potential incoming radiation
#' @param nrecday integer: frequency of the sub-daily measurements,
#'    48 for half hourly measurements
#' @return correlation coefficents at daily timescale
daily_corr <- function(x, y, Rg_pot, nrecday = 48) {
  nday = length(x) / nrecday
  dff <- data.frame(
    x = x, y = y, Rg_pot = Rg_pot,
    iday = rep(1:nday, each = nrecday)
  )
  # corf <- function(x,y) {
  #   # pearson product-moment correlation coefficient without na checking
  #   mean( (x-mean(x))*(y-mean(y)) ) / (sd(x)*sd(y))
  # }
  ans <- dff %>%
    #drop_na() %>%  # can drop entire days
    #filter(.data$Rg_pot > 0) %>% # will fail at polar night - drops all days
    # set to NA - witch will  give an NA correlation if all are NA
    mutate(x = ifelse(.data$Rg_pot > 0, .data$x, NA)) %>%
    group_by(.data$iday) %>%
    summarise(
      # warning on sd(x)==0 or sd(y)==0
      r2 = suppressWarnings(cor(.data$x,.data$y, use = "na.or.complete")^2)
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
#' @param nrecday integer scalar: number of records within one day
#' @param sigma parameter to \code{mmand::gaussianSmooth} in units
#'   number of days
#' @param GPP numeric vector of gross primary productivity (umol m-2 s-1)
#'
#' @return numeric vector of gradients, repeated nrecday to match length of GPP
#' @export
compute_GPPgrad <- function(GPP, nrecday = 48, sigma = 20.0){
  gradGPP=GPP
  gradGPP[is.na(GPP)] <- 0
  gradGPP[(GPP < -9000)] <- 0
  #GPPgrad <- rep=np.repeat(np.gradient(gaussian_filter(gradGPP.reshape(-1,nStepsPerDay).mean(axis=1),sigma=[20])),nStepsPerDay)
  gppDay <- colMeans(matrix(GPP, nrow = nrecday))
  gppDay_smooth <- gaussianSmooth(gppDay, sigma) # from mmand
  gppDay_grad <- gradient_equi(gppDay_smooth)
  GPPgrad <- rep(gppDay_grad, each = nrecday)
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


#' Fit a random forest model for predicting water use effciency (WUE)
#'
#' @param data_train trainign dataset where T = ET and hence WUE = GPP/ET
#' @param control list with entries rfseed. See \code{\link{tea_config}}.
#'
#' @return a tidymodels workflow object
#' @export
tea_fit_wue <- function(data_train, control) {
  data_train_wue <- data_train %>%
    mutate(TEA_WUE = .data$GPP/.data$ET) # assuming T = ET)
  rf_recipe <-
    recipe(
      TEA_WUE ~ Rg + Tair + RH + u + Rg_pot_daily + Rgpotgrad + year + GPPgrad +
        DWCI +
        C_Rg_ET + CSWI
      , data = data_train_wue
    )
    # ) %>%
    # step_log(Sale_Price, base = 10) %>%
    # step_other(Neighborhood, Overall_Qual, threshold = 50) %>%
    # step_novel(Neighborhood, Overall_Qual) %>%
    # step_dummy(Neighborhood, Overall_Qual)
  rf_mod <- rand_forest(trees = 100, mtry = round(11/3), min_n = 1) %>%
    set_engine("ranger", importance = "impurity", seed = control$rfseed,
               quantreg = TRUE) %>%
    set_mode("regression")
  set.seed(control$rfseed)
  rf_wf <- workflows::workflow() %>%
    add_model(rf_mod) %>%
    add_recipe(rf_recipe) %>%
    fit(data_train_wue)
}

#' Predict water use efficiency WUE, transpiration T, and evaporation E
#'
#' @param rf tidymodels workflow
#' @param data data.frame with predictors as in traning data
#' @param control see \code{\link{tea_config}}
#'
#' @return \code{data} with appended columns: \code{WUE_perc}, \code{T_perc},
#'   and \code{ET_perc}
#'   where \code{perc} is each prediction percentile corresponding to the
#'   \code{control$quantile_wue}, e.g. 75 for the 0.75 prediction quantile.
#' @export
tea_predict <- function(rf, data, control) {
  wue_pred <- pred_ranger_quantiles(rf, data, control$quantiles_wue)
  ET_pred <- compute_TandE(data$ET, wue_pred)
  ET_pred_wide <- ET_pred %>%
    select(-.data$ET) %>%
    pivot_wider(.data$id, .data$perc, values_from = !c(.data$id, .data$perc)) %>%
    select(-.data$id)
  ans <- replace_columns(data, ET_pred_wide)
}

#' Predict quantiles from workflow's ranger fitting object for new data
#'
#' The workflow need to be constructed with
#' \code{set_engine("ranger", quantreg = TRUE, ...)}.
#' The newdata is transformed/prepared by \code{\link{bake}}.
#'
#' @param rf tidymodels workflow
#' @param newdata data.frame with predictors as in traning data
#' @param quantiles numeric vector of probabilities of prediction distribution
#'
#' @return data.frame with one column for each quantile
pred_ranger_quantiles <- function(rf, newdata, quantiles = c(0.5,0.75)) {
  predict(
    rf$fit$fit$fit,
    workflows::extract_recipe(rf) %>% bake(newdata),
    type = "quantiles",
    quantiles = quantiles
  ) %>%
    chuck("predictions") %>% # pick element named predictions
    as_tibble() %>%
    set_names(quantiles)
}

#' compute T and E from given ET and quantiles of water use efficiency (WUE)
#'
#' @param ET numeric vector of evapotranspiration
#' @param wue_pred tibble with each column a prediction quantile of WUE
#'   specifying the quantile as a number (0.0..1.0) in the column name
#'
#' @return tibble in long format with columsn ET, perc (1..100), WUE, T, E
#'   column id identifies the original row-number in the wide format
#' @export
compute_TandE <- function(ET, wue_pred) {
  wue_pred %>%
    mutate(id = 1:n()) %>%
    bind_cols(ET = ET) %>%
    pivot_longer(-c(.data$ET, .data$id), "perc", values_to = "WUE",
                 names_transform = list(perc = ~round(as.numeric(.x)*100))) %>%
    mutate(
      T = .data$ET * .data$WUE,
      E = .data$ET - .data$T
    )
}






