.tmp.f <- function(){
  library(testthat)
  library(dplyr)
  library(purrr)
}

# regression tests against commit 6dbac52 of
# https://github.com/oscarperezpriego/ETpartitioningTutorial
# modified calculate_Chi_o and calculate_WUE_o to not deviding VPD in mean
# again by 10

test_that("calculate_longterm_leaf with wrong column names",{
  expect_error(
    ans <- calculate_longterm_leaf( rename(FIHyy, GPP_bla = 'GPP')), "GPP")
})

test_that("calculate_longterm_leaf",{
  lt  <- calculate_longterm_leaf(FIHyy, C = 1.189, altitude = 270)
  # regression to values from ETpartitioningTutorial
  expect_equal( lt$chi_o, 0.69, tolerance = 0.01 )
  expect_equal( lt$WUE_o, 9.03, tolerance = 0.01 )
})

test_that("compute_stomatal_conductance_jarvis",{
  tz = attr(FIHyy$timestamp, "tzone")
  dsf <-  filter(FIHyy, between(
    timestamp, ISOdatetime(2010,5,15,0,0,0,tz=tz), ISOdatetime(2010,5,16,0,0,0,tz=tz)))
  ans <- ETPart:::compute_stomatal_conductance_jarvis(
    par=c(200, 0.2, 25), Q=dsf$Q, VPD=dsf$VPD, Tair=dsf$Tair, gcmax = 1)
  # answer from previous run, might be wrong
  expect_equal( mean(ans), 0.204, tolerance = 0.01 )
})

test_that("cost_optim_opriego",{
  tz = attr(FIHyy$timestamp, "tzone")
  dsf <-  filter(FIHyy, between(
    timestamp, ISOdatetime(2010,5,15,0,0,0,tz=tz), ISOdatetime(2010,5,16,0,0,0,tz=tz)))
  dsf <- filter(dsf, GPP>0 & Q>0 & Rg>0)
  constants = etpart_constants()
  ra <- ETPart:::compute_aerodynamic_conductance(dsf$u, dsf$ustar)
  dens = ETPart:::calculate_air_density(dsf$Pair, dsf$Tair, constants) #[kg m-3]
  Mden = dens/constants$M ##<< molar air density [mol m-3].
  GPP_max = quantile(dsf$GPP, probs=c(0.90), na.rm=T)
  Dmax = mean(dsf$VPD[dsf$GPP>GPP_max]/10, na.rm=T)
  GPP_sd_threshold = pmax(dsf$GPP*0.1, dsf$GPP_sd)
  cfg = priego_config()
  #cfg = priego_config(wWUE = 0.1)
  ans <- ETPart:::cost_optim_opriego(
    par=c(200, 0.2, 25, 0.6),
    chi_o = 0.88, WUE_o= 22.5
    ,GPP = dsf$GPP
    ,GPP_sd = dsf$GPP_sd
    ,H=dsf$H
    ,VPD=dsf$VPD/10
    ,Tair=dsf$Tair
    ,Pair=dsf$Pair
    ,Q=ifelse(is.na(dsf$Q), dsf$Rg * 2, dsf$Q)
    ,Ca=dsf$Ca
    ,ustar=dsf$ustar
    ,u=dsf$u
    ,ra = ra
    ,Mden = Mden, GPP_max = GPP_max, Dmax = Dmax
    ,GPP_sd_threshold = GPP_sd_threshold
    ,VPD_plant = dsf$VPD/10, # use VPD here instead of leaf temp
    ,config = cfg
  )
  # answer from previous run, might be wrong
  #with Priego original cost function: expect_equal( ans, 5.927, tolerance = 0.001 )
  expect_equal( ans, 85.83, tolerance = 0.01 )
})

test_that("optimize_function",{
  tz = attr(FIHyy$timestamp, "tzone")
  dsf <-  filter(FIHyy, between(
    timestamp, ISOdatetime(2010,5,15,0,0,0,tz=tz), ISOdatetime(2010,5,19,0,0,0,tz=tz)))
    # mutate(
    #   GPP_sd = ifelse(is.na(.data$GPP_sd), .data$GPP*0.1,.data$GPP_sd),
    #   Q = ifelse(is.na(.data$Q)==TRUE, .data$Rg*2, .data$Q)
    # )
  cfg = priego_config(burninlength = 1200, niter = 2000, updatecov = 250)
  ans <-  ETPart:::optim_priego(dsf, chi_o = 0.69, WUE_o= 9.03, config = cfg)
  expect_equal(length(ans$paropt), 4)
  # from oscarperezpriego/ETpartitioning:
  ans$paropt
  #expect_equal( ans$paropt[1], c(a1=275), tolerance = 0.5)
  expect_equal( ans$paropt[1], c(a1=396), tolerance = 0.5)
})

.tmp.f <- function(){
  cfg = priego_config(
    burninlength = 1200, niter = 200, updatecov = 250, GPPs2_bias2obs = 0.0)
  cfg = priego_config(
    burninlength = 500, niter = 8000, updatecov = 250, GPPs2_bias2obs = 0.05)
  ans <-  ETPart:::optim_priego(dsf, chi_o = 0.69, WUE_o= 9.03, config = cfg)
  plot(ans$parMCMC)
  pairs(ans$parMCMC)
  #plot(density(-2*log(ans$parMCMC$SS)))
  plot(ans$parMCMC$SS ~ ans$parMCMC$pars[,4])
  plot(density(ans$parMCMC$pars[,4]))
  plot(ans$parMCMC$SS)
  ans$parMCMC$bestpar
  summary(ans$parMCMC)
}

test_that("predict_transpiration_opriego",{
  tz = attr(FIHyy$timestamp, "tzone")
  dsf <-  filter(FIHyy, between(
    timestamp, ISOdatetime(2010,5,15,0,0,0,tz=tz),
    ISOdatetime(2010,5,19,0,0,0,tz=tz)))
  paropt <- c(a1 = 254, Do = 0.269, Topt = 19.6, beta = 0.555)
  transpiration_mod <- predict_transpiration_opriego(
    dsf, paropt, chi_o = 0.69, WUE_o= 9.03)$T
  # regression to previous value
  expect_equal(mean(transpiration_mod), 0.041, tolerance = 0.001)
})

test_that("estimate_T_priego_5days",{
  #tmp <-  FIHyy[ FIHyy$TIMESTAMP_START>  201105150000,]
  #tmp <-  tmp[tmp$TIMESTAMP_START<  201105160000,]
  tz = attr(FIHyy$timestamp, "tzone")
  lt  <- calculate_longterm_leaf(FIHyy, altitude = 200)
  data <-  filter(FIHyy, between(
    timestamp, ISOdatetime(2010,5,15,0,0,0,tz=tz), ISOdatetime(2010,5,22,0,0,0,tz=tz)))
  data <- data %>%
    #filter(NIGHT == 0) %>%  # filter only during fitting
    # -1 to associate midnight to previous day
    mutate(cumday = ETPart:::get_cumulative_day(timestamp-1)) %>%
    filter(cumday <= 6)
  cfg = priego_config(niter = 2000, burninlength = 1000, updatecov = 250)
  cfg = priego_config(GPPs2_bias2obs = 0.08)
  dfans <- map_dfr(unique(data$cumday), function(iday){
    ETPart:::estimate_T_priego_5days(
      data, iday, lt$chi_o, lt$WUE_o, config = cfg)$data
  })
  #expect_equal(nrow(dfans),nrow(data))
  expect_equal(select(dfans, -c(T, T_sd, GPP_pred, GPP_pred_sd)), data)
  #plot(T ~ ET, data=dfans, xlab="ET (mm/hour)"); abline(0,1)
  # regression to previous value: same magnitude
  #expect_equal(coef(lm(T~ET,data=dfans))["ET"], c(ET=0.5), tolerance= 0.1)
  .tmp.f <- function(){
    plot(T ~ timestamp, dfans, type="l")
    lines(T+1.96*T_sd ~ timestamp, dfans, type="l", col="gray")
    lines(T-1.96*T_sd ~ timestamp, dfans, type="l", col="gray")
    plot(GPP_pred ~ GPP, dfans); abline(0,1)
    plot((GPP_pred - GPP) ~ GPP_sd, dfans)
    # bias ratio
    s2m = (dfans$GPP_pred - dfans$GPP)^2 - dfans$GPP_sd^2
    #plot(density(s2m, na.rm = TRUE))
    mean(s2m, na.rm = TRUE)
    mean(s2m, na.rm = TRUE) / mean(dfans$GPP_sd^2, na.rm = TRUE)
    #0.08 # add 8% bias to observation error
  }
})

test_that("estimate_T_priego",{
  skip_on_cran()
  #tmp <-  FIHyy[ FIHyy$TIMESTAMP_START>  201105150000,]
  #tmp <-  tmp[tmp$TIMESTAMP_START<  201105160000,]
  tz = attr(FIHyy$timestamp, "tzone")
  lt  <- calculate_longterm_leaf(FIHyy, altitude = 200)
  data <-  filter(FIHyy, between(
    timestamp, ISOdatetime(2010,5,15,0,0,0,tz=tz), ISOdatetime(2010,5,22,0,0,0,tz=tz)))
  res <- ETPart:::estimate_T_priego(
    data, chi_o = lt$chi_o, WUE_o = lt$WUE_o, is_verbose = FALSE)
  dfans <- res$data
  #expect_equal(nrow(dfans),nrow(data))
  expect_equal(select(dfans, -c(
    .data$T,.data$T_sd, .data$GPP_pred, .data$GPP_pred_sd, .data$cumday)), data)
  #plot(T ~ ET, data=dfans, xlab="ET (mm/hour)"); abline(0,1)
  #plot(beta ~ cumday, res$paropt)
  # regression to previous value: same magnitude
  #expect_equal(coef(lm(T~ET,data=dfans))["ET"], c(ET=0.5), tolerance= 0.1)
})



.tmp.f <- function(){
  # entire FIHyy dataset
  lt <- calculate_longterm_leaf(FIHyy, altitude = 100)
  data <- FIHyy %>%
    #filter(NIGHT == 0) %>%  # filter only during fitting
    # -1 to associate midnight to previous day
    mutate(yday = ETPart:::get_cumulative_day(timestamp-1)) %>%
    filter(between(yday, 250, 260))
  res <- ETPart:::estimate_T_priego(data, chi_o = lt$chi_o, WUE_o = lt$WUE_o)
  dfT <- resopt$data
  plot(T ~ ET, data=dfT, xlab="ET (mm/hour)"); abline(0,1)
  plot(beta ~ cumday, data=resopt$paropt)
  plot(T ~ ET, data=filter(dfT, between(cumday,4,7)), xlab="ET (mm/hour)"); abline(0,1)
  plot(T ~ ET, data=filter(dfT, between(cumday,1,3)), xlab="ET (mm/hour)"); abline(0,1)
  plot(T ~ timestamp, data=dfT, type = "l")
}









