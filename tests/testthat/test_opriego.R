.tmp.f <- function(){
  library(testthat)
  library(dplyr)
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
  lt  <- calculate_longterm_leaf(FIHyy, C = 1.189, Z = 0.27)
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
  ra <- ETPart:::compute_aerodynamic_conductance(dsf$u, dsf$ustar)
  constants = list(
    Cp = 1003.5, ##<< heat capacity [J kg-1 K-1].
    R_gas_constant = 0.287058, ##<< [J kg-1 deg K-1].
    M = 0.0289644 ##<< molar mass, [kg mol-1].
  )
  dens = ETPart:::calculate_air_density(dsf$Pair, dsf$Tair, constants) #[kg m-3]
  Mden = dens/constants$M ##<< molar air density [mol m-3].
  GPP_max = quantile(dsf$GPP, probs=c(0.90), na.rm=T)
  Dmax = mean(dsf$VPD[dsf$GPP>GPP_max]/10, na.rm=T)
  GPP_sd_threshold = pmax(dsf$GPP*0.1, dsf$GPP_sd)
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
    ,VPD_plant = dsf$VPD/10 # use VPD here instead of leaf temp
  )
  # answer from previous run, might be wrong
  expect_equal( ans, 5.927, tolerance = 0.001 )
})

test_that("optimize_function",{
  tz = attr(FIHyy$timestamp, "tzone")
  dsf <-  filter(FIHyy, between(
    timestamp, ISOdatetime(2010,5,15,0,0,0,tz=tz), ISOdatetime(2010,5,19,0,0,0,tz=tz)))
  ans <-  ETPart:::optim_priego(
    dsf, chi_o = 0.88, WUE_o= 22.5,
    config = priego_config(burninlength = 500, niter = 1000))
  expect_equal(length(ans$paropt), 4)
  # from ETpartitioning:
  expect_equal( ans$paropt[1], c(a1=295), tolerance = 40)
})

test_that("5dayloop",{
  #tmp <-  FIHyy[ FIHyy$TIMESTAMP_START>  201105150000,]
  #tmp <-  tmp[tmp$TIMESTAMP_START<  201105160000,]
  dsf <-  FIHyy[ FIHyy$TIMESTAMP_START > 201405150000,]
  dsf <-  dsf[dsf$TIMESTAMP_START <=  201406150000,]
  dsf <- dsf %>% mutate(
    datetime = BerkeleyJulianDateToPOSIXct(.data$TIMESTAMP_END, "Etc/GMT+0"))
  data <- dsf
  lt  <- calculate_longterm_leaf(FIHyy, C = 1.1, Z = 0.2)
  data <- data %>%
    #filter(NIGHT == 0) %>%  # filter only during fitting
    # -1 to associate midnight to previous day
    mutate(cumday = ETPart:::get_cumulative_day(datetime-1)) %>%
    filter(cumday <= 6)
  dfans <- map_df(unique(data$cumday), function(iday){
    ETPart:::estimate_T_priego_5days(data, iday, lt$chi_o, lt$WUE_o)
  })
  data %>% group_by(cumday) %>% mutate()
})









