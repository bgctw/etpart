.tmp.f <- function(){
  library(testthat)
  library(dplyr)
}

test_that("calculate_longterm_leaf with wrong column names",{
  expect_error(
    ans <- calculate_longterm_leaf(
      EddySample, "GPP_NT_VUT_USTAR50_bla", "VPD_F?_foobar", "TA_F",
      C = 1.1, Z = 0.2)
    ,"GPP_NT_VUT_USTAR50_bla")
})

test_that("calculate_longterm_leaf",{
  # dataTest = structure(list(Photos = c(-2.40437006950378, -2.74130988121033,
  #                                  -2.457279920578, -2.62785005569458, -2.72841000556946, -2.72528004646301
  # ), VPD = c(0.132400000095367, 0.233400011062622, 0.25220000743866,
  #            0.21470000743866, 0.225699996948242, 0.367700004577637)
  # , Tair = c(6.19799995422363,
  #            3.50399994850159, 3.53900003433228, 3.46900010108948, 4.38700008392334,
  #            6.20100021362305))
  # , .Names = c("Photos", "VPD", "Tair"), row.names = c(NA,
  #                                                      6L), class = "data.frame")
  lt  <- calculate_longterm_leaf(EddySample, C = 1.1, Z = 0.2)
  # answer from previous run, might be wrong
  # bug in orgiginal code devided VPD twice by 10
  #expect_equal( lt$chi_o, 0.88, tolerance = 0.1 )
  expect_equal( lt$chi_o, 0.67, tolerance = 0.01 )
  #expect_equal( lt$WUE_o, 24.5, tolerance = 0.1 )
  expect_equal( lt$WUE_o, 7.0, tolerance = 0.1 )
})

test_that("compute_stomatal_conductance_jarvis",{
  dsf <-  EddySample[EddySample$TIMESTAMP_START > 201405150000,]
  dsf <-  dsf[dsf$TIMESTAMP_START < 201405160000,]
  ans <- mean(ETPart:::compute_stomatal_conductance_jarvis(
    par=c(200, 0.2, 25), Q=dsf$PPFD_IN, VPD=dsf$VPD_F, Tair=dsf$TA_F, gcmax = 1))
  # answer from previous run, might be wrong
  expect_equal( ans, 0.26, tolerance = 0.1 )
})

test_that("cost_optim_opriego",{
  dsf <-  EddySample[ EddySample$TIMESTAMP_START>  201405150000,]
  dsf <-  dsf[dsf$TIMESTAMP_START<  201405160000,]
  dsf <- filter(dsf, NEE_VUT_USTAR50_JOINTUNC>0 & PPFD_IN>0 & SW_IN_F>0)
  ra <- ETPart:::compute_aerodynamic_conductance(dsf$WS_F, dsf$USTAR)
  ans <- ETPart:::cost_optim_opriego(
    par=c(200, 0.2, 25, 0.6),
    Chi_o = 0.88, WUE_o= 22.5
    ,Photos = dsf$GPP_NT_VUT_MEAN
    ,Photos_unc = dsf$NEE_VUT_USTAR50_JOINTUNC
    ,H=dsf$H_F_MDS
    ,VPD=dsf$VPD_F/10
    ,Tair=dsf$TA_F
    ,Pair=dsf$PA_F
    ,Q=ifelse(is.na(dsf$PPFD_IN), dsf$SW_IN_F * 2, dsf$PPFD_IN)
    ,Ca=dsf$CO2_F_MDS
    ,Ustar=dsf$USTAR
    ,WS=dsf$WS_F
    ,ra_c = ra$ra_c, ra_w = ra$ra_w
    )
  # answer from previous run, might be wrong
  expect_equal( mean(ans, na.rm=T), 4.1, tolerance = 0.1 )
})

test_that("optimize_function",{
  #tmp <-  EddySample[ EddySample$TIMESTAMP_START>  201105150000,]
  #tmp <-  tmp[tmp$TIMESTAMP_START<  201105160000,]
  dsf <-  EddySample[ EddySample$TIMESTAMP_START>  201405150000,]
  dsf <-  dsf[dsf$TIMESTAMP_START<  201405160000,]
  ans <-  optim_priego(par_lower= c(0,0, 10, 0)
                             ,par_upper = c(400,0.4, 30, 1)
                             ,data=dsf
                             ,Chi_o = 0.88, WUE_o= 22.5)
  expect_equal(length(ans), 4)
  # answer from previous run, might be wrong
  #expect_equal( as.numeric(ans[1]), 382.0955, tolerance = 20 )
  # from ETpartitioning:
  expect_equal( as.numeric(ans[1]), 211, tolerance = 40)
})









