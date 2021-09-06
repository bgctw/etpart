.tmp.f <- function(){
  library(testthat)
  library(dplyr)
}

# data.frame spanning more than one year
set.seed(0815)
nday = 367
isprecip = rpois(48*nday, 0.005)
.dsin = -cos(1:(nday*48)*2*pi/48)
#plot(.sin)
df <- data.frame(
  timestamp = seq(
    ISOdate(2018, 1, 1, 0, tz = "Etc/GMT-1"), by = "30 min", length.out = 48*nday)
  ,GPP = pmax(0,.dsin)
  ,precip = isprecip * 5*rlnorm(48*nday)
  ,ET = 1/24*pmax(0,.dsin)
)

.tmp.f <- function() {
  plot(ET ~ timestamp, head(df, 3*48))
  plot(precip ~ timestamp, head(df, 48*48))
}

test_that("compute_daily_GPP", {
  GPPday = compute_daily_GPP(df$GPP, df$timestamp)
  expect_equal(length(GPPday), length(df$GPP))
  expect_equal(GPPday[1], sum(df$GPP[1:48]))
  expect_equal(GPPday[48], sum(df$GPP[1:48]))
  expect_equal(GPPday[49], sum(df$GPP[48 + 1:48]))
  #
  # check not starting with midnight
  dfs <- df[2:50,]
  GPPday = compute_daily_GPP(dfs$GPP, dfs$timestamp)
  expect_equal(length(GPPday), length(dfs$GPP))
  expect_equal(GPPday[1], sum(dfs$GPP[1:47]))
  expect_equal(GPPday[47], sum(dfs$GPP[1:47]))
  expect_equal(GPPday[48], sum(dfs$GPP[48:49]))
})

test_that("compute_cswi", {
  cswi = compute_cswi(df, smax = 5)
  expect_equal(length(cswi), nrow(df))
  df2 = mutate(df, cswi = cswi)
  .tmp.f <- function(){
    plot(cswi ~ timestamp, head(df2, 48*5))
    plot(cswi ~ timestamp, df2)
  }
})



