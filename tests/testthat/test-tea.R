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
  ,Rg = 200 * pmax(0,.dsin)
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

test_that("compute_diurnal_centroid", {
  #increase ET at 10:00am at second day
  # so that dci is shifted towards morning for second day
  df2 <- df; df2$ET[48+2*10] <- 2*df2$ET[48+2*10]
  dci = compute_diurnal_centroid(df2$ET)
  expect_equal(length(dci), nday)
  expect_equal(dci[1], 12)
  expect_true(dci[2] < 12)
  #
  ndci <- compute_norm_diurnal_centroid(df2$ET, df2$Rg)
  expect_equal(ndci[1], 0)
  expect_equal(ndci[2], dci[2] - 12)
})

test_that("daily_correlation", {
  #increase ET at 10:00am at second day
  # so that dci is shifted towards morning for second day
  df2 <- df; df2$ET[48+2*10] <- 2*df2$ET[48+2*10]
  dci = ETPart:::daily_corr(df2$ET, df2$GPP, df2$Rg)
  expect_equal(length(dci), nday)
  expect_equal(dci[1], 1.0)
  expect_true(dci[2] < 1.0)
})

test_that("compute_DWCI", {
  tmp <- compute_DWCI(DETha)
})




