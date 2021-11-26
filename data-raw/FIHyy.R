#library(REddyProc)
library(dplyr)
library(bigleaf)
#library(ncdfTools) # https://github.com/bgctw/ncdfTools
library(readr)



berkeley_conv <- ETPart:::get_berkeley_conversion()
colspec <- do.call("cols_only", c(list(
  TIMESTAMP_END = col_character()),
  as.list(structure(
  rep("d", length(berkeley_conv$filecolname)),
  names = berkeley_conv$filecolname))))
df0 <- read_csv(
  file.path("tmp","FIHyy.csv"), col_types = colspec, na = c("", "NA", "-9999"))

FIHyy <- construct_data_FN15_(df0)
# df0 <- df0 %>% select(c("TIMESTAMP_END",berkeley_conv$filecolname))
# names(df0)[-1] <- berkeley_conv$colname
#
# df <- df0 %>% mutate(
#   timestamp = REddyProc::BerkeleyJulianDateToPOSIXct(df0$TIMESTAMP_END) - 15*60
#   , ET = LE.to.ET(LE, TA) * 60*30# from bigleaf convert mm/sec to mm/half-hour
#   , ET_sd = LE_RANDUNC * ET/LE
#   , NEE_sd = NEE_RANDUNC
#   , RH = bigleaf::VPD.to.rH(VPD/10, TA)
#   , quality_flag = (NEE_QC == 0) & (LE_QC == 0) &
#       is.finite(ET) & is.finite(GPP_NT)
#   )
#
# FIHyy <- df %>% select(
#   timestamp
#   , ET, GPP = GPP_DT, RH, VPD, Rg = SW_IN, Rg_pot = SW_IN_POT, Tair = TA
#   , ET_sd, NEE_sd
#   , precip = P, u = WS, quality_flag
# )

.tmp.f <- function(){
  df <- FIHyy
  nrec <- nrow(df)
  nday <- nrec/48
  df$iday <- rep(1:nday, each = 48)
  dfday <- df %>% group_by(iday) %>% summarize(
    ET = sum(ET, na.rm = FALSE)
    ,LE = sum(LE, na.rm = TRUE)
    ,timestamp = first(timestamp)
    )
  summary(dfday)
  plot(ET ~ timestamp, dfday)
  plot(LE ~ timestamp, dfday)
}

.tmp.save <- function(){
  usethis::use_data(FIHyy, overwrite = TRUE)
  saveRDS(FIHyy, file.path("tmp","FIHyy.rds"))
}
#
