library(REddyProc)
library(dplyr)
library(bigleaf)
library(ncdfTools) # https://github.com/bgctw/ncdfTools
library(readr)


df0 <- read.csv(file.path("tmp","FIHyy.csv")) %>%
df %>% df0 %>%
  select(timestamp = TIMESTAMP_END, "NEE","NEE_f","Rg",)
summary(df)

df <- df %>%
  mutate(
    ET = LE.to.ET(LE, Tair) # from bigleaf
    , precip = dfprecip$precip
  )

# partition GPP
EProc <- sEddyProc$new('DE-Tha', df,
                       c('NEE','Rg','Tair','VPD', 'Ustar','ET'))
EProc$sMDSGapFillAfterUstar('NEE', uStarTh = 0.45, uStarSuffix = "")
EProc$sSetLocationInfo(LatDeg = 51.0, LongDeg = 13.6, TimeZoneHour = 1)
EProc$sMDSGapFill('Tair', FillAll = FALSE,  minNWarnRunLength = NA)
EProc$sMDSGapFill('VPD', FillAll = FALSE,  minNWarnRunLength = NA)
EProc$sMRFluxPartition()
EProc$sMDSGapFill('ET') # gapfill ET to get model, errors, and sd
tmp <- EProc$sExportResults()
df <- df %>% mutate(
  GPP = tmp$GPP_f
  ,GPP_sd = tmp$NEE_fsd # use the uncertainty of NEE for GPP_sd
  ,NEE_fall = tmp$NEE_fall
  ,Rg_pot = tmp$PotRad
  ,ET_fall = tmp$ET_fall
  ,ET_sd = tmp$ET_fsd
  ,Tair_f = tmp$Tair_f # use gap-filled versions
  ,VPD_f = tmp$VPD_f
)

DETha <- df %>% select(
  timestamp = DateTime, ET, NEE, GPP, Tair = Tair_f, rH, VPD = VPD_f,
  Rg, Rg_pot, u = Ustar,
  GPP_sd, NEE_fall, ET_fall, ET_sd)

usethis::use_data(DETha, overwrite = TRUE)
