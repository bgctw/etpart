# Convenence functions to read standard datasets
# Will reduce support work with people redoing this.

#' @export
construct_data_FN15 = function(
  data,
  GPP = 'GPP_NT_VUT_USTAR50',
  GPP_sd = 'NEE_VUT_USTAR50_RANDUNC',
  Ca = 'CO2_F_MDS',
  #Ca_QC = 'CO2_F_MDS_QC',
  #GPP_DT = 'GPP_DT_VUT_USTAR50',
  #GPP_NT = 'GPP_NT_VUT_USTAR50',
  #G = 'G_F_MDS',
  #G_QC = 'G_F_MDS_QC',
  #H_CORR = 'H_CORR',
  H = 'H_F_MDS',
  #H_QC = 'H_F_MDS_QC',
  #H_RANDUNC = 'H_RANDUNC',
  #H_RANDUNC_N = 'H_RANDUNC_N',
  #LE_CORR = 'LE_CORR',
  LE = 'LE_F_MDS',
  LE_QC = 'LE_F_MDS_QC',
  LE_sd = 'LE_RANDUNC',
  #LE_RANDUNC_N = 'LE_RANDUNC_N',
  #LW_IN = 'LW_IN_F_MDS',
  #LW_OUT = 'LW_OUT',
  #NEE = 'NEE_VUT_USTAR50',
  NEE_QC = 'NEE_VUT_USTAR50_QC',
  #NEE_RANDUNC = 'NEE_VUT_USTAR50_RANDUNC',
  #NEE_RANDUNC_N = 'NEE_VUT_USTAR50_RANDUNC_N',
  #NETRAD = 'NETRAD',
  NIGHT = 'NIGHT',
  precip = 'P',
  Pair = 'PA',
  #PPFD_DIF = 'PPFD_DIF',
  Q = 'PPFD_IN',
  #PPFD_OUT = 'PPFD_OUT',
  #RECO_DT = 'RECO_DT_VUT_USTAR50',
  #RECO_NT = 'RECO_NT_VUT_USTAR50',
  #SWC_1 = 'SWC_F_MDS_1',
  #SWC_1_QC = 'SWC_F_MDS_1_QC',
  #SWC_2 = 'SWC_F_MDS_2',
  #SWC_2_QC = 'SWC_F_MDS_2_QC',
  Rg = 'SW_IN_F_MDS',
  #Rg_QC = 'SW_IN_F_MDS_QC',
  #SW_OUT = 'SW_OUT',
  Rg_pot = 'SW_IN_POT',
  Tair = 'TA_F_MDS',
  #TA_QC = 'TA_F_MDS_QC',
  #TS_1 = 'TS_F_MDS_1',
  #TS_2 = 'TS_F_MDS_2',
  #TS_1_QC = 'TS_F_MDS_1_QC',
  #TS_2_QC = 'TS_F_MDS_2_QC',
  ustar = 'USTAR',
  VPD = 'VPD_F_MDS',
  #VPD_QC = 'VPD_F_MDS_QC',
  u = 'WS_F'
){
  data %>%
    select(
      TIMESTAMP_END,
      GPP=GPP, GPP_sd=GPP_sd, Ca=Ca, H=H, LE=LE, LE_sd=LE_sd, NIGHT=NIGHT,
      precip=precip, Pair=Pair, Q=Q, Rg=Rg, Rg_pot=Rg_pot, Tair=Tair, ustar=ustar,
      VPD=VPD,u=u
    ) %>%
    mutate(
      timestamp = BerkeleyJulianDateToPOSIXct(TIMESTAMP_END) - 15*60
      #, ET = LE.to.ET(LE, Tair) * 60*30# from bigleaf convert mm/sec to mm/half-hour
      , ET = LE.to.ET(LE, Tair) * 60*60# from bigleaf convert mm/sec to mm/half-hour
      , ET_sd = LE_sd * ET/LE
      #, NEE_sd = NEE_RANDUNC
      , RH = bigleaf::VPD.to.rH(VPD/10, Tair)
      , isnight = (NIGHT == 1)
      , quality_flag = (NEE_QC == 0) & (LE_QC == 0) &
        is.finite(ET) & is.finite(GPP)
    ) %>%
    select(-c(NIGHT, TIMESTAMP_END)) %>%
    select( timestamp, ET, ET_sd, GPP, GPP_sd, RH, everything())
}

#' @export
POSIXctToBerkeleyJulianDate <- function(
  ### convert POSIXct to JulianDate format used in Berkeley release
  sDateTime  ##<< POSIXct vector
  ,tz = getTZone(sDateTime)
) {
  ##author<< TW,
  ##seealso<< \code{\link{BerkeleyJulianDateToPOSIXct}},
  ##details<<
  ## In the Berkeley-Release of the Fluxnet data, the time is stored as an number
  ## with base10-digits representing YYYYMMddhhmm
  charRep <- strftime(sDateTime, format = "%Y%m%d%H%M", tz = tz)
  ans <- as.numeric(charRep)
  ans
}

#' @export
BerkeleyJulianDateToPOSIXct <- function(
  ### convert JulianDate format used in Berkeley release to POSIXct
  julianDate  ##<< numeric vector representing times (see details for format)
  , tz = "UTC"  ##<< time zone used to represent the dates
  , ...    ##<< further arguments to \code{\link{strptime}}
) {
  ##author<< TW,
  ##seealso<< \code{\link{POSIXctToBerkeleyJulianDate}}
  ## \code{\link{fConvertTimeToPosix}}
  ##details<<
  ## In the Berkeley-Release of the Fluxnet data, the time is stored as an number
  ## with base10-digits representing YYYYMMddhhmm
  ans <- as.POSIXct(strptime(as.character(julianDate), "%Y%m%d%H%M", tz = tz, ...))
  ans
}

