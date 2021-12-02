# Convenence functions to read standard datasets
# Will reduce support work with people redoing this.

#' Format Fluxnet2015 data.frame
#'
#' Construct data.frame formatted for partitioning from Fluxnet2015 format
#' https://fluxnet.org/data/fluxnet2015-dataset/
#'
#' @param data data frame with column names given by the following arguments.
#' @param GPP photosynthesis data (umol CO2 m-2 s-1).
#' @param GPP_sd photosynthesis uncertainties (umol CO2 m-2 s-1).
#' @param Ca atmospheric CO2 concentration (umol Co2 mol air-1).
#' @param H sensible heat flux (W m-2).
#' @param LE Latent heat flux (W m-2)
#' @param LE_QC quality flag for LE (0 measured, 1 good quality, ...)
#' @param LE_sd uncertainty of LE
#' @param NEE_QC quality flag for NEE (0 measured, 1 good quality, ...)
#' @param NIGHT Flag indicating nighttime interval (1 for night, zero otherwise)
#' @param precip P used if measured (mm per dataset resolution: either hour or half-hour)
#' @param Pair atmospheric pressure (kPa).
#' @param Q photosynthetic active radiation (umol m-2 s-1).
#' @param Rg Shortwave radiation, incoming (W/m2)
#' @param Rg_pot potential Shortwave radiation, incoming (W/m2)
#' @param Tair air temperature (deg C).
#' @param ustar wind friction velocity (m s-1).
#' @param VPD vapor pressure deficit (hPa).
#' @param u wind velocity (m s-1).
#'
#' @return data.frame with colums as given by argument names and
#' \item{timestamp}{POSIXct time: end of averaging period}
#' \item{ET}{evapotranspiration computed from LE in mm/hour}
#' \item{ET_sd}{uncertainty of ET}
#' \item{RH}{relative humidity (1/1) computed from VPD and Tair}
#' \item{isnight}{boolean TRUE if its nighttime computed from NIGHT}
#' \item{qualtiy_flag}{
#'   (NEE_QC == 0) & (LE_QC == 0) & is.finite(ET) & is.finite(GPP)}
#'
#' @importFrom bigleaf LE.to.ET
#' @importFrom bigleaf VPD.to.rH
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
      .data$TIMESTAMP_END,
      GPP=.data$GPP, GPP_sd=.data$GPP_sd, Ca=.data$Ca, H=.data$H, LE=.data$LE,
      LE_sd=.data$LE_sd, NIGHT=.data$NIGHT,
      precip=.data$precip, Pair=.data$Pair, Q=.data$Q, Rg=.data$Rg,
      Rg_pot=.data$Rg_pot, Tair=.data$Tair, ustar=.data$ustar,
      VPD=.data$VPD,u=.data$u
    ) %>%
    mutate(
      timestamp = BerkeleyJulianDateToPOSIXct(.data$TIMESTAMP_END) - 15*60
      #, ET = LE.to.ET(LE, Tair) * 60*30# from bigleaf convert mm/sec to mm/half-hour
      , ET = LE.to.ET(.data$LE, .data$Tair) * 60*60 # convert mm/sec to mm/hour
      , ET_sd = .data$LE_sd * .data$ET/.data$LE
      #, NEE_sd = NEE_RANDUNC
      , RH = VPD.to.rH(.data$VPD/10, .data$Tair)
      , isnight = (.data$NIGHT == 1)
      , quality_flag = (.data$NEE_QC == 0) & (.data$LE_QC == 0) &
        is.finite(.data$ET) & is.finite(.data$GPP)
    ) %>%
    select(-c(.data$NIGHT, .data$TIMESTAMP_END)) %>%
    select(
      .data$timestamp, .data$ET, .data$ET_sd, .data$GPP, .data$GPP_sd, .data$RH,
      everything())
}

#' convert POSIXct to JulianDate format used in Berkeley release
#'
#' @param sDateTime POSIXct vector
#' @param tz timezone to use
#'
#' @details
#' In the Berkeley-Release of the Fluxnet data, the time is stored as an number
#' with base10-digits representing YYYYMMddhhmm
#' @seealso \code{\link{BerkeleyJulianDateToPOSIXct}}
#'
#' @export
POSIXctToBerkeleyJulianDate <- function(sDateTime, tz = get_tzone(sDateTime)) {
  charRep <- strftime(sDateTime, format = "%Y%m%d%H%M", tz = tz)
  ans <- as.numeric(charRep)
  ans
}

#' convert JulianDate format used in Berkeley release to POSIXct
#'
#' @param julianDate numeric vector representing times (see details for format)
#' @param tz time zone used to represent the dates
#' @param ... further arguments to \code{\link{strptime}}
#'
#' @details
#' In the Berkeley-Release of the Fluxnet data, the time is stored as an number
#'  with base10-digits representing YYYYMMddhhmm
#'
#' @seealso
#' \code{\link{POSIXctToBerkeleyJulianDate}}
#' @export
BerkeleyJulianDateToPOSIXct <- function(julianDate, tz = "UTC", ...) {
  ans <- as.POSIXct(strptime(as.character(julianDate), "%Y%m%d%H%M", tz = tz, ...))
  ans
}

#' Extract the timezone attribute from POSIXct with default on missing
#'
#' @param x POSIXct vector
#' @param default time zone returned, if x has not timezone associated or
#' attribute is the zero string
#'
#' @examples
#' get_tzone(as.POSIXct("2010-07-01 16:00:00", tz = "etc/GMT-1") )
#' get_tzone(as.POSIXct("2010-07-01 16:00:00") )
#' # printed with local time zone, but actually has no tz attribute
#' get_tzone(Sys.time())
#'
#' @export
get_tzone <- function(x, default = "GMT") {
  tzone <- attr(x, "tzone")
  if (length(tzone) && nzchar(tzone)) tzone else default
}

