#' Extracts the site metadata from fluxdatasite url
#'
#' @param site site name, e.g. FI-Hyy
#' @param fluxdatasite url to be appended by /site
#'
#' @return list named by first entry in html-table and entries
#'   contact.name: vector of strings, contact.email: vector of strings parsed
#'   from entry Tower.Team
#' @export
get_fluxnet_metadata <- function(site, fluxdatasite = "http://sites.fluxdata.org") {
  if (!require("rvest")) stop(
    "package rvest needs to be installed for using get_fluxnet_metadata")
  url <- paste0(fluxdatasite,"/",site)
  page <- rvest::read_html(url)
  maininfo <- html_table(page, fill = TRUE)[[1]]
  labels <- gsub(":$","",trimws(maininfo[[1]]))
  values <- trimws(maininfo[[2]])
  as.numeric.try <- function(val){
    numval = suppressWarnings(as.numeric(val))
    if (!is.na(numval)) numval else val
  }
  metadata <- structure(lapply(values, as.numeric.try), names = make.names(labels))
  towerteam0 <-  strsplit(metadata$Tower.Team,"\n")[[1]] %>% trimws()
  towerteam <- towerteam0[(sapply(towerteam0, nzchar))] # skip empty lines
  metadata$contact.name = gsub("^(.+): (.+) <(.+)>(.*)$","\\2", towerteam)
  metadata$contact.email = gsub("^(.+): (.+) <(.+)>(.*)$","\\3", towerteam)
  metadata
}

get_berkeley_conversion <- function() {
  tribble(
  ~filecolname, ~colname,
 'CO2_F_MDS', 'CO2',
 'CO2_F_MDS_QC', 'CO2_QC',
 'GPP_DT_VUT_USTAR50', 'GPP_DT',
 'GPP_NT_VUT_USTAR50', 'GPP_NT',
 'G_F_MDS', 'G',
 'G_F_MDS_QC', 'G_QC',
 'H_CORR', 'H_CORR',
 'H_F_MDS', 'H',
 'H_F_MDS_QC', 'H_QC',
 'H_RANDUNC', 'H_RANDUNC',
 'H_RANDUNC_N', 'H_RANDUNC_N',
 'LE_CORR', 'LE_CORR',
 'LE_F_MDS', 'LE',
 'LE_F_MDS_QC', 'LE_QC',
 'LE_RANDUNC', 'LE_RANDUNC',
 'LE_RANDUNC_N', 'LE_RANDUNC_N',
 'LW_IN_F_MDS', 'LW_IN',
 'LW_OUT', 'LW_OUT',
 'NEE_VUT_USTAR50', 'NEE',
 'NEE_VUT_USTAR50_QC', 'NEE_QC',
 'NEE_VUT_USTAR50_RANDUNC', 'NEE_RANDUNC',
 'NEE_VUT_USTAR50_RANDUNC_N', 'NEE_RANDUNC_N',
 'NETRAD', 'NETRAD',
 'NIGHT', 'NIGHT',
 'P', 'P',
 'PA', 'PA',
 'PPFD_DIF', 'PPFD_DIF',
 'PPFD_IN', 'PPFD_IN',
 'PPFD_OUT', 'PPFD_OUT',
 'RECO_DT_VUT_USTAR50', 'RECO_DT',
 'RECO_NT_VUT_USTAR50', 'RECO_NT',
 'SWC_F_MDS_1', 'SWC_1',
 'SWC_F_MDS_1_QC', 'SWC_1_QC',
 'SWC_F_MDS_2', 'SWC_2',
 'SWC_F_MDS_2_QC', 'SWC_2_QC',
 'SW_IN_F_MDS', 'SW_IN',
 'SW_IN_F_MDS_QC', 'SW_IN_QC',
 'SW_OUT', 'SW_OUT',
 'SW_IN_POT', 'SW_IN_POT',
 'TA_F_MDS', 'TA',
 'TA_F_MDS_QC', 'TA_QC',
 'TS_F_MDS_1', 'TS_1',
 'TS_F_MDS_2', 'TS_2',
 'TS_F_MDS_1_QC', 'TS_1_QC',
 'TS_F_MDS_2_QC', 'TS_2_QC',
 'USTAR', 'USTAR',
 'VPD_F_MDS', 'VPD',
 'VPD_F_MDS_QC', 'VPD_QC',
 'WS_F', 'WS'
  )
}
