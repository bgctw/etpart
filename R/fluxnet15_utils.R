# Convenence functions to read standard datasets
# Will reduce support work with people redoing this.

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


#' Read a file in the format of Fluxnet 2015 release
#'
#' Assigns default units to the columns and keeps variable name attributes
#' as in original file.
#'
#' @param file_path scalar string: the path to the csv file
#' @param additional_columns character vector of columns to
#'   read in addition of standard columns of \code{\link{read_from_fluxnet15}}.
#'   Can be a character vector or a object return by \code{\link{cols}}
#' @param colname_NEE name (scalar string) of column that reports NEE observations
#' @param ... further arguments to \code{\link{read_csv}}
#'
#' @examples
#'   ds_fn15 <- Example_DETha98 %>%
#'      fConvertTimeToPosix('YDH',Year = 'Year',Day = 'DoY', Hour = 'Hour') %>%
#'      dplyr::mutate(
#'         TIMESTAMP_END = POSIXctToBerkeleyJulianDate(.data$DateTime),
#'         season = factor(199801)
#'      ) %>%
#'      dplyr::rename(SW_IN = .data$Rg, TA = .data$Tair, USTAR = .data$Ustar) %>%
#'      dplyr::select(dplyr::one_of(c(
#'        "TIMESTAMP_END","NEE","SW_IN","TA","VPD","USTAR","season")))
#'   head(ds_fn15)
#'   fname <- tempfile()
#'   readr::write_csv(ds_fn15, fname)
#'
#'   # standard columns are renamed to REddyProc defaults
#'   ds_eproc <- fLoadFluxnet15(fname)
#'   head(ds_eproc)
#'   EProc <- sEddyProc$new("DE-Tha", ds_eproc)
#'   head(EProc$sExportData())
#'
#'   # Additional columns can be specified, e.g. factor column season
#'   ds_eproc <- fLoadFluxnet15(fname,
#'     additional_columns = readr::cols(season = readr::col_factor()))
#'   head(ds_eproc)
#'   EProc <- sEddyProc$new("DE-Tha", ds_eproc,
#'     c("NEE", "Rg", "Tair", "VPD", "Ustar","season"),
#'     ColNamesNonNumeric = "season"
#'     )
#'   head(EProc$sExportData())
#' @export
fLoadFluxnet15 <- function(file_path, additional_columns = character(0),
                           colname_NEE = "NEE", ...) {
  col <- col_standard <- cols_only(
    TIMESTAMP_END = col_character(),
    NEE = col_double(),
    LE = col_double(),
    H = col_double(),
    SW_IN = col_double(),
    TA = col_double(),
    TS = col_double(),
    USTAR = col_double(),
    VPD = col_double()
  )
  df_units <- tribble(
    ~varname, ~unit,
    colname_NEE, "umolm-2s-1",
    "LE", "Wm-2",
    "H", "Wm-2",
    "SW_IN", "Wm-2",
    "TA", "degC",
    "TS", "degC",
    "USTAR", "ms-1",
    "VPD", "hPa",
  )
  names(col_standard$cols)[2] <- colname_NEE
  colsInFile <- read_lines(file_path, n_max = 1L) %>% strsplit(",") %>%  "[["(1)
  col$cols <- col$cols[names(col$cols) %in% colsInFile]
  if (length(additional_columns)) {
    col_add <- if (inherits(additional_columns,"col_spec")) {
      additional_columns
    } else {
      col_add <- cols(rep(col_guess(),length(additional_columns)))
      names(col_add$cols) <- additional_columns
      col_add
    }
    col$cols <- c(col_standard$cols, col_add$cols)
  }
  # df_fn15 <- read_csv(file_path, ...)
  # df_fn15 <- read_csv(file_path, col_types = col_standard, ...)
  df_fn15 <- read_csv(file_path, col_types = col, ...)
  df_unitsin <- df_units %>% filter(.data$varname %in% names(df_fn15))
  df_fn15 <- df_fn15 %>% as.data.frame() %>%
    set_varunit_attributes(df_unitsin$varname, df_unitsin$unit)
  read_from_fluxnet15(df_fn15, colname_NEE = colname_NEE)
}

.tmp.f <- function(){
  ds <- Example_DETha98 %>%
    fConvertTimeToPosix('YDH',Year = 'Year',Day = 'DoY', Hour = 'Hour') %>%
    rename(SW_IN = .data$Rg, TA = .data$Tair, USTAR = .data$Ustar) %>%
    mutate(TIMESTAMP_END = POSIXctToBerkeleyJulianDate(.data$DateTime)) %>%
    select(one_of(c("TIMESTAMP_END","NEE","SW_IN","TA","VPD","USTAR")))
  fname <- tempfile()
  write_csv(ds, fname)
  ds_fn15 <- fLoadFluxnet15(fname)
  head(ds_fn15)
  EProc <- sEddyProc$new("DE-Tha", ds_fn15)
  head(EProc$sExportData())
}

#' extract REddyProc input columns from data.frame in Fluxnet15 format
#'
#' Column format as described at
#' https://fluxnet.org/data/fluxnet2015-dataset/fullset-data-product/
#'
#' @details If input has numeric column USTAR_QC then USTAR of records
#' with USTAR_QC > 2 are set to NA.
#'
#' @param ds data.frame with columns TIMESTAMP_END (Time YYYYMMDDHHMM),
#' NEE, LE, H, USTAR, TA, TS, VPD, SW_IN and optionally USTAR_QC
#' @param colname_NEE name (scalar string) of column that reports NEE observations
#' @return data.frame with additional columns 'DateTime', 'NEE','Ustar' and
#'   'Rg','Tair','Tsoil' if columns 'SW_IN','TA', or 'TS' are present respectively
#' @export
read_from_fluxnet15 <- function(ds, colname_NEE = "NEE"){
  ustar_qc <- if (
    ("USTAR_QC" %in% names(ds)) && is.numeric(ds$USTAR_QC)) ds$USTAR_QC else
      rep(1L, nrow(ds))
  ds_eproc <- ds %>% mutate(
    DateTime = BerkeleyJulianDateToPOSIXct(.data$TIMESTAMP_END),
    NEE = .data[[colname_NEE]],
    Ustar = ifelse(ustar_qc <= 2L, .data$USTAR, NA_real_)
  )

  if ("TA" %in% names(ds_eproc)) ds_eproc <- ds_eproc %>% mutate( Tair = .data$TA)
  if ("TS" %in% names(ds_eproc)) ds_eproc <- ds_eproc %>% mutate( Tsoil = .data$TS)
  if ("SW_IN" %in% names(ds_eproc)) ds_eproc <- ds_eproc %>% mutate( Rg = .data$SW_IN)
  ds_eproc
}



#' extract processing results in Fluxnet15 format
#'
#' extract processing results with columns corresponding to Fluxnet15 release
#'
#' @param EProc sEddyProc class with uncertainty also in meteo variables and
#' both nighttime and daytime partitioning columns present
#' @param keep_other_cols set to TRUE to report also other columns
#' @param is_export_nonfilled set to FALSE to not export columns before gapfilling
#'
#' @return data.frame with columns names of Fluxnet15. Timestamps are
#'   in ISO string format \code{\link{POSIXctToBerkeleyJulianDate}}
#' @export
extract_FN15 <- function(EProc = .self, is_export_nonfilled = TRUE, keep_other_cols = FALSE) {
  input <- bind_cols(EProc$sDATA, select(EProc$sTEMP, -1L))
  time <- EProc$sDATA$sDateTime
  timestep <- difftime(time[2],time[1], units = "hours")
  output_time <- tibble(
    TIMESTAMP_START = POSIXctToBerkeleyJulianDate(time - timestep/2),
    TIMESTAMP_END = POSIXctToBerkeleyJulianDate(time + timestep/2)
  )
  # do not import stringr for dependencies
  str_replace <- function(x,pattern,replacement) gsub(pattern, replacement, x)
  replaceFun <- function(pattern, replacement,...){
    tmp <- input %>% select(matches(pattern))
    names(tmp) <- names(tmp) %>% str_replace(pattern, replacement)
    tmp
  }
  replace_patterns_uStar <- tribble(
    ~pattern, ~replacement, ~variable, ~method,
    "^NEE_U(\\d\\d)_f$", "NEE_VUT_\\1", "NEE", "",
    "^GPP_U(\\d\\d)_f$", "GPP_NT_VUT_\\1", "GPP", "NT",
    "^GPP_DT_U(\\d\\d)$", "GPP_DT_VUT_\\1", "GPP", "DT",
    "^Reco_U(\\d\\d)$", "RECO_NT_VUT_\\1", "GPP", "NT",
    "^Reco_DT_U(\\d\\d)$", "RECO_DT_VUT_\\1", "GPP", "DT",
    "^Ustar_Thresh_U(\\d\\d)$", "USTAR_THRESHOLD_VUT_\\1", "", "",
    "^NEE_U(\\d\\d)_fqc$", "NEE_VUT_USTAR\\1_QC", "NEE", "",
    # the following extract only one column
    # by putting them as a pattern it works also if the column does not exist
    #better add fqc for all ustar "^NEE_U50_fqc$", "NEE_VUT_USTAR50_QC", "NEE", "",
    "^NEE_U50_fsd$", "NEE_VUT_USTAR50_RANDUNC", "NEE", "",
    "^NEE_U50_fnum$", "NEE_VUT_USTAR50_RANDUNC_N", "NEE", "",
    "^FP_qc$", "GPP_DT_U50_QC", "GPP", ""		##<< quality flag: 0: good parameter fit,
    ## 1: some parameters out of range, required refit,
    ## 2: next parameter estimate is more than two weeks away
  )
  output_ustar <- replace_patterns_uStar %>% pmap(replaceFun) %>% bind_cols()
  if (!ncol(output_ustar)) output_ustar <- NULL
  replace_patterns_filled <- tribble(
    ~pattern, ~replacement, ~variable, ~method,
    "^night$", "NIGHT", "", "",
    "^Rg_f$", "SW_IN_F_MDS", "Rg", "",
    "^Rg_fqc$", "SW_IN_F_MDS_QC", "Rg", "",
    "^PotRad_NEW$", "SW_IN_POT", "Rg", "",
    "^Tair_f$", "TA_F_MDS", "", "",
    "^Tair_fqc$", "TA_F_MDS_QC", "", "",
    "^VPD_f$", "VPD_F_MDS", "", "",
    "^VPD_fqc$", "VPD_F_MDS_QC", "", ""
  )
  output_filled <- replace_patterns_filled %>% pmap(replaceFun) %>% bind_cols()
  if (!ncol(output_filled)) output_filled <- NULL
  replace_patterns_orig <- tribble(
    ~pattern, ~replacement, ~variable, ~method,
    "^NEE$", "NEE_ORIG", "NEE", "",
    "^Rg$", "SW_IN", "Rg", "",
    "^Tair$", "TA", "", "Tair",
    "^Ustar$", "USTAR", "uStar", "",
    "^VPD$", "VPD", "VPD", ""
  )
  output_orig <- if (isTRUE(is_export_nonfilled)){
    replace_patterns_orig %>% pmap(replaceFun) %>% bind_cols()
  } else NULL
  if (!ncol(output_orig)) output_orig <- NULL
  #
  output <- bind_cols(output_time, output_orig, output_filled, output_ustar)
  if (isTRUE(keep_other_cols)) output <- bind_cols(
    output, input[,setdiff(names(input),names(output)),drop = FALSE])
  output
}



