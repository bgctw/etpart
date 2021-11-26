get_cumulative_day <- function(
  ### get the cumulative day from a single equidistant time series
  datetime  ##<< POSIXct vector
){
  ##note<< day is computed in associated time zone. Hence, convert UTC to local
  ## time before
  ##value<< cumulative day starting
  # only works with equidistant where no days are missing
  #cumDay <- cumsum(c(FALSE,diff(as.POSIXlt(dateHalfHour - 15L*60L)$yday) != 0))
  as.numeric(as.Date(datetime))-min(as.numeric(as.Date(datetime)), na.rm=T)+1
}
get_halfhour_ofday <- function(
  ### get the half-Hour within one day
  dateHalfHour  ##<< POSIXct vector, indicating the end
  ## a half-hour (00:30, 01:00, 01:30, ...)
){
  hour <- as.POSIXlt(dateHalfHour - 15L*60L)$hour
  min <- as.POSIXlt(dateHalfHour - 30L*60L)$min
  ##value<< 00:30 -> 1, 01:00 -> 2, ... , 23:30 -> 47, 00:00 -> 48
  hour*2L + min/30L + 1L
}

check_required_cols <- function(data, required_names) {
  iMissing <- which( !(required_names %in% names(data)))
  if (length(iMissing)) stop(
    "Need to provide columns ",
    paste(c(required_names)[iMissing], collapse = ", "))
}

replace_columns <- function(data, ..., add.unique = TRUE) {
  # from sjmisc
  # evaluate dots
  .dots <- match.call(expand.dots = FALSE)$`...`

  # if add_columns had no arguments, .dots are NULL
  # this crashes R when calling bind_cols
  if (is.null(.dots)) {
    stop("You must specify at least one more data frame with columns to add.", call. = FALSE)
  }

  # bind all data frames to one
  tmp <- dplyr::bind_cols(...)

  # check for identical column names
  data.doubles <- colnames(data) %in% colnames(tmp)
  tmp.doubles <- colnames(tmp) %in% colnames(data)

  # replace duplicate variables in "data" with duplicates from "..."
  data[, data.doubles] <- tmp[, tmp.doubles, drop = FALSE]

  # add remaining columns that were not duplicates
  if (add.unique)
    x <- dplyr::bind_cols(data, tmp[, !tmp.doubles, drop = FALSE])
  else
    x <- data

  x
}


