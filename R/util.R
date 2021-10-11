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
