#' partition ET by TEA
#'
#' The Transpiration Estimation Algorithm (TEA) partitions ET by first
#' constraining the data to dry conditions where T/ET~1
#' (see \code{\link{tea_filter}}) and then learning water use efficiency (WUE)
#' relationship with predictors from the resulting training data
#' (see \code{link{tea_fit_wue}}).
#' The learned model is then used to predict WUE, T, and ET for each
#'
#' @param data data.frame with required columns: XX
#' @param config list with configuration options see \code{\link{tea_config}}
#'
#' @return see \code{\link{tea_predict}}
#' @export
#'
#' @examples
partition_tea <- function(data,  config = tea_config()) {
  tea_checkvars(data)
  dff <- tea_filter(data)
  rf <- tea_fit_wue(dff)
  datap <- tea_predict_wue(rf, data)
}

#' configure hyperparameters of the TEA algorithm
#'
#' @param CSWIlimit numeric scalar upper threshold of conservative surface
#'     wetness index for filter
#' @param perc
#'
#' @return
#' @export
#'
#' @examples
tea_config <- function(CSWIlimit = -0.5, perc = c(0.75)) {
  list(
    CSWIlimit = CSWIlimit
    ,perc = perc
  )
}



