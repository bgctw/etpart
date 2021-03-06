% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tea.R
\name{partition_tea}
\alias{partition_tea}
\title{partition ET by TEA}
\usage{
partition_tea(data, control = tea_config())
}
\arguments{
\item{data}{data.frame with required columns: XX}

\item{control}{list with configuration options see \code{\link{tea_config}}}
}
\value{
see \code{\link{tea_predict}}, \code{data} with predictions
percentiles of WUE, E and T appended.
}
\description{
The Transpiration Estimation Algorithm (TEA) partitions ET by first
constraining the data to dry conditions where T/ET~1
(see \code{\link{tea_filter}}) and then learning water use efficiency (WUE)
relationship with predictors from the resulting training data
(see \code{link{tea_fit_wue}}).
The learned model is then used to predict WUE, T, and ET for each
}
