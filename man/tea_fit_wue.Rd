% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tea.R
\name{tea_fit_wue}
\alias{tea_fit_wue}
\title{Fit a random forest model for predicting water use effciency (WUE)}
\usage{
tea_fit_wue(data_train, control)
}
\arguments{
\item{data_train}{trainign dataset where T = ET and hence WUE = GPP/ET}

\item{control}{list with entries rfseed. See \code{\link{tea_config}}.}
}
\value{
a tidymodels workflow object
}
\description{
Fit a random forest model for predicting water use effciency (WUE)
}
