% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tea.R
\name{compute_diurnal_centroid}
\alias{compute_diurnal_centroid}
\title{Diurnal centroid of sub-daily fluxes}
\usage{
compute_diurnal_centroid(flux, nrecday = 48)
}
\arguments{
\item{flux}{numeric vector of sub-daily flux that must be continuous
and regular of full days}

\item{nrecday}{integer: frequency of the sub-daily measurements,
48 for half hourly measurements}
}
\value{
The diurnal centroid, in hours at a daily frequency
}
\description{
Calculates the daily flux weighted time of a sub-daily flux.
}
