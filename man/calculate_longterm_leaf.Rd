% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/priego.R
\name{calculate_longterm_leaf}
\alias{calculate_longterm_leaf}
\title{leaf-to-ambient CO2 ratio (chi_o) and Water use efficiency (WUE_o)}
\usage{
calculate_longterm_leaf(
  data,
  altitude,
  C = config$C,
  config = priego_config(),
  constants = etpart_constants()
)
}
\arguments{
\item{data}{Data.frame with columns
\itemize{
\item GPP: photosynthesis (umol CO2 m-2 s-1)
\item VPD: vapor pressure deficit (kPa).
\item Tair: air temperature (deg C).
}}

\item{altitude}{numeric value defining elevation (m).}

\item{C}{Empirical coeficient for C3 species.
(see Wang et al., 2017; Plant Nature)}

\item{config}{configuration (\code{\link{priego_config}})
with entry C for default of argument C}

\item{constants}{physical constants, see \code{\link{etpart_constants}}}
}
\value{
list with numeric entries:
\item{chi_o}{long-term effective "internal" leaf-to-ambient CO2 (unitless)}
\item{WUE_o}{long-term effective Water use efficiency (umolCO2 mmol-1)}
}
\description{
Calculate long-term effective "internal"
leaf-to-ambient CO2 (chi_o) and Water use efficiency (WUE_o)
}
\details{
the following metrics are calculated:

chi_o:

\deqn{logistic_chi_o = 0.0545*(Tair_g-25)-0.58*log(VPD_g)-0.0815*Z+C}
\deqn{chi_o <- exp(logistic_chi_o)/(1+exp(logistic_chi_o))}
WUE_o:

\deqn{WUE_o <- (390*(1-chi_o)*96)/(1.6*VPD_g)*0.001}

\code{Tair_g} and \code{VPD_g} are calculated based on the mean value of the
growing period.
The growing period is estimated as those periods over the 85 quantile of GPP.
}
\examples{
calculate_longterm_leaf(FIHyy, altitude=60)
}
\references{
Wang, H., I. C. Prentice, et al., (2017), Towards a universal model
for carbon dioxide uptake by plants, Nature Plants, 3(9), 734-741.

Reichstein, M., et al. (2005), On the separation of net ecosystem exchange
into assimilation and ecosystem respiration: review and improved algorithm,
Global Change Biology, 11(9), 1424-1439.
}
