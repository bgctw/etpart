#' Constants used throughout ETPart
#'
#' @param Cp heat capacity (J kg-1 K-1)
#' @param R_gas_constant (J kg-1 deg K-1)
#' @param M_air molar mass of air (kg mol-1)
#'
#' @return a list with above arguments as entries.
#' @export
etpart_constants <- function(
  Cp = 1003.5,
  R_gas_constant = 0.287058,
  M_air = 0.0289644
) {
  list(
    Cp = Cp,
    R_gas_constant = R_gas_constant,
    M_air = M_air
  )
}
