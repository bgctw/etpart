# ET partitioning according to
# Priego et al. 2018: Partitioning Eddy Covariance Water Flux Components Using
# Physiological and Micrometeorological Approaches

#' partition ET by Priego approach
#'
#' Model transpiration by optimizing short-term adaptation of ratio of leaf
#' internal to ambiant CO2 concentration.
#'   (see Perez-Priego et al., 2018).
#'
#' @param data          Data.frame with variables;
#' \itemize{
#'    \item those required by \code{\link{predict_transpiration_opriego}} and
#'    \item Rg: incoming short-wave radiation (W m-2).
#'    }
#' @param ...   further arguments passed to both
#'    \code{\link{calculate_longterm_leaf}}
#'    and  \code{\link{predict_transpiration_opriego}}
#'    such as \code{config=\link{priego_config}()} and
#'    \code{constants=\link{etpart_constants}()}
#'
#' @details
#' For each five-day window the short-term variation in leaf-internal
#' water to CO2 ratio (chi) and the modifiers for photosynthesis fitted so to
#'   minimize not only the mismatch between observed and modeled GPP but
#'   also the unit cost of transpiration.
#'   The transpiration cost is specified  by introducing a conditional factor
#'   demand (phi), which invokes the optimality hypothesis.
#'   The phi term is to be defined as the integrated cost of transpiration
#'   (i.e. transpiration_mod/photos_mod) over a time period (5 days) normalized
#'   by a factor describing the long-term effective water use efficiency (WUE_o).
#'
#' These parameters are then used to predict transpiration.
#'
#' Evaporation is computed as difference to ET from eddy-covariance.
#'
#' The multi-objective function is defined as:
#'
#'  \deqn{OF = sum((photos-photos_{mod})/photos_{unc})\^2)/n + phi}
#'
#'  where phi invokes optimality theory by minimizing the following term
#'
#'  \deqn{phi = (sum(transpiration_mod)/sum(photos_mod)*WUE_o}
#'
#' The 4 model parameters (a1, Do, Topt and beta, see Perez-Priego 2018) are
#'  \describe{
#'  \item{a1}{radiation curvature}
#'  \item{D0}{empirical coef. related to response of stomatal closure to VPD.}
#'  \item{Topt}{optimum temperature}
#'  \item{beta}{A plant state variable defining the carbon cost of water.}
#'   }
#'
#' @seealso
#' \code{\link{calculate_longterm_leaf}},
#' \code{\link{predict_transpiration_opriego}},
#' \code{\link{priego_config}},
#' \code{\link{etpart_constants}}
#'
#' @return data with updated or added columns T and E in mm/hour.
#'
#' @references
#' Perez-Priego, O., G. Katul, M. Reichstein et al. Partitioning eddy covariance
#' water flux components using physiological and micrometeorological approaches,
#' Journal of Geophysical Research: Biogeosciences. In press
#'
#' Reichstein, M., et al. (2005), On the separation of net ecosystem exchange
#' into assimilation and ecosystem respiration: review and improved algorithm,
#' Global Change Biology, 11(9), 1424-1439.
#' @export
partition_priego <- function(data, ...) {
  lt <- calculate_longterm_leaf(data, ...)
  dfT <- estimate_T_priego_5days(data, lt$chi_o, lt$WUE_o, ...)
  dfT %>% mutate(E = ET - T)
}

#' Configure parameters of the Priego transpiration partitioning
#'
#' @param par_lower     A vector containing the lower bound of the parameters
#'   (a1,Do,To,beta)
#' @param par_upper     A vector containing the upper bound of the parameters
#'   (a1,Do,To,beta)
#' @param C Empirical coeficient for C3 species (see Wang et al., 2017; Plant Nature).
#' @param niter number of iterations for the MCMC
#' @param updatecov number of iterations after which the parameter covariance
#'    matrix is (re)evaluated based on the parameters kept thus far, and used
#'    to update the MCMC jumps.
#' @param ntrydr maximal number of tries for the delayed rejection procedure
#' @param burninlength number of initial iterations to be removed from output.
#'
#' @return a list with above arguments as entries.
#' @seealso \code{\link{partition_priego}}
#' @export
priego_config <- function(
  C = 1.189,
  par_lower = c(a1=0, Do=0, To=10, beta=0),
  par_upper = c(a1=400, Do=0.4, To=30, beta=1),
  #niter = 20000,
  niter = 2000,
  updatecov = 200,
  ntrydr = 1,
  burninlength = 600
) {
  list(
    C = C,
    par_lower = par_lower,
    par_upper = par_upper,
    niter = niter,
    updatecov = updatecov,
    ntrydr = ntrydr,
    burninlength = burninlength
  )
}

#' leaf-to-ambient CO2 ratio (chi_o) and Water use efficiency (WUE_o)
#'
#' Calculate long-term effective "internal"
#'    leaf-to-ambient CO2 (chi_o) and Water use efficiency (WUE_o)
#'
#' @param data      Data.frame with columns
#' \itemize{
#'    \item GPP: photosynthesis (umol CO2 m-2 s-1)
#'    \item VPD: vapor pressure deficit (kPa).
#'    \item Tair: air temperature (deg C).
#'    }
#' @param Z         Z- numeric value defining elevation (km).
#' @param C         Empirical coeficient for C3 species.
#'   (see Wang et al., 2017; Plant Nature)
#' @param constants physical constants, see \code{\link{etpart_constants}}
#'
#' @details the following metrics are calculated:
#'
#' chi_o:
#'
#'   \deqn{logistic_chi_o = 0.0545*(Tair_g-25)-0.58*log(VPD_g)-0.0815*Z+C}
#'   \deqn{chi_o <- exp(logistic_chi_o)/(1+exp(logistic_chi_o))}
#' WUE_o:
#'
#'   \deqn{WUE_o <- (390*(1-chi_o)*96)/(1.6*VPD_g)*0.001}
#'
#' \code{Tair_g} and \code{VPD_g} are calculated based on the mean value of the growing period.
#' The growing period is estimated as those periods over the 85 quantile of GPP.
#'
#' @return list with numeric entries:
#' \item{chi_o}{long-term effective "internal" leaf-to-ambient CO2 (unitless)}
#' \item{WUE_o}{long-term effective Water use efficiency (umolCO2 mmol-1)}
#'
#' @references
#' Wang, H., I. C. Prentice, et al., (2017), Towards a universal model
#' for carbon dioxide uptake by plants, Nature Plants, 3(9), 734-741.
#'
#' Reichstein, M., et al. (2005), On the separation of net ecosystem exchange
#' into assimilation and ecosystem respiration: review and improved algorithm,
#' Global Change Biology, 11(9), 1424-1439.
#'
#' @examples
#' calculate_longterm_leaf(FIHyy, Z=0.060)
#' @export
#' @importFrom stats quantile
calculate_longterm_leaf <- function(
  data
  ,altitude
  ,C = config$C
  ,config = priego_config()
  ,constants = etpart_constants()
) {
  check_required_cols(data, c("GPP","VPD","Tair"))
  dsagg <- data %>%
    mutate(VPD_kPa = VPD/10) %>%  # Converting VPD units (hPa -> kPa)
    # Defining optimal growing period according to quantiles of photosynthesis
    filter(GPP >  quantile(GPP, probs = 0.85, na.rm=T)) %>%
    summarize(
      Tair = mean(Tair, na.rm = TRUE), VPD = mean(VPD_kPa, na.rm = TRUE))
  logistic_chi_o = 0.0545*(dsagg$Tair-25)-0.58*log(dsagg$VPD)-0.0815*altitude+C
  chi_o <- exp(logistic_chi_o)/(1+exp(logistic_chi_o)) # Longterm effective CiCa
  ## 0.001 to convert umol/mol into umol/mmol
  WUE_o <- (390*(1-chi_o)*96)/(1.6*dsagg$VPD)*0.001
  list(chi_o = chi_o, WUE_o = WUE_o)
}

estimate_T_priego_5days <- function(
  data, iday, chi_o, WUE_o, config = priego_config(), ...
) {
  ds5 <- filter(data, between(cumday, iday - 2, iday + 2))
  popt <- ds5 %>%
    # better do inside optim_priego
    # mutate(
    #   GPP_sd = ifelse(is.na(.data$GPP_sd), .data$GPP*0.1,.data$GPP_sd),
    #   Q = ifelse(is.na(.data$Q)==TRUE, .data$Rg*2,  .data$Q)
    # ) %>%
    optim_priego(chi_o, WUE_o, config=config, ...)
  ans <- ds5 %>%
    filter(cumday == iday) %>%
    mutate(T = predict_transpiration_opriego(.data, popt$paropt, chi_o, WUE_o, ...))
}


#' @importFrom stats median
#' @importFrom FME Latinhyper
#' @importFrom  FME modMCMC
optim_priego <- function(data, chi_o, WUE_o
     ,config = priego_config(), constants = etpart_constants()
){
  dsf = data %>%
    # Rejecting bad data and filtering for daytime data
    filter(!isnight & GPP>0 & Q>0 & Rg>0) %>%
    mutate(
      VPD_kPa = VPD/10, # Convert VPD units (hPa -> kPa)
      #If PAR is not provided we use SW_in instead as an approximation of PAR
      # *2: conversion factor between W m2 to umol m-2 s-1
      Q = ifelse(is.na(Q), Rg*2, Q),
      #landa = (3147.5-2.37*(Tair+273.15))*1000 # Latent heat of evaporation [J kg-1]
      GPP_sd = ifelse(is.na(.data$GPP_sd), .data$GPP*0.1, .data$GPP_sd),
    )
  pars <- Latinhyper(cbind(config$par_lower,config$par_upper),num = 1)
  # try computing all that does not depend on parameters once outside cost
  ra <- compute_aerodynamic_conductance(dsf$u, dsf$ustar)
  dens = calculate_air_density(dsf$Pair, dsf$Tair, constants) # [kg m-3].
  Mden = dens/constants$M_air ##<< molar air density [mol m-3].
  GPP_max = quantile(dsf$GPP, probs=c(0.90), na.rm=T)
  Dmax = mean(dsf$VPD_kPa[dsf$GPP>GPP_max], na.rm=T)
  GPP_sd_threshold = pmax(dsf$GPP*0.1, dsf$GPP_sd)
  # VPD_plant = estimate_VPD_plant(
  #   H=H, Tair=Tair, Pair=Pair, VPD=VPD_kPa, ra=ra, constants=constants, dens=dens)
  VPD_plant = dsf$VPD_kPa
  min.RSS <- function(p) cost_optim_opriego(
    p, chi_o, WUE_o,
    GPP = dsf$GPP,
    GPP_sd = dsf$GPP_sd,
    H = dsf$H,
    VPD = dsf$VPD_kPa,
    Tair = dsf$Tair,
    Pair = dsf$Pair,
    Q = dsf$Q,
    Ca = dsf$Ca,
    ustar = dsf$ustar,
    u = dsf$u,
    constants = NULL, # actually never used because of precomputations
    ra = ra, VPD_plant = VPD_plant,
    Mden = Mden, GPP_max = GPP_max, Dmax = Dmax,
    GPP_sd_threshold = GPP_sd_threshold
  )
  rss <- min.RSS(pars)
  parMCMC <- try(modMCMC(f = min.RSS
                         ,p = as.numeric(pars)
                         ,niter = config$niter
                         ,updatecov = config$updatecov, ntrydr = config$ntrydr
                         ,lower = config$par_lower , upper = config$par_upper
                         ,burninlength = config$burninlength))
  out <- summary(parMCMC)
  paropt <- c(
    a1 = out[1,1],
    Do = out[1,2],
    Topt = out[1,3],
    beta = out[1,4]
  )
  list(paropt = paropt, parMCMC = parMCMC)
}

calculate_air_density <- function(Pair, Tair, constants=etpart_constants()) {
  dens = Pair/(constants$R_gas_constant*(Tair+273.15))
}

compute_aerodynamic_conductance <- function(u,ustar) {
  #--  Aerodynamic conductance
  ra_m <- u/ustar^2 ##<< aerodynamic resistance to momentum transfer.
  ra_b <- 6.2*ustar^-0.67 ##<< aerodynamic resistance to heat transfer.
  ra <- ra_m+ra_b ##<< Monteith and Unsworth [2013]
  ra_w <- ra_m+2*(1.05/0.71/1.57)^(2/3)*ra_b # originally by Hicks et al., 1987.
  ra_c <- ra_m+2*(1.05/0.71)^(2/3)*ra_b
  list(ra_w = ra_w, ra_c = ra_c, ra = ra)
}

cost_optim_opriego <- function(par, chi_o, WUE_o,
  GPP, GPP_sd, H, VPD, Tair, Pair, Q, Ca, ustar, u,
  # need supply the arguments below for performance, see optim_priego
  constants,
  ra, VPD_plant,
  Mden, GPP_max, Dmax, GPP_sd_threshold
  ) {
  # note, that VPD and VPD_Plant need to be specified in kPa (not hPa)
  pred <- predict_T_GPP_opriego(
    par, chi_o, WUE_o,
    GPP, GPP_sd, H, VPD, Tair, Pair, Q, Ca, ustar, u,
    constants,
    estimate_VPD_eddy,
    ra, dens=NULL, VPD_plant,
    Mden, GPP_max, Dmax, GPP_sd_threshold
  )
  WaterCost_i <- sum(pred$T, na.rm=T)/sum(pred$GPP, na.rm=T)
  Phi <- WaterCost_i*WUE_o
  FO <- sum(((pred$GPP-GPP)/GPP_sd_threshold)^2, na.rm=T)/length(GPP)
  FO+Phi
}

predict_T_GPP_opriego <- function(
  par, chi_o, WUE_o,
  GPP, GPP_sd, H, VPD, Tair, Pair, Q, Ca, ustar, u,
  # the following intermediates are independent of par and maybe precomputed
  constants = etpart_constants(),
  festimate_VPD = estimate_VPD_plant,
  ra = compute_aerodynamic_conductance(u, ustar),
  dens = calculate_air_density(Pair, Tair, constants), # [kg m-3].
  VPD_plant = festimate_VPD(
    H=H, Tair=Tair, Pair=Pair, VPD=VPD, ra=ra, constants=constants, dens=dens),
  Mden = dens/constants$M, ##<< molar air density [mol m-3].
  GPP_max = quantile(GPP, probs=c(0.90), na.rm=T),
  Dmax = mean(VPD[GPP>GPP_max], na.rm=T),
  GPP_sd_threshold = ifelse(
    GPP*0.1 > GPP, GPP*0.1, GPP_sd)
) {
  beta <- par[4]
  g_bulk <- estimate_canopy_conductances(
    par, ra, chi_o, Ca,
    Tair=Tair, VPD=VPD, Q=Q, Dmax=Dmax, GPP_max=GPP_max, Mden=Mden)
  chi = chi_o*(1/(1+beta*VPD_plant^0.5))
  transpiration_mod <- g_bulk$gw_bulk*VPD_plant/Pair*1000 ##<< [mmol H2O m-2 s-1]
  GPP_mod <- g_bulk$gc_bulk*Ca*(1-chi)
  # T in mm/m2/s GPP in GPP_mod in ?umol CO2 m-2 s-1
  list(T = transpiration_mod, GPP = GPP_mod)
}

estimate_canopy_conductances <- function(
  par, ra, chi_o, Ca, Tair, VPD, Q, Dmax, GPP_max, Mden) {
  a1 <- par[1]
  D0 <- par[2]
  Topt <-  par[3]
  beta <- par[4]
  #-- Defining optimum parameters
  Chimax <- chi_o*(1/(1+beta*Dmax^0.5))
  # We assume that a max conductance is achieved at GPP_max
  # under optimum conditions.
  gcmax <- median(GPP_max/(Mden*Ca*(1-Chimax)), na.rm=T)
  #--  Calculating canopy stomatal conductance to CO2 [m s-1]
  gc_mod <- compute_stomatal_conductance_jarvis(par=par,Q,VPD,Tair,gcmax)
  gw_mod <- 1.6*gc_mod ## leaf canopy conductance to water vapor [m s-1]
  #--  Calculating "bulk" surface conductance
  gc_bulk <- Mden/(1/(gc_mod)+ra$ra_c) # bulk canopy conductance CO2 [mol m-2 s-1]
  gw_bulk <- Mden/(1/(gw_mod)+ra$ra_w) # bulk canopy conductance water[mol m-2 s-1]
  list(gc_bulk = gc_bulk, gw_bulk = gw_bulk)
}


estimate_VPD_plant <- function(
  H, Tair, Pair, VPD, ra, constants,
  dens = Pair/(constants$R_gas_constant*(Tair+273.15)) # Air density [kg m-3].
) {
  #--  plant temperature (Tplant) and plant to air vapor pressure deficit (VPD_plant)
  # Approximation of a surface temperature as canopy temperature (deg C)
  Tplant <- (H*ra$ra/(constants$Cp*dens))+Tair
  # saturated vapor pressure deficit at the plant surface.
  es_plant <- 0.61078*exp((17.269*Tplant)/(237.3+ Tplant))
  # saturated vapor pressure deficit at the plant surface.
  es_air <- 0.61078*exp((17.269*Tair)/(237.3+ Tair))
  # atmospheric vapor pressure [kPa]
  ea <- es_air-VPD
  # Plant-to-air vapor pressure deficit [kPa]
  VPD_plant <- es_plant-ea
}
estimate_VPD_eddy <- function(
  H, Tair, Pair, VPD, ra, constants,
  dens = Pair/(constants$R_gas_constant*(Tair+273.15)) # Air density [kg m-3].
) { VPD }


#' Canopy stomatal conductance by Jarvis.
#'
#' Calculate canopy stomatal conductance using Jarvis's approach.
#'
#' @param par       Set of parameter for the respective sensitivity function.
#' @param Q         Vector containing time series of
#'   photosynthetic active radiation (umol CO2 m-2 s-1)
#' @param VPD       Vector containing time series of
#'   vapor pressure deficit (kPa).
#' @param Tair      Vector containing time series of
#'   air temperature (deg C).
#' @param gcmax     Empirical parameter defining
#'   maximum conductance (m s-1 or mol m-2 s-1).
#'
#' @details the following metrics are calculated:
#'
#'\deqn{compute_stomatal_conductance_jarvis <- gcmax*f(Q)*f(VPD)*f(Tair)}
#'
#'
#' @return a numeric value:
#' canopy leaf conductance (units refer to that given by gmax)
#'
#' @references
#' Perez-Priego, O., G. Katul, M. Reichstein et al. Partitioning eddy covariance
#' water flux components using physiological and micrometeorological approaches,
#' Journal of Geophysical Research: Biogeosciences. In press
#'
#' @seealso
#'   \code{\link{partition_priego}},
#'
#' @examples
#'  ## Selecting a single day (e.g. 15-05-2010)
#' tmp <- FIHyy[(FIHyy$timestamp > ISOdatetime(2010,5,15,0,0,0,tz=tz)) &
#'              (FIHyy$timestamp <= ISOdatetime(2010,5,16,0,0,0,tz=tz)),]
#' gc = compute_stomatal_conductance_jarvis(par=c(200, 0.2, 25)
#' ,Q = tmp$Q
#' ,VPD = tmp$VPD/10 # convert hPa to kPa
#' ,Tair = tmp$Tair
#' ,gcmax = 1)
#' plot(gc ~ tmp$timestamp)
#' @export
compute_stomatal_conductance_jarvis <- function(par, Q, VPD, Tair, gcmax) {
  #-- Fitting parameters
  a1 <- par[1] ##<< radiation curvature
  D0 <- par[2] # empirical coef. related to response of stomatal closure to VPD.
  Topt <-  par[3] ##<< optimum temperature
  # light response curve
  FQ <- (Q)/(Q+a1) ##<<  See Baldocchi et al., 1991 AFM and Jarvis 1976
  # VPD sensitivity
  Fd <- exp(-D0*VPD)
  # Optimum Temperature curve
  Tl <- 0 ##<<  minimum temperature [deg C]
  Th <- 50 ##<<  max temperature [deg C]
  b4 <- (Th-Topt)/(Th-Tl) ##<< parameter
  b3 <- 1/((Topt-Tl)*(Th-Topt)^b4) ##<< parameter
  Ftemp <- b3*(Tair-Tl)*(Th-Tair)^b4
  ## Calculating conductance
  sens_fun <- FQ*Fd*Ftemp
  # scaling between 0 and 1.
  sensitivity_function_scaled <- sens_fun/max(sens_fun, na.rm=T)
  gcmax*sensitivity_function_scaled
}

#' Compute transpiration given parameters
#'
#' @param par optimized parameters see \code{\link{partition_priego}}
#' @param data          Data.frame with columns
#' \itemize{
#'   \item GPP:     photosynthesis data (umol CO2 m-2 s-1).
#'   \item GPP_sd:  photosynthesis uncertainties (umol CO2 m-2 s-1).
#'   \item H:       sensible heat flux (W m-2).
#'   \item VPD:     vapor pressure deficit (hPa).
#'   \item Tair:    air temperature (deg C).
#'   \item Pair:    atmospheric pressure (kPa).
#'   \item Q:       photosynthetic active radiation (umol m-2 s-1).
#'   \item Ca:      atmospheric CO2 concentration (umol Co2 mol air-1).
#'   \item ustar:   wind friction velocity (m s-1).
#'   \item u:       wind velocity (m s-1).
#'   }
#' @param chi_o         Long-term effective chi
#' @param WUE_o         Long-term effective WUE
#' @param ...           Further arguments to \code{predict_T_GPP_opriego}
#'    such as \code{constants} (see \code{\link{etpart_constants}})
#'    or a non-default \code{VPD_plant} or
#'    precomputed intermediates.
#' @return vector of estimated transpiration for each record in data
#'   in mm/hour (kg m-2 hour-1)
#' @seealso
#'   \code{\link{partition_priego}},
#'   \code{\link{compute_stomatal_conductance_jarvis}}
#' @export
predict_transpiration_opriego <- function(
  data, par, chi_o, WUE_o, ...
) {
  check_required_cols(
    data, c("GPP","GPP_sd","H","VPD","Tair","Pair","Q","Ca","ustar","u"))
  pred <- predict_T_GPP_opriego(
    par, chi_o, WUE_o,
    GPP=data$GPP, GPP_sd=data$GPP_sd, H=data$H,
    VPD=data$VPD/10, Tair=data$Tair, Pair=data$Pair,
    Q=data$Q, Ca=data$Ca, ustar=data$ustar, u=data$u,
    ...
  )
  pred$T * (18.01528/1e6)*3600  # from mmol m-2 s-2 to mm per hour
}




