# ET partitioning according to
# Priego et al. 2018: Partitioning Eddy Covariance Water Flux Components Using
# Physiological and Micrometeorological Approaches


#' @title "internal" leaf-to-ambient CO2 (chi_o) and Water use efficiency (WUE_o)
#'
#' @description Calculation of long-term effective "internal"
#'    leaf-to-ambient CO2 (chi_o) and Water use efficiency (WUE_o)
#'
#' @param data      Data.frame or matrix containing all required variables.
#' @param ColPhotos column name of numeric vector containing time series of
#'   photosynthesis data (umol CO2 m-2 s-1)
#' @param ColVPD    column name of numeric vector containing time series of
#'   vapor pressure deficit (kPa).
#' @param ColTair   column name of numeric vector containing time series of
#'   air temperature (deg C).
#' @param Z         Z- numeric value defining elevation (km).
#' @param C         Empirical coeficient for C3 species.
#'   (see Wang et al., 2017; Plant Nature)

#' @export
#' @importFrom stats quantile
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
#' Tair_g and VPD_g are calculated based on the mean value of the growing period.
#' The growing period is estimated as those periods over the 85 quantile of Photos.
#'
#' @return list with numeric entries:
#' \item{chi_o}{long-term effective "internal" leaf-to-ambient CO2 (unitless)}
#' \item{WUE_o}{long-term effective "Water use efficiency (umolCO2 mmol-1)}
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
#' calculate_chi_o(EddySample)
calculate_longterm_leaf <- function(
  data
  ,ColPhotos = "GPP_NT_VUT_MEAN"
  ,ColVPD = "VPD_F"
  ,ColTair = "TA_F"
  ,C = 1.189
  ,Z = 0.27
) {
  iMissing <- which( !(c(ColPhotos, ColVPD, ColTair) %in% names(data)))
  if (length(iMissing)) stop(
    "Need to provide columns ",
    paste(c(ColPhotos, ColVPD, ColTair)[iMissing], collapse = ", "))
  dsagg <- data %>%
    select(Photos = ColPhotos, VPD = ColVPD, Tair = ColTair) %>%
    mutate(VPD_kPa = VPD/10) %>%  # Converting VPD units (hPa -> kPa)
    # Defining optimal growing period according to quantiles of photosynthesis
    filter(Photos >  quantile(Photos, probs = 0.85, na.rm=T)) %>%
    summarize(
      Tair_g = mean(Tair, na.rm = TRUE), VPD_g = mean(VPD_kPa, na.rm = TRUE))
  logistic_chi_o = 0.0545*(dsagg$Tair_g-25)-0.58*log(dsagg$VPD_g)-0.0815*Z+C
  chi_o <- exp(logistic_chi_o)/(1+exp(logistic_chi_o)) # Longterm effective CiCa
  ## 0.001 to convert umol/mol into umol/mmol
  WUE_o <- (390*(1-chi_o)*96)/(1.6*dsagg$VPD_g)*0.001
  list(chi_o = chi_o, WUE_o = WUE_o)
}



# ' # Find the optimal paramerers
#' @title  Model optimization routine
#'
#' @description Routine to estimate optimal parameters of a photosynthesis model
#'   using a multi-constraint Markov Chain Monte Carlo (MCMC)
#'   (see Perez-Priego et al., 2018).
#'
#' @param par_lower     A vector containing the lower bound of the parameters
#'   (a1,Do,To,beta)
#' @param par_upper     A vector containing the upper bound of the parameters
#'   (a1,Do,To,beta)
#' @param data          Data.frame or matrix containing all required variables.
#' @param ColPhotos     Column name of numeric vector containing time series of
#'   photosynthesis data (umol CO2 m-2 s-1).
#' @param ColPhotos_unc Column name of numeric vector containing time series of
#'   photosynthesis uncertainties (umol CO2 m-2 s-1).
#' @param ColH          Column name of numeric vector containing time series of
#'   sensible heat flux (W m-2).
#' @param ColVPD        Column name of numeric vector containing time series of
#'  vapor pressure deficit (hPa).
#' @param ColTair       Column name of numeric vector containing time series of
#'   air temperature (deg C).
#' @param ColPair       Column name of numeric vector containing time series of
#'   atmospheric pressure (kPa).
#' @param ColQ          Column name of numeric vector containing time series of
#'   photosynthetic active radiation (umol m-2 s-1).
#' @param ColCa         Column name of numeric vector containing time series of
#'   atmospheric CO2 concentration (umol Co2 mol air-1).
#' @param ColUstar      Column name of numeric vector containing time series of
#'   wind friction velocity (m s-1).
#' @param ColWS         Column name of numeric vector containing time series of
#'   wind velocity (m s-1).
#' @param ColSW_in      Column name of numeric vector containing time series of
#'   incoming short-wave radiation (W m-2).
#' @param Chi_o         Long-term effective chi
#' @param WUE_o         Long-term effective WUE
#'
#'
#' @export
#' @importFrom stats median
#' @importFrom FME Latinhyper
#' @importFrom  FME modMCMC
#'
#' @details the multi-objective function is defined as:
#'
#'
#'            \deqn{OF <- sum((photos-photosy_mod)/photos_unc)^2)/n + phi}
#'
#'
#'          where phi invokes optimality theory by minimizing the following term
#'
#'          \deqn{phi <- (sum(transpiration_mod)/sum(photos_mod)*WUE_o}
#'
#' @note The 4 model parameters (a1, Do, Topt and beta, see Perez-Priego 2018)
#'   are estimated using a multi-constraint Markov Chain Monte Carlo (MCMC).
#'   The objective function (OF) is to find those numerical solutions that
#'   minimize not only the mismatch between observed and modeled Photos but
#'   also the unit cost of transpiration by introducing a conditional factor
#'   demand (phi), which invokes the optimality hypothesis.
#'   The phi term is to be defined as the integrated cost of transpiration
#'   (i.e. transpiration_mod/photos_mod) over a time period (5 days) normalized
#'   by a factor describing the long-term effective water use efficiency (WUE_o).
#' @return a numeric vector containing 4 optimal parameters (a1,Do,To,beta):
#'  \item{a1}{radiation curvature}
#'  \item{D0}{empirical coef. related to response of stomatal closure to VPD.}
#'  \item{Topt}{optimum temperature}
#'  \item{beta}{A plant state variable defining the carbon cost of water.}
#'
#' @references
#' Perez-Priego, O., G. Katul, M. Reichstein et al. Partitioning eddy covariance
#' water flux components using physiological and micrometeorological approaches,
#' Journal of Geophysical Research: Biogeosciences. In press
#'
#' Reichstein, M., et al. (2005), On the separation of net ecosystem exchange
#' into assimilation and ecosystem respiration: review and improved algorithm,
#' Global Change Biology, 11(9), 1424-1439.
#'
#' @examples
#'  ## Selecting a single day (e.g. 15-05-2011)
#'  tmp <-  EddySample[ EddySample$TIMESTAMP_START>  201105150000,]
#'  tmp <-  tmp[tmp$TIMESTAMP_START<  201105160000,]
#'  ## Defining parameter values
#'
#'  optimal_parameters(par_lower = c(0, 0, 10, 0)
#'                     ,par_upper = c(400,0.4, 30, 1)
#'                    ,data = tmp
#'                    ,Chi_o = 0.88
#'                    ,WUE_o = 24.25)
optim_priego <- function(par_lower, par_upper, data, Chi_o, WUE_o
     ,ColPhotos = "GPP_NT_VUT_MEAN"
     ,ColPhotos_unc = "NEE_VUT_USTAR50_JOINTUNC"
     ,ColH = "H_F_MDS"
     ,ColVPD = "VPD_F"
     ,ColTair = "TA_F"
     ,ColPair = "PA_F"
     ,ColQ = "PPFD_IN"
     ,ColCa = "CO2_F_MDS"
     ,ColUstar = "USTAR"
     ,ColWS = "WS_F"
     ,ColSW_in = "SW_IN_F"
){
  dsf = data %>%
    rename(Photos = ColPhotos, Photos_unc = ColPhotos_unc, H = ColH, VPD = ColVPD,
           Tair = ColTair, Pair = ColPair, Q = ColQ, Q_in = ColSW_in, Ca = ColCa,
           Ustar = ColUstar, WS = ColWS) %>%
    filter(Photos>0 & Q>0 & Q_in>0) %>% # Rejecting bad data and filtering for daytime data
    mutate(
      VPD = VPD/10, # Convert VPD units (hPa -> kPa)
      Q_in = Q_in * 2, # convertion factor between W m2 to umol m-2 s-1
      #If PAR is not provided we use SW_in instead as an approximation of PAR
      Q = ifelse(is.na(Q), Q_in, Q),
      landa = (3147.5-2.37*(Tair+273.15))*1000 # Latent heat of evaporisation [J kg-1]
    )
  pars <- Latinhyper(cbind(par_lower,par_upper),num = 1)
  # try computing all that does not depend on parameters once outside cost
  ra <- compute_aerodynamic_conductance(dsf$WS, dsf$Ustar)
  constants = list(
    Cp = 1003.5, ##<< heat capacity [J kg-1 K-1].
    R_gas_constant = 0.287058, ##<< [J kg-1 deg K-1].
    M = 0.0289644 ##<< molar mass, [kg mol-1].
  )
  # Airdensity [kg m-3].
  dens = dsf$Pair/(constants$R_gas_constant*(dsf$Tair+273.15))
  Mden = dens/constants$M ##<< molar air density [mol m-3].
  Photos_max = quantile(dsf$Photos, probs=c(0.90), na.rm=T)
  Dmax = mean(dsf$VPD[dsf$Photos>Photos_max], na.rm=T)
  Photos_unc_threshold = ifelse(
    dsf$Photos*0.1 > dsf$Photos_unc, dsf$Photos*0.1, dsf$Photos_unc)
  min.RSS <- function(p) cost_optim_opriego(
    p, Chi_o, WUE_o,
    Photos = dsf$Photos,
    Photos_unc = dsf$Photos_unc,
    H = dsf$H,
    VPD = dsf$VPD,
    Tair = dsf$Tair,
    Pair = dsf$Pair,
    Q = dsf$Q,
    Ca = dsf$Ca,
    landa = dsf$landa,
    Ustar = dsf$Ustar,
    WS = dsf$WS,
    ra_c = ra$ra_c, ra_w = ra$ra_w,
    dens = dens, Mden = Mden, Photos_max = Photos_max, Dmax = Dmax,
    Photos_unc_threshold = Photos_unc_threshold
  )
  rss <- min.RSS(pars)
  parMCMC <- try(modMCMC(f = min.RSS
                         ,p = as.numeric(pars)
                         ,niter = 20000
                         ,updatecov = 500, ntrydr = 3
                         ,lower = par_lower , upper = par_upper
                         ,burninlength = 10000))
  out <- summary(parMCMC)
  paropt <- c(
    a1 = out[1,1],
    Do = out[1,2],
    Topt = out[1,3],
    beta = out[1,4]
  )
}

compute_aerodynamic_conductance <- function(WS,Ustar) {
  #--  Aerodynamic conductance
  ra_m <- WS/Ustar^2 ##<< aerodynamic resistance to momentum transfer.
  ra_b <- 6.2*Ustar^-0.67 ##<< aerodynamic resistance to heat transfer.
  ra <- ra_m+ra_b ##<< Monteith and Unsworth [2013]
  ra_w <- ra_m+2*(1.05/0.71/1.57)^(2/3)*ra_b # originally by Hicks et al., 1987.
  ra_c <- ra_m+2*(1.05/0.71)^(2/3)*ra_b
  list(ra_w = ra_w, ra_c = ra_c)
}

cost_optim_opriego <- function(par, Chi_o, WUE_o,
  Photos, Photos_unc, H, VPD, Tair, Pair, Q, Ca, landa, Ustar, WS,
  ra_w, ra_c,
  constants = list(
    Cp = 1003.5, ##<< heat capacity [J kg-1 K-1].
    R_gas_constant = 0.287058, ##<< [J kg-1 deg K-1].
    M = 0.0289644 ##<< molar mass, [kg mol-1].
  ),
  # may supply the arguments below for performance
  dens = Pair/(constants$R_gas_constant*(Tair+273.15)), # Air density [kg m-3].
  Mden = dens/constants$M, ##<< molar air density [mol m-3].
  Photos_max = quantile(Photos, probs=c(0.90), na.rm=T),
  Dmax = mean(VPD[Photos>Photos_max], na.rm=T),
  # We define the threshold of the uncertainties as the 1% of the magnitude of
  # bthe observation.
  Photos_unc_threshold = ifelse(Photos*0.1 > Photos_unc, Photos*0.1, Photos_unc)
  ) {
  # need correct column names and VPD in kPa, nonmissing Q in umol m-2 s-1
  # need filtered data
  # Chi_o is calculated using calculate_longerterm
  #-- Fitting parameters
  a1 <- par[1]
  D0 <- par[2]
  Topt <-  par[3]
  beta <- par[4]
  #--  plant temperature (Tplant) and plant to air vapor pressure deficit (VPD_plant)
  # Tplant <- (H*ra/(constants$Cp*dens))+Tair ##<< Approaximation of a surface temperature as canopy temperature (deg C)
  # es_plant <- 0.61078*exp((17.269*Tplant)/(237.3+ Tplant)) ##<< saturated vapor presure deficit at the plant surface.
  # es_air <- 0.61078*exp((17.269*Tair)/(237.3+ Tair)) ##<< saturated vapor presure deficit at the plant surface.
  # ea <- es_air-VPD ##<< atmospheric vapor pressure [kPa]
  # VPD_plant <- es_plant-ea ##<< Plant-to-air vapor pressure deficit [kPa]
  VPD_plant <- VPD
  #-- Defining optimum parameters
  Chimax <- Chi_o*(1/(1+beta*Dmax^0.5))
  # We assume that a max conductance is achieved at Photos_max
  # under optimum conditions.
  gcmax <- median(Photos_max/(Mden*Ca*(1-Chimax)), na.rm=T)
  #--  Calculating canopy stomatal conductance to CO2 [m s-1]
  gc_mod <- compute_stomatal_conductance_jarvis(par=par,Q,VPD,Tair,gcmax)
  gw_mod <- 1.6*gc_mod ## leaf canopy conductance to water vapor [m s-1]
  #--  Calculating "bulk" surface conductance
  gc_bulk <- Mden/(1/(gc_mod)+ra_c) # bulk canopy conductance CO2 [mol m-2 s-1]
  gw_bulk <- Mden/(1/(gw_mod)+ra_w) # bulk canopy conductance water[mol m-2 s-1]
  #-- Estimating Chi
  Chi <- Chi_o*(1/(1+beta*VPD_plant^0.5))
  #-- Estimating photosynthesis and transpiration rates
  Photos_mod <- gc_bulk*Ca*(1-Chi)
  transpiration_mod <- gw_bulk*VPD_plant/Pair*1000 ##<< [mmol H2O m-2 s-1]
  ##-- Objective function
  WaterCost_i <- sum(transpiration_mod, na.rm=T)/sum(Photos_mod, na.rm=T)
  Phi <- WaterCost_i*WUE_o
  FO <- sum(((Photos_mod-Photos)/Photos_unc_threshold)^2, na.rm=T)/length(Photos_mod)
  # FO <- sum((Photos_mod-Photos)^2, na.rm=T)/length(Photos_mod)
  ## Summed cost
  FO+Phi
}

#' @title Calculate stomatal conductance using Jarvis's approach.
#'
#' @description Calculate canopy stomatal conductance using Jarvis's approach.
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
#' @export
#' @details the following metrics are calculated:
#'
#'            \deqn{compute_stomatal_conductance_jarvis <- gcmax*f(Q)*f(VPD)*f(Tair)}
#'
#'
#' @return a numeric value:
#'         \item{compute_stomatal_conductance_jarvis}{canopy leaf conductance (untis refer to that given by gmax)}
#'
#'
#' @references Perez-Priego, O., G. Katul, M. Reichstein et al. Partitioning eddy covariance
#'             water flux components using physiological and micrometeorological approaches,
#'             Journal of Geophysical Research: Biogeosciences. In press
#'
#'
#' @examples
#'  ## Selecting a single day (e.g. 15-05-2011)
#'  tmp <-  EddySample[ EddySample$TIMESTAMP_START>  201105150000,]
#'  tmp <-  tmp[tmp$TIMESTAMP_START<  201105160000,]
#' compute_stomatal_conductance_jarvis(par=c(200, 0.2, 25)
#' ,Q = tmp$PPFD_IN * 2 #(convert to mumuol)
#' ,VPD = tmp$VPD_F
#' ,Tair = tmp$TA_F
#' ,gcmax = 1)
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




