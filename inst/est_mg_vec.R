## script contains functions that enable retrieval of soil moisture
## from brightness temperatures
## fresnel eqs: reflectivity
## mironov eqs: real and imaginary dielectric
## and additional helper functions for the larger functions
#get predicted brightness temperature

# load the colliander version of tau omega
source(here("R", "colliandertb.R"))

# Helpers -----------------------------------------------------------------

tau_omega <- function(r, rhfac, gamma, Tair, Tsoil = NULL, omega, twotemps = F, split = F, noT = F, rtms = F) {
  if (twotemps == T & is.null(Tsoil)) {
    stop("Tsoil must be provided when twotemps=T; see tau_omega for more")
  }

  if (split) {
    tbs <- (1 - rhfac * r) * gamma * Tsoil
    tbc <- (1 - omega) * (1 - gamma) * (1 + ((rhfac * r) * gamma)) * Tair # am I missing the double bounce here?
    # why is there only one roughness?
    # gamma is the vegetation transmissivity
    return(list(tb = tbs + tbc, tbs = tbs, tbc = tbc))

  } else if (twotemps == T) {

    return((1 - rhfac * r) * gamma * Tsoil + (1 - omega) * (1 - gamma) * (1 + ((rhfac * r) * gamma)) * Tair) # plot r and gamma

  } else if (rtms) {

    tb <- RTMS_tb(r = r, rhfac = rhfac, gamma = gamma, Tsoil = Tsoil, Tveg = Tair, omega = omega)[[4]]
    return(tb)

  } else if (noT == TRUE) {
    r_rough <- rhfac * r
    return((1 - r_rough) * gamma + (1 - omega) * (1 - gamma) * (1 + (r_rough * gamma)))

  } else {

    return((1 - rhfac * r) * gamma * Tair + (1 - omega) * (1 - gamma) * (1 + ((rhfac * r) * gamma)) * Tair)
  }
}


mironov <- function(f, # frequency
                    mv, # soil moisture volume
                    cf) { # clay fraction
  eps_0 = 8.854 * 10 ^ (-12)
  eps_winf = 4.9
  fHz = f

  #Initializing the GRMDM spectroscopic parameters with clay fraction
  #RI & NAC of dry soils
  znd = 1.634 - 0.539 * cf + 0.2748 * cf ^ 2
  zkd = 0.03952 - 0.04038 * cf

  #Maximum bound water fraction
  zxmvt = 0.02863 + 0.30673 * cf
  #zxmvt = 0.02863 + 0.0030673 * cf

  #Bound water parameters
  zep0b = 79.8 - 85.4 * cf + 32.7 * cf ^ 2
  ztaub = 1.062 * 10 ^ (-11) + 3.450 * 10 ^ (-12) * cf
  zsigmab = 0.3112 + 0.467 * cf

  #Unbound (free) water parameters
  zep0u = 100
  ztauu = 8.5 * 10 ^ (-12)
  zsigmau = 0.3631 + 1.217 * cf

  #Computation of epsilon water (bound & unbound)
  zcxb = (zep0b - eps_winf) / (1 + (2 * pi * fHz * ztaub) ^ 2)
  zepwbx = eps_winf + zcxb
  zepwby = zcxb * (2 * pi * fHz * ztaub) + zsigmab / (2 * pi * eps_0 * fHz)
  zcxu = (zep0u - eps_winf) / (1 + (2 * pi * fHz * ztauu) ^ 2)
  zepwux = eps_winf + zcxu
  zepwuy = zcxu * (2 * pi * fHz * ztauu) + zsigmau / (2 * pi * eps_0 * fHz)

  #Computation of refractive index of water (bound & unbound)
  znb = sqrt(sqrt(zepwbx ^ 2 + zepwby ^ 2) + zepwbx) / sqrt(2)
  zkb = sqrt(sqrt(zepwbx ^ 2 + zepwby ^ 2) - zepwbx) / sqrt(2)
  znu = sqrt(sqrt(zepwux ^ 2 + zepwuy ^ 2) + zepwux) / sqrt(2)
  zku = sqrt(sqrt(zepwux ^ 2 + zepwuy ^ 2) - zepwux) / sqrt(2)

  #Computation of soil refractive index (nm & km): xmv can be a vector
  zxmvt2 = min(mv, zxmvt)
  zflag = mv >= zxmvt
  znm = znd + (znb - 1) * zxmvt2 + (znu - 1) * (mv - zxmvt) * zflag
  zkm = zkd + zkb * zxmvt2 + zku * (mv - zxmvt) * zflag

  #Computation of soil dielectric constant:
  zepmx = znm ^ 2 - zkm ^ 2
  zepmy = znm * zkm * 2
  #eps = zepmx + i*zepmy
  return(list(zepmx, zepmy))

}

fresnelr <- function(eps, theta){

  # theta is converted to radians then cos(theta) and sin2 theta are calculated
  thetar = theta*pi/180
  costheta = cos(thetar)
  sin2theta = sin(thetar)^2

  # calculate fresnel h and v polarization reflectivity
  fH = abs((costheta - sqrt(eps - sin2theta)) / (costheta + sqrt(eps - sin2theta)))^2
  fV = abs((eps*costheta - sqrt(eps - sin2theta)) / (eps*costheta + sqrt(eps - sin2theta)))^2

  return(list(fH=fH,
              fV=fV))
}

mpdi_transfer <- function(TbV, TbH, fH, fV, omega, inc_angle = 40){
  # the equations here are from Meesters et al., 2005 IEEE Geoscience and remote sensing letters
  # variable names they use:
  # emissivities (eH and eV) are functions of the soil dielectric constant (k) and the incidence angle (u)
  # replaced the emissivity terms to make more sense with the other functions.
  eH <- 1-fH
  eV <- 1-fV
  theta <- inc_angle * pi / 180
  MPDI <- (TbV - TbH)/(TbV + TbH) # equation 5

  d <- (1/2) * (omega / (1-omega)) # equation 9

  a_ku <- (1/2) * (((eV-eH)/MPDI) - (eV-eH)) # equation 8

  ad <- a_ku * d # making equation 13 easier to troubleshoot
  tau <- cos(theta) * log(ad + sqrt(ad^2 + a_ku + 1)) # equation 13
  gamma <- exp(-1 * tau/cos(theta)) # equation 3

  return(gamma)
}



# predictions  ------------------------------------------------------------
pred_H_V_tb<- function(Tb_H, Tb_V,#measured brightness temps
                         fH, fV, # Fresnel reflectivity (h and v pol)
                         gamma,
                         Tair, Tsoil,
                         omega,
                         rhfac,
                        twotemps, split = F, rtms=F){

  #take a vector or single value of possible reflectivity values over a realistic range of eps
  #predict the brightness temperatures for the array

  # predicted h pol brightness temperature
  pred_T_H = tau_omega(r = fH,
                       rhfac = rhfac,
                       gamma = gamma,
                       Tair = Tair, Tsoil= Tsoil,
                       omega = omega, twotemps=twotemps, split = split, rtms=rtms)
  # predicted v pol brightness temperature
  pred_T_V = tau_omega(r = fV,
                       rhfac = rhfac,
                       gamma = gamma,
                       Tair = Tair, Tsoil= Tsoil,
                       omega = omega, twotemps=twotemps, split = split, rtms=rtms)

  if(split == T){
    pred_T_H = pred_T_H$tb; pred_T_V = pred_T_V$tb

    tbs_h = pred_T_H$tbs; tbc_h = pred_T_H$tbc
    tbs_v = pred_T_V$tbs; tbc_v = pred_T_V$tbc

    #compute the squared differences between predicted and observed Tb
    costfun =(pred_T_H - Tb_H)^2 + (pred_T_V - Tb_V)^2 #(pred_T_H - Tb_H)^2 +
    tb_output <- list(cf = costfun, # return the cost function value
                      predTBH = pred_T_H, # return brightness temperatures
                      predTBV = pred_T_V,
                      tb_soil_h = tbs_h, tb_can_h = tbc_h,
                      tb_soil_v = tbs_v, tb_can_v = tbc_v)
  }else{

  #compute the squared differences between predicted and observed Tb
  tbh_res = (pred_T_H - Tb_H)^2
  tbv_res = (pred_T_V - Tb_V)^2
  costfun = tbh_res + tbv_res
  #microwave polarization difference index
  mpdi = (pred_T_V - pred_T_H)/(pred_T_V + pred_T_H)

  tb_output <- list(cf = costfun, # return the cost function value
                 predTBH = pred_T_H, # return brightness temperatures
                 predTBV = pred_T_V,
                 tbh_cf = tbh_res,
                 tbv_cf = tbv_res,
                 mpdi = mpdi,
                 mpdi_cf = mpdi+costfun)
  }

  return(tb_output)

}


# retrieval of sm and vod -------------------------------------------------
# TbH = tb_nowet_interval$tbh[1]
# TbV = tb_nowet_interval$tbv[1]
# Tair = tb_nowet_interval$phys_temp[1]
# Tsoil = tb_nowet_interval$mean_soil_temp[1]
# omega = omegaInput
# roughness = roughnessInput
# vod = seq(0, 3, by=0.001)
# twotemps = T
# split=F
# mpdi=T
# smc <- mean_smc[1]

get_sm_vod <- function(smc,
                      TbH, TbV, Tair, Tsoil = NULL,
                      vod = seq(0, 3, length.out=10),
                      omega = 0.05, inc_angle = 40,
                      cf = 0.232, roughness = 0.16,# lambda = 0, cf.out = F,
                      twotemps=F, split=F, rtms=FALSE, mpdi=FALSE){

  ##EPSILON##
  # first, create vectors to hold the epsilon values calculated from mironov
  eps_real = 0 * smc
  eps_im = 0 * smc

  # for each index of epsilon determine the real and imaginary component of
  # epsilon
  for (i in seq_along(eps_real)){
    eps_i = mironov(1.4e9, smc[i], cf)
    eps_real[i] = eps_i[[1]]
    eps_im[i] = eps_i[[2]] * -1

  }

  # epsilon is a complex number with real and imaginary components
  eps = complex(real = eps_real, imaginary = eps_im)

  # theta is converted to radians then cos(theta) and sin2 theta are calculated
  theta = inc_angle * pi / 180
  mycostheta = cos(theta)
  mysin2theta = sin(theta)^2

  # calculate fresnel h and v polarization reflectivity
  fH = abs((mycostheta - sqrt(eps - mysin2theta)) / (mycostheta + sqrt(eps - mysin2theta)))^2
  fV = abs((eps*mycostheta - sqrt(eps - mysin2theta)) / (eps*mycostheta + sqrt(eps - mysin2theta)))^2

  ##VOD AND GAMMA
  vod = vod # input vector of potential VOD values (potentially unnecessary)

  if(mpdi) {
    print("HI")
    gamma = mpdi_transfer(TbV, TbH, fH, fV, omega = omega, inc_angle = inc_angle)

  }else if(rtms){
    gamma=exp(-1 * vod)# calculate gamma

  }else{
    gamma=exp(-1 * vod / mycostheta)# calculate gamma
  }

  omega = omega # omega from tau-omega
  rhfac = exp(-roughness*mycostheta)#^2roughness factor applied to reflectivities

  ##ESTIMATING COST FUNCTION FOR ALL GAMMA/VOD AND SOIL MOISTURE VALUES##
  eps_baseline <- vector() # vector to store cost function values

  # matrix to hold cost functions values for each epsilon and gamma
  out_cf <- matrix(nrow=length(eps), ncol=length(gamma))
  out_tbh <- matrix(nrow=length(eps), ncol=length(gamma))
  out_tbv <- matrix(nrow=length(eps), ncol=length(gamma))

  # return residuals
  out_tbh_cf <- matrix(nrow=length(eps), ncol=length(gamma))
  out_tbv_cf <- matrix(nrow=length(eps), ncol=length(gamma))

  # return mpdi
  out_mpdi <- matrix(nrow=length(eps), ncol=length(gamma))
  # for each value of epsilon and each value of gamma
  # determine cost function minimizing the error in
  # both H and V brightness temperatures
  if (split) {

    out_tbs_h <- matrix(nrow = length(eps), ncol = length(gamma))
    out_tbc_h <- matrix(nrow = length(eps), ncol = length(gamma))

    out_tbs_v <- matrix(nrow = length(eps), ncol = length(gamma))
    out_tbc_v <- matrix(nrow = length(eps), ncol = length(gamma))

  }

  for (e in seq_along(fH)){
    for (g in seq_along(gamma)){
      # input each value of gamma and epsilon into temp retrieval
      eps_baseline = pred_H_V_tb(
                                 Tb_H = TbH, Tb_V = TbV,
                                 fH = fH[e], fV = fV[e],
                                 rhfac = rhfac,
                                 Tair = Tair,
                                 Tsoil = Tsoil,
                                 omega = omega,
                                 gamma = gamma[g],
                                 twotemps = twotemps,
                                 split = split, rtms=rtms)
      # export the cost function value for each gamma and epsilon
      # into matrix

      out_cf[e, g] <- eps_baseline[[1]]
      out_tbh[e, g] <- eps_baseline[[2]]
      out_tbv[e, g] <- eps_baseline[[3]]
      out_tbh_cf[e, g] <- eps_baseline[["tbh_cf"]]
      out_tbv_cf[e, g] <- eps_baseline[["tbv_cf"]]
      out_mpdi[e, g] <- eps_baseline[["mpdi"]]

      if(split){

        out_tbs_h[e,g] <- eps_baseline[[4]]
        out_tbc_h[e,g] <- eps_baseline[[5]]
        out_tbs_v[e,g] <- eps_baseline[[6]]
        out_tbc_v[e,g] <- eps_baseline[[7]]

      }

    }
  }

  # find the indices of the smallest smallest cost function value
  # returns a matrix dim of (1, 2)
  if(mpdi == F){
  min_index <- which(out_cf == min(out_cf, na.rm = TRUE), arr.ind = TRUE)
  # take gamma with the lowest cost function and calculate the VOD
  vod_out <- vod[min_index[2]]

  }else{
    min_index <- which(out_mpdi == min(out_mpdi, na.rm = TRUE), arr.ind = TRUE)
    vod_out <- vod[min_index[2]]
    }
  # # regularize the vod such that it is close to the mean value of the GNSS VOD
  # out_cf_regularized <- out_cf + lambda*(vod_out - 0.68)^2
  # min_index_regularized <- which(out_cf_regularized == min(out_cf_regularized, na.rm = TRUE),
  #                                arr.ind = TRUE)
  # vod_out_regularized <- vod[min_index_regularized[2]]

  if (split) {
    output <- list(
      min_index,
      costfx = out_cf[min_index[1], min_index[2]],
      epsilon = eps[min_index[1]],
      pred_tbh = out_tbh[min_index[1], min_index[2]][1],
      pred_tbv = out_tbv[min_index[1], min_index[2]][1],
      tb_soil_h = out_tbs_h[min_index[1], min_index[2]][1],
      tb_can_h = out_tbc_h[min_index[1], min_index[2]][1],
      tb_soil_v = out_tbs_v[min_index[1], min_index[2]][1],
      tb_can_v = out_tbc_v[min_index[1], min_index[2]][1],
      tbh_cf = out_tbh_cf[min_index[1], min_index[2]],
      tbv_cf = out_tbv_cf[min_index[1], min_index[2]],
      soil_moisture = smc[min_index[1]],
      vod = vod_out#,
      # vod_reg = vod_out_regularized,
      # cost_fx_reg = out_cf_regularized[min_index[1], min_index[2]]
      #cost_mat = out_cf
    )
  } else {

    output <- list(
      min_index,
      costfx = out_cf[min_index[1], min_index[2]],
      epsilon = eps[min_index[1]],
      pred_tbh = out_tbh[min_index[1], min_index[2]][1],
      pred_tbv = out_tbv[min_index[1], min_index[2]][1],
      soil_moisture = smc[min_index[1]],
      tbh_cf = out_tbh_cf[min_index[1], min_index[2]],
      tbv_cf = out_tbv_cf[min_index[1], min_index[2]],
      vod = vod_out,
      mpdi = out_mpdi[min_index[1], min_index[2]]#,
      # vod_reg = vod_out_regularized,
      # cost_fx_reg = out_cf_regularized[min_index[1], min_index[2]]
      #cost_mat = out_cf
    )
  }

  return(output)

}

# run the retrieval for various scenarios using the above functions.

ret_vod_dual <- function (TbH, TbV,# inputs of brightness temperatures
                          sm,
                          vod = seq(0, 3, by=0.01), # soil moisture values and VOD search values
                          Tair, # physical temperature
                          Tsoil = NULL, # soil temperature
                          omega= 0.05, # scattering albedo
                          roughness= 0.15, inc_angle= 40,
                          ret_fx = get_sm_vod, # function choice default is dual
                          twotemps=TRUE, mpdi= F,#lambda = 0,
                          split=F, rtms = F){ # are you using two temps (Ts and Tc)?

  if(twotemps==T & is.null(Tsoil)){stop("Tsoil must be provided when twotemps=T; ret_vod_dual")}
  if(rtms == T){warning("Note:rtms=T so watch out if that was not intended")}

  vod_xpol <- numeric(length(TbH))
  cost_xpol <- numeric(length(TbH))
  sm_xpol <- numeric(length(TbH))
  # vod_reg_xpol <- numeric(length(TbH))
  # cost_reg_xpol <- numeric(length(TbH))
  reflec_xpol <- numeric(length(TbH))
  tbpred_vpol <- numeric(length(TbH))
  tbpred_hpol <- numeric(length(TbH))
  tbv_cost <- numeric(length(TbH))
  tbh_cost <- numeric(length(TbH))
  mpdi_xpol <- numeric(length(TbH))

  if(length(omega)==1){o <- rep(omega, length(TbH))}else{o<-omega}

  for (i in seq_along(TbH)){

    est <- ret_fx(smc = sm[i], TbH=TbH[i], TbV=TbV[i],
                  vod = vod,
                  Tair = Tair[i], Tsoil=Tsoil[i],
                  omega = o[i], roughness=roughness, inc_angle = inc_angle,
                  twotemps=twotemps, rtms = rtms, mpdi = mpdi)

    vod_xpol[i]<- est$vod
    cost_xpol[i] <- est$costfx
    sm_xpol[i] <- est$soil_moisture
    tbpred_vpol[i] <- est$pred_tbv
    tbpred_hpol[i] <- est$pred_tbh
    tbh_cost[i] <- est$tbh_cf
    tbv_cost[i] <- est$tbv_cf
    mpdi_xpol[i] <- est$mpdi
    # vod_reg_xpol[i] <- est$vod_reg
    # cost_reg_xpol[i] <- est$cost_fx_reg


  }
  if(!identical(sm_xpol,sm)){
    warning("'Estimated' soil moisture values do not match input.")
  }

  return(list(vod = vod_xpol, sm = sm_xpol,
              pred_tbh1 = tbpred_hpol, pred_tbv1=tbpred_vpol,
              costfx = cost_xpol,

              # vod_reg = vod_reg_xpol,
              # costfx_reg=cost_reg_xpol,
              tbh_cost = tbh_cost, tbv_cost=tbv_cost))
}


## file management

load_data_retrieval <- \(){
  require(here)
  retrieval_output_file_names <- list.files(here('outputs'), pattern= "retrievals_output")

  retrieval_data_file_names <- list.files(here("data"), pattern = "data_for_retrieval")

  most_recent_retrieval_output <- retrieval_output_file_names[length(retrieval_output_file_names)]

  most_recent_retrieval_data <- retrieval_data_file_names[length(retrieval_data_file_names)]

  load(here("data", most_recent_retrieval_data),envir=globalenv())
  load(here("outputs", most_recent_retrieval_output),envir=globalenv())

  cat(paste("The data and retrieval files were loaded into the environment.\n\n most recent data file:",most_recent_retrieval_data,"\n most recent retrieval file:",most_recent_retrieval_output,"\n")
      )
}

# convert minutes to time
minute_to_time <- function(minute_of_day) {

  # Round the minute to the nearest 30
  rounded_minutes <- round(minute_of_day / 30) * 30

  hours <- floor(rounded_minutes / 60)
  minutes <- rounded_minutes %% 60

  # Format the output as HH:MM
  return(sprintf("%02d:%02d", hours, minutes))
}
