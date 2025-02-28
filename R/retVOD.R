# run the retrieval for various scenarios using the above functions.

#' Title
#'
#' @param TbH
#' @param TbV
#' @param sm
#' @param vod
#' @param Tair
#' @param Tsoil
#' @param omega
#' @param roughness
#' @param inc_angle
#' @param ret_fx
#' @param twotemps
#' @param mpdi
#' @param split
#' @param rtms
#'
#' @returns
#' @export
#'
#' @examples
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
