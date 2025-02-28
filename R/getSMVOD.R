#' Title
#'
#' @param smc
#' @param TbH
#' @param TbV
#' @param Tair
#' @param Tsoil
#' @param vod
#' @param omega
#' @param inc_angle
#' @param cf
#' @param roughness
#' @param twotemps
#' @param split
#' @param rtms
#' @param mpdi
#'
#' @returns
#' @export
#'
#' @examples
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
