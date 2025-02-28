#' Predict brightness temperatures
#'
#' @param TbH Observed horizontally polarized brightness temperature
#' @param TbV Observed vertically polarized brightness temperature
#' @param fH Estimated horizontal fresnel reflectivity
#' @param fV Estimated vertical fresnel reflectivity
#' @param gamma
#' @param Tair
#' @param Tsoil
#' @param omega
#' @param rhfac
#'
#' @return Vertical and horizontally polarized brightness temperatures
#' @export

estTb <- function(
    tbH, tbV,
    fH, fV,
    gamma,
    Tair, Tsoil,
    omega,
    rhfac) {
  # take a vector or single value of possible reflectivity values over a realistic range of eps

  # estimate h pol brightness temperature
  pred_T_H <- tau_omega(
    r = fH,
    rhfac = rhfac,
    gamma = gamma,
    Tair = Tair, Tsoil = Tsoil,
    omega = omega, twotemps = twotemps, split = split, rtms = rtms
  )
  # predicted v pol brightness temperature
  pred_T_V <- tau_omega(
    r = fV,
    rhfac = rhfac,
    gamma = gamma,
    Tair = Tair, Tsoil = Tsoil,
    omega = omega, twotemps = twotemps, split = split, rtms = rtms
  )

    pred_T_H <- pred_T_H$tb
    pred_T_V <- pred_T_V$tb

    tbs_h <- pred_T_H$tbs
    tbc_h <- pred_T_H$tbc
    tbs_v <- pred_T_V$tbs
    tbc_v <- pred_T_V$tbc

    # compute the squared differences between predicted and observed Tb
    costfun <- (pred_T_H - Tb_H)^2 + (pred_T_V - Tb_V)^2 # (pred_T_H - Tb_H)^2 +
    tb_output <- list(
      cf = costfun, # return the cost function value
      predTBH = pred_T_H, # return brightness temperatures
      predTBV = pred_T_V,
      tb_soil_h = tbs_h, tb_can_h = tbc_h,
      tb_soil_v = tbs_v, tb_can_v = tbc_v
    )

  return(tb_output)
}
