#' Predict brightness temperatures
#'
#' @param tbH Observed horizontally polarized brightness temperature
#' @param tbV Observed vertically polarized brightness temperature
#' @param fH Estimated horizontal fresnel reflectivity
#' @param fV Estimated vertical fresnel reflectivity
#' @param gamma Transfer function
#' @param Tair Air temperature
#' @param Tsoil Soil temperature
#' @param omega Canopy opacity
#' @param rhfac Roughness factor
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

  # Estimate brightness temperature for H and V polarizations
  pred_tbH <- tau_omega(r = fH, rhfac = rhfac, gamma = gamma, Tair = Tair, Tsoil = Tsoil, omega = omega)$totaltb
  pred_tbV <- tau_omega(r = fV, rhfac = rhfac, gamma = gamma, Tair = Tair, Tsoil = Tsoil, omega = omega)$totaltb

  # Compute the squared difference between predicted and observed brightness temperatures
  res <- tbResiduals(tbH = tbH, tbV = tbV, tbHpred = pred_tbH, tbVpred = pred_tbV)

  # Return the cost function value and predicted brightness temperatures
  return(structure(list(res, pred_tbH, pred_tbV),
                   .Names = c("residuals", "pred_tbH", "pred_tbV"))
  )
}
