#' Tau-omega model
#'
#' @param r Fresnel reflectivity
#' @param rhfac Roughness factor to apply to `r`
#' @param gamma Vegetation transmissivity transfer function. VOD divided by the cosine of incidence angle
#' @param Tair Air temperature in deg Kelvin
#' @param Tsoil Soil temperature in deg Kelvin
#' @param omega Impediance of canopy?
#' @param noT Do you want to only calculate the emissitivity?
# @param rtms
#'
#' @returns Vector of brightness temperatures
#' @export
#'
tau_omega <- function(r, rhfac, gamma, Tair, Tsoil, omega, noT = F) {

  # calculate rough fresnel reflectivity
  r_rough <- rhfac * r
  # calculate the soil-associated brightness temperature
  tbs <- (1 - rhfac * r) * gamma * Tsoil
  # calculate the canopy-associated brightness temperature
  tbc <- (1 - omega) * (1 - gamma) * (1 + ((rhfac * r) * gamma)) * Tair
  # calculate the total brightness temperature
  totaltb <- tbs + tbc
  # return the soil and canopy brightness temperatures and the total brightness temperature
  invisible(structure(
    list(tbs, tbc, totaltb),
    names = c("tbs", "tbc", "totaltb"),
    inputs = list(
      r = r,
      rhfac = rhfac,
      gamma = gamma,
      Tsoil = Tsoil,
      Tair = Tair,
      omega = omega
    ),
    class = "tb"
  ))
  #
#   if (rtms) {
#
#     tb <- RTMS_tb(r = r, rhfac = rhfac, gamma = gamma, Tsoil = Tsoil, Tveg = Tair, omega = omega)[[4]]
#     return(tb)
#
#   } else if (noT == TRUE) {
#     r_rough <- rhfac * r
#     return((1 - r_rough) * gamma + (1 - omega) * (1 - gamma) * (1 + (r_rough * gamma)))
#   }
}


#' @export
print.tb <- function(x, ...){

  cat("Soil Brightness Temperature:", x$tbs, "\n")
  cat("Canopy Brightness Temperature:", x$tbc, "\n")
  cat("Total:", x$totaltb, "\n")

  cat("\n To see the inputs, check attr(obj,'inputs')")

}
