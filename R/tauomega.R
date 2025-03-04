#' Tau-omega model
#'
#' @param r Fresnel reflectivity
#' @param rhfac Roughness factor to apply to `r`
#' @param gamma Vegetation transmissivity transfer function. VOD divided by the cosine of incidence angle
#' @param Tair Air temperature in deg Kelvin
#' @param Tsoil Soil temperature in deg Kelvin
#' @param omega Impediance of canopy?
#' @param emiss Do you want to only calculate the emissitivity?
#'
#' @returns Vector of brightness temperatures
#' @export
#'
tau_omega <- function(r, rhfac, gamma, Tair=NULL, Tsoil=NULL, omega, emiss = F) {

if(is.null(Tair)|is.null(Tsoil)){
  stop("Tair or Tsoil not provided.")
}

if(emiss == F){
  # calculate rough fresnel reflectivity
  r_rough <- rhfac * r
  # calculate the soil-associated brightness temperature
  tbs <- (1 - r_rough) * gamma * Tsoil
  # calculate the canopy-associated brightness temperature
  tbc <- (1 - omega) * (1 - gamma) * (1 + (r_rough * gamma)) * Tair
  # calculate the total brightness temperature
  totaltb <- tbs + tbc

  # return the soil and canopy brightness temperatures and the total brightness temperature
  to_out <- structure(
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
  )

  return(to_out)

}else{

    r_rough <- rhfac * r
    emissivity <- (1 - r_rough) * gamma + (1 - omega) * (1 - gamma) * (1 + (r_rough * gamma))
    # return the soil and canopy brightness temperatures and the total brightness temperature
    to_out<-structure(
      list(emissivity),
      .Names = "emissivity",
      inputs = list(
        r = r,
        rhfac = rhfac,
        gamma = gamma,
        Tsoil = Tsoil,
        Tair = Tair,
        omega = omega
      )
    )
    return(to_out)
  }

}

#' @export
print.tb <- function(x, ...){

  cat("Soil Brightness Temperature:", x$tbs, "\n")
  cat("Canopy Brightness Temperature:", x$tbc, "\n")
  cat("Total:", x$totaltb, "\n")

  cat("\n To see the inputs, check attr(obj,'inputs')")

}
