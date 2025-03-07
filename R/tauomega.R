#' Tau-omega model
#'
#' @param r_rough Rough surface reflectivity (i.e., smooth reflectivity * rhfac)
#' @param gamma Vegetation transmissivity transfer function. VOD divided by the cosine of incidence angle
#' @param Tair Air temperature in deg Kelvin
#' @param Tsoil Soil temperature in deg Kelvin
#' @param omega Scattering albedo
#' @param emiss Do you want to only calculate the emissitivity?
#'
#'
#' @details The reflectivities should be accounting to the roughness of the surface before input into this function. See `fresnelr.R`for more details
#'          for the way that the roughness factor was calculated and applied.
#' @returns Vector of brightness temperatures
#' @export
#'
tau_omega <- function(r_rough, gamma, Tair = NULL, Tsoil = NULL, omega, emiss = F) {
  if (emiss == F) {
    if (is.null(Tair) | is.null(Tsoil)) {
      stop("Tair or Tsoil not provided.")
    }

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
        r_rough = r_rough,
        gamma = gamma,
        Tsoil = Tsoil,
        Tair = Tair,
        omega = omega
      ),
      class = "tb"
    )

    return(to_out)

  } else {
    emissivity <- (1 - r_rough) * gamma + (1 - omega) * (1 - gamma) * (1 + (r_rough * gamma))
    # return the soil and canopy brightness temperatures and the total brightness temperature
    to_out <- structure(
      list(emissivity),
      .Names = "emissivity",
      inputs = list(
        r_rough = r_rough,
        gamma = gamma,
        omega = omega
      )
    )
    return(to_out)
  }
}

#' @export
print.tb <- function(x, ...) {
  cat("Soil Brightness Temperature:", x$tbs, "\n")
  cat("Canopy Brightness Temperature:", x$tbc, "\n")
  cat("Total:", x$totaltb, "\n")

  cat("\n To see the inputs, check attr(obj,'inputs')")
}
