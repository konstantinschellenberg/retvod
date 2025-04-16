#' Fresnel's reflectivity calculation
#'
#' @param eps_real Real part of the complex dielectric constant
#' @param eps_imag Imaginary part of the complex dielectric constant
#' @param theta Incidence angle in degrees
#' @param h Soil surface roughness
#'
#' @return List of horizontal and vertical polarization reflectivity
#' @export
#'
fresnelr <- function(eps_real, eps_imag, theta, h){

  # if (eps_real<0|eps_imag==0){
  #   stop("fresnelr: Invalid epsilon value.")
  # }
  #
  # stopifnot(theta>0)# no negative angles

  # concatenate real and imaginary parts of epsilon
  eps <- complex(real=eps_real, imaginary=eps_imag)

  # theta is converted to radians
  thetar = theta*pi/180

  # calculate cos(theta) and sin2(theta)
  costheta = cos(thetar)
  sin2theta = sin(thetar)^2
  # calculate fresnel h and v polarization reflectivity
  fH = abs((costheta - sqrt(eps - sin2theta)) / (costheta + sqrt(eps - sin2theta)))^2
  fV = abs((eps*costheta - sqrt(eps - sin2theta)) / (eps*costheta + sqrt(eps - sin2theta)))^2
  # calculate the roughness factor to apply to smooth reflectivities
  rhfac <- exp(-h * costheta) # could be squared here
  fH_rough <- fH * rhfac
  fV_rough <- fV * rhfac

  #return the rough fresnel reflectivities
  return(list(fH=fH_rough,
              fV=fV_rough))
}
