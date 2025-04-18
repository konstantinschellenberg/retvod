#' Fresnel's reflectivity calculation
#'
#' @param eps Dielectric
#' @param theta Incidence angle in degrees
#' @param h Soil surface roughness
#'
#' @return List of horizontal and vertical polarization reflectivity
#' @export
#'
#'
fresnelr <- function(eps, theta, h){

  if (Re(eps)<0|Im(eps)==0){
    stop("fresnelr: Invalid epsilon value.")
  }

  stopifnot(theta>0)# no negative angles

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
