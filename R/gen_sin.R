#' Generate sinusoid
#'
#' @param len length of the output vector
#' @param rangeL lower bound to values
#' @param rangeH upper bound to values

gen_sin <- function(len, rangeL, rangeH){

  amplitude = (rangeH-rangeL)/2
  offset = (rangeH + rangeL)/2

  # Generate sinusoidal values
  sinusoid <- offset + amplitude * sin(seq(0, 2 * pi, length.out = len))

  return(sinusoid)
}
