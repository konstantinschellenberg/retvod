#' Solve for soil moisture and VOD
#'
#' @param smc Soil moisture
#' @param vod Test VOD values
#' @param tbH Horizontal brightness temperature
#' @param tbV Vertical brightness temperature
#' @param Tair Air temperature
#' @param Tsoil Soil temperature
#' @param omega Scattering albedo
#' @param inc_angle Incidence angle
#' @param clay_frac Clay fraction
#' @param roughness Roughness estimate
#'
#' @description This function uses the a range of input soil moisture and vod values to solve for the best VOD value given the air, soil, and observed brightness temperatures.
#' @details The function returns the predicted brightness temperatures (i.e., using tau-omega), residuals for each polarization, and the predicted soil moisture and VOD. Mironov is used to determine the dielectric constant for a given soil moisture and has a set frequency of 1.4e9Hz (for L-band). Also the clay fraction at MOFLUX is 23.2%.
#' @return List of predicted brightness temperatures, soil moisture and VOD:
#'
#' @export
#'
solveSmVod <- function(smc,
                       vod,
                       tbH, tbV,
                       Tair, Tsoil,
                       omega, inc_angle,
                       clay_frac = 0.232, roughness, mat=F) {

  ## calculate epsilon (dielectric) for each value of soil moisture
  eps_list <- sapply(smc, \(s) mironov(1.4e9, s, clay_frac)$dielectric)

  ## calculate fresnel h and v polarization reflectivities
  reflecs <- sapply(eps_list, \(e) fresnelr(eps = e, theta = inc_angle), simplify = TRUE)
  fH <- reflecs["fH", ] |> unlist() # Extract Horizontal Reflectivity
  fV <- reflecs["fV", ] |> unlist() # Extract Vertical Reflectivity

  ## calculate gamma
  cosTheta <- cos(inc_angle * (pi / 180))
  gamma <- exp(-1 * (vod / cosTheta))

  ## apply roughness correction
  rhfac <- exp(-roughness * cosTheta) # could be squared here

  ## initialize output matrices
  num_eps <- length(eps_list) # number of dielectric values
  num_vod <- length(vod) # number of test VOD values

  # array to hold residuals and predictions for each epsilon and gamma
  results <- array(NA,
    dim = c(num_eps, num_vod, 5), # rows, columns, calc values
    dimnames = list(NULL, NULL, c("pred_tbH", "pred_tbV", "cf_total", "cf_tbH", "cf_tbV"))
  )

  ## Compute cost function for all combinations of epsilon and VOD
  for (e in seq_along(fH)) {
    for (g in seq_along(gamma)) {
      result <- estTb(
        tbH = tbH, tbV = tbV,
        fH = fH[e], fV = fV[e], gamma = gamma[g],
        rhfac = rhfac,
        Tair = Tair,
        Tsoil = Tsoil,
        omega = omega
      )
      # store results
      results[e, g, "pred_tbH"] <- result$pred_tbH
      results[e, g, "pred_tbV"] <- result$pred_tbV

      results[e, g, "cf_total"] <- result$residuals$totaltb
      results[e, g, "cf_tbH"] <- result$residuals$tbH
      results[e, g, "cf_tbV"] <- result$residuals$tbV
    }
  }

  min_index <- which(results[, , "cf_total"] == min(results[, , "cf_total"],
                                                    na.rm = TRUE),

    arr.ind = TRUE
  )

  vod_solution <- vod[min_index[2]]

  output <- list(
    min_cf_index = min_index,
    cf_tb = results[min_index[1], min_index[2], "cf_total"]|>unname(),
    epsilon = eps_list[min_index[1]],
    pred_tbH = results[min_index[1], min_index[2], "pred_tbH"]|>unname(),
    pred_tbV = results[min_index[1], min_index[2], "pred_tbV"]|>unname(),
    cf_tbH = results[min_index[1], min_index[2], "cf_tbH"]|>unname(),
    cf_tbV = results[min_index[1], min_index[2], "cf_tbV"]|>unname(),
    sm_est = smc[min_index[1]],
    vod_est = vod_solution,
    cf_mat = if (mat == T) results[, , "cf_total"]
  )
  return(structure(output,
                   flag = if (nrow(min_index>1)) "Tie for lowest residuals found")
  )
}
