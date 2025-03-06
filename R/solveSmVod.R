#' Solve for soil moisture and VOD
#'
#' @param reflecs H and V pol reflectivities
#' @param gamma Gamma estimates for each test VOD value
#' @param tbH Horizontal brightness temperature
#' @param tbV Vertical brightness temperature
#' @param Tair Air temperature
#' @param Tsoil Soil temperature
#' @param omega Scattering albedo
#' @param mat TRUE returns the cost function matrix
#'
#' @description This function uses the a range of input soil moisture and vod values to solve for the best VOD value given the air, soil, and observed brightness temperatures.
#' @details The function returns the predicted brightness temperatures (i.e., using tau-omega), residuals for each polarization, and the predicted soil moisture and VOD. Mironov is used to determine the dielectric constant for a given soil moisture and has a set frequency of 1.4e9Hz (for L-band). Also the clay fraction at MOFLUX is 23.2%.
#' @return List of predicted brightness temperatures, soil moisture and VOD:
#'
#' @export
#'
solveSmVod <- function(reflecs,
                       gamma,
                       tbH, tbV,
                       Tair, Tsoil,
                       omega, mat=F) {
  ## initialize output matrices
  num_r <- length(reflecs) # number of dielectric values
  num_gamma <- length(gamma) # number of test VOD values

  # array to hold residuals and predictions for each epsilon and gamma
  results <- array(NA,
    dim = c(num_eps, num_vod, 5), # rows, columns, calc values
    dimnames = list(smc, vod, c("pred_tbH", "pred_tbV", "cf_total", "cf_tbH", "cf_tbV"))
  )

  ## Compute cost function for all combinations of epsilon and VOD
  for (e in seq_along(smc)) {
    for (g in seq_along(vod)) {
      result <- estTb(
        tbH = tbH, tbV = tbV,
        fH = reflecs[[e]][["fH"]], fV = reflecs[[e]][["fH"]],
        gamma = gamma[g],
        Tair = Tair,
        Tsoil = Tsoil,
        omega = omega
      )
      # store results
      results[e,g, ] <- c(result$pred_tbH, result$pred_tbV, result$residuals$totaltb,
                          result$residuals$tbH, result$residuals$tbV)
    }
  }

  min_index <- which(results[, , "cf_total"] == min(results[, , "cf_total"],
                                                    na.rm = TRUE),

    arr.ind = TRUE
  )

  if(length(smc)==1){
    best_row <- 1
    best_col <- min_index
  }else if (length(vod)==1){
    best_row <- min_index
    best_col <- 1
  }else{
    best_row <- min_index[1]
    best_col <- min_index[2]
  }

  output <- list(
    min_cf_index = min_index,
    cf_tb = results[best_row, best_col, "cf_total"]|>unname(),
    epsilon = eps_list[best_row],
    pred_tbH = results[best_row, best_col, "pred_tbH"]|>unname(),
    pred_tbV = results[best_row, best_col, "pred_tbV"]|>unname(),
    cf_tbH = results[best_row, best_col, "cf_tbH"]|>unname(),
    cf_tbV = results[best_row, best_col, "cf_tbV"]|>unname(),
    sm_est = smc[best_row],
    vod_est = vod[best_col],
    cf_mat = if (mat == T) results[, , "cf_total"] else NA
  )
  return(structure(output,
                   flag = if (length(best_row)>1) "Tie for lowest residuals found")
  )
}
