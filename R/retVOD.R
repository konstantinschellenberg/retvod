#' Retrieve VOD for all values
#'
#' @param tbH Horizontal brightness temperatures
#' @param tbV Vertical brightness temperatures
#' @param smc Soil moisture
#' @param vod VOD estimate
#' @param Tair Air temperature
#' @param Tsoil Soil temperature
#' @param omega Canopy opacity, scattering albedo
#' @param roughness Soil roughness estimate
#' @param inc_angle Incidence angle
#'
#' @return Retrieved VOD and auxillary information
#' @export
#'
#' @import cli

retVOD <- function(tbH, tbV,
                   smc,
                   vod,
                   Tair,
                   Tsoil,
                   omega,
                   roughness,
                   inc_angle) {

  if (length(tbH)!= length(tbV)|| length(tbH)!= length(smc)) {
    stop("tbH, tbV, and smc lengths differ.")
  }

  # Prepare omega vector (if only one value, repeat for each TbH)
  o <- if (length(omega) == 1) rep_len(omega, length(tbH)) else omega

  # initialize lists for the results
  results <- vector("list", length(tbH))
  names(results) <- c(
    "vodEst", "cfEst", "smEst",
     "tbHpred","tbVpred",
    "tbHcost", "tbVcost"
  )

  # for each brightness temperature retrieve VOD and soil moisture
  cli::cli_progress_bar("Retrieving VOD values...", total = length(tbH))
  for (i in seq_along(tbH)) {
    est <- solveSmVod(
      smc = smc[i], tbH = tbH[i], tbV = tbV[i],
      vod = vod,
      Tair = Tair[i], Tsoil = Tsoil[i],
      omega = o[i], roughness = roughness, inc_angle = inc_angle
    )

    # return elements from retrieval

    results$vodEst[i] <- est$vod_est
    results$cfEst[i] <- est$cf_tb
    results$smEst[i] <- est$sm_est
    results$tbVpred[i] <- est$pred_tbV
    results$tbHpred[i] <- est$pred_tbH
    results$tbHcost[i] <- est$cf_tbH
    results$tbVcost[i] <- est$cf_tbV

    cli_progress_update()  # Update progress bar
  }

  results$rmse <- mean(sqrt(results$cfEst), na.rm=T)

  cli_progress_done()  # Complete progress bar
  cli_alert_success("VOD retrieval complete!")

  if (!identical(results$smEst, smc)) {
    warning("'Estimated' soil moisture values do not match input.")
  }

  return(structure(results,
    creation_time = Sys.time(),
    inputs = list(h=tbH, v=tbV, sm = smc, omega = omega, rough = roughness,angle =inc_angle),
    class = "retVOD"
  ))
}


#' @export
#' @importFrom graphics hist par abline
#' @importFrom grDevices dev.new
plot.retVOD <- function(x, ...){

  dev.new(width = 10, height = 6)

  inputs <- attr(x, "inputs")
  par(mfrow = c(3,3), mar = c(5, 5, 1, 1))

  plot(x$vodEst, main = "Retrieved Vegetation Optical Depth",
       ylab = "VOD",
       xlab = "Index")

  plot(x$cfEst, main = "Total Tb Residuals",
       ylab = "Brightness Temps (K^2)",
       xlab = "Index")

  hist(x$cfEst, main = "Histogram Tb Residuals", xlab= "Tb Residuals")

  plot(x$tbHcost, main = "TbH Residuals (K^2)",
       ylab = "Brightness Temps (K^2)",
       xlab = "Index")

  plot(x$tbVcost, main = "TbV Residuals (K^2)",
       ylab = "Brightness Temps (K^2)",
       xlab = "Index")

  plot(x$tbHpred ~ x$tbHcost, main = "Predicted TbH ~ TbH Residuals",
       ylab = "Brightness Temps (K)",
       xlab = "Tb H Residuals (K^2)")
  abline(a = 0, b = -1, col = "red")
  plot(x$tbVpred ~ x$tbVcost, main = "Predicted TbV ~ TbV Residuals",
       ylab = "Brightness Temps (K)",
       xlab = "Tb V Residuals (K^2)")
  abline(a = 0, b = -1, col = "red")

  plot(inputs$h ~ x$tbHcost, main = "Observed TbH ~ TbH Residuals",
       ylab = "Brightness Temps (K)",
       xlab = "Tb H Residuals (K^2)")
  abline(a = 0, b = 1, col = "red")
  plot(inputs$v ~ x$tbVcost, main = "Observed TbV ~ TbV Residuals",
       ylab = "Brightness Temps (K)",
       xlab = "Tb V Residuals (K^2)")
  abline(a = 0, b = 1, col = "red")

}
