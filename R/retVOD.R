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
  }

  if (!identical(results$smEst, smc)) {
    warning("'Estimated' soil moisture values do not match input.")
  }

  return(structure(results,
    creation_time = Sys.time()
  ))
}
