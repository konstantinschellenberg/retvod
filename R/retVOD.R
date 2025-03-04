#' Retrieve VOD for all values
#'
#' @param tbH Horizontal brightness temperatures
#' @param tbV Vertical brightness temperatures
#' @param sm Soil moisture
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

  lengths_equal <- all(sapply(list(tbH, tbV, smc), length) == length(tbH))

  if (!lengths_equal) {
    stop("tbH, tbV, and smc lengths differ.")
  }

  # Prepare omega vector (if only one value, repeat for each TbH)
  o <- if (length(omega) == 1) rep(omega, length(tbH)) else omega

  # initialize lists for the results
  results <- vector("list", length(tbH))
  names(results) <- c(
    "vod_xpol", "cost_xpol", "sm_xpol",
    "tbpred_vpol", "tbpred_hpol",
    "tbh_cost", "tbv_cost"
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
    results$vod_xpol[i] <- est$vod_est
    results$cost_xpol[i] <- est$cf_tb
    results$sm_xpol[i] <- est$sm_est
    results$tbpred_vpol[i] <- est$pred_tbV
    results$tbpred_hpol[i] <- est$pred_tbH
    results$tbh_cost[i] <- est$cf_tbH
    results$tbv_cost[i] <- est$cf_tbV
  }

  if (!identical(results$sm_xpol, smc)) {
    warning("'Estimated' soil moisture values do not match input.")
  }

  return(structure(results,
    creation_time = Sys.time()
  ))
}
