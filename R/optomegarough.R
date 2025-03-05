#' Optimize omega and roughness
#'
#' @param df A data.frame
#' @param range_omega Range of test omega values
#' @param range_roughness Range of test roughness values
#' @param precision Numeric indicating increment step of the grid search
#' @return List of all residuals and the best residuals
#' @export
#'
#' @importFrom tidyr pivot_longer
#' @importFrom tidyselect contains

optOR <- function(df, range_omega, range_roughness, precision) {
  optOut <- vector("list", length(range_omega) * length(range_roughness))
  names(optOut) <- character(length(optOut))

  index <- 1
  for (o in seq_along(range_omega)) {
    for (r in seq_along(range_roughness)) {
      retrieval <- retVOD(
        tbH = df[["tbh"]], tbV = df[["tbv"]],
        smc = df[["mean_smc"]], vod = seq(0,3, by = precision),
        Tair = df[["phys_temp"]],
        Tsoil = df[["mean_soil_temp"]],
        omega = range_omega[o],
        roughness = range_roughness[r],
        inc_angle = 40
      )

      name <- paste("omega", round(o, 5), "rough", r, sep = "_")
      print(name)
      optOut[[index]] <- sum(retrieval$costfx, na.rm = T)
      names(optOut)[index] <- name
      index <- index + 1
    }
  }

  # Convert results list to dataframe
  optimization_df <- do.call(cbind, optOut) |>
    as.data.frame() |>
    tidyr::pivot_longer(cols = tidyselect::contains("omega"), names_to = "omega_roughness", values_to = "cf")

  # Find best parameters
  best_params <- optimization_df[which.min(optimization_df$cf), ]

  return(structure(list(results = optimization_df, best = best_params, pres = precision),
         creation_time = Sys.time()))
}

