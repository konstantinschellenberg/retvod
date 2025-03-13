#' Calculate residuals for predicted brightness temperatures
#'
#' @param tbH Obs.H polarization brightness temperature
#' @param tbV Obs.V polarization brightness temperature
#' @param tbHpred Pred. H polarization brightness temperature
#' @param tbVpred Pred. V polarization brightness temperature
#'
#' @return Residuals for each predicted brightness temperature
#'
#' @export

tbResiduals <- function(tbH, tbV, tbHpred, tbVpred) {
  if(!length(tbH)==length(tbHpred)){
    stop("tbH and tbHpred have different lengths.")
  }
  if(!length(tbV)==length(tbVpred)){
    stop("tbV and tbVpred have different lengths.")
  }
  if(identical(tbH,tbHpred)|identical(tbV,tbVpred)){
    warning("Input vectors are identical.")
  }
  # compute the squared differences between predicted and observed Tb
  tbH_res <- (tbHpred - tbH)^2
  tbV_res <- (tbVpred - tbV)^2
  # cost function is the sum of the two square differences
  totaltb_res <- tbH_res + tbV_res
  rse <- sqrt(totaltb_res)
  # return the residuals for each polarization, the total, and sqrt(total)
  res <- structure(
    list(
      tbH_res,
      tbV_res,
      totaltb_res,
      rse
    ),
    .Names = c("tbH", "tbV", "totaltb", "rootSE"),
    class = c("list", "res")
  )
  return(res)
}
