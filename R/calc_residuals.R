tbResiduals <- function(TbH, TbV, TbHpred, TbVpred) {
  # compute the squared differences between predicted and observed Tb
  tbh_res <- (TbHpred - TbH)^2
  tbv_res <- (TbVpred - TbV)^2
  # cost function is the sum of the two square differences
  cf <- tbh_res + tbv_res
  rse <- sqrt(costfun)
  # return the residuals for each polarization, the total, and sqrt(total)
  res <- structure(
    list(
      tbh_res,
      tbv_res,
      cf,
      rse
    ),
    .Names = c("TbH", "TbV", "costfun", "rootSE")
  )
  return(res)
}
