test_that("All residuals are returned", {

H <- seq(200,250, length.out = 10)
V <- seq(260,280, length.out = 10)

pH <- seq(200,250, length.out = 10) + 1
pV <- seq(260,280, length.out = 10) + 2

res_test <- tbResiduals(tbH = H, tbV = V, tbHpred = pH, tbVpred = pV)

expect_true(all(res_test$tbH==1))
expect_true(all(res_test$tbV==4))
expect_true(all(res_test$totaltb==5))
expect_true(all(res_test$rootSE==sqrt(res_test$totaltb)))

})


test_that("Error handling works", {
  H <- seq(200,250, length.out = 10)
  V <- seq(260,280, length.out = 10)

  pH <- seq(200,250, length.out = 10) + 1
  pV <- seq(260,280, length.out = 10) + 2

  expect_error(tbResiduals(tbH = H, tbV = V, tbHpred = c(pH,1), tbVpred = pV),"tbH and tbHpred have different lengths.")
  expect_error(tbResiduals(tbH = H, tbV = V, tbHpred = pH, tbVpred = c(pV,1)),"tbV and tbVpred have different lengths.")
  expect_warning(tbResiduals(tbH = H, tbV = V, tbHpred = H, tbVpred = pV),"Input vectors are identical.")

})
