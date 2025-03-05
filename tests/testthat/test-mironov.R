test_that("Check output", {
  # check that the output is of type complex
  # i.e., eps is of type complex
  test_miro <- mironov(f = 1.4e9, smv = 0.1, cf = 0.232)
  expect_type(test_miro$dielectric, "complex")
  # check that combined output and the inputs are the same
  expect_equal(test_miro$real, Re(test_miro$dielectric))
  expect_equal(test_miro$imaginary, -Im(test_miro$dielectric))
  # check that the output value is consistent
  epsilon <- 4.87186642-0.4439285i
  expect_equal(test_miro$dielectric,epsilon)
})

test_that("Check input", {
  f <- 1
  smv <- 2
  cf <- 3

  suppressWarnings({
    expect_error(mironov(f = f, smv = 0.1, cf = 0.1))
    #expect_error(mironov(f = 1.4e10, smv = smv, cf = 0.1))
    expect_error(mironov(f = 1.4e10, smv = 0.1, cf = cf))
  })
})
