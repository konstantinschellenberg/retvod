test_that("Check output", {

  #Check that you get the expected value
  test_eps <- mironov(1.4e9, 0.4, 0.232)

  reflecs <- fresnelr(test_eps$real, test_eps$imaginary, 40, h=0.1)

  # Expected values for known inputs (replace with actual known correct values)
  thetar = 40*pi/180
  h=0.1
  rhfac = exp(-h * cos(thetar)) # could be squared here
  hnum = cos(thetar) - sqrt(test_eps$dielectric - sin(thetar)^2)
  hdem = cos(thetar) + sqrt(test_eps$dielectric - sin(thetar)^2)

  vnum = test_eps$dielectric*cos(thetar) - sqrt(test_eps$dielectric - sin(thetar)^2)
  vdem = test_eps$dielectric*cos(thetar) + sqrt(test_eps$dielectric - sin(thetar)^2)

  expected_fH <- abs(hnum/hdem)^2 * rhfac
  expected_fV <- abs(vnum/vdem)^2 * rhfac

  expect_true(reflecs$fH > reflecs$fV)
  expect_true(reflecs$fH==expected_fH)
  expect_true(reflecs$fV==expected_fV)
})

test_that("Function handles invalid inputs", {

  eps <- complex(real=10, imaginary=6)
  # Negative epsilon values
  expect_error(fresnelr(-1, 30, h=0.1))

  # Zero epsilon values (invalid as refractive index cannot be zero)
  expect_error(fresnelr(0, 30, h = 0.1))

  # Negative angle no bueno
  expect_error(fresnelr(eps, -40, h=0.1))
})

test_that("Total reflectivity and transmissivity should sum to 1", {
  eps <- mironov(1.4e9, 0.4, 0.232)$dielectric  # Example refractive index
  theta <- 40   # Example angle (degrees)
  h = 0.1
  # Call the fresnelr function
  result <- fresnelr(eps, theta, h)

  fH <- result$fH
  fV <- result$fV

  # Total reflectivity (can also include polarization weightings)
  total_reflectivity <- fH + fV

  # Total transmissivity is simply 1 - total reflectivity
  total_transmissivity <- 1 - total_reflectivity

   # Check if total reflectivity + transmissivity is approximately 1
  expect_equal(total_reflectivity + total_transmissivity, 1, info = paste("Total reflectivity and transmissivity do not sum to 1 at angle", theta))
  # tolerance = 1e-3,
})
