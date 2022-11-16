test_that("median value is correct",{
  expect_equal(mediantest(c(1,3),c(5,-1),paired=TRUE)$estimate,0)
  expect_equal(mediantest(c(1,3))$estimate,2)
})




