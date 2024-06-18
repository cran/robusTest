set.seed(1)
n<-40
x<-rnorm(n)
y<-x^2+0.3*rnorm(n)
robust_res=indeptest(x,y)

test_that("indeptest on Simu1",{
  expect_equal(robust_res$statistic,stat_indeptest(x,y))
  expect_equal(round(robust_res$statistic,4),0.8696)
  expect_equal(round(robust_res$p.value,2),0.01)
  #expect_equal(round(robust_res$p.value,2),round(indeptest(x,y,N=50000)$p.value,2))
})

x1<-c(0.2, 0.3, 0.1, 0.4)
y1<-c(0.5, 0.4, 0.05, 0.2)

x2<-c(0.1, 0.3, 0.5, 0.3)
y2<-c(0.5, 0.45, 0.5, 0.25)

test_that("warning for indeptest when there are ties",{
  expect_warning(indeptest(x1,y1))
  expect_warning(indeptest(x2,y1))
  expect_warning(indeptest(x1,y2))
  expect_warning(indeptest(x2,y2))
  expect_silent(indeptest(x2,y2,ties.break = "random"))
})

x <- c(1:150)*2
y <- rep(c(-0.95,  0.73,  1.79, -0.64, -0.50),30)
y <-round(y+seq(-0.1,0.1,length.out=150),4)

result=indeptest(x,y)

test_that("Independent test on deterministic (n=150) sample",{
  expect_equal(round(result$statistic,4),0.6695)
  #expect_equal(round(result$p.value,3),0.214)
  expect_true(abs(0.2144-result$p.value)<=0.02)
})

# n=40
# power_res=rep(NA,1000)
# for (i in 1:1000)
# {
#   x<-rnorm(n)
#   y<-x^2+0.3*rnorm(n)
#   power_res[i]=indeptest(x,y)$p.value
# }
#
# reject_res=rep(NA,1000)
# for (i in 1:1000)
# {
#   x<-rnorm(n)
#   y<-rnorm(n)
#   reject_res[i]=indeptest(x,y)$p.value
# }
# skip_on_cran()
# test_that("indeptest has the correct rejection rate in a simulation setting",{
#   expect_gt(mean(reject_res<=0.05),0.03)
#   expect_lt(mean(reject_res<=0.05),0.07)
# })

# # set.seed(1)
# # n<-4000
# # x<-rnorm(n)
# # y<-0.01*x^2+rnorm(n)
# # result=indeptest(x,y)
# # result$statistic #0.5973543
# # result$p.value#0.476434
# #
# # set.seed(1)
# # n<-40
# # x<-rnorm(n)
# # y<-0.1*x^2+rnorm(n)
# # result=indeptest(x,y)
# # result$statistic #0.54154
# # result$p.value#0.472538
# #
# # set.seed(1)
# # n<-80
# # x<-rnorm(n)
# # y<-0.1*x^2+rnorm(n)
# # result=indeptest(x,y)
# # result$statistic #0.5757875
# # result$p.value# 0.425226
# #
# #
# # set.seed(1)
# # n<-130
# # x<-rnorm(n)
# # y<-0.1*x^2+rnorm(n)
# # result=indeptest(x,y)
# # result$statistic #0.6672388
# # result$p.value# 0.21383
# #
# # set.seed(1)
# # n<-140
# # x<-rnorm(n)
# # y<-0.5*x^2+rnorm(n)
# # result=indeptest(x,y)
# # result$statistic #0.5433134
# # result$p.value# 0.582986
#
# #
# #
# #
# #
