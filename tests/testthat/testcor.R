test_that("Pearson correlation value, statistic and pvalue are correct on Evans dataset",{
  expect_equal(with(Evans,cortest(CHL[CDH==1],DBP[CDH==1]))$estimate,as.numeric(with(Evans,cor.test(CHL[CDH==1],DBP[CDH==1]))$estimate))
  expect_equal(round(with(Evans,cortest(CHL[CDH==1],DBP[CDH==1]))$p.value,3),0.017)
  expect_equal(round(with(Evans,cortest(CHL[CDH==1],DBP[CDH==1]))$statistic,4),2.4126)
})

test_that("Kendall correlation value, statistic and pvalue are correct",{
  expect_warning(with(Evans,cortest(CHL[CDH==1],DBP[CDH==1],method="kendall")))
})

x <- c(-1.80592115,-0.04479747,-1.23845984, 0.06650748,
       0.69532038, -0.85064100, -0.52527723, -0.57914233,
       0.19006873, 0.79010406)

epsi <- c(1.8250331, -0.4214741, -0.1426509,  1.2405776,
          0.4579464, 1.7526687, 0.2209095, -1.9728033,
          -1.3070188, 2.5863419)
y <- x^2+0.3*epsi

result_Pearson <- cortest(x,y)

test_that("Pearson correlation value, statistic and pvalue are correct on deterministic (small) sample",{
  expect_equal(round(result_Pearson$estimate,4),-0.5946)
  expect_equal(round(result_Pearson$statistic,4),-1.3177)
  expect_equal(round(result_Pearson$p.value,3),0.285)
})

x <- c(1:150)*2
y <- rep(c(-0.95,  0.73,  1.79, -0.64, -0.50),30)
result_Pearson2 <- cortest(x,y)

test_that("Pearson correlation value, statistic and pvalue are correct on deterministic (large) sample",{
  expect_equal(round(result_Pearson2$estimate,4),-0.0021)
  expect_equal(round(result_Pearson2$statistic,4),-0.0259)
  expect_equal(round(result_Pearson2$p.value,4),0.9794)
})

y <-round(y+seq(-0.1,0.1,length.out=150),4)
result_Spearman <- cortest(x,y,method="spearman")

# test_that("Spearman correlation value, statistic and pvalue are correct on deterministic (large) sample",{
#   expect_equal(round(result_Spearman$estimate,4),0.2216)
#   #expect_equal(round(result_Spearman$statistic,4),2.8494)
#   expect_true(abs(2.849407-result_Spearman$statistic)<=0.01)
#   expect_true(abs(0.004380078-result_Spearman$p.value)<=0.1)
#   #expect_equal(round(result_Spearman$p.value,4),0.0044)
# })

result_Kendall <- cortest(x,y,method="kendall")

test_that("Kendall correlation value, statistic and pvalue are correct on deterministic (large) sample",{
  expect_equal(round(result_Kendall$estimate,4),0.2064)
  expect_equal(round(result_Kendall$statistic,4),3.8891)
  expect_equal(round(result_Kendall$p.value,7),0.0001006)
})


# set.seed(1)
# x <- Evans$CHL[Evans$CDH==1]+runif(sum(Evans$CDH==1),0,0.01)
# set.seed(1)
# y <- Evans$DBP[Evans$CDH==1]+runif(sum(Evans$CDH==1),0,0.0001)
# result=cortest(x,y,method="spearman")
#
# test_that("Spearman returns the same as the R version",{
#   expect_warning(with(Evans,cortest(CHL[CDH==1],DBP[CDH==1],method="spearman")))
#   expect_equal(round(result$statistic,2),2.57)
#   expect_equal(round(result$p.value,2),0.01)
#   expect_equal(round(result$estimate,2),0.26)
# })


# set.seed(1)
# n<-4000
# x<-rnorm(n)
# y<-0.01*x^2+rnorm(n)
# #result=cortest(x,y)
# #result$statistic #-0.200384
# #result$p.value#0.8411905
# test_that("Pearson pvalue is correct on simulation 1",{
#   expect_equal(round(cortest(x,y)$p.value,2),0.84)
# })
#
#
#
# set.seed(1)
# n<-40
# x<-rnorm(n)
# y<-0.1*x^2+rnorm(n)
# # result=cortest(x,y)
# # result$statistic # 1.289596
# # result$p.value# 0.22207
#
# test_that("Pearson pvalue is correct on simulation 2",{
#   expect_equal(round(cortest(x,y)$p.value,2),0.22)
# })
#
# set.seed(1)
# n<-80
# x<-rnorm(n)
# y<-0.1*x^2+rnorm(n)
# # result=cortest(x,y)
# # result$statistic #-2.100994
# # result$p.value# 0.03901
#
# test_that("Pearson pvalue is correct on simulation 3",{
#   expect_equal(round(cortest(x,y)$p.value,2),0.04)
# })
#
# set.seed(1)
# n<-130
# x<-rnorm(n)
# y<-0.1*x^2+rnorm(n)
# # result=cortest(x,y)
# # result$statistic #1.161425
# # result$p.value# 0.2476303
#
# test_that("Pearson pvalue is correct on simulation 3",{
#   expect_equal(round(cortest(x,y)$p.value,2),0.25)
# })
#
# set.seed(1)
# n<-140
# x<-rnorm(n)
# y<-0.5*x^2+rnorm(n)
# # result=cortest(x,y)
# # result$statistic #0.3470794
# # result$p.value# 0.7290604
#
# test_that("Pearson pvalue is correct on simulation 3",{
#   expect_equal(round(cortest(x,y)$p.value,2),0.73)
# })
#
