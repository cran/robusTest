data_rec=read.table("Recensement.txt",header=TRUE)
dataM=data_rec[(data_rec$SEXE=="M"),]

data3=data_rec[data_rec$CATEGORIE %in% c(3,4,5,10),]

test_that("vartest on Recensement data",{
  expect_equal(round(vartest(dataM$SAL_HOR~dataM$SYNDICAT)$p.value,4),0.0107)
  expect_equal(round(vartest(data3$SAL_HOR~data3$CATEGORIE)$p.value,3),0.041)
})


x <- c(1:150)*2
y <- rep(c(-0.95,  0.73,  1.79, -0.64, -0.50),30)
y <-y+seq(-0.1,0.1,length.out=150)
xbis <- x/sd(x)
result <- vartest(xbis,y)

result_stat <- round(result$statistic,4) #carefull, statistic has name "t"
names(result_stat) <- NULL

test_that("vartest statistic and pvalue are correct on deterministic (large) sample",{
  expect_equal(result_stat,-0.5778)
  expect_equal(round(result$p.value,4),0.5639)
})
