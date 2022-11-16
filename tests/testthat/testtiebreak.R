x <- c(1,2,2,3,4,5,5,5,7)
xbreak=tiebreak(x)

y <- c(4,9,12,11,2,10)
xy_break=tiebreak(x,y)


test_that("tiebreak works well (on fixed samples)",{
  expect_equal(sum(duplicated(xbreak)),0)
  expect_equal(sum(xy_break$x%in%xy_break$y),0)
  expect_equal(length(tiebreak(x,y,nb_break=TRUE)),3)
  expect_equal(length(tiebreak(x,nb_break=TRUE)),2)
})

x<-sample(1:10,size=30,replace=TRUE)
y<-sample(1:10,size=30,replace=TRUE)
res_tab=table(x)

test_that("tiebreak works well (on random samples)",{
  expect_equal(tiebreak(x,nb_break=TRUE)$nb_break,sum(res_tab[unlist(res_tab)!=1]-1))
  expect_equal(tiebreak(x,y,nb_break=TRUE)$nb_break,sum(x%in%y))
})

