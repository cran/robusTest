#' Bidimensional Empirical Cumulative Distribution Function
#'
#' Compute the empirical cumulative distribution function for a bivariate continuous distribution.
#' @param x,y the two continuous variables. Must be of same length.
#' @details The bidimensional e.c.d.f. (empirical cumulative distribution function) Fn is a step function with jumps i/n at observation values, where i is
#' the number of tied observations at that value.
#'
#' For observations (x1,y1), ..., (x_n,y_n), Fn is defined as
#' \deqn{Fn(t1,t2) = 1/n sum_{i=1}^n Indicator(xi<=t1,yi<=t2)}
#'
#' @return The result is returned as a matrix of dimension (n*n) where the entry (i,j) corresponds to Fn(xi,yj), i=1, ...,n,
#' j=1, ...,n.
#' @note Missing values are removed such that if a value of \code{x} (resp. \code{y}) is missing then the corresponding
#' values of both \code{x} and \code{y} are removed. The bidimensional e.c.d.f. is then computed on the remaining elements.
#' @seealso \code{\link{indeptest}}; the \code{bivariate} package also provides plots for the
#' bidimensional e.c.d.f.
#' @export
#' @examples
#' #Simulated data #1
#' x<-c(0.2, 0.3, 0.1, 0.4)
#' y<-c(0.5, 0.4, 0.05, 0.2)
#' ecdf2D(x,y)
#'
#' #Simulated data #2
#' n<-40
#' x<-rnorm(n)
#' y<-x^2+0.3*rnorm(n)
#' plot(x,y)
#' ecdf2D(x,y)

ecdf2D <- function(x,y)
{
  if (length(x)!=length(y)) stop("'x' and 'y' must have the same length")
  if (sum(is.na(x))!=0)
  {
    na.ind=which(is.na(x))
    x<-x[-na.ind];y<-y[-na.ind]
  }
  if (sum(is.na(y))!=0)
  {
    na.ind=which(is.na(y))
    x<-x[-na.ind];y<-y[-na.ind]
  }
  n<-length(x)
  sort_x <- order(x)
  sort_y <- order(y)
  result=matrix(0,n,n)
  for(i in 1:(n))
  {
    cpt=0
    for(j in 1:(n))
    {
      if(y[sort_x[j]] <= y[sort_y[i]]){cpt<-cpt+1}
      result[i,j] <-cpt/n
    }
  }
  #truc=list(result=result)
  #class(truc)<-"ecdf2D"
  #return(truc)
  return(list(ecdf=result,x=x,y=y))
}



