#' Test for the median in the one and two-sample paired tests
#'
#' In the one sample case, test if the median of the random variable is equal to 0. In the paired two-sample case, test if the median of the
#' difference between the two random variables is equal to 0. The function also returns the confidence interval for the median of the
#' random variable in the one-sample case or for the median of the difference between the two random variables in the two-sample case.
#' @param x,y two continuous variables.
#' @param alternative indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less".
#' @param paired a logical value. If \code{paired=TRUE}, the user must provide values for \code{x} and \code{y}
#' and the paired test is implemented. If \code{paired=FALSE}, only \code{x} must be provided.
#' @param conf.level confidence level for the confidence interval of the median.
#' @details The null hypothesis for the one sample median test is: H0 Med(X)=0 where Med represents the median.
#' The alternative is specified by the \code{alternative} argument: "\code{greater}" means that Med(X)>0, "\code{less}"
#'  means that Med(X)<0 and "\code{two.sided}" means that Med(X) is different from 0. The null hypothesis for the paired median test is: H0 Med(X-Y)=0. Both tests are asymptotically
#'  calibrated in the sense that the rejection probability under the null hypothesis is asymptotically equal to the level of the test. The
#'  test and the confidence interval are constructed based on asymptotic results on the rank statistics for the uniform distribution, as described in the book
#'  \emph{Asymptotic Statistics} from A. W. Van der Vaart, 1998.
#' @note The paired median test can be implemented either by providing the variables \code{x} and \code{y} with \code{paired=TRUE} or by just providing
#' one vector equal to the difference between \code{x} and \code{y} with \code{paired=FALSE}.
#' @return Returns the result of the test with its corresponding p-value and the value of median of X or of X-Y. The confidence interval is also displayed
#' but only when \code{alternative=two.sided}.
#' @keywords test
#' @seealso \code{\link{cortest}}, \code{\link{indeptest}}, \code{\link{vartest}}, \code{\link{wilcoxtest}}.
#' @author See \emph{Asymptotic Statistics}.
#' A. W. Van der Vaart, 1998.
#' @export
#' @examples
#' #Simulations
#' x <- rexp(200)-log(2)
#' mediantest(x)
#' x <- rexp(200)-log(2)+0.2
#' mediantest(x,alternative="greater")
#'
#' n=100 #sample size
#' M=1000 #number of replications
#' res1=res2=rep(NA,M)
#' testone=function(n){
#' D=rchisq(n,df=4)-qchisq(df=4, p=0.5)
#' list(test1=mediantest(D)$p.value,test2=binom.test(sum(D>0),n)$p.value)
#' } #test2 is the sign test.
#' for (i in 1:M)
#' {
#' result=testone(n)
#' res1[i]=result$test1
#' res2[i]=result$test2
#' }
#' mean(res1<0.05) #0.048
#' mean(res2<0.05) # 0.04

mediantest <- function(x,y=NULL,alternative="two.sided",paired=FALSE,conf.level=0.95) {UseMethod("mediantest")}
#' @export
mediantest.default=function(x,y=NULL,alternative="two.sided",paired=FALSE,conf.level=0.95)
{
  alpha=1-conf.level
  if (paired==TRUE)
  {
    if (is.null(y)){ stop("'y' is missing for paired test")}
    if (is.null(x)){ stop("'x' is missing for paired test")}
    if (length(x)!=length(y)){stop("'x' and 'y' should have the same size for the paired median test")}
    #Perform the paired two sample test
    X <- x-y
    medX <- stats::median(X)
    pairedTest<-TRUE
  }
  if (paired==FALSE)
  {
    if (is.null(x)) stop("'x' is missing")
    if (is.null(y)==FALSE) stop("there should be no 'y' for the one sample median test")
    X<-x
    medX <- stats::median(X)
    pairedTest<-FALSE
  }
  n <- length(X)
  y<-stats::ecdf(X)
  k=n*y(0)+1
  xsort <- sort(x)
  CIl <- xsort[ceiling(-stats::qnorm(1-alpha/2)*sqrt(n)/2+n/2)]
  CIr <- xsort[ceiling(stats::qnorm(1-alpha/2)*sqrt(n)/2+n/2)]
  if (alternative=="two.sided" | alternative=="t"){
    alpha1<-2*(1-stats::pnorm((2*k-n-0.5)/sqrt(n)))
    alpha2<-2*(1-stats::pnorm((n-2*k+0.5)/sqrt(n)))
    Pval <- min(alpha1, alpha2)}
  if (alternative=="less"| alternative=="l"){
    Pval <- 1-stats::pnorm((2*k-n-1.8)/sqrt(n))}
  if (alternative=="greater"| alternative=="g"){
    Pval <- 1-stats::pnorm((n-2*k+2.1)/sqrt(n))}
  result <- list(p.value=Pval,CI=c(CIl,CIr),conf.level=conf.level,estimate=medX, alternative=alternative,pairedTest=pairedTest)
  class(result)<-"mediantest"
  return(result)
}

#' @export
print.mediantest <- function(x, ...)
{
  medval=x$estimate
  if (x$pairedTest==FALSE){
    cat("\nMedian test for the on sample case\n\n")
    cat(paste("p-value = ",round(x$p.value,4),"\n",sep= ""))
    if (x$alternative=="two.sided" | x$alternative=="t"){
      cat("alternative hypothesis: median (X) is not equal to zero\n")
      cat(paste(x$conf.level*100," % asymptotic confidence interval for the median:","\n",round(x$CI[1],4),"  ",round(x$CI[2],4),"\n",sep=""))
      }
    if (x$alternative=="less" | x$alternative=="l"){
      cat("alternative hypothesis: median (X) is negative\n")}
    if (x$alternative=="greater" | x$alternative=="g"){
      cat("alternative hypothesis: median (X) is positive\n")}
    cat("median estimate:\n")
    print(medval)
  }
  if (x$pairedTest==TRUE){
    cat("\nMedian test for the two sample case\n\n")
    cat(paste("p-value = ",round(x$p.value,4),"\n",sep= ""))
    if (x$alternative=="two.sided" | x$alternative=="t"){
      cat("alternative hypothesis: median (X-Y) is not equal to zero\n")
      cat(paste(x$conf.level*100," % asymptotic confidence interval for median (X-Y):","\n",round(x$CI[1],4),"  ",round(x$CI[2],4),"\n",sep=""))
      }
    if (x$alternative=="less" | x$alternative=="l"){
      cat("alternative hypothesis: median (X-Y) is negative\n")}
    if (x$alternative=="greater" | x$alternative=="g"){
      cat("alternative hypothesis: median (X-Y) is positive\n")}
    cat("median estimate of the difference:\n")
    print(medval)
  }
}




