#' Calibrated variance test between two or more independent samples
#'
#' Tests the equality of variance for continuous independent samples using Welch's ANOVA on the square of the centered variables.
#' @param x,y the two continuous variables for the two samples variance test.
#' @param ... other parameters such as \code{alternative} which indicates the alternative hypothesis
#' as one of "two.sided", "greater" or "less", or \code{data} and \code{formula} to include a dataset with a formula
#' from the command \code{vartest(x~y,data=dataset)} with dataset the name of the data.
#' @details For the two sample variance test, the null hypothesis is: H0 var(X)=var(Y) where var represents the variance.
#' The alternative is specified by the \code{alternative} argument.
#' The test uses Welch's ANOVA procedure on the variables \eqn{(x_i-mean(x))^2}, \eqn{(y_i-mean(y))^2} and is asymptotically well calibrated in the sense that the rejection probability under the null hypothesis is asymptotically equal to the level of the test.
#' The test can be applied to more than two groups in which case it tests if all groups have the same variance.
#' @return Returns the result of the test with its corresponding p-value and the value of the test statistic. In the two sample case, the function
#' also returns the confidence interval for the difference of the two variances.
#' @note The function can also be called using formulas: type \code{vartest(x~y,data)} with x the quantitative variable
#' and y a factor variable with two or more levels. For more than two groups, only the formula is valid.
#' @keywords test
#' @seealso \code{\link{cortest}}, \code{\link{indeptest}}, \code{\link{mediantest}}, \code{\link{wilcoxtest}}.
#' @export
#' @examples
#' #Application on the Evans dataset
#' data(Evans)
#' #Description of this dataset is available in the lbreg package
#' with(Evans,var.test(CHL[CDH==0],CHL[CDH==1]))
#' with(Evans,vartest(CHL[CDH==0],CHL[CDH==1]))
#' vartest(CHL~CDH,data=Evans) #using formulas
#'
#' #Similar pvalues between var.test and vartest
#'
#' \donttest{
#' #Simulated data
#' n=60 #sample size
#' M=10000 #number of replications
#' testone=function(n){
#' X=rnorm(n,0,1)
#' Y=rchisq(2*n,df=2)/2
#' list(test1=vartest(X,Y)$p.value,test2=var.test(X,Y)$p.value) #var.test is the standard Fisher test
#' }
#' res1=res2=rep(NA,M)
#' # Replications to check if the corrected Fisher test and the standard test are well calibrated
#' for (i in 1:M)
#' {
#' result=testone(n)
#' res1[i]=result$test1
#' res2[i]=result$test2
#' }
#' mean(res1<0.05)  #0.0515
#' mean(res2<0.05)  #0.1509}
#'

vartest <- function(x,y,...) {UseMethod("vartest")}#alternative="two.sided",conf.level=0.95,data=list(),...
#' @export
vartest.default=function(x,y,alternative="two.sided",conf.level=0.95,...)
{
  X2 <- (x-mean(x))^2
  Y2 <- (y-mean(y))^2
  ttest<-stats::t.test(X2,Y2,conf.level = conf.level)
  result <- list(statistic=ttest$statistic, p.value=ttest$p.value,CI=c(ttest$conf.int),conf.level=conf.level, estimate=ttest$estimate,alternative=alternative)
  class(result)<-"vartest"
  return(result)
}

#' @export
vartest.formula <- function(x,y,data=list(),...)#,alternative="two.sided"
{
  formula<-x
  mf <- stats::model.frame(formula=formula, data=data)
  #Fact<-as.factor(Fact)
  #reg<-lm(model.response(mf)~Fact)
  Fact<-mf[,2]
  reg<-stats::lm(formula,data=data)
  vtest<-stats::oneway.test((reg$residuals)^2~Fact)
  #return(oneway.test((reg$residuals)^2~Fact)$p.value)
  result <- list(p.value=vtest$p.value)
  class(result)<-"vartestform"
  return(result)
}
#' @export
print.vartest <- function(x, ...)
{
  varval=c(round(x$estimate[1],6),round(x$estimate[2],6))
  names(varval)<-c("variance of x", "variance of y") #ratio of variances"
  cat("\nCorrected F test to compare two variances\n\n")
  cat(paste("F = ", round(x$statistic,4), ", " , "p-value = ",round(x$p.value,4),"\n",sep= ""))
  if (x$alternative=="two.sided" | x$alternative=="t"){
    cat("alternative hypothesis: true difference of variances is not equal to 0\n")}#ratio
  if (x$alternative=="less" | x$alternative=="l"){
    cat("alternative hypothesis: true difference of variances is less than 0\n")}
  if (x$alternative=="greater" | x$alternative=="g"){
    cat("alternative hypothesis: true difference of variances is greater than 0\n")}
  cat(paste(x$conf.level*100," % percent confidence interval:","\n",round(x$CI[1],4),"  ",round(x$CI[2],4),"\n",sep=""))
  cat("sample estimates:\n")
  #cat("ratio of variances\n")
  print(varval)
}

#' @export
print.vartestform <- function(x, ...)
{
  cat("\nCorrected F test to compare variances\n\n")
  cat(paste("p-value = ",round(x$p.value,4),"\n",sep= ""))
  cat("alternative hypothesis: all the variances are not equal\n")
}

