#' Calibrated Wilcoxon rank sum and signed rank tests
#'
#' Compares the distribution between two random variables by testing if one variable tends to take larger (or smaller) values than the other. The test
#' works for independent and paired variables by using corrected versions of the Wilcoxon (or equivalently Mann-Whitney in the two-sample case) for one and two-sample tests.
#' @param x,y two continuous variables.
#' @param alternative indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less".
#' @param ties.break the method used to break ties in case there are ties in the x or y vectors. Can be \code{"none"} or \code{"random"}.
#' @param paired a logical value. If \code{paired=TRUE}, you must provide values for \code{x} and \code{y} (of same length)
#' and the paired test is implemented. If \code{paired=FALSE}, the paired test is implemented when \code{y} is null and
#' only \code{x} is provided and the two sample test (for independent variables) is implemented when both \code{x} and \code{y} are provided.
#' @param ... it is possible to use a formula with or without specifying a dataset
#' from the commands \code{wilcoxtest(x~y)} or \code{wilcoxtest(x~y,data=dataset)} with dataset the name of the data.
#' @details For two independent samples, the null hypothesis for the corrected Wilcoxon (Mann-Whitney) test is: H0 Med(X-Y)=0 where Med represents the median.
#' The alternative is specified by the \code{alternative} argument: "\code{greater}" means that Med(X-Y)>0 and "\code{less}"
#'  means that Med(X-Y)<0. The null hypothesis for the paired Wilcoxon test is: H0 Med(D1+D2)=0 where D1 is the difference
#'  between X1 and Y1 taken on the same pair (same with D2 on a different pair). Both tests are asymptotically well calibrated in the sense that the rejection probability under the
#'  null hypothesis is asymptotically equal to the level of the test.
#' @note The function can also be called using formulas: type \code{wilcoxtest(x~y,data)} with x the quantitative variable
#' and y a factor variable with two levels. The option \code{ties.break} handles ties in the Wilcoxon test. If \code{ties.break="none"} the ties are ignored, if
#' \code{ties.break="random"} they are randomly broken. For the Wilcoxon rank sum test the ties between the \code{x} and \code{y} are
#' detected and broken (but the ties inside the \code{x} and \code{y} vectors are not changed). For the signed rank test, the ties in the
#' vector \code{x-y} (or in the \code{x} vector in case \code{y=NULL}) are randomly broken.
#' @return Returns the result of the test with its corresponding p-value and the value of the test statistic.
#' @keywords test
#' @seealso \code{\link{cortest}}, \code{\link{indeptest}}, \code{\link{mediantest}}, \code{\link{vartest}}.
#' @export
#' @examples
#' #Application on the Evans dataset
#' data(Evans)
#' #Description of this dataset is available in the lbreg package
#' with(Evans,wilcox.test(CHL[CDH==0],CHL[CDH==1]))
#' with(Evans,wilcoxtest(CHL[CDH==0],CHL[CDH==1]))
#' wilcoxtest(CHL~CDH,data=Evans) #using formulas
#' wilcoxtest(CHL~CDH,data=Evans,ties.break="random")
#' #the same test where ties are randomly broken
#'
#' \donttest{
#' #For independent samples
#' n=100 #sample size
#' M=1000 #number of replications
#' testone=function(n){
#' X=runif(n,-0.5,0.5)
#' Y=rnorm(3*n,0,0.04)
#' list(test1=wilcoxtest(X,Y)$p.value,test2=wilcox.test(X,Y)$p.value)
#' #wilcox.test is the standard Wilcoxon test
#' }
#'
#' #Simulation under the null hypothesis
#' #(note that P(X>Y)=0.5)
#' #Takes a few seconds to run
#' res1=res2=rep(NA,M)
#' for (i in 1:M)
#' {
#' result=testone(n)
#' res1[i]=result$test1
#' res2[i]=result$test2
#' }
#' mean(res1<0.05)
#' mean(res2<0.05)}
#' \donttest{
#' #For paired samples
#' #We use the value of the median of a Gamma distributed variable with shape
#' #parameter equal to 1/5 and scale parameter equal to 1. This value is
#' #computed from the command qgamma(shape=1/5, scale=1, 0.5)
#' n=100 #sample size
#' M=1000 #number of replications
#' testone=function(n){
#' D=rgamma(n,shape=1/10,scale=1)-qgamma(shape=1/5, scale=1, 0.5)/2
#' list(test1=wilcoxtest(D,ties.break = "random")$p.value,test2=wilcox.test(D)$p.value)
#' #wilcox.test is the standard paired Wilcoxon test
#' }
#' #Simulation under the null hypothesis
#' #(note that Med(D_1+D_2)=0)
#' #Takes a few seconds to run
#' for (i in 1:M)
#' {
#' result=testone(n)
#' res1[i]=result$test1
#' res2[i]=result$test2
#' }
#' mean(res1<0.05)
#' mean(res2<0.05)}

wilcoxtest <- function(x,y=NULL,alternative="two.sided",ties.break="none",paired=FALSE,...) {UseMethod("wilcoxtest")}
#' @export
wilcoxtest.default=function(x,y=NULL,alternative="two.sided",ties.break="none",paired=FALSE,...)
{
  Message=FALSE
  if (paired==TRUE)
  {
    if (is.null(y)){ stop("'y' is missing for paired test")}
    if (is.null(x)){ stop("'x' is missing for paired test")}
    #Perform the paired two sample test
    X <- x-y
    n <- length(X)
    dupliX=duplicated(X)
    nb_dupliX=sum(dupliX)
    if (nb_dupliX!=0){
      #if ((length(X)!=length(unique(X)))){
      if (ties.break=="none") {
        warning("The data contains ties! Use ties.break='random'")}
      if (ties.break=="random") {
        X[dupliX]=X[dupliX]+runif(nb_dupliX,-0.00001,0.00001)
        #Xsort=sort(X,index.return=TRUE)
        #if (sum(diff(Xsort$x)==0)>0) {
        #  index=which(diff(Xsort$x)==0)#which value should be changed in the ordered sample
        #  X[Xsort$ix[index]]<-X[Xsort$ix[index]]+runif(length(Xsort$ix[index]),-0.00001,0.00001)
        #}
        Message=TRUE
      }
    }
    R <- array(0,dim=c(n,n))
    Diag <- vector(mode = "numeric", length = n)
    for(i in 1:n)
    {
      for(j in 1:n)
      {
        R[i,j]<-(X[j]+X[i]>0)
        Diag[i]=R[i,i]
      }
    }
    H <- apply(R,1,mean)
    V <- var(H)
    Tn <- sqrt(n)*((sum(R)-sum(Diag))/(n*(n-1))-0.5)/(2*sqrt(V))
    pairedTest<-TRUE
  }
  if (paired==FALSE)
  {
    if (is.null(x)) stop("'x' is missing")
    if (is.null(y))
    {
      #Perform the paired two sample test with X being the difference between the variables in the same pair
      n <- length(x)
      #if ((length(X)!=length(unique(X)))){
      dupliX=duplicated(x)
      nb_dupliX=sum(dupliX)
      if (nb_dupliX!=0){
        if (ties.break=="none") {
          warning("The data contains ties!")}
        if (ties.break=="random") {
          x[dupliX]=x[dupliX]+runif(nb_dupliX,-0.00001,0.00001)
          #Xsort=sort(X,index.return=TRUE)
          #if (sum(diff(Xsort$x)==0)>0) {
          #  index=which(diff(Xsort$x)==0)#which value should be changed in the ordered sample
          #  X[Xsort$ix[index]]<-X[Xsort$ix[index]]+runif(length(Xsort$ix[index]),-0.00001,0.00001)
          #}
          Message=TRUE
        }
      }
      R <- array(0,dim=c(n,n))
      Diag <- vector(mode = "numeric", length = n)
      for(i in 1:n)
      {
        for(j in 1:n)
        {
          R[i,j]<-(x[j]+x[i]>0)
          Diag[i]=R[i,i]
        }
      }
      H <- apply(R,1,mean)
      V <- var(H)
      Tn <- sqrt(n)*((sum(R)-sum(Diag))/(n*(n-1))-0.5)/(2*sqrt(V))
      pairedTest<-TRUE
    } else {
      #Perform the two sample test with X and Y being two independent variables
      n <- length(x)
      m <- length(y)
      ties=x%in%y
      if (sum(ties)!=0){
        if (ties.break=="none") {
          warning("The data contains ties between the two vectors! Use ties.break='random'")}
        if (ties.break=="random") {
          x[ties] <- x[ties]+runif(sum(ties),-0.00001,0.00001)
          Message=TRUE
        }
      }
      R <- array(0,dim=c(n,m))
      for(i in 1:n)
      {
        for(j in 1:m)
        {
          R[i,j]<-(y[j]>x[i])
        }
      }
      H <- apply(R,1,mean)
      G <- apply(R,2,mean)
      V <- var(H)/n + var(G)/m
      Tn <- (mean(R)-0.5)/sqrt(V)
      pairedTest<-FALSE
    }
  }
  if (alternative=="two.sided" | alternative=="t"){
    Pval <- 2*(1-pnorm(abs(Tn)))}
  if (alternative=="less"| alternative=="l"){
    Pval <- pnorm(Tn)}
  if (alternative=="greater"| alternative=="g"){
    Pval <- 1-pnorm(Tn)}
  result <- list(statistic=Tn, p.value=Pval, alternative=alternative,pairedTest=pairedTest,Message=Message)
  class(result)<-"testW"
  return(result)
}

#' @export
wilcoxtest.formula=function(x,y=NULL,alternative="two.sided",ties.break="none",paired=FALSE,data=list(),...)
  #wilcoxtest.formula=function(formula,data=list(),alternative="two.sided",ties.break="none",paired=FALSE)
{
  formula<-x
  mf <- stats::model.frame(formula=formula, data=data)
  response <- attr(attr(mf, "terms"), "response")
  Fact<-factor(mf[[-response]])
  DATA <- stats::setNames(split(mf[[response]], Fact), c("x", "y"))
  result <- wilcoxtest.default(DATA[[1]],DATA[[2]],alternative=alternative,ties.break = ties.break,paired=paired)
  return(result)
}

#' @export
print.testW <- function(x, ...)
{
  if (x$pairedTest==FALSE){
    cat("\nCorrected Wilcoxon rank sum test\n\n")
    cat(paste("W = ", round(x$statistic,4), ", " , "p-value = ",round(x$p.value,4),"\n",sep= ""))
    if (x$alternative=="two.sided" | x$alternative=="t"){
      cat("alternative hypothesis: median (X-Y) is not equal to zero\n")}#X and Y tend to take different values
    if (x$alternative=="less" | x$alternative=="l"){
      cat("alternative hypothesis: median (X-Y) is negative\n")}# X tends to be smaller than Y
    if (x$alternative=="greater" | x$alternative=="g"){
      cat("alternative hypothesis: median (X-Y) is positive\n")}}
  if (x$pairedTest==TRUE){
    cat("\nCorrected Wilcoxon signed rank test\n\n")
    cat(paste("W = ", round(x$statistic,4), ", " , "p-value = ",round(x$p.value,4),"\n",sep= ""))
    if (x$alternative=="two.sided" | x$alternative=="t"){
      cat("alternative hypothesis: median (D1+D2) is not equal to zero\n")}#X and Y tend to take different values
    if (x$alternative=="less" | x$alternative=="l"){
      cat("alternative hypothesis: median (D1+D2) is negative\n")}# X tends to be smaller than Y
    if (x$alternative=="greater" | x$alternative=="g"){
      cat("alternative hypothesis: median (D1+D2) is positive\n")}}#X tends to be larger than Y
  if (x$Message==TRUE) {
    cat("\nTies were detected in the dataset and they were randomly broken")
  }
}




