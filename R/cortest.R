#' Calibrated tests for correlation between paired samples
#'
#' Tests the association/correlation for continuous paired samples using corrected versions of Pearson's, Kendall's and Spearman's correlation tests. These three tests are asymptotically well calibrated.
#' @param x,y the two continuous variables. Must be of same length.
#' @param alternative indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less".
#' @param method a character string indicating which test to implement. Can be \code{"pearson"}, \code{"kendall"} or \code{"spearman"}.
#' @param ties.break the method used to break ties in case there are ties in the x or y vectors. Can be \code{"none"} or \code{"random"}.
#' @param conf.level confidence level for the confidence interval of the correlation coefficient. It is used only for the Pearson's correlation test.
#' @details Three tests are implemented. The null hypothesis for the corrected Pearson test is: H0 Cor(X,Y)=0 where Cor represents Pearson's correlation coefficient.
#' The alternative is specified by the \code{alternative} argument. The null hypothesis for the corrected Kendall test is: H0 tau=0 where tau represents Kendall's tau coefficient.
#' The null hypothesis for the corrected Spearman test is: H0 rho=0 where rho represents Spearman's rho coefficient.
#' All tests are asymptotically well calibrated in the sense that the rejection probability under the null hypothesis is asymptotically equal to the level of the test.
#' For the Pearson test, the exact distribution of the test statistic under the Gaussian case has been tabulated for n<130. For n>=130, the Student distribution with n-2 degrees of
#' freedom is used.
#'
#' When Pearson's correlation test is used, a confidence interval for Pearson's correlation coefficient is also returned. This confidence interval has been implemented
#' from the delta-method. It should be noted that this method is asymptotic and can display very narrow intervals for small sample sizes and thus can suffer from
#' low coverage probabilities. We therefore recommend to use confidence intervals for Pearson's correlation coefficient only when n is at least larger than 100.
#'
#' The Kendall and Spearman correlation tests are not valid in the presence of ties in the \code{x} or \code{y} vector. If ties occur, they should be broken using the
#' \code{ties.break} option. Note that they can also be broken outside \code{cortest} using the function \code{tiebreak}.
#' @return Returns the result of the test with its corresponding p-value, the value of the test statistic and the estimated value of Pearson's correlation coefficient,
#' Kendall's tau or Spearman's rho. For Pearson's correlation test an asymptotic confidence interval for the correlation coefficient is also returned.
#' @note The option \code{ties.break} handles ties for both Kendall's and Spearman's test. If \code{ties.break="none"} the ties are ignored, if \code{ties.break="random"} they are randomly broken.
#' Note that only ties inside each vector are broken (but not ties between the two vectors).
#' @keywords test
#' @seealso \code{\link{vartest}}, \code{\link{indeptest}}, \code{\link{mediantest}}, \code{\link{wilcoxtest}}, \code{\link{tiebreak}}.
#' @importFrom stats cov pnorm pt qnorm qt runif sd var
#' @export
#' @examples
#' #Application on the Evans dataset
#' #Description of this dataset is available in the lbreg package
#' data(Evans)
#' with(Evans,cor.test(CHL[CDH==1],DBP[CDH==1]))
#' with(Evans,cortest(CHL[CDH==1],DBP[CDH==1]))
#' #The pvalues are very different!
#'
#' with(Evans,cortest(CHL[CDH==1],DBP[CDH==1],method="kendall",ties.break="random"))
#' with(Evans,cortest(CHL[CDH==1],DBP[CDH==1],method="spearman",ties.break="random"))
#'
#' #We use the function tiebreak to remove ties and compare the results from cor.test with cortest
#' X=tiebreak(Evans$CHL[Evans$CDH==1])
#' Y=tiebreak(Evans$DBP[Evans$CDH==1])
#' cor.test(X,Y,method="kendall")
#' cortest(X,Y,method="kendall")
#' cor.test(X,Y,method="spearman")
#' cortest(X,Y,method="spearman")
#'
#' \donttest{
#' #Simulated data
#' n=100 #sample size
#' M=10000 #number of replications
#' testone=function(n){
#' X=rnorm(n,0,1)
#' epsi=rnorm(n,0,1)
#' Y=X^2+0.3*epsi
#' list(test1=cortest(X,Y)$p.value,test2=cor.test(X,Y)$p.value) #cor.test is the standard Pearson test
#' }
#' res1=res2=rep(NA,M)
#' # Replications in order to check if the the corrected Pearson test and
#' # the standard test are well calibrated
#' for (i in 1:M)
#' {
#' result=testone(n)
#' res1[i]=result$test1
#' res2[i]=result$test2
#' }
#' mean(res1<0.05)
#' mean(res2<0.05)}
#'
#' \donttest{
#' #Replications with Kendall's test (may take a few minutes to run)
#' M=500
#' testone=function(n){
#' X=rnorm(n,0,1)
#' epsi=rnorm(n,0,1)
#' Y=X^2+0.3*epsi
#' list(test1=cortest(X,Y)$p.value,test2=cor.test(X,Y)$p.value,
#' test3=cortest(X,Y,method="kendall")$p.value,
#' test4=cor.test(X,Y,method="kendall")$p.value,
#' test5=cortest(X,Y,method="spearman")$p.value,
#' test6=cor.test(X,Y,method="spearman")$p.value)
#' #cor.test is the standard Pearson, Kendall or Spearman correlation test
#' }
#' res1=res2=res3=res4=res5=res6=rep(NA,M)
#' # Replications to check if the tests are well calibrated
#' for (i in 1:M)
#' {
#' result=testone(n)
#' res1[i]=result$test1
#' res2[i]=result$test2
#' res3[i]=result$test3
#' res4[i]=result$test4
#' res5[i]=result$test5
#' res6[i]=result$test6
#' }
#' mean(res1<0.05)
#' mean(res2<0.05)
#' mean(res3<0.05)
#' mean(res4<0.05)
#' mean(res5<0.05)
#' mean(res6<0.05)}

cortest <- function(x,y,alternative="two.sided",method="pearson",ties.break="none",conf.level=0.95) {UseMethod("cortest")}
#' @export
cortest.default=function(x,y,alternative="two.sided",method="pearson",ties.break="none",conf.level=0.95)#,ties.break="random"
{
  if (length(x)!=length(y)) stop("'x' and 'y' must have the same length")
  n <- length(x)
  if (n<=2) stop("lengths of 'x' and 'y'  must be greater than 2")
  Message=FALSE
  alpha=1-conf.level
  if (method=="pearson"){
    x<-(x-mean(x))/sd(x)
    y<-(y-mean(y))/sd(y)
    estimate=stats::cor(x,y)
    R <- x*y
    num=sum(R)
    deno <- sqrt(sum(R^2)-((sum(R))^2)/n)
    Tn <- num/deno
    #Construction of CI using the delta-method
    C11=var(x*y)
    C22=var(x^2)
    C33=var(y^2)
    C12=cov(x^2,x*y)
    C13=cov(y^2,x*y)
    C23=cov(y^2,x^2)
    varlim=C11+C22*(estimate^2)/4+C33*(estimate)^2/4-C12*estimate-C13*estimate+C23*estimate^2/2
    if (alternative=="two.sided" | alternative=="t"){
      Pval<-pval_pear_alt_two(Tn,n)
      #Pval <- 2*(1-pt(abs(Tn),n-2))
      CIl <- estimate-qt(1-alpha/2,n-2)*sqrt(varlim/n)
      CIr <- estimate+qt(1-alpha/2,n-2)*sqrt(varlim/n)
      if (CIl<(-1)) { CIl<-(-1)}
      if (CIr>1) { CIr<-1}
    }
    #CIl <- (-qt(1-alpha/2,n-2)*deno+num)/cor_denom #does not work!!!
    #CIr <- (qt(1-alpha/2,n-2)*deno+num)/cor_denom} #does not work!!!
    if (alternative=="less"| alternative=="l"){
      Pval<-pval_pear_alt_less(Tn,n)
      #Pval <- pt(Tn,n-2)
      CIl <- -1
      CIr <- estimate+qt(1-alpha,n-2)*sqrt(varlim/n)
      # CIr <- (qt(1-alpha,n-2)*deno+num)/cor_denom #does not work!!!
    }
    if (alternative=="greater"| alternative=="g"){
      Pval<-pval_pear_alt_great(Tn,n)
      #Pval <- 1-pt(Tn,n-2)
      CIl <- estimate-qt(1-alpha,n-2)*sqrt(varlim/n)
      #CIl <- (qt(alpha,n-2)*deno+num)/cor_denom #does not work!!!
      CIr <- 1}
  }
  if (method=="kendall"){
    duplix=duplicated(x)
    nb_duplix=sum(duplix)
    dupliy=duplicated(y)
    nb_dupliy=sum(dupliy)
    if((nb_dupliy+nb_duplix)!=0){
      #if ((length(x)!=length(unique(x)))|(length(Y)!=length(unique(Y)))) {
      if (ties.break=="none") {
        warning("The data contains ties! Use ties.break='random'")}
      if (ties.break=="random") {
        Message=TRUE
        if (nb_duplix!=0){
          x[duplix]=x[duplix]+runif(nb_duplix,-0.00001,0.00001)}
        if (nb_dupliy!=0){
          y[dupliy]=y[dupliy]+runif(nb_dupliy,-0.00001,0.00001)}
      }
      # xsort=sort(X,index.return=TRUE)
      # if (sum(diff(Xsort$x)==0)>0) {
      #   index=which(diff(Xsort$x)==0)#which value should be changed in the ordered sample
      #   X[Xsort$ix[index]]<-X[Xsort$ix[index]]+runif(length(Xsort$ix[index]),-0.00001,0.00001)
      # }
      # Ysort=sort(Y,index.return=TRUE)
      # if (sum(diff(Ysort$x)==0)>0) {
      #   index=which(diff(Ysort$x)==0)#which value should be changed in the ordered sample
      #   Y[Ysort$ix[index]]<-Y[Ysort$ix[index]]+runif(length(Ysort$ix[index]),-0.00001,0.00001)
      # }
      # #X=X+runif(length(X),-0.00001,0.00001) #break all values
      # #Y=Y+runif(length(Y),-0.00001,0.00001)
    }
    R <- array(0,dim=c(n,n))
    S <- array(0,dim=c(n,n))
    for(i in 1:n)
    {
      for(j in 1:n)
      {
        R[i,j]<-((x[j]-x[i])*(y[j]-y[i]))>0
        S[i,j]<-(x[j]>x[i]) & (y[j]>y[i])
      }
    }
    H <- apply(S,1,mean)+apply(S,2,mean)
    V <- var(H)
    estimate=2*(sum(R)/(n*(n-1))-0.5)
    Tn <- sqrt(n)*estimate/(4*sqrt(V))
    if (alternative=="two.sided" | alternative=="t"){
      Pval <- 2*(1-pnorm(abs(Tn)))
      CIl <- -qnorm(1-alpha/2)*(4*sqrt(V))/sqrt(n)+estimate
      CIr <- qnorm(1-alpha/2)*(4*sqrt(V))/sqrt(n)+estimate}
    if (alternative=="less"| alternative=="l"){
      Pval <- pnorm(Tn)
      CIl <- -1
      CIr <- qnorm(1-alpha)*(4*sqrt(V))/sqrt(n)+estimate}
    if (alternative=="greater"| alternative=="g"){
      Pval <- 1-pnorm(Tn)
      CIl <- qnorm(alpha)*(4*sqrt(V))/sqrt(n)+estimate
      CIr <- 1}
  }
  if (method=="spearman"){
    duplix=duplicated(x)
    nb_duplix=sum(duplix)
    dupliy=duplicated(y)
    nb_dupliy=sum(dupliy)
    if((nb_dupliy+nb_duplix)!=0){
      #if ((length(x)!=length(unique(x)))|(length(Y)!=length(unique(Y)))) {
      if (ties.break=="none") {
        warning("The data contains ties! Use ties.break='random'")}
      if (ties.break=="random") {
        Message=TRUE
        if (nb_duplix!=0){
          x[duplix]=x[duplix]+runif(nb_duplix,-0.00001,0.00001)}
        if (nb_dupliy!=0){
          y[dupliy]=y[dupliy]+runif(nb_dupliy,-0.00001,0.00001)}
      }
    }
    cppRet = spearmanCore(x, y);
    H = cppRet$H;
    sumR = cppRet$sumR;
    V<-16*var(H)
    estimate=3*sumR/(n^3-n)
    Tn<-sqrt(n)*sumR/(n*(n-1)*(n-2)*sqrt(V))
    if (alternative=="two.sided" | alternative=="t"){
      Pval <- 2*(1-pnorm(abs(Tn)))}
    if (alternative=="less"| alternative=="l"){
      Pval <- pnorm(Tn)}
    if (alternative=="greater"| alternative=="g"){
      Pval <- 1-pnorm(Tn)}
    CIl=CIr=NULL
  }
  result <- list(statistic=Tn, p.value=Pval,CI=c(CIl,CIr),conf.level=conf.level, estimate=estimate,alternative=alternative,method=method,Message=Message)
  class(result)<-"test"
  return(result)
}

#' @export
print.test <- function(x, ...)
{
  corval=x$estimate
  if (x$method=="pearson"){
    names(corval)<-"cor"
    cat("\nCorrected Pearson correlation test\n\n")
    if (round(x$p.value,4)==0){
      cat(paste("t = ", round(x$statistic,4), ", " , "p-value <1e-4","\n",sep= ""))
    } else {
    cat(paste("t = ", round(x$statistic,4), ", " , "p-value = ",round(x$p.value,4),"\n",sep= ""))}
    if (x$alternative=="two.sided" | x$alternative=="t"){
      cat("alternative hypothesis: true correlation is not equal to 0\n")
      #print(x$CI)
    }
    if (x$alternative=="less" | x$alternative=="l"){
      cat("alternative hypothesis: true correlation is less than 0\n")
    }
    if (x$alternative=="greater" | x$alternative=="g"){
      cat("alternative hypothesis: true correlation is greater than 0\n")
    }
    cat(paste(x$conf.level*100," % asymptotic confidence interval for the correlation coefficient:","\n",round(x$CI[1],4),"  ",round(x$CI[2],4),"\n",sep=""))
    cat("sample estimates:\n")
    print(corval)
  }
  if (x$method=="kendall"){
    names(corval)<-"tau"
    cat("\nCorrected Kendall correlation test\n\n")
    if (round(x$p.value,4)==0){
      cat(paste("t = ", round(x$statistic,4), ", " , "p-value <1e-4","\n",sep= ""))
    } else {
    cat(paste("t = ", round(x$statistic,4), ", " , "p-value = ",round(x$p.value,4),"\n",sep= ""))}
    if (x$alternative=="two.sided" | x$alternative=="t"){
      cat("alternative hypothesis: true tau is not equal to 0\n")}
    if (x$alternative=="less" | x$alternative=="l"){
      cat("alternative hypothesis: true tau is less than 0\n")}
    if (x$alternative=="greater" | x$alternative=="g"){
      cat("alternative hypothesis: true tau is greater than 0\n")}
    cat(paste(x$conf.level*100," % asymptotic confidence interval:","\n",round(x$CI[1],4),"  ",round(x$CI[2],4),"\n",sep=""))
    cat("sample estimates:\n")
    print(corval)
    if (x$Message==TRUE) {
      cat("\nTies were detected in the dataset and they were randomly broken")
    }
  }
  if (x$method=="spearman"){
    names(corval)<-"rho"
    cat("\nCorrected Spearman correlation test\n\n")
    if (round(x$p.value,4)==0){
      cat(paste("t = ", round(x$statistic,4), ", " , "p-value <1e-4","\n",sep= ""))
    } else {
    cat(paste("S = ", round(x$statistic,4), ", " , "p-value = ",round(x$p.value,4),"\n",sep= ""))}
    if (x$alternative=="two.sided" | x$alternative=="t"){
      cat("alternative hypothesis: true rho is not equal to 0\n")}
    if (x$alternative=="less" | x$alternative=="l"){
      cat("alternative hypothesis: true rho is less than 0\n")}
    if (x$alternative=="greater" | x$alternative=="g"){
      cat("alternative hypothesis: true rho is greater than 0\n")}
    cat("sample estimates:\n")
    print(corval)
    if (x$Message==TRUE) {
      cat("\nTies were detected in the dataset and they were randomly broken")
    }
  }
}

pval_pear_alt_two<-function(Tn,n)
{
  if (n<=129){
    y1<-(1:(2e5))/(2e5)
    y1<-y1[c(seq(10,40000,by=40),seq(40000,160000,by=600),seq(160040,2e5,by=40))]
    #y1<-y1[seq(1,2e5,by=6)]
    x1<-robust_Pearson_table[[n]]
    funstep<-stats::stepfun(x1,c(0,y1))
    Pval<-1-funstep(abs(Tn))+funstep(-abs(Tn))
    #2*(1-funstep(abs(Tn)))
  } else {
    Pval <- 2*(1-pt(abs(Tn),n-2))
  }
  return(Pval)
}

pval_pear_alt_less<-function(Tn,n)
{
  if (n<=129){
    y1<-(1:(2e5))/(2e5)
    y1<-y1[c(seq(10,40000,by=40),seq(40000,160000,by=600),seq(160040,2e5,by=40))]
    #y1<-y1[seq(1,2e5,by=6)]
    x1<-robust_Pearson_table[[n]]
    funstep<-stats::stepfun(x1,c(0,y1))
    Pval<-funstep(Tn)
  } else {
    Pval <- pt(Tn,n-2)
  }
  return(Pval)
}

pval_pear_alt_great<-function(Tn,n)
{
  if (n<=129){
    y1<-(1:(2e5))/(2e5)
    y1<-y1[c(seq(10,40000,by=40),seq(40000,160000,by=600),seq(160040,2e5,by=40))]
    #y1<-y1[seq(1,2e5,by=6)]
    x1<-robust_Pearson_table[[n]]
    funstep<-stats::stepfun(x1,c(0,y1))
    Pval<-1-funstep(Tn)
  } else {
    Pval <- 1-pt(Tn,n-2)
  }
  return(Pval)
}
