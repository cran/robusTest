#' Robust independence test for two continuous variables of Kolmogorov-Smirnov's type
#'
#' Test the independence between two continuous variables based on the maximum distance
#' between the joint empirical cumulative distribution function and the product of the marginal
#' empirical cumulative distribution functions.
#' @param x,y the two continuous variables. Must be of same length.
#' @param N the number of Monte-Carlo replications if simu=TRUE.
#' @param simu if TRUE a Monte-Carlo simulation with \code{N} replications is used to determine the
#' distribution of the test statistic under the null hypothesis. If FALSE, pre computed tables are used (see Details
#' for more information).
#' @param ties.break the method used to break ties in case there are ties in the x or y vectors. Can be \code{"none"},
#' \code{"random"} or \code{"rep_random"}.
#' @param nb_tiebreak the number of repetition for breaking the ties when \code{ties.break="rep_random"}.
#' @details For two continuous variables, \code{indeptest} tests H0 X and Y are independent
#' against H1 X and Y are not independent.
#'
#' For observations (x1,y1), ..., (x_n,y_n), the bivariate e.c.d.f.
#' (empirical cumulative distribution function) Fn is defined as:
#' \deqn{Fn(t1,t2) = sum_{i=1}^n Indicator(xi<=t1,yi<=t2)/n.}
#'
#' Let Fn(t1) and Fn(t2) be the marginals e.c.d.f. The test statistic is defined as:
#' \deqn{n^(1/2) sup_{t1,t2} |Fn(t1,t2)-Fn(t1)*Fn(t2)|.}
#'
#'Under H0 the test statistic is distribution free and is equivalent to
#'the same test statistic computed for two independent continuous uniform variables in \eqn{[0,1]},
#'where the supremum is taken for t1,t2 in \eqn{[0,1]}. Using this result, the distribution of the test
#'statistic is obtained using Monte-Carlo simulations. The user can either use the argument simu=TRUE to
#'perform the Monte-Carlo simulation (with N the number of replications) or simply use the available tables
#'by choosing simu=FALSE. In the latter case, the exact distribution is estimated for n=1, ...,150. For \eqn{151<=n<=175}, the
#'distribution with n=150 is used. For \eqn{176<=n<=250}, the distribution with n=200 is used.
#'For \eqn{251<=n<=400}, the distribution with n=300 is used. For \eqn{401<=n<=750}, the distribution with n=500 is used.
#'For \eqn{n>=751}, the distribution with n=1000 is used. Those tables were computed using 2e^5 replications in Monte-Carlo simulations.
#' @return Returns the result of the test with its corresponding p-value and the value of the test statistic.
#' @note Only a two sided alternative is possible with this test. Missing values are removed such that if a value
#' of \code{x} (resp. \code{y}) is missing then the corresponding
#' values of both \code{x} and \code{y} are removed. The test is then implemented on the remaining elements. If \code{ties.break="none"} the ties are ignored, putting
#' mass (number of ties)/n at tied observations in the computation of the empirical cumulative distribution functions.
#' If \code{ties.break="random"} they are randomly broken. If \code{ties.break="rep_random"} they are randomly broken \code{nb_tiebreak}
#' times where \code{nb_tiebreak} is a parameter of the function. In that case, the test statistic and the p values are computed by taking
#' the average over all replications.
#'
#' This function is implemented using the Rcpp package.
#' @keywords test
#' @seealso \code{\link{cortest}}, \code{\link{vartest}}, \code{\link{mediantest}}, \code{\link{wilcoxtest}}.
#' See also the \code{hoeffd} function in the \code{Hmisc} package for the Hoeffding test.
#' @author See \emph{Distribution Free Tests of Independence Based on the Sample Distribution Function}.
#' J. R. Blum, J. Kiefer and M. Rosenblatt, 1961.
#' @export
#' @examples
#' #Simulated data 1
#' x<-c(0.2, 0.3, 0.1, 0.4)
#' y<-c(0.5, 0.4, 0.05, 0.2)
#' indeptest(x,y)
#'
#' #Simulated data 2
#' n<-40 #sample size
#' x<-rnorm(n)
#' y<-x^2+0.3*rnorm(n)
#' plot(x,y)
#' indeptest(x,y)
#'
#' #Application on the Evans dataset
#' #Description of this dataset is available in the lbreg package
#' data(Evans)
#' with(Evans,plot(CHL[CDH==1],DBP[CDH==1]))
#' with(Evans,cor.test(CHL[CDH==1],DBP[CDH==1])) #the standard Pearson test
#' with(Evans,cortest(CHL[CDH==1],DBP[CDH==1])) #the robust Pearson test
#' with(Evans,indeptest(CHL[CDH==1],DBP[CDH==1])) #the robust independence test
#' #The robust tests give very different pvalues than the standard Pearson test!
#'
#' #Breaking the ties
#' #The ties are broken once
#' with(Evans,indeptest(CHL[CDH==1],DBP[CDH==1],ties.break="random"))
#' #The ties are broken repeatedly and the average of the test statistics and p.values
#' #are computed
#' with(Evans,indeptest(CHL[CDH==1],DBP[CDH==1],ties.break="rep_random",nb_tiebreak=100))

#' @export
indeptest <- function(x,y,N=50000,simu=FALSE,ties.break="none",nb_tiebreak=100) {UseMethod("indeptest")}
#' @export
indeptest.default<-function (x,y,N=50000,simu=FALSE,ties.break="none",nb_tiebreak=100){
  if (length(x)!=length(y)) stop("'x' and 'y' must have the same length")
  Message=FALSE
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
  n <- length(x)
  if (n<3) stop("lengths of 'x' and 'y' must be greater than 2")
  if (ties.break=="rep_random"){
    dupliX=duplicated(x)
    nb_dupliX=sum(dupliX)
    dupliY=duplicated(y)
    nb_dupliY=sum(dupliY)
    ties=x%in%y
    if ((nb_dupliX+nb_dupliY)!=0 | sum(ties)!=0){
      Message=TRUE
      pval_vect<-stat_vect<-rep(NA,nb_tiebreak)
      for (j in 1:nb_tiebreak)
      {
        newx=x
        newy=y
        if (nb_dupliX!=0){
          newx=x[dupliX]+runif(nb_dupliX,-0.00001,0.00001)
          newy=y}
        if (nb_dupliY!=0){
          newy[dupliY]=y[dupliY]+runif(nb_dupliY,-0.00001,0.00001)
          newx=x}
        if (sum(ties)!=0){
          newx[ties] <- x[ties]+runif(sum(ties),-0.00001,0.00001)
          newy=y}
        result=pval_comput(newx,newy,simu=simu,N=N)
        pval_vect[j]<-result$Pval
        stat_vect[j]<-result$Tn
        Pval=mean(pval_vect)
        Tn=mean(stat_vect)
      }
    } else {
      Pval<-pval_comput(x,y,simu=simu,N=N)
      }
    } else {
  dupliX=duplicated(x)
  nb_dupliX=sum(dupliX)
  dupliY=duplicated(y)
  nb_dupliY=sum(dupliY)
  ties=x%in%y
  if ((nb_dupliX+nb_dupliY)!=0 | sum(ties)!=0){
    if (ties.break=="none") {
      warning("The data contains ties! Use ties.break='random'")}
    if (ties.break=="random") {
      Message=TRUE
      if (nb_dupliX!=0){
        x[dupliX]=x[dupliX]+runif(nb_dupliX,-0.00001,0.00001)}
      if (nb_dupliY!=0){
        y[dupliY]=y[dupliY]+runif(nb_dupliY,-0.00001,0.00001)}
    if (sum(ties)!=0){
    x[ties] <- x[ties]+runif(sum(ties),-0.00001,0.00001)}
    }
  }
  result=pval_comput(x,y,simu=simu,N=N)
  Pval<-result$Pval
  Tn<-result$Tn
  }
  #Pval<-1-ecdf_fun(Tn)
  result <- list(statistic=Tn, p.value=Pval,message=Message)
  class(result)<-"indeptest"
  return(result)
}
#' @export
print.indeptest <- function(x, ...)
{
  cat("\nRobust independence test for two continuous variables\n\n")
  if (round(x$p.value,4)==0){
    cat(paste("t = ", round(x$statistic,4), ", " , "p-value <1e-4","\n",sep= ""))
  } else {
  cat(paste("t = ", round(x$statistic,4), ", " , "p-value = ",round(x$p.value,4),"\n",sep= ""))}
  if (x$message==TRUE) {
    cat("\nTies were detected in the dataset and they were randomly broken")
  }
}

pval_comput<-function(x,y,simu=FALSE,N=50000){
  Tn<-stat_indeptest(x,y)
  n<-length(x)
  if (simu==TRUE){
    ecdf_fun<-simulecdf(n,N)
    Pval<-1-ecdf_fun(Tn)
  } else {
    #y1<-(1:(5e5))/(5e5)
    #data(ecdf10.Rdata, envir=environment())#paste(ecdf,n,.Rdata,sep="")
    #load(paste("ecdf",n,".Rdata",sep=""))#Tables/
    if (3<=n & n<=150)
    {
      x1<-robust_table_indep[[n]]#$x
      y1 <- (c(1:5e5)/5e5)[c(seq(20,1e5,by=2000),seq(1e5+200,5e5,by=200))]
      #y1<-robust_table_indep[[n]]$y
      ##y1<-robust_table[[n]]$y
      ###load(system.file(paste("data_tables/ecdf",n,".RData",sep=""),package="testRcpp"))
      ###data(list=paste("ecdf",n,sep=""))
    } else {
      if (151<=n & n<=175)
      {
        x1<-robust_table_indep[[150]]#$x
        y1 <- (c(1:5e5)/5e5)[c(seq(20,1e5,by=2000),seq(1e5+200,5e5,by=200))]
        # y1<-robust_table_indep[[150]]$y
        # #y1<-robust_table[[150]]$y
        # ##load(system.file("data_tables/ecdf150.RData",package="testRcpp"))
        # ##data(ecdf150)
      } else {
        if (176<=n & n<=250)
        {
          x1<-robust_table_indep[[200]]#robust_table_indep[[151]]$x
          y1 <- (c(1:5e5)/5e5)[c(seq(20,1e5,by=2000),seq(1e5+200,5e5,by=200))]
          # y1<-robust_table_indep[[151]]$y
          # #y1<-robust_table[[151]]$y
          # ##load(system.file("data_tables/ecdf200.RData",package="testRcpp"))
          # ##data(ecdf200)
        } else {
          if (251<=n & n<=400)
          {
            x1<-robust_table_indep[[300]]#robust_table_indep[[152]]$x
            y1 <- (c(1:5e5)/5e5)[c(seq(20,1e5,by=2000),seq(1e5+200,5e5,by=200))]
            # y1<-robust_table_indep[[152]]$y
            # #y1<-robust_table[[152]]$y
            # ##load(system.file("data_tables/ecdf300.RData",package="testRcpp"))
            # ##data(ecdf300)
          } else {
            if (401<=n & n<=750)
            {
              x1<-robust_table_indep[[500]]#robust_table_indep[[153]]$x
              y1 <- (c(1:5e5)/5e5)[c(seq(20,1e5,by=2000),seq(1e5+200,5e5,by=200))]
              # y1<-robust_table_indep[[153]]$y
              # #y1<-robust_table[[153]]$y
              # ##load(system.file("data_tables/ecdf500.RData",package="testRcpp"))
              # ##data(ecdf500)
            } else {
              if (751<=n)
              {
                x1<-robust_table_indep[[1000]]#robust_table_indep[[154]]$x
                y1 <- (c(1:5e5)/5e5)[c(seq(20,1e5,by=2000),seq(1e5+200,5e5,by=200))]
                # y1<-robust_table_indep[[154]]$y
                # #y1<-robust_table[[154]]$y
                # ##load(system.file("data_tables/ecdf1000.RData",package="testRcpp"))
                # ##data(ecdf1000)
              }
            }
          }
        }
      }
    }
    funstep<-stats::stepfun(x1,c(0,y1))
    Pval<-1-funstep(Tn)
  }
  return(list(Tn=Tn,Pval=Pval))
}



