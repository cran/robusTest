#' Break the ties in a given vector or between two vectors
#'
#' If the vector contains ties (either inside a single or between two vectors), the function breaks them using a random perturbation.
#' @param x,y the variables containing ties.
#' @param nb_break if TRUE return also the number of values that have been broken
#' @details If \code{y=NULL} the function detects the ties in the vector \code{x}. A uniform variable on the interval \eqn{[-e^(-5),e^(-5)]} is added
#' to the value of all the ties but one in the vector \code{x}. If \code{y} is also provided, the function detects the ties between
#' \code{x} and \code{y} and break them (only in the \code{x} vector) by adding a uniform variable on the interval \eqn{[-e^(-5),e^(-5)]} to these values.
#' If \code{nb_break=TRUE} the result is returned as a list that also includes the number of values that have been broken.
#' @return After breaking the ties, either returns the \code{x} variable as a vector or returns a list containing the \code{x} and \code{y}
#' variables. If \code{nb_break=TRUE}, the number of values that have been broken is also an element of the list.
#' @export
#' @examples
#' x <- c(1,2,2,3,4,5,5,5,7)
#' xbreak=tiebreak(x)
#' xbreak
#' #a uniform value has been added to the second, sixth and seventh value of x.
#' tiebreak(x,nb_break=TRUE) #3 values have been broken in x
#' sum(duplicated(xbreak))#check if the breaking procedure has worked.
#' y <- c(4,9,12,11,2,10)
#' xy_break=tiebreak(x,y)
#' xy_break$x
#' xy_break$y #a uniform value has been added to the second, third and fifth value of x.
#' xy_break$x%in%xy_break$y #check that no values for xbreak can be found in ybreak
#' tiebreak(x,y,nb_break=TRUE) #also returns the number of broken values

#' @export
tiebreak=function(x,y=NULL,nb_break=FALSE){
  if (is.null(y)){
    dupli=duplicated(x)
    nb_dupli=sum(dupli)
    if (nb_dupli==0) {
      warning("The data does not contain ties")
    } else {
      x[dupli]=x[dupli]+runif(nb_dupli,-0.00001,0.00001)}
    #Xsort=sort(X,index.return=TRUE)
    # index=which(diff(Xsort$x)==0)#which value should be changed in the ordered sample
    #X[Xsort$ix[index]]<-X[Xsort$ix[index]]+runif(length(Xsort$ix[index]),-0.00001,0.00001)
    if (nb_break==TRUE){
      return(list(x=x,nb_break=nb_dupli))
    } else {
    return(x)}
  } else {
    if (is.null(x)){ stop("a value for 'x' need to be provided")}
    ties=x%in%y
    if (sum(ties)==0){warning("The data does not contain ties")} else {
      x[ties] <- x[ties]+runif(sum(ties),-0.00001,0.00001)
    }
    if (nb_break==TRUE){
      return(list(x=x,y=y,nb_break=sum(ties)))
    } else {
    return(list(x=x,y=y))}
  }
}

