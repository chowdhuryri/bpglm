#' Function to calculate mu1 and mu2
#' @param X1    : Covariates X1
#' @param X2    : Covariats X2 same as X1
#' @param icob1 : Beta estimates corresponding to y1
#' @param icob2 : Beta estimates corresponding to y2
#' @param mtype : model type
#' @return return estimated values of Lambdas
link_func<-function(X1,X2,icob1,icob2,mtype)
{
      lambda1 <- exp(X1%*%icob1)
      lambda2 <- exp(X2%*%icob2)

      Lambdas<-cbind(lambda1,lambda2)
      Lambdas
}
