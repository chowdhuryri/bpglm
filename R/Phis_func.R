#' Function to calculate mu1 and mu2
#' @param y1    : Outcome variabe y1
#' @param y2    : Outcome variable y2
#' @param X1    : Covariates X1
#' @param X2    : Covariats X2 same as X1
#' @param icob1 : Beta estimates corresponding to y1
#' @param icob2 : Beta estimates corresponding to y2
#' @param tobs : Total number of observations n
#' @param nvar1 : Number of covariates ncol(X1)
#' @param nvar2 : Number of covariates ncol(X2)
#' @param nvar : Number of covariates ncol(X1+X2)
#' @param mtype : model type
#' @return return estimated values of phi1 and phi2
#'
phis_func<-function(y1,y2,X1,X2,icob1,icob2,tobs,nvar1,nvar2,nvar,mtype)
{
      xb1<-X1%*%icob1
      xb2<-X2%*%icob2

      lambda1 <- exp(xb1)

      if (mtype==1){

            lambda2 <- exp(xb2)

      } else if (mtype==2){

            lambda2<-as.vector( exp(xb1) * exp(xb2))

      } else if (mtype==3){

            lambda2 <- exp(xb2)

      } else if (mtype==4){

            lambda2 <- exp(xb2)

      } else if (mtype==5){

            lambda2 <- exp(xb2)

      } else if (mtype==6){

            lambda2<-as.vector( exp(xb1) * exp(xb2))
      }

      phi1<-(1/(tobs-nvar1))*sum((y1-lambda1)^2/lambda1)
      phi2<-(1/(tobs-nvar1))*sum((y2-lambda2)^2/lambda2)
      Phis<-cbind(phi1,phi2)
      Phis
}
