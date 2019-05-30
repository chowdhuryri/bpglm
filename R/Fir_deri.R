#' Function to calculate First derivatives
#' @param y1    : Outcome variabe y1
#' @param y2    : Outcome variable y2
#' @param X1    : Covariates X1
#' @param X2    : Covariats X2 same as X1
#' @param Lambdas : Estimated mu1 and mu2 from model
#' @param phi1 : Estimated phi1
#' @param phi2 : Estimated phi2
#' @param tobs : Total number of observations n
#' @param nvar : Number of covariates ncol(X1)
#' @param mtype : model type
#' @importFrom gsl gamma_inc
#'
fderiv<-function(y1,y2,X1,X2,Lambdas, phi1,phi2,tobs,nvar,mtype){

      lambda1    <- Lambdas[,1]
      lambda2    <- Lambdas[,2]

      if (mtype==1){

            ut1   <- as.numeric(y1-lambda1) * X1
            uth1  <- matrix(apply(ut1,2,sum),nrow=1)

            ut2   <- as.numeric(y2-lambda2) * X2
            uth2  <- matrix(apply(ut2,2,sum),nrow=1)

      } else if (mtype==2){

            ut1   <- as.numeric(y1-lambda1) * X1
            uth1  <- matrix(apply(ut1,2,sum),nrow=1)

            ut2   <- as.numeric(y2-y1*lambda2) * X2
            uth2  <- matrix(apply(ut2,2,sum),nrow=1)

      } else if (mtype==3){

            ut1  <-as.numeric(y1 - lambda1*exp(lambda1)/(exp(lambda1) - 1)) * X1
            uth1 <- matrix(apply(ut1,2,sum),nrow=1)

            ut2  <-as.numeric(y2 - lambda2*exp(lambda2)/(exp(lambda2) - 1)) * X2
            uth2 <- matrix(apply(ut2,2,sum),nrow=1)

      } else if (mtype==4){

            ut1  <-as.numeric(y1 - lambda1 *exp(lambda1)/(exp(lambda1) - 1)) * X1
            uth1 <- matrix(apply(ut1,2,sum),nrow=1)

            ut2  <-as.numeric(y2-y1*lambda2*exp(y1*lambda2)/(exp(y1*lambda2)-1)) * X2
            uth2 <- matrix(apply(ut2,2,sum),nrow=1)

      } else if (mtype==5){

            k1<-max(y1)
            k2<-max(y2)

            ut1<- as.numeric(y1 - lambda1 + 1 / gamma_inc(k1 + 1, lambda1) * lambda1 * lambda1 ^ k1 * exp(-lambda1)) * X1
            uth1 <- matrix(apply(ut1,2,sum),nrow=1)

            ut2<- as.numeric(y2 - lambda2 + 1 / gamma_inc(k2 + 1, lambda2) * lambda2 * lambda2 ^ k2 * exp(-lambda2)) * X2
            uth2 <- matrix(apply(ut2,2,sum),nrow=1)

      } else if (mtype==6){

            k1<-max(y1)
            k2<-max(y2)

            ut1<- as.numeric(y1 - lambda1 + 1 / gamma_inc(k1 + 1, lambda1) * lambda1 * lambda1 ^ k1 * exp(-lambda1)) * X1
            uth1 <- matrix(apply(ut1,2,sum),nrow=1)

            ut2 <- -X2 * as.numeric((gamma_inc(k2 + 1,lambda2 * y1) * lambda2 * y1 - exp(-lambda2 * y1 + lambda2) * y1 * (lambda2 * y1) ^ k2
                   - gamma_inc(k2 + 1,lambda2 * y1) * y2) / gamma_inc(k2 + 1,lambda2 * y1))
            uth2 <- matrix(apply(ut2,2,sum),nrow=1)

      }

      uth<-cbind(uth1,uth2)
      uth
}
