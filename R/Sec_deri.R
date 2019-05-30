#' Function to calculate Second derivatives
#' @param y1    : Outcome variabe y1
#' @param y2    : Outcome variable y2
#' @param X1    : Covariates X1
#' @param X2    : Covariats X2 same as X1
#' @param Lambdas : Estimated mu1 and mu2 from model
#' @param phi1 : Estimated phi1
#' @param phi2 : Estimated phi2
#' @param tobs : Total number of observations n
#' @param nvar1 : Number of covariates ncol(X1)
#' @param nvar2 : Number of covariates ncol(X2)
#' @param nvar : Number of covariates ncol(X1+X2)
#' @param mtype : model type
#' @return ith second derivative matrix
#'
#' @importFrom magic adiag
#' @importFrom gsl gamma_inc
#'
sderiv<-function(y1,y2,X1,X2,Lambdas,phi1,phi2,tobs,nvar1,nvar2,nvar,mtype){
  AA   <- matrix(0,nrow=nvar1,ncol=nvar1) # 0,nrow=nvar,ncol=nvar1
  BB   <- matrix(0,nrow=nvar2,ncol=nvar2) # 0,nrow=nvar,ncol=nvar2
#  it1  <- matrix(0,nrow=nvar1,ncol=nvar1) # 0,nrow=nvar,ncol=nvar
#  it2  <- matrix(0,nrow=nvar2,ncol=nvar2) # 0,nrow=nvar,ncol=nvar
  ith  <- matrix(0,nrow=nvar,ncol=nvar)
  lambda1  <- Lambdas[,1]
  lambda2  <- Lambdas[,2]

  if(mtype==1){

    y1aa <- -1 * lambda1
    y2aa <- -1 * lambda2

  } else if(mtype==2){

    y1aa <-  -1 * lambda1
    y2aa <-  as.numeric(-1 * lambda2 * y1)

  } else if (mtype==3){

    y1aa <- (-lambda1 * exp(lambda1) / (exp(lambda1) - 1) - lambda1 ^ 2
             * exp(lambda1) / (exp(lambda1) - 1) + lambda1 ^ 2 * exp(lambda1)^2
             / (exp(lambda1) - 1) ^ 2)
    y2aa <- (-lambda2 * exp(lambda2) / (exp(lambda2) - 1) - lambda2 ^ 2
             * exp(lambda2) / (exp(lambda2) - 1) + lambda2 ^ 2 * exp(lambda2)^2
             / (exp(lambda2) - 1) ^ 2)

  } else if (mtype==4){

      y1aa <- (-lambda1 * exp(lambda1) / (exp(lambda1) - 1) - lambda1 ^ 2
               * exp(lambda1) / (exp(lambda1) - 1) + lambda1 ^ 2 * exp(lambda1)
               ^ 2 / (exp(lambda1) - 1) ^ 2)

      y2aa <- as.numeric(-y1 * lambda2 * exp(y1 * lambda2) / (exp(y1 * lambda2) - 1)
              - y1 ^ 2 * lambda2 ^ 2 * exp(y1 * lambda2) / (exp(y1 * lambda2)-1)
              + y1 ^ 2 * lambda2 ^ 2 * exp(y1 * lambda2) ^ 2 / (exp(y1 * lambda2) - 1) ^ 2)

  } else if (mtype==5){

        k1<-max(y1)
        k2<-max(y2)
        y1aa <- (-lambda1 - exp(-lambda1) * (-lambda1 ^ (k1 + 1)
                * gamma_inc(k1 + 1, lambda1) * k1 + lambda1 ^ (k1 + 2)
                * gamma_inc(k1 + 1, lambda1) - lambda1 ^ (k1 + 1)
                * gamma_inc(k1 + 1, lambda1) - lambda1 ^ (2 + 2 * k1)
                * exp(-lambda1)) / gamma_inc(k1 + 1, lambda1) ^ 2)

        y2aa <- (-lambda2 - exp(-lambda2) * (-lambda2 ^ (k2 + 1)
                * gamma_inc(k2 + 1, lambda2) * k2 + lambda2 ^ (k2 + 2)
                * gamma_inc(k2 + 1, lambda2) - lambda2 ^ (k2 + 1)
                * gamma_inc(k2 + 1, lambda2) - lambda2 ^ (2 + 2 * k2)
                * exp(-lambda2)) / gamma_inc(k2 + 1, lambda2) ^ 2)

  } else if (mtype==6){

        k1<-max(y1)
        k2<-max(y2)

        y1aa <- (-lambda1 - exp(-lambda1) * (-lambda1 ^ (k1 + 1)
                * gamma_inc(k1 + 1, lambda1) * k1 + lambda1 ^ (k1 + 2)
                * gamma_inc(k1 + 1, lambda1) - lambda1 ^ (k1 + 1)
                * gamma_inc(k1 + 1, lambda1) - lambda1 ^ (2 + 2 * k1)
                * exp(-lambda1)) / gamma_inc(k1 + 1, lambda1) ^ 2)

        y2aa <- as.numeric(-lambda2 * y1 + exp(-lambda2 * y1) * lambda2 * y1 *
                (-(lambda2 * y1) ^ k2 * lambda2 * y1 * gamma_inc(k2 + 1, lambda2 * y1)
                 + exp(-lambda2 * y1) * (lambda2 * y1) ^ (2 * k2) * lambda2 * y1
                 + (lambda2 * y1) ^ k2 * gamma_inc(k2 + 1, lambda2 * y1) * k2
                 + (lambda2 * y1) ^ k2 * gamma_inc(k2 + 1, lambda2 * y1))
                / gamma_inc(k2 + 1, lambda2 * y1) ^ 2)

 }

    #  if(nvar>2){
      if(nvar1>1){
            AA<-matrix(apply((t(apply(X1, 1, function(x) (x %*% t(x))))*y1aa),
                             2,sum),nrow=nvar1,ncol=nvar1)
      }
      if(nvar2>1){
            BB<-matrix(apply((t(apply(X2, 1, function(x) (x %*% t(x))))*y2aa),
                             2,sum),nrow=nvar2,ncol=nvar2)
      }

  #    } else if(nvar==2){
      if(nvar1==1){
            AA<-matrix(apply((t(apply(X1, 1, function(x) (x %*% t(x))))*y1aa),
                             1,sum),nrow=nvar1,ncol=nvar1)
      }
      if(nvar2==1){
            BB<-matrix(apply((t(apply(X2, 1, function(x) (x %*% t(x))))*y2aa),
                             1,sum),nrow=nvar2,ncol=nvar2)
      }
  #    }

      ith  <- adiag(AA,BB)
      ith
}
