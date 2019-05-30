#' Function to calculate First derivatives
#' @param y1    : Outcome variabe y1
#' @param y2    : Outcome variable y2
#' @param Lambdas : Estimated mu1 and mu2 from model
#' @param phi1 : Estimated phi1
#' @param phi2 : Estimated phi2
#' @param mtype : model type
#' @return returns estimated values llc - log likelihood
#' @importFrom gsl gamma_inc

loglik_func<-function(y1,y2,Lambdas,phi1,phi2,mtype){
      llc     <- 0
      y1logl  <- 0
      y2logl  <- 0
      lambda1 <- Lambdas[,1]
      lambda2 <- Lambdas[,2]

      if (mtype==1){

            y1logll <- y1*log(lambda1)-lambda1-log(factorial(y1))
            y2logll <- y2*log(lambda2)-lambda2-log(factorial(y2))
            y1logl <-sum(y1logll)
            y2logl <-sum(y2logll)

      } else if (mtype==2){

            pm <-ifelse(y1==0,log(1.00000001),log(y1))
            y1logll <- y1*log(lambda1)-lambda1-log(factorial(y1))
            y2logll <- y2*log(lambda2)-lambda2*y1+y2*pm-log(factorial(y2))
            y1logl <-sum(y1logll)
            y2logl <-sum(y2logll)

      } else if (mtype==3){

            pm      <-ifelse(y1==0,log(1.00000001),log(y1))
            y1logll <- y1 * log(lambda1) - log(factorial(y1)) - log(exp(lambda1) - 1)
            y2logll <- y2 * log(lambda2) - log(factorial(y2)) - log(exp(lambda2) - 1)
            y1logl  <-sum(y1logll)
            y2logl  <-sum(y2logll)

      } else if (mtype==4){

            pm      <-ifelse(y1==0,log(1.00000001),log(y1))
            y1logll <-  y1 * log(lambda1) - log(factorial(y1)) - log(exp(lambda1) - 1)
            y2logll <- y2 * log(lambda2) + y2 * pm - log(factorial(y2)) - log(exp(y1 * lambda2) - 1)
            y1logl  <-sum(y1logll)
            y2logl  <-sum(y2logll)

      } else if (mtype==5){

            k1<-max(y1)
            k2<-max(y2)

            c1 <- gamma(k1 + 1) / gamma_inc(k1 + 1, lambda1)
            c2 <- gamma(k2 + 1) / gamma_inc(k2 + 1, lambda2)

            pm <-ifelse(y1==0,log(1.00000001),log(y1))

            y1logll <- y1*log(lambda1)-lambda1-log(factorial(y1)) +log(c1)
            y2logll <- y2*log(lambda2)-lambda2-log(factorial(y2)) +log(c2)
            y1logl  <-sum(y1logll)
            y2logl  <-sum(y2logll)

      } else if (mtype==6){

            k1<-max(y1)
            k2<-max(y2)

            c1 <- gamma(k1 + 1) / gamma_inc(k1 + 1, lambda1)
            c2 <- gamma(k2 + 1) / gamma_inc(k2 + 1, lambda2 * y1)

            pm <-ifelse(y1==0,log(1.00000001),log(y1))

            y1logll <- y1*log(lambda1)-lambda1-log(factorial(y1)) +log(c1)
            y2logll <- y2*log(lambda2)-lambda2*y1+y2*pm-log(factorial(y2))+log(c2)
            y1logl  <-sum(y1logll)
            y2logl  <-sum(y2logll)

      }

  llc<-cbind(y1logl+y2logl,y1logl,y2logl)
  llc1<-list(llc,y1logll,y2logll)
  return(llc1)
}


