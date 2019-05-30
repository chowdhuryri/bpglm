#' This is a function calculate deviance statistics.
#' @param y1    : Outcome variabe y1
#' @param y2    : Outcome variable y2
#' @param Lambdas : Estimated mu1 and mu2 from model
#' @param phi1 : Estimated phi1
#' @param phi2 : Estimated phi2
#' @param mtype :model type
#' @return Provide deviance
#' @importFrom gsl gamma_inc
#'
deva_func<-function(y1,y2,Lambdas,phi1,phi2,mtype){
      lambda1 <- Lambdas[,1]
      lambda2 <- Lambdas[,2]

      if (mtype==1){

            A<-ifelse(y1==0,0,log(y1/lambda1))
            A1<-y1*A
            B<-y1-lambda1
            C<-ifelse(y2==0,0,log(y2/lambda2))
            C1<-y2*C
            D<-y2-lambda2
          #  devia<-2*sum((A1-B+C1-D))
            devia1<-2*sum(A1-B)
            devia2<-2*sum(C1-D)

      } else if (mtype==2){

            A<-ifelse(y1==0,0,log(y1/lambda1))
            A1<-y1*A
            B<-y2*ifelse(y2==0,0,log(y2/(lambda1*lambda2)))
            D<-y1-lambda1
            E<-(y2-((lambda1*lambda2)/lambda1)*y1)
           # devia<-2*sum(A1+B-y2*A-D-E)
            devia1<-2*sum(A1-D)
            devia2<-2*sum(B-y2*A-E)

      } else if (mtype==3){
            mu1    <-(lambda1*exp(lambda1)/(exp(lambda1)-1))
            mu2    <-(lambda2*exp(lambda2)/(exp(lambda2)-1))

            devia1 <-2*sum(y1*log(y1/mu1) - log((exp(y1)-1)/(exp(mu1)-1)))
            devia2 <-2*sum(y2*log(y2/mu2) - log((exp(y2)-1)/(exp(mu2)-1)))

      } else if (mtype==4){

            mu1    <-(lambda1*exp(lambda1)/(exp(lambda1)-1))
            mu2    <-(lambda2*y1*exp(lambda2*y1))/(exp(lambda2*y1)-1)

            devia1<-2*sum(y1*log(y1/mu1) - log((exp(y1)-1)/(exp(mu1)-1)))
            devia2<-2*sum(y2*log(y2/mu2) - log((exp(y2)-1)/(exp(mu2/mu1*y1)-1)))

      } else if (mtype==5){

            k1<-max(y1)
            k2<-max(y2)
            mu1    <-lambda1 * gamma(k1 + 1) * k1 / gamma_inc(k1 + 1,lambda1)
            mu2    <-lambda2 * gamma(k2 + 1) * k2 / gamma_inc(k2 + 1,lambda2)

            c1 <- gamma(k1 + 1) / gamma_inc(k1 + 1, y1)
            c2 <- gamma(k2 + 1) / gamma_inc(k2 + 1, y2)
            c1d <- gamma(k1 + 1) / gamma_inc(k1 + 1, mu1)
            c2d <- gamma(k2 + 1) / gamma_inc(k2 + 1, mu2)

            A<-ifelse(y1==0,0,log(y1/mu1))
            A1<-y1*A
            B<-y1-mu1
            C<-ifelse(y2==0,0,log(y2/lambda2)) #mu2 replaced bylambda2
            C1<-y2*C
            D<-y2-lambda2 #mu2 replaced bylambda2

            devia1<-2*sum(A1-B + log(c1) -log(c2))
            devia2<-2*sum(C1-D + log(c1d) -log(c2d))

      } else if (mtype==6){
            icc<-1
       if(icc==0){
             devia1<-0
             devia2<-0
       } else{
            k1<-max(y1)
            k2<-max(y2)
            mu1    <-lambda1 * gamma(k1 + 1) * k1 / gamma_inc(k1 + 1,lambda1)
            mu2    <-as.numeric(lambda2 * y1 *gamma(k2 + 1) * k2 / gamma_inc(k2 + 1,lambda2*y1))

            c1 <- gamma(k1 + 1) / gamma_inc(k1 + 1, y1)
            c2 <- gamma(k2 + 1) / gamma_inc(k2 + 1, y2 * y1)
            c1d <- gamma(k1 + 1) / gamma_inc(k1 + 1, mu1)
            c2d <- gamma(k2 + 1) / gamma_inc(k2 + 1, mu2 * y1)

            A<-ifelse(y1==0,0,log(y1/mu1))
            A1<-y1*A
            B<-y2*ifelse(y2==0,0,log(y2/mu2))
            D<-y1-mu1
            E<-(y2-(mu2/mu1)*y1)

            devia1<-2*sum(A1-D + log(c1) -log(c2))
            devia2<-2*sum(B-y2*A-E + log(c1d) -log(c2d))
            }
      }

      devia<-cbind(devia1+devia2,devia1,devia2)

      return(devia)
}
