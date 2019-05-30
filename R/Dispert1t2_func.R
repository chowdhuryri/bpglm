#' This is a function to calculate General test for over dispersion
#' An extension of Hosmer-Lemeshow test for Poisson-Poisson
#' @param y1    : Outcome variabe y1
#' @param y2    : Outcome variable y2
#' @param Lambdas : Lambda estimate from model
#' @param p : number of parameters
#' @param mtype : model type
#' @return Provides T1 statistic with 2*G d.f.
#' @importFrom stats pchisq
#' @importFrom stats aggregate
#' @importFrom gsl gamma_inc
#' @importFrom hypergeo genhypergeo

disper_fun<-function(y1,y2,Lambdas,p,mtype)  {

      y1<-as.matrix(y1,ncol=1)
      y2<-as.matrix(y2,ncol=1)

      lambda1<-as.matrix(Lambdas[,1],ncol=1)
      lambda2<-as.matrix(Lambdas[,2],ncol=1)

      l2y1<-lambda2*y1
      k1<-max(y1)
      k2<-max(y2)

      if (mtype==1){

            meany1 <-  lambda1
            meany2 <-  lambda2
            vary1  <- meany1
            vary2  <- meany2

            Vmuk1  <- lambda1
            Vmuk2  <- meany2

      } else if (mtype==2){

            meany1 <-  lambda1
            meany2 <-  lambda1*lambda2
            vary1  <- meany1
            vary2  <- meany2

            Vmuk1  <- lambda1
            Vmuk2  <- meany2

      } else if (mtype==3){

            meany1    <-(lambda1*exp(lambda1)/(exp(lambda1)-1))
            meany2    <-(lambda2*exp(lambda2)/(exp(lambda2)-1))

            vary1     <-(lambda1*exp(lambda1)/(exp(lambda1)-1))*(1-(lambda1/(exp(lambda1)-1)))
            vary2     <-(lambda2*exp(lambda2)/(exp(lambda2)-1))*(1-(lambda2/(exp(lambda2)-1)))

            #deltak1<- meany1 - lambda1
            #deltak2<- meany2 - lambda2
           # Vmuk1  <- lambda1-deltak1*(meany1-1)
            #Vmuk2  <- lambda2-deltak2*(meany2-1)


            Vmuk1  <-meany1
            Vmuk2  <-meany2

      } else if (mtype==4){

            meany1    <-(lambda1*exp(lambda1)/(exp(lambda1)-1))
            meany2    <-((lambda2*y1*exp(lambda2*y1))/(exp(lambda2*y1)-1))

            vary1     <-(lambda1*exp(lambda1)/(exp(lambda1)-1))*(1-(lambda1/(exp(lambda1)-1)))
            vary2     <-((lambda2*y1*exp(lambda2*y1))/(exp(lambda2*y1)-1))*(1-(lambda2*y1)/((exp(lambda2*y1)-1)))

            #deltak1<- meany1 - lambda1
            #deltak2<- meany2 - lambda2
            #Vmuk1  <- lambda1-deltak1*(meany1-1)
            #Vmuk2  <- lambda2-deltak2*(meany2-1)

            Vmuk1  <-meany1
            Vmuk2  <-meany2

       } else if (mtype==5){

            meany1 <- gamma_inc(k1,lambda1) * lambda1 * k1 / gamma_inc(k1+1,lambda1)
            meany2 <- gamma_inc(k2,lambda2) * lambda2 * k2 / gamma_inc(k2+1,lambda2)

            vary1  <- -gamma_inc(k1,lambda1) ^ 2 * lambda1 ^ 2 * k1 ^ 2/gamma_inc(k1 + 1,lambda1) ^ 2 + 1 / gamma_inc(k1 + 1,lambda1)* gamma(k1 + 1) * (lambda1 * (lambda1 + 1) - (k1 + 1) ^ 2
                      * lambda1 ^ (k1 + 1) / factorial(k1 + 1) / exp(lambda1)
                      * genhypergeo(c(1,k1 + 2),c(k1 + 1,k1 + 1), lambda1))

            vary2  <- -gamma_inc(k2,lambda2) ^ 2 * lambda2 ^ 2 * k2 ^ 2/gamma_inc(k2 + 1,lambda2) ^ 2 + 1 / gamma_inc(k2 + 1,lambda2)* gamma(k2 + 1) * (lambda2 * (lambda2 + 1) - (k2 + 1) ^ 2
                      * lambda2 ^ (k2 + 1) / factorial(k2 + 1) / exp(lambda2)
                      * genhypergeo(c(1,k2 + 2),c(k2 + 1,k2 + 1), lambda2))

            deltak1<- meany1 - lambda1
            deltak2<- meany2 - lambda2

            Vmuk1  <- lambda1-deltak1*(meany1-k1-1)
            Vmuk2  <- lambda2-deltak2*(meany2-k2-1)

      } else if (mtype==6){

            meany1 <- gamma_inc(k1,lambda1) * lambda1 * k1 / gamma_inc(k1+1,lambda1)
            meany2 <- lambda2 * y1 * gamma_inc(k2,l2y1) * k2 / gamma_inc(k2+1,l2y1)

            vary1  <- -gamma_inc(k1,lambda1) ^ 2 * lambda1 ^ 2 * k1 ^ 2 /gamma_inc(k1 + 1,lambda1) ^ 2 + 1 /gamma_inc(k1 + 1,lambda1)* gamma(k1 + 1) * (lambda1 * (lambda1 + 1) - (k1 + 1) ^ 2
                      * lambda1 ^ (k1 + 1) / factorial(k1 + 1) / exp(lambda1)
                      * genhypergeo(c(1,k1 + 2), c(k1 + 1,k1 + 1), lambda1))

            vary2  <- -lambda2 ^ 2 * y1 ^ 2 * gamma_inc(k2,l2y1) ^ 2 * k2 ^ 2/gamma_inc(k2+1,l2y1) ^ 2 + gamma(k2 + 1)/gamma_inc(k2+1,l2y1)* (lambda2 ^ 2 * y1 ^ 2 + lambda2 * y1 - (k2 + 1) ^ 2
                      / exp(lambda2 * y1) * (lambda2 * y1) ^ (k2 + 1) /
                      factorial(k2 + 1) * genhypergeo(c(1,k2 + 2), c(k2 + 1,k2 + 1), l2y1))

                        deltak1<- meany1 - lambda1
                        deltak2<- meany2 - lambda2
                        Vmuk1  <- lambda1-deltak1*(meany1-k1-1)
                        Vmuk2  <- lambda2-deltak2*(meany2-k2-1)

      } else if (mtype==-2){

            meany1 <-exp(lambda1) * gamma_inc(k1, lambda1) * k1 * lambda1 / (exp(lambda1)
                                                                             * gamma_inc(k1 + 1, lambda1) - gamma(k1 + 1))
            meany2 <-exp(lambda2) * gamma_inc(k2, lambda2) * k2 * lambda2 / (exp(lambda2) * gamma_inc(k2 + 1, lambda2)
                                                                             - gamma(k2 + 1))

            vary1  <--(exp(2 * lambda1) * gamma_inc(k1, lambda1) ^ 2 * k1 ^ 2 * lambda1 ^ 2 - gamma_inc(k1 + 1, lambda1)
                       * gamma(k1 + 1)* exp(2 * lambda1) * lambda1 ^ 2 - gamma_inc(k1 + 1, lambda1) * gamma(k1 + 1)
                       * exp(2 * lambda1) * lambda1 + gamma_inc(k1 + 1, lambda1) * genhypergeo(c(1,k1 + 2), c(k1 + 1,k1 + 1), lambda1)
                       * lambda1 ^ (k1 + 1) * exp(lambda1)* k1 + gamma(k1 + 1) ^ 2 * exp(lambda1) * lambda1 ^ 2
                       + gamma_inc(k1 + 1, lambda1)* genhypergeo(c(1,k1 + 2), c(k1 + 1,k1 + 1), lambda1) * lambda1 ^ (k1 + 1)
                       * exp(lambda1) + gamma(k1 + 1) ^ 2 * exp(lambda1) * lambda1 - gamma(k1 + 1)
                       * genhypergeo(c(1,k1 + 2), c(k1 + 1,k1 + 1), lambda1) * lambda1 ^ (k1 + 1) * k1 - gamma(k1 + 1)
                       * genhypergeo(c(1,k1 + 2), c(k1 + 1,k1 + 1), lambda1) * lambda1 ^ (k1 + 1)) / (exp(lambda1)
                                                                                                      * gamma_inc(k1 + 1, lambda1) - gamma(k1 + 1)) ^ 2
            vary2  <--(exp(2 * lambda2) * gamma_inc(k2, lambda2) ^ 2 * k2 ^ 2 * lambda2 ^ 2 - gamma_inc(k2 + 1, lambda2)
                       * gamma(k2 + 1) * exp(2 * lambda2) * lambda2 ^ 2 - gamma_inc(k2 + 1, lambda2) * gamma(k2 + 1)
                       * exp(2 * lambda2) * lambda2 + gamma_inc(k2 + 1, lambda2) * genhypergeo(c(1,k2 + 2), c(k2 + 1,k2 + 1), lambda2)
                       * lambda2 ^ (k2 + 1) * exp(lambda2)* k2 + gamma(k2 + 1) ^ 2 * exp(lambda2) * lambda2 ^ 2
                       + gamma_inc(k2 + 1, lambda2)* genhypergeo(c(1,k2 + 2), c(k2 + 1,k2 + 1), lambda2) * lambda2 ^ (k2 + 1)
                       * exp(lambda2) + gamma(k2 + 1) ^ 2 * exp(lambda2) * lambda2 - gamma(k2 + 1) * genhypergeo(c(1,k2 + 2),
                        c(k2 + 1,k2 + 1), lambda2) * lambda2 ^ (k2 + 1)* k2 - gamma(k2 + 1) * genhypergeo(c(1,k2 + 2),
                        c(k2 + 1,k2 + 1), lambda2) * lambda2 ^ (k2 + 1)) / (exp(lambda2)* gamma_inc(k2 + 1, lambda2)
                        - gamma(k2 + 1)) ^ 2

            deltak1<- meany1 - lambda1
            deltak2<- meany2 - lambda2
            Vmuk1  <- lambda1-deltak1*(meany1-k1-1)
            Vmuk2  <- lambda2-deltak2*(meany2-k1-1)

      } else if (mtype==-3){

            meany1    <-exp(lambda1) * gamma_inc(k1, lambda1) * k1 * lambda1 / (exp(lambda1) * gamma_inc(k1 + 1, lambda1)
                        - gamma(k1+ 1))
            meany2    <-lambda2 * y1 * exp(lambda2 * y1) * k2 * (k2 + 1) * gamma_inc(k2, lambda2 * y1) / (exp(lambda2 * y1) *
                        gamma_inc(k2 + 1, lambda2 * y1) * k2 + exp(lambda2 * y1) * gamma_inc(k2 + 1, lambda2 * y1) - gamma(k2 + 2))

            vary1     <--(exp(2 * lambda1) * gamma_inc(k1, lambda1) ^ 2 * k1 ^ 2 * lambda1 ^ 2 - gamma_inc(k1 + 1, lambda1)
                          * gamma(k1 + 1)* exp(2 * lambda1) * lambda1 ^ 2 - gamma_inc(k1 + 1, lambda1) * gamma(k1 + 1)
                          * exp(2 * lambda1) * lambda1 + gamma_inc(k1 + 1, lambda1) * genhypergeo(c(1,k1 + 2), c(k1 + 1,k1 + 1),
                            lambda1) * lambda1 ^ (k1 + 1) * exp(lambda1)* k1 + gamma(k1 + 1) ^ 2 * exp(lambda1) * lambda1 ^ 2
                          + gamma_inc(k1 + 1, lambda1)* genhypergeo(c(1,k1 + 2), c(k1 + 1,k1 + 1), lambda1) * lambda1 ^ (k1 + 1)
                          * exp(lambda1) + gamma(k1 + 1) ^ 2
                          * exp(lambda1) * lambda1 - gamma(k1 + 1) * genhypergeo(c(1,k1 + 2), c(k1 + 1,k1 + 1), lambda1) *
                                lambda1 ^ (k1 + 1)* k1 - gamma(k1 + 1) * genhypergeo(c(1,k1 + 2), c(k1 + 1,k1 + 1), lambda1) *
                                lambda1 ^ (k1 + 1)) / (exp(lambda1)* gamma_inc(k1 + 1, lambda1) - gamma(k1 + 1)) ^ 2

            vary2    <-lambda2 * y1 * (-exp(2 * lambda2 * y1) * k2 ^ 4 * lambda2 * y1 * gamma_inc(k2, lambda2 * y1) ^ 2 - 2
                                       * exp(2 * lambda2 * y1) * k2 ^ 3 * lambda2 * y1 * gamma_inc(k2, lambda2 * y1) ^ 2
                                       - exp(2 * lambda2 * y1)* k2 ^ 2 * lambda2 * y1 * gamma_inc(k2, lambda2 * y1) ^ 2
                                       - genhypergeo(c(1,k2 + 2), c(k2 + 1,k2 + 1), lambda2 * y1)* exp(lambda2 * y1)
                                       * (lambda2 * y1) ^ k2 * gamma_inc(k2 + 1, lambda2 * y1) * k2 ^ 3 + gamma(k2 + 2)
                                       * exp(2 * lambda2 * y1) * gamma_inc(k2 + 1, lambda2 * y1) * k2 * lambda2 * y1 - 3
                                       * genhypergeo(c(1,k2 + 2), c(k2 + 1,k2 + 1), lambda2 * y1) * exp(lambda2 * y1)
                                       * (lambda2 * y1) ^ k2 * gamma_inc(k2 + 1, lambda2 * y1) * k2 ^ 2 + gamma(k2 + 2)
                                       * exp(2 * lambda2 * y1) * gamma_inc(k2 + 1, lambda2 * y1) * lambda2 * y1
                                       + genhypergeo(c(1,k2 + 2), c(k2 + 1,k2 + 1), lambda2 * y1) * gamma(k2 + 2)
                                       * (lambda2 * y1) ^ k2 * k2 ^ 2 - 3 * genhypergeo(c(1,k2 + 2), c(k2 + 1,k2 + 1),
                                          lambda2 * y1) * exp(lambda2 * y1) * (lambda2 * y1) ^ k2 * gamma_inc(k2 + 1, lambda2 * y1)
                                       * k2 - gamma(k2 + 2) ^ 2 * exp(lambda2 * y1) * lambda2 * y1 + gamma(k2 + 2)
                                       * exp(2 * lambda2 * y1) * gamma_inc(k2 + 1, lambda2 * y1) * k2 + 2
                                       * genhypergeo(c(1,k2 + 2), c(k2 + 1,k2 + 1), lambda2 * y1) * gamma(k2 + 2)
                                       * (lambda2 * y1) ^ k2 * k2 - genhypergeo(c(1,k2 + 2), c(k2 + 1,k2 + 1), lambda2 * y1)
                                       * exp(lambda2 * y1) * (lambda2 * y1) ^ k2 * gamma_inc(k2 + 1, lambda2 * y1) + gamma(k2 + 2)
                                       * exp(2 * lambda2 * y1) * gamma_inc(k2 + 1, lambda2 * y1)
                                       + genhypergeo(c(1,k2 + 2), c(k2 + 1,k2 + 1), lambda2 * y1) * gamma(k2 + 2)
                                       * (lambda2 * y1) ^ k2 - gamma(k2 + 2) ^ 2 * exp(lambda2 * y1)) / (-exp(lambda2 * y1)
                                       * gamma_inc(k2 + 1, lambda2 * y1) * k2 - exp(lambda2 * y1)
                                       * gamma_inc(k2 + 1, lambda2 * y1) + gamma(k2 + 2)) ^ 2

            deltak1<- meany1 - lambda1
            deltak2<- meany2 - lambda2
            Vmuk1  <- lambda1-deltak1*(meany1-k1-1)
            Vmuk2  <- lambda2-deltak2*(meany2-k2-1)

          }

      phi1<-(sum((y1-meany1)^2/Vmuk1))/(nrow(y1)-p)
      phi2<-(sum((y2-meany2)^2/Vmuk2))/(nrow(y2)-p)

      #cat("Phi1 = ",phi1,"\n")
     # cat("Phi2 = ",phi2,"\n")

      DDD<-cbind(y1,y2)
      y1.zero <- which(DDD[,1]==0,arr.ind = T)
      DDD[y1.zero,1]<-1

      ttd<-cbind(DDD[,1:2],meany1,meany2,vary1,vary2)

      y1s1<-aggregate(ttd[,1],as.data.frame(ttd[,1]),sum)
      muhy1<-aggregate(ttd[,3],as.data.frame(ttd[,1]),mean)
      vahy1<-aggregate(ttd[,5],as.data.frame(ttd[,1]),mean)

      y2by1<-aggregate(ttd[,2],as.data.frame(ttd[,1]),mean)
      muhy2<-aggregate(ttd[,4],as.data.frame(ttd[,1]),mean)
      vahy2<-aggregate(ttd[,6],as.data.frame(ttd[,1]),mean)

      if(y1s1[1,1]==0){
            chi1<-(y1s1[1,1]-muhy1[1,2])* (1/vahy1[1,2]) *(y1s1[1,1]-muhy1[1,2])
            chi2<-(y1s1[1,1]-muhy1[1,2])* ((1/vahy1[1,2])*phi1) *(y1s1[1,1]-muhy1[1,2])
            y1s1<-y1s1[2:nrow(y1s1),]
            muhy1<-muhy1[2:nrow(muhy1),]
            vahy1<-vahy1[2:nrow(vahy1),]
            y2by1<-y2by1[2:nrow(y2by1),]
            muhy2<-muhy2[2:nrow(muhy2),]
            vahy2<-vahy2[2:nrow(vahy2),]
      }
      matrow<-length(y1s1[,1])

      chidT1<-0
      chidT2<-0

      for (i in 1:matrow){

            ST1<-0
            Ut1<-matrix(0,2,1)
            It1<-matrix(0,2,2)
            ST2<-0
            It2<-matrix(0,2,2)

            Ut1[1,1]<-y1s1[i,1]-muhy1[i,2]
            Ut1[2,1]<-y2by1[i,2]-muhy2[i,2]
            It1[1,1]<-vahy1[i,2]
            It1[1,2]<-0
            It1[2,1]<-0
            It1[2,2]<-vahy2[i,2]

            It2[1,1]<-vahy1[i,2] *phi1
            It2[1,2]<-0
            It2[2,1]<-0
            It2[2,2]<-vahy2[i,2]  *phi2

            ST1<-t(Ut1)%*%solve(It1)%*%Ut1
            ST2<-t(Ut1)%*%solve(It2)%*%Ut1
            chidT1<-chidT1+ST1
            chidT2<-chidT2+ST2
      }
      tdf<-2*matrow
      if(exists("chi1")){
            tdf<-tdf+1
            chidT1<-chidT1+chi1
            chidT2<-chidT2+chi2
      }

      pchidT1 <-pchisq(chidT1, tdf, ncp = 0, lower.tail = FALSE, log.p = FALSE)
      pchitabT1<-data.frame(ChiSquare=chidT1,DF=tdf,p_value=pchidT1)
      pchidT2 <-pchisq(chidT2, tdf, ncp = 0, lower.tail = FALSE, log.p = FALSE)
      pchitabT2<-data.frame(ChiSquare=chidT2,DF=tdf,p_value=pchidT2)
      rownames(pchitabT1)<-"T1 (GOF)" # This is for GOF No phi multiplied
      rownames(pchitabT2)<-"T2(Dispersion Test)"
      t1t2<-rbind(pchitabT1, pchitabT2)

      #a1<-mean(meany1)
      #a2<-mean(vary1)
      #a3<-mean(meany2)
      #a4<-mean(vary2)
      #cat("Muhat1=",a1,"\n")
      #cat("Variance1=",a2,"\n")
      #cat("Muhat2=",a3,"\n")
      #cat("Variance2=",a4,"\n")
      #print(t1t2)
      return(t1t2)

}

