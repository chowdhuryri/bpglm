
is.int <- function(x, tol = .Machine$double.eps^0.5)  {
      abs(x - round(x)) < tol}

tollf<-function(icob,icob1){
      icob<-matrix(icob,ncol=1)
      icob1<-matrix(icob1,ncol=1)
      icd<-icob-icob1
      ttcon<-sqrt((t(icd)%*%icd)/(t(icob)%*%icob))
      return(ttcon)
}

#' This is a function to perform Chi-square test.
#' Takes the input of Log-Likelihood from Null Model and Alternative Model.
#' @param NullM  Reduced model fit results
#' @param AlterM  Full model fit results
#' @return Provides Chisquare test and p-values
#'
#' @export
#' @importFrom stats pchisq

ChiNA<-function(NullM,AlterM){
      MyChi <--2*(NullM$logLik[1])+2*(AlterM$logLik[1])
      tdf   <-(AlterM$nvar)
      pchiV <-pchisq(MyChi, tdf, ncp = 0, lower.tail = FALSE, log.p = FALSE)
      Chitab<-data.frame(ChiSquare=MyChi,DF=tdf,p_value=pchiV)
      cat("Likelihood Ratio Test: Null Model vs. Alternative Model","\n")
      cat("----------------------------------------------------------", "\n")
      print(Chitab,digits=8)
      cat("----------------------------------------------------------", "\n")

}

#' This is a function to perform Chi-square test.
#' Takes the input of Log-Likelihood from Full Model and Reduced Model.
#' @param Redm  Reduced model fit results
#' @param Fullm  Full model fit results
#' @return Provides Chisquare test and p-values
#'
#' @export
#' @importFrom stats pchisq

ChiRF<-function(Redm,Fullm){
      MyChi <--2*(Redm$logLik[1])+2*(Fullm$logLik[1])
      tdf   <- Fullm$nvar - Redm$nvar
      pchiV <-pchisq(MyChi, tdf, ncp = 0, lower.tail = FALSE, log.p = FALSE)
      Chitab<-data.frame(ChiSquare=MyChi,DF=tdf,p_value=pchiV)
      cat("Likelihood Ratio Test: Reduced Model vs. Full Model","\n")
      cat("----------------------------------------------------------", "\n")
      print(Chitab,digits=8,row.names=FALSE)
      cat("----------------------------------------------------------", "\n")

}

#' This is a function to perform Chi-square test for goodness of fit. Based on
#' the predicted probability.
#' @param y1    : Outcome variabe y1
#' @param y2    : Outcome variable y2
#' @param Lambdas : Estimated mu1 and mu2.
#' @param nvar : total number of parameter in bivariate model (including cons.)
#' @param mtype : model type
#' @return Provides Chisquare test and p-values
#' @importFrom stats pchisq

PearChiPred<-function(y1,y2,Lambdas,nvar,mtype){
      mydt<-cbind(data.frame(cbind(y1,y2)),Lambdas)
      cname<-c(names(mydt)[1:2],"lambda1","lambda2")

      colnames(mydt)  <- cname
      lam1_tot        <- sum(mydt$lambda1)
      lam2_tot        <- sum(mydt$lambda2)
      y1_split        <- split(mydt,list(mydt[,1]))
      y2_split        <- split(mydt,list(mydt[,2]))
      y1_y2_split     <- split(mydt[,c(1:2,4)],list(mydt[,1],mydt[,2]))
      y1_mar_prob     <- matrix(sapply(y1_split, function(x) sum(x$lambda1))/lam1_tot,ncol=1)
      y2_mar_prob     <- matrix(sapply(y2_split, function(x) sum(x$lambda2))/lam2_tot,nrow=1)
      y1_con_sum      <- matrix(sapply(y1_y2_split, function(x) sum(x$lambda2)),nrow=length(y1_mar_prob))
      y1_con_mar_sum  <- matrix(apply(y1_con_sum,1,sum),nrow=dim(y1_con_sum)[1],ncol=dim(y1_con_sum)[2])
      y1_con_prob     <- y1_con_sum/y1_con_mar_sum
      y1_mar_prob_n   <- matrix(y1_mar_prob,nrow=dim(y1_mar_prob)[1],ncol=dim(y1_con_mar_sum)[2])

      if (mtype==2 | mtype==4 | mtype==6) {
            y1y2_joint_prob <-  y1_con_prob*y1_mar_prob_n
      } else if (mtype==1 | mtype==3 | mtype==5) {
            y1y2_joint_prob <-  y1_mar_prob %*% y2_mar_prob
      }

      y1y2_pred_count <- y1y2_joint_prob*nrow(mydt)
      y1y2_obs_count  <- as.matrix(table(mydt[,1],mydt[,2]))
      y1y2_obs_prob   <- y1y2_obs_count/nrow(mydt)
      rownames(y1_con_prob)     <- rownames(y1y2_obs_count)
      colnames(y1_con_prob)     <- colnames(y1y2_obs_count)
      rownames(y1y2_joint_prob) <- rownames(y1y2_obs_count)
      colnames(y1y2_joint_prob) <- colnames(y1y2_obs_count)
      rownames(y1y2_pred_count) <- rownames(y1y2_obs_count)
      colnames(y1y2_pred_count) <- colnames(y1y2_obs_count)

      chisqA1<- (y1y2_obs_count-y1y2_pred_count)^2/y1y2_pred_count
      chisqA1<- sum(chisqA1,na.rm=TRUE)
      #chigfdf<- (nrow(y1y2_pred_count)-1)*(ncol(y1y2_pred_count)-1)
      chigfdf<- nrow(y1y2_pred_count)*ncol(y1y2_pred_count)-nvar
      chisqA1_prob<-pchisq(chisqA1, chigfdf, ncp = 0, lower.tail = FALSE, log.p = FALSE)

      chiq_res<-data.frame("Chi-square" = chisqA1, "D.F" = chigfdf, "p-value" = chisqA1_prob)
      list_res<-list("ConProbPred" = y1_con_prob, "JointProbPred"=y1y2_joint_prob,
                     "JointProbObs"=y1y2_obs_prob, "ObsCount"= y1y2_obs_count,
                     "CountPred"=y1y2_pred_count, "ProbY1Pred"=y1_mar_prob,
                     "ProbY2Pred"=y2_mar_prob, "GOF.Chisquare" =chiq_res)
      return(list_res)
}

#' This is a function to calculate the z-test for y1 and y2 for dispersion (Impirical). Based on observed y1 and y2.
#' @param y1    : Outcome variabe y1
#' @param y2    : Outcome variable y2
#' @param tobs  : Total number of observations n
#' @return Provides  z-value and p-value for test for dispersion
#' @importFrom stats pnorm
#' @importFrom stats var

dpztest<-function(y1,y2,tobs){
      z1<-(mean(y1)-var(y1))/(sqrt(mean(y1))/sqrt(tobs))
      z2<-(mean(y2)-var(y2))/(sqrt(mean(y2))/sqrt(tobs))
      zz<-rbind(z1,z2)
      zpro <- 2 * pnorm(-abs(zz))
      zdes1<-data.frame(Z_value=zz,p_value=zpro)
      rownames(zdes1)<-c("Z(Y1) Marginal 2-tail","Z(Y2) Marginal 2-tail")
      return(zdes1)
} # end function


#' This is a function to perform Chi-square test.
#' Input saved results from Model 1 and Model 2.
#' @param Mod1  Model 1 fit results
#' @param Mod2  Model 2 fit results
#' @return ProvidesVoung Z-value and p-value
#'
#' @export
#' @importFrom stats pnorm
vountest<-function(Mod1, Mod2){
      outmat<-data.frame(matrix(0,nrow=Mod1$N,ncol=2))

      for (mo in 1: 2){
            if (mo==1){
                  y1logl   <- Mod1$logliky1
                  y2logl   <- Mod1$logliky2
                  phi1    <- Mod1$Phi1
                  phi2    <- Mod1$Phi2
                  y1      <- Mod1$y1
                  y2      <- Mod1$y2
            } else if (mo==2){
                  y1logl  <- Mod2$logliky1
                  y2logl  <- Mod2$logliky2
                  phi1    <- Mod2$Phi1
                  phi2    <- Mod2$Phi2
                  y1      <- Mod2$y1
                  y2      <- Mod2$y2
            }

            outmat[,mo]<-y1logl+y2logl
      } #mo

      outmat$mi <- outmat[,1]-outmat[,2]
      outmat$mid<- outmat$mi-(Mod1$nvar-Mod2$nvar)*log(Mod1$N)/(2*Mod1$N)
      sm        <- (sum(outmat$mi^2)-(sum(outmat$mi))^2/Mod1$N)/Mod1$N
      V1        <- mean(outmat$mi)*sqrt(Mod1$N)/sqrt(sm)
      V2        <- mean(outmat$mid)*sqrt(Mod1$N)/sqrt(sm)
      VV<-rbind(V1,V2)
      np1<-rbind(Mod1$nvar,Mod1$nvar)
      np2<-rbind(Mod2$nvar,Mod2$nvar)
      #vpro <- 2 * pnorm(-abs(VV))
      pcgt<-sum(as.numeric(outmat$mi>0))
      pclt<-sum(as.numeric(outmat$mi<0))
      cat("Positive mi>0 = ",pcgt,"\n")
      cat("Negative mi<0 = ",pclt,"\n")
      vpro <- pnorm(abs(VV),lower.tail=FALSE)
      vdes1<-data.frame(Z_value=VV,p=np1,q=np2,p_value=vpro)
      rownames(vdes1)<-c("V1 Un adj. one-tail","V1 Adj. one-tail")

      return(vdes1)
}

