#' This is main function to estimate the regression coefficients from BIVPOI
#' @param Y1Y2 :  Data frame (y1,y2)
#' @param X1 : Data frame for X1
#' @param X2 : Data frame for X2
#' @param mxit :  Maximum number of iterations
#' @param icob :  Vector of initial parameter estimates
#' @param mtype : Model type: 1= BIVPOI marginal/marginal; 2= BIVPOI marginal/conditional
#' 3 = ZTBIVPOI marginal/marginal; 4= ZTBIVPOI marginal/conditional
#' 5 = RTBIVPOI marginal/marginal; 6= RTBIVPOI marginal/conditional
#' @param ppy : Default value 1 to print results, 0 to suppress print
#' @param presc :  Precision of estimates  default value is 100000
#' @param phi :  Default "TRUE" Phi's will be estimated; 1 assumes Phis as 1
#' @return Returns a list of estimates and model statistics
#' @export

bpglm<-function(Y1Y2,X1=NULL,X2=NULL, mxit=150, icob=NULL, mtype=1,ppy=1,
                presc=100000,phi=NULL){

      if(is.null(X1)==TRUE & is.null(X2)==FALSE){
            X2<-as.matrix(X1,nrow=dim(X1)[1],ncol=dim(X1)[2])
      }

      if(is.null(X1)==TRUE){
            X1<-matrix(1,nrow=dim(Y1Y2)[1],ncol=1)
            colnames(X1)<-"Y1:Constant"
      } else {
            cons1<-matrix(1,nrow=dim(Y1Y2)[1],ncol=1)
            colnames(cons1)<-"Y1:Constant"
            X1   <-as.matrix(cbind(cons1,X1),nrow=dim(Y1Y2)[1])
      }
      nvar1 <- ncol(X1)
      rrna1 <-colnames(X1)

      if(is.null(X2)==TRUE){
            X2<-matrix(1,nrow=dim(Y1Y2)[1],ncol=1)
            colnames(X2)<-"Y2:Constant"
      } else {
            cons2<-matrix(1,nrow=dim(Y1Y2)[1],ncol=1)
            colnames(cons2)<-"Y2:Constant"
            X2   <-as.matrix(cbind(cons2,X2),nrow=dim(Y1Y2)[1])
      }

      nvar2 <- ncol(X2)
      rrna2 <-colnames(X2)
      nvar  <-nvar1+nvar2
      rrn<-c(rrna1,rrna2)
      dmat<-cbind(Y1Y2,X1,X2)
#print(head(dmat))
      dmat  <- poichkfun(dmat, mxit = 150, mtype)

      conv  <-0  # Not Converged
      tobs  <-nrow(dmat)
      y1n   <-names(dmat)[1]
      y2n   <-names(dmat)[2]
      y1    <-as.matrix(dmat[,1])
      y2    <-as.matrix(dmat[,2])
      X1    <-as.matrix(dmat[,3:(nvar1+2)])
      X2    <-as.matrix(dmat[,(nvar1+3):(nvar1+nvar2+2)])
#print(head(X1))
#print(head(X2))

      if (dim(dmat)[2] >4) {
            if(ppy==1){
                  cat("----- Bi-variate Poisson Regression -----","\n")
                  cat("----- Model with covariates --------------------","\n")
            }
      } else if(dim(dmat)[2] ==4) {
            if(ppy==1){
                  cat("----- Bi-variate Poisson Regression ------","\n")
                  cat("----- Constant only model.....No covariates -----","\n")
            }
      }

      if(is.null(icob)==TRUE){
            icob<-matrix(.01,nrow=nvar,ncol=1)
      }

      # loop for itaration

      for (ita in 1:mxit) {

            nicob1<-matrix(icob[1:nvar1],ncol=1)
            nicob2<-matrix(icob[(nvar1+1):(nvar)],ncol=1)
            Lambdas<-link_func(X1,X2,nicob1,nicob2,mtype)

            ## Phis
            if(is.null(phi)==TRUE){
                  Phis<-phis_func(y1,y2,X1,X2,nicob1,nicob2,tobs,nvar1,nvar2,nvar,mtype)
                  phi1<-Phis[1]
                  phi2<-Phis[2]
            } else{
                  phi1<-1
                  phi2<-1
            }
            ## Log Likelihood
            LLC<-loglik_func(y1,y2,Lambdas,phi1,phi2,mtype)
            y1logll<-LLC[[2]]
            y2logll<-LLC[[3]]
            LLC    <-LLC[[1]]
            cat("Iteration = ", ita, "\n")
            cat("Log Likelihood = ", LLC[1], "\n")

            # Deviance
            deva <-deva_func(y1,y2,Lambdas,phi1,phi2,mtype)

            ## Score vector
            uth  <- fderiv(y1,y2,X1,X2,Lambdas,phi1,phi2,tobs,nvar,mtype)

            ## Information Matrix
            ith  <- sderiv(y1,y2,X1,X2,Lambdas,phi1,phi2,tobs,nvar1,nvar2,nvar,mtype)
            ithus<-ith*(-1)
            ithinv<-solve(ithus)

            ## Beta estimates

            icob1<-icob - solve(ith) %*% t(uth)

            icob1se<-matrix(sqrt(diag(ithinv)),ncol=1)
            icobT <- icob1/icob1se
            icobTp <- 2*(1-pt(abs(icobT),tobs-nvar))

            scoreT<-cbind(ithinv,icob1)
            AIC<-(-2*(LLC))+cbind(nvar*2,nvar1*2,nvar2*2)
            AICC<-(-2*(LLC))+cbind(nvar*2*(tobs/(tobs-nvar*2-1)),nvar1*2*(tobs/(tobs-nvar1*2-1)),nvar2*2*(tobs/(tobs-nvar2*2-1)))
            BIC<-(-2*(LLC))+cbind(nvar*log(tobs),nvar1*log(tobs),nvar2*log(tobs))

            ## Dispersion adjuested s.d and p-values
            ph1     <-matrix(phi1,nrow=nvar1,ncol=1)
            ph2     <-matrix(phi2,nrow=nvar2,ncol=1)
            ph      <-rbind(ph1,ph2)
            icobsed <-matrix(diag(ithinv),ncol=1)
            icobsed <-sqrt(icobsed*ph)
            icobTd  <- icob1/icobsed
            icobTpd <- 2*(1-pt(abs(icobTd),tobs-nvar))

            ## Condition for end of itaration
            coc <- sum(abs(icob1-icob))
            cob <- sum(abs(icob1/presc))
            ttcon<-tollf(icob,icob1)

            icob<-icob1

            ## End of Newton Raphson
           # cat("Phi1 = ",phi1,"Phi2 = ",phi2,"\n")
           # print(data.frame(Var.Names=rrn,Coeff=icob1,s.e=icob1se,t_value=icobT,p_value=icobTp,row.names = NULL))
            if (ita == mxit & coc > cob)
                  stop("Function did not converge.... \n
                       Try by increasing maximum number of itarations....")

            if((coc<=cob)) {
               restab<-data.frame(Var.Names=rrn,Coeff=icob1,s.e=icob1se,t.value=icobT,p.value=icobTp,Adj.s.e=icobsed,Adj.p.value=icobTpd,row.names = NULL)
                  if(ppy==1){
                       cat("","\n")

        if(mtype==1){cat("-------- BIVP INDEPENDENT (Marginal/Marginal) MODEL------------","\n")}
        if(mtype==2){cat("-------- BIVP DEPENDENT (Marginal/Conditional) MODEL--------------------","\n")}
        if(mtype==3){cat("-------- ZTBIVP INDEPENDENT (Marginal/Marginal) MODEL---------","\n")}
        if(mtype==4){cat("-------- ZTBIVP DEPENDENT (Marginal/Conditional) MODEL---------","\n")}
        if(mtype==5){cat("-------- RTBIVP INDEPENDENT (Marginal/Marginal)----","\n")}
        if(mtype==6){cat("------- RTBIVP DEPENDENT (Marginal/Conditional)----","\n")}

                        cat("","\n")
                        cat("Log Likelihood =", LLC[1], "\n")
                        cat("AIC =", AIC[1], "\n")
                        cat("AICC =", AICC[1], "\n")
                        cat("BIC =", BIC[1], "\n")
                        cat("Deviance =", deva[1], "\n")
                        cat("Phi1 = ",phi1,"Phi2 = ",phi2,"\n")
                        cat("","\n")
                        cat("Iteration = ", ita, "\n")
                        cat("Parameter Estimates", "\n")
                        cat("----------------------------------------------------------------------", "\n")
                        is.num <- sapply(restab, is.numeric)
                        restab[is.num] <- lapply(restab[is.num], round, 6)
                        print(restab,row.names=FALSE)
                        cat("----------------------------------------------------------------------", "\n")
                        cat("\n")
                        cat("Overdispersion: Z-test for y1 & y2 =","\n" )
                        zzt<-dpztest(y1,y2,tobs)
                        is.num <- sapply(zzt, is.numeric)
                        zzt[is.num] <- lapply(zzt[is.num], round, 5)
                        print(zzt)
                        cat("\n")
                      #  cat("GOF test:","\n")
                       # print(gofres$GOF.Chisquare,row.names=FALSE)
                       # cat("\n")
                        #cat("GOF T1:","\n")
                        # cat("Good-ness-of-fit (T1) Overdispersion (T2):","\n")
                        # goft1<-t1t2(y1,y2,Lambdas,phi1,phi2,mtype)
                        # print(goft1)
                        # cat("\n")

                         cat("Good-ness-of-fit (T1) Overdispersion (T2):","\n")
                         goft2<-disper_fun(y1,y2,Lambdas,nvar,mtype)
                         print(goft2)
                         cat("\n")

                         cat("Pearson Chisquare GOF using predicted probability","\n")
                         gofres<-PearChiPred(y1,y2,Lambdas,nvar,mtype)
                         print(gofres$GOF.Chisquare,row.names=FALSE)
                         cat("\n")

                        if (is.null(phi)==TRUE) {mp<-nvar+1
                        } else {mp<-nvar}
                  } #ppy
                  conv<-1
                  listall<-list(restab,ithinv,ith,uth,LLC,tobs,nvar1,nvar2,nvar,AIC,BIC,AICC,Lambdas,phi1,phi2,mtype,conv,ttcon,gofres,
                                zzt,goft2,deva,mp,y1,y2,y1logll,y2logll)
                  names(listall)<-c("coeff","ithetainv","ith","utheta","logLik","N","nvar1","nvar2","nvar","AIC","BIC","AICC","Lambdas","Phi1",
                                    "Phi2","mtype","convg","Control","GOF.Chi","Disp.Ztest","T1","Deviance","nparam","y1","y2",
                                    "logliky1","logliky2")
                  mxit<-mxit
                  if(ppy==1){
                        print("Function Converged.....")
                  }
                  break
            }

      }   ## end loop for itaration"
      return(listall)
}    #    end of function
