#' Function to check function arguments
#' @param dmat :  Data frame (y1,y2,X1,X2)
#' @param mxit :  Maximum number of iterations
#' @param mtype : Model type: 1= Double Poisson with proposed NC; 2= Efron's Double Poisson Model-1; 3= Efron's Double Poisson Model-2
#' @return returns the data matrix after casewise deletion incase of missing values
#' #@importFrom("stats", "complete.cases", "na.omit")
#' @importFrom stats complete.cases
#' @importFrom stats na.omit
poichkfun<-function(dmat, mxit=150, mtype){

      dmat<-as.matrix(dmat)
      ##### Start function argument checks ######
      if(missing(dmat))
            stop("Provide data file...")

    #  if((ncol(dmat) %% 2 == 0)=="FALSE")
     #       stop("Data file should have even number of columns...
      #           Ex.: Y1,Y2,X11,X12,...,X1p,X21,X22,...,X2p.")

      if(nrow(dmat[!complete.cases(dmat),])>0){
            warning("....Missing value(s) present in data file....
                    ....Number of rows deletd with missing values is = ",nrow(dmat[!complete.cases(dmat),]))
            dmat<-na.omit(dmat)
      }

      if(missing(mtype))
            stop("Need to specify the model type....Valid types are
                 1 = Bivariate Poisson Model (Marginal/Marginal)
                 2 = Bivariate Poisson Model (Marginal/Conditional)
                 3 = Bivariate BIVZTP  Model (Marginal/Marginal)
                 4 = Bivariate BIVZTP  Model (Marginal/Conditional)
                 5 = Bivariate BIVRTP  Model (Marginal/Marginal)
                 6 = Bivariate BIVRTP  Model (Marginal/Conditional)")

      if(mtype< 1 | mtype >6)
            stop("Valid types are....
                 1 = Bivariate Poisson Model (Marginal/Marginal)
                 2 = Bivariate Poisson Model (Marginal/Conditional)
                 3 = Bivariate BIVZTP  Model (Marginal/Marginal)
                 4 = Bivariate BIVZTP  Model (Marginal/Conditional)
                 5 = Bivariate BIVRTP  Model (Marginal/Marginal)
                 6 = Bivariate BIVRTP  Model (Marginal/Conditional)")

      if(is.int(mtype)=="FALSE")
            stop("Need to specify the model type....Valid types are
                 1 = Bivariate Poisson Model (Marginal/Marginal)
                 2 = Bivariate Poisson Model (Marginal/Conditional)
                 3 = Bivariate BIVZTP  Model (Marginal/Marginal)
                 4 = Bivariate BIVZTP  Model (Marginal/Conditional)
                 5 = Bivariate BIVRTP  Model (Marginal/Marginal)
                 6 = Bivariate BIVRTP  Model (Marginal/Conditional)")

      return(dmat)
}

