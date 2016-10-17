






##' Amount of cloud estimatives functions.
##' @param data a data frame with all atmospherics variables
##' @return a vector with cloud amount estimatives
##' @export
    CQB <- function(data){
        a <- maxlim(with(data,0.34^2 + 4 * 0.458 * (0.803-K)),max_ = Inf)
        maxlim( ( 0.34-sqrt(a) ) / (-2 * 0.458))
    }  # Quadratic regression of Black(1956)
    
##' Amount of cloud estimatives functions.
##' @param data a data frame with all atmospherics variables
##' @return a vector with cloud amount estimatives  
##' @export
    CKC <- function(data){
        with(data,maxlim( (4/3*(1-K))^(1/3.4) ))   
        } # Kasten & Czeplack (1980)
    
##' Amount of cloud estimatives functions.
##' 
##' @param data a data frame with all atmospherics variables
##' 
##' @return a vector with cloud amount estimatives 
##' @export
    CCB <- function(data){ 
        with(data,maxlim( 2.33 - 3.33*K  ) ) 
        }      # Campbell (1985)
    
##' Amount of cloud estimatives functions.
##' @param data a data frame with all atmospherics variables
##' @param alt Site sea level heigth
##' @return a vector with cloud amount estimatives
##' @export
    CKZ <- function(data, alt = 88.){ 
            a <- 1./( 0.78*exp(-0.00085*alt))
          with(data,maxlim(  sqrt( (1-K)*a ) ))
    }                                                     # Konzelmann (1994)
    
##' Amount of cloud estimatives functions.
##' @param data a data frame with all atmospherics variables
##' @return a vector with cloud amount estimatives
##' @export
    CWU <- function(data){
        num <- with(data,252.7 - (Rg * 60*60*24/(4.19*10000)))
        den <- with(data,0.695*(Rpot * 60*60*24/(4.19*10000)))
        maxlim( 1+ (num/den) )
    } # Weishampel and Urban (1996)

##' Amount of cloud estimatives functions.
##' @param data a data frame with all atmospherics variables
##' @return a vector with cloud amount estimatives 
##' @export
    CJG <- function(data){ 
        with(data,maxlim( ifelse(K < 0.9 , 1.1-K,2*(1-K)) ))
        } # Jedge (2006)

##
#/////////////////////////////////////////////////////////////////////////////////////////////////
#                                        PARAMETRIZAÇÔES DE 
#----                                       EMISSIVIDADE 

##' Emissivity from atmosphere
##' @param data Data frame with all atmospherics variables
##' @param func Function for amount of cloud 
##' @param coef1,coef2,coef3 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @export
    EAN <- function(data,func, 
                    coef1 = 0.83, coef2 = 0.18, coef3 = 0.067, 
                    adjust = FALSE){ 
            
            sigma <- 5.67051*10^(-8)
            
            if(adjust){
                
                tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( coef1 - coef2*( 10^(-coef3*es))  ),
                                data = data, na.action = "na.exclude",
                                start = list(coef1 = coef1,coef2=coef2,coef3=coef3) )
                
                new.coefs <- coef(tmp.nls) 
                new.emiss <- do.call(EAN,as.list(modifyList(formals(EAN),
                                                    c(list(data = data, func = func),
                                                      as.list(new.coefs)))))
           
                
                return(list(emiss = new.emiss, coefs = new.coefs))
                
            } else {
                
                return( with(data,maxlim( coef1 - coef2*(10^(-coef3*es)) )) )
                
            }
            
    }   ## Angstrom (1915)

    
##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @export
    EBR <- function(data,func,
                    coef1 = 0.51, coef2 = 0.066,
                    adjust = FALSE){
       
        sigma <- 5.67051*10^(-8)
        
        if(adjust){
            
            tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( coef1 + coef2*sqrt(es)  ),
                            data = data, na.action = "na.exclude",
                            start = list(coef1 = coef1,coef2=coef2) )
            
            new.coefs <- coef(tmp.nls) 
            new.emiss <- do.call(EBR,as.list(modifyList(formals(EBR),
                                                        c(list(data = data, func = func),
                                                          as.list(new.coefs)))))
            
            return(list(emiss = new.emiss, coefs = new.coefs))
            
            
        } else {
        
            return(with(data,maxlim( coef1 + coef2*sqrt(es) )) )
            
        }
        
    }    ## Brunt (1932)

##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @export
    ESW <- function(data,func,
                    coef1 = 0.0000092,
                    adjust = FALSE){ 
        
        sigma <- 5.67051*10^(-8)
        
        if(adjust){
            
            tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim(coef1*Ta^2),
                            data = data, na.action = "na.exclude",
                            start = list(coef1 = coef1) )
            
            new.coefs <- coef(tmp.nls) 
            new.emiss <- do.call(ESW,as.list(modifyList(formals(ESW),
                                                        c(list(data = data, func = func),
                                                          as.list(new.coefs)))))
            
            return(list(emiss = new.emiss, coefs = new.coefs))
            
        } else {
            
            return( with(data,maxlim(coef1*Ta^2)) )    
            
        }
        
    }     ## Swinbank (1963)

##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @export
    EIJ <- function(data,func, 
                    coef1 = 0.261, coef2 = -0.000777,
                    adjust = FALSE){
       
         sigma <- 5.67051*10^(-8)
        
        if(adjust){
            
            tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim(1 - coef1 * exp(coef2 * (273 - Ta)^2)) ,
                            data = data, na.action = "na.exclude",
                            start = list(coef1 = coef1, coef2 = coef2) )
            
            new.coefs <- coef(tmp.nls) 
            new.emiss <- do.call(EIJ,as.list(modifyList(formals(EIJ),
                                                        c(list(data = data, func = func),
                                                          as.list(new.coefs)))))
            
            return(list(emiss = new.emiss, coefs = new.coefs))
            
        } else {
            return(with(data,maxlim(1 - coef1 * exp(coef2 * (273 - Ta)^2)))     )
        }
        
    } ## Idso & Jackson (1969)
    # NOTE: (Ta - 273.15) OU (273 - Ta)

##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud
##' @param coef1,coef2 Scheme coeficients  
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @export
    EBT <- function(data,func,
                    coef1 = 1.24, coef2 = 1/7, 
                    adjust = FALSE){ 
        
        sigma <- 5.67051*10^(-8)
        
        if(adjust){
            
            tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( coef1*(es/Ta)^(coef2) ) ,
                            data = data, na.action = "na.exclude",
                            start = list(coef1 = coef1, coef2 = coef2) )
            
            new.coefs <- coef(tmp.nls) 
            new.emiss <- do.call(EBT,as.list(modifyList(formals(EBT),
                                                        c(list(data = data, func = func),
                                                          as.list(new.coefs)))))
            
            return(list(emiss = new.emiss, coefs = new.coefs))
            
        } else {
            
            return(with(data,maxlim(coef1*(es/Ta)^(coef2) ))    )
            
        }
        
    }   ## Brutsaert (1934)
    

##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @export
    EID <- function(data,func,
                    coef1 = 0.7, coef2 = 0.0000595,
                    adjust = FALSE){ 
        
        sigma <- 5.67051*10^(-8)
        
        if(adjust){
            
            tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( coef1 + coef2*es*exp(1500/Ta) )  ,
                            data = data, na.action = "na.exclude",
                            start = list(coef1 = coef1, coef2 = coef2) )
            
            new.coefs <- coef(tmp.nls) 
            new.emiss <- do.call(EID,as.list(modifyList(formals(EID),
                                                        c(list(data = data, func = func),
                                                          as.list(new.coefs)))))
            
            return(list(emiss = new.emiss, coefs = new.coefs))
            
        } else {
            
            return( with(data,maxlim( coef1 + coef2*es*exp(1500/Ta) ))  )
            
        }
        
        
        }   ## Idso (1981)

##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2,coef3 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @export
    EKZ <- function(data,func,
                    coef1 = 0.23, coef2 = 0.484, coef3 = 1/8,
                    adjust = FALSE) {
        
        sigma <- 5.67051*10^(-8)
        
        if(adjust){
            
            tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( coef1 + coef2*(es/Ta)^(coef3) )  ,
                            data = data, na.action = "na.exclude",
                            start = list(coef1 = coef1, coef2 = coef2) )
            
            new.coefs <- coef(tmp.nls) 
            new.emiss <- do.call(EKZ,as.list(modifyList(formals(EKZ),
                                                        c(list(data = data, func = func),
                                                          as.list(new.coefs)))))
            
            return(list(emiss = new.emiss, coefs = new.coefs))
            
        } else {
            
            return( with(data,maxlim( coef1 + coef2*(es/Ta)^(coef3) ))  )
            
        }
        
        
        }   ## Konzelmann (1994)

##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2,coef3 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @export
    EPR <- function(data,func,
                    coef1 = 1, coef2 = 1.2, coef3 = 3,
                    adjust = FALSE) {
        
        sigma <- 5.67051*10^(-8)
        
        if(adjust){
            
            tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( 1 - ((coef1+46.5*es/Ta) * exp(-sqrt(coef2+coef3*46.5*es/Ta))) )  ,
                            data = data, na.action = "na.exclude",
                            start = list(coef1 = coef1,  coef3 = coef3 ) )
            
            new.coefs <- coef(tmp.nls) 
            new.emiss <- do.call(EPR,as.list(modifyList(formals(EPR),
                                                        c(list(data = data, func = func),
                                                          as.list(new.coefs)))))
            
            return(list(emiss = new.emiss, coefs = new.coefs))
            
        } else {
            
            return( with(data,maxlim(1 - ((coef1+46.5*es/Ta) * exp(-sqrt(coef2+coef3*46.5*es/Ta)))))  )
            
        }
        
        
        } ## Prata (1996)

#/////////////////////////////////////////////////////////////////////////////////////////////////
# Parametrizações de ...
#---- EMISSIVIDADE EFETIVA COM INDICE DE ATENUAÇÃO

##' Effective emissivity from atmosphere with cloud atenuation
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2,coef3,coef4 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @export
    ALH <- function(data,func, 
                    coef1 = 1.18, coef2 = 1.24, coef3 = -.34, coef4 = 1.37,
                    adjust=FALSE){
        
        if(adjust){
            
            emiss <- EBT(data = data,func = func, adjust = TRUE)
            sigma <- 5.67051*10^(-8)
            # emiss$emiss <- do.call("EBT",
            #                        args = as.list(modifyList(formals(EBT), 
            #                                                  c(list(data = data, func = "-"), 
            #                                                    as.list(emiss$coefs)))
            #                                       )
            #                        )
            #         
            suppressWarnings(
                tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( (coef1/coef2) * emiss$emiss *  (coef3*K+coef4) ) ,
                            data = data, #na.action = "na.exclude",
                            start = list(coef1 = coef1, coef3 = coef3 ) )
            )
            new.coefs <- coef(tmp.nls) 
            new.emiss <- predict(tmp.nls)
            
            return(list(emiss = new.emiss, coefs = new.coefs))
            
        } else {
            
            return(maxlim(coef1/coef2*EBT(data) * with(data,(coef3*K+coef4)) ))
            
        }
                
        
        
        } ## Lhomme (2007)
    
    
##' Effective emissivity from atmosphere with cloud atenuation
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @return a vector with emissivity estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
##' @import stats
##' @import utils
##' @export
    ABM <- function(data,func,
                    coef1 = 0.3,
                    adjust = FALSE){ 
        
        if(adjust){
            
            sigma <- 5.67051*10^(-8)
            
            emiss <- EID(data = data,func = func , adjust = TRUE)
            # emiss$emiss <- do.call("EID",
            #                        args = as.list(modifyList(formals(EID), 
            #                                                  c(list(data = data, func = "-"), 
            #                                                    as.list(emiss$coefs)))
            #                        )
            # )
            
            suppressWarnings(
                tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( emiss$emiss * (1+coef1*(1-data$K)^2 ) ) ,
                                data = data, 
                                start = list(coef1 = coef1) )
            )
            new.coefs <- coef(tmp.nls) 
            new.emiss <- predict(tmp.nls)
            
            return(list(emiss = new.emiss, coefs = new.coefs))
            
        } else {
            return( maxlim( EID(data) *with(data,(1+coef1*(1-K)^2 ))    ) )
        }
        
        
        
    } ## Stockli (2007)
    
##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @export
    AGB <- function(data,func,
                    coef1 = 0.84, coef2 = 21.,
                    adjust = FALSE){

        sigma <- 5.67051*10^(-8)
        
        if(adjust){
            
            tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( (coef1*(rh-68))/(sigma*Ta^4) + (1.- coef2*K/Ta)^4 ) ,
                            data = data, na.action = "na.exclude",
                            start = list(coef1 = coef1, coef2 = coef2) )
            
            new.coefs <- coef(tmp.nls) 
            new.emiss <- predict(tmp.nls)
            
            return(list(emiss = new.emiss, coefs = new.coefs))
            
        } else {
        
            a <- with(data,(coef1*(rh-68))/(sigma*Ta^4) )
            b <- with(data,(1.- coef2*K/Ta)^4)
            return(  maxlim(a+b)    )
        }
        
    } ## Gabathuler (2001)
    
#/////////////////////////////////////////////////////////////////////////////////////////////////
# Parametrizações de ...
#---- EMISSIVIDADE EFETIVA COM COBERTURA DE NUVENS
#    
    
##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @export
    FAN <- function(data,func,
                    coef1 = 0.22,
                    adjust = FALSE){
        
      C <- do.call(func , args = list(data = data)) 
      
      sigma <- 5.67051*10^(-8)
      
      if(adjust){
          
          emiss <- EAN(data = data,func=func, adjust = TRUE)
          # emiss$emiss <- do.call("EAN",
          #                        args = as.list(modifyList(formals(EAN),
          #                                                  c(list(data = data, func = "-"),
          #                                                    as.list(emiss$coefs)))
          #                        )
          # )

          suppressWarnings(
              tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( emiss$emiss * (1+coef1*C) ) ,
                              data = data,# na.action = "na.exclude",
                              start = list(coef1 = coef1) )
          )
          new.coefs <- coef(tmp.nls) 
          new.emiss <- predict(tmp.nls)
          
          return(list(emiss = new.emiss, coefs = new.coefs))
          
      } else {
          
          maxlim( EAN(data)*(1+coef1*C)  )    
          
      }
      
    }   ## Angstrom (1915)

##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @export
    FBR <- function(data,func,
                    coef1 = 0.22,
                    adjust = FALSE){
        
        C <- do.call(func , args = list(data = data)) 
        
        sigma <- 5.67051*10^(-8)
        
        if(adjust){
            
            emiss <- EBR(data = data,func=func, adjust = TRUE)
            # emiss$emiss <- do.call("EBR",
            #                        args = as.list(modifyList(formals(EBR), 
            #                                                  c(list(data = data, func = "-"), 
            #                                                    as.list(emiss$coefs)))
            #                        )
            # )
            
            suppressWarnings(
                tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( emiss$emiss * (1+coef1*C) ) ,
                                data = data, #na.action = "na.exclude",
                                start = list(coef1 = coef1) )
            )
            new.coefs <- coef(tmp.nls) 
            new.emiss <- predict(tmp.nls)
            
            return(list(emiss = new.emiss, coefs = new.coefs))
            
        } else {
            
            return(   maxlim( EBR(data)*(1+coef1*C)  )    )
            
        }
        
    }                                                  ## Brutsaert (1982) ou (1975)**

##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @export
    FHY <- function(data,func,
                    coef1 = 0.69, coef2 = 0.979,
                    adjust = FALSE){
        
        C <- do.call(func , args = list(data = data)) 
        
        sigma <- 5.67051*10^(-8)
        
        if(adjust){
            
            suppressWarnings(
                tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( coef1 *(1-C^6) + coef2 * C^4 ) ,
                                data = data, #na.action = "na.exclude",
                                start = list(coef1 = coef1, coef2 = coef2) )
            )
            new.coefs <- coef(tmp.nls) 
            new.emiss <- do.call(FHY,as.list(modifyList(formals(FHY),
                                                        c(list(data = data, func = func),
                                                          as.list(new.coefs)))))
            
            return(list(emiss = new.emiss, coefs = new.coefs))
            
        } else {
            
            return(maxlim( coef1 *(1-C^6) + coef2 * C^4))
        }
        
    }     ## HYBRID (2009)
    
##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @export
    FKZ <- function(data,func, 
                    coef1 = 0.963,
                    adjust = FALSE){
        C <- do.call(func , args = list(data = data)) 
        
        sigma <- 5.67051*10^(-8)
        
        if(adjust){
            
            emiss <- EKZ(data = data,func = func, adjust = TRUE)
            # emiss$emiss <- do.call("EKZ",
            #                        args = as.list(modifyList(formals(EKZ), 
            #                                                  c(list(data = data, func = "-"), 
            #                                                    as.list(emiss$coefs)))
            #                        )
            # )
            
            suppressWarnings(
                tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( emiss$emiss *(1.0-C^3)+ coef1 * C^3 ) ,
                                data = data,# na.action = "na.exclude",
                                start = list(coef1 = coef1) )
            )
            new.coefs <- coef(tmp.nls) 
            new.emiss <- predict(tmp.nls)
            
            return(list(emiss = new.emiss, coefs = new.coefs))
            
            
        } else {
            return(   EKZ(data)*(1.0-C^3)+ coef1 * C^3    )
        }
        
    }      ## Konzelmann (1994)
    
##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @export
    FIJ <- function(data,func,
                    coef1 = 0.22,
                    adjust = FALSE){
        
        C <- do.call(func , args = list(data = data)) 
        
        sigma <- 5.67051*10^(-8)
        
        if(adjust){
            
            emiss <- EIJ(data = data, func = func, adjust = TRUE)
            # emiss$emiss <- do.call("EIJ",
            #                        args = as.list(modifyList(formals(EIJ), 
            #                                                  c(list(data = data, func = "-"), 
            #                                                    as.list(emiss$coefs)))
            #                        )
            # )
            # 
            suppressWarnings(
                tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( emiss$emiss *(1+coef1*C^2) ) ,
                                data = data,#na.action = "na.exclude",
                                start = list(coef1 = coef1) )
            )
            new.coefs <- coef(tmp.nls) 
            new.emiss <- predict(tmp.nls)
            
            return(list(emiss = new.emiss, coefs = new.coefs))
            
        } else {
            
            return(  EIJ(data)*(1+coef1*C^2) )    
            
        }
        
    }                                                   
    
#############################
### OTHERS FUNCTIONS
### 

##' Incoming solar radiation atenuattion (K)
##' 
##' @param Rg global radiation serie.
##' @param dates A vector with dates of time serie, same lenght than Rg
##' @param lon longitude from local analisys 
##' @param lat latitude from local analisys
##' @param timezone Local time diference with GMT (-1 for fluxes measurement)
##' @return Vector with K index columns
##' @author Roilan Hernandez
##' @importFrom stats setNames
##' @importFrom dplyr %>% mutate select 
##' @importFrom plyr . 
##' @export
    kloudines <- function(Rg, dates,lon=-53.76,lat=-29.72,timezone=-4){
  
        if(lon==-53.76 & lat==-29.72) 
            warning("Latitude e Longitude de Santa Maria, RS, Brazil",
                    call. = TRUE,immediate. = TRUE)
        
        Rpot <- PotRad(date = dates,lon = lon, lat = lat,timezone = timezone)
        
        K <- Rg/Rpot
        K <- ifelse(is.infinite(K),0.0, K ) 
        K <- maxlim(K)
        K <- ifelse(is.na(K), 0.0,K)
 
        return(K)
    }


################################
#####     GARBAGE         ######

    
    # plot(with(data,Li/(sigma*Ta^4) ), pch = 19, cex =0.3, col = "red")
    # points(new.Li, pch = 19, cex =0.3)
    # points(with(data, maxlim( 0.937097 - 0.167061*( 10^(-0.041560*es))) ), col = "green", pch = 19, cex =0.3)
    # 
    # data$new.Li <- predict(tmp.nls) * sigma * data$Ta^4 
    # data$old.Li <- with(data,maxlim( coef1 - coef2*( 10^(-coef3*es))) * sigma * Ta^4 )
    # 
    # timePlot(data, c("Li","new.Li", "old.Li"), group = TRUE, lty = 1, avg.time = "day", lwd =2)
    # 
    # scatterPlot(data, x = "Li", y = "new.Li", linear = TRUE, mod.line = TRUE, avg.time = "day")
    # 
    # 
    # str(new.Li); str(with(data,Li/(sigma*Ta^4) ))
    
    
    # aa <- EAN_gc(data = data,func = "-")
    # bb <- EAN_gc(data = data,func = "-",coef1 = 0.93709740, coef2 = 0.16706067 ,coef3 = 0.04155977,adjust = TRUE)
    # plot(aa, pch = 19, cex = 0.3,col = "red", ylim = c(0.5,1.0))
    # points(bb$emiss, pch = 19, cex = 0.3,col = "blue", ylim = c(0.5,1.0))
    # 
    
