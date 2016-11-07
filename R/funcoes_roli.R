






##' Amount of cloud estimatives functions.
##' @param data a data frame with all atmospherics variables
##' @return a vector with cloud amount estimatives
##' @export
##' @references Black JN (1956) The distribution of solar radiation over the Earth's surface. Arch Meteor Geophy B 7:165–189
##' 
    CQB <- function(data){
        a <- maxlim(with(data,0.34^2 + 4 * 0.458 * (0.803-K)),max_ = Inf)
        a <- ifelse(is.infinite(a),NA,a)
        maxlim( ( 0.34-sqrt(a) ) / (-2 * 0.458))
    }  # Quadratic regression of Black(1956)
    
##' Amount of cloud estimatives functions.
##' @param data a data frame with all atmospherics variables
##' @return a vector with cloud amount estimatives  
##' @export
##' @references Kasten F, Czeplak G (1980) Solar and terrestrial radiation dependent on the amount and type of cloud. Sol Energy 24:177–189
    CKC <- function(data){
        with(data,maxlim( (4/3*(1-K))^(1/3.4) ))   
        } # Kasten & Czeplack (1980)
    
##' Amount of cloud estimatives functions.
##' 
##' @param data a data frame with all atmospherics variables
##' 
##' @return a vector with cloud amount estimatives 
##' @export
##' @references Campbell GS (1985) Soil physics with BASIC. Elsevier, Amsterdam
    CCB <- function(data){ 
        with(data,maxlim( 2.33 - 3.33*K  ) ) 
        }      # Campbell (1985)
    
##' Amount of cloud estimatives functions.
##' @param data a data frame with all atmospherics variables
##' @param alt Site sea level heigth
##' @return a vector with cloud amount estimatives
##' @export
##' @references Konzelmann T, van de Wal RSW, Greuell W, Bintanja R, Henneken  EAC, Abe-Ouchi A 
##' (1994) Parameterization of global and longwave incoming radiation for the Greenland Ice Sheet. 
##' Glob Planet Chang 9:143–164
    CKZ <- function(data, alt = 88.){ 
            a <- 1./( 0.78*exp(-0.00085*alt))
          with(data,maxlim(  sqrt( (1-K)*a ) ))
    }     # Konzelmann (1994)
    
##' Amount of cloud estimatives functions.
##' @param data a data frame with all atmospherics variables
##' @return a vector with cloud amount estimatives
##' @export
##' @references Weishampel JF, Urban DL (1996) Coupling a spatially-explicit forest
##' gap model with a 3-D solar routine to simulate latitudinal effects.
##' Ecol Model 86:101–111
    CWU <- function(data){
        num <- with(data,252.7 - (Rg * 60*60*24/(4.19*10000)))
        den <- with(data,0.695*(Rpot * 60*60*24/(4.19*10000)))
        maxlim( 1+ (num/den) )
    } # Weishampel and Urban (1996)

##' Amount of cloud estimatives functions.
##' @param data a data frame with all atmospherics variables
##' @return a vector with cloud amount estimatives 
##' @export
##' @references Jegede, O. O., Ogolo, E.O., Aregbesola, T.O. (2006) Estimating net
##' radiation using routine meteorological data at a tropical location
##' in Nigeria. Int J Sustain Energy 25:107–115
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
##' @references Angstrom, A. (1915) A study of the radiation of the atmosphere.
##' Smithsonian Miscellaneous Collections 65(3)
    EAN <- function(data,func, 
                    coef1 = 0.83, coef2 = 0.18, coef3 = 0.067, 
                    adjust = FALSE){ 
            
            sigma <- 5.67051*10^(-8)
            
            if(adjust){
                
                
                tmp.nls <-
                    nls( Li/(sigma*Ta^4) ~ maxlim( coef1 - coef2*( 10^(-coef3*es))  ),
                                data = data, #na.action = "na.exclude",
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
##' @references Brunt D (1932) Notes on radiation in the atmosphere I. Q J Roy
##' Meteor Soc 58:389–420
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
##' @references Swinbank WC (1963) Long-wave radiation from clear skies. Q J Roy
##' Meteor Soc 89:339–348
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
##' @references Idso SB, Jackson RD (1969) Thermal radiation from the atmosphere.
##' J Geophys Res 74:5397–5403
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
##' @references Brutsaert W (1975) On a derivable formula for long-wave radiation
##' from clear skies. Water Resour Res 11:742–744
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
##' @references Idso SB (1981) A set of equations for full spectrum and 8- to 14-um
##'  and 10.5- to 12.5-um thermal radiation from cloudless skies.
##'  Water Resour Res 17:295–304
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
##' @references Konzelmann T, van de Wal RSW, Greuell W, Bintanja R, Henneken
##' EAC, Abe-Ouchi A (1994) Parameterization of global and
##' longwave incoming radiation for the Greenland Ice Sheet. Glob
##' Planet Chang 9:143–164
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
##' @references Prata AJ (1996) A new long-wave formula for estimating downward clearsky
##' radiation at the surface. Q J Roy Meteor Soc 122:1127–1151
    EPR <- function(data,func,
                    coef1 = 1, coef2 = 1.2, coef3 = 3,
                    adjust = FALSE) {
        
        sigma <- 5.67051*10^(-8)
        
        if(adjust){
            
            tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( 1 - ((coef1+46.5*es/Ta) * 
                                                               exp(-sqrt(coef2+coef3*46.5*es/Ta))) )  ,
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

    
    
##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @export
##' @references Niemala (2001); TODO 
ENM <- function(data,func,
                coef1 = 0.72, coef2 = 0.0009,
                adjust = FALSE) {
    
    sigma <- 5.67051*10^(-8)
    
    if(adjust){
        
        tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( coef1 + sign(es- 20.) *(coef2 * (es - 20.)))  ,
                        data = data, #na.action = "na.exclude",
                        start = list(coef1 = coef1,  coef2 = coef2 ) )
        
        new.coefs <- coef(tmp.nls) 
        new.emiss <- do.call(ENM,as.list(modifyList(formals(ENM),
                                                    c(list(data = data, func = func),
                                                      as.list(new.coefs)))))
        
        return(list(emiss = new.emiss, coefs = new.coefs))
        
    } else {
        
        return( with(data,maxlim( coef1 + sign(es- 20.)  * (coef2 * (es - 20.0)))  ))
        
    }

} ## Niemala (2001)



##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2,coef3 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @export
##' @references Aimi, D. (2017); TODO 
EAI <- function(data,func,
                coef1 = 0.52843, coef2 = -0.00820,coef3 = 24242.65010,
                adjust = FALSE) {
    
    sigma <- 5.67051*10^(-8)
    
    if(adjust){
        
        tmp.nls <- nls( Li/(sigma*Ta^4) ~ 
                            maxlim( (coef1 + ( coef2 * (es - 20.0)) + 
                                         coef3 / (sigma * Ta^4 * ifelse(Cover != "CN", Ta,es) )) 
                                    * exp(es/Ta) ), #
                        data = data,
                        start = list(coef1 = coef1,  coef2 = coef2, coef3 = coef3) )
        
        new.coefs <- coef(tmp.nls) 
        new.emiss <- do.call(EAI,as.list(modifyList(formals(EAI),
                                                    c(list(data = data, func = func),
                                                      as.list(new.coefs)))))
        
        return(list(emiss = new.emiss, coefs = new.coefs))
        
    } else {
        
        return( with(data,
                     maxlim( (coef1 + ( coef2 * (es - 20.0)) + 
                                  coef3 / (sigma * Ta^4 * ifelse(Cover != "CN", Ta,es) ))  
                             * exp(es/Ta) ) # 
                     )
                )
        
    }
    
} ## Aimi (2017)


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
##' @references Lhomme JP, Vacher JJ, Rocheteau A (2007) Estimating downward
##' long-wave radiation on the Andean Altiplano. Agr For Meteorol
##' 145:139–148
    ALH <- function(data,func, 
                    coef1 = 1.18, coef2 = 1.24, coef3 = -.34, coef4 = 1.37,
                    adjust=FALSE){
        
        if(adjust){
            
            emiss <- EBT(data = data,func = func, adjust = TRUE)
            sigma <- 5.67051*10^(-8)
             
            
            data$emiss <- emiss$emiss
            data$K[is.na(data$emiss)] <- NA
            
            suppressWarnings(
                tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( (coef1/coef2) * emiss *  (coef3*K+coef4) ) ,
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
##' @references Stockli R (2007) LBA-MIP driver data gap filling algorithms.
##' Unpublished. http://www.climatemodeling.org/lba-mip/LBAmipDriverDataFillingMethods.pdf
##' Accessed 7 May 2011
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
            data$emiss <- emiss$emiss
            data$K[is.na(data$emiss)] <- NA
            
            suppressWarnings(
                tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( emiss * (1+coef1*(1-K)^2 ) ) ,
                                data = data, na.action = "na.exclude",
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
##' @references Gabathuler M, Marty CA, Hanselmann KW (2001) Parameterization
##' of incoming longwave radiation in high-mountain environments.
##' Phys Geogr 22:99–114
    AGB <- function(data,func,
                    coef1 = 0.84, coef2 = 21.,
                    adjust = FALSE){

        sigma <- 5.67051*10^(-8)
        
        if(adjust){
            
            tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( (coef1*(rh-68))/(sigma*Ta^4) + 
                                                          (1.- coef2*K/Ta)^4 ) ,
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
##' @references Angstrom A (1915) A study of the radiation of the atmosphere.
##' Smithsonian Miscellaneous Collections 65(3)
    FAN <- function(data,func,
                    coef1 = 0.22,
                    adjust = FALSE){
        
      C <- do.call(func , args = list(data = data)) 
      
      sigma <- 5.67051*10^(-8)
      
      if(adjust){
          
          emiss <- EAN(data = data,func=func, adjust = TRUE)
       
          data$emiss <- emiss$emiss
          C[is.na(data$emiss)] <- NA

          suppressWarnings(
              tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( emiss * (1+coef1*C) ) ,
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
##' @references Brutsaert W (1975) On a derivable formula for long-wave radiation
##' from clear skies. Water Resour Res 11:742–744
    FBR <- function(data,func,
                    coef1 = 0.22,
                    adjust = FALSE){
        
        C <- do.call(func , args = list(data = data)) 
        
        sigma <- 5.67051*10^(-8)
        
        if(adjust){
            
            emiss <- EBR(data = data,func=func, adjust = TRUE)
       
            data$emiss <- emiss$emiss
            C[is.na(data$emiss)] <- NA
            
            suppressWarnings(
                tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( emiss * (1+coef1*C) ) ,
                                data = data, #na.action = "na.exclude",
                                start = list(coef1 = coef1) )
            )
            new.coefs <- coef(tmp.nls) 
            new.emiss <- predict(tmp.nls)
            
            return(list(emiss = new.emiss, coefs = new.coefs))
            
        } else {
            
            return(   maxlim( EBR(data)*(1+coef1*C)  )    )
            
        }
        
    }   ## Brutsaert (1982) ou (1975)**

##' Effective emissivity from atmosphere with cloud atenuation.
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @export
##' @references Friend AD, Stevens AK, Knox RG, Cannell MGR (1997) A processbased,
##' terrestrial biosphere model of ecosystem dynamics
##' (Hybrid v3.0). Ecol Model 95:249–287
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
##' @references Konzelmann T, van de Wal RSW, Greuell W, Bintanja R, Henneken
##' EAC, Abe-Ouchi A (1994) Parameterization of global and
##' longwave incoming radiation for the Greenland Ice Sheet. Glob
##' Planet Chang 9:143–164
    FKZ <- function(data,func, 
                    coef1 = 0.963,
                    adjust = FALSE){
        C <- do.call(func , args = list(data = data)) 
        
        sigma <- 5.67051*10^(-8)
        
        if(adjust){
            
            emiss <- EKZ(data = data,func = func, adjust = TRUE)
            
            data$emiss <- emiss$emiss
            C[is.na(data$emiss)] <- NA
            
            suppressWarnings(
                tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( emiss *(1.0-C^3)+ coef1 * C^3 ) ,
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
##' @references Idso SB, Jackson RD (1969) Thermal radiation from the atmosphere. 
##' J Geophys Res 74:5397–5403
    FIJ <- function(data,func,
                    coef1 = 0.22,
                    adjust = FALSE){
        
        C <- do.call(func , args = list(data = data)) 
        
        sigma <- 5.67051*10^(-8)
        
        if(adjust){
            
            emiss <- EIJ(data = data, func = func, adjust = TRUE)
            
            data$emiss <- emiss$emiss
            C[is.na(data$emiss)] <- NA
            
            suppressWarnings(
                tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( emiss *(1+coef1*C^2) ) ,
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
    # points(with(data, maxlim( 0.937097 - 0.167061*( 10^(-0.041560*es))) ),
    # col = "green", pch = 19, cex =0.3)
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
    # bb <- EAN_gc(data = data,func = "-",
    # coef1 = 0.93709740, coef2 = 0.16706067 ,coef3 = 0.04155977,adjust = TRUE)
    # plot(aa, pch = 19, cex = 0.3,col = "red", ylim = c(0.5,1.0))
    # points(bb$emiss, pch = 19, cex = 0.3,col = "blue", ylim = c(0.5,1.0))
    # 
    #####
    ## FOR OPTIMIZATION OF INITIAL PARAMETERS
    
    # coef1.0 <- with(data, min(Li/(sigma*Ta^4),na.rm = TRUE)) * 0.5
    # 
    # st <- 
    # lm(log10(Li/(sigma*Ta^4) - coef1.0 ) ~  (es)  ,
    #       data = data[complete.cases(data),] 
    #       # ,start = list(coef1 = coef1,coef2=coef2,coef3=coef3), 
    #       # ,trace = TRUE
    #    ) %>% 
    # coef() %>% as.numeric
    # new.start <- list(coef1 = coef1.0, coef2 = 10^(st[1]), coef3 = st[2])
