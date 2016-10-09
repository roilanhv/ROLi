
##' Amount of cloud estimatives functions.
##' @param data a data frame with all atmospherics variables
##' @return a vector with cloud amount estimatives
##' @export
    CQB <- function(data){
        a <- maxlim(with(data,0.34^2 + 4 * 0.458 * (0.803-K)),max_ = Inf)
        maxlim( ( 0.34-sqrt(a) ) / (-2 * 0.458))
    }  # Quadratic regression of Black(1956)
    
##' Amount of cloud estimatives functions.
##' 
##' @param data a data frame with all atmospherics variables
##' 
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
##' 
##' @param data a data frame with all atmospherics variables
##' @param alt Site sea level heigth
##' @return a vector with cloud amount estimatives
##' @export
    CKZ <- function(data, alt = 88.){ 
            a <- 1./( 0.78*exp(-0.00085*alt))
          with(data,maxlim(  sqrt( (1-K)*a ) ))
    }                                                     # Konzelmann (1994)
    
##' Amount of cloud estimatives functions.
##' 
##' @param data a data frame with all atmospherics variables
##' 
##' @return a vector with cloud amount estimatives
##' @export
    CWU <- function(data){
        num <- with(data,252.7 - (Rg * 60*60*24/(4.19*10000)))
        den <- with(data,0.695*(Rpot * 60*60*24/(4.19*10000)))
        maxlim( 1+ (num/den) )
    } # Weishampel and Urban (1996)

##' Amount of cloud estimatives functions.
##' 
##' @param data a data frame with all atmospherics variables
##' 
##' @return a vector with cloud amount estimatives 
##' @export
    CJG <- function(data){ 
        with(data,maxlim( ifelse(K < 0.9 , 1.1-K,2*(1-K)) ))
        } # Jedge (2006)

##
#/////////////////////////////////////////////////////////////////////////////////////////////////
# Parametrizações de ...
#---- EMISSIVIDADE 

##' Emissivity from atmosphere
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @export
    EAN <- function(data,func){ with(data,maxlim( 0.83 - 0.18*(10^(-0.067*es)) )) }                                   
    ## Angstrom (1915)

##' Emissivity from atmosphere
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @export
    EBR <- function(data,func){ with(data,maxlim( 0.51 + 0.066*sqrt(es) )) }                                         
    ## Brunt (1932)

##' Emissivity from atmosphere
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @export
    ESW <- function(data,func){ 
        with(data,maxlim(0.0000092*Ta^2)) 
        }     ## Swinbank (1963)

##' Emissivity from atmosphere
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @export
    EIJ <- function(data,func){ 
        with(data,maxlim(1 - 0.261 * exp(-0.000777 * (273 - Ta)^2))) 
        } ## Idso & Jackson (1969)
        # NOTE: (Ta - 273.15) OU (273 - Ta)

##' Emissivity from atmosphere
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @export
    EBT <- function(data,func){ 
        with(data,maxlim(1.24*(es/Ta)^(1/7) ))
        }   ## Brutsaert (1934)
    

##' Emissivity from atmosphere
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @export
    EID <- function(data,func){ 
        with(data,maxlim(0.7 + 0.0000595*es*exp(1500/Ta) )) 
        }   ## Idso (1981)

##' Emissivity from atmosphere
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @export
    EKZ <- function(data,func) {
        with(data,maxlim( 0.23 + 0.484*(es/Ta)^(1/8) )) 
        }   ## Konzelmann (1994)

##' Emissivity from atmosphere
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @export
    EPR <- function(data,func){ 
        with(data,maxlim(1 - ((1+46.5*es/Ta) * exp(-sqrt(1.2+3*46.5*es/Ta)))))
        } ## Prata (1996)

#/////////////////////////////////////////////////////////////////////////////////////////////////
# Parametrizações de ...
#---- EMISSIVIDADE EFETIVA COM INDICE DE ATENUAÇÃO

##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @export
    ALH <- function(data,func){ 
        1.18/1.24*EBT(data) * with(data,(-.34*K+1.37)) 
        } ## Lhomme (2007)
    
##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
##' @export
    ABM <- function(data,func){ 
        EID(data)*with(data,(1+0.3*(1-K)^2 ))
        } ## Stockli (2007)
    
##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @export
    AGB <- function(data,func){
        sigma <- 5.67051*10^(-8)
        a <- with(data,(0.84*(rh-68))/(sigma*Ta^4) )
        b <- with(data,(1- 21*K/Ta)^4)
        maxlim(a+b)
    } ## Gabathuler (2001)
    
#/////////////////////////////////////////////////////////////////////////////////////////////////
# Parametrizações de ...
#---- EMISSIVIDADE EFETIVA COM COBERTURA DE NUVENS
#    
    
##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @export
    FAN <- function(data,func){
      C <- do.call(func , args = list(data = data)) 
      maxlim( EAN(data)*(1+0.22*C)  )
    }   ## Angstrom (1915)

##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @export
    FBR <- function(data,func){ 
      C <- do.call(func , args = list(data = data)) 
      maxlim( EBR(data) * (1 + 0.22*C) ) 
      }                                                  ## Brutsaert (1982) ou (1975)**

##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @export
        # FMU <- function()
    FHY <- function(data,func){ 
      C <- do.call(func , args = list(data = data)) 
      maxlim(0.69*(1-C^6) + 0.979*C^4)
      }                                                 ## HYBRID (2009)
    
##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @export
    FKZ <- function(data,func){
        C <- do.call(func , args = list(data = data)) 
        EKZ(data)*(1.0-C^3)+0.963*C^3
    }                                                   ## Konzelmann (1994)
    
##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @export
    FIJ <- function(data,func){
        C <- do.call(func , args = list(data = data)) 
        EIJ(data)*(1+0.22*C^2)
    }                                                   
    

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




    
