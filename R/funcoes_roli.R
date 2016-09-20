
##' Amount of cloud estimatives functions.
##' 
##' @param data a data frame with all atmospherics variables
##' 
##' @return a vector with cloud amount estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
    CQB <- function(data){
        a <- maxlim(with(data,0.34^2 + 4 * 0.458 * (0.803-K)),max_ = Inf)
        maxlim( ( 0.34-sqrt(a) ) / (-2 * 0.458))
    }  # Quadratic regression of Black(1956)
    
##' Amount of cloud estimatives functions.
##' 
##' @param data a data frame with all atmospherics variables
##' 
##' @return a vector with cloud amount estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch    
    CKC <- function(data){ with(data,maxlim( (4/3*(1-K))^(1/3.4) ))   } # Kasten & Czeplack (1980)
    
##' Amount of cloud estimatives functions.
##' 
##' @param data a data frame with all atmospherics variables
##' 
##' @return a vector with cloud amount estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch    
    CCB <- function(data){ with(data,maxlim( 2.33 - 3.33*K  ) )  }      # Campbell (1985)
    
##' Amount of cloud estimatives functions.
##' 
##' @param data a data frame with all atmospherics variables
##' @param alt Site sea level heigth
##' @return a vector with cloud amount estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
    CKZ <- function(data, alt = 88.){ 
            a <- 1./( 0.78*exp(-0.00085*alt))
          with(data,maxlim(  sqrt( (1-K)*a ) ))
    }                                                     # Konzelmann (1994)
    
##' Amount of cloud estimatives functions.
##' 
##' @param data a data frame with all atmospherics variables
##' 
##' @return a vector with cloud amount estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
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
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch    
    CJG <- function(data){ with(data,maxlim( ifelse(K < 0.9 , 1.1-K,2*(1-K)) ))} # Jedge (2006)

##
#/////////////////////////////////////////////////////////////////////////////////////////////////
# Parametrizações de ...
#---- EMISSIVIDADE 

##' Emissivity from atmosphere
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
    EAN <- function(data,func){ with(data,maxlim( 0.83 - 0.18*(10^(-0.067*es)) )) }                                   ## Angstrom (1915)

##' Emissivity from atmosphere
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
    EBR <- function(data,func){ with(data,maxlim( 0.51 + 0.066*sqrt(es) )) }                                          ## Brunt (1932)

##' Emissivity from atmosphere
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
    ESW <- function(data,func){ 
        with(data,maxlim(0.0000092*Ta^2)) 
        }     ## Swinbank (1963)

##' Emissivity from atmosphere
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
    EIJ <- function(data,func){ 
        with(data,maxlim(1 - 0.261 * exp(-0.000777 * (Ta - 273.15)^2))) 
        } ## Idso & Jackson (1969)
        # NOTE: (Ta - 273.15) OU (273 - Ta)

##' Emissivity from atmosphere
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
    EBT <- function(data,func){ 
        with(data,maxlim(1.24*(es/Ta)^(1/7) ))
        }   ## Brutsaert (1934)
    

##' Emissivity from atmosphere
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
    EID <- function(data,func){ 
        with(data,maxlim(0.7 + 0.0000595*es*exp(1500/Ta) )) 
        }   ## Idso (1981)

##' Emissivity from atmosphere
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
    EKZ <- function(data,func) {
        with(data,maxlim( 0.23 + 0.484*(es/Ta)^(1/8) )) 
        }   ## Konzelmann (1994)

##' Emissivity from atmosphere
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
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
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
    ALH <- function(data,func){ 
        1.18/1.24*EBT(data) * with(data,(-.34*K+1.37)) 
        } ## Lhomme (2007)
    
##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
    ABM <- function(data,func){ 
        EID(data)*with(data,(1+0.3*(1-K)^2 ))
        } ## Stockli (2007)
    
##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
    AGB <- function(data,func){
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
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
    FAN <- function(data,func){
      C <- do.call(func , args = list(data = data)) 
      maxlim( EAN(data)*(1+0.22*C)  )
    }   ## Angstrom (1915)

##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
    FBR <- function(data,func){ 
      C <- do.call(func , args = list(data = data)) 
      maxlim( EBR(data) * (1 + 0.22*C) ) 
      }                                                  ## Brutsaert (1982) ou (1975)**

##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
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
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
    FKZ <- function(data,func){
        C <- do.call(func , args = list(data = data)) 
        EKZ(data)*(1.0-C^3)+0.963*C^3
    }                                                   ## Konzelmann (1994)
    
##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @return a vector with emissivity estimatives
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
    FIJ <- function(data,func){
        C <- do.call(func , args = list(data = data)) 
        EIJ(data)*(1+0.22*C^2)
    }                                                   
    
    
    
#
#/////////////////////////////////////////////////////////////////////////////////////////////////
#  FUNÇÕES NECESSÁRIAS
#
    maxlim <- function(i,max_=1,min_=0){ 
        sapply(i,function(i) min(max(i,min_),max_) ) 
        }
   
##' Incoming solar radiation atenuattion (K)
##' 
##' @param Rg global radiation serie.
##' @param dates A vector with dates of time serie, same lenght than Rg
##' @param lon longitude from local analisys 
##' @param lat latitude from local analisys
##' @param timezone Local time diference with GMT (-1 for fluxes measurement)
##' @return Vector with K index columns
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
##' @importFrom stats setNames
##' @importFrom dplyr %>% mutate select 
##' @importFrom plyr . 
##' @importFrom REddyProc fCalcPotRadiation
    kloudines <- function(Rg, dates,lon=-53.76,lat=-29.72,timezone=-4){
  
        if(lon==-53.76 & lat==-29.72) 
            warning("Latitude e Longitude de Santa Maria",call. = TRUE,immediate. = TRUE)
        
        Rpot <- fCalcPotRadiation(DoY.V.n = format(dates,"%j") %>% as.numeric,
                                  Hour.V.n = format(dates,"%H") %>% as.numeric,
                                  Lat_deg.n = lat,
                                  Long_deg.n = lon,
                                  TimeZone_h.n = timezone,
                                  useSolartime.b = TRUE)
        
        K <- Rg/Rpot %>% 
                ifelse(is.infinite(.),0.0, . ) %>%
                maxlim(.) %>%
                ifelse(is.na(.), 0.0,.)
 
        return(K)
    }

##' Potencial Radiation from date vector
##' 
##' @param date A vector with data
##' @param lon Longitude from local analisys 
##' @param lat Latitude from local analisys
##' @param timezone Local time diference with GMT (-1 for fluxes measurement)
##' @return Vector with
##' @author Roilan Hernandez, Guilherme Goergen and Jonatan Tatsch
##' @importFrom dplyr %>% 
##' @importFrom REddyProc fCalcPotRadiation
    Rg_Rpot <- function(date,lon=-53.76,lat=-29.72,timezone=-4){
        
        if(lon==-53.76 & lat==-29.72) 
            warning("Latitude e Longitude de Santa Maria",immediate. = TRUE)
        
        fCalcPotRadiation(DoY.V.n = format(date,"%j") %>% as.numeric,
                          Hour.V.n = format(date,"%H") %>% as.numeric,
                          Lat_deg.n = lat,
                          Long_deg.n = lon,
                          TimeZone_h.n = timezone,
                          useSolartime.b = TRUE
        )
    }

##' Downward Longwave Radiation from a emissivity and a cloud cover schemes
##' 
##' @param data Data frame with column of date (date), temperature (Ta), 
##' partial vapor pressure (es), potencial radiation (Rpot), relative humidity (rh), 
##' atennuation index (K)
##' @param E_fun Emissivity scheme
##' @param C_fun Cloud cover scheme
##' @return Vector with Downward Longwave Radiation time series
    roli_i <- function(data,E_fun = "FHY",C_fun = "CQB"){
        emis_ <- do.call(E_fun,list(data=data,func = C_fun)) 
        roli_est <-  with(data, emis_*sigma*Ta^4)   
        out_rol <- data.frame(ROL = roli_est) 
        names(out_rol) <- paste(E_fun,C_fun,sep = "_")
        out_rol
    }

##' Function for day and nigth identification
##' @param date A vector with dates
##' @param lon Longitude from local analisys 
##' @param lat Latitude from local analisys
##' @param timezone Local time diference with GMT (+1 for fluxes measurement)
##' @importFrom dplyr %>% 
##' @importFrom REddyProc fCalcPotRadiation
##' @return Vector with "day"/"night" string
    to.daylight <- 
        function(date,lon = -47.63,lat = -21.61, timezone = -4){

            
            if(lon == -47.63 & lat == -21.61 & timezone == -3) 
                message("Warning: This data are only valid for PdG site. RHV")
            
            date1 = format(date, "%j") %>% as.numeric()
            date2 = format(date, "%H") %>% as.numeric()
            
            ifelse(fCalcPotRadiation(DoY.V.n = date1,
                                     Hour.V.n = date2,  
                                     Lat_deg.n = lat,
                                     Long_deg.n = lon,
                                     TimeZone_h.n = timezone,    
                                     useSolartime.b = TRUE ) > 0 , 
                   "day","night")
            
        }
    
##' Function for evaluation of all parameterizations 
##' @param data_ A data frame with all atmospherics variables
##' @param lon Longitude from local analisys 
##' @param lat Latitude from local analisys
##' @param timezone Local time diference with GMT (-1 for fluxes measurement)
##' @importFrom dplyr %>% filter mutate arrange_ bind_cols
##' @importFrom openair cutData 
##' @importFrom magrittr set_names %<>%
##' @importFrom hydroGOF gof
##' @importFrom tidyr gather separate spread_
##' @return List with statistical error information between predicted and observed Li 
    eval.params <- function(data_,lon=-53.76,lat=-29.72,timezone=-4){
        
        data_ <-  
            cutData(x = data_,type = "season",hemisphere = "southern") %>% 
            mutate(daytime = to.daylight(date,lon=lon,lat=lat,timezone=timezone))
        
        lapply(with(data_, unique(season)) %>% as.vector, function(j){ 
   
            season.data <-  data_ %>%
                filter(season == j)
            
            lapply(with(season.data,unique(daytime)), function(k){
                
                in.data <- season.data %>% filter_("daytime" == k) 
                
                estats.roli <- 
                    lapply(with(in.data,unique(params)), function(i){ # i = "FBM_CQB"
                        
                        tdy.roli.filt <- 
                            in.data %>%
                            filter_("params" == i)
                        
                        gof.data <- 
                            gof(sim = with(tdy.roli.filt,value),
                                obs = with(tdy.roli.filt,Li),
                                na.rm = TRUE) %>% 
                            as.data.frame() %>%
                            set_names(i) 
                        
                    }) %>% bind_cols() 
                
                estats.roli %<>% mutate(stats = rownames(gof(1:10,10:1)) )
                
                estats.roli %<>% 
                    gather(params,value,-stats) %>%
                    arrange_("stats")
                
                estats.roli.arrange <- 
                    estats.roli %>% 
                    separate(params,sep = "_",into = c("emis","aten")) %>% 
                    spread_("emis","value")
                
                estats.roli.arrange
                
            }) %>% set_names(with(data_,unique(daytime)))
        }) %>% set_names(with(data_,unique(season)) %>% as.vector)
    }
    

   select_stats <- function(roli_list,idx = "RMSE"){
        
        stats.rmse <- NULL
        
        for(i in 1:4){
            for(j in 1:2){
                season <- names(roli_list)[i]
                day.tim <- names(roli_list[[i]])[j]
                stats.rmse[[season]][[day.tim]] <- 
                    roli_list[[season]][[day.tim]] %>% filter_("stats" == idx)
                
            }
        }
        stats.rmse
    }
    
