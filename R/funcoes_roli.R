




##
#/////////////////////////////////////////////////////////////////////////////////////////
#    PARAMETRIZAÇÔES DE  COBERTURA DE NUVENS

##' Amount of cloud estimatives functions.
##' @param data a data frame with all atmospherics variables
##' @return a vector with cloud amount estimatives
##' @export
##' @references Black JN (1956) The distribution of solar radiation over the Earth's surface.
##' Arch Meteor Geophy B 7:165–189
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
##' @references Kasten, F.  Czeplak, G. (1980) Solar and terrestrial radiation
##' dependent on the amount and type of cloud. Sol Energy 24:177–189
CKC <- function(data){
    with(data,maxlim( (4/3*(1-K))^(1/3.4) ))   
} # Kasten & Czeplack (1980)
    
##' Amount of cloud estimatives functions.
##' 
##' @param data a data frame with all atmospherics variables
##' 
##' @return a vector with cloud amount estimatives 
##' @export
##' @references Campbell, G.S. (1985) Soil physics with BASIC. Elsevier, Amsterdam
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

##' Amount of cloud estimatives functions.
##' @param data a data frame with all atmospherics variables
##' @return a vector with cloud amount estimatives 
##' @export
##' @references Stockli R (2007) LBA-MIP driver data gap filling algorithms.
##' Unpublished. http://www.climatemodeling.org/lba-mip/LBAmipDriverDataFillingMethods.pdf
##' Accessed 7 May 2011
CLM <- function(data){ 
    with(data,maxlim( 1-K ))
} # LBA_MIP  (2006)

##' Amount of cloud estimatives functions.
##' @param data a data frame with all atmospherics variables
##' @return a vector with cloud amount estimatives 
##' @export
##' @references Flerchinger, G. N. (2009) Comparison of algoritmhs for incoming 
##' atmospheric long-wave radiation, Water Resources Res., 45, W0342.  
CFG <- function(data){ 
    with(data,maxlim( ((1-K/max(K,na.rm = TRUE))/0.75)^(1./3.4) ))
} # 



##
#/////////////////////////////////////////////////////////////////////////////////////////
                                # PARAMETRIZAÇÔES DE 
                                #    EMISSIVIDADE 
Ofun <- function(Li,Ta){
    sigma <- 5.67051*10^(-8)
    maxlim(Li/(sigma*Ta^4))
}


##' Emissivity from atmosphere
##' @param data Data frame with all atmospherics variables
##' @param func Function for amount of cloud 
##' @param coef1,coef2,coef3,coef4,coef5 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @param forced Forced adjust, default TRUE
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @importFrom  minpack.lm nls.lm nls.lm.control
##' @export
##' @references Angstrom, A. (1915) A study of the radiation of the atmosphere.
##' Smithsonian Miscellaneous Collections 65(3)
EAN <- function(data,
                func = "-", 
                coef1 = 0.83, 
                coef2 = 0.18, 
                coef3 = 0.067, 
                coef4 = 0.22,
                coef5 = 1.0,
                adjust = FALSE,
                forced = TRUE){ 
    
    if(func != "-"){
        data$cp <- do.call(func , args = list(data = data)) 
        start.coefs <- c(c1 = coef1, c2 = coef2, c3 = coef3,ct = coef4, ce = coef5)
    } else {
        data$cp <- 0
        start.coefs <- c(c1 = coef1, c2 = coef2, c3 = coef3)
    }
    
    
    if(adjust ){
        
        if(func != "-"){
            
            nls.out <- try(nls( Ofun(Li,Ta) ~ ean(es,Ta,rh,cp,c1,c2,c3,ct,ce),
                                data = data,  start = start.coefs ),silent = TRUE)
            
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - ean(es,Ta,rh,cp,c1,c2,c3,ct,ce))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
            
        } else {
            
            nls.out <- try(nls( Ofun(Li,Ta) ~ ean(es,Ta,rh,cp,c1,c2,c3),
                                data = data, 
                                start = start.coefs ),silent = TRUE)
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - ean(es,Ta,rh,cp,c1,c2,c3))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 1000))
                }
            }
        }   
        
        if(class(nls.out) != "try-error"){
            
            new.coefs <- coef(nls.out) %>%
                setNames(paste0("coef",1:length(.)))
            
            new.emiss <-
                with(data = data, 
                     run_fun(E_fun = EAN,
                             data = data, 
                             func = func,
                             new.coefs = new.coefs) )
            
            return(list(emiss = new.emiss, coefs = new.coefs))
            
        } else {
            
            return(list(emiss = NA, coefs = NA))
            
        }
    }  else {
        
        return(with(data, ean(es,Ta,rh,cp, c1=coef1,c2 = coef2, c3= coef3,ct= coef4,ce = coef5) ))
        
    }
    
}   ## Angstrom (1915)

ean <- function(es, Ta,rh,cp,
                c1 = 0.83, 
                c2 = 0.18, 
                c3 = 0.067, 
                ct = 0.22,
                ce = 1.0){
    
    maxlim( (c1 - c2*(10^(-c3*es)))*(1+ ct * cp^ce) )
    
}

    
##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2,coef3,coef4 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @param forced Forced adjust, default TRUE
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @importFrom  minpack.lm nls.lm nls.lm.control
##' @export
##' @references Brunt D (1932) Notes on radiation in the atmosphere I. Q J Roy
##' Meteor Soc 58:389–420
EBR <- function(data,
                func = "-",
                coef1 = 0.51, 
                coef2 = 0.066,
                coef3 = 0.22,
                coef4 = 1.0,
                adjust = FALSE,
                forced = TRUE){
    
    
    if(func != "-"){
        data$cp <- do.call(func , args = list(data = data)) 
        start.coefs <- c(c1 = coef1, c2 = coef2,ct = coef3, ce = coef4)
    } else { 
        data$cp <- 0
        start.coefs <- c(c1 = coef1, c2 = coef2) #*
    }
    
    if(adjust){
        
        if(func != "-"){
            
            nls.out <- try(nls( Ofun(Li,Ta) ~ ebr(es,Ta,rh,cp,c1,c2,ct,ce),
                                data = data,  start = start.coefs ),silent = TRUE)
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - ebr(es,Ta,rh,cp,c1,c2,ct,ce))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
        } else {
            nls.out <- try(nls( Ofun(Li,Ta) ~ ebr(es,Ta,rh,cp,c1,c2),
                                data = data, 
                                start = start.coefs ),silent = TRUE)
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - ebr(es,Ta,rh,cp,c1,c2))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
        }   
        
        if(class(nls.out) != "try-error"){
            new.coefs <- coef(nls.out) %>%
                setNames(paste0("coef",1:length(.)))
            
            
            new.emiss <-
                with(data = data, 
                     run_fun(E_fun = EBR, ####
                             data = data, 
                             func = func,
                             new.coefs = new.coefs) )
            
            return(list(emiss = new.emiss, coefs = new.coefs)) 
            
        } else { 
            return(list(emiss = NA, coefs = NA))
        }
        
    } else {
        
        return(with(data,ebr(es,Ta,rh,cp,c1 = coef1,c2 = coef2, ct = coef3, ce = coef4)) )
        
    }
    
}    ## Brunt (1932)

ebr <- function(es,Ta,rh,cp,
                c1 = 0.51, 
                c2 = 0.066,
                ct = 0.22,
                ce = 1.0){
    
    maxlim( (c1 + c2*sqrt(es))*(1+ct*cp^ce) )
    
}

##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud
##' @param coef1,coef2,coef3,coef4 Scheme coeficients  
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @param forced Forced adjust, default TRUE
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @importFrom  minpack.lm nls.lm nls.lm.control
##' @export
##' @references Brutsaert W (1975) On a derivable formula for long-wave radiation
##' from clear skies. Water Resour Res 11:742–744
EBT <- function(data,
                func = "-",
                coef1 = 1.24, 
                coef2 = 1/7,
                coef3 = 0.22,
                coef4 = 1.0,
                adjust = FALSE,
                forced = TRUE){ 
    
    
    if(func != "-"){
        data$cp <- do.call(func , args = list(data = data)) 
        start.coefs <- c(c1 = coef1,c2=coef2,ct = coef3, ce = coef4)
    } else { 
        data$cp <- 0
        start.coefs <- c(c1 = coef1,c2=coef2)
    }
    
    if(adjust){
        
        if(func != "-"){
            
            nls.out <- try(nls( Ofun(Li,Ta) ~ ebt(es,Ta,rh,cp,c1,c2,ct,ce),
                                data = data,  start = start.coefs ),silent = TRUE)
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - ebt(es,Ta,rh,cp,c1,c2,ct,ce))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
        } else {
            nls.out <- try(nls( Ofun(Li,Ta) ~ ebt(es,Ta,rh,cp,c1,c2),
                                data = data, 
                                start = start.coefs ),silent = TRUE)
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - ebt(es,Ta,rh,cp,c1,c2))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
            
        }   
        
        if(class(nls.out) != "try-error"){
            new.coefs <- coef(nls.out) %>%
                setNames(paste0("coef",1:length(.)))
            
            new.emiss <-
                with(data = data, 
                     run_fun(E_fun = EBT, ####
                             data = data, 
                             func = func,
                             new.coefs = new.coefs) )
            
            return(list(emiss = new.emiss, coefs = new.coefs)) 
        } else {
            return(list(emiss = NA, coefs = NA)) 
        }
    } else {
        return(with(data,ebt(es,Ta,rh,cp,c1 = coef1,c2= coef2,ct=coef3,ce=coef4))    )
    }
    
}   ## Brutsaert (1934)

ebt <- function(es,Ta,rh,cp,
                c1 = 1.24, 
                c2 = 1/7,
                ct = 0.22,
                ce = 1.0){
    
    maxlim(c1*(es/Ta)^(c2)*(1.0+ ct*cp^ce) )
    
}


##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2,coef3,coef4,coef5 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @param forced Forced adjust, default TRUE
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @importFrom  minpack.lm nls.lm nls.lm.control
##' @export
##' @references Dilley , A. C. (1998) Estimating downward clear-sky long-wave 
##' irradiance at the surface from screen temperature and precpitable water. 
##' Q. J. R. Meteorol. Soc., 96, 313-319.
EDO <- function(data,
                func = "-",
                coef1 = 59.38, 
                coef2 = 113.7, 
                coef3 = 96.96,
                coef4 = 0.22, 
                coef5 = 1.0, 
                adjust = FALSE,
                forced = TRUE){ 
    
    if(func != "-"){
        data$cp <- do.call(func , args = list(data = data)) 
        start.coefs <- c(c1 = coef1, c2 = coef2, c3 = coef3, ct = coef4, ce = coef5)
    } else {
        data$cp <- 0
        start.coefs <- c(c1 = coef1, c2 = coef2, c3 = coef3)
    }
    
    if(adjust){
        # backup >>
        # tmp.nls <- nls( Li ~ ((1.0+coef4*cp^coef5)*(  coef1+  
        #                       coef2*(Ta/273.15)^6 + coef3 * sqrt((465/25)*(es/Ta))))  ,
        #                 data = data , 
        #                 start = start.coefs)
        
        if(func != "-"){
            
            nls.out <- try(nls( Li ~ edo(es,Ta,rh,cp,c1,c2,c3,ct,ce),
                                data = data,  start = start.coefs ),silent = TRUE)
            if(forced){    
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Li - edo(es,Ta,rh,cp,c1,c2,c3,ct,ce))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
        } else {
            nls.out <- try(nls(  Li ~ edo(es,Ta,rh,cp,c1,c2,c3),
                                 data = data, 
                                 start = start.coefs ),silent = TRUE)
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Li - edo(es,Ta,rh,cp,c1,c2,c3))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
            
        }   
        
        if(class(nls.out) != "try-error"){
            new.coefs <- coef(nls.out) %>%
                setNames(paste0("coef",1:length(.)))
            
            new.emiss <-
                with(data = data, 
                     run_fun(E_fun = EDO, ####
                             data = data, 
                             func = func,
                             new.coefs = new.coefs) )
            
            return(list(emiss = new.emiss, coefs = new.coefs)) 
        } else {
            return(list(emiss = NA, coefs = NA)) 
        }
    } else {
        
        sigma <- 5.67051*10^(-8)
        return( with(data, edo(es,Ta,rh,cp,c1 = coef1,c2 = coef2,c3 = coef3,ct = coef4,ce = coef5)/(sigma*Ta^4)))    
        
    }
    
}     ## Dilley (1963)

edo <- function(es, Ta,rh,cp,
                c1 = 0.83, 
                c2 = 0.18, 
                c3 = 0.067, 
                ct = 0.22,
                ce = 1.0){
    
    # sigma <- 5.67051*10^(-8)
    
    (c1 + c2 * (Ta/273.15)^6 + c3*sqrt((465/25)*(es/Ta))) * (1.0+ct*cp^ce)
    
}


##' Effective emissivity from atmosphere with cloud atenuation
##' 
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2,coef3,coef4 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @param forced Forced adjust, default TRUE
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @importFrom  minpack.lm nls.lm nls.lm.control
##' @export
##' @references Gabathuler M, Marty CA, Hanselmann KW (2001) Parameterization
##' of incoming longwave radiation in high-mountain environments.
##' Phys Geogr 22,99–114
EGB <- function(data,
                func = "-",
                coef1 = 0.84, 
                coef2 = 21.,
                coef3 = 0.22, 
                coef4 = 1.,
                adjust = FALSE,
                forced = TRUE){ 
    
    
    if(func != "-"){
        data$cp <- do.call(func , args = list(data = data)) 
        start.coefs <- list(c1 = coef1, c2 = coef2, ct = coef3, ce = coef4)
    } else { 
        data$cp <- 0
        start.coefs <- list(c1 = coef1, c2 = coef2)
    }
    
    if(adjust){
        
        if(func != "-"){
            
            nls.out <- try(nls( Ofun(Li,Ta) ~ egb(es,Ta,rh,cp,K,c1,c2,ct,ce),
                                data = data,  start = start.coefs ),silent = TRUE)
            
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - egb(es,Ta,rh,cp,K,c1,c2,ct,ce))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
        } else {
            
            nls.out <- try(nls( Ofun(Li,Ta) ~ egb(es,Ta,rh,cp,K,c1,c2),
                                data = data, 
                                start = start.coefs ),silent = TRUE)
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - egb(es,Ta,rh,cp,K,c1,c2))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
        }   
        
        if(class(nls.out) != 'try-error'){
            new.coefs <- coef(nls.out) %>%
                setNames(paste0("coef",1:length(.)))
            
            new.emiss <-
                with(data = data, 
                     run_fun(E_fun = EGB, ####
                             data = data, 
                             func = func,
                             new.coefs = new.coefs) )
            
            return(list(emiss = new.emiss, coefs = new.coefs))     
        } else {
            return(list(emiss = NA, coefs = NA))     
        }
    } else {
        
        return(with(data,egb(es,Ta,rh,cp,K,c1 = coef1,c2=coef2,ct=coef3,ce=coef4)))
        
    }
    
} ## Gabathuler (2001)

egb <- function(es,Ta,rh,cp,K,
                c1 = 0.51, 
                c2 = 0.066,
                ct = 0.22,
                ce = 1.0){
    sigma <- 5.67051*10^(-8)
    a <- (c1*(rh-68))/(sigma*Ta^4)
    b <- (1.- c2*K/Ta)^4
    return( maxlim( (a+b)*(1.+ct*cp^ce)) )
}

##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud
##' @param coef1,coef2,coef3,coef4,coef5 Scheme coeficients  
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @param forced Forced adjust, default TRUE
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @importFrom  minpack.lm nls.lm nls.lm.control
##' @export
##' @references Garratt, J.A. (1992), Extreme maximun land surface temperatures,
##' J. Appl. Meteorol., 31, 1096-1105.
EGR <- function(data,func = "-",
                coef1 = 0.79,
                coef2 = 0.17,
                coef3 = 0.096,
                coef4 = 0.22,
                coef5 = 1.0,
                adjust = FALSE,
                forced = TRUE){ 
    
    
    if(func != "-"){
        data$cp <- do.call(func , args = list(data = data)) 
        start.coefs <- c(c1 = coef1, c2 = coef2, c3 = coef3,ct = coef4, ce = coef5) 
    } else { 
        data$cp <- 0
        start.coefs <- c(c1 = coef1, c2 = coef2, c3 = coef3) 
    }
    
    if(adjust){
        
        if(func != "-"){
            
            nls.out <- try(nls( Ofun(Li,Ta) ~ egr(es,Ta,rh,cp,c1,c2,c3,ct,ce),
                                data = data,  start = start.coefs ),silent = TRUE)
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - egr(es,Ta,rh,cp,c1,c2,c3,ct,ce))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
        } else {
            nls.out <- try(nls( Ofun(Li,Ta) ~ egr(es,Ta,rh,cp,c1,c2,c3),
                                data = data, 
                                start = start.coefs ),silent = TRUE)
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - egr(es,Ta,rh,cp,c1,c2,c3))
                        out[!is.na(out)]  
                    }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
        }   
        
        if(class(nls.out) != "try-error"){
            new.coefs <- coef(nls.out) %>%
                setNames(paste0("coef",1:length(.)))
            
            new.emiss <-
                with(data = data, 
                     run_fun(E_fun = EGR, ####
                             data = data, 
                             func = func,
                             new.coefs = new.coefs) )
            
            return(list(emiss = new.emiss, coefs = new.coefs))     
        } else {
            return(list(emiss = NA, coefs = NA))     
        }
    } else {
        return(with(data,egr(es,Ta,rh,cp,c1=coef1,c2=coef2,c3=coef3,ct=coef4,ce=coef5)))
    }
    
}   ## Garratt (1992)

egr <- function(es,Ta,rh,cp, 
                c1 = 0.79,
                c2 = 0.17,
                c3 = 0.096,
                ct = 0.22,
                ce = 1.0){
    
    maxlim( (c1 - c2 * exp(-c3*es) )*(1.0 + ct*cp^ce) )
    
}

##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2,coef3,coef4 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @param forced Forced adjust, default TRUE
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @importFrom  minpack.lm nls.lm nls.lm.control
##' @export
##' @references Idso SB, Jackson RD (1969) Thermal radiation from the atmosphere.
##' J Geophys Res 74:5397–5403
EIJ <- function(data,func = "-", 
                coef1 = 0.261, 
                coef2 = 0.0007,
                coef3 = 0.22, 
                coef4 = 1.0,
                adjust = FALSE,
                forced = TRUE){
    
    
    if(func != "-"){
        data$cp <- do.call(func , args = list(data = data)) 
        start.coefs <- c(c1 = coef1, c2 = coef2,  ct=coef3,ce=coef4)
    } else { 
        data$cp <- 0
        start.coefs <- c(c1 = coef1, c2 = coef2)
    }
    
    if(adjust){
        
        if(func != "-"){
            
            nls.out <- try(nls( Ofun(Li,Ta) ~ eij(es,Ta,rh,cp,c1,c2,ct,ce),
                                data = data,  start = start.coefs ),silent = TRUE)
            
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - eij(es,Ta,rh,cp,c1,c2,ct,ce))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
        } else {
            nls.out <- try(nls( Ofun(Li,Ta) ~ eij(es,Ta,rh,cp,c1,c2),
                                data = data, 
                                start = start.coefs ),silent = TRUE)
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - eij(es,Ta,rh,cp,c1,c2))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
        }   
        
        if(class(nls.out) != "try-error"){
            new.coefs <- coef(nls.out) %>%
                setNames(paste0("coef",1:length(.)))
            
            new.emiss <-
                with(data = data, 
                     run_fun(E_fun = EIJ, ####
                             data = data, 
                             func = func,
                             new.coefs = new.coefs) )
            
            return(list(emiss = new.emiss, coefs = new.coefs)) 
        } else {
            return(list(emiss = NA, coefs = NA)) 
        }
    } else {
        return(with(data,eij(es,Ta,rh,cp,c1=coef1,c2=coef2,ct=coef3,ce=coef4)) )
    }
    
} ## Idso & Jackson (1969)

eij <- function(es,Ta,rh,cp,
                c1 = 0.261, 
                c2 = 0.0007,
                ct = 0.22,
                ce = 1.0){
    
    maxlim((1.-c1*exp(c2*(273.15-Ta)^2))*(1.+ct*cp^ce))
    
}


##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2,coef3,coef4 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @param forced Forced adjust, default TRUE
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @importFrom  minpack.lm nls.lm nls.lm.control
##' @export
##' @references Idso SB (1981) A set of equations for full spectrum and 8- to 14-um
##'  and 10.5- to 12.5-um thermal radiation from cloudless skies.
##'  Water Resour Res 17:295–304
EID <- function(data,func = "-",
                coef1 = 0.7,
                coef2 = 0.0000595,
                coef3 = 0.22, 
                coef4 = 1.0,
                adjust = FALSE,
                forced = TRUE){ 
    
    
    
    if(func != "-"){
        data$cp <- do.call(func , args = list(data = data)) 
        start.coefs <- c(c1 = coef1, c2 = coef2,ct = coef3, ce = coef4)
    } else { 
        data$cp <- 0
        start.coefs <- c(c1 = coef1, c2 = coef2)
    }
    
    if(adjust){
        
        if(func != "-"){
            
            nls.out <- try(nls( Ofun(Li,Ta) ~ eid(es,Ta,rh,cp,c1,c2,ct,ce),
                                data = data,  start = start.coefs ),silent = TRUE)
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - eid(es,Ta,rh,cp,c1,c2,ct,ce))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
        } else {
            nls.out <- try(nls( Ofun(Li,Ta) ~ eid(es,Ta,rh,cp,c1,c2),
                                data = data, 
                                start = start.coefs ),silent = TRUE)
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - eid(es,Ta,rh,cp,c1,c2))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
        }   
        
        if(class(nls.out) != "try-error"){
            new.coefs <- coef(nls.out) %>%
                setNames(paste0("coef",1:length(.)))
            
            new.emiss <-
                with(data = data, 
                     run_fun(E_fun = EID, ####
                             data = data, 
                             func = func,
                             new.coefs = new.coefs) )
            
            return(list(emiss = new.emiss, coefs = new.coefs)) 
        } else {
            return(list(emiss = NA, coefs = NA)) 
        }
    } else {
        
        return( with(data,eid(es,Ta,rh,cp,c1=coef1,c2=coef2,ct=coef3,ce=coef4))  )
        
    }
    
    
}   ## Idso (1981)

eid <- function(es,Ta,rh,cp,
                c1 = 0.7,
                c2 = 0.0000595,
                ct = 0.22,
                ce = 1.0){
    
    maxlim( (c1 + c2*es*exp(1500/Ta))*(1.+ct*cp^ce) )
    
}


##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2,coef3,coef4,coef5 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @param forced Forced adjust, default TRUE
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @importFrom  minpack.lm nls.lm nls.lm.control
##' @export
##' @references Konzelmann T, van de Wal RSW, Greuell W, Bintanja R, Henneken
##' EAC, Abe-Ouchi A (1994) Parameterization of global and
##' longwave incoming radiation for the Greenland Ice Sheet. Glob
##' Planet Chang 9:143–164
EKZ <- function(data,func = "-",
                coef1 = 0.23, 
                coef2 = 0.484,
                coef3 = 1/8,
                coef4 = 0.22,
                coef5 = 1.,
                adjust = FALSE,
                forced = TRUE) {
    
    
    if(func != "-"){
        data$cp <- do.call(func , args = list(data = data)) 
        start.coefs <- c(c1 = coef1, c2 = coef2, c3 = coef3,  ct = coef4, ce = coef5)
    } else { 
        data$cp <- 0
        start.coefs <- c(c1 = coef1, c2 = coef2, c3 = coef3)
    }
    
    
    if(adjust){
        
        if(func != "-"){
            
            nls.out <- try(nls( Ofun(Li,Ta) ~ ekz(es,Ta,rh,cp,c1,c2,c3,ct,ce),
                                data = data,  start = start.coefs ),silent = TRUE)
            if(forced){
            if(class(nls.out) == "try-error"){
                
                resEOfun <- function(par , idata = data) {
                    idata <- cbind(idata,data.frame(t(par)))
                    out <- with(idata, Ofun(Li,Ta) - ekz(es,Ta,rh,cp,c1,c2,c3,ct,ce))
                    out[!is.na(out)]  }
                
                nls.out <- 
                    nls.lm(fn = resEOfun,
                           par = start.coefs,
                           idata = data, 
                           control = nls.lm.control(nprint = 1,maxiter = 75))
            }
            }
        } else {
            nls.out <- try(nls( Ofun(Li,Ta) ~ ekz(es,Ta,rh,cp,c1,c2,c3),
                                data = data, 
                                start = start.coefs ),silent = TRUE)
            if(forced){
            if(class(nls.out) == "try-error"){
                
                resEOfun <- function(par , idata = data) {
                    idata <- cbind(idata,data.frame(t(par)))
                    out <- with(idata, Ofun(Li,Ta) - ekz(es,Ta,rh,cp,c1,c2,c3))
                    out[!is.na(out)]  }
                
                nls.out <- 
                    nls.lm(fn = resEOfun,
                           par = start.coefs,
                           idata = data, 
                           control = nls.lm.control(nprint = 1,maxiter = 75))
            }
                }
        }   
        
        if(class(nls.out) != "try-error"){
        new.coefs <- coef(nls.out) %>%
            setNames(paste0("coef",1:length(.)))
        
        new.emiss <-
            with(data = data, 
                 run_fun(E_fun = EKZ, ####
                         data = data, 
                         func = func,
                         new.coefs = new.coefs) )
        
        return(list(emiss = new.emiss, coefs = new.coefs))   
        } else {
            return(list(emiss = NA, coefs = NA))   
        }
    } else {
        
        return( with(data,ekz(es,Ta,rh,cp,c1=coef1,c2=coef2,c3=coef3,ct=coef4,ce=coef5))  )
        
    }
    
    
}   ## Konzelmann (1994)

ekz <- function(es,Ta,rh,cp,
                c1 = 0.23, 
                c2 = 0.484,
                c3 = 1/8,
                ct = 0.22,
                ce = 1.0){
    
  maxlim( (c1 + c2*(es/Ta)^(c3))*(1.+ct*cp^ce) )
    
}



##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2,coef3,coef4,coef5 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @param forced Forced adjust, default TRUE
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @importFrom  minpack.lm nls.lm nls.lm.control
##' @export
##' @references Niemala, S. (2001) Comparison of surface radiative flux 
##' parameterizations part I: Longwave radiation. Atmos. Res., 58, 1-18.
ENM <- function(data,func = "-",
                coef1 = 0.72,
                coef2 = 0.09,
                coef3 = 0.76,
                coef4 = 0.22,
                coef5 = 1.0,
                adjust = FALSE,
                forced = TRUE){ 
    
    if(func != "-"){
        data$cp <- do.call(func , args = list(data = data)) 
        start.coefs <- c(c1 = coef1, c2 = coef2,c3 = coef3,ct = coef4, ce =coef5) 
    } else { 
        data$cp <- 0
        start.coefs <- c(c1 = coef1, c2 = coef2,c3 = coef3)
    }
    
    if(adjust){
        
        if(func != "-"){
            
            nls.out <- try(nls( Ofun(Li,Ta) ~ enm(es,Ta,rh,cp,c1,c2,c3,ct,ce),
                                data = data,  start = start.coefs ),silent = TRUE)
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - enm(es,Ta,rh,cp,c1,c2,c3,ct,ce))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
        } else {
            nls.out <- try(nls( Ofun(Li,Ta) ~ enm(es,Ta,rh,cp,c1,c2,c3),
                                data = data, 
                                start = start.coefs ),silent = TRUE)
            
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - enm(es,Ta,rh,cp,c1,c2,c3))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
        }   
        
        if(class(nls.out) != "try-error"){
            new.coefs <- coef(nls.out) %>%
                setNames(paste0("coef",1:length(.)))
            
            new.emiss <-
                with(data = data, 
                     run_fun(E_fun = ENM, ####
                             data = data, 
                             func = func,
                             new.coefs = new.coefs) )
            
            return(list(emiss = new.emiss, coefs = new.coefs))
            
        } else {
            return(list(emiss = NA, coefs = NA)) 
        }
        
    } else {
        
        return(with(data,enm(es,Ta,rh,cp,c1=coef1,c2=coef2,c3=coef3,ct=coef4,ce=coef5)))
        
    }
    
}   ## Niemala (2001)

enm <- function(es,Ta,rh,cp,
                c1 = 0.72,
                c2 = 0.09,
                c3 = 0.76,
                ct = 0.22,
                ce = 1.0){
    
    maxlim((c1 + sign(es-20.)*ifelse(sign(es-20.) > 0,c2,c3)*(es-20.))*(1.+ct*cp^ce))
    
}

##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2,coef3,coef4,coef5 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @param forced Forced adjust, default TRUE
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @importFrom  minpack.lm nls.lm nls.lm.control
##' @export
##' @references Prata AJ (1996) A new long-wave formula for estimating downward clearsky
##' radiation at the surface. Q J Roy Meteor Soc 122:1127–1151
EPR <- function(data,func = "-",
                coef1 = 1, 
                coef2 = 1.2, 
                coef3 = 3,
                coef4 = 0.22, 
                coef5 =1.,
                adjust = FALSE,
                forced = TRUE) {
    
    if(func != "-"){
        data$cp <- do.call(func , args = list(data = data)) 
        start.coefs <- c(c1 = coef1,  c3 = coef3,ct = coef4, ce = coef5)
    } else { 
        data$cp <- 0
        start.coefs <- c(c1 = coef1,  c3 = coef3)
    }
    
    if(adjust){
        
        if(func != "-"){
            
            nls.out <- try(nls( Ofun(Li,Ta) ~ epr(es,Ta,rh,cp,c1,c3,ct,ce),
                                data = data,  start = start.coefs ),silent = TRUE)
            
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - epr(es,Ta,rh,cp,c1,c3,ct,ce))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
        } else {
            nls.out <- try(nls( Ofun(Li,Ta) ~ epr(es,Ta,rh,cp,c1,c3),
                                data = data, 
                                start = start.coefs ),silent = TRUE)
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - epr(es,Ta,rh,cp,c1,c3))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
        }   
        
        if(class(nls.out) != "try-error"){
            new.coefs <- coef(nls.out) %>%
                setNames(paste0("coef",1:length(.)))
            
            new.emiss <-
                with(data = data, 
                     run_fun(E_fun = EPR, ####
                             data = data, 
                             func = func,
                             new.coefs = new.coefs) )
            
            return(list(emiss = new.emiss, coefs = new.coefs))     
        } else {
            return(list(emiss = NA, coefs = NA))     
        }
    } else {
        
        return( with(data, epr(es,Ta,rh,cp,c1=coef1,c2=coef2,c3=coef3,ct=coef4,ce=coef5)) )
        
    }
    
} ## Prata (1996)


epr <- function(es,Ta,rh,cp,
                c1 = 1, 
                c2 = 1.2,
                c3 = 3,
                ct = 0.22,
                ce = 1.0){
    
    maxlim( (1 - ((c1+46.5*es/Ta) *exp(-sqrt(c2+c3*46.5*es/Ta))))*(1.0+ct*cp^ce))
}


##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2,coef3,coef4 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @param forced Forced adjust, default TRUE
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @importFrom  minpack.lm nls.lm nls.lm.control
##' @export
##' @references Sutterlund, D. R. (1979) An improved equation for estimating longwave 
##' radiation from the atmosphere, Water Res. Res., 15, 1649-1650.
EST <- function(data,
                func = "-",
                coef1 = 1.08,
                coef3 = 0.22, 
                coef4 = 1.0, 
                adjust = FALSE,
                forced = TRUE){ 
    
    
    if(func != "-"){
        data$cp <- do.call(func , args = list(data = data)) 
        start.coefs <- c(c1 = coef1, ct = coef3, ce= coef4)
    } else { 
        data$cp <- 0
        start.coefs <- c(c1 = coef1)
    }
    
    if(adjust){
        
        if(func != "-"){
            
            nls.out <- try(nls( Ofun(Li,Ta) ~ est(es,Ta,rh,cp,c1,ct,ce),
                                data = data,  start = start.coefs ),silent = TRUE)
            
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - est(es,Ta,rh,cp,c1,ct,ce))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
        } else {
            nls.out <- try(nls( Ofun(Li,Ta) ~ est(es,Ta,rh,cp,c1),
                                data = data, 
                                start = start.coefs ),silent = TRUE)
            
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - est(es,Ta,rh,cp,c1))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
            
        }   
        
        if(class(nls.out) == "try-error"){
            new.coefs <- coef(nls.out) %>%
                setNames(paste0("coef",1:length(.)))
            
            new.emiss <-
                with(data = data, 
                     run_fun(E_fun = EST, ####
                             data = data, 
                             func = func,
                             new.coefs = new.coefs) )
            
            return(list(emiss = new.emiss, coefs = new.coefs)) 
        } else {
            return(list(emiss = NA, coefs = NA)) 
        }
    } else {
        
        return( with(data,est(es,Ta,rh,cp,c1 = coef1, ct = coef3, ce = coef4)) )    
        
    }
    
}     ## Swinbank (1963)

est <- function(es,Ta,rh,cp,
                c1 = 1.08, 
                ct = 0.22,
                ce = 1.0){
    
    maxlim(c1*(1.-exp(-es^(Ta/2016)))*(1.0+ct*cp^ce))
    
}


##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2,coef3 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @param forced Forced adjust, default TRUE
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @importFrom  minpack.lm nls.lm nls.lm.control
##' @export
##' @references Swinbank WC (1963) Long-wave radiation from clear skies. Q J Roy
##' Meteor Soc 89:339–348
ESW <- function(data,
                func = "-",
                coef1 = 0.0000092,
                coef2 = 0.22,
                coef3 = 1.0, 
                adjust = FALSE,
                forced = TRUE){ 
    
    
    if(func != "-"){
        data$cp <- do.call(func , args = list(data = data)) 
        start.coefs <- c(c1 = coef1, ct = coef2, ce = coef3)
    } else {
        data$cp <- 0
        start.coefs <- c(c1 = coef1)
    }
    
    if(adjust ){
        
        if(func != "-"){
            
            nls.out <- try(nls( Ofun(Li,Ta) ~ esw(es,Ta,rh,cp,c1,ct,ce),
                                data = data,  start = start.coefs ),silent = TRUE)
            if(forced){
            if(class(nls.out) == "try-error"){
                
                resEOfun <- function(par , idata = data) {
                    idata <- cbind(idata,data.frame(t(par)))
                    out <- with(idata, Ofun(Li,Ta) - esw(es,Ta,rh,cp,c1,ct,ce))
                    out[!is.na(out)]  }
                
                nls.out <- 
                    nls.lm(fn = resEOfun,
                           par = start.coefs,
                           idata = data, 
                           control = nls.lm.control(nprint = 1,maxiter = 75))
            }
            }
        } else {
            nls.out <- try(nls( Ofun(Li,Ta) ~ esw(es,Ta,rh,cp,c1),
                                data = data, 
                                start = start.coefs ),silent = TRUE)
            if(forced){
            if(class(nls.out) == "try-error"){
                
                resEOfun <- function(par , idata = data) {
                    idata <- cbind(idata,data.frame(t(par)))
                    out <- with(idata, Ofun(Li,Ta) - esw(es,Ta,rh,cp,c1))
                    out[!is.na(out)]  }
                
                nls.out <- 
                    nls.lm(fn = resEOfun,
                           par = start.coefs,
                           idata = data, 
                           control = nls.lm.control(nprint = 1,maxiter = 75))
            }}
        }   
        
        if(class(nls.out) != "try-error"){
        new.coefs <- coef(nls.out) %>%
            setNames(paste0("coef",1:length(.)))
        
        new.emiss <-
            with(data = data, 
                 run_fun(E_fun = ESW, ####
                         data = data, 
                         func = func,
                         new.coefs = new.coefs) )
        
        return(list(emiss = new.emiss, coefs = new.coefs))
        } else {
            return(list(emiss = NA, coefs = NA))
        }
    } else {
        
        return( with(data,esw(es,Ta,rh,cp,c1 = coef1,ct = coef2,ce = coef3)) )    
        
    }
    
}     ## Swinbank (1963)
  
esw <- function(es,Ta,rh,cp,
                c1 = 0.0000092,
                ct = 0.22,
                ce = 1.0){
    
    maxlim(c1*Ta^2*(1.0+ct*cp^ce))
    
}


##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2,coef3,coef4,coef5 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @param forced Forced adjust, default TRUE
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @importFrom  minpack.lm nls.lm
##' @export
##' @references Aimi, D. (2017); TODO 
EAI <- function(data,func = "-",
                coef1 = 0.52843, 
                coef2 = -0.00820,
                coef3 = 24242.65010,
                coef4 = 0.22, 
                coef5 = 1.,
                adjust = FALSE,
                forced = TRUE) {
    
    if(func != "-"){
        data$cp <- do.call(func , args = list(data = data)) 
        start.coefs <- list(c1 = coef1,  c2 = coef2, c3 = coef3,
                            ct = coef4, ce = coef5)
    } else {
        data$cp <- 0
        start.coefs <- list(c1 = coef1,  c2 = coef2, c3 = coef3)
    }
    
    if(adjust ){
        if(func != "-"){
            
            nls.out <- try(nls( Ofun(Li,Ta) ~ eai(es,Ta,rh,cp,c1,c2,c3,ct,ce),
                                data = data, start = start.coefs ),silent = TRUE)
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - eai(es,Ta,rh,cp,c1,c2,c3,ct,ce))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
        } else {
            
            nls.out <- try(nls( Ofun(Li,Ta) ~ eai(es,Ta,rh,cp,c1,c2,c3),
                                data = data, start = start.coefs ),
                           silent = TRUE)
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - eai(es,Ta,rh,cp,c1,c2,c3))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }}
        }
        
        if(class(nls.out) != "try-error"){
            new.coefs <- coef(nls.out) %>%
                setNames(paste0("coef",1:length(.)))
            
            new.emiss <-
                with(data = data, 
                     run_fun(E_fun = EAI, ####
                             data = data, 
                             func = func,
                             new.coefs = new.coefs) )
            
            return(list(emiss = new.emiss, coefs = new.coefs))
        } else {
            return(list(emiss = NA, coefs = NA))
        }
        
    } else {
        
        return( with(data,eai(es,Ta,rh,cp,c1= coef1,c2= coef2,c3= coef3,ct= coef4,ce= coef5))) # 
        
    }
    
} ## Aimi (2017)

eai <- function(es,Ta,rh,cp,
                c1 = 0.52843, 
                c2 = -0.00820,
                c3 = 24242.65010,
                # c3 = 1.5, 
                ct = 0.22, 
                ce = 1.){
    sigma <- 5.67051*10^(-8)
    maxlim((c1 + ( c2 * (es - 20.0) ) + c3 / (sigma * Ta^5) )*(1.0+ct*cp^ce))
    # maxlim((c1 + ( c2 * (es - 20) ) + c3 * (273.15/Ta)^5 )*(1.0+ct*cp^ce))
    
}



##' Emissivity from atmosphere
##' @param data a data frame with all atmospherics variables
##' @param func a function for amount of cloud 
##' @param coef1,coef2,coef3,coef4,coef5 Scheme coeficients 
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @param forced Forced adjust, default TRUE
##' @return a vector with emissivity estimatives
##' @import stats
##' @import utils
##' @importFrom  minpack.lm nls.lm
##' @export
##' @references Aimi, D. (2017); TODO 
EAM <- function(data,func = "-",
                coef1 = 0.52843, 
                coef2 = 0.0820,
                coef3 = 0.028,
                coef4 = 0.22, 
                coef5 = 1.,
                adjust = FALSE,
                forced = TRUE) {
    
    if(func != "-"){
        data$cp <- do.call(func , args = list(data = data)) 
        start.coefs <- list(c1 = coef1,  c2 = coef2, c3 = coef3,
                            ct = coef4, ce = coef5)
    } else {
        data$cp <- 0
        start.coefs <- list(c1 = coef1,  c2 = coef2, c3 = coef3)
    }
    
    if(adjust ){
        if(func != "-"){
            
            nls.out <- try(nls( Ofun(Li,Ta) ~ eam(es,Ta,rh,cp,c1,c2,c3,ct,ce),
                                data = data, start = start.coefs ),silent = TRUE)
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - eam(es,Ta,rh,cp,c1,c2,c3,ct,ce))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }
            }
        } else {
            
            nls.out <- try(nls( Ofun(Li,Ta) ~ eam(es,Ta,rh,cp,c1,c2,c3),
                                data = data, start = start.coefs ),
                           silent = TRUE)
            if(forced){
                if(class(nls.out) == "try-error"){
                    
                    resEOfun <- function(par , idata = data) {
                        idata <- cbind(idata,data.frame(t(par)))
                        out <- with(idata, Ofun(Li,Ta) - eam(es,Ta,rh,cp,c1,c2,c3))
                        out[!is.na(out)]  }
                    
                    nls.out <- 
                        nls.lm(fn = resEOfun,
                               par = start.coefs,
                               idata = data, 
                               control = nls.lm.control(nprint = 1,maxiter = 75))
                }}
        }
        
        if(class(nls.out) != "try-error"){
            new.coefs <- coef(nls.out) %>%
                setNames(paste0("coef",1:length(.)))
            
            new.emiss <-
                with(data = data, 
                     run_fun(E_fun = EAM, ####
                             data = data, 
                             func = func,
                             new.coefs = new.coefs) )
            
            return(list(emiss = new.emiss, coefs = new.coefs))
        } else {
            return(list(emiss = NA, coefs = NA))
        }
        
    } else {
        
        return( with(data,eam(es,Ta,rh,cp,c1= coef1,c2= coef2,c3= coef3, ct= coef4,ce= coef5))) # 
        
    }
    
} ## Aimi (2017)



eam <- function(es,Ta,rh,cp,
                c1 = 0.52843, 
                c2 = 0.820,
                c3 = 0.028,
                ct = 0.22, 
                ce = 1.){
    
    maxlim(c1*(Ta - 273.16)/es +  c2 * rh/100  + c3 * (Ta/273.16) )*(1.0+ct*cp^ce)   ####         !!!!
    
}

#/////////////////////////////////////////////////////////////////////////////////////////////////
# Parametrizações de ...
#---- EMISSIVIDADE EFETIVA COM INDICE DE ATENUAÇÃO
# 
# ##' Effective emissivity from atmosphere with cloud atenuation
# ##' @param data a data frame with all atmospherics variables
# ##' @param func a function for amount of cloud 
# ##' @param coef1,coef2,coef3,coef4 Scheme coeficients 
# ##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
# ##' @return a vector with emissivity estimatives
# ##' @import stats
# ##' @import utils
# ##' @export
# ##' @references Lhomme JP, Vacher JJ, Rocheteau A (2007) Estimating downward
# ##' long-wave radiation on the Andean Altiplano. Agr For Meteorol
# ##' 145:139–148
#     ALH <- function(data,func, 
#                     coef1 = 1.18, coef2 = 1.24, coef3 = -.34, coef4 = 1.37,
#                     adjust=FALSE){
#         
#         if(adjust){
#             
#             emiss <- EBT(data = data,func = func, adjust = TRUE)
#             sigma <- 5.67051*10^(-8)
#              
#             
#             data$emiss <- emiss$emiss
#             data$K[is.na(data$emiss)] <- NA
#             
#             suppressWarnings(
#                 tmp.nls <- nls( Li/(sigma*Ta^4) ~ 
#                                     maxlim( (coef1/coef2) * emiss *  (coef3*K+coef4) ) ,
#                             data = data, 
#                             start = list(coef1 = coef1, coef3 = coef3 ) )
#             )
#             new.coefs <- coef(tmp.nls) 
#             new.emiss <- predict(tmp.nls)
#             
#             return(list(emiss = new.emiss, coefs = new.coefs))
#             
#         } else {
#             
#             return(maxlim(coef1/coef2*EBT(data) * with(data,(coef3*K+coef4)) ))
#             
#         }
#                 
#         
#         
#         } ## Lhomme (2007)
#     
#     
    
#/////////////////////////////////////////////////////////////////////////////////////////////////
# Parametrizações de ...
#---- EMISSIVIDADE EFETIVA COM COBERTURA DE NUVENS
#    
#     
# ##' Effective emissivity from atmosphere with cloud atenuation
# ##' 
# ##' @param data a data frame with all atmospherics variables
# ##' @param func a function for amount of cloud 
# ##' @param coef1 Scheme coeficients 
# ##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
# ##' @return a vector with emissivity estimatives
# ##' @import stats
# ##' @import utils
# ##' @export
# ##' @references Angstrom A (1915) A study of the radiation of the atmosphere.
# ##' Smithsonian Miscellaneous Collections 65(3)
#     FAN <- function(data,func,
#                     coef1 = 0.22,
#                     adjust = FALSE){
#         
#       C <- do.call(func , args = list(data = data)) 
#       
#       sigma <- 5.67051*10^(-8)
#       
#       if(adjust){
#           
#           emiss <- EAN(data = data,func=func, adjust = TRUE)
#        
#           data$emiss <- emiss$emiss
#           C[is.na(data$emiss)] <- NA
# 
#           suppressWarnings(
#               tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( emiss * (1+coef1*C) ) ,
#                               data = data,# na.action = "na.exclude",
#                               start = list(coef1 = coef1) )
#           )
#           new.coefs <- coef(tmp.nls) 
#           new.emiss <- predict(tmp.nls)
#           
#           return(list(emiss = new.emiss, coefs = new.coefs))
#           
#       } else {
#           
#           maxlim( EAN(data)*(1+coef1*C)  )    
#           
#       }
#       
#     }   ## Angstrom (1915)
# 
# ##' Effective emissivity from atmosphere with cloud atenuation
# ##' 
# ##' @param data a data frame with all atmospherics variables
# ##' @param func a function for amount of cloud 
# ##' @param coef1 Scheme coeficients 
# ##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
# ##' @return a vector with emissivity estimatives
# ##' @import stats
# ##' @import utils
# ##' @export
# ##' @references Brutsaert W (1975) On a derivable formula for long-wave radiation
# ##' from clear skies. Water Resour Res 11:742–744
#     FBR <- function(data,func,
#                     coef1 = 0.22,
#                     adjust = FALSE){
#         
#         C <- do.call(func , args = list(data = data)) 
#         
#         sigma <- 5.67051*10^(-8)
#         
#         if(adjust){
#             
#             emiss <- EBR(data = data,func=func, adjust = TRUE)
#        
#             data$emiss <- emiss$emiss
#             C[is.na(data$emiss)] <- NA
#             
#             suppressWarnings(
#                 tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( emiss * (1+coef1*C) ) ,
#                                 data = data, #na.action = "na.exclude",
#                                 start = list(coef1 = coef1) )
#             )
#             new.coefs <- coef(tmp.nls) 
#             new.emiss <- predict(tmp.nls)
#             
#             return(list(emiss = new.emiss, coefs = new.coefs))
#             
#         } else {
#             
#             return(   maxlim( EBR(data)*(1+coef1*C)  )    )
#             
#         }
#         
#     }   ## Brutsaert (1982) ou (1975)**
# 
# ##' Effective emissivity from atmosphere with cloud atenuation.
# ##' @param data a data frame with all atmospherics variables
# ##' @param func a function for amount of cloud 
# ##' @param coef1,coef2 Scheme coeficients 
# ##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
# ##' @return a vector with emissivity estimatives
# ##' @import stats
# ##' @import utils
# ##' @export
# ##' @references Friend AD, Stevens AK, Knox RG, Cannell MGR (1997) A processbased,
# ##' terrestrial biosphere model of ecosystem dynamics
# ##' (Hybrid v3.0). Ecol Model 95:249–287
#     FHY <- function(data,func,
#                     coef1 = 0.69, coef2 = 0.979,
#                     adjust = FALSE){
#         
#         C <- do.call(func , args = list(data = data)) 
#         
#         sigma <- 5.67051*10^(-8)
#         
#         if(adjust){
#            
#             
#             suppressWarnings(
#                 tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( coef1 *(1-C^6) + coef2 * C^4 ) ,
#                                 data = data, #na.action = "na.exclude",
#                                 start = list(coef1 = coef1, coef2 = coef2) )
#             )
#             new.coefs <- coef(tmp.nls) 
#             new.emiss <- do.call(FHY,as.list(modifyList(formals(FHY),
#                                                         c(list(data = data, func = func),
#                                                           as.list(new.coefs)))))
#             
#             return(list(emiss = new.emiss, coefs = new.coefs))
#             
#         } else {
#             
#             return(maxlim( coef1 *(1-C^6) + coef2 * C^4))
#         }
#         
#     }     ## HYBRID (2009)
# 
# 
# ##' Effective emissivity from atmosphere with cloud atenuation
# ##' 
# ##' @param data a data frame with all atmospherics variables
# ##' @param func a function for amount of cloud 
# ##' @param coef1 Scheme coeficients 
# ##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
# ##' @return a vector with emissivity estimatives
# ##' @import stats
# ##' @import utils
# ##' @export
# ##' @references Konzelmann T, van de Wal RSW, Greuell W, Bintanja R, Henneken
# ##' EAC, Abe-Ouchi A (1994) Parameterization of global and
# ##' longwave incoming radiation for the Greenland Ice Sheet. Glob
# ##' Planet Chang 9:143–164
#     FKZ <- function(data,func, 
#                     coef1 = 0.963,
#                     adjust = FALSE){
#         C <- do.call(func , args = list(data = data)) 
#         
#         sigma <- 5.67051*10^(-8)
#         
#         if(adjust){
#             
#             emiss <- EKZ(data = data,func = func, adjust = TRUE)
#             
#             data$emiss <- emiss$emiss
#             C[is.na(data$emiss)] <- NA
#             
#             suppressWarnings(
#                 tmp.nls <- nls(Li/(sigma*Ta^4) ~ maxlim( emiss *(1.0-C^3)+ coef1 * C^3 ) ,
#                                 data = data,# na.action = "na.exclude",
#                                 start = list(coef1 = coef1) )
#             )
#             new.coefs <- coef(tmp.nls) 
#             new.emiss <- predict(tmp.nls)
#             
#             return(list(emiss = new.emiss, coefs = new.coefs))
#             
#             
#         } else {
#             return(   EKZ(data)*(1.0-C^3)+ coef1 * C^3    )
#         }
#         
#     }      ## Konzelmann (1994)
#     
# 
# ##' Effective emissivity from atmosphere with cloud atenuation
# ##' 
# ##' @param data a data frame with all atmospherics variables
# ##' @param func a function for amount of cloud 
# ##' @param coef1 Scheme coeficients 
# ##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
# ##' @return a vector with emissivity estimatives
# ##' @import stats
# ##' @import utils
# ##' @export
# ##' @references Idso SB, Jackson RD (1969) Thermal radiation from the atmosphere. 
# ##' J Geophys Res 74:5397–5403
#     FIJ <- function(data,func,
#                     coef1 = 0.22,
#                     adjust = FALSE){
#         
#         C <- do.call(func , args = list(data = data)) 
#         
#         sigma <- 5.67051*10^(-8)
#         
#         if(adjust){
#             
#             emiss <- EIJ(data = data, func = func, adjust = TRUE)
#             
#             data$emiss <- emiss$emiss
#             C[is.na(data$emiss)] <- NA
#             
#             suppressWarnings(
#                 tmp.nls <- nls( Li/(sigma*Ta^4) ~ maxlim( emiss *(1+coef1*C^2) ) ,
#                                 data = data,#na.action = "na.exclude",
#                                 start = list(coef1 = coef1) )
#             )
#             new.coefs <- coef(tmp.nls) 
#             new.emiss <- predict(tmp.nls)
#             
#             return(list(emiss = new.emiss, coefs = new.coefs))
#             
#         } else {
#             
#             return(  EIJ(data)*(1+coef1*C^2) )    
#             
#         }
#         
#     }                                                   
#     
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

    
##' Create from hourly clearness index a mean daily clearness index and 
##' smooth clearness index
##' 
##' @param data_ Dataframe with at least date and Rg column.
##' @param lon longitude from local analisys 
##' @param lat latitude from local analisys
##' @param timezone local time diference with GMT (-1 for fluxes measurement)
##' @param window_size lenght of window to smoothing  
##' @return The same input dataframe with K_hourly, K_mean (daily mean), 
##' K_smooth (smooth K_mean index in specific temporal window) with K index columns
##' @author Roilan Hernandez
##' @importFrom stats setNames
##' @importFrom dplyr %>% mutate select group_by left_join
##' @importFrom plyr . 
##' @importFrom zoo rollmean
##' @export    
correct.clearness.index <- function(data_,
                                    lon, lat, timezone,
                                    window_size = 6){
    if(!("Rpot" %in% names(data_))) {
        data_ <- 
            data_ %>%
            mutate(Rpot =  PotRad(date = date,lon = lon, lat = lat,timezone = timezone))
    }
    
    K_mean <- 
        data_ %>%
        select(date, Rg, Rpot) %>%
        subset(Rpot > 100) %>%
        group_by(day = as.Date(date)) %>%
        summarise(K_mean = mean(Rg,na.rm = TRUE)/mean(Rpot,na.rm = TRUE)) %>%
        select(day,K_mean)
    
    to_K <- 
        data_ %>% 
        mutate(day = as.Date(date)) %>%
        left_join(., K_mean,by = "day") %>% 
        subset(!is.na(K_mean)) %>% 
        mutate(K = rollmean(K_mean,k = window_size,fill = mean(K_mean),
                                 align = "center",na.pad = TRUE)) %>% 
        select(date,K_mean,K)
    
    data_ %>% 
        left_join(., to_K,by = "date") %>%
        mutate(K_hourly = kloudines(dates = date, Rg,lon=lon, lat=lat, timezone=timezone) ) 
    
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
