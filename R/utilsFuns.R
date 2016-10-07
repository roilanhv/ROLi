
#  FUNÇÕES NECESSÁRIAS
#
maxlim <- function(i,max_=1,min_=0){ 
    sapply(i,function(i) min(max(i,min_),max_) ) 
}



fCalcPotRadiation <- function (DoY.V.n, Hour.V.n, Lat_deg.n, Long_deg.n, TimeZone_h.n, 
                               useSolartime.b = TRUE) 
{
    SolElev_rad.V.n <- fCalcSunPosition(DoY.V.n, Hour.V.n, Lat_deg.n, 
                                        Long_deg.n, TimeZone_h.n, useSolartime.b = useSolartime.b)$SolElev
    ExtRadiation.V.n <- fCalcExtRadiation(DoY.V.n)
    PotRadiation.V.n <- ifelse(SolElev_rad.V.n <= 0, 0, ExtRadiation.V.n * 
                                   sin(SolElev_rad.V.n))
    attr(PotRadiation.V.n, "varnames") <- "PotRad"
    attr(PotRadiation.V.n, "units") <- attr(ExtRadiation.V.n, 
                                            "units")
    PotRadiation.V.n
}



fCalcExtRadiation <- function (DoY.V.n) 
{
    FracYear_rad.V.n <- 2 * pi * (DoY.V.n - 1)/365.24
    SolarIrr_Wm2.c <- 1366.1
    ExtRadiation.V.n <- SolarIrr_Wm2.c * (1.00011 + 0.034221 * 
                                              cos(FracYear_rad.V.n) + 0.00128 * sin(FracYear_rad.V.n) + 
                                              0.000719 * cos(2 * FracYear_rad.V.n) + 7.7e-05 * sin(2 * 
                                                                                                       FracYear_rad.V.n))
    attr(ExtRadiation.V.n, "varnames") <- "ExtRad"
    attr(ExtRadiation.V.n, "units") <- "W_m-2"
    ExtRadiation.V.n
}



fCalcSunPosition <- function (DoY.V.n, Hour.V.n, Lat_deg.n, Long_deg.n, TimeZone_h.n, 
                              useSolartime.b = TRUE) 
{
    FracYear_rad.V.n <- 2 * pi * (DoY.V.n - 1)/365.24
    EqTime_h.V.n <- (0.0072 * cos(FracYear_rad.V.n) - 0.0528 * 
                         cos(2 * FracYear_rad.V.n) - 0.0012 * cos(3 * FracYear_rad.V.n) - 
                         0.1229 * sin(FracYear_rad.V.n) - 0.1565 * sin(2 * FracYear_rad.V.n) - 
                         0.0041 * sin(3 * FracYear_rad.V.n))
    LocTime_h.V.n <- (Long_deg.n/15 - TimeZone_h.n)
    SolTime_h.V.n <- if (useSolartime.b) {
        Hour.V.n + LocTime_h.V.n + EqTime_h.V.n
    }
    else {
        warning("Solar position calculated without correction for local time and equation of time.")
        Hour.V.n
    }
    SolTime_rad.V.n <- (SolTime_h.V.n - 12) * pi/12
    SolTime_rad.V.n <- ifelse(SolTime_rad.V.n < -pi, SolTime_rad.V.n + 
                                  2 * pi, SolTime_rad.V.n)
    attr(SolTime_h.V.n, "varnames") <- "SolTime"
    attr(SolTime_h.V.n, "units") <- "hour"
    SolDecl_rad.V.n <- ((0.33281 - 22.984 * cos(FracYear_rad.V.n) - 
                             0.3499 * cos(2 * FracYear_rad.V.n) - 0.1398 * cos(3 * 
                                                                                   FracYear_rad.V.n) + 3.7872 * sin(FracYear_rad.V.n) + 
                             0.03205 * sin(2 * FracYear_rad.V.n) + 0.07187 * sin(3 * 
                                                                                     FracYear_rad.V.n))/180 * pi)
    attr(SolDecl_rad.V.n, "varnames") <- "SolDecl"
    attr(SolDecl_rad.V.n, "units") <- "rad"
    SolElev_rad.V.n <- asin(sin(SolDecl_rad.V.n) * sin(Lat_deg.n/180 * 
                                                           pi) + cos(SolDecl_rad.V.n) * cos(Lat_deg.n/180 * pi) * 
                                cos(SolTime_rad.V.n))
    attr(SolElev_rad.V.n, "varnames") <- "SolElev"
    attr(SolElev_rad.V.n, "units") <- "rad"
    SolAzim_cos.V.n <- ((cos(SolDecl_rad.V.n) * cos(SolTime_rad.V.n) - 
                             sin(SolElev_rad.V.n) * cos(Lat_deg.n/180 * pi))/(sin(Lat_deg.n/180 * 
                                                                                      pi) * cos(SolElev_rad.V.n)))
    SolAzim_cos.V.n[SolAzim_cos.V.n > +1] <- 1
    SolAzim_cos.V.n[SolAzim_cos.V.n < -1] <- 1
    SolAzim_rad.V.n <- acos(SolAzim_cos.V.n)
    SolAzim_rad.V.n <- ifelse(SolTime_rad.V.n < 0, pi - SolAzim_rad.V.n, 
                              pi + SolAzim_rad.V.n)
    attr(SolAzim_cos.V.n, "varnames") <- "SolAzim"
    attr(SolAzim_cos.V.n, "units") <- "rad"
    SolPosition.L <- list(SolTime = SolTime_h.V.n, SolDecl = SolDecl_rad.V.n, 
                          SolElev = SolElev_rad.V.n, SolAzim = SolAzim_rad.V.n)
    SolPosition.L
}
