
##' Downward Longwave Radiation from a emissivity and a cloud cover schemes
##' 
##' @param data Data frame with column of date (date), temperature (Ta), 
##' partial vapor pressure (es), potencial radiation (Rpot), relative humidity (rh), 
##' atennuation index (K)
##' @param E_fun Emissivity scheme
##' @param C_fun Cloud cover scheme
##' @return Vector with Downward Longwave Radiation time series
##' @examples 
##' # Downward longwave for Santa Maria site, January and July of  2014
##' Li_sm <- get.Li(data = data2,E_fun = "FHY",C_fun = "CQB")
##' head(Li_sm)
##' summary(Li_sm)
##' 
get.Li <- function(data,E_fun = "FHY",C_fun = "CQB"){
    sigma <- 5.67051*10^(-8) # W m^(-2) T^(-4)
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
##' @param timezone Local time diference with GMT (-1 for fluxes measurement)
##' @importFrom dplyr %>% 
##' @return Vector with "day"/"night" string
##' @author Roilan Hernandez
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

##' Potencial Radiation from date vector
##' @param date A vector with data
##' @param lon Longitude from local analisys 
##' @param lat Latitude from local analisys
##' @param timezone Local time diference with GMT (-1 for fluxes measurement)
##' @return Vector with potential radiation for given dates and coordinates.
##' @importFrom dplyr %>% 
##' @author Roilan Hernandez
##' 
PotRad <-  function(date,lon=-53.76,lat=-29.72,timezone=-4){
    
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



# ##' Function for evaluation of all parameterizations
# ##' @param data_ A data frame with all atmospherics variables
# ##' @param lon Longitude from local analisys
# ##' @param lat Latitude from local analisys
# ##' @param timezone Local time diference with GMT (-1 for fluxes measurement)
# ##' @importFrom dplyr %>% filter_ filter mutate arrange_ bind_cols
# ##' @importFrom openair cutData
# ##' @importFrom magrittr set_names %<>%
# ##' @importFrom hydroGOF gof
# ##' @importFrom tidyr gather separate spread_
# ##' @return List with statistical error information between predicted and observed Li
# eval.params <- function(data_,lon=-53.76,lat=-29.72,timezone=-4){
# 
#     data_ <-
#         cutData(x = data_,type = "season",hemisphere = "southern") %>%
#         mutate(daytime = to.daylight(date,lon=lon,lat=lat,timezone=timezone))
# 
#     lapply(with(data_, unique(season)) %>% as.vector, function(j){
# 
#         season.data <-  data_ %>%
#             filter(season == j)
# 
#         lapply(with(season.data,unique(daytime)), function(k){
# 
#             in.data <-
#                 season.data %>%
#                 filter_("daytime" == k)
# 
#             estats.roli <-
#                 lapply(with(in.data,unique(params)), function(i){ # i = "FBM_CQB"
# 
#                     tdy.roli.filt <-
#                         in.data %>%
#                         filter_("params" == i)
# 
#                     gof.data <-
#                         gof(sim = with(tdy.roli.filt,value),
#                             obs = with(tdy.roli.filt,Li),
#                             na.rm = TRUE) %>%
#                         as.data.frame() %>%
#                         set_names(i)
# 
#                 }) %>% bind_cols()
# 
#             estats.roli %<>% mutate(stats = rownames(gof(1:10,10:1)) )
# 
#             estats.roli %<>%
#                 gather(params,value,-stats) %>%
#                 arrange_("stats")
# 
#             estats.roli.arrange <-
#                 estats.roli %>%
#                 separate(params,sep = "_",into = c("emis","aten")) %>%
#                 spread_("emis","value")
# 
#             estats.roli.arrange
# 
#         }) %>% set_names(with(data_,unique(daytime)))
#     }) %>% set_names(with(data_,unique(season)) %>% as.vector)
# }

# 
#    select_stats <- function(roli_list,idx = "RMSE"){
#         
#         stats.rmse <- NULL
#         
#         for(i in 1:4){
#             for(j in 1:2){
#                 season <- names(roli_list)[i]
#                 day.tim <- names(roli_list[[i]])[j]
#                 stats.rmse[[season]][[day.tim]] <- 
#                     roli_list[[season]][[day.tim]] %>% filter_("stats" == idx)
#                 
#             }
#         }
#         stats.rmse
#     }
#     