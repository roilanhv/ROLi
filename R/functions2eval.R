
##' Downward Longwave Radiation from a emissivity and a cloud cover schemes
##' 
##' @param data_ Data frame with column of date (date), temperature (Ta), 
##' partial vapor pressure (es), potencial radiation (Rpot), relative humidity (rh), 
##' atennuation index (K)
##' @param E_fun Emissivity scheme
##' @param C_fun Cloud cover scheme
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @param log_file path to a log file
##' @return Vector with Downward Longwave Radiation time series
##' @examples 
##' # Downward longwave for Santa Maria site, January and July of  2014
##' # without adjusting
##' Li_sm <- get.Li(data = data2,E_fun = "EAN",C_fun = "CQB")
##' head(Li_sm)
##' summary(Li_sm)
##' @export
##' @importFrom readr write_lines
get.Li <- function(data_,
                   E_fun = "EAN",C_fun = "CQB",
                   adjust = FALSE,
                   log_file = NULL){
    
    sigma <- 5.67051*10^(-8) # W m^(-2) T^(-4)
        message("Combining emissivity: ", E_fun, ", whith cloud from: ", C_fun)
    
    emis_ <- try(do.call(E_fun,list(data=data_,func = C_fun, adjust = adjust)),silent = TRUE)
    
    if(adjust){
        
        if(class(emis_) == "try-error") {
            emis_ <- list(emiss = NA,coefs = NA)
            roli_est <-  with(data_, emis_$emiss*sigma*Ta^4)
        } else {  roli_est <-  with(data_, emis_$emiss*sigma*Ta^4) }
        
    } else {
        if(class(emis_) == "try-error") emis_ <- NA
        roli_est <-  with(data_, emis_*sigma*Ta^4)   
    }
    
    if(!is.null(log_file)){
        # write_lines(x = paste0("\n Emissivity: ", E_fun, ", Cloud Cover: ", C_fun), path = "./AJUST.log", append = TRUE )
        # write_lines(x = paste0("\t",names(emis_$coefs),collapse = "\t\t"), path = "./AJUST.log", append = TRUE)
        write_lines(x = paste(E_fun,C_fun, paste0(emis_$coefs %>% round(7),collapse = " ")), path = log_file, append = TRUE)
    }
    
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
##' @export
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
##' @export
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

#' Function for statistical analise
#' 
#' Recommended use of function from hydroGOF package (ex., rmse, mae, pbias, me, nrmse, rSD...)
#' @param data_li Data frame with observed and calculated  downward longwave radiation series.
#' @param statistic Statistical index to be calculated between observation and each simulation, "rmse" default.
#' @param avg.time Temporal resolution in analysis, "hourly" for hourly observations or "daily" for convertion to mean daily values. 
#' @importFrom tidyr gather
#' @importFrom dplyr rename rename_ group_by summarise select summarise_each funs select_ bind_cols
#' @import hydroGOF
#' @examples
#' LiStats <- CalcStats(data_li = data_li, statistic = "pbias",avg.time = "daily")
#' LiStats
#' @export
CalcStats <- function(data_li, 
                      statistic = c("rmse","pbias"),
                      avg.time = "hourly"){
    
     Li <- params <- value <- obs <- sim <- NULL
    
    if(avg.time == "daily") {
        
        result <- 
            data_li %>% 
            gather( params, value, -Li, -date) %>%
            group_by(params,day = as.Date(date)) %>%
            summarise_each(funs(mean)) %>%
            select(date, Li, params, value) %>%
            rename(obs =Li, sim = value)
        
    } else {
 
    result <- 
        data_li %>%
        gather( params, value, -Li, -date) %>%
        rename_(obs = "Li", sim = "value") %>% 
        group_by(params) 
    }
    
    schems <- data.frame(schemes = unique(result$params))
    
    result <-
    lapply(1:length(statistic), function(i){

            result %>%
            summarise(RESULT = do.call(statistic[i], list(obs = obs, sim = sim, na.rm = TRUE))) %>%
            setNames(c("schemes",statistic[i])) %>%
            select_(statistic[i])
        
    }) %>% bind_cols()
    
    bind_cols(schems,result)
    
}

#' Function to get all scheme calculation, adjust = TRUE make a adjusting NLS.
#' @param data Data frame with atmospheric variables
#' @param Ovrcst_sch Schemes for cloud cover index
#' @param Emiss_sch Schemes of atmosphere emissivity 
#' @param adjust FALSE, TRUE for NLS adjusting 
#' @param log_file path to a log file
#' @return Data frame with Li observed and all combintions of schemes for calculations of Li
#' @author Roilan Hernandez
#' @export
#' @importFrom dplyr bind_cols 
#' @importFrom readr write_file write_lines
get.AllSchems <- function(data,
                          Ovrcst_sch = c("CQB","CKC","CCB","CKZ","CWU","CJG", "CLM", "CFG"),
                          Emiss_sch  = c("EAN","EBR","EDO","EGR","EIJ","EID","EKZ","ENM",
                                         "EPR","EST","ESW","EAI"),
                          adjust = FALSE, 
                          log_file = NULL){
    
    roli_comb <- rbind(expand.grid(Emiss_sch,"-"), expand.grid(Emiss_sch, Ovrcst_sch) )
    
    if(!is.null(log_file)){
    write_file(x = "\t START A NEW AJUST TASK", path = log_file, append = FALSE)
    write_lines(x = paste("\t",Sys.time(), "\n"), path = log_file, append = TRUE)
    write_lines(x = paste0("\t",length(Emiss_sch), " EMISSIVITY SCHEMES: ", 
                           paste0(Emiss_sch,collapse = "-") ),
                path = log_file, append = TRUE )
    write_lines(x = paste0("\t",length(Ovrcst_sch), " CLOUD COVER SCHEMES: ", 
                           paste0(Ovrcst_sch,collapse = "-") ),
                path = log_file, append = TRUE )
    write_lines(x = "EMISSIVITY CLOUD COEF1 COEF2 COEF3 COEF4 COEF5",
                path = log_file, append = TRUE)
    }
    
    Li.sims <- 
        lapply(1:nrow(roli_comb), function(i){
            
            tmp.Li <- 
            get.Li(data_ = data,
                   E_fun = roli_comb[i,1] %>% as.character,
                   C_fun = roli_comb[i,2] %>% as.character,
                   adjust = adjust,
                   log_file = log_file)
            
            return(tmp.Li)
            
        }) %>% bind_cols()
    
    names(Li.sims) <- gsub("_-","",names(Li.sims))
    
    cbind(data.frame(date = data$date), Li.sims)
}


#' Split data and Li series by time factor (ex., "daytime", "season", "cover")
#' @param data_ Data frame with column date (POSIXct), Obs and all series for each scheme.
#' @param split_class Type of cut data analises (daytime and season by default, other types needs as a factor column)
#' @param lat,lon,timezone Latitude, longitude e timezone of observations local.
#' @param round Digits for round
#' @importFrom openair cutData
#' @importFrom tidyr gather_
#' @importFrom dplyr rename mutate %>% group_by_ summarise ungroup
#' @import hydroGOF
#' @export
split_stats <- function(data_ ,
                        split_class = c("daytime","season","Cover","ID"),
                        lon = -53.18,lat = -29.71,timezone = -3, round = 3){
                         
    hemisphere <- ifelse(lat < 0.0, "southern", "northern")

    if(!("Obs" %in% names(data_))){ return(message("Column Obs isn't in the input data"))}

    output <-
    data_ %>%
        mutate(daytime = to.daylight(date,lon = lon,lat = lat,timezone = timezone)) #%>%
        cutData(type = "season",hemisphere = hemisphere) %>%
        gather_(key_col = "params", value_col = "Sim",
                gather_cols = names(data_)[!names(data_) %in% c("date",split_class,"Obs")]) %>%
        group_by_(.dots = c("params",split_class)) %>%
        summarise(RMSE = rmse(obs = Obs, sim = Sim, na.rm = TRUE) %>% round(round),
                  NRMSE = nrmse(obs = Obs, sim = Sim, na.rm = TRUE) %>% round(round),
                  rSD = (1- sd(Sim,na.rm = TRUE)/sd(Obs,na.rm = TRUE))  %>% round(round),
                  MAE = mae(obs = Obs, sim = Sim, na.rm = TRUE)%>% round(round),
                  PBIAS = pbias(obs = Obs, sim = Sim, na.rm = TRUE)%>% round(round),
                  NSE = NSE(obs = Obs, sim = Sim, na.rm = TRUE)%>% round(round),
                  R2 = cor(Sim, Obs, method = "pearson", use = "pairwise.complete.obs")^2 %>%
                      round(round),
                  NObs = sum(is.na(Sim),!is.na(Sim)) ,
                  PNAs = floor( (sum(is.na(Sim))/NObs)*100) ) %>%
        ungroup()

    return(output)
    }

# ##' Function for evaluation of all parameterizations
# ##' @param data_ A data frame with all atmospherics variables and simulated series
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