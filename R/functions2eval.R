




##' Downward Longwave Radiation from a emissivity and a cloud cover schemes
##' 
##' @param data_ Data frame with column of date (date), temperature (Ta), 
##' partial vapor pressure (es), potencial radiation (Rpot), relative humidity (rh), 
##' atennuation index (K)
##' @param E_fun Emissivity scheme
##' @param C_fun Cloud cover scheme
##' @param adjust FALSE, TRUE if nonlinear least square adjusting wanted
##' @param method "non-linear" (default) for Non linear Least Square adjust, 
##' "montecarlo" for MonteCarlo optimization. Later be usseful when a NLS can't 
##' adjust observed data allowing optimization.
##' @param nsample population number evaluated in each iteration 
##' (only when method = "montecarlo").
##' @param max_iter maximun number of iterations (only when method = "montecarlo").
##' @param stats statistical function to be minimized (only when method = "montecarlo"),
##' NOTE: the best result should be 0.0 (ex., if stats = r (correlation), then transform to 
##' rMod = 1.0 - r, so the best result is when r== 1.0, so rMod == 0.0)
##' @param log_file path to a log file
##' @return Vector with Downward Longwave Radiation time series
##' @export
##' @import readr
get.Li <- function(data_,
                   E_fun = "EAN",
                   C_fun = "-",
                   adjust = FALSE,
                   method = c("non-linear","montecarlo"),
                   nsample = 1000,
                   max_iter = 10,
                   stats = "rmse",
                   log_file = NULL){
    
    sigma <- 5.67051*10^(-8) # W m^(-2) T^(-4)
    
        message("Combining emissivity: ", E_fun, ", whith cloud from: ", C_fun)
    
        if(adjust){
            message("      Making adjusting.")
            
            if( method[1] == "non-linear" ){
                message("           Non-Linear Adjusting")
                emis_ <- try( do.call(E_fun,list(data=data_, func = C_fun, adjust = adjust)),
                              silent = TRUE)
                
            } else if ( method[1] == "montecarlo" ){ 
                
                message("           MonteCarlo Adjusting")
                emis_ <- do.call(E_fun,list(data=data_, 
                                            func = C_fun, 
                                            adjust = adjust,
                                            nsample = nsample,
                                            method = "montecarlo",
                                            max_iter = max_iter,
                                            stats = stats))
                
            }
            
            
            if(class(emis_) == "try-error" & !is.na(method[2] == "montecarlo") ){
                
                message("           Error in NLS, passing to MonteCarlo approach >")
                emis_ <- do.call(E_fun,list(data=data_, 
                                            func = C_fun, 
                                            adjust = adjust,
                                            nsample = nsample,
                                            method = "montecarlo",
                                            max_iter = max_iter,
                                            stats = stats))
                
            }
            
        } else {
            
            emis_ <- try( do.call(E_fun,list(data=data_, func = C_fun)), silent = TRUE)
            
        }
        
        
        if(class(emis_) == "try-error") {
            emis_ <- list(emiss = NA,coefs = NA)
            roli_est <-  with(data_, emis_$emiss*sigma*Ta^4)
        } else {  
            roli_est <-  with(data_, emis_*sigma*Ta^4) 
        }
    #     
    # } else {
    #     if(class(emis_) == "try-error") emis_ <- NA
    #     roli_est <-  with(data_, emis_*sigma*Ta^4)   
    # }
    
    if(!is.null(log_file)){
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

##' Function to get all scheme calculation, adjust = TRUE make a adjusting NLS.
##' @param data Data frame with atmospheric variables
##' @param Ovrcst_sch Schemes for cloud cover index
##' @param Emiss_sch Schemes of atmosphere emissivity
##' @param adjust FALSE, TRUE for NLS adjusting
##' @param log_file path to a log file
##' @param method "non-linear" (default) for Non linear Least Square adjust,
##' "montecarlo" for MonteCarlo optimization. Later be usseful when a NLS can't 
##' adjust observed data allowing optimization.
##' @param nsample population number evaluated in each iteration 
##' (only when method = "montecarlo").
##' @param max_iter maximun number of iterations (only when method = "montecarlo").
##' @param stats statistical function to be minimized (only when method = "montecarlo"),
##' NOTE: the best result should be 0.0 (ex., if stats = r (correlation), then transform to 
##' rMod = 1.0 - r, so the best result is when r== 1.0, so rMod == 0.0)
#' @return Data frame with Li observed and all combintions of schemes for calculations of Li
#' @author Roilan Hernandez
#' @export
#' @importFrom dplyr bind_cols 
#' @import readr
get.AllSchems <- function(data,
                          Ovrcst_sch = c("CQB","CKC","CCB","CKZ","CWU","CJG", "CLM", "CFG"),
                          Emiss_sch  = c("EAN","EBR","EDO","EGR","EIJ","EID","EKZ","ENM",
                                         "EPR","EST","ESW"),
                          adjust = FALSE,
                          method = c("non-linear","montecarlo"),
                          nsample = 1000,
                          max_iter = 10,
                          stats = "rmse",
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
            # i = 1
            tmp.Li <- 
            get.Li(data_ = data,
                   E_fun = roli_comb[i,1] %>% as.character,
                   C_fun = roli_comb[i,2] %>% as.character,
                   adjust = adjust,
                   method = method,
                   nsample = nsample,
                   max_iter = max_iter,
                   stats = stats,
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
        mutate(daytime = to.daylight(date,lon = lon,lat = lat,timezone = timezone)) %>%
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



#' Calculate statistical index 
#' @param data_ Data frame with column Li observed and all series for each scheme.
#' @param statistic Statistical index to calculate
#' @param round Digits for round
#' @importFrom tidyr gather separate spread
#' @importFrom dplyr rename mutate %>% group_by_ summarise ungroup select
#' @import hydroGOF
#' @export
table_stats <- function(data_, statistic = "rmse", round = 2){
    statistic <- statistic[1]
    
    llply(data_,function(site){ 
        site %>%
            select(Li,starts_with("E",ignore.case = FALSE)) %>%
            gather(SCHEMES, Sim, -Li) %>%
            group_by(SCHEMES) %>%
            summarise(RESULT = do.call(statistic, list(obs = Li, sim = Sim, na.rm = TRUE)) %>% round(round)) %>% 
            separate(SCHEMES,into = c("Emissivity","CloudCover"),sep = "_") %>%
            mutate(CloudCover = ifelse(is.na(CloudCover),"-",CloudCover)) %>%
            spread(CloudCover,RESULT) 
    },.progress = "text")
    
}

##' Calculate emissivity
##' @param E_fun qkbf
##' @param data ADKVN
##' @param func SLKDJhvb
##' @param new.coefs qeohv
##' @export
run_fun <- function(E_fun,data, func,new.coefs){
    do.call(E_fun,
            as.list(modifyList(formals(E_fun), 
                               c(list(data = data, func = func), 
                                 as.list(new.coefs)))))
}


##' Calculate statistics
##' @param stats qkbf
##' @param obs ADKVN
##' @param sim SLKDJhvb
##' @import hydroGOF
##' @export
run_stats <- function(stats,obs,sim){
    do.call(stats,
            as.list(modifyList(formals(stats), 
                               c(list(sim = sim,
                                      obs = obs, 
                                      na.rm = TRUE)
                                 )
                               )
                    )
            )
}




##' MonteCarlo simulation for optimization
##' 
##' Make a MonteCarlo simulation from initial value on each schemes' coeficients. Initial sample 
##' has a 'nsample'*10 population number for each parameter in scheme, determined by a LHS 
##' (Latin Hypercubic Sample) with frontiers +- 5 times the value of coeficient:
##' $$ A_i = [a_i-5a_i;a_i+5a_i] $$
##' The MonteCarlo design follows a similar way in genetic algorithms, where a initial population 
##' of value for  parameter produce a L_i serie that is evaluated with choosen estatistic (stats).
##' Always two equal dimenson population are compared and choosen combinations of parameters that 
##' produce a 'stats' minor than a o.5 quantile. That way ensures that best combinations be conserved
##' if its performance is better than a half of total combinations 'stats' value. Finally a coeficients 
##' combination with best performance prevales and is considered the optimization. Furthermore, 
##' a 'nsample', 'max_iter' and 'stats' will define a result.  
##' @param data Data frame with all atmospherics variables
##' @param E_fun Function for emissivity 
##' @param func Function from cloud cover calculation
##' @param coefs Coeficients in E_fun function
##' @param nsample Number of population individuous
##' @param max_iter Maximun number of iterations
##' @param stats Statistical index to be minimized
##' @importFrom  parallel mclapply detectCores
##' @importFrom  dplyr bind_rows arrange filter_
##' @import stats
##' @import utils
##' @export
##' @return Vector with best combination of parameters
MonteCarlo <- function(data,
                       E_fun,
                       func,
                       coefs,
                       nsample ,
                       max_iter ,
                       stats){
    
    sigmaSB <- 5.67051 * 10^(-8)
    
    params <- 
        LHSU(xmax = coefs + coefs*5,
             xmin = coefs - coefs*5,
             nsample = nsample*10)  %>%
        as.data.frame() %>%
        setNames(names(coefs))                                                       ### **
    
    params[1,] <- coefs
    
    population_father <- NULL
    
        cat(" <")
    
    for(iter in 1:max_iter){
        cat("-")
        population_child <- 
            mclapply(1:length(params[[1]]),                                              ### **
                     function(i){
                         # i = 100
                         emiss <- with(data = data, 
                                       run_fun(E_fun = E_fun,
                                               data = data, 
                                               func = func,
                                               new.coefs = params[i,]) )
                         
                         stats_emiss <- run_stats(stats = stats,
                                                  sim = with(data, emiss * (sigmaSB*Ta^4)) %>% 
                                                      as.numeric(),
                                                  obs = with(data, Li) %>% as.numeric())
                         
                         return(cbind(params[i,], stats_var = stats_emiss) )
                    
                     },
                     mc.cores = max( detectCores()-1,1),                       ### **
                     mc.preschedule = FALSE) %>% 
            bind_rows()  ### **
        
        if(!is.null(population_father)){ 
            
            population_child <- 
                bind_rows(population_child,                                     ### **
                          population_father) %>%
                mutate(q.5 = quantile(stats_var , probs = 0.5, na.rm = TRUE)) %>%
                subset(stats_var < q.5) %>%               ### **       ### **
                arrange(stats_var) %>%
                select(-q.5)                                   ### **
            
        }
        
        params <- 
            LHSU(xmax = apply(population_child, 2, max, na.rm = TRUE)[1:length(coefs)],                                     ### **
                 xmin = apply(population_child, 2, min, na.rm = TRUE)[1:length(coefs)],                                     ### **
                 nsample = nsample)  %>%
            as.data.frame() %>%
            setNames(names(coefs))                                     ### **
        
        population_father <- population_child
        
    } ## acaba o for 
    cat(">\n")
    return(population_father[1,1:length(coefs)])
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