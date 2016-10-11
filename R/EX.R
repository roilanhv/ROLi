# require(dplyr)
# require(magrittr)
# require(readr)
# require(tibble)
# require(openair)
# require(REddyProc)
# data2 <-
# read_rds("~/Dropbox/NEW_WORKS/EPGMET_2016/Dados/Dados_SM_Nov2013Set2015_1hora_posFBE.rds") %>%
#     as.tbl() %>%
#     dplyr::select(date, Rg_fill, Tar_fill,UR_fill,Li) %>%
#     rename(Rg = Rg_fill,
#            Ta = Tar_fill,
#            rh = UR_fill) %>%
#     # selectByDate(month = c(1,7)) %>%
#     mutate(Rpot = PotRad(date,lon=-53.76,lat=-29.72 )) %>%
#     mutate(es = (rh/100) * (0.6112 * exp((17.67*Ta)/(243.5+Ta))) *10) %>%
#     mutate(Ta = Ta+273.15) %>%
#     mutate(K = kloudines(Rg = Rg,dates = date,lon=-53.76,lat=-29.72 ))
# 
# data2
 # save(data2,file = "./data/sm_data.rda",compress = "xz")
# 
