
#  FUNÇÕES NECESSÁRIAS
#
maxlim <- function(i,max_=1,min_=0){ 
    sapply(i,function(i) min(max(i,min_),max_) ) 
}
