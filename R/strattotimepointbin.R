strattotimepointbin <-
function(x,binborder,depoval,hiatuslist=list(),unit="sediment per time"){
  stopifnot(all(is.finite(x)),all(is.finite(binborder)),all(is.finite(depoval)),is.unsorted(binborder,strictly = TRUE)==FALSE,is.list(hiatuslist),length(binborder)==(length(depoval)+1),min(depoval)>0)
  if(any(unit=="sediment per time",unit=='time per sediment')==FALSE){
    stop("Error: Incompatible unit (either \"sediment per time\" or \"time per sediment\")")
  }
  #check input for hiatuslist
  if (length(hiatuslist)>0 ){
    if ( any(sapply(hiatuslist,length)!=2,unlist(sapply(hiatuslist,function(x) tail(x,1)))<=0) ){ #do all entries of the list have 2 components (1 for strat. height, 1 for duration?, and are the durations positive?)
      stop("Incompatible input format of hiatuses. Please check help page")
    }
    hiatheight=unlist(sapply(hiatuslist,function(x) head(x,1))) #get stratigraphic height of all hiatuses
    hiatuslist=hiatuslist[hiatheight<max(binborder) & hiatheight>min(binborder)]
  }
  usedunit=unit
  usedhiatuslist=hiatuslist
  
  ll=pointtransform(points=x,xdep=binborder,ydep=depoval,direction='height to time',depositionmodel='binned deposition rate',hiatuslist=usedhiatuslist,unit=usedunit)
  return(ll[which(names(ll)!="report")])
}
