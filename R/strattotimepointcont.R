strattotimepointcont <-
function(x,xdep,ydep,hiatuslist=list(),unit="sediment per time"){
  stopifnot(all(is.finite(x)),all(is.finite(xdep)),all(is.finite(ydep)),is.unsorted(xdep,strictly = TRUE)==FALSE,is.list(hiatuslist),length(xdep)==length(ydep),min(ydep)>0)
  if(any(unit=="sediment per time",unit=='time per sediment')==FALSE){
    stop("Error: Incompatible unit (either \"sediment per time\" or \"time per sediment\")")
  }
  #check input for hiatuslist
  if (length(hiatuslist)>0 ){
    if ( any(sapply(hiatuslist,length)!=2,unlist(sapply(hiatuslist,function(x) tail(x,1)))<=0) ){ #do all entries of the list have 2 components (1 for strat. height, 1 for duration?, and are the durations positive?)
      stop("Incompatible input format of hiatuses. Please check help page")
    }
    hiatheight=unlist(sapply(hiatuslist,function(x) head(x,1))) #get stratigraphic height of all hiatuses
    hiatuslist=hiatuslist[hiatheight<max(xdep)&hiatheight>min(xdep)]
  }
  usedunit=unit
  usedhiatuslist=hiatuslist
  
  ll=pointtransform(points=x,xdep=xdep,ydep=ydep,direction='height to time',depositionmodel='piecewise linear deposition rate',hiatuslist=usedhiatuslist,unit=usedunit)
  return(ll[which(names(ll)!="report")])
}
