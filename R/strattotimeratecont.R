strattotimeratecont <-
function(xdep,ydep,xsig,ysig,pos=NULL,hiatuslist=list(),unit="sediment per time"){
  xpat=xsig
  ypat=ysig
  stopifnot(is.list(hiatuslist),all(is.finite(xdep)),all(is.finite(ydep)),all(is.finite(xpat)),all(is.finite(ypat)),is.unsorted(xdep,strictly = TRUE)==FALSE,is.unsorted(xpat,strictly = TRUE)==FALSE,min(ypat)>=0,min(ydep)>0,length(xdep)>=2,any(is.null(pos),all(is.finite(pos))),length(xdep)==length(ydep),length(xpat)==length(ypat))
  lims=c(max(min(xdep),min(xpat)),min(max(xdep),max(xpat))) #limits where both signal and deporate are defined
  
  if (length(hiatuslist)>0 ){
    if ( any(sapply(hiatuslist,length)!=2,unlist(sapply(hiatuslist,function(x) tail(x,1)))<=0) ){ #do all entries of the list have 2 components (1 for strat. height, 1 for duration?)
      stop("Incompatible input format of hiatuses. Please check help page")
    }
    hiatheight=unlist(sapply(hiatuslist,function(x) head(x,1))) #get stratigraphic height of all hiatuses
    #remove hiatuses that are defined outside the stratigraphic height
    hiatuslist=hiatuslist[ hiatheight > lims[1]  &  hiatheight < lims[2] ]
  }
  if(any(unit=="sediment per time",unit=='time per sediment')==FALSE){
    stop("Error: Incompatible unit (either \"sediment per time\" or \"time per sediment\")")
  }
  ll=patterntransform(xdep,ydep,xpat,ypat,direction='height to time',depositionmodel='piecewise linear deposition rate',patternmode='piecewise linear',pos=pos,unit=unit,hiatuslist = hiatuslist)
  outlist=list(age=ll$time,val=ll$val)
  return(outlist)
}
