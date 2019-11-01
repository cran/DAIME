timetostratratebin <-
function(binborder,depoval,signalval,pos=NULL){
  xdep=binborder
  ydep=depoval
  xpat=binborder
  ypat=signalval
  stopifnot(all(is.finite(xdep)),all(is.finite(ydep)),all(is.finite(xpat)),all(is.finite(ypat)),is.unsorted(xdep,strictly = TRUE)==FALSE,is.unsorted(xpat,strictly = TRUE)==FALSE,min(ypat)>=0,length(xdep)>=2,length(xpat)>=2,any(is.null(pos),all(is.finite(pos))),length(xdep)==(length(ydep)+1),length(xpat)==(length(ypat)+1))
  ll=patterntransform(xdep,ydep,xpat,ypat,direction='time to height',depositionmodel='binned deposition rate',patternmode='binned',pos=pos)
  return(ll)          
}
