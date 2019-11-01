timetostratratecont <-
function(xdep,ydep,xsig,ysig,pos=NULL){
  xpat=xsig
  ypat=ysig
  stopifnot(all(is.finite(xdep)),all(is.finite(ydep)),all(is.finite(xpat)),all(is.finite(ypat)),is.unsorted(xdep,strictly = TRUE)==FALSE,is.unsorted(xpat,strictly = TRUE)==FALSE,min(ypat)>=0,length(xdep)>=2,length(xpat)>=2,any(is.null(pos),all(is.finite(pos))),length(xdep)==length(ydep),length(xpat)==length(ypat))
  
  ll=patterntransform(xdep,ydep,xpat,ypat,direction='time to height',depositionmodel='piecewise linear deposition rate',patternmode='piecewise linear',pos=pos)

return(ll)
}
