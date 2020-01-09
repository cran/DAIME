timetostratpointcont <-
function(x,xdep,ydep){
  stopifnot(all(is.finite(x)),all(is.finite(xdep)),all(is.finite(ydep)),is.unsorted(xdep,strictly = TRUE)==FALSE,length(xdep)==length(ydep))
  ll=pointtransform(points=x,xdep=xdep,ydep=ydep,direction='time to height',depositionmodel='piecewise linear deposition rate')
  return(ll[which(names(ll)!="report")])
}
