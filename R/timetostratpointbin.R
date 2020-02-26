timetostratpointbin <-
function(x,binborder,depoval){
  stopifnot(all(is.finite(x)),all(is.finite(binborder)),all(is.finite(depoval)),is.unsorted(binborder,strictly = TRUE)==FALSE,length(binborder)==(length(depoval)+1))
  ll=pointtransform(points=x,xdep=binborder,ydep=depoval,direction='time to height',depositionmodel='binned deposition rate')
  outlist=list(height=ll$height,age=ll$time)
  return(outlist)
}
