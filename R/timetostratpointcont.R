timetostratpointcont <-
function(x,xdep,ydep){
  #determine all relevant y values of the depositon rate
  resl=approx(xdep,ydep,xout=sort(unique(c(xdep,x))),yleft=0,yright=0)
  #integrate deposition rate using the trapezoidal rule
  I0 = c(0,cumsum((0.5*(c(0,resl[[2]])+c(resl[[2]],0))*(c(resl[[1]],0)-c(0,resl[[1]])))[2:length(resl[[1]])]))                              
  #determine transformed points in the section
  return(approx(resl[[1]],I0,xout=x,yleft=0,yright=0)[[2]])
}
