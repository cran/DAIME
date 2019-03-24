timetostratratecont <-
function(xdep,ydep,xsig,ysig,pos=NULL){
  if(is.null(pos)==TRUE){pos=xsig}
  #determine all relevant y values of the depositon rate
  resl=approx(xdep,ydep,xout=sort(unique(c(xdep,xsig,pos))),yleft=0,yright=0)
  #integrate deposition rate using the trapezoidal rule
  I0 = c(0,cumsum((0.5*(c(0,resl[[2]])+c(resl[[2]],0))*(c(resl[[1]],0)-c(0,resl[[1]])))[2:length(resl[[1]])]))                              
  #determine where the transformed rate function is evaluated
  xtrans=approx(resl[[1]],I0,xout=pos)[[2]] 
  #get values of the transformed rate function
  ytrans=approx(xsig,ysig,xout=pos,yleft=0,yright=0)[[2]]/(approx(xdep,ydep,xout=pos,yleft=0,yright=0)[[2]]) 
return(list(x=xtrans,y=replace(ytrans,is.finite(ytrans)==FALSE,0)))
}
