timetostratratecont <-
function(xdep,ydep,xsig,ysig,pos=NULL){
  #check input
  {
    stopifnot(is.unsorted(xdep,strictly=TRUE)==FALSE,is.unsorted(xsig,strictly=TRUE)==FALSE,length(xdep)==length(ydep),length(xsig)==length(ysig),length(xdep)==length(xsig),min(ysig)>=0)
  }
  #preparing input data
  {
    lims=c(max(min(xdep),min(xsig)),min(max(xdep),max(xsig))) #limits where both signal and deporate are defined
    if(is.null(pos)){ #default setting for pos: if no points for evaluation are given
      pos=seq(from=lims[1],to=lims[2],length.out=min(1000,10*length(c(xdep,xsig))))
    }
    #remove points for which both the deposition rate and the signal are undefined
    relevantpoints=(pos>=lims[1] & pos<=lims[2])
    xt=pos[relevantpoints]
    #initialize storage for potential hiatuses
    hiatii=numeric()
    #determine all relevant function values of the depositon rate
    xvals=sort(unique(c(xdep,xt,xsig)))
    yvals=approx(xdep,ydep,xout=xvals,yleft=0,yright=0)[[2]]
  }
  #build age model by integrating over the deposition rate using the trapezoidal rule
  {
    intvals = c(0,cumsum((0.5*(c(0,yvals)+c(yvals,0))*(c(xvals,0)-c(0,xvals)))[2:length(xvals)]))                              
  }
  #if sediment is removed/sediment accumulation is stagnating, adjust age model:
  {
    if (min(ydep)<=0){ 
      #remove sections where sedimetn is removed
      ##define starting values of running indexes for the removal of sediment
      xvals=c(min(xvals),min(xvals),xvals) #duplicate values to avoid index malfunction later
      intvals=c(min(intvals),min(intvals),intvals) #duplicate values to avoid index malfunction later
      highpoint=intvals[length(intvals)] #height of last point deposited 
      highindex=length(intvals) #index of last deposited point
      while(highindex>2){   ##### procedure is runnign "backwards", e.g. from latest deposited sedimetn to earliest deposited sedimetn
        if(intvals[highindex-1] < highpoint){ #if sediment is only accumulated
          highindex=highindex-1 #change running indexes
          highpoint=intvals[highindex]
        }
        else if (intvals[highindex-1]==highpoint){ #if deposition stagnates
          highindex=highindex-1
          hiatii=c(hiatii, highpoint)
        }
        else if(intvals[highindex-1] > highpoint & intvals[highindex-2]>= highpoint ){ #if sedimetn is removed or not accumulated in the two prior steps
          hiatii=c(hiatii,highpoint) #store information about location of erosion 
          highindex=highindex-1 #change running indices
          intvals[highindex]=highpoint #adjust sedimetn to the highpoint, i.e. the highest point that will not be removed after the hiatus /erode
        }
        else if(intvals[highindex-1]>highpoint & intvals[highindex-2]<highpoint){ #if transition of sediment removal to accumulation
          x1=xvals[highindex-2]
          x2=xvals[highindex-1]
          y=approx(xdep,ydep,xout=c(x1,x2))[[2]]
          if(y[1]==y[2]){
            m=(intvals[highindex-1]-intvals[highindex-2])/(x2-x1)
            c=intvals[highindex-2]-m*x1
            xnew=(highpoint-c)/m
          }
          else{
            a=0.5*(y[2]-y[1])/(x2-x1)
            b=y[1]-2*a*x1
            c=-a*x2^2-b*x2+intvals[highindex-1]
            xnew=(-b+sqrt(b^2-4*a*(c-highpoint)))/(2*a)
            if(xnew < x1 | xnew>x2){
              xnew=(-b-sqrt(b^2-4*a*(c-highpoint)))/(2*a)
            }
          }
          intvals[highindex-1]=highpoint 
          intvals=c(intvals[1:(highindex-2)],highpoint,intvals[(highindex-1):length(intvals)])
          xvals=c(xvals[1:(highindex-2)],xnew,xvals[(highindex-1):length(xvals)]) 
          highindex=highindex-1
          hiatii=c(hiatii,highpoint)
        }
      }
      intvals=intvals[-c(1,2)] #remove additional values attached in the beginning
      xvals=xvals[-c(1,2)]
    }
  }
  #transform points from time into height
  {
    xtrans=approx(xvals,intvals,xout=xt)[[2]] #transform time into stratigraphic height 
    ytrans=approx(xsig,ysig,xout=xt)[[2]]/approx(xdep,ydep,xout=xt)[[2]] #values of the stratigraphic rate at the transformed points
    #remove points that coincide with a hiatus
    if(length(hiatii)>0){
      hiatii=unique(hiatii)
      for(i in 1:length(hiatii)){
        ytrans[xtrans==hiatii[i]]=NA
        xtrans[xtrans==hiatii[i]]=NA
      }
    }
  }
  #adjust output size to input size
  {
    xtrans=replace(rep(NA,length(relevantpoints)),relevantpoints,xtrans)
    ytrans=replace(rep(NA,length(relevantpoints)),relevantpoints,ytrans)
  }
return(list(height=xtrans,val=ytrans))
}
