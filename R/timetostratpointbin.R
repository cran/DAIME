timetostratpointbin <-
function(x,binborder,depoval){
  #check input
  {
    stopifnot(is.unsorted(binborder,strictly=TRUE)==FALSE,(length(binborder)-1)==length(depoval))
  }
  #preparing input data
  {
    lims=c(min(binborder),max(binborder)) #limits where both signal and deporate are defined
    #remove points for which both the deposition rate and the signal are undefined
    relevantpoints=(x>=lims[1] & x<=lims[2])
    xt=x[relevantpoints]
    #initialize storage for potential hiatuses
    hiatii=numeric()
    #determine all relevant x values
    xvals=binborder
  }
  #build age model by integrating over the deposition rate
  {
    intvals=cumsum(c(0,diff(xvals)*depoval))
  }
  #if sediment is removed/sediment accumulation is stagnating, adjust age model:
  {
    if(min(depoval)<=0){
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
        else if(intvals[highindex-1]>highpoint & intvals[highindex-2]<highpoint){ #case 4: transition from removal to accumulation
          m=(intvals[highindex-1]-intvals[highindex-2])/(xvals[highindex-1]-xvals[highindex-2]) #determien at which value the sedimetn crosses the value "highpoint"
          xnew=(highpoint-(intvals[highindex-2]-m*xvals[highindex-2]))/m
          intvals[highindex-1]=highpoint #remove sedimetn
          xvals=c(xvals[1:(highindex-2)],xnew,xvals[(highindex-1):length(xvals)]) #insert new bin boundary at the value where the sediemtn is taking on the value highpoint
          intvals=c(intvals[1:(highindex-2)],highpoint,intvals[(highindex-1):length(intvals)])
          hiatii=c(hiatii,highpoint)
          highindex=highindex-1
        }
      }
      xvals=xvals[-c(1,2)]
      intvals=intvals[-c(1,2)]
    }
  }
  #transform points from time into height
  {
    xtrans=approx(xvals,intvals,xout=xt)[[2]] #transform time into stratigraphic heigh
    #remove points that coincide with a hiatus
    if(length(hiatii)>0){
      hiatii=unique(hiatii)
      for(i in 1:length(hiatii)){
        xtrans[xtrans==hiatii[i]]=NA
      }
    }
  }
  #adjust output to input size
  {
    xtrans=replace(rep(NA,length(relevantpoints)),relevantpoints,xtrans)
  }
  return(list(height=xtrans,time=x))
}
