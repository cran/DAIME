pointtransform=function(points,xdep,ydep,direction='time to height',depositionmodel='piecewise linear deposition rate',hiatuslist=list(),unit='time per sediment',timetype='time'){
  #### Rough input check: is.finiteness, is.ordered, and length####
    if(!all(is.finite(c(points,ydep,xdep)))){
      stop("Values in points, xdep, and ydep need to be finite")
    }
    if (length(xdep)<2){
      stop("Degenerate input: xdep needs to contain at least two elements")
    }
    if (timetype!='time' & timetype!='age'){
      stop("timetype needs to be either \"time\" or \"age\" ")
    }
    if (timetype=='age'){
      if(direction=='time to height'){
        points=-points
        xdep=-xdep
      }
      if(direction=="height to time"&depositionmodel=='age model'){
        ydep=-ydep
      }
    }
    if(any(diff(xdep)<=0)){
      stop("Need strictly increasing vector for xdep")
    }
    lims=c(min(xdep),max(xdep))
    relevantpoints=(points>=lims[1] & points<=lims[2])
    xt=points[relevantpoints]
    hiatii=numeric()
    
  #### Check mode of deposition model and compatibility of inputs ####
    if (depositionmodel=='binned deposition rate'){
      if(length(xdep)!=(length(ydep)+1)){
        stop("For binned deposition rates, length(xdep)==length(ydep)+1 is required")
      }
    }
    else if (depositionmodel=='age model' | depositionmodel=='piecewise linear deposition rate'){
      if(length(xdep)!=length(ydep)){
        stop("Age models and piecewise linear deposition rates require length(xdep)==length(ydep)")
      }
    }
    else{
      stop("Incompatible mode for deposition model. Set depositionmodel to either \"piecewise linear deposition rate\",\"binned deposition rate\", or to \"age model\". ")
    }
    
  #### Check direction of transformation ####
    if (direction=='height to time'){
      if (depositionmodel=='age model'){
        if( any(diff(ydep)<=0) ){
          stop("Ambiguous age model for the transformation from height to time, need strictly increasing ydep. Please use the input option \"hiatuslist\" to define hiatuses")
        }
      }
      else{
        if (any(ydep<=0) ){
          stop("Negative or zero deposition rates are impossible in the section, need strictly positive values of ydep. Please use the input option \"hiatuslist\" to define hiatuses when the transformation is from height to time")
        }
        if(unit=='sediment per time'){
          ydep=ydep^{-1}
        }
        else if (unit!='time per sediment'){
          stop("Incompatible unit. Set unit to either \"sediment per time\" or to \"time per sediment\"")
        }
      }
    }
    else if (direction=='time to height'){
      #everything is fine
    }
    else{
      stop("Incompatible direction of transformation. Set direction to either \"time to height\" or to \"height to time\"")
    }
  
  #### Check Hiatuslist ####
  if(direction=='height to time'){
    if(!is.list(hiatuslist)){
      stop("hiatuslist needs be a list")
    }
    if (length(hiatuslist)>0 ){
      if ( any(sapply(hiatuslist,length)!=2) ){ #do all entries of the list have 2 components (1 for strat. height, 1 for duration?)
        stop("Elements of hiatuslist need to be of length two: first entry corresponding to stratigraphic height, second to duration of hiatus")
      }
      if(!all(sapply(hiatuslist,function(x) all(is.finite(x))))){
        stop("Non-finite entries in hiatuslist")
      }
      if(any(unlist(sapply(hiatuslist,function(x) tail(x,1)))<=0)){
        stop("Second entries of the vectors in hiatuslist need to be strictly positive, since hiatuses of negative duration make no sense!")
      }
      hiatheight=unlist(sapply(hiatuslist,function(x) head(x,1))) #get stratigraphic height of all hiatuses
      hiatdur=unlist(sapply(hiatuslist,function(x) tail(x,1)))
      relhiat=unique(hiatheight[hiatheight > lims[1]  &  hiatheight < lims[2]]) #relevant hiatuses
      hiatuslist=list()
      #remove hiatuses that are defined outside the stratigraphic height and combine those with identical heights
      if(length(relhiat)!=0){
        for(i in 1:length(relhiat)){
          hiatuslist[[i]]=c(relhiat[i],sum(hiatdur[hiatheight==relhiat[i]]))
        }
      }
    }
  }
  
  #### Build age models for transformation ####
    if(depositionmodel=='binned deposition rate'){
      xvals=xdep
      intvals=cumsum(c(0,diff(xvals)*ydep))
    }
    else if(depositionmodel=='piecewise linear deposition rate'){
      #determine all relevant function values of the depositon rate
      xvals=sort(unique(c(xdep,xt)))
      yvals=approx(xdep,ydep,xout=xvals,yleft=0,yright=0)$y
      #build age model by integrating over the deposition rate using the trapezoidal rule
      intvals = c(0,cumsum((0.5*(c(0,yvals)+c(yvals,0))*(c(xvals,0)-c(0,xvals)))[2:length(xvals)])) 
    }
    else if(depositionmodel=='age model'){
      xvals=xdep
      intvals=ydep
    }
  
  #### Insert hiatuses for transformation from time to height ####
  if(direction == 'time to height'){
    erosionpresent=any((depositionmodel=='piecewise linear deposition rate'&min(ydep)<=0),(depositionmodel=='binned deposition rate'&(min(ydep)<=0)),(depositionmodel=='age model' & (min(diff(c(ydep,(max(ydep)+1))))<=0)))
    #error in next line: no non-missing arguments to min; returning Inf DOUBLECHECK
    if (erosionpresent){ 
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
          if(depositionmodel=='age model' | depositionmodel=='binned deposition rate'){
            m=(intvals[highindex-1]-intvals[highindex-2])/(x2-x1)
            c=intvals[highindex-2]-m*x1
            xnew=(highpoint-c)/m  
          }
          else{ #if depomodel is piecewise linear
            y=approx(xdep,ydep,xout=c(x1,x2))$y 
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
  #### Insert hiatused for transformation from height to time ####
  if(direction == 'height to time'){
    #insert hiatuses 
    if(length(hiatuslist)>0){
      hiatheight=unlist(sapply(hiatuslist,function(x) head(x,1)))
      for (i in 1:length(hiatuslist)){ #insert hiatuses into age model
        #get time at which the hiatus begins
        jumpval=approx(xvals,intvals,xout=hiatheight[i],ties="ordered")$y 
        #split age model: one part below the hiatus to the beginning of the hiatus, one part from the end of the hiatus to the top. the second part is shifted by the duration of the hiatus
        intvals=c(intvals[ xvals<hiatheight[i] ],jumpval,jumpval+hiatuslist[[i]][2],intvals[ xvals>hiatheight[i] ]+hiatuslist[[i]][2] ) 
        #adjust x values to match the new age model
        xvals=c(xvals[ xvals<hiatheight[i] ],hiatheight[i],hiatheight[i],xvals[ xvals>hiatheight[i]])
      }
    }
  }

  #### Transform points ####
  xtrans=approx(xvals,intvals,xout=xt,ties="ordered")$y #transform points 
  if (direction=='time to height'){ #Remove values coinciding with hiatuses
    if(length(hiatii)>0){
      hiatii=unique(hiatii)
      for(i in 1:length(hiatii)){
        xtrans[xtrans==hiatii[i]]=NA
      }
    }
    height=replace(rep(NA,length(points)),relevantpoints,xtrans)
    age=points
  }
  else if (direction =='height to time'){
    #remove values that are located right on a hiatus
    if(length(hiatuslist)>0){
      for(j in 1:length(hiatheight)){
        xtrans[xt==hiatheight[j]]=NA
      }
    }
    age=replace(rep(NA,length(points)),relevantpoints,xtrans)
    height=points
  }
  #### Write summary ####
  summarysentence=paste("Transformed",as.character(length(points)),"point(s) from ",direction,"using a",depositionmodel, "as deposition model.")  
  if(direction=='time to height'){
    summarysentence=paste(summarysentence,as.character(length(hiatii[hiatii!=min(xdep)])),"hiatus(es) in the section.")
  }
  if(direction=="height to time"){
  summarysentence=paste(summarysentence,as.character(length(hiatuslist)),"hiatus(es) in the section.")
    if(depositionmodel !='age model'){
      if(unit=='sediment per time'){
        summarysentence=paste(summarysentence,"Interpreted deposition rate as deposition rate at a height in the section having the units",unit)
      }
      if(unit=='time per sediment'){
        summarysentence=paste(summarysentence,"Interpreted deposition rate as inverse deposition rate having the units",unit)
      }
    }
  }
  #### Output ####
  if (timetype=='age'){
    outlist=list(age=-age,height=height,report=paste(summarysentence, "Results given in age."))
  }
  if (timetype=='time'){
    outlist=list(time=age,height=height,report=paste(summarysentence, "Results given in time."))
  }
    
  return(outlist)
}