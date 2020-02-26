patterntransform=function(xdep,ydep,xpat,ypat,direction='time to height',depositionmodel='piecewise linear deposition rate',patternmode='piecewise linear',pos=NULL,hiatuslist=list(),unit='time per sediment',timetype='time'){
  #### Check input: finite, length xdep/xpat, increasing, define/check pos ####
  if(!all(is.finite(c(xpat,ypat,ydep,xdep)))){
    stop("Values in xpat, ypat, xdep, and ydep need to be finite")
  }
  if (length(xdep)<2 | length(xpat)<2){
    stop("Degenerate input: both xdep and xpat need to contain at least two elements")
  }
  
  if (timetype!='time' & timetype!='age'){
    stop("timetype needs to be either \"time\" or \"age\" ")
  }
  if (timetype=='age'){
    if(direction=='time to height'){
      xpat=-xpat
      xdep=-xdep
    }
    if(direction=="height to time"&depositionmodel=='age model'){
      ydep=-ydep
    }
  }
  
  if(any(c(diff(xdep),diff(xpat))<=0)){
    stop("Need strictly increasing vectors for xdep and xpat")
  }
  lims=c(max(min(xdep),min(xpat)),min(max(xdep),max(xpat)))
  if( is.null(pos) ){
    pos=seq(lims[1],lims[2],length.out = min(1000,10*length(c(xdep,xpat))))
  }
  if(!all(is.finite(pos))){
    stop("Need finite values for pos")
  }
  relevantpoints=(pos>=lims[1] & pos<=lims[2])
  xt=pos[relevantpoints]
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
    
  #### Check Patterns ####
  if (patternmode=='piecewise linear'){
    if(length(xpat)!=length(ypat)){
      stop("Piecewise linear patterns require length(xpat)==length(ypat)")
    }
  }
  else if(patternmode=='binned'){
    if(length(xpat)!=(length(ypat)+1)){
      stop("Binned patterns require length(xpat)==(length(ypat)+1)")
    }
  }
  else{
    stop("Incompatible mode for input of patterns. Set patternmode to either \"piecewise linear\" or \"binned\" )")
  }
  if(any(ypat<0)){
    stop("Patterns need to be positive, so ypat>=0 is required")
  }
  #### Check Hiatuslist ####
  if(direction=='height to time'){
    if(!is.list(hiatuslist)){
      stop("hiatuslist should be a list")
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
  xvals=sort(unique(c(xt,xdep,xpat)))
  if(depositionmodel=='binned deposition rate'){
    ydep2=cumsum(c(0,diff(xdep)*ydep))
    intvals=approx(xdep,ydep2,xout=xvals,yleft=0,yright=0)$y
  }
  else if(depositionmodel=='piecewise linear deposition rate'){
    #determine all relevant function values of the depositon rate
    ydep2=approx(xdep,ydep,xout=xvals,yleft=0,yright=0)$y
    #build age model by integrating over the deposition rate using the trapezoidal rule
    intvals = c(0,cumsum((0.5*(c(0,ydep2)+c(ydep2,0))*(c(xvals,0)-c(0,xvals)))[2:length(xvals)])) 
  }
  else if(depositionmodel=='age model'){
    intvals=approx(xdep,ydep,xout=xvals,yleft=0,yright=0)$y
  }
  #### Insert hiatuses for transformation from time to height ####
  if(direction == 'time to height'){
    erosionpresent=any((depositionmodel=='piecewise linear deposition rate'&min(ydep)<=0),(depositionmodel=='binned deposition rate'&(min(ydep)<=0)),(depositionmodel=='age model' & (min(diff(c(ydep,(max(ydep)+1))))<=0)))
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
  #### Insert hiatuses for transformation from height to time ####
  if(direction == 'height to time'){
    #insert hiatuses (if present)
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

  #### transform points ####
  xtrans=approx(xvals,intvals,xout=xt,ties="ordered")$y #transform points 
  if (direction=='time to height'){
    if(length(hiatii)>0){
      hiatii=unique(hiatii)
      for(i in 1:length(hiatii)){
        xtrans[xtrans==hiatii[i]]=NA
      }
    }
  }
  if (direction =='height to time'){
    #remove values that are located right on a hiatus
    if(length(hiatuslist)>0){
      for(j in 1:length(hiatheight)){
        xtrans[xt==hiatheight[j]]=NA
      }
    }
  }
  #### Calculate Values ####
  if(patternmode=='binned'){
    numerator=approx(xpat,c(ypat,ypat[length(ypat)]),xout=xt,method='constant')$y
  }
  else if(patternmode=='piecewise linear'){
    numerator=approx(xpat,ypat,xout=xt)$y
  }
  if (depositionmodel=='age model'){
    grad=diff(ydep)/diff(xdep)
    denominator=approx(xdep,c(grad,grad[length(grad)]),xout=xt,method='constant')$y
  }
  else if(depositionmodel=='binned deposition rate'){
    denominator=approx(xdep,c(ydep,ydep[length(ydep)]),xout=xt,method='constant')$y
  }
  else if(depositionmodel=='piecewise linear deposition rate'){
    denominator=approx(xdep,ydep,xout=xt)$y
  }
  ytrans=numerator/denominator
  ytrans[is.na(xtrans)]=NA
  
  #### Write summary ####
  summarysentence=paste("Transformed",patternmode,"pattern from ",direction,"using a",depositionmodel, "as deposition model.")  
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
  if(direction=='time to height'){
    if (timetype=='time'){
      outlist=list(height=xtrans,val=ytrans,report=paste(summarysentence,"Interpreted input as time"))
    }
    if (timetype=='age'){
      outlist=list(height=xtrans,val=ytrans,report=paste(summarysentence,"Interpreted input as age"))
    }
  }
  else if(direction=='height to time'){
    if (timetype=='age'){
      if (depositionmodel=='age model'){
        outlist=list(age=-xtrans,val=ytrans,report=paste(summarysentence, "Results given in age."))
      }
      else{
        outlist=list(age=-(xtrans-max(intvals)),val=ytrans,report=paste(summarysentence, "Results given in age."))
      }
      
    }
    if (timetype=='time'){
      outlist=list(time=xtrans,val=ytrans,report=paste(summarysentence, "Results given in time."))
    }
  }
  return(outlist)
}