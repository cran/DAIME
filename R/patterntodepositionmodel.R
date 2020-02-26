patterntodepositionmodel=function(xheight,yheight,xage=NULL,yage=NULL,heightmode='piecewise linear',agemode='piecewise linear',atheight=NULL,atage=NULL,rescalefor=1,timetype='time'){
  #### Check input and set default settings ####
  #initialize default settings for xage and yage: constant rate of 1 over the unit interval
  defaultage=FALSE
  if (timetype!='time' & timetype!='age'){
    stop("timetype needs to be either \"time\" or \"age\" ")
  }
  if (all(is.null(xage),is.null(yage))){
    if(timetype=='time'){
      xage=c(0,1)
    }
    if (timetype=='age'){
      xage=c(1,0)
    }
    yage=c(1,1)
    agemode='piecewise linear'
    defaultage=TRUE
  }
  ##catch input errors##
  if (xor(is.null(xage),is.null(yage))){ 
    stop("Either (1) enter input for both \"xage\" and \"yage\" to define temporal pattern or (2) skip input for both of them to use the default setting of a constant temporal pattern.")
  }
  if (timetype=='age'){
    xage=-xage
  }
  #check input: strictly positive values for both patterns, defined on strictly increasing x-values/binborders
  if(!all(is.finite(c(xage,yage,xheight,yheight)))){
    stop("Need finite inputs for xheight and yheight (and for xage and yage if used)")
  }
  #check for positivity
  if(any(c(yheight,yage)<=0) ){
    stop("Patterns must be positive. Need strictly positive values for yheight (and for yage if used)")
  }
  #check for correct length of xage and xheight
  if(length(xheight)<2 | length(xage)<2){
    stop("Degenerate input. Need length(xheight)>=2  (and length(xage)>=2 if xage is used)")
  }
  #check for strictly increasing values
  if(any(c(diff(xage),diff(xheight))<=0)){
    stop("Need strictly increasing vectors for xheight (and for xage if used)")
  }
  #### Define points in time/height where age model will be defined ####
  if(is.null(atheight)){
    atheight=seq(min(xheight),max(xheight),length.out=100)
  }
  else{
    if(!all(is.finite(atheight))){
      stop("Need finite values for atheight")
    }
  }
  relheights=atheight[atheight>=min(xheight) & atheight <= max(xheight)]
  
  if(is.null(atage)){
    atage=seq(min(xage),max(xage),length.out=100)
  }
  else{
    if(!all(is.finite(atage))){
      stop("Need finite values for atage")
    }
  }
  relages=atage[atage>=min(xage) & atage <= max(xage)]

  #### check and prepare inputs for the stratigraphic pattern ####
  heightvals=unique(sort(c(xheight,relheights)))
  if (heightmode=='piecewise linear'){
    if(length(xheight)!=length(yheight)){
      stop("For piecewise linear stratigraphic patterns, length(xheight) must match length(yheight)")
    }
    yheight=approx(xheight,yheight,xout=heightvals)$y
    xheight=heightvals
    intvals.height = c(0,cumsum((0.5*(c(0,yheight)+c(yheight,0))*(c(xheight,0)-c(0,xheight)))[2:length(xheight)]))
  }
  else if(heightmode=='binned'){
    if(length(xheight)!=(length(yheight)+1)){
      stop("For binned stratigraphic patterns, length(xheight) must match length(yheight)+1")
    }
    intvals.height=cumsum(c(0,diff(xheight)*yheight))
    intvals.height=approx(xheight,intvals.height,xout=heightvals)$y
    xheight=heightvals
  }
  else{
    stop("Incompatible mode for the stratigraphic pattern: set heightmode to either \"piecewise linear\" or \"binned\". ")
  }

  
  #### check and prepare inputs for the temporal pattern ####
  agevals=unique(sort(c(xage,relages)))
  if (agemode=='piecewise linear'){
    if(length(xage)!=length(yage)){
      stop("For piecewise linear temporal patterns, length(xage) must match length(yage)")
    }
    yage=approx(xage,yage,xout=agevals)$y
    xage=agevals
    intvals.age = c(0,cumsum((0.5*(c(0,yage)+c(yage,0))*(c(xage,0)-c(0,xage)))[2:length(xage)]))
  }
  else if(agemode=='binned'){
    if(length(xage)!=(length(yage)+1)){
      stop("For binned temporal patterns, length(xage) must match length(yage)+1")
    }
    intvals.age=cumsum(c(0,diff(xage)*yage))
    intvals.age=approx(xage,intvals.age,xout=agevals)$y
    xage=agevals
  }
  else{
    stop("Incompatible mode for the temporal pattern: Set agemode to either \"piecewise linear\" or \"binned\". ")
  }

  
  ####rescale temporal and stratigraphic rates ####
  if (rescalefor=='temporal pattern'){
    rescaleby=max(intvals.height)
  }
  else if (rescalefor=='stratigraphic pattern'){
    rescaleby=max(intvals.height)
  }
  else{
    if(any(!is.numeric(rescalefor),is.infinite(rescalefor),rescalefor<=0)){
      stop("Please use either a strictly positive numeric value or \"temporal pattern\" or \"stratigraphic pattern\" for the input option \"rescalefor\"")
    }
    rescaleby=rescalefor
  }
  yage=rescaleby*yage/max(intvals.age)
  intvals.age=rescaleby*intvals.age/max(intvals.age)
  yheight=rescaleby*yheight/max(intvals.height)
  intvals.height=rescaleby*intvals.height/max(intvals.height)
  
  #### Define Auxiliary Function ####
  #subroutine invint: determine preimage of values under given a piecewise quadratic function
  #given a piecewise linear function pilinfun=approxfun(x,y)
  #invint determines at what values the integral over pilinfun takes on the values fvals
  #I are the values of the integral over pilinfun from -infty to x
  invint=function(x,y,I,fvals){
    a=0.5*diff(y)/diff(x)             #leading coefficients for the piecewise quadratic function
    b=y[-length(x)]-2*a*x[-length(x)]                     #coefficients of linear parts for the piecewise quadratic function
    c=-a*x[-length(x)]^2-b*x[-length(x)]+I[-length(x)]  #constant coefficients for the piecewise quadratic function
    sieve <- rep(0,length(fvals))
    for (j in 1:(length(I)-1)){                           #determine in which interval of the piecewise quadratic function the j-th inverted entry in fvals will be in.
      sieve = sieve + as.numeric( fvals >= I[j] )
    }
    ct = c[sieve] - fvals;                       #determine the coefficients (for each entry of fvals) for the quadratic equations to solve
    at = a[sieve]
    bt = b[sieve]
    root1 = (-bt+sqrt(bt^2-4*at*ct))/(2*at)        #solve quadratic equations
    root2 = (-bt-sqrt(bt^2-4*at*ct))/(2*at)
    mysample = root1
    for (k in (1:length(root1))[at!=0]){           #blblb
      if(root2[k] >= x[sieve[k]] & root2[k] <= x[sieve[k]+1] ){
        mysample[k] = root2[k]
      }
    }
    mysample[at==0] = -ct[at==0]/bt[at==0]         #when an interval of the piecewise quadratic function is linear, the inversion can be done analytically
    return(mysample)
  }
  
  #### get values that need to be added ####
  if (heightmode=='binned'){
    xadd.height=approx(intvals.height,xheight,xout=intvals.age)$y
  }
  else if(heightmode=='piecewise linear'){
    xadd.height=invint(xheight,yheight,intvals.height,intvals.age)
  }
  if(agemode=='binned'){
    xadd.age=approx(intvals.age,xage,xout=intvals.height)$y
  }
  else if (agemode=='piecewise linear'){
    xadd.age=invint(xage,yage,intvals.age,intvals.height)
  }
  #### sort results ####
  xheight.end=sort(c(xheight,xadd.height))
  xage.end=sort(c(xage,xadd.age))
  keepme=(!duplicated(xheight.end)) & (!duplicated(xage.end))
  age=xage.end[keepme]
  height=xheight.end[keepme]
  
  if (is.null(atheight) & is.null(atage)){
    if(!is.null(xage)){
      age2=approx(height,age,xout=height,yright=max(age),yleft=min(age))
      height2=approx(age,height,xout=age,yright=max(height),yleft=min(height))
      ageout=sort(c(age2$y,height2$x))
      heightout=sort(c(age2$x,height2$y))
    }
    else{
      heightout=xheight
      ageout=approx(height,age,xout=height,yright=max(age),yleft=min(age))$y
    }
  }
else{
  age2=approx(height,age,xout=atheight,yright=max(age),yleft=min(age))
  height2=approx(age,height,xout=atage,yright=max(height),yleft=min(height))
  ageout=sort(c(age2$y,height2$x))
  heightout=sort(c(age2$x,height2$y))
}
  #### Write summary ####
  summarysentence=paste("Generated age model based on the dilution/condensation of a",agemode,"temporal pattern relative to a",heightmode,"stratigraphic pattern.")
  if(defaultage){
    summarysentence=paste(summarysentence,"Using the default option of a temporal pattern that is constant for a duration of one time unit.")
  }
  #### Output ####
  keepme=(!duplicated(heightout)) & (!duplicated(ageout))
  if(timetype=='time'){
    outlist=list(time=ageout[keepme],height=heightout[keepme],report=paste(summarysentence, "Results given in time."))
  }
  if(timetype=='age'){
    outlist=list(age=ageout[keepme],height=heightout[keepme],report=paste(summarysentence, "Results given in age."))
  }
  return(outlist)
}
