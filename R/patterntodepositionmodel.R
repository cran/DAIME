patterntodepositionmodel=function(xheight,yheight,xage=NULL,yage=NULL,heightmode='piecewise linear',agemode='piecewise linear'){
  #Check input and set default settings
  {
    #initialize default settings for xage and yage: constant rate of 1 over the unit interval
    if (all(is.null(xage),is.null(yage))){
      xage=c(0,1)
      yage=c(1,1)
    }
    #catch input error
    if (xor(is.null(xage),is.null(yage))){ 
      stop("Error. Please either (1) enter input for both \"xage\" and \"yage\" to define a pattern in time or (2) skip input for both of them to use the default setting of a constant pattern in time.")
    }
    #check input: strictly positive values for both patterns, defined on strictly increasing x-values/binborders
    stopifnot(all(is.finite(xage)),all(is.finite(yage)),all(is.finite(xheight)),all(is.finite(yheight)),min(yheight)>0,min(yage)>0,is.unsorted(xage,strictly = TRUE)==FALSE,is.unsorted(xheight,strictly = TRUE)==FALSE,length(xheight)>=2)
  }
  #prepare inputs for the pattern in stratigraphic height
  {
    heightvals=unique(sort(c(xheight,seq(min(xheight),max(xheight),length.out=100))))
    if (heightmode=='piecewise linear'){
      stopifnot(length(xheight)==length(yheight))
      yheight=approx(xheight,yheight,xout=heightvals)$y
      xheight=heightvals
      intvals.height = c(0,cumsum((0.5*(c(0,yheight)+c(yheight,0))*(c(xheight,0)-c(0,xheight)))[2:length(xheight)]))
    }
    else if(heightmode=='binned'){
      stopifnot(length(xheight)==(length(yheight)+1))
      intvals.height=cumsum(c(0,diff(xheight)*yheight))
      intvals.height=approx(xheight,intvals.height,xout=heightvals)$y
      xheight=heightvals
    }
    else{
      stop("Please enter a compatible mode for the pattern in stratigraphic height: Either \"piecewise linear\" or \"binned\". ")
    }
    #rescale the rates to have a input volume of 1
    yheight=yheight/max(intvals.height)
    intvals.height=intvals.height/max(intvals.height)
  }
  #check and prepare inputs for the pattern in time
  {
    agevals=unique(sort(c(xage,seq(min(xage),max(xage),length.out=100))))
    if (agemode=='piecewise linear'){
      stopifnot(length(xage)==length(yage))
      yage=approx(xage,yage,xout=agevals)$y
      xage=agevals
      intvals.age = c(0,cumsum((0.5*(c(0,yage)+c(yage,0))*(c(xage,0)-c(0,xage)))[2:length(xage)]))
    }
    else if(agemode=='binned'){
      stopifnot(length(xage)==(length(yage)+1))
      intvals.age=cumsum(c(0,diff(xage)*yage))
      intvals.age=approx(xage,intvals.age,xout=agevals)$y
      xage=agevals
    }
    else{
      stop("Please enter a compatible mode for the pattern in time: Either \"piecewise linear\" or \"binned\". ")
    }
    #rescale the rates to have a input volume of 1
    yage=yage/max(intvals.age)
    intvals.age=intvals.age/max(intvals.age)
  }
  {
    ###subroutine invint: determine preimage of values under given a piecewise quadratic function
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
  }
  #get x values that need to be added
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
  #sort results
  xheight.end=sort(c(xheight,xadd.height))
  xage.end=sort(c(xage,xadd.age))
  keepme=(!duplicated(xheight.end)) & (!duplicated(xage.end))
  return(list(age=xage.end[keepme],height=xheight.end[keepme]))
}
