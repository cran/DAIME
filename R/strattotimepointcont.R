strattotimepointcont <-
function(x,xdep,ydep,hiatuslist=list(),unit="sediment per time"){
  #check input
  {
    stopifnot(is.list(hiatuslist),is.unsorted(xdep,strictly=TRUE)==FALSE,min(ydep)>0,length(xdep)==length(ydep))
    #check input format of hiatuses
    if (length(hiatuslist)>0 ){
      if ( any(sapply(hiatuslist,length)!=2) ){ #do all entries of the list have 2 components (1 for strat. height, 1 for duration?)
        stop("Incompatible input format of hiatuses. Please check help page")
      }
      hiatheight=unlist(sapply(hiatuslist,function(x) head(x,1))) #get stratigraphic height of all hiatuses
      if ( max(hiatheight)>=max(xdep)|min(hiatheight)<=min(xdep)) { #remove hiatuses that are defined outside the stratigraphic height
        hiatuslist=hiatuslist[hiatheight<max(xdep)&hiatheight>min(xdep)]
      }
    }
  }
  
  #preparing input data
  {
    lims=c(min(xdep),max(xdep)) #limits where both signal and deporate are defined
    #remove points for which both the deposition rate and the signal are undefined
    relevantpoints=(x>=lims[1] & x<=lims[2])
    xt=x[relevantpoints]
    #adjust deposition rate dependent on input unit
    if (unit=="sediment per time"){
      ydep=ydep^(-1)
      }
    else if (unit=="time per sediment"){
      ydep=ydep
      }
    else{
      stop("Error: Incompatible unit (either \"sediment per time\" or \"time per sediment\")")
    }  
    #determine all relevant x and y values of the deposiiton rate for later evaluation
    xvals=sort(unique(c(xdep,xt)))
    yvals=approx(xdep,ydep,xout=xvals,yright=0,yleft=0)[[2]]
  }
  #build age model by integrating over the inverse deposition rate using the trapezoidal rule
  {
    intvals = c(0,cumsum((0.5*(c(0,yvals )+c(yvals,0))*(c(xvals,0)-c(0,xvals)))[2:length(xvals)]))                             
  }
  #in case of hiatuses, modify age model
  {
    if (length(hiatuslist)>0 ){
      hiatheight=unlist(sapply(hiatuslist,function(x) head(x,1))) #get stratigraphic height of all hiatuses
      for (i in 1:length(hiatuslist)){ #insert hiatuses into age model
        #get time at which the hiatus begins
        jumpval=approx(xvals,intvals,xout=hiatheight[i],ties="ordered")[[2]] 
        #split age model: one part below the hiatus to the beginning of the hiatus, one part from the end of the hiatus to the top. the second part is shifted by the duration of the hiatus
        intvals=c(intvals[ xvals<hiatheight[i] ],jumpval,jumpval+hiatuslist[[i]][2],intvals[ xvals>hiatheight[i] ]+hiatuslist[[i]][2] ) 
        #adjust x values to match the new age model
        xvals=c(xvals[ xvals<hiatheight[i] ],hiatheight[i],hiatheight[i],xvals[ xvals>hiatheight[i]])
      }
    }
  }
  #transform stratigraphic locations into time
  {
    xtrans=approx(xvals,intvals,xout=xt,ties="ordered")[[2]] #transform stratigraphic heights into time 
    #remove values that are located right on a hiatus
    if(length(hiatuslist)>0){
      for(j in 1:length(hiatheight)){
        xtrans[xt==hiatheight[j]]=NA
      }
    }
  }
  #adjust output size to input size
  {
    xtrans=replace(rep(NA,length(relevantpoints)),relevantpoints,xtrans)
  }
  return(list(time=xtrans,height=x))
}
