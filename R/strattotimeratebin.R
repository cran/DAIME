strattotimeratebin <-
function(binborder,depoval,signalval,pos=NULL,hiatuslist=list(),unit="sediment per time"){
  #check input
  {
    stopifnot(is.list(hiatuslist),is.unsorted(binborder,strictly=TRUE)==FALSE,min(depoval)>0,min(signalval)>=0,(length(binborder)-1)==length(depoval),(length(binborder)-1)==length(signalval))
    #check input format of hiatuses
    if (length(hiatuslist)>0 ){
      if ( any(sapply(hiatuslist,length)!=2) ){ #do all entries of the list have 2 components (1 for strat. height, 1 for duration?)
        stop("Incompatible input format of hiatuses. Please check help page")
      }
      hiatheight=unlist(sapply(hiatuslist,function(x) head(x,1))) #get stratigraphic height of all hiatuses
      if ( max(hiatheight)>=max(binborder)|min(hiatheight)<=min(binborder)) { #remove hiatuses that are defined outside the stratigraphic height
        hiatuslist=hiatuslist[hiatheight<max(binborder)&hiatheight>min(binborder)]
      }
    }
  }
  
  #preparing input data
  {
    lims=c(min(binborder),max(binborder)) #limits where both signal and deporate are defined
    if(is.null(pos)){ #default setting for pos: if no points for evaluation are given
      pos=seq(from=lims[1],to=lims[2],length.out=min(1000,20*length(binborder)))
    }
    #remove points for which the deposition rate is undefined
    relevantpoints=(pos>=lims[1] & pos<=lims[2])
    xt=pos[relevantpoints]
    #adjust deposition rate dependent on input unit
    if (unit=="sediment per time"){
      depoval=depoval^(-1)
    }
    else if (unit=="time per sediment"){
      depoval=depoval
    }
    else{
      stop("error: incompatible unit (either \"sediment per time\" or \"time per sediment\")")
    }
    xvals=binborder
  }
  
  #build age model by integrating over the inverse deposition rate
  {
    intvals=cumsum(c(0,diff(xvals)*depoval))
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
  #transform stratigraphic rate into time
  {
    xtrans=approx(xvals,intvals,xout=xt,ties="ordered")[[2]] #transform stratigraphic heights into time 
    ytrans=approx(binborder,c(signalval,tail(signalval,1)),method='constant',xout=xt)[[2]]/approx(binborder,c(depoval,tail(depoval,1)),method='constant',xout=xt)[[2]] #values of the temporal rate at the points in time
    #remove values that are located right on a hiatus
    if(length(hiatuslist)>0){
      for(j in 1:length(hiatheight)){
        xtrans[xt==hiatheight[j]]=NA
        ytrans[xt==hiatheight[j]]=NA
      }
    }
  }
  #adjust output to input size
  {
    xtrans=replace(rep(NA,length(relevantpoints)),relevantpoints,xtrans)
    ytrans=replace(rep(NA,length(relevantpoints)),relevantpoints,ytrans)
  }
return(list(time=xtrans,val=ytrans))
}
