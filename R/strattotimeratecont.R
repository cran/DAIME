strattotimeratecont <-
function(xdep,ydep,xsig,ysig,pos=NULL,hiatuslist=list(),unit="sediment per time"){
  if(is.null(pos)){pos=xsig}

  #determine all relevant y values of the depositon rate
  resl=approx(xdep,ydep,xout=sort(unique(c(xdep,xsig,pos))),yleft=0,yright=0)
  if (unit=="sediment per time"){
    #get integrated inverse deposition rate using the trapezoidal rule
    I0 = c(0,cumsum((0.5*(c(0,resl[[2]]^(-1))+c(resl[[2]]^(-1),0))*(c(resl[[1]],0)-c(0,resl[[1]])))[2:length(resl[[1]])]))                              
  }
  else if (unit=="time per sediment"){
    #get integrated inverse deposition rate using the trapezoidal rule
    I0 = c(0,cumsum((0.5*(c(0,resl[[2]] )+c(resl[[2]],0))*(c(resl[[1]],0)-c(0,resl[[1]])))[2:length(resl[[1]])]))                             
  }
  else{stop("error: incompatible unit (either \"sediment per time\" or \"time per sediment\")")}                         

  if (is.list(hiatuslist)==FALSE){
    stop("Incompatible input format. Please use a list for hiatii.")
  }
  if (length(hiatuslist)>0 ){ #if there are hiatii
    if ( all(sapply(hiatuslist,length)==2) ){
      for (i in 1:length(hiatuslist)){
        jumpval=approx(resl[[1]],I0,xout=hiatuslist[[i]][1],ties="ordered")[[2]] 
        I0=c(I0[ resl[[1]]<hiatuslist[[i]][1] ],jumpval,jumpval+hiatuslist[[i]][2],I0[ resl[[1]]>hiatuslist[[i]][1] ]+hiatuslist[[i]][2] ) 
        resl[[1]]=c(resl[[1]][ resl[[1]]<hiatuslist[[i]][1] ],hiatuslist[[i]][1],hiatuslist[[i]][1],resl[[1]][ resl[[1]]>hiatuslist[[i]][1] ])
      }
    }
    else{ stop("Incompatible input format of hiatuses. Please check help page")}
  }

xtrans=approx(resl[[1]],I0,xout=pos,ties='ordered',yleft=0,yright=0)[[2]]
ytrans=approx(xsig,ysig,yleft=0,yright=0,xout=pos)[[2]]/approx(xdep,ydep,yleft=0,yright=0,xout=pos)[[2]]

return(list(x=xtrans,y=replace(ytrans,is.finite(ytrans)==FALSE,0)))
}
