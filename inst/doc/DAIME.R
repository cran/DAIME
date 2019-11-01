## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=4,
  fig.height=4 
)

## ----fig1----------------------------------------------------------------

par(mar=c(4,4,1,0),mgp=c(2.1,0.4,0))
agemodelage=seq(0,3,length.out = 100)
agemodelheight=splinefunH(x=c(0,1,3),y=c(0,0.8,2),m=c(0.2,2,0.3))(agemodelage)
plot(agemodelage,agemodelheight,type='l',xlab='Time/Age \n ( agemodelage )',ylab='Stratigraphic Height \n ( agemodelheight )',main='Age Model',lwd=4)

