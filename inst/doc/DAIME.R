## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=4,
  fig.height=4 
)

## ----fig1---------------------------------------------------------------------
require(DAIME)
par(mar=c(4,4,1,0),mgp=c(2.1,0.4,0))
agemodelage=seq(0,3,length.out = 100)
agemodelheight=splinefunH(x=c(0,1,3),y=c(0,0.8,2),m=c(0.2,2,0.3))(agemodelage)
plot(agemodelage,agemodelheight,type='l',xlab='Time/Age \n ( agemodelage )',ylab='Stratigraphic Height \n ( agemodelheight )',main='Age Model',lwd=4)

## -----------------------------------------------------------------------------
require(DAIME)
my.agemodel=SeymourIslandAgeModels$C #Select one of the age models
#plot age model
plot(my.agemodel$age,my.agemodel$height,type='l',xlab='Age [ma]',ylab = 'Height [m]',xlim=rev(range(my.agemodel$age)))
#Select extinction hypotheses with intermediate influence of the Deccan trapps
temp.pat=KPgLastOccurrences$Intermediate.Deccan
temp.pat$parameters #check description of extinction hypotheses
#Plot extinction hypotheses (=last occurrences)
plot(temp.pat$age,temp.pat$val,type='l',xlab='Age [ma]',ylab='Intensity [Last occurrences per Myr]',xlim=rev(range(temp.pat$age)))

## -----------------------------------------------------------------------------
strat.pat=patterntransform(xdep=my.agemodel$age,ydep=my.agemodel$height,xpat=temp.pat$age,ypat=temp.pat$val,direction = 'time to height',depositionmodel = 'age model',patternmode = 'piecewise linear',timetype = 'age')
strat.pat$report #What was done?
#Plot last occurrences in stratigraphic height
plot(strat.pat$height,strat.pat$val,type='l',xlab='Stratigraphic Height [m]',ylab='Intensity [Last occurrences per m]')

## -----------------------------------------------------------------------------
samplingbins=SeymourIslandBins$Macellari.1984.Section.D. #Bins are loaded with DAIME package
trans.bins=pointtransform(points=samplingbins,xdep=my.agemodel$height,ydep=my.agemodel$age,direction='height to time',depositionmodel = 'age model',timetype = 'age')
trans.bins$age #ages corresponding to the limits of the sampling bins
#Note that points that are outside the deposition model are returned as NA !
abs(diff(trans.bins$age)) #Duration of bins in Myr. The bin before the K/Pg contains almost 0.15 Myr

## -----------------------------------------------------------------------------
#Sampling relevant sampling bins from Witts et al. (2016)
WittsBins=SeymourIslandBins$Witts.et.al.2016.Section.A 
WittsBins=WittsBins[WittsBins>=min(my.agemodel$height) &WittsBins<=max(my.agemodel$height)]
#No. of fossil finds of Maorites densicostatus in section A from WItts et al. (2016)
Maorites.densicostatus=c(7,7,29,1,0,1,6) 
#Fossils per meter
Foss.Rate=Maorites.densicostatus/diff(WittsBins)
#Plot fossils per meter
plot(approx(x=WittsBins,y=c(Foss.Rate,0),method='constant',xout=my.agemodel$height),type='l',xlab='Stratigraphic Height [m]',ylab = 'Fossil Abundance [foss. per m]',main='Maorites densicostatus')

## -----------------------------------------------------------------------------
temp.pat=patterntransform(xdep=my.agemodel$height,ydep=my.agemodel$age,xpat=WittsBins,ypat=Foss.Rate,depositionmodel = 'age model',direction='height to time',patternmode = 'binned',timetype = 'age')
plot(temp.pat$age,temp.pat$val,type='l',xlim = rev(range(my.agemodel$age)),xlab='Age [ma]',ylab='Fossil abundance [Fossils per Myr]', main='Maorites densicostatus')

## -----------------------------------------------------------------------------
palyn=ColDePalluel$palynomorphs
geochem=ColDePalluel$geochemistry

## -----------------------------------------------------------------------------
my.agemodel=patterntodepositionmodel(xheight=palyn$Depth..m.,yheight = palyn$Pollen.sac....,timetype = 'time')
plot(my.agemodel$time,my.agemodel$height,type='l',xlab='Relative Time',ylab='Stratigraphic Height [m]')
lines(range(my.agemodel$time),range(my.agemodel$height),lty=3)
legend('topleft',legend='Age Model',lty=1,lwd=1)

## -----------------------------------------------------------------------------
#Transform dinoflagellate cysts into relative time
temp.pat=patterntransform(xdep=my.agemodel$height,my.agemodel$time,xpat=palyn$Depth..m.,ypat=palyn$Dinofl.cyst....,direction='height to time',timetype = 'time',depositionmodel = 'age model')
plot(palyn$Depth..m.,palyn$Dinofl.cyst....,type='l',xlab='Stratigraphic Height [m]',ylab='Dinoflagellate cysts [counts]')
plot(temp.pat$time,temp.pat$val,xlab='Relative Time',ylab='Dinoflagellate cysts [counts]',type='l')

## -----------------------------------------------------------------------------
#Time of spike
spike.time=temp.pat$time[which(temp.pat$val==max(temp.pat$val))]
#Transform time of spike into height
poi=pointtransform(points=spike.time,xdep=my.agemodel$time,ydep=my.agemodel$height,direction='time to height',depositionmodel = 'age model')
#stratigraphic height of spike
poi$height
#get lithology closest to the spike
geochem$Lithology[which(abs(geochem$Depth..m.-poi$height)==min(abs(geochem$Depth..m.-poi$height)))]

## -----------------------------------------------------------------------------
#get times of deposition for geochemnistry data
sampletimes.geochem=pointtransform(points=geochem$Depth..m.,xdep=my.agemodel$height,my.agemodel$time,direction='height to time',depositionmodel = 'age model',timetype = 'time')

#transform delta13 C into time
plot(sampletimes.geochem$height, geochem$delta13C.carb....PDB.,type='l',xlab='Stratigraphic Height [m]', ylab = "delts13C [per mil]")
plot(sampletimes.geochem$time, geochem$delta13C.carb....PDB.,type='l',xlab = 'Relative Time',ylab='delta13C [per mil]')


