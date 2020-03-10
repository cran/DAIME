## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----fig.show='hold', fig.height=3.2, fig.width=3.2, echo=FALSE---------------
agemodelage=seq(0,3,length.out = 100)
agemodelheight=splinefunH(x=c(0,1,3),y=c(0,0.8,2),m=c(0.2,2,0.3))(agemodelage)
plot(agemodelage,agemodelheight,
     type='l',
     xlab='Time',
     ylab='Stratigraphic Height',
     main='Age Model',
     lwd=4,mgp=c(2.2,0.5,0))
plot(agemodelage,agemodelheight,
     type='l',
     xlab='Age',
     ylab='Stratigraphic Height',
     main='Age Model',
     lwd=4,mgp=c(2.2,0.5,0),xaxt='n')
axis(side=1,labels=3:0,at=0:3)


## ---- echo=FALSE--------------------------------------------------------------
xdep=c(0,1,2)
ydep=c(1,0.2,2)
plot(xdep,ydep,
     xlab='Time',ylab = 'Sediment deposited per time unit',
     main='Deposition Rate',
     type='l',ylim=c(0,2))

## ---- echo=FALSE--------------------------------------------------------------
xpat=seq(0,1,length.out = 100)
ypat=splinefunH(x=c(0,0.4,1),y=c(0.2,0.8,0.5),m=c(0.2,2,0.3))(xpat)
plot(xpat,ypat,type='l',xlab='',ylab='',main='Stratigraphic Pattern',lwd=4,ylim=c(0,2))
mtext(side=2,text='Fossils per Stratigraphic Height',line=2)
mtext(side=1,text='Stratigraphic Height',line=3)

## ----fig.show='hold', fig.height=3.2, fig.width=3.2, echo=FALSE---------------
require(DAIME)
xdep=seq(0,1,length.out = 100)
ydep=splinefunH(x=c(0,0.5,1),y=c(0.1,0.8,2),m=c(0.2,2,0.3))(xdep)
plot(xdep,ydep,
     xlab='Stratigraphic Height',
     ylab='Deposition Rate \n [sediment per unit of time]',
     main='Deposition Rate',
     lwd=4,ylim=c(0,2),type='l',mgp=c(2.2,0.5,0))
temp.pat=patterntransform(xdep,ydep,
                          xpat,ypat,
                          direction = 'height to time',
                          unit='sediment per time')
plot(temp.pat$time,temp.pat$val,
     xlab='Time',ylab='Fossils per unit of time',
     main='Temporal Pattern',
     type='l',lwd=4,mgp=c(1.5,0.5,0))


## ----fig.show='hold', fig.height=3.2, fig.width=3.2, echo=FALSE---------------
xdep=seq(0,1,length.out = 100)
ydep=splinefunH(x=c(0,0.5,1),y=c(0.1,0.8,2),m=c(0.2,2,0.3))(xdep)
plot(xdep,ydep,
     xlab='Stratigraphic Height',
     ylab='Inverse Deposition Rate\n[time per unit of sediment]',
     main='Inverse\nDeposition Rate',
     lwd=4,ylim=c(0,2),type='l',mgp=c(2.2,0.5,0))
temp.pat=patterntransform(xdep,ydep,
                          xpat,ypat,
                          direction = 'height to time',
                          unit='time per sediment')
plot(temp.pat$time,temp.pat$val,
     xlab='Time',ylab='Fossils per unit of time',
     main='Temporal Pattern',
     type='l',lwd=4,mgp=c(1.5,0.5,0),
     ylim=c(0,max(temp.pat$val)))


## ----  fig.height=5, fig.width=5----------------------------------------------
require(DAIME)
 #Check help page for data sources and information on the age models:
?SeymourIslandAgeModels

my.agemodel=SeymourIslandAgeModels$C #Select one of the age models
#plot age model
plot(my.agemodel$age,my.agemodel$height,
     xlab='Age [Ma]',ylab = 'Height [m]',
     xlim=rev(range(my.agemodel$age)),type='l')
## Get extinction hypotheses
#Check help page for data sources and information on the last occurrences:
?KPgLastOccurrences

#Select extinction hypotheses with intermediate influence of the Deccan trapps
temp.pat=KPgLastOccurrences$Intermediate.Deccan
temp.pat$parameters #check description of extinction hypotheses
#Plot extinction hypotheses (=last occurrences)
plot(temp.pat$age,temp.pat$val,
     xlab='Age [Ma]',ylab='Intensity [Last occurrences per Myr]',
     xlim=rev(range(temp.pat$age)),type='l')

## ----  fig.height=5, fig.width=5----------------------------------------------
strat.pat=patterntransform(xdep=my.agemodel$age,ydep=my.agemodel$height,
                           xpat=temp.pat$age,ypat=temp.pat$val,
                           direction = 'time to height',
                           depositionmodel = 'age model',
                           patternmode = 'piecewise linear',
                           timetype = 'age')
strat.pat$report #What was done?
#Plot last occurrences in stratigraphic height
plot(strat.pat$height,strat.pat$val,
     xlab='Stratigraphic height [m]',ylab='Intensity [Last occurrences per m]',
     type='l')

## ----  fig.height=5, fig.width=5----------------------------------------------
##Sampling bins are automatically loaded with the DAIME package
#Check help page for data sources and information on the sampling bins:
?SeymourIslandBins

samplingbins=SeymourIslandBins$Macellari.1984.Section.D. 
trans.bins=pointtransform(points=samplingbins,
                          xdep=my.agemodel$height,ydep=my.agemodel$age,
                          direction='height to time',
                          depositionmodel = 'age model',
                          timetype = 'age')
trans.bins$age #ages corresponding to the limits of the sampling bins
#Note that points that are outside the deposition model are returned as NA !
abs(diff(trans.bins$age)) #Duration of bins in Myr.
#The bin before the K/Pg contains almost 0.15 Myr

## ----  fig.height=5, fig.width=5----------------------------------------------
#Sampling relevant sampling bins from Witts et al. (2016)
WittsBins=SeymourIslandBins$Witts.et.al.2016.Section.A 
WittsBins=WittsBins[WittsBins>=min(my.agemodel$height) &WittsBins<=max(my.agemodel$height)]
#No. of fossil finds of Maorites densicostatus in section A from WItts et al. (2016)
Maorites.densicostatus=c(7,7,29,1,0,1,6) 
#Fossils per meter
Foss.Rate=Maorites.densicostatus/diff(WittsBins)
#Plot fossils per meter
plot(approx(x=WittsBins,y=c(Foss.Rate,0),method='constant',xout=my.agemodel$height),
     xlab='Stratigraphic height [m]',
     ylab = 'Fossil abundance [Fossils per m]',main='Maorites densicostatus',
     type='l')

## ----  fig.height=5, fig.width=5----------------------------------------------
temp.pat=patterntransform(xdep=my.agemodel$height,ydep=my.agemodel$age,
                          xpat=WittsBins,ypat=Foss.Rate,
                          depositionmodel = 'age model',direction='height to time',
                          patternmode = 'binned',
                          timetype = 'age')
plot(temp.pat$age,temp.pat$val,
     xlab='Age [ma]',ylab='Fossil abundance [Fossils per Myr]',
     main='Maorites densicostatus', xlim = rev(range(my.agemodel$age)),
     type='l')

## -----------------------------------------------------------------------------
#Check help page for data sources and information on the data:
?ColDePalluel
#Use the following abbreviations
palyn=ColDePalluel$palynomorphs
geochem=ColDePalluel$geochemistry

## ----  fig.height=5, fig.width=5----------------------------------------------
my.agemodel=patterntodepositionmodel(xheight=palyn$Depth..m.,
                                     yheight = palyn$Pollen.sac....,
                                     timetype = 'time')
plot(my.agemodel$time,my.agemodel$height,
     xlab='Relative Time',ylab='Stratigraphic Height [m]',
     type='l')
lines(range(my.agemodel$time),range(my.agemodel$height),lty=3)
legend('topleft',legend='Age Model',lty=1,lwd=1)

## ----  fig.height=5, fig.width=5----------------------------------------------
#Transform dinoflagellate cysts into relative time
temp.pat=patterntransform(xdep=my.agemodel$height,ydep=my.agemodel$time,
                          xpat=palyn$Depth..m.,ypat=palyn$Dinofl.cyst....,
                          direction='height to time',timetype = 'time',
                          depositionmodel = 'age model')
plot(palyn$Depth..m.,palyn$Dinofl.cyst....,
     xlab='Stratigraphic Height [m]',ylab='Dinoflagellate cysts [counts]',
     type='l')
plot(temp.pat$time,temp.pat$val,
     xlab='Relative Time', ylab='Dinoflagellate cysts [counts]',
     type='l')

## ----  fig.height=5, fig.width=5----------------------------------------------
#Time of spike
spike.time=temp.pat$time[which(temp.pat$val==max(temp.pat$val))]
#Transform time of spike into height
poi=pointtransform(points=spike.time,
                   xdep=my.agemodel$time,ydep=my.agemodel$height,
                   direction='time to height',
                   depositionmodel = 'age model')
#stratigraphic height of spike
poi$height
#get lithology closest to the spike
geochem$Lithology[which(abs(geochem$Depth..m.-poi$height)==min(abs(geochem$Depth..m.-poi$height)))]

## ----  fig.height=5, fig.width=5----------------------------------------------
#get times of deposition for geochemnistry data
sampletimes.geochem=pointtransform(points=geochem$Depth..m.,xdep=my.agemodel$height,
                                   ydep=my.agemodel$time,
                                   direction='height to time',
                                   depositionmodel = 'age model',
                                   timetype = 'time')

#transform delta13 C into time
plot(sampletimes.geochem$height, geochem$delta13C.carb....PDB.,
     xlab='Stratigraphic Height [m]', ylab = "delta13C [per mil]",
     type='l')
plot(sampletimes.geochem$time, geochem$delta13C.carb....PDB.,
     xlab = 'Relative Time',ylab='delta13C [per mil]',
     type='l')


## ----  fig.height=5, fig.width=5----------------------------------------------
#Check help page for data sources and information on the data:
?LakeSuperior
#The stratigraphic pattern is given by the measured excess lead
plot(LakeSuperior$BH03_3$Depth..m.,LakeSuperior$BH03_3$X210Pb.xs..Bq.kg.,
     xlab='Core depth [m]',ylab='Excess lead [Bq/kg]',
     main='Stratigraphic Pattern')
#The temporal pattern is excess lead content in the sediment
#this will follow an exponential decay
dec.const=log(2)/22.3 #decay constant of 210Pb
temp.pat.x=seq(0,300,by=0.1)
temp.pat.y=exp(-dec.const*temp.pat.x)
plot(temp.pat.x,temp.pat.y,
     xlab = 'Sediment Age [years]',
     ylab='Proportion of initial Excess lead remaining in sediment',
     main='Temporal pattern',
     type='l')

## ----  fig.height=5, fig.width=5----------------------------------------------
agemodel=patterntodepositionmodel(xheight=LakeSuperior$BH03_3$Depth..m.,
                                  yheight = LakeSuperior$BH03_3$X210Pb.xs..Bq.kg.,
                                  xage = temp.pat.x,yage=temp.pat.y,
                                  timetype = 'time',
                                  rescalefor = 'stratigraphic pattern')
#note that "timetype" is set to 'time', and not 'age'.
#This is since although we are interested in age, it is examined from youngest to oldest (in reverse order)
plot(agemodel$time,agemodel$height,
     xlab='Years before present',
     ylab='Core depth [m]',
     main='Age Model',
     lwd=3,type='l')
lines(x=1000*LakeSuperior$BH03_3$Age.model..ka., #Age model from O\'Beirne et al. (2017) is given in kyr
      y=LakeSuperior$BH03_3$Depth..m.,
      lty=3)
legend('bottomright',
       lwd=c(3,1),lty=c(1,3),
       legend=c('DAIME','O\'Beirne et al. (2017)'))

## ----  fig.height=5, fig.width=5----------------------------------------------
#Assume a changing input of 210Pb:
Pb.input=splinefunH(x=c(0,89,300),y=c(1,5,1),m=c(0.1,-0.1,-0.01))(temp.pat.x)
plot(temp.pat.x,Pb.input,
     xlab='Years before present', ylab='210Pb input',
     ylim=c(0,6),type='l')
#multiplying the input with the curve from the radioactive decay yields the temporal pattern:
temp.pat.y2=Pb.input*temp.pat.y
plot(temp.pat.x,temp.pat.y2,
     xlab='Sediment Age',ylab='210 Pb content of Sediment',
     main='Temporal pattern',
     type='l')
#This pattern can now again be used to create an age model
agemodel=patterntodepositionmodel(xheight=LakeSuperior$BH03_3$Depth..m.,
                                  yheight = LakeSuperior$BH03_3$X210Pb.xs..Bq.kg.,
                                  xage = temp.pat.x,
                                  yage=temp.pat.y2,
                                  timetype = 'time',
                                  rescalefor = 'stratigraphic pattern')
#
plot(agemodel$time,agemodel$height,
     xlab='Years before present', ylab='Core depth [m]',
     main='Age Model',
     lwd=3, xlim=c(0,150),type='l')

