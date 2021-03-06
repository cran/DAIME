\name{patterntransform}
\alias{patterntransform}

\title{
Transform patterns between stratigraphic height and time/age
}
\description{
This function (1) takes temporal patterns and determines the corresponding stratigraphic patterns or (2) takes stratigraphic patterns and reconstructs the underlying temporal patterns.

These patterns can for example be rates of (first/last) fossil occurrences, rates of morphological change, or input rates into the sediment, e.g. of geochemical proxies.
}
\usage{
patterntransform ( xdep , ydep , xpat , ypat ,  direction = 'time to height' ,
depositionmodel = 'piecewise linear deposition rate' ,  patternmode = 'piecewise linear' ,
pos=NULL , hiatuslist=list() ,unit = 'time per sediment',timetype='time' )
}

\arguments{
\item{xdep}{
Vector of strictly increasing real numbers
}
  \item{ydep}{
Vector of real numbers. \code{xdep} and \code{ydep} describe the deposition model, and how they are interpreted depends on the option \code{depositionmodel} (and when \code{direction='height to time'}, also on the option \code{unit}). By default, \code{xdep} and \code{ydep} describe a piecewise linear deposition rate, i.e. \code{deporate=approxfun(xdep,ydep)}. For more options on how to describe deposition models, see \code{depositionmodel} below
}
\item{xpat}{
Vector of strictly increasing real numbers
}
\item{ypat}{
Vector containing positive real numbers. \code{xpat} and \code{ypat} describe the pattern to be transformed. By default, it is taken to be piecewise linear and is accordingly given by \code{pattern=approxfun(xpat,ypat)}. The option \code{patternmode} allows to change to binned patterns (see below). Whether the input pattern is a temporal or a stratigraphic pattern is determined by \code{direction}
}
  \item{direction}{OPTIONAL, default is \code{'time to height'}. Determines the direction of the transformation. Either \code{'time to height'} (in which case the input pattern described by \code{xpat} and \code{ypat} is interpreted as a temporal pattern and will be transformed into a stratigraphic pattern) or \code{'height to time'} (in which case the input pattern described by \code{xpat} and \code{ypat} is interpreted as a stratigraphic pattern and will be transformed into a temporal pattern) 
}
  \item{depositionmodel}{OPTIONAL, default is \code{'piecewise linear deposition rate'}. Either \code{'piecewise linear deposition rate'}, \code{'binned deposition rate'}, or \code{'age model'}. Determines how the input of \code{xdep} and \code{ydep} is interpreted. Dependent on the choice of \code{direction} and \code{depositionmodel}, different restrictions on the input of \code{ydep} exist, for details and examples see the vignette (available via \code{vignette('DAIME') or the error messages}
}
  \item{patternmode}{OPTIONAL, default is \code{'piecewise linear'}. Either \code{'piecewise linear'} or \code{'binned'}. Determines whether \code{xpat} and \code{ypat} are taken as a piecewise linear or a binned description of the pattern. In the first case, \code{length(xpat)} needs to match \code{length(ypat)}, and the pattern is given by \code{approxfun(xpat,ypat)}. In the second case, \code{length(xpat)} needs to match \code{length(ypat)+1}, and the value of the pattern between the \code{xpat[i]} and \code{xpat[i+1]} is given by \code{ypat[i]}
}
  \item{pos}{OPTIONAL, default is \code{NULL}. A vector of points that will be transformed and then used determine the value of the pattern at these points. In the default setting, the number and location of points is chosen automatically
}
  \item{hiatuslist}{
OPTIONAL, default is an empty list. List of hiatuses to be included into the transformation, only used when \code{direction='height to time'}, e.g. the transformation is from stratigraphic height to time. A list, containing vectors of length two as elements. Every element in the list corresponds to a hiatus: \code{hiatuslist[[i]][1]} is the stratigraphic height of the i-th hiatus, and \code{hiatuslist[[i]][2]} is its duration
}
  \item{unit}{
OPTIONAL, default is 'sediment per time'. Only used when \code{direction='height to time'}, e.g. the transformation is from stratigraphic height to time, and \code{deposition model} is \code{'piecewise linear deposition rate'} or \code{'binned deposition rate'}. Either \code{'sediment per time'} or \code{'time per sediment'}. Determines the unit of the (binned or piecewise linear) deposition rate in the sense that if \code{unit='sediment per time'}, the deposition rate describes the deposition rate with which a given point in the section was deposited. If \code{unit='time per sediment'}, the deposition rate describes the time needed to deposit one unit of sediment
}
\item{timetype}{OPTIONAL, default is "time". Either "time" or "age", determines whether input/output will be interpreted/given as time or age
}
}

\value{
Returns a list with two elements. If \code{direction='time to height'}, these elements are
\item{height }{vector containing stratigraphic heights}
\item{val }{vector containing the values of the stratigraphic pattern (=transformed temporal pattern) at the stratigraphic heights given by \code{height}}
\item{report }{A short summary of the task performed}

If \code{direction='height to time'} and \code{timetype='time'}, these elements are
\item{time }{vector containing points in time}
\item{val }{vector containing the values of the temporal pattern (=transformed stratigraphic pattern) at the points in time given by \code{time}}
\item{report }{A short summary of the task performed}

If \code{direction='height to time'} and \code{timetype='age'}, these elements are
\item{age }{vector containing ages}
\item{val }{vector containing the values of the temporal pattern (=transformed stratigraphic pattern) at the ages given by \code{age}}
\item{report }{A short summary of the task performed}

Output of \code{NA} in \code{age},\code{time}, \code{height}, or \code{val} indicates that some values coincide with a hiatus or intervals where the deposition model is undefined.
}

\references{
Hohmann, Niklas. 2018. Quantifying the Effects of Changing Deposition Rates and Hiatii on the Stratigraphic Distribution of Fossils. <doi:10.13140/RG.2.2.23372.51841>
}

\author{
Niklas Hohmann
}

\seealso{
For an overview of the functions in the DAIME package and examples using stratigraphic data see the vignette (available via \code{vignette('DAIME')} )

\code{\link{pointtransform}} for the transformation of points

\code{\link{patterntodepositionmodel}} to create deposition models based in sedimentary dilution/condensation of known patterns

\code{\link{strattotimeratebin}}, \code{\link{strattotimeratecont}}, \code{\link{timetostratratebin}}, and  \code{\link{timetostratratecont}} for diverse (outdated) wrappers of \code{patterntransform}. 
}

\examples{
### Example 1: transform patterns from time into the section
##define (piecewise linear) deposition rate
my.xdep=seq(0,12,length.out=100)
my.ydep=splinefunH(x=c(0,2,4,6,8,10,12),y=c(1,5,6,1,0.5,1,6),m=c(0,1.5,-0.5,-0.5,0,0.5,0))(my.xdep)
#Plot deposition rate
plot(my.xdep,my.ydep,type='l',main='Deposition Rate',xlab='Time', ylab='Sediment per Time Unit',
  ylim=c(0,max(my.ydep))) 
##define (piecewise linear) temporal pattern, here fossil occurrences
patternname='Fossil Occurrences'
my.temp.pat.x= seq(from=min(my.xdep),to=max(my.xdep),length.out=100)
my.temp.pat.y=splinefunH(x=c(0,4,12),y=c(0.5,2,0.5),m=c(0,0,0))(my.temp.pat.x)  

plot(my.temp.pat.x,my.temp.pat.y,type='l',main=paste(patternname,'in time'),
  xlab='Time',ylab=patternname,ylim=c(0,max(my.temp.pat.y)))
##transform temporal pattern into stratigraphic pattern
strat.pat=patterntransform(xdep=my.xdep,ydep=my.ydep,xpat=my.temp.pat.x,ypat=my.temp.pat.y,
  direction = 'time to height' , depositionmodel = 'piecewise linear deposition rate' ,
  patternmode = 'piecewise linear')
strat.pat$report
my.strat.pat.x=strat.pat$height
my.strat.pat.y=strat.pat$val
plot(my.strat.pat.x,my.strat.pat.y,type='l',xlab='Stratigraphic Height',
ylab=patternname, main=paste(patternname,' in the section'))

## Variation: modify deposition rate to include erosion
my.xdep2=my.xdep
my.ydep2=splinefunH(x=c(0,2,4,6,8,10,12),y=c(1,5,6,-2,5,1,6),m=c(0,1.5,-0.5,-0.5,0,0.5,0))(my.xdep2)
#Plot deposition rate
plot(my.xdep2,my.ydep2,type='l',main='Deposition Rate',xlab='Time',ylab='Sediment per Time Unit')
lines(range(my.xdep2),c(0,0))
strat.pat2=patterntransform(xdep=my.xdep2,ydep=my.ydep2,xpat=my.temp.pat.x,ypat=my.temp.pat.y,
  direction = 'time to height' , depositionmodel = 'piecewise linear deposition rate' ,
  patternmode = 'piecewise linear')
strat.pat2$report
my.strat.pat.x2=strat.pat2$height
my.strat.pat.y2=strat.pat2$val
plot(my.strat.pat.x2,my.strat.pat.y2,type='l',xlab='Stratigraphic Height',
  ylab=patternname,main=paste(patternname,' in the section'),ylim=c(0,1))
#the spike is because the deposition rate is very small when it transitions from 
#negative to positive, generating a punctual extreme condensation

### Example 2: Transform patterns from the section into time
patternname='Fossil Occurrences'
#piecewise linear deposition rate
my.xdep=seq(0,12,length.out=100)
my.ydep=splinefunH(x=c(0,2,4,6,8,10,12),y=c(1,5,6,1,0.5,1,6),m=c(0,1.5,-0.5,-0.5,0,0.5,0))(my.xdep)
## Same pattern, this time interpreted as stratigraphic pattern
my.strat.pat.x= seq(from=min(my.xdep),to=max(my.xdep),length.out=100)
my.strat.pat.y=splinefunH(x=c(0,4,12),y=c(0.5,2,0.5),m=c(0,0,0))(my.strat.pat.x)
plot(my.strat.pat.x,my.strat.pat.y,type='l',main=paste(patternname,'in the section'),
  xlab='Stratigraphic Height',ylab=patternname,ylim=c(0,max(my.strat.pat.y)))

##interpret deposition rate as inverse deposition rate with the unit "time per sediment"
usedunit='time per sediment'
plot(my.xdep,my.ydep,type='l',main=paste('Deposition Rate as',usedunit),
  xlab='Stratigraphic Height',ylab=usedunit,ylim=c(0,max(my.ydep)))
#transform pattern
temp.pat1=patterntransform(xdep=my.xdep,ydep=my.ydep,xpat=my.strat.pat.x,ypat=my.strat.pat.y,
  direction = 'height to time' , depositionmodel = 'piecewise linear deposition rate' ,
  patternmode = 'piecewise linear', unit=usedunit)
temp.pat1$report
#plot result
plot(temp.pat1$time,temp.pat1$val,type='l',xlab='Time',
  ylab=patternname,main=paste(patternname, 'in time'),ylim=c(0,max(temp.pat1$val)))
  
##take the deposition rate as ordinary deposition rate with unit "sediment per time"
usedunit='sediment per time'
plot(my.xdep,my.ydep,type='l',main=paste('Deposition Rate as',usedunit),
xlab='Stratigraphic Height',ylab=usedunit,ylim=c(0,max(my.ydep)))
#transform pattern
temp.pat2=patterntransform(xdep=my.xdep,ydep=my.ydep,xpat=my.strat.pat.x,ypat=my.strat.pat.y,
  direction = 'height to time' ,depositionmodel = 'piecewise linear deposition rate' ,
  patternmode = 'piecewise linear', unit=usedunit)
temp.pat2$report
#plot result
plot(temp.pat2$time,temp.pat2$val,type='l',xlab='Time',
  ylab=patternname,main=paste(patternname, 'in time'),ylim=c(0,max(temp.pat2$val)))
#Note the big difference in the resulting patterns in time depending on the unit used.

}

\keyword{ Deposition rate }
\keyword{ Isotope }