\name{patterntransform}
\alias{patterntransform}

\title{
Transform patterns between stratigraphic height and time
}
\description{
This function (1) takes temporal patterns and determines the corresponding stratigraphic patterns or (2) takes stratigraphic patterns and reconstructs the underlying temporal patterns.

These patterns can for example be rates of (first/last) fossil occurrences, rates of morphological change, or input rates into the sediment, e.g. of geochemical proxies.
}
\usage{
patterntransform ( xdep , ydep , xpat , ypat ,  direction = 'time to height' ,
depositionmodel = 'piecewise linear deposition rate' ,  patternmode = 'piecewise linear' ,
pos=NULL , hiatuslist=list() ,unit = 'time per sediment' )
}

\arguments{
\item{xdep}{
Vector of strictly increasing real numbers
}
  \item{ydep}{
Vector of real numbers. \code{xdep} and \code{ydep} describe the deposition model, and how they are interpreted depends on the option \code{depositionmodel} (and when \code{direction='height to time'}, also on the option \code{unit}). By default, \code{xdep} and \code{ydep} describe a piecewise linear deposition rate, i.e. \code{deporate=approxfun(xdep,ydep)}. For more options on how to describe deposition models, see the input \code{depositionmodel} below
}
\item{xpat}{
Vector containing strictly increasing real numbers
}
\item{ypat}{
Vector containing positive real numbers. \code{xpat} and \code{ypat} describe the pattern to be transformed. By default, it is taken to be piecewise linear and is accordingly given by \code{pattern=approxfun(xpat,ypat)}. The option \code{patternmode} allows to change to binned patterns (see below). Whether the input pattern is a temporal or a stratigraphic pattern is determined by \code{direction}
}
  \item{direction}{OPTIONAL, default is \code{'time to height'}. Determines the direction of the transformation. Either \code{'time to height'} (in which case the input pattern described by \code{xpat} and \code{ypat} is interpreted as a temporal pattern and will be transformed into a stratigraphic pattern) or \code{'height to time'} (in which case the input pattern is described by \code{xpat} and \code{ypat} is interpreted as a stratigraphic pattern and will be transformed into a temporal pattern) 
}
  \item{depositionmodel}{OPTIONAL, default is \code{'piecewise linear deposition rate'}. Either \code{'piecewise linear deposition rate'}, \code{'binned deposition rate'}, or \code{'age model'}. Determines how the input of \code{xdep} and \code{ydep} is interpreted. Dependent on the choice of \code{direction} and \code{depositionmodel}, different restrictions on the input of \code{ydep} exist, for details and examples see the vignette (available via \code{vignette('DAIME')}
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
}

\value{
Returns a list with two elements. If \code{direction='time to height'}, these elements are
\item{height }{vector containing stratigraphic heights}
\item{val }{vector containing the values of the stratigraphic pattern (=transformed temporal pattern) at the stratigraphic heights given by \code{height}}

If \code{direction='height to time'}, these elements are
\item{age }{vector containing points in time}
\item{val }{vector containing the values of the temporal pattern (=transformed stratigraphic pattern) at the points in time given by \code{age}}

Output of \code{NA} in \code{age}, \code{height}, or \code{val} indicates that some values coincide with a hiatus or intervals where the deposition rate is undefined.
}

\references{
Hohmann, Niklas. 2018. Quantifying the Effects of Changing Deposition Rates and Hiatii on the Stratigraphic Distribution of Fossils. <doi:10.13140/RG.2.2.23372.51841>
}

\author{
Niklas Hohmann
}

\seealso{
For an overview of the functions in the DAIME package, see its vignette (available via \code{vignette('DAIME')} )

\code{\link{pointtransform}} for the transformation of points

\code{\link{patterntodepositionmodel}} to create deposition models based in sedimentary dilution/condensation of known patterns

\code{\link{strattotimeratebin}}, \code{\link{strattotimeratecont}}, \code{\link{timetostratratebin}}, and  \code{\link{timetostratratecont}} for diverse wrappers of \code{patterntransform}. 
}

\examples{
##define (piecewise linear) deposition rate
xdep=seq(0,12,length.out=100)
ydep=splinefunH(x=c(0,2,4,6,8,10,12),y=c(1,5,6,1,0.5,1,6),m=c(0,1.5,-0.5,-0.5,0,0.5,0))(xdep)
#Plot deposition rate
plot(xdep,ydep,type='l',main='Deposition Rate',xlab='Time', ylab='Sediment per Time Unit') 
##define (piecewise linear) pattern, here fossil occurrences
patternname='Fossil Occurrences'
xpat= seq(from=min(xdep),to=max(xdep),length.out=100)
ypat=splinefunH(x=c(0,4,12),y=c(0.5,2,0.5),m=c(0,0,0))(xpat) #function values of the signal
plot(xpat,ypat,type='l',main=paste(patternname,'in time'),
  xlab='Time',ylab=patternname)

### Example 1: transform patterns from time into the section
reslist=patterntransform(xdep,ydep,xpat,ypat,direction = 'time to height' , 
depositionmodel = 'piecewise linear deposition rate' , patternmode = 'piecewise linear')
plot(reslist$height,reslist$val,type='l',xlab='Stratigraphic Height',
ylab=patternname, main=paste(patternname,' in the section'))

## modify deposition rate to include erosion
ydep=splinefunH(x=c(0,2,4,6,8,10,12),y=c(1,5,6,-2,5,1,6),m=c(0,1.5,-0.5,-0.5,0,0.5,0))(xdep)
#Plot deposition rate
plot(xdep,ydep,type='l',main='Deposition Rate',xlab='Time',ylab='Sediment per Time Unit')
lines(c(0,12),c(0,0))
reslist=patterntransform(xdep,ydep,xpat,ypat,direction = 'time to height' , 
depositionmodel = 'piecewise linear deposition rate' ,  patternmode = 'piecewise linear')
plot(reslist$height,reslist$val,type='l',xlab='Stratigraphic Height',
  ylab=patternname,main=paste(patternname,' in the section'),ylim=c(0,1))
#the spike is because the deposition rate is very small when it transitions from 
#negative to positive, generating a punctual extreme condensation

### Example 2: Transform patterns from the section into time
## Same pattern, this time interpreted as stratigraphic pattern
xpat= seq(from=min(xdep),to=max(xdep),length.out=100)
ypat=splinefunH(x=c(0,4,12),y=c(0.5,2,0.5),m=c(0,0,0))(xpat) #function values of the signal
plot(xpat,ypat,type='l',main=paste(patternname,'in the section'),
  xlab='Stratigraphic Height',ylab=patternname)
#piecewise linear deposition rate
xdep=seq(0,12,length.out=100)
ydep=splinefunH(x=c(0,2,4,6,8,10,12),y=c(1,5,6,1,0.5,1,6),m=c(0,1.5,-0.5,-0.5,0,0.5,0))(xdep)
##unit of the deposition rate
usedunit='time per sediment'
plot(xdep,ydep,type='l',main=paste('Deposition Rate as',usedunit),xlab='Stratigraphic Height',
  ylab=usedunit)

##transform
reslist=patterntransform(xdep,ydep,xpat,ypat,direction = 'height to time' , 
depositionmodel = 'piecewise linear deposition rate' , 
patternmode = 'piecewise linear', unit=usedunit)
##plot result
plot(reslist$age,reslist$val,type='l',xlab='Time',
  ylab=patternname,main=paste(patternname, 'in time'))
  
##use different unit for transformation
usedunit='sediment per time'
plot(xdep,ydep,type='l',main=paste('Deposition Rate as',usedunit),
xlab='Stratigraphic Height',ylab=usedunit)
##transform
reslist=patterntransform(xdep,ydep,xpat,ypat,direction = 'height to time' ,
depositionmodel = 'piecewise linear deposition rate' ,
patternmode = 'piecewise linear', unit=usedunit)
##plot result
plot(reslist$age,reslist$val,type='l',xlab='Time',
  ylab=patternname,main=paste(patternname, 'in time'))
#Note the big difference in the resulting patterns in time depending on the unit used.
}

\keyword{ Deposition rate }
\keyword{ Isotope }