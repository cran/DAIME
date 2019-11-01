\name{strattotimeratecont}
\alias{strattotimeratecont}

\title{
Transform Patterns from Stratigraphic Height into Time (Continuous)
}

\description{
This function transforms a given continuous stratigraphic pattern into the corresponding temporal pattern. The transformation is based on the age model derived from the given deposition rate. Essentially a wrapper of \code{patterntransform(...,direction='height to time',depositionmodel='piecewise linear deposition rate',patternmode='piecewise linear')} kept for backwards compatibility of the package.
}

\usage{
strattotimeratecont(xdep, ydep, xsig, ysig, pos = NULL,   hiatuslist = list(),
  unit = "sediment per time")
}

\arguments{
  \item{xdep}{
Vector of strictly increasing real numbers 
}
  \item{ydep}{
Vector of the same length of \code{xdep}, containing strictly positive real numbers. \code{xdep} and \code{ydep} define the deposition rate in the sense that \code{depositionrate=approxfun(xdep,ydep)}. The units of te deposition rate are by default \code{"sediment per time"} and can be changed using the optional input \code{unit}
}
  \item{xsig}{
Vector of strictly increasing real numbers 
}
  \item{ysig}{
Vector of the same length of \code{xsig}, containing positive real numbers. \code{xsig} and \code{ysig} define the stratigraphic pattern in the sense that \code{stratigraphicrate=approxfun(xsig,ysig)}.
}
  \item{pos}{
OPTIONAL. Vector of points in stratigrapic height that will be transformed into time, and determine where the temporal pattern is evaluated. Default setting is a vector spanning the interval of interest with evenly spaced points (no. is chosen by the code)
}
  \item{hiatuslist}{
OPTIONAL, default is an empty list. A list of hiatuses to insert. Every entry in the list corresponds to a hiatus, which is given as a vector of length two. The first entry of the vector is the stratigraphic height of the hiatus, the second entry is the duration of the hiatus.
}
  \item{unit}{
OPTIONAL, default is "sediment per time". Either "sediment per time" or "time per sediment". Determines the unit that is assumed for the deposition rate.
}
}

\value{
Returns a list containing:
\item{age }{vector containing times}
\item{val }{vector containign the values of the temporal pattern at the times given by the first entry of the list}
The temporal pattern can then be approximated by \code{approxfun(output$age,output$val)}. Output of \code{NA} indicates that some values coincide with a hiatus or are located at places where the deposition rate/stratigraphic pattern are undefined.
}

\references{
Hohmann, Niklas. 2018. Quantifying the Effects of Changing Deposition Rates and Hiatii on the Stratigraphic Distribution of Fossils. <doi:10.13140/RG.2.2.23372.51841>
}

\author{
Niklas Hohmann
}

\seealso{
\code{\link{patterntransform}} for the more general function underlying this wrapper

For an overview of the functions in the DAIME package, see its vignette (available via \code{vignette("DAIME")} )
}
\examples{
##define deposition rate and stratigraphic pattern
xdep=seq(0,12,length.out=100)
ydep=splinefunH(x=c(0,2,4,6,8,10,12),y=c(1,5,6,1,0.5,1,6),m=c(0,1.5,-0.5,-0.5,0,0.5,0))(xdep)
usedunit="sediment per time" #unit of deposition rate is sediment per time unit (default setting)
#Plot deposition rate
plot(xdep,ydep,type='l',main='Deposition Rate',xlab='Stratigraphic Height',
  ylab=usedunit)
#define a pattern
xsig= seq(from=min(xdep),to=max(xdep),length.out=100)
ysig=splinefunH(x=c(0,4,12),y=c(0.5,2,0.5),m=c(0,0,0))(xsig) #function values of the signal
plot(xsig,ysig,type='l',main='Stratigraphic Pattern',
  xlab='Stratigraphic Height',ylab='Intensity')
  
##transform stratigraphic pattern into temporal pattern
usedunit="sediment per time" #unit of deposition rate is sediment per time unit (default setting)
reslist=strattotimeratecont(xdep,ydep,xsig,ysig,unit=usedunit)
#plot temporal pattern (transformed stratigrapic pattern)
plot(reslist$age,reslist$val,type='l',xlab='Time',ylab='Intensity',
  main=paste('deposition rate interpreted as',usedunit))
#now using same input, but different interpretation of units
usedunit="time per sediment" #now interpret the deposition rate in different units
reslist=strattotimeratecont(xdep,ydep,xsig,ysig,unit=usedunit)
#plot temporal pattern (transformed stratigrapic pattern)
plot(reslist$age,reslist$val,type='l',xlab='Time',ylab='Intensity temporal pattern',
  main=paste('deposition rate interpreted as',usedunit))
#note how different the results look!!

#Insert hiatus
stratigraphicheight=5 #strat. height of the hiatus
duration=4 #duration of the hiatus
hiatuslist=list(c(stratigraphicheight,duration)) #required input format for hiatuses
reslist=strattotimeratecont(xdep,ydep,xsig,ysig,hiatuslist=hiatuslist) #unit is back to default
plot(reslist$age,reslist$val,xlab='Time',ylab='Intensity',main='Temporal Pattern')
#hiatus is approx between 2 and 5
}

\keyword{ Deposition rate }
\keyword{ Stratigraphic rate }
\keyword{ strattotime }
\keyword{ Stratigraphy }
