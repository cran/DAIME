\name{strattotimeratecont}
\alias{strattotimeratecont}

\title{
Transform Rates from Stratigraphic Height into Time (Continuous)
}

\description{
This function transforms a given continuous stratigraphic rate/signal into the corresponding temporal rate/signal. The transformation is based on the age model constructed from the deposition rate.

This can for example be used to reconstruct the patterns of (first/last) fossil occurrences  in time, given a given depositional environment.
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
Vector of the same length of \code{xsig}, containing positive real numbers. \code{xsig} and \code{ysig} define the stratigraphic rate in the sense that \code{stratigraphicrate=approxfun(xsig,ysig)}.
}
  \item{pos}{
OPTIONAL. Vector of points in stratigrapic height that will be transformed into time, and determine where the temporal rate is evaluated. Default setting is a vector spanning the interval of interest with evenly spaced points (no. is chosen by the code)
}
  \item{hiatuslist}{
OPTIONAL, default is an empty list. A list of hiatuses to insert. Every entry in the list corresponds to a hiatus, which is given as a vector of length two. The first entry of the vector is the stratigraphic height of the hiatus, the second entry is the duration of the hiatus.
}
  \item{unit}{
OPTIONAL, default is "sediment per time". Either "sediment per time" or "time per sediment". Determines the unit that is assumed for the deposition rate.
}
}

\value{
Returns a list containing
\item{time }{vector containing times}
\item{val }{vector containign the values of the temporal rate/signal at the times given by the first entry of the list}
The temporal rate is then approximated by \code{approxfun(output$time,output$val)}. Output of \code{NA} indicates that some values coincide with a hiatus or are located at places where the deposition rate/stratigraphic rate are undefined.
}

\references{
Hohmann, Niklas. 2018. Quantifying the Effects of Changing Deposition Rates and Hiatii on the Stratigraphic Distribution of Fossils. <doi:10.13140/RG.2.2.23372.51841>
}

\author{
Niklas Hohmann
}

\seealso{
\code{\link{timetostratratecont}} for the transformation of rates from time into stratigraphic height,
\code{\link{strattotimepointcont}} for the transformation of single samples or isotope ratios from stratigraphic height into time.

For an overview of the functions in the DAIME package, see its vignette (available via \code{vignette("DAIME")} )
}
\examples{
##define deposition rate and stratigraphic rate
xdep=seq(0,12,length.out=100)
ydep=splinefunH(x=c(0,2,4,6,8,10,12),y=c(1,5,6,1,0.5,1,6),m=c(0,1.5,-0.5,-0.5,0,0.5,0))(xdep)
usedunit="sediment per time" #unit of deposition rate is sediment per time unit (default setting)
#Plot deposition rate
plot(xdep,ydep,type='l',main='Deposition Rate',xlab='Stratigraphic Height',
  ylab=usedunit)
#define a rate ("signal")
xsig= seq(from=min(xdep),to=max(xdep),length.out=100)
ysig=splinefunH(x=c(0,4,12),y=c(0.5,2,0.5),m=c(0,0,0))(xsig) #function values of the signal
plot(xsig,ysig,type='l',main='Stratigraphic Signal/Rate',
  xlab='Stratigraphic Height',ylab='Stratigraphic Rate')
  
##transform stratigraphic rate into temporal rate
usedunit="sediment per time" #unit of deposition rate is sediment per time unit (default setting)
reslist=strattotimeratecont(xdep,ydep,xsig,ysig,unit=usedunit)
#plot temporal rate (transformed stratigrapic rate)
plot(reslist$time,reslist$val,type='l',xlab='Time',ylab='temporal rate',
  main=paste('deposition rate interpreted as',usedunit))
#now using same input, but different interpretation of units
usedunit="time per sediment" #now interpret the deposition rate in different units
reslist=strattotimeratecont(xdep,ydep,xsig,ysig,unit=usedunit)
#plot temporal rate (transformed stratigrapic rate)
plot(reslist$time,reslist$val,type='l',xlab='Time',ylab='temporal rate',
  main=paste('deposition rate interpreted as',usedunit))
#note how different the results look!!

#Insert hiatus
stratigraphicheight=5 #strat. height of the hiatus
duration=4 #duration of the hiatus
hiatuslist=list(c(stratigraphicheight,duration)) #required input format for hiatuses
reslist=strattotimeratecont(xdep,ydep,xsig,ysig,hiatuslist=hiatuslist) #unit is back to default
plot(reslist$time,reslist$val,xlab='Time',ylab='Intensity',main='Temporal Rate/Signal')
#hiatus is approx between 2 and 5
}

\keyword{ Deposition rate }
\keyword{ Stratigraphic rate }
\keyword{ strattotime }
\keyword{ Stratigraphy }
