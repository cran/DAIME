\name{strattotimeratecont}
\alias{strattotimeratecont}

\title{
Transform Rates from Stratigraphic Height into Time
}
\description{
This function takes a deposition rate and transforms a given stratigraphic rate/signal from stratigraphic height into time.

This can for example be used to reconstruct the patterns of (first/last) fossil occurrences  in time, given  a given depositional environment.
}
\usage{
strattotimeratecont(xdep, ydep, xsig, ysig, pos = NULL,
  hiatuslist = list(), unit = "sediment per time")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xdep}{
vector of real numbers 
}
  \item{ydep}{
vector of the same length of \code{xdep}, containing positive real numbers. \code{xdep} and \code{ydep} define the deposition rate in the sense that \code{depositionrate=approxfun(xdep,ydep)}. The units of te deposition rate are by default \code{"sediment per time"} and can be changed using the optional input \code{unit}
}
  \item{xsig}{
vector of real numbers
}
  \item{ysig}{
vector of the same length of \code{xsig}, containing positive real numbers. \code{xsig} and \code{ysig} define the stratigraphic rate in the sense that \code{stratigraphicrate=approxfun(xsig,ysig)}.
}
  \item{pos}{
OPTIONAL, default is \code{NULL}. At what height is the stratigraphic rate transformed? (by default \code{pos=xsig} is used)
}
  \item{hiatuslist}{
OPTIONAL, default is an empty list. A list of hiatuses to insert. Every entry in the list corresponds to a hiatus, which is given as a vector of length two. The first entry of the vector is the stratigraphic height of the hiatus, the second entry is the duration of the hiatus.
}
  \item{unit}{
OPTIONAL, default is "sediment per time". Either "sediment per time" or "time per sediment". Determines the unit of the deposition rate.
}
}

\value{
Returns a list containing
\item{x }{vector containing times (of deposition of the points in \code{xdep} (\code{pos} if used))}
\item{y }{vector containing the values of the temporal rate/signal at \code{x}}
\code{x} and \code{y} define the temporal rate/signal in the sense that \code{temporalrate=approxfun(x,y)}
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
}
\examples{
#define deposition rate
xdep=seq(0,12,length.out=100)
ydep=splinefunH(x=c(0,2,4,6,8,10,12),y=c(1,5,6,1,0.5,1,6),m=c(0,1.5,-0.5,-0.5,0,0.5,0))(xdep)
usedunit="sediment per time" #unit of deposition rate is sediment per time unit (default setting)
plot(xdep,ydep,type='l',main='Deposition Rate',xlab='Stratigraphic Height',
  ylab='Sediment per time unit')
#Plot deposition rate

#define a rate ("signal")
xsig= seq(0,12,length.out=100)
ysig=splinefunH(x=c(0,4,12),y=c(0.5,2,0.5),m=c(0,0,0))(xsig)
plot(xsig,ysig,type='l',main='Stratigraphic Signal/Rate',
  xlab='Stratigraphic Height',ylab='intensity')

resultlist=strattotimeratecont(xdep,ydep,xsig,ysig,unit=usedunit)
plot(resultlist,type='l',xlab='Time',ylab='Intensity',main='Temporal Rate/Signal')

usedunit="time per sediment" #now interpret the deposition rate in different units
resultlist=strattotimeratecont(xdep,ydep,xsig,ysig,unit=usedunit)
plot(resultlist,type='l',xlab='Time',ylab='Intensity',main='Temporal Rate/Signal')
}

\keyword{ DepAIM }
\keyword{ Deposition rate }
\keyword{ Stratigraphic rate }
\keyword{ strattotime }
