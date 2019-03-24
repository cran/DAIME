\name{strattotimepointcont}
\alias{strattotimepointcont}

\title{
Transform Points and Isotope Ratios from Stratigraphic Height into Time
}
\description{
This function takes a deposition rate and determines at what time a given stratigraphic height was deposited.

This can be used to create age models, determine the time of deposition of single samples in the outcrop, and transform isotope ratios from stratigraphic height into time (see examples).
}
\usage{
strattotimepointcont(x, xdep, ydep, hiatuslist = list(), unit = "sediment per time")
}

\arguments{
  \item{x}{
  vector containing stratigraphic heights for which the time of deposition is determined
}
  \item{xdep}{
vector of real numbers 
}
  \item{ydep}{
vector of the same length of \code{xdep}, containing positive real numbers. \code{xdep} and \code{ydep} define the deposition rate in the sense that \code{depositionrate=approxfun(xdep,ydep)}. The units of te deposition rate are by default \code{"sediment per time"} and can be changed using the optional input \code{unit}.
}
  \item{hiatuslist}{
OPTIONAL, default is an empty list. A list of hiatuses to insert. Every entry in the list corresponds to a hiatus, which is given as a vector of length two. The first entry of the vector is the stratigraphic height of the hiatus, the second entry is the duration of the hiatus.
}
  \item{unit}{
OPTIONAL, default is "sediment per time". Either "sediment per time" or "time per sediment". Determines the unit of the deposition rate
}
}

\value{
A vector of the same length as \code{x}, containing the times at which the points in the section given by  \code{x} are deposited.
}
\references{
Hohmann, Niklas. 2018. Quantifying the Effects of Changing Deposition Rates and Hiatii on the Stratigraphic Distribution of Fossils. <doi:10.13140/RG.2.2.23372.51841>
}
\author{
Niklas Hohmann
}

\seealso{
\code{\link{strattotimeratecont}} for the transformation of rates from stratigraphic height into time,
\code{\link{timetostratpointcont}} for the transformation of single samples or isotope ratios from time into stratigraphic height.
}
\examples{
#define deposition rate
xdep=seq(0,12,length.out=100)
ydep=splinefunH(x=c(0,2,4,6,8,10,12),y=c(1,5,6,1,0.5,1,6),m=c(0,1.5,-0.5,-0.5,0,0.5,0))(xdep)
#unit of deposition rate is sediment per time unit (default setting)
usedunit="sediment per time"
#Plot deposition rate
plot(xdep,ydep,type='l',main='Deposition Rate',xlab='Stratigraphic Height',ylab=usedunit) 

#at what time was the point with tratigraphic height 9 deposited?
stratheight=9
strattotimepointcont(stratheight,xdep,ydep) 

#get the age model
stratheight=seq(min(xdep),max(xdep),length.out=100)
time=strattotimepointcont(stratheight,xdep,ydep)
#plot age model
plot(time,stratheight,type='l',xlab='Time',ylab='Stratigraphic Height',main='Age model')

samplelocation=sort(runif(20,min=0,max=12))
isotoperatio=sin(samplelocation)*rnorm(length(samplelocation))
plot(samplelocation,isotoperatio,type='l',xlab='Stratigraphic Height',ylab='Isotope Ratio')

#transform only (!) sample locations, NOT values
sampletime=strattotimepointcont(samplelocation,xdep,ydep)
plot(sampletime,isotoperatio,type='l',xlab='Time',ylab='Isotope Ratio')
}

\keyword{ DepAIM }
\keyword{ Deposition rate }
\keyword{ Isotope }
\keyword{ strattotime }
