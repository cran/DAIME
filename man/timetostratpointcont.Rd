\name{timetostratpointcont}
\alias{timetostratpointcont}

\title{
Transform Points and Isotope Ratios from Time into Stratigraphic Height
}

\description{
This function takes a deposition rate and determines what stratigraphic height was deposited at a given time.

This can be used to create age models, determine the location of single samples in the outcrop that were deposited at a given time, and transform isotope ratios from time into stratigraphic height (see examples).
}

\usage{
timetostratpointcont(x, xdep, ydep)
}

\arguments{
  \item{x}{
vector containing the points in time whose position in stratigraphic height are to be determined
}
  \item{xdep}{
vector of real numbers 
}
  \item{ydep}{
vector of the same length of \code{xdep}, containing positive real numbers. \code{xdep} and \code{ydep} define the deposition rate in the sense that \code{depositionrate=approxfun(xdep,ydep)}
}
}

\value{
A vector of the same length as \code{x}, containing the stratigraphic heights that are deposited at the times given in \code{x}.
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
#Plot deposition rate
plot(xdep,ydep,type='l',main='Deposition Rate',xlab='Time',ylab='Deposition Rate') 

#at what stratigraphic height is the point which was deposited after 9 time units?
time=9
timetostratpointcont(time,xdep,ydep) 

#get the age model
time=seq(min(xdep),max(xdep),length.out=100)
stratheight=timetostratpointcont(time,xdep,ydep)
#plot age model
plot(time,stratheight,type='l',xlab='Time',ylab='Stratigraphic Height',main='Age model')

#transform a fake isotope ratio curve
##create fake ratios and sample locations
sampletime=sort(runif(20,min=0,max=12))
isotoperatio=sin(sampletime)*rnorm(length(sampletime))
plot(sampletime,isotoperatio,type='l',xlab='Time',ylab='Isotope Ratio')

#transform only (!) sample times, NOT values
samplelocation=timetostratpointcont(sampletime,xdep,ydep)
#this is the resulting isotope ratio curve
plot(samplelocation,isotoperatio,type='l',xlab='Time',ylab='Isotope Ratio')
}

\keyword{ DepAIM }
\keyword{ Deposition rate }
\keyword{ Isotope }
\keyword{ timetostrat }
