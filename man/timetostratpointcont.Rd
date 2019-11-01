\name{timetostratpointcont}
\alias{timetostratpointcont}

\title{
Transform Points and Isotope Ratios from Time into Stratigraphic Height (Continuous)
}

\description{
This function takes a continuous deposition rate to determined what stratigraphic height of a section was deposited at a given time. Essentially a wrapper of \code{pointtransform(...,direction='time to height',depositionmodel='piecewise linear deposition rate')} kept for backwards compatibility of the package.
}

\usage{
timetostratpointcont(x, xdep, ydep)
}

\arguments{
\item{x}{
Vector containing the times whose stratigraphic heights are to be determined
}
\item{xdep}{
Vector of strictly increasing real numbers 
}
  \item{ydep}{
Vector of the same length of \code{xdep}. \code{xdep} and \code{ydep} define the deposition rate in the sense that \code{depositionrate=approxfun(xdep,ydep)}.
}
}

\value{
Returns a list containing
\item{height }{vector containing the stratigraphic heights that were deposited at the times given by the input \code{x}}
\item{age }{vector containing the ages at which the stratigraphic heights given by \code{height} were deposited. Essentially a duplicate of the input \code{x}}
Output of \code{NA} in \code{height} indicates that some values coincide with a hiatus or are located at places where the deposition rate is undefined.
}

\references{
Hohmann, Niklas. 2018. Quantifying the Effects of Changing Deposition Rates and Hiatii on the Stratigraphic Distribution of Fossils. <doi:10.13140/RG.2.2.23372.51841>
}

\author{
Niklas Hohmann
}

\seealso{
\code{\link{pointtransform}} for the more general function underlying this wrapper

For an overview of the functions in the DAIME package, see its vignette (available via \code{vignette("DAIME")} )
}
\examples{
##define deposition rate
xdep=seq(0,12,length.out=100)
ydep=splinefunH(x=c(0,2,4,6,8,10,12),y=c(1,5,6,1,0.5,1,6),m=c(0,1.5,-0.5,-0.5,0,0.5,0))(xdep)
#Plot deposition rate
plot(xdep,ydep,type='l',main='Deposition Rate',xlab='Time',ylab='Deposition Rate') 

##at what stratigraphic height can an object be found that was deposited in the 
#sediment after 9 time units?
timetostratpointcont(9,xdep,ydep) 

##create age model
#points that will be transformed into stratigraphic height
time=seq(min(xdep),max(xdep),length.out=100) 
reslist=timetostratpointcont(time,xdep,ydep)
#plot age model
plot(reslist$height,reslist$age,type='l',ylab='Time',xlab='Stratigraphic Height',main='Age model')
#age model with Removal of sediment
#define deposition rate with negative deposition rate, e.g. removal of sediment
xdep=seq(0,12,length.out=100)
ydep=splinefunH(x=c(0,2,4,6,8,10,12),y=c(1,5,6,0,-2,1,6),m=c(0,1.5,-0.5,-0.5,0,0.5,0))(xdep)
#Plot deposition rate
plot(xdep,ydep,main='Deposition Rate',xlab='Time',ylab='Deposition Rate')
#sediment removal is done automatically. e.g. create age model as before:
reslist=timetostratpointcont(time,xdep,ydep)
#plot age model
plot(reslist$height,reslist$age,type='l',ylab='Time',xlab='Stratigraphic Height',main='Age model')


#transform isotope ratio curves
ydep=splinefunH(x=c(0,2,4,6,8,10,12),y=c(1,5,6,1,0.5,1,6),m=c(0,1.5,-0.5,-0.5,0,0.5,0))(xdep)
##create fake ratios and sample locations
sampletime=sort(runif(20,min=min(xdep),max=max(xdep))) #times where the samples were taken
isotoperatio=sin(sampletime)*rnorm(length(sampletime)) #isotope ratios
plot(sampletime,isotoperatio,type='l',xlab='Time',ylab='Isotope Ratio')
#transform only (!) sample times, NOT values
reslist=timetostratpointcont(sampletime,xdep,ydep)
#this is the resulting isotope ratio curve in stratigraphic height
plot(reslist$height,isotoperatio,type='l',xlab='Stratigraphic Height',ylab='Isotope Ratio')
}

\keyword{ DAIME }
\keyword{ Deposition rate }
\keyword{ Isotope }
\keyword{ timetostrat }
\keyword{ Stratigraphy }
