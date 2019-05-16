\name{timetostratratecont}
\alias{timetostratratecont}

\title{
Transform Rates from Time into Stratigraphic Height (Continuous)
}

\description{
This function transforms a given continuous temporal rate/signal into the corresponding stratigraphic rate/signal. The transformation is based on the age model constructed from the deposition rate.

This can for example be used to predict the patterns of (first/last) fossil occurrences in the outcrop under a given depositional environment.
}

\usage{
timetostratratecont(xdep, ydep, xsig, ysig, pos = NULL)
}

\arguments{
  \item{xdep}{
Vector of strictly increasing real numbers 
}
  \item{ydep}{
Vector of the same length of \code{xdep}, containing strictly positive real numbers. \code{xdep} and \code{ydep} define the deposition rate in the sense that \code{depositionrate=approxfun(xdep,ydep)}
}
  \item{xsig}{
Vector of strictly increasing real numbers 
}
  \item{ysig}{
Vector of the same length of \code{xsig}, containing positive real numbers. \code{xsig} and \code{ysig} define the temporal rate in the sense that \code{stratigraphicrate=approxfun(xsig,ysig)}.
}
  \item{pos}{
OPTIONAL. Vector of points in time height that will be transformed into stratigraphic height, and determine where the stratigraphic rate is evaluated. Default setting is a vector spanning the interval of interest with evenly spaced points (no. is chosen by the code)
}
}

\value{
Returns a list containing
\item{x }{vector containing stratigraphic heights (that are deposited at times \code{xdep} (\code{pos} if used))}
\item{y }{vector containign the values of the stratigraphic rate/signal at \code{x}}
\code{x} and \code{y} define the stratigraphic rate/signal in the sense that \code{stratsig=approxfun(x,y)}
}

\references{
Hohmann, Niklas. 2018. Quantifying the Effects of Changing Deposition Rates and Hiatii on the Stratigraphic Distribution of Fossils. <doi:10.13140/RG.2.2.23372.51841>
}

\author{
Niklas Hohmann
}

\seealso{
\code{\link{strattotimeratecont}} for the transformation of rates from stratigraphic height into time,
\code{\link{timetostratpointcont}} for the transformation of sigle samples or isotope ratios from time into stratigraphic height.

For an overview of the functions in the DAIME package, see its vignette (available via \code{vignette("DAIME")} )
}

\examples{
##define deposition rate
xdep=seq(0,12,length.out=100)
ydep=splinefunH(x=c(0,2,4,6,8,10,12),y=c(1,5,6,1,0.5,1,6),m=c(0,1.5,-0.5,-0.5,0,0.5,0))(xdep)
#Plot deposition rate
plot(xdep,ydep,type='l',main='Deposition Rate',xlab='Time',ylab='Sediment per Time Unit')

#define a temporal rate ("signal")
xsig= seq(0,12,length.out=100)
ysig=splinefunH(x=c(0,4,12),y=c(0.5,2,0.5),m=c(0,0,0))(xsig)
plot(xsig,ysig,type='l',main='Temporal Signal/Rate',xlab='Time',ylab='Signal Intensity')

##transform temporal rate into stratigraphic rate (signal observable in the section)
reslist=timetostratratecont(xdep,ydep,xsig,ysig)
plot(reslist$height,reslist$val,type='l',xlab='Stratigraphic Height',
  ylab='Signal Intensity',main='Stratigraphic Rate/Signal')
  
##with removal of sediment
ydep=splinefunH(x=c(0,2,4,6,8,10,12),y=c(1,5,6,-2,5,1,6),m=c(0,1.5,-0.5,-0.5,0,0.5,0))(xdep)
#Plot deposition rate
plot(xdep,ydep,type='l',main='Deposition Rate',xlab='Time',ylab='Sediment per Time Unit')
reslist=timetostratratecont(xdep,ydep,xsig,ysig)
plot(reslist$height,reslist$val,type='l',xlab='Stratigraphic Height',
  ylab='Signal Intensity',main='Stratigraphic Rate/Signal',ylim=c(0,1))
#the spike is because the deposition rate is very small when it transitions from 
#negative to positive, generating a punctual extreme condensation
  
}

\keyword{ DAIME}
\keyword{ Deposition rate }
\keyword{ temporal rate }
\keyword{ timetostrat }
