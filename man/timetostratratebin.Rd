\name{timetostratratebin}
\alias{timetostratratebin}

\title{
Transform Rates from Time into Stratigraphic Height (Binned)
}

\description{
This function transforms a given binned temporal rate/signal into the corresponding stratigraphic rate/signal. The transformation is based on the age model constructed from the deposition rate.

This can for example be used to predict the patterns of (first/last) fossil occurrences in the outcrop under a given depositional environment.
}

\usage{
timetostratratebin(binborder, depoval, signalval, pos = NULL)
}

\arguments{
  \item{binborder}{
Vector of strictly increasing numerical values. Defines the borders of the bins in time
}
  \item{depoval}{
Strictly positive vector of length \code{length(binborder)-1}. Defines the deposition rate in the bins, i.e. \code{depoval[i]} is the deposition rate in the time interval between \code{binborder[i]} and  \code{binborder[i+1]}
}
  \item{signalval}{
Positive vector of length \code{length(binborder)-1}. Defines the temporal rate/signal in the bins, i.e. \code{signalval[i]} is the temporal rate between the times \code{binborder[i]} and  \code{binborder[i+1]}
}
  \item{pos}{
OPTIONAL. Vector of points in time that will be transformed into stratigraphic rate, and where the stratigraphic rate will evaluated. Default setting is a vector spanning the interval of interest with evenly spaced points (no. is chosen by the code)
}
}

\value{
A list containing the following entries
\item{x }{x values, corresponding to the locations of \code{pos}, transformed into stratigraphic height }
\item{y }{function values of the transformed rate at x}
The results return NA if the values coincide with a hiatus or both deposition rate and temporal rate are undefined
}

\references{
Hohmann, Niklas. 2018. Quantifying the Effects of Changing Deposition Rates and Hiatii on the Stratigraphic Distribution of Fossils. <doi:10.13140/RG.2.2.23372.51841>
}

\author{
Niklas Hohmann
}

\seealso{
For an overview of the functions in the DAIME package, see its vignette (available via \code{vignette("DAIME")} )
}

\examples{
## Define deposition rate
binborder=1:6 #temporal bins for the deposition rate
depoval=c(5,4,3,1,2) #deposition rate in the bins
#plot deposition rate
depositionrate=approxfun(binborder,c(depoval,tail(depoval,1)),method="constant",yleft=NA,yright=NA)
plot(depositionrate(seq(from=min(binborder),to=max(binborder),length.out=100)),xlab='time',
  ylab='deposition rate', main='Deposition rate')
## Define temporal rate
signalval=c(1,0,5,2,1)
#plot temporal rate
temporalrate=approxfun(binborder,c(signalval,tail(signalval,1)),method="constant",
  yleft=NA,yright=NA)
plot(temporalrate(seq(from=min(binborder),to=max(binborder),length.out=100)),xlab='time',
  ylab='temporal rate', main='Temporal rate')

## Transform temporal rate into stratigraphic rate
reslist=timetostratratebin(binborder,depoval,signalval)
#plot resulting stratigraphic rate
plot(reslist$height,reslist$val,xlab='Stratigraphic Height',ylab='Intensity',
  main='Stratigraphic Rate')
## With removal of sediment (hiatus)
depoval=c(5,4,-3,1,2) #erosion rate is 3 in the middle time bin
reslist=timetostratratebin(binborder,depoval,signalval)
#plot resulting stratigraphic rate
plot(reslist$height,reslist$val,xlab='Stratigraphic Height',ylab='Intensity',
  main='Stratigraphic Rate')
}

\keyword{ Deposition rate }
\keyword{ Stratigraphy }
