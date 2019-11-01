\name{timetostratratebin}
\alias{timetostratratebin}

\title{
Transform Patterns from Time into Stratigraphic Height (Binned)
}

\description{
This function transforms a given binned temporal pattern into the corresponding stratigraphic pattern. The transformation is based on the age model derived from the given deposition rate. Essentially a wrapper of \code{patterntransform(...,direction='time to height',depositionmodel='binned deposition rate',patternmode='binned')} kept for backwards compatibility of the package.
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
Positive vector of length \code{length(binborder)-1}. Defines the temporal pattern in the bins, i.e. \code{signalval[i]} is the value of the temporal pattern between the times \code{binborder[i]} and  \code{binborder[i+1]}
}
  \item{pos}{
OPTIONAL. Vector of points in time that will be transformed into stratigraphic height, and where the stratigraphic pattern will evaluated. Default setting is a vector spanning the interval of interest with evenly spaced points (no. is chosen by the code)
}
}

\value{
A list containing the following entries:
\item{height }{vector containing stratigraphic heights}
\item{val }{vector containing the values of the stratigraphic pattern at the height given by \code{height}}
NA is returned if the values coincide with a hiatus or when the deposition rate/temporal pattern are undefined.
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
## Define deposition rate
binborder=1:6 #temporal bins for the deposition rate
depoval=c(5,4,3,1,2) #deposition rate in the bins
#plot deposition rate
depositionrate=approxfun(binborder,c(depoval,tail(depoval,1)),method="constant",yleft=NA,yright=NA)
plot(depositionrate(seq(from=min(binborder),to=max(binborder),length.out=100)),xlab='time',
  ylab='deposition rate', main='Deposition rate')
## Define temporal pattern
signalval=c(1,0,5,2,1)
#plot temporal pattern
temporalpat=approxfun(binborder,c(signalval,tail(signalval,1)),method="constant",
  yleft=NA,yright=NA)
plot(temporalpat(seq(from=min(binborder),to=max(binborder),length.out=100)),xlab='time',
  ylab='Intensity', main='Temporal Pattern')

## Transform temporal pattern into stratigraphic pattern
reslist=timetostratratebin(binborder,depoval,signalval)
#plot resulting stratigraphic pattern
plot(reslist$height,reslist$val,xlab='Stratigraphic Height',ylab='Intensity',
  main='Stratigraphic Pattern')
## With removal of sediment (hiatus)
depoval=c(5,4,-3,1,2) #erosion rate is 3 in the middle time bin
reslist=timetostratratebin(binborder,depoval,signalval)
#plot resulting stratigraphic pattern
plot(reslist$height,reslist$val,xlab='Stratigraphic Height',ylab='Intensity',
  main='Stratigraphic Pattern')
}

\keyword{ Deposition rate }
\keyword{ Stratigraphy }
