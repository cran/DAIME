\name{timetostratpointbin}
\alias{timetostratpointbin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform Points and Isotope Ratios from Time into Stratigraphic Height (Binned)
}

\description{
This function takes a binned deposition rate to determined what stratigraphic height of a section was deposited at a given time. Essentially a wrapper of \code{pointtransform(...,direction='time to height',depositionmodel='binned deposition rate')} kept for backwards compatibility of the package.
}

\usage{
timetostratpointbin(x, binborder, depoval)
}

\arguments{
  \item{x}{
Vector containing the times whose stratigraphic heights are to be determined
}
  \item{binborder}{
Vector of strictly increasing numerical values. Defines the borders of the bins in time
}
  \item{depoval}{
Vector of length \code{length(binborder)-1}. Defines the deposition rate in the bins, i.e. \code{depoval[i]} is the deposition rate in the time interval between      
\code{binborder[i]} and  \code{binborder[i+1]}
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
##Define deposition rate
binborder=1:6 #temporal bins for the deposition rate
depoval=c(5,4,3,1,2) #deposition rate in the bins
#plot deposition rate
depositionrate=approxfun(binborder,c(depoval,tail(depoval,1)),method="constant",yleft=NA,yright=NA)
plot(depositionrate(seq(from=min(binborder),to=max(binborder),length.out=100)),xlab='time',
  ylab='deposition rate', main='Deposition rate')

##at what stratigraphic height can an object be found that was deposited in the 
##sediment after 5 time units?
timetostratpointbin(5,binborder,depoval) 

##create age model
#points that will be transformed into stratigraphic height
time=seq(from=min(binborder),to=max(binborder),length.out=100) 
reslist=timetostratpointbin(time,binborder,depoval)
#plot age model
plot(reslist$height,reslist$age,type='l',ylab='Time',xlab='Stratigraphic Height',main='Age model')

## Age model with removal of sediment (hiatus)
depoval=c(5,4,-3,1,2) #in the midle time bin, erosion rate is 3
reslist=timetostratpointbin(time,binborder,depoval)
#plot age model. the gap represents the hiatus
plot(reslist$height,reslist$age,type='l',ylab='Time',xlab='Stratigraphic Height',
  main='Age model with erosion')
#A object deposited in the sediment after 3.5 time units is destroyed due to the hiatus:
timetostratpointbin(3.5,binborder,depoval) 

##transform isotope ratio curves
depoval=c(5,4,2,1,0.1)
#create fake ratios and sample locations
sampletime=sort(runif(20,min=min(binborder),max=max(binborder))) #times where the samples were taken
isotoperatio=sin(sampletime)*rnorm(length(sampletime)) #isotope ratios
plot(sampletime,isotoperatio,type='l',xlab='Time',ylab='Isotope Ratio')
#!!transform only (!) sample times, NOT isotope values!!
reslist=timetostratpointbin(sampletime,binborder,depoval)
#this is the resulting isotope ratio curve in stratigraphic height
plot(reslist$height,isotoperatio,type='l',xlab='Stratigraphic Height',ylab='Isotope Ratio')

}

\keyword{ Deposition rate }
\keyword{ Stratigraphy }
