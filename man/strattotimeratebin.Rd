\name{strattotimeratebin}
\alias{strattotimeratebin}

\title{
Transform Rates from Stratigraphic Height into Time (Binned)
}

\description{
This function transforms a given binned stratigraphic rate/signal into the corresponding temporal rate/signal. The transformation is based on the age model constructed from the deposition rate.

This can for example be used to reconstruct the patterns of (first/last) fossil occurrences  in time, given a given depositional environment.
}

\usage{
strattotimeratebin(binborder, depoval, signalval, pos = NULL, hiatuslist = list(), 
  unit = "sediment per time")
}

\arguments{
  \item{binborder}{
Vector of strictly increasing numerical values. Defines the borders of the bins in stratigraphic height
}
  \item{depoval}{
Strictly positive vector of length \code{length(binborder)-1}. Defines the deposition rate in the bins, i.e. \code{depoval[i]} is the deposition rate between the stratigraphic heights \code{binborder[i]} and  \code{binborder[i+1]}
}
  \item{signalval}{
Positive vector of length \code{length(binborder)-1}. Defines the stratigraphic rate/signal in the bins, i.e. \code{signalval[i]} is the stratigraphic between the stratigraphic heights \code{binborder[i]} and  \code{binborder[i+1]}
}
  \item{pos}{
OPTIONAL. Vector of points in stratigrapic height that will be transformed into time, and where  the temporal rate will evaluated. Default setting is a vector spanning the interval of interest with evenly spaced points (no. is chosen by the code)
}
  \item{hiatuslist}{
OPTIONAL, default is an empty list. A list of hiatuses to insert. Every entry in the list corresponds to a hiatus, which is given as a vector of length two. The first entry of the vector is the stratigraphic height of the hiatus, the second entry is the duration of the hiatus.
}
  \item{unit}{
OPTIONAL, default is "sediment per time". Either "sediment per time" or "time per sediment". Determines the unit that is assumed for the input given by \code{depoval}.
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
For an overview of the functions in the DAIME package, see its vignette (available via \code{vignette("DAIME")} )
}

\examples{
##define deposition rate and stratigraphic rate
binborder=1:6 #bins in tratigraphic height
depoval=c(1,4,0.1,2,2)
usedunit="sediment per time" #unit of deposition rate is sediment per time unit (default setting)
#Plot deposition rate
plot(approx(binborder,c(depoval,tail(depoval,1)),method='constant',
  xout=seq(min(binborder),max(binborder),length.out=100)),type='l',main='Deposition Rate',
  xlab='Stratigraphic Height',ylab=usedunit,ylim=c(0,max(depoval)))
#define stratigraphic rate
signalval=c(runif(5))
#plot stratigraphic rate
plot(approx(binborder,c(signalval,tail(signalval,1)),method='constant', 
  xout=seq(min(binborder),max(binborder),length.out=100)),type='l',main='Deposition Rate',
  xlab='Stratigraphic Height', ylab='Stratigraphic Rate',ylim=c(0,max(signalval)))
  
##transform stratigraphic rate into temporal rate
usedunit="sediment per time" #unit of deposition rate is sediment per time unit (default setting)
reslist=strattotimeratebin(binborder,depoval, signalval,unit=usedunit)
#plot results
plot(reslist$time,reslist$val,xlab='time',ylab='temporal rate',
  main=paste('depoval interpreted as',usedunit))
usedunit="time per sediment" #use other input interpretation
reslist=strattotimeratebin(binborder,depoval, signalval,unit=usedunit)
#note how different the results look!!
plot(reslist$time,reslist$val,xlab='time',ylab='temporal rate',
  main=paste('depoval interpreted as',usedunit))
  
##insert a hiatus
stratigraphicheight=3.5 #strat. height of the hiatus
duration=2 #duration of the hiatus
hiatuslist=list(c(stratigraphicheight,duration)) #required input format for hiatuses
#usedunit is back to default setting!
reslist=strattotimeratebin(binborder,depoval, signalval,hiatuslist=hiatuslist)
#the hiatus corresponds to the gap  in the middle of the temporal rate
plot(reslist$time,reslist$val,xlab='time',ylab='temporal rate') 

}

\keyword{ Deposition rate }
\keyword{ Stratigraphy }
