\name{strattotimepointbin}
\alias{strattotimepointbin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transform Points and Isotope Ratios from Stratigraphic Height into Time (Binned)
}
\description{
This function takes a binned deposition rate and reconstructs at what time a given stratigraphic height was deposited. Essentially a wrapper of \code{pointtransform(...,direction='height to time',depositionmodel='binned deposition rate')} kept for backwards compatibility of the package. It is recommended to use \code{pointtransform} instead of this function.
}
\usage{
strattotimepointbin(x, binborder, depoval, hiatuslist = list(),
  unit = "sediment per time")
}

\arguments{
  \item{x}{
Vector containing the stratigraphic heights whose time of deposition are to be determined
}
  \item{binborder}{
Vector of strictly increasing numerical values. Defines the borders of the bins in stratigraphic height
}
  \item{depoval}{
Strictly positive vector of length \code{length(binborder)-1}. Defines the deposition rate in the bins, i.e. \code{depoval[i]} is the deposition rate between the stratigraphic heights \code{binborder[i]} and  \code{binborder[i+1]}
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
\item{age }{vector containing ages of deposition of the stratigraphic heights given by the input \code{x}}
\item{height }{vector containing the heights that were deposited at the times given by \code{age}. Essentially a duplicate of the input \code{x}}
Output of \code{NA} in \code{age} indicates that some values coincide with a hiatus or are located at places where the deposition rate is undefined.
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
binborder=1:6 #bins in tratigraphic height
depoval=c(1,4,5,2,2)
usedunit="sediment per time" #unit of deposition rate is sediment per time unit (default setting)
#Plot deposition rate
plot(approx(binborder,c(depoval,tail(depoval,1)),method='constant',xout=seq(min(binborder),
  max(binborder),length.out=100)),type='l',main='Deposition Rate',xlab='Stratigraphic Height',
  ylab=usedunit,ylim=c(0,max(depoval)))
  
##at what time was the point with tratigraphic height 3.5 deposited?
strattotimepointbin(3.5,binborder,depoval) 

##create an age model
stratheight=seq(min(binborder),max(binborder),length.out=100)
reslist=strattotimepointbin(stratheight,binborder,depoval,unit=usedunit)
#plot age model
usedunit="sediment per time"
plot(reslist$age,reslist$height,xlab='Time',ylab='Stratigraphic Height',
  main=paste('Age model with unit',usedunit))
#create age model but with other units for sedimentn input
usedunit='time per sediment'
reslist=strattotimepointbin(stratheight,binborder,depoval,unit=usedunit)
#plot age model (note the difference this setting makes)
plot(reslist$age,reslist$height,xlab='Time',ylab='Stratigraphic Height',
  main=paste('Age model with unit',usedunit))

##create age model with a hiatus
stratigraphicheight=5 #strat. height of the hiatus
duration=2 #duration of the hiatus
hiatuslist=list(c(stratigraphicheight,duration)) #required input format for hiatuses
reslist=strattotimepointbin(stratheight,binborder,depoval,hiatuslist=hiatuslist)
#!using default setting for units again!
plot(reslist$age,reslist$height,xlab='Time',ylab='stratigraphic height') 

##Transform isotope ratios
depoval=c(1,4,0.1,2,2) 
#create fake oxygen isotope curves
samplelocation=sort(runif(20,min=min(binborder),max=max(binborder))) #where the samples are taken
isotoperatio=sin(samplelocation)*rnorm(length(samplelocation)) #isotope ratios of the samples
plot(samplelocation,isotoperatio,type='l',xlab='Stratigraphic Height',ylab='Isotope Ratio')
#transform only (!) sample locations, NOT values
reslist=strattotimepointbin(samplelocation,binborder,depoval)
#Isotope ratios in time 
plot(reslist$age,isotoperatio,type='l',xlab='Time',ylab='Isotope Ratio')

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Deposition rate }
\keyword{ Stratigraphy }
\keyword{ Isotope }
