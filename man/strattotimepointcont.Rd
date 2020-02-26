\name{strattotimepointcont}
\alias{strattotimepointcont}

\title{
Transform Points and Isotope Ratios from Stratigraphic Height into Time (Continuos)
}
\description{
This function takes a continuous deposition rate and reconstructs at what time a given stratigraphic height was deposited. Essentially a wrapper of \code{pointtransform(...,direction='height to time',depositionmodel='piecewise linear deposition rate')} kept for backwards compatibility of the package. It is recommended to use \code{pointtransform} instead of this function.
}
\usage{
strattotimepointcont(x, xdep, ydep, hiatuslist = list(), unit = "sediment per time")
}

\arguments{
  \item{x}{
Vector containing the stratigraphic heights whose time of deposition are to be determined
}
\item{xdep}{
Vector of strictly increasing real numbers 
}
  \item{ydep}{
Vector of the same length of \code{xdep}, containing strictly positive real numbers. \code{xdep} and \code{ydep} define the deposition rate in the sense that \code{depositionrate=approxfun(xdep,ydep)}. The units of te deposition rate are by default \code{"sediment per time"} and can be changed using the optional input \code{unit}
}
  \item{hiatuslist}{
OPTIONAL, default is an empty list. A list of hiatuses to insert. Every entry in the list corresponds to a hiatus, which is given as a vector of length two. The first entry of the vector is the stratigraphic height of the hiatus, the second entry is the duration of the hiatus.
}
  \item{unit}{
OPTIONAL, default is "sediment per time". Either "sediment per time" or "time per sediment". Determines the unit that is assumed for the deposition rate.
}
}

\value{
Returns a list containing
\item{age }{vector containing times of deposition of the stratigraphic heights given by the input \code{x}}
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
xdep=seq(0,12,length.out=100)
ydep=splinefunH(x=c(0,2,4,6,8,10,12),y=c(1,5,6,1,0.5,1,6),m=c(0,1.5,-0.5,-0.5,0,0.5,0))(xdep)
#unit of deposition rate is sediment per time unit (default setting)
usedunit="sediment per time" #unit of deposition rate is sediment per time unit (default setting)
#Plot deposition rate
plot(xdep,ydep,type='l',main='Deposition Rate',xlab='Stratigraphic Height',ylab=usedunit) 

##at what time was the point with tratigraphic height 9 deposited?
strattotimepointcont(9,xdep,ydep) 

##create an age model
stratheight=seq(min(xdep),max(xdep),length.out=100)
usedunit="sediment per time"
reslist=strattotimepointcont(stratheight,xdep,ydep,unit=usedunit)
#plot age model
plot(reslist$age,reslist$height,xlab='Time',ylab='Stratigraphic Height',
  main=paste('Age model with unit',usedunit))
#create age model but with other units for sedimentn input
usedunit='time per sediment'
reslist=strattotimepointcont(stratheight,xdep,ydep,unit=usedunit)
#plot age model (note the difference this setting makes)
plot(reslist$age,reslist$height,xlab='Time',ylab='Stratigraphic Height',
  main=paste('Age model with unit',usedunit))

##create age model with a hiatus
stratigraphicheight=5 #strat. height of the hiatus
duration=2 #duration of the hiatus
hiatuslist=list(c(stratigraphicheight,duration)) #required input format for hiatuses
reslist=strattotimepointcont(stratheight,xdep,ydep,hiatuslist=hiatuslist)
#!using default setting for units again!
plot(reslist$age,reslist$height,xlab='Time',ylab='stratigraphic height')

##Transform isotope ratios
#create fake oxygen isotope curves
samplelocation=sort(runif(20,min=min(xdep),max=max(xdep))) #where the samples are taken
isotoperatio=sin(samplelocation)*rnorm(length(samplelocation)) #isotope ratios of the samples
plot(samplelocation,isotoperatio,type='l',xlab='Stratigraphic Height',ylab='Isotope Ratio')
#transform only (!) sample locations, NOT values
reslist=strattotimepointcont(samplelocation,xdep,ydep)
#Isotope ratios in time 
plot(reslist$age,isotoperatio,type='l',xlab='Time',ylab='Isotope Ratio')
}

\keyword{ Deposition rate }
\keyword{ Isotope }
\keyword{ strattotime }
