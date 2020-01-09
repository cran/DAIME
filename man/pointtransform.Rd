\name{pointtransform}
\alias{pointtransform}

\title{
Transform points between stratigraphic height and time
}
\description{
This function (1) takes times and determines the stratigraphic heights that were deposited at said times or (2) takes stratigraphic heights and determines their time of deposition.

This can be used to (1) reconstruct the age/stratigraphic height of single objects placed in the sediment (2) construct age models from deposition rates and (3) transform (isotope) ratios or percentages from stratigraphic height into time and vice versa (see examples).
}
\usage{
pointtransform( points , xdep , ydep ,
  direction = 'time to height' , depositionmodel = 'piecewise linear deposition rate' ,
  hiatuslist = list() , unit = 'time per sediment' )
}

\arguments{
  \item{points}{
Vector containing points that are transformed. Whether the entries in \code{points} are interpreted as stratigraphic heights or as times depends on the option \code{direction}.
}
\item{xdep}{
Vector of strictly increasing real numbers 
}
  \item{ydep}{
Vector of real numbers. \code{xdep} and \code{ydep} describe the deposition model, and how they are interpreted depends on the option \code{depositionmodel} (and when \code{direction='height to time'}, also on the option \code{unit}). By default, \code{xdep} and \code{ydep} describe a piecewise linear deposition rate, i.e. \code{deporate=approxfun(xdep,ydep)}. For more options on how to describe deposition models, see \code{depositionmodel}
}
  \item{direction}{OPTIONAL, default is \code{'time to height'}. Determines the direction of the transformation. Either \code{'time to height'} (in which case the input given by \code{points} corresponds to times) or \code{'height to time'} (in which case the input of \code{points} corresponds to stratigraphic heights)
}
  \item{depositionmodel}{OPTIONAL, default is \code{'piecewise linear deposition rate'}. Either \code{'piecewise linear deposition rate'}, \code{'binned deposition rate'}, or \code{'age model'}. Determines how the input of \code{xdep} and \code{ydep} is interpreted. Dependent on the choice of \code{direction} and \code{depositionmodel}, different restrictions on the input of \code{ydep} exist, for details and examples see the vignette (available via \code{vignette('DAIME') or the error messages}
}
  \item{hiatuslist}{
OPTIONAL, default is an empty list. List of hiatuses to be included into the transformation, only used when \code{direction='height to time'}, e.g. the transformation is from stratigraphic height to time. A list, containing vectors of length two as elements. Every element in the list corresponds to a hiatus: \code{hiatuslist[[i]][1]} is the stratigraphic height of the i-th hiatus, and \code{hiatuslist[[i]][2]} is its duration
}
  \item{unit}{
OPTIONAL, default is 'sediment per time'. Only used when \code{direction='height to time'}, e.g. the transformation is from stratigraphic height to time, and \code{depositionmodel} is \code{'piecewise linear deposition rate'} or \code{'binned deposition rate'}. Either \code{'sediment per time'} or \code{'time per sediment'}. Determines the unit of the (binned or piecewise linear) deposition rate. If \code{unit='sediment per time'}, the deposition rate describes the deposition rate with which a given point in the section was deposited. If \code{unit='time per sediment'}, the deposition rate describes the time needed to deposit one unit of sediment (also know as the inverse deposition rate)
}
}

\value{
Returns a list containing:
\item{age }{Vector containing the times at which the stratigraphic heights in \code{height} were deposited}
\item{height }{Vector containing stratigraphic heights that were deposited at the times given in \code{age}}
\item{report }{A short summary of the task performed}

If \code{direction='time to height'}, \code{age} is a duplicate of the input \code{points}, and \code{height} contains the stratigraphic height that were deposited at the times given by \code{points}.

If \code{direction='height to time'}, \code{height} is a duplicate of the input \code{points}, and \code{age} contains the times at which the stratigraphic height given by \code{points} were deposited.

Output of \code{NA} in a vector indicates that the transformation was unsuccessful. This happens when entries of \code{points} coincide with a hiatus or are not in the domain of definition of the deposition model.
}

\references{
Hohmann, Niklas. 2018. Quantifying the Effects of Changing Deposition Rates and Hiatii on the Stratigraphic Distribution of Fossils. <doi:10.13140/RG.2.2.23372.51841>
}

\author{
Niklas Hohmann
}

\seealso{
For an overview of the functions in the DAIME package, see its vignette (available via \code{vignette('DAIME')} )

\code{\link{patterntransform}} for the transformation of input rates (stratigraphic patterns into temporal patterns and vice versa)

\code{\link{patterntodepositionmodel}} to create deposition models based on sedimentary dilution/condensation of known patterns

\code{\link{strattotimepointbin}}, \code{\link{strattotimepointcont}}, \code{\link{timetostratpointbin}}, and  \code{\link{timetostratpointcont}} for diverse wrappers of \code{pointtransform}. 
}
\examples{
##define deposition rate
xdep=seq(0,12,length.out=100)
ydep=splinefunH(x=c(0,2,4,6,8,10,12),y=c(1,5,6,1,0.5,1,6),m=c(0,1.5,-0.5,-0.5,0,0.5,0))(xdep)
#unit of deposition rate is sediment per time unit (default setting)
usedunit='sediment per time' #unit of deposition rate is sediment per time unit (default setting)
#Plot deposition rate (in stratigraphic height)
plot(xdep,ydep,type='l',main='Deposition Rate',xlab='Stratigraphic Height',ylab=usedunit,
  ylim=c(0,max(ydep)))

### Example 1: Determine age/stratigraphic height of single objects
##at what time was the object found at stratigraphic height 9 deposited?
pointtransform(9,xdep,ydep,direction='height to time', unit=usedunit) 
#change unit used
usedunit='time per sediment'
pointtransform(9,xdep,ydep,direction='height to time', unit=usedunit) 
#note how different the results are!
##Now, take deposition rate as deposition rate in time
##at what stratigraphic height will an object appear that was deposited at time 5?
pointtransform(5,xdep,ydep,direction='time to height') 

### Example 2: Create Age model based on a deposition rate
##create an age model. Essentially transform many points, which then approximate the age model
stratheight=seq(min(xdep),max(xdep),length.out=1000) #many points to approx. age model
usedunit='sediment per time'
reslist=pointtransform(stratheight,xdep,ydep,direction='height to time',unit=usedunit)
reslist$report
#plot age model
plot(reslist$age,reslist$height,xlab='Time',ylab='Stratigraphic Height',
  main=paste('Age model based on deposition rate with unit',usedunit))
#create age model but with other units for sedimentn input
usedunit='time per sediment'
reslist=pointtransform(stratheight,xdep,ydep,direction='height to time',unit=usedunit)
reslist$report
#plot age model (note the difference the setting of unit makes in terms of time 
#required to deposit the section!)
plot(reslist$age,reslist$height,xlab='Time',ylab='Stratigraphic Height',
  main=paste('Age model based on deposition rate with unit',usedunit))

##create age model with a hiatus 1: height to time
stratigraphicheight=5 #strat. height of the hiatus
duration=10 #duration of the hiatus
myhiatuslist=list(c(stratigraphicheight,duration)) #required input format for hiatuses
reslist=pointtransform(stratheight,xdep,ydep,direction='height to time',hiatuslist=myhiatuslist)
reslist$report
#!using default setting for unit (sediment per time) again!
plot(reslist$age,reslist$height,xlab='Time',ylab='stratigraphic height')
#the gap corresponds to the hiatus

##create age model with a hiatus 2: time to height
xdep=c(0,6,8,12)
ydep=c(0,6,5,12)
plot(xdep,ydep,type='l',main='Age Model, not eroded',xlab='time')
reslist=pointtransform(stratheight,xdep,ydep,direction='time to height',depositionmodel='age model')
reslist$report
plot(reslist$age,reslist$height,xlab='Time',ylab='stratigraphic height'
  ,main='Age model, eroded (with hiatus)')

###Example 3: Transform (isotope) ratios
##define deposition rate
xdep=seq(0,12,length.out=100)
ydep=splinefunH(x=c(0,2,4,6,8,10,12),y=c(1,5,6,1,0.5,1,6),m=c(0,1.5,-0.5,-0.5,0,0.5,0))(xdep)
#create fake (oxygen) isotope curves
samplelocation=seq(0,12,length.out=20) #where the samples are taken
isotoperatio=(-1)^(0:19) +1
plot(xdep,ydep,type='l',ylim=c(0,7),xlab='Stratigraphic Height',ylab='')
lines(samplelocation,isotoperatio,type='l',lwd=4)
legend('topleft',lwd=c(1,4),legend=c('Deposition rate','(Isotope) ratio'))
#transform only (!) sample locations, NOT values 
#again using the default setting for unit
reslist=pointtransform(samplelocation,xdep,ydep,direction='height to time')
#Isotope ratios in time 
plot(reslist$age,isotoperatio,type='l',xlab='Time',ylab='Isotope Ratio')
}

\keyword{ Deposition Rate }
\keyword{ Isotope }