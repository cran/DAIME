\name{patterntodepositionmodel}
\alias{patterntodepositionmodel}

\title{
Use sedimentary dilution/condensation of known patterns to create age model
}

\description{
Takes a pair of a stratigraphic and a temporal pattern to determine the deposition model that transforms them into each other. The deposition model is returned as an age model.
}

\usage{
patterntodepositionmodel(xheight, yheight, 
    xage = NULL, yage = NULL,
    heightmode = 'piecewise linear', agemode = 'piecewise linear',
    atheight=NULL,atage=NULL,rescalefor=1,timetype='time')
}

\arguments{
\item{xheight}{A vector of strictly increasing numbers
}
\item{yheight}{A vector of strictly positive numbers. \code{xheight} and \code{yheight} define the stratigraphic pattern. By default, it is assumed to be piecewise linear and is accordingly given by \code{approxfun(xheight,yheight)}. This can be changed to a binned representation using the option \code{heightmode} (see below)
}
\item{xage}{OPTIONAL, default is \code{NULL}
}
\item{yage}{OPTIONAL, default is \code{NULL}. \code{xage} and \code{yage} can be used to define the temporal pattern used for the reconstruction of the deposition model. If no input is handed over (default setting), it is assumed that the temporal pattern is constant over a time interval of duration 1. Using the option \code{agemode} allows to define both piecewise linear and binned temporal patterns (see below)
}
\item{heightmode}{OPTIONAL, default is \code{'piecewise linear'}. Either \code{'piecewise linear'} or \code{'binned'}. Determines whether \code{xheight} and \code{yheight} are taken as a piecewise linear or a binned description of the stratigraphic pattern. In the first case, \code{length(xheight)} needs to match \code{length(yheight)}, and the stratigraphic pattern is given by \code{approxfun(xheight,yheight)}. In the second case, \code{length(xheight)} needs to match \code{length(yheight)+1}, and the value of the stratigraphic pattern between the stratigraphic heights \code{xheight[i]} and \code{xheight[i+1]} is given by \code{yheight[i]}
}
\item{agemode}{OPTIONAL, default is \code{'piecewise linear'}. Either \code{'piecewise linear'} or \code{'binned'}. Determines whether \code{xage} and \code{yage} are taken as a piecewise linear or a binned description of the temporal pattern. In the first case, \code{length(xage)} needs to match \code{length(yage)}, and the temporal pattern is given by \code{approxfun(xage,yage)}. In the second case, \code{length(xage)} needs to match \code{length(yage)+1}, and the value of the temporal pattern between the ages \code{xage[i]} and \code{xage[i+1]} is given by \code{yage[i]}
}
\item{atheight}{OPTIONAL, default is \code{NULL}. Stratigraphic heights that should be included in the age model. By default, the function chooses those heights automatically
}
\item{atage}{OPTIONAL, default is \code{NULL}. Ages that should be included in the age model. By default, the function chooses those ages automatically
}
\item{rescalefor}{OPTIONAL, default is 1. Either a strictly positive number, 'temporal pattern', or 'stratigraphic pattern'. Determines what the total input into the sediment is.}
\item{timetype}{OPTIONAL, default is "time". Either "time" or "age", determines whether input/output will be interpreted/given as time or age
}
}

\value{
If \code{timetype='time'}, a list containing the following entries:
\item{time }{Vector of times}
\item{height }{Vector of stratigraphic heights}
\item{report }{A short summary of the task performed}

If \code{timetype='age'}, a list containing the following entries:
\item{age }{Vector of ages}
\item{height }{Vector of stratigraphic heights}
\item{report }{A short summary of the task performed}

\code{age}/\code{time} and \code{height} form the age model in the sense that the age/time at which \code{height[i]} was deposited is given by \code{age[i]}/\code{time[i]}. Conversely at the age/time \code{age[i]}/\code{time[i]}, \code{height[i]} was deposited
}

\references{
Hohmann, Niklas. 2018. Quantifying the Effects of Changing Deposition Rates and Hiatii on the Stratigraphic Distribution of Fossils. <doi:10.13140/RG.2.2.23372.51841>
}

\author{
Niklas Hohmann
}

\seealso{
\code{\link{pointtransform}} for the transformations of points using age models

\code{\link{patterntransform}} for the transformations of temporal and stratigraphic patterns (such as input rates into the sediment) using age models

For an overview of the functions in the DAIME package and examples using stratigraphic data see the vignette (available via \code{vignette('DAIME')} )
}

\examples{
### Reconstruct deposition model based on condensation/dilution of a constant input
#assume a constant input of pollen (this is an arbitrary choice) into the sediment through time
#this is the observed stratigraphic pattern
my.xheight=seq(0,10,by=1)
my.yheight=seq(0,3,length.out=11)+rexp(11)
plot(my.xheight,my.yheight,xlab='Stratigraphic Height',ylab='Pollen Abundance',type='l',
  ylim=c(0,max(my.yheight)),main='Stratigraphic Pattern')
#reconstruct deposition model based on the assumption of constant pollen input in time
my.agemodel=patterntodepositionmodel(xheight=my.xheight,yheight=my.yheight)
my.agemodel$report
agemodelage=my.agemodel$time
agemodelheight=my.agemodel$height
plot(agemodelage,agemodelheight,type='l',
  lwd=6,main='Reconstructed Deposition Models as Age Model')
legend('topleft',legend='Age Model',lwd=6,col='black')
#approximate deposition rate (=derivative of the age model)
grad=diff(agemodelheight)/diff(agemodelage)
xbase=agemodelage[-1]
plot(xbase,grad,xlab='Time',ylab='Deposition Rate',
  main='Reconstructed Deposition Model as Deposition Rate',type='l',lwd=6,ylim=c(0,max(grad)))

#now assume pollen input into the sediment is decreasing and lasts for 2 time units
my.xage=c(0,2)
my.yage=c(5,1)
plot(my.xage,my.yage,type='l',xlab='Time',ylab='Pollen Input',ylim=c(0,max(my.yage)),
  lwd=6,main='Temporal Pattern')
#reconstruct age model based on these updated assumptions
my.agemodel2=patterntodepositionmodel(xheight=my.xheight,yheight=my.yheight,
    xage=my.xage,yage=my.yage)
my.agemodel2$report
agemodelage2=my.agemodel2$time
agemodelheight2=my.agemodel2$height
plot(agemodelage2,agemodelheight2,type='l',lwd=6, 
  main='Reconstructed Deposition Model as Age Model')
legend('topleft',legend='Age Model',lwd=6,col='black')

#if a pattern is given as bins, use the option 'heightmode' or 'agemode'
#define stratigraphic pattern as binned function
my.xheight3=seq(0,10,length.out=11) #borders of the bins used to define the stratigraphic pattern
my.yheight3=seq(0,3,length.out=10)+rexp(10) #note that xheight has one element more than yheight
barplot(my.yheight3,width=diff(my.xheight3),ylab='Pollen Abundance',xlab='Stratigraphic Height',
  space=0,main='Stratigraphic Pattern')
#reconstruct age model as in the first example, but with the binned pollen observation 
#  in stratigraphic height
#Note that pollen input in time is again assumed to be constant 
#(default setting if no xage and yage given)
my.agemodel3=patterntodepositionmodel(xheight=my.xheight3,yheight=my.yheight3,heightmode='binned')
my.agemodel3$report
agemodelage3=my.agemodel3$time
agemodelheight3=my.agemodel3$height
plot(agemodelage3,agemodelheight3,type='l',lwd=6,
  main='Reconstructed Deposition Model as Age Model')
legend('topleft',legend='Age Model',lwd=6,col='black')
#The same option is also available for the temporal pattern (no example given)

}

\keyword{ Deposition rate }
\keyword{ Stratigraphy }
