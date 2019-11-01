\name{patterntodepositionmodel}
\alias{patterntodepositionmodel}

\title{
Use sedimentary dilution/condensation of known patterns to create age models
}

\description{
Takes patterns (in stratigraphic height and/or time) to determine the deposition model that transforms them into each other. The deposition model is returned as an age model.
}

\usage{
patterntodepositionmodel(xheight, yheight, 
    xage = NULL, yage = NULL,
    heightmode = 'piecewise linear', agemode = 'piecewise linear')
}

\arguments{
\item{xheight}{A vector of strictly increasing numbers
}
\item{yheight}{A vector of strictly positive numbers. \code{xheight} and \code{yheight} define the stratigraphic pattern. By default, it is assumed to be piecewise linear and is accordingly given by \code{approxfun(xheight,yheight)}. This can be changed using the option \code{heightmode} ( see below)
}
\item{xage}{OPTIONAL, default is \code{NULL}
}
\item{yage}{OPTIONAL, default is \code{NULL}. \code{xage} and \code{yage} can be used to define the temporal pattern used for the reconstruction of the deposition model. If no input is handed over (default setting), it is assumed that the temporal pattern is constant over a time interval of duration 1. Using the option \code{agemode} allows to define both piecewise linear and binned temporal patterns (see below)
}
\item{heightmode}{OPTIONAL, default is \code{'piecewise linear'}. Either \code{'piecewise linear'} or \code{'binned'}. Determines whether \code{xheight} and \code{yheight} are taken as a piecewise linear or a binned description of the stratigraphic pattern. In the first case, \code{length(xheight)} needs to match \code{length(yheight)}, and the stratigraphic pattern is given by \code{approxfun(xheight,yheight)}. In the second case, \code{length(xheight)} needs to match \code{length(yheight)+1}, and the value of the stratigraphic pattern between the stratigraphic heights \code{xheight[i]} and \code{xheight[i+1]} is given by \code{yheight[i]}
}
\item{agemode}{OPTIONAL, default is \code{'piecewise linear'}. Either \code{'piecewise linear'} or \code{'binned'}. Determines whether \code{xage} and \code{yage} are taken as a piecewise linear or a binned description of the temporal pattern. In the first case, \code{length(xage)} needs to match \code{length(yage)}, and the temporal pattern is given by \code{approxfun(xage,yage)}. In the second case, \code{length(xage)} needs to match \code{length(yage)+1}, and the value of the temporal pattern between the ages \code{xage[i]} and \code{xage[i+1]} is given by \code{yage[i]}.}
}

\value{
A list containing the following entries:
\item{age }{Vector of ages}
\item{height }{Vector of stratigraphic heights}
They form the age model in the sense that the age at which \code{height[i]} was deposited is given by \code{age[i]}. Conversely at the time \code{age[i]}, \code{height[i]} was deposited.
}

\references{
Hohmann, Niklas. 2018. Quantifying the Effects of Changing Deposition Rates and Hiatii on the Stratigraphic Distribution of Fossils. <doi:10.13140/RG.2.2.23372.51841>
}

\author{
Niklas Hohmann
}

\seealso{
\code{\link{pointtransform}} for the transformations of points using age models

\code{\link{patterntransform}} for the transformations of patterns (such as input rates into the sediment) using age models

For an overview of the functions in the DAIME package, see its vignette (available via \code{vignette('DAIME')} )
}

\examples{
### Reconstruct deposition model based on condensation/dilution of a constant input
#assume a constant input of pollen (this is an arbitrary choice) into the sediment through time
#this is the observed stratigraphic pattern
xheight=seq(0,10,by=1)
yheight=seq(0,3,length.out=11)+rexp(11)
plot(xheight,yheight,xlab='Stratigraphic Height',ylab='Pollen Abundance',type='l',
  ylim=c(0,max(yheight)),main='Stratigraphic Pattern')
#reconstruct deposition model based on the assumption of constant pollen input in time
resultlist=patterntodepositionmodel(xheight,yheight)
plot(resultlist$age,resultlist$height,xlab='Age',ylab='Stratigraphic Height',type='l',
  lwd=6,main='Reconstructed Deposition Models as Age Model')
legend('topleft',legend='Age Model',lwd=6,col='black')
#approximate deposition rate (=derivative of the age model)
grad=diff(resultlist$height)/diff(resultlist$age)
xbase=resultlist$age[-1]
plot(xbase,grad,xlab='Time',ylab='Deposition Rate',
  main='Reconstructed Deposition Model as Deposition Rate',type='l',lwd=6,ylim=c(0,max(grad)))

#now assume pollen input into the sediment is decreasing and lasts for 2 time units
xage=c(0,2)
yage=c(5,1)
plot(xage,yage,type='l',xlab='Time',ylab='Pollen Input',ylim=c(0,max(yage)),
  lwd=6,main='Temporal Pattern')
#reconstruct age model based on these updated assumptions
resultlist=patterntodepositionmodel(xheight,yheight,xage,yage)
plot(resultlist$age,resultlist$height,xlab='Age',ylab='Stratigraphic Height',type='l',lwd=6, 
  main='Reconstructed Deposition Model as Age Model')
legend('topleft',legend='Age Model',lwd=6,col='black')

#if a pattern is given as bins, use the option 'heightmode' or 'agemode'
#define stratigraphic pattern as binned function
xheight=seq(0,10,length.out=11) #borders of the bins used to define the stratigraphic pattern
yheight=seq(0,3,length.out=10)+rexp(10) #note that xheight has one element more than yheight
barplot(yheight,width=diff(xheight),ylab='Pollen Abundance',xlab='Stratigraphic Height',
  space=0,main='Stratigraphic Pattern')
#reconstruct age model as in the first example, but with the binned pollen observation 
#  in stratigraphic height
#Note that pollen input in time is again assumed to be constant 
#(default setting if no xage and yage given)
resultlist=patterntodepositionmodel(xheight,yheight,heightmode='binned')
plot(resultlist$age,resultlist$height,xlab='Age',ylab='Stratigraphic Height',type='l',lwd=6,
  main='Reconstructed Deposition Model as Age Model')
legend('topleft',legend='Age Model',lwd=6,col='black')
#The same option is also available for the temporal pattern (no example given)
}

\keyword{ Deposition rate }
\keyword{ Stratigraphy }
