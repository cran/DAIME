\name{SeymourIslandAgeModels}
\alias{SeymourIslandAgeModels}
\docType{data}
\title{Age models for the K/Pg boundary on Seymour Island, Antarctica}
\description{
Contains Age models for K/Pg boundary on Seymour Island, Antarctica. Age models are based

*On the Stratigraphic height of the K/Pg boundary of 1007.5 m with an age of 66.04 ma (Renne et al. 2013)

*On magnetostratigraphy from Tobin et al. (2012) and ages from Gradstein et al. (2012), yielding an age of 66.389 ma at the stratigraphic height of 934.4 m.

The stratigraphic height of a low deposition interval is assumed to be at 959.5 m based on section A from Witts et al. (2016).
}
\usage{data(SeymourIslandBins)}
\format{A list with entries, each named after the age model described in Hohmann (in preparation). Each of these age models is divided into
An list containing two vecotrs 'agemodelage' and 'agemodelheight' corresponding to the time of deposition of a given stratigraphic height
A description of the paraameters used to generate the age models
}

\examples{
data("SeymourIslandAgeModels")
#Select age model B
my.agemodel=SeymourIslandAgeModels$B
#check parameters of age model
my.agemodel$parameters
plot(my.agemodel$age,my.agemodel$height,type='l')
}
\source{
Gradstein, Felix M., et al., eds. The geologic time scale 2012. elsevier, 2012.

Renne, Paul R., et al. "Time scales of critical events around the Cretaceous-Paleogene boundary." Science 339.6120 (2013): 684-687.

Tobin, Thomas S., et al. "Extinction patterns, delta18 O trends, and magnetostratigraphy from a southern high-latitude Cretaceous/Paleogene section: Links with Deccan volcanism." Palaeogeography, Palaeoclimatology, Palaeoecology 350 (2012): 180-188.

Hohmann, Niklas. 2018. Quantifying the Effects of Changing Deposition Rates and Hiatii on the Stratigraphic Distribution of Fossils. <doi:10.13140/RG.2.2.23372.51841>.

Witts, James D., et al. "Macrofossil evidence for a rapid and severe Cretaceous/Paleogene mass extinction in Antarctica." Nature communications 7 (2016): 11738.
}
\references{  None}
\keyword{dataset}
\keyword{Seymour}
\keyword{K/Pg}
\keyword{Stratigraphy}