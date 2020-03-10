\name{ColDePalluel.rdata}
\alias{ColDePalluel}
\docType{data}
\title{Microcpaleontological, palynological, and geochemical data from the Col de Palluel, SE France (Late Albanian)}
\description{
Contains data on geochemistry, palynology, planktic forams, and nannofossils from the Col de Palluel from Bornemann et al. (2005)
}
\usage{data(ColDePalluel)}
\format{A list with four entries named "geochemistry", "palynomorphs", "nannofossils", and "plankticforaminifera" containing the data given in Bornemann et al. (2005).
}

\examples{
data("ColDePalluel")
#plot carbonate content
plot(ColDePalluel$geochemistry$Depth..m.,ColDePalluel$geochemistry$CaCO3....,type='l')
}
\source{
Bornemann, Andre; Pross, Joerg; Reichelt, Kerstin; Herrle, Jens O; Hemleben, Christoph; Mutterlose, Joerg (2005): Micropaleontological investigation of Late Albanian Col de Palluel. PANGAEA, https://doi.org/10.1594/PANGAEA.737945 \href{https://creativecommons.org/licenses/by/3.0/}{CC BY 3.0}
}
\references{
Bornemann, A et al. (2005): Reconstruction of short-term palaeoceanographic changes during the formation of the Late Albian Niveau Breistroffer black shales (Oceanic Anoxic Event 1d, SE France). Journal of the Geological Society of London, 162(4), 623-639, https://doi.org/10.1144/0016-764903-171
}
\keyword{dataset}
\keyword{Albanian}
\keyword{Micropaleontology}
\keyword{Stratigraphy}