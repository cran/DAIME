\name{SeymourIslandBins.rdata}
\alias{SeymourIslandBins}
\docType{data}
\title{Sampling bins approaching and covering the K/Pg boundary on Seymour Island, Antarctica}

\description{
Contains the sampling bins from Seymour Island from
\itemize{
\item Macellari (1984), section D 
\item Witts et al. (2016), section A
}
Reference for stratigraphic heights is the composite section by Witts et al. (2016) with the K/Pg boundary at 1007.5 m. The bins from Macellari were correlated using the glauconite bed and the last ammonite occurrences at the K/PG boundary (Macellari 1984, p. 506 and 518, samples 401,402, and 403).
}
\usage{data(SeymourIslandBins)}

\format{A list with two entries: One corresponding to the bins from Macellari (1984), the other to the bins from Witts et al. (2016)
}
\examples{
data("SeymourIslandBins")
#Compare bins from Macellari and Witts before the K/Pg extinction
plot(SeymourIslandBins$Macellari.1984.Section.D.,
  rep(0,length(SeymourIslandBins$Macellari.1984.Section.D.)),
  ylim=c(0,1.1),xlab='Stratigraphic Height [m]',xlim=c(930,1007.5))
points(SeymourIslandBins$Witts.et.al.2016.Section.A,
  rep(0.2,length(SeymourIslandBins$Witts.et.al.2016.Section.A)),pch=19)
lines(c(1007.5,1007.5),c(0,2))
legend('topleft',pch=c(1,19,NA),lty=c(NA,NA,1),legend=c('Macellari','Witts','K/Pg boundary'))
}

\source{
Macellari, Carlos E. Late Cretaceous stratigraphy, sedimentology, and macropaleontology of Seymour Island, Antarctic Peninsula. Diss. The Ohio State University, 1984.

Witts, James D., et al. "Macrofossil evidence for a rapid and severe Cretaceous/Paleogene mass extinction in Antarctica." Nature communications 7 (2016): 11738.
}

\keyword{dataset}
\keyword{Seymour}
\keyword{K/Pg}
\keyword{Stratigraphy}
