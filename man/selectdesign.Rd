\name{selectdesign}

\alias{selectdesign}

\title{Selection on one assignment possibility}

\description{One assignment possibility is randomly selected from all theoretical possibilities.}

\usage{selectdesign(design, MT, limit, starts = read.table(file.choose(new=FALSE)))}

\arguments{

  \item{design}{Type of single-case design: \code{"AB"}, \code{"ABA"}, \code{"ABAB"}, \code{"CRD"} (completely randomized design), \code{"RBD"} (randomized block design), \code{"ATD"} (alternating treatments design), or \code{"MBD"} (multiple-baseline AB design).}

  \item{MT}{Measurement times: number of observations.}

  \item{limit}{For phase designs: minimum number of observations per phase. 
For alternating treatments designs: maximum number of consecutive administrations of the same condition.}

  \item{starts}{Only for multiple baseline designs: location of the file where the possible start points can be found. Default: a window pops up in which the file can be selected.}

}

\details{

When using the default \code{starts} argument, only the \code{design} argument is required for multiple baseline designs. In this default version a window pops up in which is asked in what file the possible start points can be found. In this startpoint file, each row should contain all possibilities for one unit, separated by a tab. The rows and columns should not be labeled.

For multiple baseline designs, a possible combination of start points for each unit is returned. For all other designs, a possible sequence of conditions is returned (e.g., "A" "A" "A" "A" "B" "B" "B"). 

}

\references{

Bulte, I., & Onghena, P. (2008). An R package for single-case randomization tests. Behavior Research Methods, 40, 467-478.

Bulte, I., & Onghena, P. (2009). Randomization tests for multiple baseline designs: An extension of the SCRT-R package. Behavior Research Methods, 41, 477-485.

\url{http://ppw.kuleuven.be/english/research/mesrg}

}

\seealso{

\code{\link{quantity}} to calculate the number of assignment possibilities.

\code{\link{assignments}} to enumerate all assignment possibilities.

}

\examples{selectdesign(design = "ABAB", MT = 24, limit = 4)}

\keyword{Single-case design}
\keyword{Randomization test}
