\name{quantity}

\alias{quantity}

\title{Number of assignment possibilities}

\description{The number of assignment possibilities for the specified design is calculated.}

\usage{quantity(design, MT, limit, starts = file.choose(new = FALSE), 
assignments = file.choose(new = FALSE))}

\arguments{

  \item{design}{Type of single-case design: \code{"AB"}, \code{"ABA"}, \code{"ABAB"}, \code{"CRD"} (completely randomized design), \code{"RBD"} (randomized block design), \code{"ATD"} (alternating treatments design), \code{"MBD"} (multiple-baseline AB design) or \code{"Custom"} (user specified design).}

  \item{MT}{Measurement times: number of observations.}

  \item{limit}{For phase designs: minimum number of observations per phase. 
For alternating treatments designs: maximum number of consecutive administrations of the same condition.}

  \item{starts}{Only for multiple baseline designs: location of the file where the possible start points can be found. Default: a window pops up in which the file can be selected.}
  
  \item{assignments}{Only for user specified designs: location of the file where all the possible assignments can be found. Default: a window pops up in which the file can be selected.}

}

\details{

For multiple baseline designs, when using the default \code{starts} argument, only the \code{design} argument is required. In this default version a window pops up in which is asked in what file the possible start points can be found. In this startpoint file, each row should contain all possibilities for one unit, separated by a tab. The rows and columns should not be labeled.

The number of possible assignments for a multiple baseline design is calculated assuming NO overlap between the possible start points of the different units (staggered administration). 

For user specified designs, when using the default \code{assignments} argument, a window pops up in which is asked in what file all the possible assignments can be found. In this file, each row should contain the sequence of conditions in one possible assignment, separated by a tab. There should be one row for every possible assignment. The rows and columns should not be labeled.

}

\references{

Bulte, I., & Onghena, P. (2008). An R package for single-case randomization tests. Behavior Research Methods, 40, 467-478.

Bulte, I., & Onghena, P. (2009). Randomization tests for multiple baseline designs: An extension of the SCRT-R package. Behavior Research Methods, 41, 477-485.

\url{http://ppw.kuleuven.be/home/english/research/mesrg}

}

\seealso{

\code{\link{assignments}} to enumerate all assignment possibilities.

\code{\link{selectdesign}} to randomly select one of the assignment possibilities.

}

\examples{quantity(design = "ABAB", MT = 24, limit = 4)}

\keyword{Single-case design}
\keyword{Randomization test}