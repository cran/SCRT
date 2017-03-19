\name{assignments}

\alias{assignments}

\title{All assignments possibilities}

\description{All assignment possibilities for the specified design are enumerated.}

\usage{assignments(design, save = "no", MT, limit, 
starts = file.choose(new = FALSE), assignments = file.choose(new = FALSE))}

\arguments{

  \item{design}{Type of single-case design: \code{"AB"}, \code{"ABA"}, \code{"ABAB"}, \code{"CRD"} (completely randomized design), \code{"RBD"} (randomized block design), \code{"ATD"} (alternating treatments design), \code{"MBD"} (multiple-baseline AB design) or \code{"Custom"} (user specified design).}

  \item{save}{Save the possible assignments to a file (\code{save="yes"}) or just see them as output in the R console (default: \code{save="no")}.}

  \item{MT}{Measurement times: number of observations.}

  \item{limit}{For phase designs: minimum number of observations per phase. 
For alternating treatments designs: maximum number of consecutive administrations of the same condition.}

  \item{starts}{Only for multiple baseline designs: location of the file where the possible start points can be found. Default: a window pops up in which the file can be selected.}
  
  \item{assignments}{Only for user specified designs: location of the file where all the possible assignments can be found. Default: a window pops up in which the file can be selected.}

}

\details{

When choosing to save the possible assignments to a file, a window will pop up (for multiple baseline designs or user specified designs this is the second pop-up window!!!) to ask where to save them. This location can be an existing file, as well as a new file that can be created by giving a file name and the extension .txt. In this latter case a confirmation is required ("The file does not exist yet. Create the file?"). 

For multiple baseline designs, when using the default \code{starts} argument, first a window pops up in which is asked in what file the possible start points can be found. In this startpoint file, each row should contain all possibilities for one unit, separated by a tab. The rows and columns should not be labeled.

For user specified designs, when using the default \code{assignments} argument, first a window pops up in which is asked in what file all the possible assignments can be found. In this file, each row should contain the sequence of conditions in one possible assignment, separated by a tab. There should be one row for every possible assignment. The rows and columns should not be labeled.

For multiple baseline designs, the possible combinations of start points for each unit are returned. For all other designs, the possible sequences of conditions are returned (e.g., "A" "A" "A" "A" "B" "B" "B"). 

}

\references{

Bulte, I., & Onghena, P. (2008). An R package for single-case randomization tests. Behavior Research Methods, 40, 467-478.

Bulte, I., & Onghena, P. (2009). Randomization tests for multiple baseline designs: An extension of the SCRT-R package. Behavior Research Methods, 41, 477-485.

\url{http://ppw.kuleuven.be/english/research/mesrg}

}

\seealso{

\code{\link{quantity}} to calculate the number of assignment possibilities.

\code{\link{selectdesign}} to randomly select one of the assignment possibilities.

}

\examples{assignments(design = "ABAB", save = "no", MT = 24, limit = 4)}

\keyword{Single-case design}
\keyword{Randomization test}
