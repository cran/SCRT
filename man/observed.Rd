\name{observed}

\alias{observed}

\title{Observed test statistic}

\description{The observed test statistic is calculated from the obtained raw data.}

\usage{observed(design, statistic, data = read.table(file.choose(new = FALSE)))}

\arguments{
  \item{design}{Type of single-case design: \code{"AB"}, \code{"ABA"}, \code{"ABAB"}, \code{"CRD"} (completely randomized design), \code{"RBD"} (randomized block design), \code{"ATD"} (alternating treatments design), \code{"MBD"} (multiple-baseline AB design) or \code{"Custom"} (user specified design).}

  \item{statistic}{Test statistic. For alternation designs, multiple-baseline designs and AB phase designs, there are 3 built-in possibilities: \code{"A-B"}, \code{"B-A"}, and \code{"|A-B|"}, which stand for the (absolute value of the) difference between condition means. 
For phase designs with more than 2 phases, 3 more built-in options are available: \code{"PA-PB"}, \code{"PB-PA"}, and \code{"|PA-PB|"} refer to the (absolute value of the) difference between the means of phase means.
Additionally, it is possible to specify a custom test statistic using the variable identifiers \code{"A"} and \code{"B"} (or in the case of phase deisgns with more than 2 phases, \code{"A1"}, \code{"B1"}, \code{"A2"}, \code{"B2"}, \code{"A"} and \code{"B"}) and any of the basic R functions. For example, \code{"abs(mean(A) - mean(B))"} can be used as a test statistic and it will be the same as using \code{"|A-B|"}.}

  \item{data}{File in which the data can be found. Default: a window pops up in which the file can be selected.}

}

\details{

When using the default \code{data} argument, a window will pop up to ask in what file the data can be found. This text file containing the data should consist of two columns for single-case phase and alternation designs: the first with the condition labels and the second with the obtained scores.

For multiple-baseline designs it should consist of these two columns for EACH unit. This way, each row represents one measurement occasion. It is important not to label the rows or columns.

Missing data should be indicated as \code{NA}. For calculations, missing data are omitted.

}

\references{

Bulte, I., & Onghena, P. (2008). An R package for single-case randomization tests. Behavior Research Methods, 40, 467-478.

Bulte, I., & Onghena, P. (2009). Randomization tests for multiple baseline designs: An extension of the SCRT-R package. Behavior Research Methods, 41, 477-485.

\url{http://ppw.kuleuven.be/home/english/research/mesrg}

}

\seealso{

\code{\link{distribution.systematic}} to generate the exhaustive randomization distribtion and 
\code{\link{pvalue.systematic}} to obtain the corresponding p-value.

\code{\link{distribution.random}} to generate the nonexhaustive randomization distribution and 
\code{\link{pvalue.random}} to obtain the corresponding p-value.

}

\examples{
data(ABAB)
observed(design = "ABAB", statistic = "PA-PB", data = ABAB)}

\keyword{Single-case design}
\keyword{Randomization test}
\keyword{Test statistic}

