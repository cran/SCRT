\name{graph1}

\alias{graph1}

\title{Graphical representation of single-case data}

\description{The observed single-case data are plotted.}

\usage{graph1(design,data=read.table(file.choose(new=FALSE)),
xlab="Measurement Times",ylab="Scores")}

\arguments{

  \item{design}{Type of single-case design: \code{"AB"}, \code{"ABA"}, \code{"ABAB"}, \code{"CRD"} (completely randomized design), \code{"RBD"} (randomized block design), \code{"ATD"} (alternating treatments design), or \code{"MBD"} (multiple-baseline AB design).}

  \item{data}{File in which the data can be found. Default: a window pops up in which the file can be selected.} 
  
  \item{xlab}{Label x axis}
  
  \item{ylab}{Label y axis}

}

\details{

When using the default \code{data} argument, a window will pop up to ask in what file the data can be found. This text file containing the data should consist of two columns for single-case phase and alternation designs: the first with the condition labels and the second with the obtained scores.
For multiple-baseline designs it should consist of these two columns for EACH unit. This way, each row represents one measurement occasion. It is important not to label the rows or columns.

For alternation designs, after the plot is drawn, the location of the legend should be indicated by a left mouse click.

}

\references{

Bulte, I., & Onghena, P. (2008). An R package for single-case randomization tests. Behavior Research Methods, 40, 467-478.

Bulte, I., & Onghena, P. (2009). Randomization tests for multiple baseline designs: An extension of the SCRT-R package. Behavior Research Methods, 41, 477-485.

Bulte, I., & Onghena, P. (in press). When the Truth Hits You Between the Eyes: A Software Tool for the Visual Analysis of Single-Case Experimental Data. Manuscript accepted for publication in Methodology. 

\url{http://ppw.kuleuven.be/english/research/mesrg}

}

\examples{
data(ABAB)
graph1(design = "ABAB", data = ABAB)}

\keyword{Single-case design}
\keyword{Graph}