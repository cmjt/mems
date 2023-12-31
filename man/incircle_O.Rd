% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mems.r
\name{incircle_O}
\alias{incircle_O}
\title{Calculate incumcircle radius given vertex A, B, C
locations of a triangle}
\usage{
incircle_O(A, B, C)
}
\arguments{
\item{A}{A numeric vector of length 2 specifying vertex location "A".}

\item{B}{A numeric vector of length 2 specifying vertex location "B".}

\item{C}{A numeric vector of length 2 specifying vertex location "C".}
}
\description{
Calculate incumcircle radius given vertex A, B, C
locations of a triangle
}
\examples{
incircle_O(c(1, 1), c(2, 1), c(3, 2.5))
}
