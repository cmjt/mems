% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mems.r
\name{metrics}
\alias{metrics}
\title{A range of triangle quality metrics}
\usage{
metrics(a, b, c)
}
\arguments{
\item{a}{A numeric triangle edgelength (opposite angle A).}

\item{b}{A numeric triangle edgelength (opposite angle B).}

\item{c}{A numeric triangle edgelength (opposite angle C).}
}
\value{
A list of triangle 'quality' metrics collated in \url{https://doi.org/10.1016/j.finel.2019.01.010}
\itemize{
\item \code{metric_1}: \eqn{\frac{4\sqrt{3}|A|}{\Sigma_{i = 1}^3 L_i^2}},
where \eqn{L_i} is the length of edge \eqn{i}.
\item \code{metric_2}: \eqn{6 \sqrt{\frac{A}{\sqrt{3}}}}.
\item \code{metric_3}:\eqn{\frac{\text{min}(L_i)}{\text{max}(L_i)}}.
\item \code{metric_4}:  \eqn{3 \frac{\text{min}(\phi_i)}{\pi}}, where
\eqn{\phi_i} is the \eqn{i^{th}} internal angle.
\item \code{metric_5}: \eqn{\frac{4\sqrt{3} A}{\text{max}(L_i)} * \Sigma_{i = 1}^3 L_i^2}.
\item \code{metric_6}: \eqn{\frac{1}{4 q_b} + q_w \frac{\sqrt{3}}{16}}, where
\eqn{q_b = \frac{\Sigma_{i = 1}^3 L_i^2}{4 \sqrt{3} A}} and
\eqn{q_m = (\frac{1}{3A})(\Sigma_{i = 1}^3 L_i)^2}.
}
}
\description{
A range of triangle quality metrics
}