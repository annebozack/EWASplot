% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/volcano_plot.R
\name{volcano_plot}
\alias{volcano_plot}
\title{Volcano plot}
\usage{
volcano_plot(
  pvector,
  effectvector,
  FDR = FALSE,
  title = NULL,
  col = "gray40",
  col.sig = "gray18",
  cex = 0.5,
  cex.sig = 0.75,
  alpha = 0.5,
  alpha.sig = 1,
  size.line.sig = 0.5
)
}
\arguments{
\item{pvector}{Vector of p-values}

\item{effectvector}{Vector of effect sizes}

\item{FDR}{Should an FDR level of significance be indicated? Defaults to FALSE.}

\item{title}{Title text.}

\item{col}{Color that will be used for probes for non-significant probes. Defaults to 'gray40'.}

\item{col.sig}{Color of FDR or Bonferroni significant probes. Defaults to 'gray18'.}

\item{cex}{Size of non-significant probes. Defaults to 0.5.}

\item{cex.sig}{Size of significant probes. Defaults to 0.5.}

\item{alpha}{Alpha of non-significant probes. Defaults to 0.5.}

\item{alpha.sig}{alpha Alpha of non-significant probes. Defaults to 1.}

\item{Size}{of line indicating FDR and Bonferroni significance. Defaults to 0.5.}
}
\value{
Volcano plot
}
\description{
This function creates a volcano plot of p-values and effect sizes using ggplot2.
}
