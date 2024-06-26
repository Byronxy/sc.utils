% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cell_markers.R
\name{CMenrich}
\alias{CMenrich}
\title{CellMarker Enrichment}
\usage{
CMenrich(gene.list, species)
}
\arguments{
\item{gene.list}{character vector with a list of \bold{Gene Symbols}}

\item{species}{Can either be \emph{human} or \emph{mouse} for the moment.}
}
\value{
A \emph{list} two elememts. First element (enrichments) is a data.frame structured as follows:
\itemize{
\item {\bold{p.value}} {P-value of the enrichment, the lower the better}
\item {\bold{overlap}} {How many common genes between your costom list and the markers of the given cell type}
\item {\bold{signature}} {Total markers of the given cell type }
\item {\bold{genes}} {Common genes}
}
The second element is a list with the actual gene names for each enrichment.
}
\description{
Given a list of markers, indentifies the putative cell type using the comphensive catlog creatd in
}
\examples{
results=CMenrich(gene.list=c('NEUROD1','ARX','CD79A'),species='human') .# some gene names that I came up with ....

}
