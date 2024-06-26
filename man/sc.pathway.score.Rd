% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_pathway_score_seurat.R
\name{sc.pathway.score}
\alias{sc.pathway.score}
\alias{sc.pathway.score,Seurat-method}
\alias{sc.pathway.score.Seurat}
\title{Run single cell pathway score pipeline}
\usage{
sc.pathway.score(object, ...)

\S4method{sc.pathway.score}{Seurat}(
  object,
  assay.use = NULL,
  slot = "counts",
  method = "ssGSEA",
  gmtFile = NULL,
  ncores = 4,
  msigdb = FALSE,
  msigdb.species = "Homo sapiens",
  metabolism.geneset = TRUE,
  metabolism.type = c("KEGG", "REACTOME")
)

sc.pathway.score.Seurat(
  object,
  assay.use = NULL,
  slot = "counts",
  method = c("ssGSEA", "gsva"),
  gmtFile = NULL,
  ncores = 4,
  msigdb = FALSE,
  msigdb.species = "Homo sapiens",
  metabolism.geneset = TRUE,
  metabolism.type = "KEGG"
)
}
\arguments{
\item{object}{Seurat object}
}
\value{
Seurat (version 3) object.
}
\description{
Run single cell pathway score pipeline

Run pathway score pipeline with seurat object
}
