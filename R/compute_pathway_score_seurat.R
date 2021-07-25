library(Seurat)

#' Run pathway score pipeline
#' @param object Pipeline object
#' @export
#' @rdname sc.pathway.score
sc.pathway.score <- function(object,...){
  UseMethod("sc.pathway.score")
}


#' @rdname sc.pathway.score
#' @return Seurat (version 3) object.
sc.pathway.score.Seurat <- function(object,
                                    assay.use = NULL,
                                    slot = "counts",
                                    method = "ssGSEA",
                                    gmtFile,
                                    ncores = 1,
                                    metabolism.type = "KEGG"
                                    ){
  assay.use <- assay.use %||% DefaultAssay(object)
  expr <- GetAssayData(object = object, slot = slot, assay = assay.use)
  expr<-data.frame(as.matrix(expr))


  signatures_KEGG_metab <- system.file("data", "KEGG_metabolism_nc.gmt", package = "sc.utils")
  signatures_REACTOME_metab <- system.file("data", "REACTOME_metabolism.gmt", package = "sc.utils")


  if (metabolism.type == "KEGG")  {gmtFile<-signatures_KEGG_metab; cat("Your choice is: KEGG\n")}
  if (metabolism.type == "REACTOME")  {gmtFile<-signatures_REACTOME_metab; cat("Your choice is: REACTOME\n")}

  #----------------calculate score--------------------------------

  # for ssGSEA or GSVA, we recommend to use the whole gene expression matrix rather than matrix only containing highly variable genes
  #ssGSEA
  if (method == "ssgsea") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile) #signature read
    gsva_es <- gsva(as.matrix(expr), geneSets, method=c("ssgsea"), kcdf=c("Poisson"), parallel.sz=ncores) #
    signature_exp<-data.frame(gsva_es)

  }

  #GSVA
  if (method == "gsva") {
    library(GSVA)
    library(GSEABase)
    geneSets <- getGmt(gmtFile) #signature read
    gsva_es <- gsva(as.matrix(expr), geneSets, method=c("gsva"), kcdf=c("Poisson"), parallel.sz=ncores) #
    signature_exp<-data.frame(gsva_es)
  }

  object@assays$pathway_score$score<-signature_exp

  return(object)
}
