#' @title Run single cell pathway score pipeline
#' @rdname sc.pathway.score
#' @param object Pipeline object
#' @export
setGeneric("sc.pathway.score", function(object,
                                    ...)
  standardGeneric("sc.pathway.score"))

#' @rdname sc.pathway.score
#' @exportMethod sc.pathway.score
setMethod("sc.pathway.score", signature(object = "Seurat"),
          function(object,
                   assay.use = NULL,
                   slot = "counts",
                   method = "ssGSEA",
                   gmtFile = NULL,
                   ncores = 4,
                   msigdb = FALSE,
                   msigdb.species = "Homo sapiens",
                   metabolism.geneset = TRUE,
                   metabolism.type = c("KEGG","REACTOME")) {
            sc.pathway.score.Seurat (object = object,
                                     assay.use = assay.use,
                                     slot = slot,
                                     method = method,
                                     gmtFile = gmtFile,
                                     ncores = ncores,
                                     msigdb = msigdb,
                                     msigdb.species = msigdb.species,
                                     metabolism.geneset = metabolism.geneset,
                                     metabolism.type = metabolism.type)
          })

#' @title Run pathway score pipeline with seurat object
#' @rdname sc.pathway.score
#' @param object Seurat object
#' @importFrom SeuratObject GetAssayData Reductions Embeddings FetchData
#' @return Seurat (version 3) object.
#' @export
#'
sc.pathway.score.Seurat <- function(object,
                   assay.use = NULL,
                   slot = "counts",
                   method = c("ssGSEA","gsva"),
                   gmtFile = NULL,
                   ncores = 4,
                   msigdb = FALSE,
                   msigdb.species = "Homo sapiens",
                   metabolism.geneset = TRUE,
                   metabolism.type = "KEGG") {

            assay.use <- assay.use %||% DefaultAssay(object)

            expr <- SeuratObject::GetAssayData(object = object, slot = slot, assay = assay.use)

            expr<-data.frame(as.matrix(expr))

            #print(expr[1:5,1:5])
            if (msigdb == TRUE) {
              library(msigdbr)
              gmt <- msigdbr(species = msigdb.species)
            }


            if (metabolism.geneset == TRUE  & is.null(gmtFile)) {
              signatures_KEGG_metab <- system.file("data", "KEGG_metabolism_nc.gmt", package = "sc.utils")
              signatures_REACTOME_metab <- system.file("data", "REACTOME_metabolism.gmt", package = "sc.utils")
              if (metabolism.type == "KEGG")  {gmtFile<-signatures_KEGG_metab; cat("Your choice is: KEGG\n")}
              if (metabolism.type == "REACTOME")  {gmtFile<-signatures_REACTOME_metab; cat("Your choice is: REACTOME\n")}
              cat(metabolism.type)
            }

            # check gmt file
            if (is.null(gmtFile)) {
              stop("gmt file should be provided! ")
            }

            #----------------calculate score--------------------------------

            # for ssGSEA or GSVA, we recommend to use the whole gene expression matrix rather than matrix only containing highly variable genes
            #ssGSEA
            if (method == "ssgsea") {
              library(GSVA)
              library(GSEABase)
              geneSets <- GSEABase::getGmt(gmtFile) #signature read
              gsva_es <- GSVA::gsva(as.matrix(expr), geneSets, method=c("ssgsea"), kcdf=c("Poisson"), parallel.sz=ncores) #
              signature_exp<-data.frame(gsva_es)

            }

            #GSVA
            if (method == "gsva") {
              library(GSVA)
              library(GSEABase)
              geneSets <- GSEABase::getGmt(gmtFile) #signature read
              gsva_es <- GSVA::gsva(as.matrix(expr), geneSets, method=c("gsva"), kcdf=c("Poisson"), parallel.sz=ncores) #
              signature_exp<-data.frame(gsva_es)
            }

            object@assays$pathway_score$score<-signature_exp

            return(object)
}

##

#' @title transform pathway score per cell to pathway score per cell type
#' @param object Seurat object from sc.pathway.score.seurat
#' @param cluster the feature to aggregate
#' @param avg_method aggregate method, mean or median
#' @return a data frame
#' @export
aggregate_pathway_score <- function(object,cluster = "cell_type",avg_method = c("mean","median")){
  #object = pbmc
  metainfo <- object@meta.data
  signature_exp <- object@assays$pathway_score$score
  signature_exp <- signature_exp %>% t() %>% as.data.frame()
  signature_exp <- cbind(metainfo[,cluster], signature_exp)
  names(signature_exp)[1] = "Cluster"
  signature_exp$Cluster = as.character(signature_exp$Cluster)

  tmp = aggregate(signature_exp[, 2:ncol(signature_exp)], list(signature_exp$Cluster), mean)

  return(tmp)
}

#' @title dotplot to visualize pathway score per cell type
#' @param df data frame from function `aggregate_pathway_score`
#' @import ggplot2 tidyr dplyr
#' @export
#' @return `ggplot2` object
dotplot_pathway_score <- function(df, topn = 5, normalization = c("zeroone","z-score")){
  tmp <- df %>% pivot_longer(!Group.1 ,names_to = "pathway", values_to = "score")

  if (normalization == "zeroone") {
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    cat("Perfrom 0 1 scaling\n")
    tmp$score = as.numeric(range01(tmp$score))
  } else if (normalization == "z-score"){
    cat("Perfrom z-score scaling\n")
    tmp = tmp %>% mutate(score = (score - mean(score))/sd(score))
  }

  tmp = tmp %>%
    dplyr::group_by_(~ pathway) %>%
    dplyr::arrange_(~ desc(score)) %>%
    dplyr::top_n(n = topn)

  ggplot(tmp, aes(x= Group.1, y= pathway, size= score, color = score)) +
    geom_point() +
    scale_color_continuous(low="red", high="blue", name = "score",
                           guide=guide_colorbar(reverse=TRUE)) +
    xlab(NULL) +
    ylab(NULL) +
    scale_size(range=c(3, 8))
}
