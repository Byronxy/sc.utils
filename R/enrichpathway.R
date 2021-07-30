#' @title use clusterprofiler to do enrichment analysis per cell type
#' @import clusterProfiler
#' @import Seurat
#' @param object Seurat object
#' @param fun fun param in clusterProfiler::compareCluster

dotplot_encich_clusterprofiler <- function(object,fun,...){
  cat("perform DEG analysis\n")
  markers <- Seurat::FindAllMarkers(object = object, only.pos = TRUE, min.pct = 0.25,
                                 thresh.use = 0.25,...)
  cat("perform id transfrom\n")
  ids=clusterProfiler::bitr(markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
  markers=merge(markers,ids,by.x='gene',by.y='SYMBOL')
  gcSample=split(markers$ENTREZID, markers$cluster)
  cat("perform enrichment\n")
  xx <- clusterProfiler::compareCluster(geneClusters = gcSample, fun=fun,...)
  print(xx)
  cat("visualization\n")
  p=dotplot(xx)

  p = p + theme(axis.text.x = element_text(angle = 45,
                                      vjust = 0.5, hjust=0.5))
  #p
  return(p)
}


