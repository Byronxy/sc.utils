
#' @title write mtx for first step in pyscenic
#' @export
#' @import SeuratObject
#' @param object `Seurat` object
#' @param path write path
#' @param slot slot for `Seurat` object
#' @param assay assay for `Seurat` object

write_pyscenic_mtx <- function(object, path =  "/export/bioinfo-team/home/xiongy/pyscenic/mtx.tsv", slot = "data", assay = "RNA"){
  mtx <- SeuratObject::GetAssayData(object = object, slot = slot, assay = assay)
  cat("Transforming\n")
  mtx <- t(as.matrix(mtx))
  cat("Writing\n")
  write.table(mtx, file = path,  sep = "\t", col.names = NA)
  cat("Finish\n")
}

#' @title get pseudocell for `Seurat` object
#' @export
#' @import SeuratObject
#' @param object `Seurat` object
#' @param pseudocell_size pseudocell size (default: 5)
#' @param group_id character, the column to aggregate in meta data
#' @param seed set.seed
get_pseudocell <- function(object = NULL,
                           slot = "data",
                           assay = "RNA",
                           group_id = "cell_type",
                           seed = 123,
                           pseudocell_size = 5){
  # refer to: https://github.com/JiaqiLiZju/CelltypeEvolution/blob/master/Pseudocell_analysis/Pseudocell_analysis_pipeline.r
  '%!in%' <- function(x,y)!('%in%'(x,y))

  mtx <- SeuratObject::GetAssayData(object = object, slot = slot, assay = assay)
  meta_data = object@meta.data
  # head(meta_data)
  if (group_id %!in% names(meta_data)) {
    stop("the group column should be in meta data")
  }
  meta_data$cell_ids = rownames(meta_data)
  meta_data <- meta_data[,c("cell_ids", group_id)]
  meta_data[,group_id] = as.character(meta_data[,group_id])
  names(meta_data) = c("cell_ids", "group_id")

  # subset into pseudocell
  set.seed(seed)

  pseudocell.size = pseudocell_size
  new_ids_list = list()
  for (i in 1:length(unique(meta_data$group_id))) {
    cluster_id = unique(meta_data$group_id)[i]
    cluster_cells <- rownames(meta_data[meta_data$group_id == cluster_id,])
    cluster_size <- length(cluster_cells)
    pseudo_ids <- floor((seq_along(cluster_cells)-1)/pseudocell.size)
    pseudo_ids <- paste0(cluster_id, "_Cell", pseudo_ids)
    names(pseudo_ids) <- sample(cluster_cells)
    new_ids_list[[i]] <- pseudo_ids
  }


  #assign new cell ids
  new_ids <- unlist(new_ids_list)
  new_ids <- as.data.frame(new_ids)
  new_ids_length <- table(new_ids)
  new_ids$cell_ids = rownames(meta_data)

  new_colnames <- rownames(new_ids)  ###add
  all.data <- mtx[,as.character(new_colnames)] ###add
  all.data <- t(all.data)###add

  #aggregate
  new.data<-aggregate(list(all.data[,1:length(all.data[1,])]),
                      list(name=new_ids[,1]),FUN=mean)
  rownames(new.data)<-new.data$name
  new.data<-new.data[,-1]

  new_ids$cell_type = gsub("[_]Cell.*$","",new_ids$new_ids)

  mtx_out = new.data
  meta_data_full_out = new_ids
  meta_data_out = meta_data_full_out %>% dplyr::select(new_ids, cell_type) %>% dplyr::distinct(.)
  rownames(meta_data_out) = meta_data_out$new_ids
  meta_data_out = meta_data_out[rownames(mtx_out),]

  out <- list(mtx_out = mtx_out,
              meta_data_full_out = meta_data_full_out,
              meta_data_out = meta_data_out)

  return(out)
}




