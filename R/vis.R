#' Function to plot postive cells or negative cells with two colors
#'
#' @param object Seurat object
#' @param gene gene name
#' @param cutoff cutoff
#' @param color_value color value  c("grey","red")
#' @param legend whether use legend, default: FALSE
#' @param reduction what is the reduction, default: umap
#'
#' @return  a ggplot2 object
#'
#'
#' @export
#'
#'
FeaturePlot_gene_pos <- function(object,
                                 gene,
                                 cutoff = 0,
                                 color_value = c("grey","red"),
                                 legend = FALSE,
                                 reduction = "umap"){
  stopifnot(length(gene) == 1)
  p = FeaturePlot(object,
                  features = gene,
                  order = T,
                  reduction = reduction)
  data = p$data
  gene_exp = as.numeric(data[,4])
  data = data %>%
    # rownames_to_column(var = "cell") %>%
    mutate(gene_pos = ifelse(gene_exp > cutoff,"Pos","Neg"))
  names(data)[1] = "x"
  names(data)[2] = "y"
  p1 <- data %>%
    ggplot(aes(x,y, color = gene_pos)) +
    geom_point(size = .2,alpha = .8)+
    theme_cowplot()+
    scale_color_manual(values = color_value) +
    # scale_color_brewer(palette = "RdYlBu")+
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    labs(x = "", y = "") +
    ggtitle(gene) +
    guides(colour = guide_legend(override.aes = list(size=3)))

  if (legend == FALSE) {
    p1 =  p1 + theme(legend.position = "none")
  }
  return(p1)

}
