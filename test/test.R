
library(Seurat)
library(sc.utils)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(tidyverse)
pbmc <- readRDS("~/project/tutorials/pbmc3k_final.rds")
pbmc$cell_type = pbmc@active.ident
pbmc <- sc.pathway.score(object = pbmc, method = "ssgsea",ncores = 4, metabolism.geneset = TRUE, metabolism.type = "KEGG")

pbmc@assays$pathway_score$score[1:5,1:5]


# vis pathway
df <- aggregate_pathway_score(pbmc)
dotplot_pathway_score(df, normalization = "zscore")


ggplot(df, aes_string(x=x, y="Description", size=size, color=colorBy)) +
  geom_point() +
  scale_color_continuous(low="red", high="blue", name = color,
                         guide=guide_colorbar(reverse=TRUE)) +
  scale_y_discrete(labels = label_func) +
  ylab(NULL) + ggtitle(title) + theme_dose(font.size) +
  scale_size(range=c(3, 8))

pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25,
                              thresh.use = 0.25)




ids=bitr(pbmc.markers$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
pbmc.markers=merge(pbmc.markers,ids,by.x='gene',by.y='SYMBOL')
gcSample=split(pbmc.markers$ENTREZID, pbmc.markers$cluster)
gcSample # entrez id , compareCluster
xx <- compareCluster(geneClusters = gcSample, fun="enricher",
                     organism="hsa", pvalueCutoff=0.05)

p=dotplot(xx)
p+ theme(axis.text.x = element_text(angle = 45,
                                    vjust = 0.5, hjust=0.5))
p

library(msigdbr)
gmt <- msigdbr(species = "Homo sapiens")
gmt2 <- gmt%>%
  dplyr::select(gs_name, entrez_gene)
gmts <- split(gmt2, gmt$gs_cat)

## TERM2GENE need a data frame
xx <- compareCluster(geneClusters = gcSample, fun="enricher",
                     TERM2GENE = gmt2)

dotplot_encich_clusterprofiler(pbmc,fun="enricher",TERM2GENE = gmt2)


data("pbmc_small")
pbmc_small
pbmc <- BuildClusterTree(object = pbmc, graph = pbmc@graphs)
Tool(object = pbmc, slot = 'BuildClusterTree')

plot(pbmc@tools$BuildClusterTree)

tree = pbmc@tools$BuildClusterTree
X <- c("red", "orange", "yellow", "green", "blue", "purple")
plot(tree,
     edge.color = sample(X, length(tree$edge)/2, replace = TRUE),
     edge.width = sample(1:10, length(tree$edge)/2, replace = TRUE))

library(ggtree)
ggtree::ggtree(tr = tree) + geom_tiplab()

out <- get_pseudocell(object = pbmc)
