library(readxl)
library(tidyverse)
library(usethis)
library(devtools)
#' gs_MT_HSP_RB_DS
#'
#' Ref: Liver tumour immune microenvironment subtypes and neutrophil heterogeneity (Supplementary Table 1)
#'
#' To avoid unexpected noise and expression artefacts by dissociation, a total of 1,514
#' genes associated with mitochondria (50 genes), heat-shock protein
#' (178 genes), ribosome (1,253 genes) and dissociation (33 genes) were
#' excluded
# data <- read_excel("C:/Users/39047/Nutstore/1/Nutstore/生物信息/单细胞/artefact_1514genes.xlsx")
# gs <- data.frame(gene = c(data$Mitochondria,data$`Heat shock protein`,data$Ribosome,data$Dissociation),
#                  gs_name = c(rep("Mitochondira",length(data$Mitochondria)),
#                               rep("Heat shock protein",length(data$`Heat shock protein`)),
#                               rep("Ribosome",length(data$Ribosome)),
#                               rep("Dissociation",length(data$Dissociation)))) %>% filter(!is.na(gene))
#
# gs_MT_HSP_RB_DS <- gs
#
#
# use_data(gs_MT_HSP_RB_DS, internal = FALSE)
# data("gs_MT_HSP_RB_DS")
