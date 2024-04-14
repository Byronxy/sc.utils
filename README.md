
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sc.utils

<!-- badges: start -->
<!-- badges: end -->

The goal of sc.utils is to integrate useful function for single cell
sequencing data analysis

## Installation

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Byronxy/sc.utils")
```

## Example

This is a basic example which shows you how to solve a common problem:

    library(sc.utils)
    ## basic example code

## Support Gene Signatures

### Cell Marker Database

| FileName            | Description                                                                         | Usage                                                                              | ID Format |
|---------------------|-------------------------------------------------------------------------------------|------------------------------------------------------------------------------------|-----------|
| All cell markers    | All cell markers of different cell types from different tissues in human and mouse. | data(“cell_markers_all_entrez”,package = “sc.utils”, envir = environment())        | EntrezID  |
| All cell markers    | All cell markers of different cell types from different tissues in human and mouse. | data(“cell_markers_all_symbal”,package = “sc.utils”, envir = environment())        | SymbolID  |
| Human cell markers  | Cell markers of different cell types from different tissues in human.               | data(“cell_markers_human_entrez”,package = “sc.utils”, envir = environment())      | EntrezID  |
| Human cell markers  | Cell markers of different cell types from different tissues in human.               | data(“cell_markers_human_symbal”,package = “sc.utils”, envir = environment())      | SymbolID  |
| Mouse cell markers  | Cell markers of different cell types from different tissues in mouse.               | data(“cell_markers_mouse_entrez”,package = “sc.utils”, envir = environment())      | EntrezID  |
| Mouse cell markers  | Cell markers of different cell types from different tissues in mouse.               | data(“cell_markers_mouse_symbal”,package = “sc.utils”, envir = environment())      | SymbolID  |
| Single cell markers | Cell markers derived from single-cell sequencing researches in human and mouse.     | data(“cell_markers_singlecell_entrez”,package = “sc.utils”, envir = environment()) | EntrezID  |
| Single cell markers | Cell markers derived from single-cell sequencing researches in human and mouse.     | data(“cell_markers_singlecell_symbal”,package = “sc.utils”, envir = environment()) | SymbolID  |

``` r
#entrez
data("cell_markers_all_entrez",package = "sc.utils", envir = environment())

data("cell_markers_human_entrez",package = "sc.utils", envir = environment())

data("cell_markers_mouse_entrez",package = "sc.utils", envir = environment())

data("cell_markers_singlecell_entrez",package = "sc.utils", envir = environment())

#symbol
data("cell_markers_all_symbal",package = "sc.utils", envir = environment())

data("cell_markers_human_symbal",package = "sc.utils", envir = environment())

data("cell_markers_mouse_symbal",package = "sc.utils", envir = environment())

data("cell_markers_singlecell_symbal",package = "sc.utils", envir = environment())
```

### In-house Pipeline Signature

#### gmt file

- gbm_single_cell_geneset.gmt

- SHH_MB_hg_signature.gmt

<!-- end list -->

``` r
library(GSEABase)
#> Loading required package: annotate
#> Loading required package: XML
#> Loading required package: graph
#> 
#> Attaching package: 'graph'
#> The following object is masked from 'package:XML':
#> 
#>     addNode
gmtFile <- paste(file.path(system.file('examples', package='sc.utils')), "SHH_MB_hg_signature.gmt", sep="/")
geneSets <- getGmt(gmtFile)
```

#### load data

``` r
#cell cycle
data("cellcycle_hg",package = "sc.utils", envir = environment())
data("cellcycle_mm",package = "sc.utils", envir = environment())

#gbm single cell signature
data("gbm_single_cell_geneset",package = "sc.utils", envir = environment())
data("gbm_single_cell_geneset_list",package = "sc.utils", envir = environment())

#shh-mg signature
data("shh_mb_hg_single_cell_geneset_list",package = "sc.utils", envir = environment())
data("shh_mb_mm_single_cell_geneset_list",package = "sc.utils", envir = environment())
```

To avoid unexpected noise and expression artefacts by dissociation, a
total of 1,514 genes associated with mitochondria (50 genes), heat-shock
protein (178 genes), ribosome (1,253 genes) and dissociation (33 genes)
were excluded.

``` r
data("gs_MT_HSP_RB_DS",package = "sc.utils", envir = environment())
```

### Useful Function

#### Data Visualization

``` r
FeaturePlot_gene_pos()
```

#### Data Transformation

#### Statistics
