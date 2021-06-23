
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sc.utils

<!-- badges: start -->
<!-- badges: end -->

The goal of sc.utils is to …

## Installation

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Byronxy/sc.utils")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(sc.utils)
## basic example code
```

## Support Gene Signatures

### Cell Marker Database

| FileName            | Description                                                                         | Usage                                                                                          | ID Format |
|---------------------|-------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------|-----------|
| All cell markers    | All cell markers of different cell types from different tissues in human and mouse. | data(“cell\_markers\_all\_entrez”,package = “enrichCellMarkers”, envir = environment())        | EntrezID  |
| All cell markers    | All cell markers of different cell types from different tissues in human and mouse. | data(“cell\_markers\_all\_symbal”,package = “enrichCellMarkers”, envir = environment())        | SymbolID  |
| Human cell markers  | Cell markers of different cell types from different tissues in human.               | data(“cell\_markers\_human\_entrez”,package = “enrichCellMarkers”, envir = environment())      | EntrezID  |
| Human cell markers  | Cell markers of different cell types from different tissues in human.               | data(“cell\_markers\_human\_symbal”,package = “enrichCellMarkers”, envir = environment())      | SymbolID  |
| Mouse cell markers  | Cell markers of different cell types from different tissues in mouse.               | data(“cell\_markers\_mouse\_entrez”,package = “enrichCellMarkers”, envir = environment())      | EntrezID  |
| Mouse cell markers  | Cell markers of different cell types from different tissues in mouse.               | data(“cell\_markers\_mouse\_symbal”,package = “enrichCellMarkers”, envir = environment())      | SymbolID  |
| Single cell markers | Cell markers derived from single-cell sequencing researches in human and mouse.     | data(“cell\_markers\_singlecell\_entrez”,package = “enrichCellMarkers”, envir = environment()) | EntrezID  |
| Single cell markers | Cell markers derived from single-cell sequencing researches in human and mouse.     | data(“cell\_markers\_singlecell\_symbal”,package = “enrichCellMarkers”, envir = environment()) | SymbolID  |

``` r
#entrez
data("cell_markers_all_entrez",package = "enrichCellMarkers", envir = environment())

data("cell_markers_human_entrez",package = "enrichCellMarkers", envir = environment())

data("cell_markers_mouse_entrez",package = "enrichCellMarkers", envir = environment())

data("cell_markers_singlecell_entrez",package = "enrichCellMarkers", envir = environment())

#symbol
data("cell_markers_all_symbal",package = "enrichCellMarkers", envir = environment())

data("cell_markers_human_symbal",package = "enrichCellMarkers", envir = environment())

data("cell_markers_mouse_symbal",package = "enrichCellMarkers", envir = environment())

data("cell_markers_singlecell_symbal",package = "enrichCellMarkers", envir = environment())
```

### In-house Pipeline Signature

#### gmt file

-   gbm\_single\_cell\_geneset.gmt

-   SHH\_MB\_hg\_signature.gmt

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
gmtFile <- paste(file.path(system.file('examples', package='enrichCellMarkers')), "SHH_MB_hg_signature.gmt", sep="/")
geneSets <- getGmt(gmtFile)
```

#### load data

``` r
#cell cycle
data("cellcycle_hg",package = "enrichCellMarkers", envir = environment())
data("cellcycle_mm",package = "enrichCellMarkers", envir = environment())

#gbm single cell signature
data("gbm_single_cell_geneset",package = "enrichCellMarkers", envir = environment())
data("gbm_single_cell_geneset_list",package = "enrichCellMarkers", envir = environment())

#shh-mg signature
data("shh_mb_hg_single_cell_geneset_list",package = "enrichCellMarkers", envir = environment())
data("shh_mb_mm_single_cell_geneset_list",package = "enrichCellMarkers", envir = environment())
```

### Useful Function

#### Data Visualization

#### Data Transformation

#### Statistics
