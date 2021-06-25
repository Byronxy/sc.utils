#' CellMarker Enrichment
#'
#' Given a list of markers, indentifies the putative cell type using the comphensive catlog creatd in
#'
#' @param gene.list character vector with a list of \bold{Gene Symbols}
#' @param species Can either be \emph{human} or \emph{mouse} for the moment.
#'
#' @return  A \emph{list} two elememts. First element (enrichments) is a data.frame structured as follows:
#' \itemize{
#' \item {\bold{p.value}} {P-value of the enrichment, the lower the better}
#' \item {\bold{overlap}} {How many common genes between your costom list and the markers of the given cell type}
#' \item {\bold{signature}} {Total markers of the given cell type }
#' \item {\bold{genes}} {Common genes}
#' }
#' The second element is a list with the actual gene names for each enrichment.
#'
#' @examples
#' results=CMenrich(gene.list=c('NEUROD1','ARX','CD79A'),species='human') .# some gene names that I came up with ....
#'
#' @export
#'
#'
CMenrich = function (gene.list,species)
{

  if (species=='human')
    data('human.CMenrich')
  else
    data('mouse.CMenrich')

  #print(gene.list)

  gene.list=intersect(gene.list,all.markers)

  #print(gene.list)

  pop.size=length(all.markers)
  samp.size=length(gene.list)

  p.val=rep(0,length(markers))
  samp.hits=rep(0,length(markers))
  pop.hits=rep(0,length(markers))
  genes=list()

  for (k in 1:length(markers))
  {
    samp.hits[k]=length(intersect(gene.list,markers[[k]]))
    pop.hits[k]=length(markers[[k]])
    p.val[k]=phyper(samp.hits[k], pop.hits[k], pop.size-pop.hits[k], samp.size,lower.tail = FALSE)+dhyper(samp.hits[k], pop.hits[k], pop.size-pop.hits[k], samp.size)
    genes[[k]]=intersect(gene.list,markers[[k]])
  }

  ix=order(p.val)

  enrich=data.frame(Cell.type=cell.type, p.value=p.val,overlap=samp.hits,signature=pop.hits)
  enrich=enrich[ix,]

  result=list(enrichments=enrich,genes=genes[ix])

}

#' Convert list to gmt file
#' @param geneSet a list object
#' @param gmt_file output file name .gmt
#'
#' @export
#'
write.gmt <- function(geneSet= NULL,gmt_file= NULL){
  if(class(geneSet) != "list"){stop("geneSet need to be a list!")}
  sink( gmt_file )
  for (i in 1:length(geneSet)){
    cat(names(geneSet)[i])
    cat('\tNA\t')
    cat(paste(geneSet[[i]],collapse = '\t'))
    cat('\n')
  }
  sink()
}


#' Covert between mouse and human symbols
#'
#' @param genesig vector, gene signature
#' @param to 'human' or 'mouse'
#' @export
convert_symbol <- function(genesig = NULL, to = "mouse"){
  if(class(genesig) != "character"){stop("genesig need to be a character!")}
  if (!requireNamespace("biomaRt")) {
    stop("Please install 'biomaRt' package firstly!")
  }
  require("biomaRt")
  #biomRt_path <- paste(file.path(system.file('examples', package='enrichCellMarkers')), "biomRt_input.Rdata", sep="/")
  data(mouse, package = "sc.utils")
  data(human, package = "sc.utils")
  if (to == "mouse") {
    genesig_res <- getLDS(attributes = c("hgnc_symbol"),filters = "hgnc_symbol",
                          values = genesig, mart = human,
                          attributesL = c("mgi_symbol"), martL = mouse,
                          uniqueRows = T)
    # genesig_res <- genesig_res$MGI.symbol
    message("convert from human to mouce")
  } else {
    genesig_res <- getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",
                          values = genesig, mart = mouse,
                          attributesL = c("hgnc_symbol"), martL = human,
                          uniqueRows = T)
    # genesig_res <- genesig_res$HGNC.symbol
    message("convert from mouce to human")
  }
  return(genesig_res)
}
