#' Regulons for Viper analysis
#' 
#' @format:
#'   A list of regulons as input to the Viper method
#' 
#' @source \url{https://www.bioconductor.org/packages/release/bioc/html/viper.html}
"ViperRegulons"


#' Pre-processed genesets from the gProfiler database (Human)
#' 
#' @format:
#'   SYMBOL: A list of gene x geneset indicator matrices in the sparseMatrix format (rownames are gene symbols)
#'   ENSEMBL: A list of gene x geneset indicator matrices in the sparseMatrix format (rownames are ENSEMBL gene IDs)
#' 
#' @source \url{https://biit.cs.ut.ee/gprofiler}
"gProfilerDB_human"


#' Pre-processed genesets from the gProfiler database (Mouse)
#' 
#' @format:
#'   SYMBOL: A list of gene x geneset indicator matrices in the sparseMatrix format (rownames are gene symbols)
#'   ENSEMBL: A list of gene x geneset indicator matrices in the sparseMatrix format (rownames are ENSEMBL gene IDs)
#' 
#' @source \url{https://biit.cs.ut.ee/gprofiler}
"gProfilerDB_mouse"


#' Pre-processed genesets from the nanoString database (Human)
#' 
#' @format:
#'   A list of gene x geneset indicator matrices in the sparseMatrix format (rownames are gene symbols)
#' 
#' @source \url{https://www.nanostring.com/}
"nanoStringDB_human"


#' Pre-processed genesets from the nanoString database (Mouse)
#' 
#' @format:
#'   A list of gene x geneset indicator matrices in the sparseMatrix format (rownames are gene symbols)
#' 
#' @source \url{https://www.nanostring.com/}
"nanoStringDB_mouse"


#' Pre-processed cell type markers from the CellMarker database (Human)
#' 
#' @format:
#'   A gene x cell-type indicator matrix in the sparseMatrix format (rownames are gene symbols)
#' 
#' @source \url{http://biocc.hrbmu.edu.cn/CellMarker/}
"CellMarkerDB_human"


#' Pre-processed cell type markers from the CellMarker database (Mouse)
#' 
#' @format:
#'   A gene x cell-type indicator matrix in the sparseMatrix format (rownames are gene symbols)
#' 
#' @source \url{http://biocc.hrbmu.edu.cn/CellMarker/}
"CellMarkerDB_mouse"



#' Curated cell type markers for hand-selected cell types (Human)
#' 
#' @format:
#'   A list of gene symbols for different tissues/cell-types
#' 
#' @source \url{http://compbio.mit.edu/ACTIONet/}
"curatedMarkers_human"



#' Reference profile for cell sorted PBMC
#' 
#' @format:
#'   SingleCellExperiment object with rowData, colData, and metadata inside
#' 
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/30726743}
"Monaco2019"


#' Enrichr functional pathways
#' 
#' @format:
#'   A list of genesets
#' 
#' @source \url{https://amp.pharm.mssm.edu/Enrichr/}
"EnrichrDB"


#' Collection of Disease-associated genesets
#' 
#' @format:
#'   A list of genesets
#' 
#' @source \url{https://amp.pharm.mssm.edu/Enrichr/}
"DiseaseDB"


#' Collection of TF-TG assocations
#' 
#' @format:
#'   A list of TF-TG associations
#' 
#' @source \url{https://amp.pharm.mssm.edu/chea3/}
"ChEA3plusDB"
