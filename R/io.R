import.sce.from.count.matrix <- function(counts, gene.names, sample_annotations = NULL, prefilter = FALSE, min.cell.frac.per.gene = 0.001, min.genes.per.cell = 300) {
	if(!is.sparseMatrix(counts)) {
		counts = as(counts, "sparseMatrix")
	}	
	rownames(counts) = gene.names    
    
    if(is.null(sample_annotations)) {
		sce <- SingleCellExperiment(assays = list(counts = counts))
    } else {
        sce <- SingleCellExperiment(assays = list(counts = counts), colData = sample_annotations)
	}
	
    if (prefilter) {
		min.cells.per.gene = round(ncol(sce) * min.cell.frac.per.gene)
        cell.counts = Matrix::rowSums(SummarizedExperiment::assays(sce)$counts > 0)
        sce = sce[cell.counts > min.cells.per.gene, ]
        
        feature.counts = Matrix::colSums(SummarizedExperiment::assays(sce)$counts > 0)
        sce = sce[, feature.counts > min.genes.per.cell]
    }
    
    return(sce)
}

import.sce.from.table <- function(fname, sep = "\t", header = TRUE, prefilter = FALSE, min.cell.frac.per.gene = 0.001, min.genes.per.cell = 300) {
    require(Matrix)
    require(SingleCellExperiment)
    
    counts = read.table(fname, header = TRUE, sep = sep, as.is = TRUE)
    
    if (!is.numeric(counts[1, 1])) {
        row.names = counts[, 1]
        counts = counts[, -1]
        rownames(counts) = row.names
    }
    
    counts = as(as.matrix(counts), "sparseMatrix")
    
    sce <- SingleCellExperiment(assays = list(counts = counts))
    
    if (prefilter) {
		min.cells.per.gene = round(ncol(sce) * min.cell.frac.per.gene)
        cell.counts = Matrix::rowSums(SummarizedExperiment::assays(sce)$counts > 0)
        sce = sce[cell.counts > min.cells.per.gene, ]
        
        feature.counts = Matrix::colSums(SummarizedExperiment::assays(sce)$counts > 0)
        sce = sce[, feature.counts > min.genes.per.cell]
    }
    
    return(sce)
}

import.sce.from.10X <- function(input_path, mtx_file = "matrix.mtx.gz", feature_annotations = "features.tsv.gz", sample_annotations = "barcodes.tsv.gz", 
    sep = "\t", header = FALSE, prefilter = FALSE, min.cell.frac.per.gene = 0.001, min.genes.per.cell = 300) {
    require(Matrix)
    require(scran)
    
    counts = readMM(paste(input_path, mtx_file, sep = "/"))
    feature_table = read.table(paste(input_path, feature_annotations, sep = "/"), header = header, sep = sep, as.is = TRUE)
    sample_annotations = read.table(paste(input_path, sample_annotations, sep = "/"), header = header, sep = sep, as.is = TRUE)

    # Handle header
    if(nrow(feature_table) == nrow(counts) + 1) {
		colnames(feature_table) = feature_table[1, ]
		feature_table = feature_table[-1, ]
	}
    if(nrow(sample_annotations) == ncol(counts) + 1) {
		colnames(sample_annotations) = sample_annotations[1, ]
		sample_annotations = sample_annotations[-1, ]
	}

    
    if(ncol(feature_table) > 2) {
		IDX = split(1:nrow(feature_table), feature_table$V3)
		expression.counts = counts[IDX$`Gene Expression`, ]
		gene.table = feature_table[IDX$`Gene Expression`, 1:2]    
		IDX = IDX[!grepl("Gene Expression", names(IDX))]
	} else {
		expression.counts = counts
		gene.table = feature_table
		IDX = list()
	}
	
    colnames(gene.table) = c("ENSEMBL", "SYMBOL")
    rownames(expression.counts) = gene.table$SYMBOL    
    colnames(expression.counts) = sample_annotations[, 1]
    
    
    if (ncol(sample_annotations) > 1) {
        sce <- SingleCellExperiment(assays = list(counts = expression.counts), colData = sample_annotations, rowData = gene.table[, 
            -3])
    } else {
        sce <- SingleCellExperiment(assays = list(counts = expression.counts), rowData = gene.table)
    }
    
    # Load additional barcoded features
    for(feature.name in names(IDX)) {
		feature.counts = counts[IDX[[feature.name]], ]
		
		row.annotations = feature_table[IDX[[feature.name]], ]
		
		rownames(feature.counts) = row.annotations[, 1]
		colnames(feature.counts) = sample_annotations[, 1]
		
		
		reducedDims(sce)[[feature.name]] = Matrix::t(feature.counts)
	}
	
    
    if (prefilter) {
		min.cells.per.gene = round(ncol(sce) * min.cell.frac.per.gene)
        cell.counts = Matrix::rowSums(SummarizedExperiment::assays(sce)$counts > 0)
        sce = sce[cell.counts > min.cells.per.gene, ]
        
        feature.counts = Matrix::colSums(SummarizedExperiment::assays(sce)$counts > 0)
        sce = sce[, feature.counts > min.genes.per.cell]
    }
    
    return(sce)
}


# https://satijalab.org/seurat/v3.0/conversion_vignette.html
import.sce.from.Seurat <- function(Seurat.obj) {
    library(Seurat)
    
    sce <- Seurat::as.SingleCellExperiment(Seurat.obj)
    
    return(sce)
}

import.sce.from.AnnData <- function(AnnData.fname, export_all = FALSE) {
    sce <- Seurat::as.SingleCellExperiment(Seurat::ReadH5AD(file = AnnData.fname))
    
    return(sce)
}

import.sce.from.CDS <- function(monocle_cds) {
    
    counts = exprs(monocle_cds)
    gene_annotations = fData(monocle_cds)
    sample_annotations = pData(monocle_cds)
    
    sce <- SingleCellExperiment(assays = list(counts = counts), colData = sample_annotations, rowData = gene_annotations)
    
    return(sce)
}


convert.sce.rownames <- function(sce, from = "ENSEMBL", to = "SYMBOL", species = "human") {
    if (species == "human") {
        library(org.Hs.eg.db)
        suppressWarnings(ids <- mapIds(org.Hs.eg.db, keys = row.names(sce), keytype = from, column = to, multiVals = "first"))
        ids[is.na(ids)] = ""
        
        rownames(sce) = ids
    } else if (species == "mouse") {
        library(org.Mm.eg.db)
        suppressWarnings(ids <- mapIds(org.Mm.eg.db, keys = row.names(sce), keytype = from, column = to, multiVals = "first"))
        
        ids[is.na(ids)] = ""
        rownames(sce) = ids
    }
    
    return(sce)
}
