import.sce.from.count.matrix <- function(counts, gene.names = NULL, prefilter = FALSE, min.cells.per.gene = 10, min.genes.per.cell = 300) {
    if (!is.null(gene.names)) {
        rownames(counts) = gene.names
    }
    
    counts = as(as.matrix(counts), "sparseMatrix")
    
    sce <- SingleCellExperiment(assays = list(counts = counts))
    
    if (prefilter) {
        cell.counts = Matrix::rowSums(sce@assays[["counts"]] > 0)
        sce = sce[cell.counts > min.cells.per.gene, ]
        
        feature.counts = Matrix::colSums(sce@assays[["counts"]] > 0)
        sce = sce[, feature.counts > min.genes.per.cell]
    }
    
    return(sce)
}

import.sce.from.table <- function(fname, sep = "\t", header = TRUE, prefilter = FALSE, min.cells.per.gene = 10, min.genes.per.cell = 300) {
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
        cell.counts = Matrix::rowSums(sce@assays[["counts"]] > 0)
        sce = sce[cell.counts > min.cells.per.gene, ]
        
        feature.counts = Matrix::colSums(sce@assays[["counts"]] > 0)
        sce = sce[, feature.counts > min.genes.per.cell]
    }
    
    return(sce)
}

import.sce.from.10X <- function(input_path, mtx_file = "matrix.mtx.gz", gene_annotations = "features.tsv.gz", sample_annotations = "barcodes.tsv.gz", 
    sep = "\t", header = FALSE, prefilter = FALSE, min.cells.per.gene = 10, min.genes.per.cell = 300) {
    require(Matrix)
    require(scran)
    
    counts = readMM(paste(input_path, mtx_file, sep = "/"))
    gene.table = read.table(paste(input_path, gene_annotations, sep = "/"), header = header, sep = sep, as.is = TRUE)
    sample_annotations = read.table(paste(input_path, sample_annotations, sep = "/"), header = header, sep = sep, as.is = TRUE)
    
    rownames(counts) = gene.table[, 1]
    colnames(counts) = sample_annotations[, 1]
    
    sce <- SingleCellExperiment(assays = list(counts = counts), colData = sample_annotations, rowData = gene.table)
    
    if (prefilter) {
        cell.counts = Matrix::rowSums(sce@assays[["counts"]] > 0)
        sce = sce[cell.counts > min.cells.per.gene, ]
        
        feature.counts = Matrix::colSums(sce@assays[["counts"]] > 0)
        sce = sce[, feature.counts > min.genes.per.cell]
    }
    
    return(sce)
}

import.sce.from.10X.TotalSeq <- function(input_path, mtx_file = "matrix.mtx.gz", feature_annotations = "features.tsv.gz", sample_annotations = "barcodes.tsv.gz", 
    sep = "\t", header = FALSE, prefilter = FALSE, min.cells.per.gene = 10, min.genes.per.cell = 300) {
    require(Matrix)
    require(scran)
    
    counts = readMM(paste(input_path, mtx_file, sep = "/"))
    feature.table = read.table(paste(input_path, feature_annotations, sep = "/"), header = header, sep = sep, as.is = TRUE)
    sample_annotations = read.table(paste(input_path, sample_annotations, sep = "/"), header = header, sep = sep, as.is = TRUE)
    
    
    IDX = split(1:nrow(feature.table), feature.table$V3)
    
    expression.counts = counts[IDX$`Gene Expression`, ]
    
    gene.table = feature.table[IDX$`Gene Expression`, 1:2]
    colnames(gene.table) = c("ENSEMBL", "SYMBOL")
    rownames(expression.counts) = gene.table$SYMBOL
    colnames(expression.counts) = sample_annotations$V1
    
    if (ncol(sample_annotations) > 1) {
        sce <- SingleCellExperiment(assays = list(counts = expression.counts), colData = sample_annotations, rowData = gene.table[, 
            -3])
    } else {
        sce <- SingleCellExperiment(assays = list(counts = expression.counts), rowData = gene.table)
    }
    
    ABC.counts = counts[IDX$`Antibody Capture`, ]
    ABC.table = feature.table[IDX$`Antibody Capture`, ]
    rownames(ABC.counts) = ABC.table$V1
    colnames(ABC.counts) = sample_annotations$V1
    
    
    sce@reducedDims[["ABC"]] = Matrix::t(ABC.counts)
    
    
    if (prefilter) {
        cell.counts = Matrix::rowSums(sce@assays[["counts"]] > 0)
        sce = sce[cell.counts > min.cells.per.gene, ]
        
        feature.counts = Matrix::colSums(sce@assays[["counts"]] > 0)
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


rename.sce.rows <- function(sce, from = "ENSEMBL", to = "SYMBOL", species = "human") {
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
