reduce.sce <- function(sce, norm.method = "default", reduced_dim = 50, max.iter = 5) {
    require(scran)
    require(ACTIONet)
    
    if (is.null(rownames(sce.norm))) {
        rownames(sce.norm) = sapply(1:nrow(sce), function(i) sprintf("Gene%d", i))
    }
    if (!("cell.hashtag" %in% colnames(colData(sce.norm)))) {
        print("tagging cells")
        time.tag = Sys.time()
        
        h = hashid_settings(salt = as.character(time.tag), min_length = 8)
        cell.hashtags = sapply(1:ncol(sce.norm), function(i) ENCODE(i, h))
        sce.norm$cell.hashtag = cell.hashtags
        
        colData(sce.norm)$original.colnames = colnames(sce.norm)
        colnames(sce.norm) = sce.norm$cell.hashtag
        
        metadata(sce.norm)$tagging.time = time.tag
    }

    if (!("logcounts" %in% names(sce@assays))) {
        print("Normalizing sce object")
        
        sce.norm = normalize.sce(sce, norm.method)
		rownames(sce.norm@assays[["counts"]]) = rownames(sce.norm@assays[["logcounts"]]) = rownames(sce.norm)
		colnames(sce.norm@assays[["counts"]]) = colnames(sce.norm@assays[["logcounts"]]) = colnames(sce.norm)
    } else {
        sce.norm = sce
		rownames(sce.norm@assays[["logcounts"]]) = rownames(sce.norm)
		colnames(sce.norm@assays[["logcounts"]]) = colnames(sce.norm)
    }
    
    # Make sure rownames are consistent Check if it doesn't have rownames
    
    
    print("Running main reduction")
    suppressWarnings({
        reduction.out = reduceGeneExpression(as(sce.norm@assays[["logcounts"]], "sparseMatrix"), reduced_dim = reduced_dim, method = 1, 
            iters = max.iter)
    })
    
    S_r = t(reduction.out$S_r)
    rownames(S_r) = colnames(sce.norm)
    colnames(S_r) = sapply(1:ncol(S_r), function(i) sprintf("PC%d", i))
    
    SingleCellExperiment::reducedDim(sce.norm, "S_r") <- S_r
    
    V = reduction.out$V
    rownames(V) = rownames(sce.norm)
    colnames(V) = sapply(1:dim(V)[2], function(i) sprintf("PC%d", i))
    
    X = rowData(sce.norm)
    PC.idx = -grep("^PC", colnames(X))
    if (length(PC.idx) > 0) 
        X = X[, PC.idx]
    rowData(sce.norm) = cbind(V, X)
    
    metadata(sce.norm)$reduction.time = Sys.time()
    return(sce.norm)
}

batch.correct.sce.Harmony <- function(sce, batch.vec) {
    require(harmony)
    sce@reducedDims$S_r = harmony::HarmonyMatrix(sce@reducedDims$S_r, batch.vec, do_pca = FALSE)
    return(sce)
}


reduce.and.batch.correct.sce.Harmony <- function(sce, batch.vec = NULL, norm.method = "default", reduced_dim = 50, max.iter = 5) {
    if (is.null(batch.vec)) {
        print("You need to provide the batch vector/attr")
        return(sce)
    }
    
    sce = reduce.sce(sce, reduced_dim = reduced_dim, max.iter = max.iter, norm.method = norm.method)
    batch.correct.sce.Harmony(sce, batch.vec)
}

reduce.and.batch.correct.sce.MNN <- function(sce, batch.vec = NULL, norm.method = "scran", reduced_dim = 50, MNN.k = 20) {
    require(scran)
    require(ACTIONet)
    
    if (is.null(batch.vec)) {
        print("You need to provide the batch vector/attr")
        return(sce)
    }
    
    sce = clearSpikes(sce)
    counts(sce) = as(counts(sce), "sparseMatrix")
    
    
    IDX = split(1:ncol(sce), batch.vec)
    
    sce.list = lapply(IDX, function(idx) suppressWarnings(normalize.sce(sce[, idx], norm.method)))
    sce.list.norm = do.call(scran::multiBatchNorm, sce.list)
    
    
    # Sort based on 'complexity'
    perm = order(sapply(sce.list.norm, function(sce) dim(sce)[2]), decreasing = TRUE)
    sce.list.norm = sce.list.norm[perm]
    
    
    sce.norm = do.call(SingleCellExperiment::cbind, sce.list.norm)
    sce.norm@assays[["logcounts"]] = as(sce.norm@assays[["logcounts"]], "sparseMatrix")
    
    metadata(sce.norm)$normalization.method = "multiBatchNorm"
    metadata(sce.norm)$normalization.time = Sys.time()
    
    
    # Make sure rownames are consistent Check if it doesn't have rownames
    if (is.null(rownames(sce.norm))) {
        rownames(sce.norm) = sapply(1:nrow(sce), function(i) sprintf("Gene%d", i))
    }
    rownames(sce.norm@assays[["logcounts"]]) = rownames(counts(sce.norm)) = rownames(sce.norm)
    
    if (!("cell.hashtag" %in% colnames(colData(sce.norm)))) {
        time.tag = Sys.time()
        
        h = hashid_settings(salt = as.character(time.tag), min_length = 8)
        cell.hashtags = sapply(1:ncol(sce), function(i) ENCODE(i, h))
        sce.norm$cell.hashtag = cell.hashtags
        
        metadata(sce.norm)$tagging.time = time.tag
    }
    if (is.null(colnames(sce.norm))) {
        # Check if it doesn't have colnames
        colnames(sce.norm) = sce.norm$cell.hashtag
    }
    colnames(sce.norm@assays[["logcounts"]]) = colnames(counts(sce.norm)) = colnames(sce.norm)
    
    
    set.seed(0)
    mnn.out <- do.call(scran::fastMNN, c(sce.list.norm, list(k = MNN.k, d = reduced_dim, auto.order = FALSE, approximate = TRUE, cos.norm = FALSE)))
    
    S_r = mnn.out$corrected
    rownames(S_r) = colnames(sce.norm)
    colnames(S_r) = sapply(1:ncol(S_r), function(i) sprintf("PC%d", i))
    
    
    SingleCellExperiment::reducedDim(sce.norm, "S_r") <- S_r
    
    V = mnn.out$rotation
    rownames(V) = rownames(sce.norm)
    colnames(V) = sapply(1:dim(V)[2], function(i) sprintf("PC%d", i))
    
    colnames(V) = sapply(1:dim(V)[2], function(i) sprintf("PC%d", i))
    
    X = rowData(sce.norm)
    PC.idx = grep("^PC", colnames(X))
    if (length(PC.idx) > 0) 
        X = X[, -PC.idx]
    
    rowData(sce.norm) = cbind(V, X)
    
    metadata(sce.norm)$reduction.time = Sys.time()
    return(sce.norm)
}


