reduce.sce <- function(sce, norm.method = "default", reduced_dim = 50,  reduction_name = "S_r", max.iter = 5, return_V = FALSE, BPPARAM = SerialParam()) {
    require(ACTIONet)

    sce.norm = sce

    if (is.null(rownames(sce.norm))) {
        rownames(sce.norm) = sapply(1:nrow(sce.norm), function(i) sprintf("Gene%d", i))
    }

    if (!("logcounts" %in% names(SummarizedExperiment::assays(sce.norm)))) {
        print("Normalizing sce object")

        sce.norm = normalize.sce(sce.norm, norm.method, BPPARAM)
        rownames(SummarizedExperiment::assays(sce.norm)[["counts"]]) = rownames(SummarizedExperiment::assays(sce.norm)[["logcounts"]]) = rownames(sce.norm)
        colnames(SummarizedExperiment::assays(sce.norm)[["counts"]]) = colnames(SummarizedExperiment::assays(sce.norm)[["logcounts"]]) = colnames(sce.norm)
    } else {
      rownames(SummarizedExperiment::assays(sce.norm)[["logcounts"]]) = rownames(sce.norm)
      colnames(SummarizedExperiment::assays(sce.norm)[["logcounts"]]) = colnames(sce.norm)
    }

    print("Running main reduction")
    suppressWarnings({
        reduction.out = reduceGeneExpression(as(SummarizedExperiment::assays(sce.norm)[["logcounts"]], "sparseMatrix"), reduced_dim = reduced_dim, method = 1,
            iters = max.iter)
    })

    S_r = t(reduction.out$S_r)
    rownames(S_r) = colnames(sce.norm)
    colnames(S_r) = sapply(1:ncol(S_r), function(i) sprintf("PC%d", i))

    SingleCellExperiment::reducedDims(sce.norm)[[reduction_name]] <- S_r

    if(return_V){
    	V = reduction.out$V
    	colnames(V) = sapply(1:dim(V)[2], function(i) sprintf('PC%d', i))
      rowData(sce.norm)[["rotation"]] = V
    	# X = rowData(sce)
    	# PC.idx = -grep('^PC', colnames(X))
    	# if(length(PC.idx) > 0)
    	# X = X[, PC.idx]
    	# rowData(sce) = cbind(V, X)
    }

    return(sce.norm)
}

batch.correct.sce.Harmony <- function(sce, batch.vec, reduction_name = "S_r", max_iter = 10) {
    require(harmony)
    reducedDims(sce)[[reduction_name]] <- harmony::HarmonyMatrix(reducedDims(sce)[[reduction_name]], batch.vec, do_pca = FALSE, max.iter.harmony = max_iter)
    return(sce)
}

reduce.and.batch.correct.sce.Harmony <- function(sce, batch.vec = NULL, norm.method = "default", reduced_dim = 50, max.iter = 10, BPPARAM = SerialParam()) {
	if( !("harmony" %in% rownames(installed.packages())) ) {
		message("You need to install harmony (https://github.com/immunogenomics/harmony) first for batch-correction.")
		return
	} else {
		library(harmony)
	}

    if (is.null(batch.vec)) {
        print("You need to provide the batch vector/attr")
        return(sce)
    }

    sce = reduce.sce(sce, reduced_dim = reduced_dim, max.iter = max.iter, norm.method = norm.method, BPPARAM = BPPARAM)
    sce = batch.correct.sce.Harmony(sce, batch.vec)
    return(sce)
}

reduce.and.batch.correct.sce.fastMNN <- function(sce, batch.attr = NULL, reduced_dim = 50, MNN.k = 20, return_V = FALSE, reduction_name = "S_r", BPPARAM = SerialParam()) {
    require(scran)
    require(ACTIONet)

    assays(sce)[["counts"]] = as(assays(sce)[["counts"]], 'sparseMatrix')


    if( length(batch.attr) ==  1) {
      IDX = split(1:dim(sce)[2], colData(sce)[[batch.attr]])
    }
    else {
      IDX = split(1:dim(sce)[2], batch.attr)
    }
BPPARAM
    sce.list = lapply(IDX, function(idx) computeSumFactors(sce[, idx], BPPARAM = BPPARAM ) )
    sce.list.norm = do.call(batchelor::multiBatchNorm, list(sce.list, BPPARAM = BPPARAM))


    # Sort based on 'complexity'
    merge_order = order(sapply(sce.list.norm, function(sce) dim(sce)[2]), decreasing = TRUE)


    sce.norm = do.call(SingleCellExperiment::cbind, sce.list.norm)
    SummarizedExperiment::assays(sce.norm)[["logcounts"]] = as(SummarizedExperiment::assays(sce.norm)[["logcounts"]], "sparseMatrix")

    set.seed(0)
    mnn.out <- do.call(batchelor::fastMNN, c(sce.list.norm, list(k = 20, d = 50, auto.merge = FALSE, merge.order = merge_order, cos.norm = FALSE, assay.type = "logcounts", BPPARAM = BPPARAM)))

    S_r = reducedDims(mnn.out)[["corrected"]]
    rownames(S_r) = colnames(sce.norm)
    colnames(S_r) = sapply(1:ncol(S_r), function(i) sprintf("PC%d", i))


    SingleCellExperiment::reducedDims(sce.norm)[[reduction_name]] <- S_r

    if(return_V){
      V = rowData(mnn.out)[["rotation"]]
      colnames(V) = sapply(1:dim(V)[2], function(i) sprintf("PC%d", i))
      rowData(sce.norm)[["rotation"]] = V
    }

    return(sce.norm)
  }

reduce.and.batch.correct.sce.fastMNN.old <- function(sce, batch.attr = NULL, reduced_dim = 50, MNN.k = 20, reduction_name = "S_r", return_V = FALSE, BPPARAM = SerialParam()) {
    require(scran)
    require(ACTIONet)

    if (is.null(batch.attr)) {
        print("You need to provide the batch vector/attr")
        return(sce)
    }

    sce = clearSpikes(sce)
    assays(sce)[["counts"]] = as(assays(sce)[["counts"]], 'sparseMatrix')


    if( length(batch.attr) ==  1) {
      IDX = split(1:dim(sce)[2], colData(sce)[[batch.attr]])
    }
    else {
      IDX = split(1:dim(sce)[2], batch.attr)
    }

    sce.list = lapply(IDX, function(idx) computeSumFactors(sce[, idx], BPPARAM = BPPARAM ) )
    sce.list.norm = do.call(scran::multiBatchNorm, sce.list)


    # Sort based on 'complexity'
    perm = order(sapply(sce.list.norm, function(sce) dim(sce)[2]), decreasing = TRUE)
    sce.list.norm = sce.list.norm[perm]


    sce.norm = do.call(SingleCellExperiment::cbind, sce.list.norm)
    SummarizedExperiment::assays(sce.norm)[["logcounts"]] = as(SummarizedExperiment::assays(sce.norm)[["logcounts"]], "sparseMatrix")

    set.seed(0)
    mnn.out <- do.call(scran::fastMNN, c(sce.list.norm, list(k = MNN.k, d = reduced_dim, auto.order = FALSE, approximate = TRUE, cos.norm = FALSE, BPPARAM = BPPARAM)))

    S_r = mnn.out$corrected
    rownames(S_r) = colnames(sce.norm)
    colnames(S_r) = sapply(1:ncol(S_r), function(i) sprintf("PC%d", i))


    SingleCellExperiment::reducedDims(sce.norm)[[reduction_name]] <- S_r

    if(return_V){
      V = mnn.out$rotation
      rownames(V) = rownames(sce.norm)
      colnames(V) = sapply(1:dim(V)[2], function(i) sprintf("PC%d", i))

      X = rowData(sce.norm)
      PC.idx = grep("^PC", colnames(X))
      if (length(PC.idx) > 0){
  			X = X[, -PC.idx]
  		}
      rowData(sce.norm) = cbind(V, X)
    }

    return(sce.norm)
	}
