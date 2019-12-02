run.ACTIONet <- function(sce, k_max = 20, layout.compactness = 50, thread_no = 8, epsilon = 3, LC = 1, arch.specificity.z = -1, core.z = 3, 
    sce.data.attr = "logcounts", sym_method = "AND", scale.initial.coordinates = TRUE, reduction_slot = "S_r", batch = NULL, batch.correction.rounds = 3, 
    batch.lambda = 1, k_min = 2, n_epochs = 500, compute.core = F, compute.signature = T, specificity.mode = "sparse") {
    require(Matrix)
    require(igraph)
    require(ACTIONet)
    
    if (!(sce.data.attr %in% names(SummarizedExperiment::assays(sce)))) {
        R.utils::printf("Attribute %s is not an assay of the input sce\n", sce.data.attr)
        return()
    }
    
    # Run ACTION
    if (!is.null(batch)) {
        batch.vec = as.numeric(factor(batch))
        ACTION.out = runACTION_withBatch(t(reducedDims(sce)[[reduction_slot]]), batch.vec, k_min = k_min, k_max = k_max, max_correction_rounds = batch.correction.rounds, 
            lambda = batch.lambda, numThreads = thread_no)
        
    } else {
        ACTION.out = runACTION(t(reducedDims(sce)[[reduction_slot]]), k_min = k_min, k_max = k_max, thread_no = thread_no)
    }
    
    # Reconstruct archetypes in the original space
    reconstruct.out = reconstructArchetypes(as(SummarizedExperiment::assays(sce)[[sce.data.attr]], "sparseMatrix"), ACTION.out$C, ACTION.out$H, z_threshold = arch.specificity.z)
    rownames(reconstruct.out$archetype_profile) = rownames(sce)
    
    # Build ACTIONet
    set.seed(0)
    build.out = buildAdaptiveACTIONet(H_stacked = reconstruct.out$H_stacked, thread_no = thread_no, LC = LC, epsilon = epsilon, sym_method = sym_method, 
        auto_adjust_LC = F)
    
    # Layout ACTIONet
    if (scale.initial.coordinates == TRUE) {
        initial.coordinates = t(scale(reducedDims(sce)[[reduction_slot]]))
    } else {
        initial.coordinates = t(reducedDims(sce)[[reduction_slot]])
    }
    
    vis.out = layoutACTIONet(build.out$ACTIONet, S_r = initial.coordinates, compactness_level = layout.compactness, n_epochs = n_epochs)
    
    # Construct igraph object
    ACTIONet = graph_from_adjacency_matrix(build.out$ACTIONet, mode = "undirected", weighted = TRUE)
    coor = vis.out$coordinates
    coor3D = vis.out$coordinates_3D
    V(ACTIONet)$x = coor[, 1]
    V(ACTIONet)$y = coor[, 2]
    V(ACTIONet)$x3D = coor3D[, 1]
    V(ACTIONet)$y3D = coor3D[, 2]
    V(ACTIONet)$z3D = coor3D[, 3]
    V(ACTIONet)$color = rgb(vis.out$colors)
    
    
    arch.Lab = t(reconstruct.out$C_stacked) %*% grDevices::convertColor(color = vis.out$colors, from = "sRGB", to = "Lab")
    arch.colors = rgb(grDevices::convertColor(color = arch.Lab, from = "Lab", to = "sRGB"))
    arch.coordinates = t(reconstruct.out$C_stacked) %*% vis.out$coordinates
    arch.coordinates_3D = t(reconstruct.out$C_stacked) %*% vis.out$coordinates_3D
    arch.vis.out = list(colors = arch.colors, coordinates = arch.coordinates, coordinates_3D = arch.coordinates_3D)
    
    ACTIONet.out = list(ACTION.out = ACTION.out, reconstruct.out = reconstruct.out, build.out = build.out, vis.out = vis.out, ACTIONet = ACTIONet, 
        arch.vis.out = arch.vis.out)
    
    # Add signature profile
    if(compute.signature == TRUE) {
		signature.profile = construct.archetype.signature.profile(sce = sce, ACTIONet.out = ACTIONet.out, reduction_slot = reduction_slot)
		ACTIONet.out$signature.profile = signature.profile    
	}
	
    # Old method of computing core archetypes -obsolete
    if(compute.core == TRUE) {
		core.out = identify.core.archetypes(ACTIONet.out, core.z)
		ACTIONet.out$core.out = core.out
	}
    
    ACTIONet.out$archetype.differential.signature = compute.archetype.feature.specificity(ACTIONet.out, sce, mode = specificity.mode, sce.data.attr = sce.data.attr)
    
    ACTIONet.out = add.archetype.labels(ACTIONet.out, k_min = k_min, k_max = k_max)
    
    print("ready to unify")
    
    ACTIONet.out$unification.out = unify.cell.states(ACTIONet.out, sce, reduction_slot = reduction_slot, sce.data.attr = sce.data.attr)    
    
    if( ('cell.hashtag' %in% names(colData(sce))) ) {
		cells = sce$cell.hashtag
	} else if ( !is.null(colnames(sce)) ) {
		cells = colnames(sce)
	} else {
		cells = as.character(1:ncol(sce))
	}

	# Compute overall connectivity of nodes
	connectivity = compute.cell.connectivity(ACTIONet.out)
    V(ACTIONet.out$ACTIONet)$connectivity = connectivity

	
    ACTIONet.out$log = list(genes = rownames(sce), cells = cells, time = Sys.time())
    
    return(ACTIONet.out)
}

reconstruct.ACTIONet <- function(ACTIONet.out, sce, compactness_level = 50, thread_no = 8, epsilon = 3, LC = 1, auto_adjust_LC = FALSE, 
    n_epochs = 500, sym_method = "AND", scale.initial.coordinates = TRUE, reduction_slot = "S_r") {
    # Build ACTIONet
    set.seed(0)
    
    if (!(reduction_slot %in% names(reducedDims(sce)))) {
        R.utils::printf("%s is not in ReducedDims of sce\n", reduction_slot)
        return()
    }
    build.out = buildAdaptiveACTIONet(H_stacked = ACTIONet.out$reconstruct.out$H_stacked, thread_no = thread_no, LC = LC, auto_adjust_LC = auto_adjust_LC, 
        epsilon = epsilon, sym_method = sym_method)
    ACTIONet.out$build.out = build.out
    
    
    # Layout ACTIONet
    if (scale.initial.coordinates == TRUE) {
        initial.coordinates = t(scale(reducedDims(sce)[[reduction_slot]]))
    } else {
        initial.coordinates = t(reducedDims(sce)[[reduction_slot]])
    }
    
    
    vis.out = layoutACTIONet(build.out$ACTIONet, S_r = initial.coordinates, compactness_level = compactness_level, n_epochs = n_epochs)
    ACTIONet.out$vis.out = vis.out
    
    # Construct igraph object
    ACTIONet = graph_from_adjacency_matrix(build.out$ACTIONet, mode = "undirected", weighted = TRUE)
    coor = vis.out$coordinates
    coor3D = vis.out$coordinates_3D
    V(ACTIONet)$x = coor[, 1]
    V(ACTIONet)$y = coor[, 2]
    V(ACTIONet)$x3D = coor3D[, 1]
    V(ACTIONet)$y3D = coor3D[, 2]
    V(ACTIONet)$z3D = coor3D[, 3]
    V(ACTIONet)$color = rgb(vis.out$colors)
    ACTIONet.out$ACTIONet = ACTIONet
    
    arch.Lab = t(ACTIONet.out$reconstruct.out$C_stacked) %*% grDevices::convertColor(color = vis.out$colors, from = "sRGB", to = "Lab")
    arch.colors = rgb(grDevices::convertColor(color = arch.Lab, from = "Lab", to = "sRGB"))
    arch.coordinates = t(ACTIONet.out$reconstruct.out$C_stacked) %*% vis.out$coordinates
    arch.coordinates_3D = t(ACTIONet.out$reconstruct.out$C_stacked) %*% vis.out$coordinates_3D
    arch.vis.out = list(colors = arch.colors, coordinates = arch.coordinates, coordinates_3D = arch.coordinates_3D)
    ACTIONet.out$arch.vis.out = arch.vis.out
    
    
    return(ACTIONet.out)
}

rerun.layout <- function(ACTIONet.out, sce, compactness_level = 50, thread_no = 8, n_epochs = 500, scale.initial.coordinates = TRUE, 
    reduction_slot = "S_r") {
    # Build ACTIONet
    set.seed(0)
    
    build.out = ACTIONet.out$build.out
    
    # Layout ACTIONet
    if (scale.initial.coordinates == TRUE) {
        initial.coordinates = t(scale(reducedDims(sce)[[reduction_slot]]))
    } else {
        initial.coordinates = t(reducedDims(sce)[[reduction_slot]])
    }
    
    vis.out = layoutACTIONet(build.out$ACTIONet, S_r = initial.coordinates, compactness_level = compactness_level, n_epochs = n_epochs)
    ACTIONet.out$vis.out = vis.out
    
    # Construct igraph object
    ACTIONet = graph_from_adjacency_matrix(build.out$ACTIONet, mode = "undirected", weighted = TRUE)
    coor = vis.out$coordinates
    coor3D = vis.out$coordinates_3D
    V(ACTIONet)$x = coor[, 1]
    V(ACTIONet)$y = coor[, 2]
    V(ACTIONet)$x3D = coor3D[, 1]
    V(ACTIONet)$y3D = coor3D[, 2]
    V(ACTIONet)$z3D = coor3D[, 3]
    V(ACTIONet)$color = rgb(vis.out$colors)
    ACTIONet.out$ACTIONet = ACTIONet
    
    arch.Lab = t(ACTIONet.out$reconstruct.out$C_stacked) %*% grDevices::convertColor(color = vis.out$colors, from = "sRGB", to = "Lab")
    arch.colors = rgb(grDevices::convertColor(color = arch.Lab, from = "Lab", to = "sRGB"))
    arch.coordinates = t(ACTIONet.out$reconstruct.out$C_stacked) %*% vis.out$coordinates
    arch.coordinates_3D = t(ACTIONet.out$reconstruct.out$C_stacked) %*% vis.out$coordinates_3D
    arch.vis.out = list(colors = arch.colors, coordinates = arch.coordinates, coordinates_3D = arch.coordinates_3D)
    ACTIONet.out$arch.vis.out = arch.vis.out
    
    return(ACTIONet.out)
}

rehash.cells <- function(sce, salt = NULL) {
    if (is.null(salt)) {
        salt = as.character(Sys.time())
    }
    
    h = hashid_settings(salt = salt, min_length = 8)
    cell.hashtags = sapply(1:ncol(sce), function(i) ENCODE(i, h))
    sce$cell.hashtag = cell.hashtags
    
    metadata(sce)$tagging.time = time.tag
    
    return(sce)
}

remove.cells <- function(ACTIONet.out, filtered.cells, force = TRUE) {
    require(igraph)
    
    if(is.character(filtered.cells)) {
		filtered.cells = match(intersect(filtered.cells, ACTIONet.out$log$cells), ACTIONet.out$log$cells)
	}
    
    core.cells = which(rowSums(ACTIONet.out$reconstruct.out$C_stacked) > 0)
    if (!force) 
        selected.cells = sort(unique(union(core.cells, setdiff(1:dim(ACTIONet.out$reconstruct.out$C_stacked)[1], filtered.cells)))) 
	else selected.cells = sort(unique(setdiff(1:dim(ACTIONet.out$reconstruct.out$C_stacked)[1], filtered.cells)))
    
    ACTIONet.out.pruned = ACTIONet.out
    
    ## Update ACTION.out
    C_trace = lapply(ACTIONet.out.pruned$ACTION.out$C, function(C) {
        if (is.null(C)) 
            return(C)
        
        subC = C[selected.cells, ]
        cs = Matrix::colSums(subC)
        cs[cs == 0] = 1
        subC = scale(subC, center = FALSE, scale = cs)
        return(subC)
    })
    ACTIONet.out.pruned$ACTION.out$C = C_trace
    
    H_trace = lapply(ACTIONet.out.pruned$ACTION.out$H, function(H) {
        if (is.null(H)) 
            return(H)
        
        subH = H[, selected.cells]
        return(subH)
    })
    ACTIONet.out.pruned$ACTION.out$H = H_trace
    
    
    ## Update reconstruct.out
    C_stacked.pruned = ACTIONet.out.pruned$reconstruct.out$C_stacked[selected.cells, ]
    cs = Matrix::colSums(C_stacked.pruned)
    cs[cs == 0] = 1
    C_stacked.pruned = scale(C_stacked.pruned, center = FALSE, scale = cs)
    ACTIONet.out.pruned$reconstruct.out$C_stacked = C_stacked.pruned
    
    
    H_stacked.pruned = ACTIONet.out.pruned$reconstruct.out$H_stacked[, selected.cells]
    ACTIONet.out.pruned$reconstruct.out$H_stacked = H_stacked.pruned
    
    ## Update build.out
    ACTIONet.out$build.out$ACTIONet = ACTIONet.out$build.out$ACTIONet[selected.cells, selected.cells]
    ACTIONet.out$build.out$ACTIONet_asym = ACTIONet.out$build.out$ACTIONet_asym[selected.cells, selected.cells]
    
    # Update ACTIONet
    ACTIONet.pruned = induced.subgraph(ACTIONet.out.pruned$ACTIONet, V(ACTIONet.out.pruned$ACTIONet)[selected.cells])
    ACTIONet.out.pruned$ACTIONet = ACTIONet.pruned
    
    
    # Update vis.out
    ACTIONet.out.pruned$vis.out$coordinates = ACTIONet.out.pruned$vis.out$coordinates[selected.cells, ]
    ACTIONet.out.pruned$vis.out$coordinates_3D = ACTIONet.out.pruned$vis.out$coordinates_3D[selected.cells, ]
    ACTIONet.out.pruned$vis.out$colors = ACTIONet.out.pruned$vis.out$colors[selected.cells, ]
    
    
    ACTIONet.out.pruned$build.out$ACTIONet = ACTIONet.out.pruned$build.out$ACTIONet[selected.cells, selected.cells]
    ACTIONet.out.pruned$build.out$ACTIONet_asym = ACTIONet.out.pruned$build.out$ACTIONet_asym[selected.cells, selected.cells]
    
    
    ## Update core.out
    ACTIONet.out.pruned$core.out$H = ACTIONet.out.pruned$core.out$H[, selected.cells]
    
    # Update unification.out
    ACTIONet.out.pruned$unification.out$C.core = ACTIONet.out.pruned$unification.out$C.core[selected.cells, ]
    cs = Matrix::colSums(ACTIONet.out.pruned$unification.out$C.core)
    cs[cs == 0] = 1
    ACTIONet.out.pruned$unification.out$C.core = scale(ACTIONet.out.pruned$unification.out$C.core, center = F, scale = cs)
    
    ACTIONet.out.pruned$unification.out$H.core = ACTIONet.out.pruned$unification.out$H.core[, selected.cells]
    ACTIONet.out.pruned$unification.out$assignments.core = ACTIONet.out.pruned$unification.out$assignments.core[selected.cells]
    ACTIONet.out.pruned$unification.out$assignments.confidence.core = ACTIONet.out.pruned$unification.out$assignments.confidence.core[selected.cells]

	# Update existing cell annotations
	for(annotation.name in names(ACTIONet.out.pruned$annotations)) {
		X = ACTIONet.out.pruned$annotations[[annotation.name]]
		
		X$Labels = X$Labels[selected.cells]
		X$Labels.confidence = X$Labels.confidence[selected.cells]
		X$cells = selected.cells
		
		ACTIONet.out.pruned$annotations[[annotation.name]] = X
		
		if(! is.null(X$highlight) ) {
			ACTIONet.out.pruned = highlight.annotations(ACTIONet.out.pruned, annotation.name)
		}
	}
    
    ## Update ACTIONet log
    ACTIONet.out.pruned$log$cells = ACTIONet.out.pruned$log$cells[selected.cells]
    
    
    return(ACTIONet.out.pruned)
}

cluster.ACTIONet.using.decomposition <- function(ACTIONet.out, annotation.name = NULL) {
	if(! ('unification.out' %in% names(ACTIONet.out)) ) {
		print("unification.out is not in ACTIONet.out. Please run unify.cell.states() first.")
		return(ACTIONet.out)
	}
	clusters = as.numeric(ACTIONet.out$unification.out$assignments.core)
    names(clusters) = paste("Cluster", as.character(clusters), sep = " ")
	
	if(! ('annotations' %in% names(ACTIONet.out)) ) {
		ACTIONet.out$annotations = list()
	}

	time.stamp = as.character(Sys.time())
	if(is.null(annotation.name)) {
		annotation.name = sprintf('Clustering_%s', time.stamp)
	}
	h = hashid_settings(salt = time.stamp, min_length = 8)
	annotation.hashtag = ENCODE(length(ACTIONet.out$annotations)+1, h)
	
	res = list(Labels = clusters, Labels.confidence = NULL, DE.profile = NULL, highlight = NULL, cells = ACTIONet.out$log$cells, time.stamp = time.stamp, annotation.name = annotation.hashtag, type = "cluster.ACTIONet.using.decomposition")
	
	cmd = sprintf("ACTIONet.out$annotations$\"%s\" = res", annotation.name)	
	eval(parse(text=cmd))
		
    return(ACTIONet.out)
}

cluster.ACTIONet <- function(ACTIONet.out, annotation.name = NULL, resolution_parameter = 0.5, arch.init = TRUE) {
    if (arch.init == TRUE) {
		if(! ('unification.out' %in% names(ACTIONet.out)) ) {
			print("unification.out is not in ACTIONet.out. Please run unify.cell.states() first.")
			return(ACTIONet.out)
		}
        print("Perform archetype-based initialization")
		initial.clusters = ACTIONet.out$unification.out$assignments.core
        clusters = as.numeric(unsigned_cluster(ACTIONet.out$build.out$ACTIONet, resolution_parameter, 0, initial.clusters))
    } else {
        print("Perform default initialization")
        clusters = as.numeric(unsigned_cluster(ACTIONet.out$build.out$ACTIONet, resolution_parameter, 0))
    }    
    names(clusters) = paste("Cluster", as.character(clusters), sep = " ")
	
	if(! ('annotations' %in% names(ACTIONet.out)) ) {
		ACTIONet.out$annotations = list()
	}

	time.stamp = as.character(Sys.time())
	if(is.null(annotation.name)) {
		annotation.name = sprintf('Clustering_%s', time.stamp)
	}
	h = hashid_settings(salt = time.stamp, min_length = 8)
	annotation.hashtag = ENCODE(length(ACTIONet.out$annotations)+1, h)
	
	res = list(Labels = clusters, Labels.confidence = NULL, DE.profile = NULL, highlight = NULL, cells = ACTIONet.out$log$cells, time.stamp = time.stamp, annotation.name = annotation.hashtag, type = "cluster.ACTIONet")
	
	cmd = sprintf("ACTIONet.out$annotations$\"%s\" = res", annotation.name)	
	eval(parse(text=cmd))

		
    return(ACTIONet.out)
}


prune.cell.scores <- function(ACTIONet.out, scores, alpha_val = 0.9, transform = FALSE, deg.scale = FALSE) {
    if (is.igraph(ACTIONet.out)) 
        ACTIONet = ACTIONet.out else ACTIONet = ACTIONet.out$ACTIONet
    
    if (length(V(ACTIONet)) != length(scores)) {
        R.utils::printf("Number of cells scores (%d) doesn't match the number of vertices in the ACTIONet (%d)\n", length(scores), 
            length(V(ACTIONet)))
    }
    G = as(get.adjacency(ACTIONet, attr = "weight"), "dgTMatrix")
    
    if (transform == TRUE) 
        scores = exp(scores)
    
    scores[scores < 0] = 0
    u = as.matrix(scores/sum(scores))
    if (alpha_val != 0) {
        smoothed.scores = batchPR(G, u, alpha_val, 1)
        x = smoothed.scores[, 1]
    } else {
        x = u
    }
    
    
    if (deg.scale == TRUE) {
        degs = Matrix::colSums(G)
        s = x/degs
        s[degs == 0] = 0
    } else {
        s = x
    }
    cond = sweepcut(G, s)
    
    selected.cells = order(s, decreasing = TRUE)[1:which.min(cond)]
    pruned.cells = setdiff(1:length(x), selected.cells)
    x[pruned.cells] = 0
    
    
    return(x)
}

extract.archetype.associated.cells <- function(ACTIONet.out, archetype.no, alpha_val = 0.9) {
    h = as.numeric(ACTIONet.out$reconstruct.out$H_stacked[archetype.no, ])
    
    h = exp(h/mean(h))
    
    x = prune.cell.scores(ACTIONet.out, h, alpha_val)
    
    idx = which(x > 0)
    
    return(ACTIONet.out$log$cells[idx])
}


impute.genes.using.ACTIONet <- function(ACTIONet.out, sce, genes, alpha_val = 0.9, thread_no = 8, prune = FALSE, rescale = FALSE, expr.slot = "logcounts") {
    require(igraph)
    
    genes = unique(genes)
    
    if (is.igraph(ACTIONet.out)) 
        ACTIONet = ACTIONet.out else ACTIONet = ACTIONet.out$ACTIONet
    
    if (length(V(ACTIONet)) != dim(sce)[2]) {
        R.utils::printf("Number of cells in the input sce (%d) doesn't match the number of vertices in the ACTIONet (%d)\n", dim(sce)[2], 
            length(V(ACTIONet)))
    }
    G = as(get.adjacency(ACTIONet, attr = "weight"), "dgTMatrix")
    
    
    matched.genes = intersect(genes, rownames(sce))
    matched.idx = match(matched.genes, rownames(sce))
    
    # Smooth/impute gene expressions
    if (!(expr.slot %in% names(SummarizedExperiment::assays(sce)))) {
        R.utils::printf("%s is not in assays of sce\n", expr.slot)
    }
    
    if (length(matched.idx) > 1) {
        raw.gene.expression = Matrix::t(as(SummarizedExperiment::assays(sce)[[expr.slot]][matched.idx, ], "dgTMatrix"))
        U = raw.gene.expression
        U[U < 0] = 0
        cs = Matrix::colSums(U)
        U = as.matrix(Matrix::sparseMatrix(i = U@i + 1, j = U@j + 1, x = U@x/cs[U@j + 1], dims = dim(U)))
        U = as.matrix(U[, cs > 0])
        gg = matched.genes[cs > 0]
    } else {
        raw.gene.expression = matrix(SummarizedExperiment::assays(sce)[[expr.slot]][matched.idx, ])
        U = raw.gene.expression/sum(raw.gene.expression)
        gg = matched.genes
    }
    
	imputed.gene.expression = batchPR(G, U, alpha_val, thread_no)
    
    imputed.gene.expression[is.na(imputed.gene.expression)] = 0
    
    # Prune values
    if (prune == TRUE) {
        imputed.gene.expression = apply(imputed.gene.expression, 2, function(x) {
            cond = sweepcut(ACTIONet.out$build.out$ACTIONet, x)
            idx = which.min(cond)
            
            perm = order(x, decreasing = TRUE)
            x[perm[(idx + 1):length(x)]] = 0
            
            return(x)
        })
    }
    
    # rescale
    if (rescale) {
        imputed.gene.expression = sapply(1:dim(imputed.gene.expression)[2], function(col) {
            x = raw.gene.expression[, col]
            y = imputed.gene.expression[, col]
            
            x.Q = quantile(x, 1)
            y.Q = quantile(y, 1)
            
            if (y.Q == 0) {
                return(array(0, length(x)))
            }
            
            y = y * x.Q/y.Q
            
            y[y > max(x)] = max(x)
            
            return(y)
        })
    }
    
    colnames(imputed.gene.expression) = gg
    
    return(imputed.gene.expression)
}

impute.genes.using.archetype <- function(ACTIONet.out, genes, prune = FALSE) {
    require(igraph)
    
    genes = unique(genes)

	profile = 
    
    matched.genes = intersect(genes, rownames(sce))
    matched.idx = match(matched.genes, rownames(sce))
    
    # Smooth/impute gene expressions
    if (!(expr.slot %in% names(SummarizedExperiment::assays(sce)))) {
        R.utils::printf("%s is not in assays of sce\n", expr.slot)
    }
    
    if (length(matched.idx) > 1) {
        raw.gene.expression = Matrix::t(as(SummarizedExperiment::assays(sce)[[expr.slot]][matched.idx, ], "dgTMatrix"))
        U = raw.gene.expression
        U[U < 0] = 0
        cs = Matrix::colSums(U)
        U = as.matrix(Matrix::sparseMatrix(i = U@i + 1, j = U@j + 1, x = U@x/cs[U@j + 1], dims = dim(U)))
        U = as.matrix(U[, cs > 0])
        gg = matched.genes[cs > 0]
    } else {
        raw.gene.expression = matrix(SummarizedExperiment::assays(sce)[[expr.slot]][matched.idx, ])
        U = raw.gene.expression/sum(raw.gene.expression)
        gg = matched.genes
    }
    
	imputed.gene.expression = batchPR(G, U, alpha_val, thread_no)
    
    imputed.gene.expression[is.na(imputed.gene.expression)] = 0
    
    # Prune values
    if (prune == TRUE) {
        imputed.gene.expression = apply(imputed.gene.expression, 2, function(x) {
            cond = sweepcut(ACTIONet.out$build.out$ACTIONet, x)
            idx = which.min(cond)
            
            perm = order(x, decreasing = TRUE)
            x[perm[(idx + 1):length(x)]] = 0
            
            return(x)
        })
    }
    
    # rescale
    if (rescale) {
        imputed.gene.expression = sapply(1:dim(imputed.gene.expression)[2], function(col) {
            x = raw.gene.expression[, col]
            y = imputed.gene.expression[, col]
            
            x.Q = quantile(x, 1)
            y.Q = quantile(y, 1)
            
            if (y.Q == 0) {
                return(array(0, length(x)))
            }
            
            y = y * x.Q/y.Q
            
            y[y > max(x)] = max(x)
            
            return(y)
        })
    }
    
    colnames(imputed.gene.expression) = gg
    
    return(imputed.gene.expression)
}

assess.label.local.enrichment <- function(P, Labels) {
	if( is.null(names(Labels)) ){
		names(Labels) = as.character(Labels)
	}
    counts = table(Labels)
    p = counts/sum(counts)
	Annot = names(Labels)[match(as.numeric(names(counts)), Labels)]
    
    X = sapply(names(p), function(label) {
        x = as.numeric(Matrix::sparseVector(x = 1, i = which(Labels == label), length = length(Labels)))
    })
	colnames(X) = Annot
	
    Exp = array(1, nrow(P)) %*% t(p)
    Obs = as(P %*% X, "dgTMatrix")
    
    # Need to rescale due to missing values within the neighborhood
    rs = Matrix::rowSums(Obs)
    Obs = Matrix::sparseMatrix(i = Obs@i + 1, j = Obs@j + 1, x = Obs@x/rs[Obs@i + 1], dims = dim(Obs))
    
    Lambda = Obs - Exp
    
    
    w2 = Matrix::rowSums(P^2)
    Nu = w2 %*% t(p)
    
    a = as.numeric(qlcMatrix::rowMax(P)) %*% t(array(1, length(p)))
    
    
    logPval = (Lambda^2)/(2 * (Nu + (a * Lambda)/3))
    logPval[Lambda < 0] = 0
    logPval[is.na(logPval)] = 0
    
    logPval = as.matrix(logPval)
    
    colnames(logPval) = Annot
    

    max.idx = apply(logPval, 1, which.max)
    updated.Labels = as.numeric(names(p))[max.idx]
    names(updated.Labels) = Annot[max.idx]
    
    updated.Labels.conf = apply(logPval, 1, max)
    
    res = list(Labels = updated.Labels, Labels.confidence = updated.Labels.conf, Enrichment = logPval)
    
    return(res)
}

infer.missing.cell.annotations <- function(ACTIONet.out, annotation.in, annotation.out, double.stochastic = FALSE, max_iter = 3, adjust.levels = T, highlight = T) {
 	Adj = get.adjacency(ACTIONet.out$ACTIONet, attr = "weight")
    A = as(Adj, "dgTMatrix")
    eps = 1e-16
    rs = Matrix::rowSums(A)
    P = sparseMatrix(i = A@i + 1, j = A@j + 1, x = A@x/rs[A@i + 1], dims = dim(A))
    if (double.stochastic == TRUE) {
        w = sqrt(Matrix::colSums(P) + eps)
        W = P %*% Matrix::Diagonal(x = 1/w, n = length(w))
        P = W %*% Matrix::t(W)
    }
    
    
	Labels = preprocess.labels(ACTIONet.out, annotation.in)			
	if(is.null(Labels)) {
		return(ACTIONet.out)
	}
	
	Annot = sort(unique(Labels))
	idx = match(Annot, Labels)
	names(Annot) = names(Labels)[idx]
    
    na.mask = is.na(Labels)
    
    i = 1
    while (sum(na.mask) > 0) {
        R.utils::printf("iter %d\n", i)
    	
        new.Labels = assess.label.local.enrichment(P, Labels)
        
        mask = na.mask & (new.Labels$Labels.confidence > 3 + log(length(Labels)))
        Labels[mask] = new.Labels$Labels[mask]

        na.mask = is.na(Labels)
        if (i == max_iter) 
            break
        
        i = i + 1
    }
    new.Labels = assess.label.local.enrichment(P, Labels)
    Labels[na.mask] = new.Labels$Labels[na.mask]
	Labels.conf = new.Labels$Labels.confidence
    
    updated.Labels = as.numeric(Labels)
    names(updated.Labels) = names(Annot)[match(Labels, Annot)]
    if(adjust.levels == T) {
		Labels = reannotate.labels(ACTIONet.out, updated.Labels)
	} else {
		Labels = updated.Labels
	}
    

   	if(! ('annotations' %in% names(ACTIONet.out)) ) {
		ACTIONet.out$annotations = list()
	}

	time.stamp = as.character(Sys.time())
	if(is.null(annotation.out)) {
		annotation.out = sprintf('InferredMissingLabels_%s', time.stamp)
	}
	h = hashid_settings(salt = time.stamp, min_length = 8)
	annotation.hashtag = ENCODE(length(ACTIONet.out$annotations)+1, h)
	
	res = list(Labels = Labels, Labels.confidence = Labels.conf, DE.profile = NULL, highlight = NULL, cells = ACTIONet.out$log$cells, time.stamp = time.stamp, annotation.name = annotation.hashtag, type = "cluster.ACTIONet")
	
	cmd = sprintf("ACTIONet.out$annotations$\"%s\" = res", annotation.out)	
	eval(parse(text=cmd))    

	if( highlight == T ) {
		print("Adding annotation highlights")
		ACTIONet.out = highlight.annotations(ACTIONet.out, annotation.name = annotation.out)			
	}
	
	return(ACTIONet.out)
}


correct.cell.annotations <- function(ACTIONet.out, annotation.in, annotation.out, LFR.threshold = 2, double.stochastic = FALSE, max_iter = 3, adjust.levels = T, min.cells = 30, highlight = T) {
 	Adj = get.adjacency(ACTIONet.out$ACTIONet, attr = "weight")
    A = as(Adj, "dgTMatrix")
    eps = 1e-16
    rs = Matrix::rowSums(A)
    P = sparseMatrix(i = A@i + 1, j = A@j + 1, x = A@x/rs[A@i + 1], dims = dim(A))
    if (double.stochastic == TRUE) {
        w = sqrt(Matrix::colSums(P) + eps)
        W = P %*% Matrix::Diagonal(x = 1/w, n = length(w))
        P = W %*% Matrix::t(W)
    }
    
    
	Labels = preprocess.labels(ACTIONet.out, annotation.in)			
	if(is.null(Labels)) {
		return(ACTIONet.out)
	}
	
	# Prunes "trivial" annotations and merges them to larger ones
	counts = table(Labels)
	Labels[Labels %in% as.numeric(names(counts)[counts < min.cells])] = NA
	mask = is.na(Labels)
	if(sum(mask) > 0) {
		Labels[mask] = NA
		tmp.label = paste(annotation.in, "pruned", sep = "_")
		ACTIONet.out = infer.missing.cell.annotations(ACTIONet.out, annotation.in = Labels, annotation.out = tmp.label, adjust.levels = T, highligh = F)
		Labels = preprocess.labels(ACTIONet.out, tmp.label)			
		ACTIONet.out$annotations = ACTIONet.out$annotations[-which(names(ACTIONet.out$annotations) == tmp.label)]	
	}
	
	Annot = sort(unique(Labels))
	idx = match(Annot, Labels)
	names(Annot) = names(Labels)[idx]
        
    
	for(i in 1:max_iter) {    
        R.utils::printf("iter %d\n", i)
    	
        new.Labels = assess.label.local.enrichment(P, Labels)
        Enrichment = new.Labels$Enrichment
        curr.enrichment = sapply(1:nrow(Enrichment), function(k) Enrichment[k, names(Labels)[k]])
        
        Diff.LFR = log2((new.Labels$Labels.confidence / curr.enrichment))
        Diff.LFR[is.na(Diff.LFR)] = 0

        Labels[Diff.LFR > LFR.threshold] = new.Labels$Labels[Diff.LFR > LFR.threshold]
        names(Labels) = Annot[Labels]
    }
	Labels.conf = new.Labels$Labels.confidence
    updated.Labels = as.numeric(Labels)
    names(updated.Labels) = names(Annot)[match(Labels, Annot)]
    
    if(adjust.levels == T) {
		Labels = reannotate.labels(ACTIONet.out, updated.Labels)
	} else {
		Labels = updated.Labels
	}
    
   	if(! ('annotations' %in% names(ACTIONet.out)) ) {
		ACTIONet.out$annotations = list()
	}

	time.stamp = as.character(Sys.time())
	if(is.null(annotation.out)) {
		annotation.out = sprintf('InferredMissingLabels_%s', time.stamp)
	}
	h = hashid_settings(salt = time.stamp, min_length = 8)
	annotation.hashtag = ENCODE(length(ACTIONet.out$annotations)+1, h)
	
	res = list(Labels = Labels, Labels.confidence = Labels.conf, DE.profile = NULL, highlight = NULL, cells = ACTIONet.out$log$cells, time.stamp = time.stamp, annotation.name = annotation.hashtag, type = "cluster.ACTIONet")
	
	cmd = sprintf("ACTIONet.out$annotations$\"%s\" = res", annotation.out)	
	eval(parse(text=cmd))    

	if( highlight == T ) {
		print("Adding annotation highlights")
		ACTIONet.out = highlight.annotations(ACTIONet.out, annotation.name = annotation.out)			
	}
		
	return(ACTIONet.out)
}
