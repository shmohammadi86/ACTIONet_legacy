construct.sparse.backbone <- function(ACTIONet.out, reduction.slot = "S_r", stretch.factor = 10) {
    core.index = ACTIONet.out$core.out$core.archs
    
    # Construct core-backbone X = ACTIONet.out$reconstruct.out$C_stacked[, core.index]
    X = t(ACTIONet.out$reconstruct.out$H_stacked[core.index, ])
    
    mat = t(reducedDims(sce)[[reduction.slot]]) %*% X
    
    # mat = t(ACTIONet.out_ACTION$reconstruct.out$H_stacked[ACTIONet.out$core.out$core.archs, ]) mat = orthoProject(mat,
    # Matrix::rowMeans(mat))
    
    backbone = cor(mat)
    
    backbone[backbone < 0] = 0
    diag(backbone) = 0
    
    backbone.graph = graph_from_adjacency_matrix(backbone, mode = "undirected", weighted = TRUE)
    
    # Construct t-spanner
    t = (2 * stretch.factor - 1)
    
    d = 1 - E(backbone.graph)$weight
    EL = get.edgelist(backbone.graph, names = FALSE)
    perm = order(d, decreasing = FALSE)
    
    backbone.graph.sparse = delete.edges(backbone.graph, E(backbone.graph))
    for (i in 1:length(d)) {
        u = EL[perm[i], 1]
        v = EL[perm[i], 2]
        sp = distances(backbone.graph.sparse, v = u, to = v)[1, 1]
        
        if (sp > t * d[perm[i]]) {
            backbone.graph.sparse = add.edges(backbone.graph.sparse, EL[perm[i], ], attr = list(weight = 1 - d[perm[i]]))
        }
    }
    
    backbone.graph = backbone.graph.sparse
    
    return(backbone.graph)
}

identify.core.archetypes <- function(ACTIONet.out, pruning.zscore.threshold = 3) {
    require(igraph)
    require(ACTIONet)
    
    mergeTree = mergeArchetypes(C_stacked = ACTIONet.out$reconstruct.out$C_stacked, H_stacked = ACTIONet.out$reconstruct.out$H_stacked)
    
    A = matrix(as.numeric(mergeTree < pruning.zscore.threshold & mergeTree != 0), ncol = ncol(mergeTree))
    A = A + t(A)
    merged.graph = graph_from_adjacency_matrix(A, weighted = TRUE, mode = "undirected")
    
    comps = components(merged.graph)
    core.archs = sapply(unique(comps$membership), function(c) min(which(comps$membership == c)))
    
    out.list = list(core.archs = core.archs, arch.membership = comps$membership, mergeTree = mergeTree)
    
    return(out.list)
}

normalize.adj <- function(Adj, double.stochastic = T) {
    A = as(Adj, "dgTMatrix")
    eps = 1e-16
    rs = Matrix::rowSums(A)
    P = sparseMatrix(i = A@i + 1, j = A@j + 1, x = A@x/rs[A@i + 1], dims = dim(A))
    if (double.stochastic == TRUE) {
        w = sqrt(Matrix::colSums(P) + eps)
        W = P %*% Matrix::Diagonal(x = 1/w, n = length(w))
        P = W %*% Matrix::t(W)
    }
    
    return(P)
}

plot.backbone.heatmap <- function(backbone, annotations.df, resolution = 0.5, CPal = "d3") {
    require(ComplexHeatmap)
    require(RColorBrewer)
    
    W = as.matrix(get.adjacency(backbone, attr = "weight"))
    
    
    set.seed(0)
    modules = unsigned_cluster(as(W, "sparseMatrix"), resolution_parameter = resolution)
    annotations.df$Module = modules
    
    annotations = colnames(annotations.df)
    Pals = lapply(annotations, function(annotation) {
        L = annotations.df[, annotation]
        if (!is.factor(L)) {
            L = factor(L, levels = sort(unique(L)))
        }
        
        all.levels = levels(L)
        Pal = ggpubr::get_palette(CPal, length(levels(L)))
        names(Pal) = levels(L)
        return(Pal)
    })
    names(Pals) = annotations
    
    ha_row = HeatmapAnnotation(df = annotations.df, col = Pals, show_legend = FALSE, show_annotation_name = FALSE, which = "row")
    
    ha_column = HeatmapAnnotation(df = annotations.df, col = Pals, show_legend = TRUE, show_annotation_name = TRUE, which = "column")
    
    gradPal = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
    
    ht = ComplexHeatmap::Heatmap(W, row_names_gp = gpar(fontsize = 0), column_names_gp = gpar(fontsize = 0), top_annotation = ha_column, 
        name = "z-score", left_annotation = ha_row, col = gradPal, raster_quality = 100)
    
    clusters = as.numeric(as.factor(annotations.df$Module))
    sapply(colnames(annotations.df)[1:(ncol(annotations.df) - 1)], function(attr) {
        labels = as.numeric(as.factor(annotations.df[, attr]))
        NMI = ClusterR::external_validation(labels, clusters)
        ARI = ClusterR::external_validation(labels, clusters)
        R.utils::printf("%s: NMI = %.2f, ARI = %.2f\n", attr, NMI, ARI)
    })
    
    return(ht)
}

unify.cell.states <- function(ACTIONet.out, sce, reduction_slot = "S_r", min.cor = 0.9, resolution = 1, alpha_val = 0.5, min.cells = 5, 
    thread_no = 8, sce.data.attr = "logcounts") {
    G = ACTIONet.out$build.out$ACTIONet
    
    print("Construct dependency map of cell states")
    S_r = t(sce@reducedDims[[reduction_slot]])
    C = ACTIONet.out$reconstruct.out$C_stacked
    W_r = S_r %*% C
    
    selected.states = which(min.cells <= Matrix::colSums(C > 0))
    
    A = as.matrix(cor(W_r[, selected.states]))
    diag(A) = 0
    A[A < min.cor] = 0
    
    
    print("Identify equivalent cell state classes")
    set.seed(0)
    modules = unsigned_cluster(as(A, "sparseMatrix"), resolution_parameter = resolution)
    
    print("Construct reduced multi-resolution cell state decomposition")
    IDX = split(1:length(modules), modules)
    C.core = sapply(IDX, function(idx) {
        if (length(idx) == 1) 
            c = as.numeric(C[, selected.states[idx]]) else {
            c = as.numeric(Matrix::rowMeans(C[, selected.states[idx]]))
        }
        return(c)
    })
    W.core = S_r %*% C.core
    H.core = runsimplexRegression(W.core, S_r)
    cellstates.core = assay(sce)[[sce.data.attr]] %*% C.core
    
    ## DE stuff
    print("Compute differential expression")
    
    ### Official DE
    DE.core = assess.feature.specificity(sce, H.core, sce.data.attr = sce.data.attr)
    
    print("Assign cells to cell states")
    ## Cell assignments
    assignments.core = apply(H.core, 2, which.max)
    assignments.confidence.core = apply(H.core, 2, max)
	names(assignments.core) = names(assignments.confidence.core) = ACTIONet.out$log$cells    

    
    ## Identify anchor cells
    print("Identify anchor cells")
    anchor.cells.core = apply(H.core, 1, which.max)

    
    ## Groups of archetypes
    equivalent.cellstates.core = sapply(split(1:length(modules), modules), function(idx) selected.states[idx])

	## Assign colors to archetypes
	Pal.Lab = t(C.core) %*% grDevices::convertColor(color = ACTIONet.out$vis.out$colors, from = "sRGB", to = "Lab")
	Pal = rgb(grDevices::convertColor(color = Pal.Lab, from = "Lab", to = "sRGB"))

    
    res = list(C.core = C.core, W.core = W.core, H.core = H.core, cellstates.core = cellstates.core, 
        DE.core = DE.core, assignments.core = assignments.core, assignments.confidence.core = assignments.confidence.core, anchor.cells.core = anchor.cells.core, 
        equivalent.classes.core = equivalent.cellstates.core, dependency.graph = A, selected.states = selected.states, Pal = Pal)

    
    return(res)
}

plot.ACTIONet.cell.state.map <- function(ACTIONet.out, sce, archetype.labels = NA, transparency.attr = NA, trans.fact = 1, stretch.factor = 10, 
    CPal = ACTIONet.color.bank, cex = 2, node.scale.factor = 3, reduction.slot = "S_r") {
    # Plot main ACTIONet first
    sketch.graph = ACTIONet.out$ACTIONet
    sketch.graph = delete.edges(sketch.graph, E(sketch.graph))
    nV = length(V(sketch.graph))
    
    
    vCol = colorspace::lighten("black", 0.97)
    vCol.border = colorspace::lighten("black", 0.9)
    
    if (is.numeric(transparency.attr)) {
        z = (transparency.attr - median(transparency.attr))/mad(transparency.attr)
        beta = 1/(1 + exp(-trans.fact * (z)))
        beta = beta^2
        
        vCol = scales::alpha(vCol, beta)
        vCol.border = scales::alpha(vCol.border, beta)
        
    }
    
    
    
    V(sketch.graph)$color = vCol
    V(sketch.graph)$frame.color = vCol.border
    
    
    
    
    
    V(sketch.graph)$size = cex
    
    coor = cbind(V(sketch.graph)$x, V(sketch.graph)$y)
    
    plot(sketch.graph, vertex.label = NA, layout = coor)
    
    
    # Now overlay the core backbone connectome on top
    core.index = ACTIONet.out$core.out$core.archs
    core.coor = ACTIONet.out$arch.vis.out$coordinates[core.index, ]
    
    Annot = NA
    if (is.numeric(archetype.labels)) {
        archetype.labels = as.character(archetype.labels)
        archetype.labels = factor(archetype.labels, levels = sort(unique(archetype.labels)))
        Annot = levels(archetype.labels)
        
        if (1 < length(CPal)) {
            Pal = CPal[1:length(Annot)]
        } else {
            Pal = ggpubr::get_palette(CPal, length(Annot))
        }
        names(Pal) = Annot
    } else if (is.factor(archetype.labels)) {
        Annot = levels(archetype.labels)
        if (1 < length(CPal)) {
            Pal = CPal[1:length(Annot)]
        } else {
            Pal = ggpubr::get_palette(CPal, length(Annot))
        }
        names(Pal) = Annot
    } else if (is.character(archetype.labels)) {
        archetype.labels = factor(archetype.labels, levels = sort(unique(archetype.labels)))
        Annot = levels(archetype.labels)
        if (1 < length(CPal)) {
            Pal = CPal[1:length(Annot)]
        } else {
            Pal = ggpubr::get_palette(CPal, length(Annot))
        }
        names(Pal) = Annot
    }
    if (length(Annot) > 1) {
        vCol = Pal[archetype.labels]
    } else {
        vCol = ACTIONet.out$arch.vis.out$colors
    }
    
    core.col = vCol[core.index]
    
    cs = as.numeric(table(ACTIONet.out$core.out$arch.membership))
    core.scale.factor = cs/max(cs)
    core.size = cex * (1 + node.scale.factor * core.scale.factor)
    
    
    backbone.graph = construct.sparse.backbone(ACTIONet.out, reduction.slot, stretch.factor)
    
    
    w.scaled = E(backbone.graph)$weight
    w.scaled = (w.scaled/max(w.scaled))
    
    E(backbone.graph)$width = 0.5 + 2.5 * w.scaled
    E(backbone.graph)$color = ggplot2::alpha("black", 0.25 + 0.35 * w.scaled)
    E(backbone.graph)$arrow.size = 0.1
    
    V(backbone.graph)$color = core.col
    V(backbone.graph)$size = core.size
    
    HSV = rgb2hsv(col2rgb(V(backbone.graph)$color))
    HSV[3, ] = HSV[3, ] * 0.8
    V(backbone.graph)$frame.color = apply(HSV, 2, function(v) do.call(hsv, as.list(v)))
    
    
    plot(backbone.graph, vertex.label = NA, layout = core.coor, add = T)
    if (!is.na(archetype.labels)) {
        legend("bottomright", legend = Annot, fill = Pal, cex = 0.5)
    }
}


plot.ACTIONet.cell.state.view <- function(ACTIONet.out, sce, archetype.labels = NA, transparency.attr = NA, trans.fact = 1, CPal = ACTIONet.color.bank, 
    cex = 2, node.scale.factor = 3, reduction.slot = "S_r") {
    # Plot main ACTIONet first
    sketch.graph = ACTIONet.out$ACTIONet
    sketch.graph = delete.edges(sketch.graph, E(sketch.graph))
    nV = length(V(sketch.graph))
    
    
    vCol = colorspace::lighten("black", 0.97)
    vCol.border = colorspace::lighten("black", 0.9)
    
    if (is.numeric(transparency.attr)) {
        z = (transparency.attr - median(transparency.attr))/mad(transparency.attr)
        beta = 1/(1 + exp(-trans.fact * (z)))
        beta = beta^2
        
        vCol = scales::alpha(vCol, beta)
        vCol.border = scales::alpha(vCol.border, beta)
        
    }
    
    
    
    V(sketch.graph)$color = vCol
    V(sketch.graph)$frame.color = vCol.border
    
    
    
    
    
    V(sketch.graph)$size = cex
    
    coor = cbind(V(sketch.graph)$x, V(sketch.graph)$y)
    
    plot(sketch.graph, vertex.label = NA, layout = coor)
    
    
    # Now overlay the core backbone connectome on top
    core.index = ACTIONet.out$core.out$core.archs
    core.coor = ACTIONet.out$arch.vis.out$coordinates[core.index, ]
    
    Annot = NA
    if (is.numeric(archetype.labels)) {
        archetype.labels = as.character(archetype.labels)
        archetype.labels = factor(archetype.labels, levels = sort(unique(archetype.labels)))
        Annot = levels(archetype.labels)
        
        if (1 < length(CPal)) {
            Pal = CPal[1:length(Annot)]
        } else {
            Pal = ggpubr::get_palette(CPal, length(Annot))
        }
        names(Pal) = Annot
    } else if (is.factor(archetype.labels)) {
        Annot = levels(archetype.labels)
        if (1 < length(CPal)) {
            Pal = CPal[1:length(Annot)]
        } else {
            Pal = ggpubr::get_palette(CPal, length(Annot))
        }
        names(Pal) = Annot
    } else if (is.character(archetype.labels)) {
        archetype.labels = factor(archetype.labels, levels = sort(unique(archetype.labels)))
        Annot = levels(archetype.labels)
        if (1 < length(CPal)) {
            Pal = CPal[1:length(Annot)]
        } else {
            Pal = ggpubr::get_palette(CPal, length(Annot))
        }
        names(Pal) = Annot
    }
    if (length(Annot) > 1) {
        vCol = Pal[archetype.labels]
    } else {
        vCol = ACTIONet.out$arch.vis.out$colors
    }
    
    core.col = vCol[core.index]
    
    cs = as.numeric(table(ACTIONet.out$core.out$arch.membership))
    core.scale.factor = cs/max(cs)
    core.size = cex * (1 + node.scale.factor * core.scale.factor)
    
    
    
    backbone.graph = graph.empty(n = length(core.size), directed = FALSE)
    V(backbone.graph)$color = core.col
    V(backbone.graph)$size = core.size
    
    HSV = rgb2hsv(col2rgb(V(backbone.graph)$color))
    HSV[3, ] = HSV[3, ] * 0.8
    V(backbone.graph)$frame.color = apply(HSV, 2, function(v) do.call(hsv, as.list(v)))
    
    
    plot(backbone.graph, vertex.label = NA, layout = core.coor, add = T)
    if (!is.na(archetype.labels)) {
        legend("bottomright", legend = Annot, fill = Pal, cex = 0.5)
    }
}


