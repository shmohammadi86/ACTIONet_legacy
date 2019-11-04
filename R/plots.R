ACTIONet.color.bank = c("#1f78b4", "#33a02c", "#e31a1c", "#6a3d9a", "#d95f02", "#e7298a", "#feaf16", "#a6761d", "#1b9e77", "#888888", "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#006E71", "#000080", "#8C564BFF", "#800000", "#e6194b", "#ffe119", "#AA4488")

# From: https://github.com/r3fang/SnapATAC/blob/master/R/plottings-utilities.R
ACTIONet.color.bank1 = c("#E31A1C", "#FFD700", "#771122", "#777711", "#1F78B4", "#68228B", "#AAAA44", "#60CC52", "#771155", "#DDDD77", 
    "#774411", "#AA7744", "#AA4455", "#117744", "#000080", "#44AA77", "#AA4488", "#DDAA77", "#D9D9D9", "#BC80BD", "#FFED6F", "#7FC97F", 
    "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", 
    "#E6AB02", "#A6761D", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", 
    "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", 
    "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FFFF33", "#A65628", 
    "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", 
    "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5")

# From: https://github.com/brianhie/scanorama/blob/master/scanorama/scanorama.py
ACTIONet.color.bank2 = c("#377eb8", "#ff7f00", "#4daf4a", "#f781bf", "#a65628", "#984ea3", "#999999", "#e41a1c", "#dede00", "#ffe119", 
    "#e6194b", "#ffbea3", "#911eb4", "#46f0f0", "#f032e6", "#d2f53c", "#008080", "#e6beff", "#aa6e28", "#800000", "#aaffc3", "#808000", 
    "#ffd8b1", "#000080", "#808080", "#fabebe", "#a3f4ff")

# From: https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors/
ACTIONet.color.bank3 = c("#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231", "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe", 
    "#008080", "#e6beff", "#9a6324", "#fffac8", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075", "#808080", "#ffffff", "#000000")


plot.ACTIONet <- function(ACTIONet.out, labels = NULL, transparency.attr = NULL, trans.z.threshold = -1, trans.fact = 2, 
	node.size = 1, CPal = ACTIONet.color.bank, add.text = FALSE, text.halo.width = 0.1, label.text.size = 0.8, 
    suppress.legend = FALSE, legend.pos = "bottomright", add.states = F, title = "") {
    
    node.size = node.size * 0.5
    
    if (is.igraph(ACTIONet.out)) 
        ACTIONet = ACTIONet.out else ACTIONet = ACTIONet.out$ACTIONet
    
    coors = cbind(V(ACTIONet)$x, V(ACTIONet)$y)
    
	labels = preprocess.labels(ACTIONet.out, labels)
	if(is.null(labels)) {
        vCol = V(ACTIONet)$color
        Annot = NULL
	} else {
		Annot = names(labels)[match(sort(unique(labels)), labels)]		
		if(length(CPal) > 1) {
            Pal = CPal[1:length(Annot)]			
		} else {
            Pal = ggpubr::get_palette(CPal, length(Annot))
		}
        names(Pal) = Annot
        vCol = Pal[names(labels)]
	}

    if (!is.null(transparency.attr)) {
        z = (transparency.attr - median(transparency.attr))/mad(transparency.attr)
        beta = 1/(1 + exp(-trans.fact * (z - trans.z.threshold)))
        beta[z > trans.z.threshold] = 1
        beta = beta^trans.fact
        
        vCol.border = scales::alpha(colorspace::darken(vCol, 0.5), beta)
        vCol = scales::alpha(vCol, beta)
    } else {
        vCol.border = colorspace::darken(vCol, 0.5)
    }

    graphics::plot(coors[, c(1, 2)], pch = 21, cex = node.size, bg = vCol, col = vCol.border, axes = FALSE, xlab = "", ylab = "", main = title)

	if(add.states == T) {
		par(new=TRUE)
		
		M = as(ACTIONet.out$unification.out$C.core, 'sparseMatrix')
		cs = Matrix::colSums(M)
		M = scale(M, center = FALSE, scale = cs)
		
	    cell.Lab = grDevices::convertColor(color = t(col2rgb(vCol)/256), from = "sRGB", to = "Lab")	    
	    core.Lab = t(t(cell.Lab) %*% M)
	    core.colors = rgb(grDevices::convertColor(color = core.Lab, from = "Lab", to = "sRGB"))
		core.colors = colorspace::lighten(core.colors, 0.1)
		
		core.coors = t(t(ACTIONet.out$vis.out$coordinates) %*% M)		
	    
	    graphics::plot(core.coors, pch = 25, cex = 4*node.size, bg = core.colors, col = "#eeeeee", axes = FALSE, xlab = "", ylab = "")
	}    
    
    
    if ( add.text == T & (!is.null(Annot)) ) {
        centroids = t(sapply(Annot, function(l) {
            idx = which(names(labels) == l)
            if(length(idx) == 1) {
				return(as.numeric(coors[idx, ]))
			} 
            sub.coors = coors[idx, ]
            D = as.matrix(dist(sub.coors))
            stats = Matrix::rowMeans(D)
            min.idx = which.min(stats)
			anchor.coor = as.numeric(sub.coors[min.idx, ])            
			return(anchor.coor)
        }))
        textHalo(x = centroids[, 1], y = centroids[, 2], labels = Annot, col = colorspace::darken(Pal, 0.5), bg = "#eeeeee", r = text.halo.width, cex = label.text.size)
    }
    
    if ( (suppress.legend == FALSE) & !is.null(Annot) ) {
        legend(legend.pos, legend = Annot, fill = Pal, cex = 0.5, bty = "n")
    }    
}


plot.ACTIONet.3D <- function(ACTIONet.out, labels = NULL, transparency.attr = NULL, trans.z.threshold = -1, trans.fact = 2, node.size = 1, CPal = ACTIONet.color.bank) {
    require(ggplot2)
    require(ggpubr)
    require(threejs)
    
    if (is.igraph(ACTIONet.out)) 
        ACTIONet = ACTIONet.out else ACTIONet = ACTIONet.out$ACTIONet
    
    nV = length(V(ACTIONet))
    coor = cbind(V(ACTIONet)$x3D, V(ACTIONet)$y3D, V(ACTIONet)$z3D)
    
    node.size = node.size * 0.2
	labels = preprocess.labels(ACTIONet.out, labels)
	if(is.null(labels)) {
        vCol = V(ACTIONet)$color
        add.text = F
        suppress.legend = T
	} else {
		Annot = names(labels)[match(sort(unique(labels)), labels)]		
		if(is.vector(CPal)) {
            Pal = CPal[1:length(Annot)]			
		} else {
            Pal = ggpubr::get_palette(CPal, length(Annot))
		}
        names(Pal) = Annot
        vCol = Pal[names(labels)]
	}

    if (!is.null(transparency.attr)) {
        z = (transparency.attr - median(transparency.attr))/mad(transparency.attr)
        beta = 1/(1 + exp(-trans.fact * (z - trans.z.threshold)))
        beta[z > trans.z.threshold] = 1
        beta = beta^trans.fact
        
        vCol.border = scales::alpha(colorspace::darken(vCol, 0.5), beta)
        vCol = scales::alpha(vCol, beta)
    } else {
        vCol.border = colorspace::darken(vCol, 0.5)
    }
        
    scatterplot3js(x = coor[, 1], y = coor[, 2], z = coor[, 3], axis.scales = FALSE, size = node.size, axis = F, grid = F, color = vCol, 
        stroke = vCol.border, bg = "black")
}



plot.ACTIONet.gene.view <- function(ACTIONet.out, top.genes = 5, CPal = NULL, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH") {
    if( !('unification.out' %in% names(ACTIONet.out)) ) {
		print('Error in plot.ACTIONet.gene.view: "unification.out" is not in ACTIONet.out. Please run unify.cell.states() first.')
		return()
	}
	gene.enrichment.table = as.matrix(ACTIONet.out$unification.out$DE.core@assays[["significance"]])
	
	filtered.rows = grep(blacklist.pattern, rownames(gene.enrichment.table))
	if(length(filtered.rows) > 0)
		gene.enrichment.table = gene.enrichment.table[-filtered.rows, ]
	
	GT = apply(gene.enrichment.table, 2, function(x) rownames(gene.enrichment.table)[order(x, decreasing = T)[1:top.genes]])
	selected.genes = sort(unique(as.character(GT)))
	
	M = Matrix::t(as(ACTIONet.out$unification.out$H.core, 'sparseMatrix'))
	cs = Matrix::colSums(M)
	M = scale(M, center = FALSE, scale = cs)
	
	core.coors = t(t(ACTIONet.out$vis.out$coordinates) %*% M)
	X = t(gene.enrichment.table[selected.genes, ])
	cs = colSums(X)
	cs[cs == 0] = 1
	X = scale(X, center = F, scale = cs)
	gene.coors = t(X) %*% core.coors

    if (is.null(CPal)) {
        Pal = ACTIONet.out$unification.out$Pal
    } else {
    	if(length(CPal) == 1) {
            Pal = ggpubr::get_palette(CPal, length(ACTIONet.out$unification.out$Pal))
    	} else {
            Pal = CPal[1:length(ACTIONet.out$unification.out$Pal)]
    	}
    }
    
    core.Lab = grDevices::convertColor(color = t(col2rgb(Pal)/256), from = "sRGB", to = "Lab")
    gene.color.Lab = t(X) %*% core.Lab
    gene.colors = rgb(grDevices::convertColor(color = gene.color.Lab, from = "Lab", to = "sRGB"))
    names(gene.colors) = selected.genes


    genes.df = data.frame(gene = selected.genes, x = gene.coors[, 1], y = gene.coors[, 2])
    
    require(ggrepel)
    require(ggplot2)
    p <- ggplot(genes.df, aes(x, y, label = gene, color = gene)) + scale_colour_manual(values = gene.colors) + geom_point(show.legend = FALSE) + geom_label_repel(show.legend = FALSE, force = 5) + theme_void()
    
    plot(p)
}


plot.ACTIONet.interactive <- function(ACTIONet.out, labels = NULL, transparency.attr = NULL, trans.z.threshold = -1, trans.fact = 2, 
	node.size = 1, CPal = ACTIONet.color.bank, enrichment.table = NULL, top.features = 7, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH", threeD = FALSE, title = "ACTIONet") {
    require(plotly)
    require(ACTIONet)
    
    if (is.igraph(ACTIONet.out)) 
        ACTIONet = ACTIONet.out else ACTIONet = ACTIONet.out$ACTIONet
    
    nV = length(V(ACTIONet))
    node.size = node.size * 5
    
	labels = preprocess.labels(ACTIONet.out, labels)
	if(is.null(labels)) {
        vCol = V(ACTIONet)$color
        Annot = NULL
	} else {
		Annot = names(labels)[match(sort(unique(labels)), labels)]		
		if(length(CPal) > 1) {
            Pal = CPal[1:length(Annot)]			
		} else {
            Pal = ggpubr::get_palette(CPal, length(Annot))
		}
        names(Pal) = Annot
        vCol = Pal[names(labels)]
	}

    if (!is.null(transparency.attr)) {
        z = (transparency.attr - median(transparency.attr))/mad(transparency.attr)
        beta = 1/(1 + exp(-trans.fact * (z - trans.z.threshold)))
        beta[z > trans.z.threshold] = 1
        beta = beta^trans.fact
        
        vCol.border = scales::alpha(colorspace::darken(vCol, 0.5), beta)
        vCol = scales::alpha(vCol, beta)
    } else {
        vCol.border = colorspace::darken(vCol, 0.5)
    }

    
	
	if( !is.null(enrichment.table) ) {
		if(ncol(enrichment.table) == nV) {
			cell.scores = Matrix::t(enrichment.table)
		} else if( (nrow(enrichment.table) != nV) ) {
			H = ACTIONet.out$unification.out$H.core
			if( (nrow(enrichment.table) == nrow(H)) | (ncol(enrichment.table) == nrow(H)) ) {
				cell.scores = map.cell.scores.from.archetype.enrichment(ACTIONet.out, enrichment.table)				
			} else {
				cell.scores = NULL
			}
		} else {
			cell.scores = enrichment.table
		}
	} else {
		if( ('unification.out' %in% names(ACTIONet.out)) ) {
			temp.enrichment.table = as.matrix(ACTIONet.out$unification.out$DE.core@assays[["significance"]])			
			if( !is.null(row.names(temp.enrichment.table)) ) {
				filtered.rows = grep(blacklist.pattern, rownames(temp.enrichment.table))
				if(length(filtered.rows) > 0)
					enrichment.table = temp.enrichment.table[-filtered.rows, ]
				else
					enrichment.table = temp.enrichment.table

				GT = apply(enrichment.table, 2, function(x) rownames(enrichment.table)[order(x, decreasing = T)[1:min(100, nrow(enrichment.table))]])
				selected.features = sort(unique(as.character(GT)))
				
				cell.scores = t(enrichment.table[selected.features, ] %*% ACTIONet.out$unification.out$H.core)
			}
		} else {
			cell.scores = NULL
		}	
	}

    if ( !is.null(cell.scores) ) {
		selected.features = colnames(cell.scores)
		node.annotations = apply(cell.scores, 1, function(x) paste(selected.features[order(x, decreasing = T)[1:top.features]], collapse = '\n'))
		# node.annotations = sapply(1:length(ACTIONet.out$log$cells), function(i) sprintf('(%s) %s', ACTIONet.out$log$cells[[i]], node.annotations[[i]]) )
	} else {
		# node.annotations = sapply(1:length(ACTIONet.out$log$cells), function(i) sprintf('(%s)', ACTIONet.out$log$cells[[i]]) )
		node.annotations = rep('', nV)
	}
    
    # Setup visualization parameters
    sketch.graph = ACTIONet.out$ACTIONet
    sketch.graph = delete.edges(sketch.graph, E(sketch.graph))
    
    node.data <- get.data.frame(sketch.graph, what = "vertices")
    edge.data <- get.data.frame(sketch.graph, what = "edges")
    
    Nv <- dim(node.data)[1]
    Ne <- dim(edge.data)[1]
    
    edge_shapes <- list()
    
    # Adjust parameters
    node.data$size = node.size
    
    axis <- list(title = "", showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
    
    node.data$vCol = vCol
    
    if(threeD == TRUE) {
		if(is.null(Annot)) {
		    node.data$vCol = vCol
			network <- plot_ly(node.data, x = ~x3D, y = ~y3D, z = ~z3D, opacity = 1, marker = list(color = ~vCol, size = ~size, opacity = 1, alpha = 1, line = list(width = 0.1 * node.size, alpha = 0.5, color = ~vCol.border)), text = node.annotations, mode = "markers", hoverinfo = "text", type = "scatter3d", showlegend = FALSE)		
			p <- plotly::layout(network, title = title, shapes = edge_shapes, xaxis = axis, yaxis = axis)
		} else {
		    node.data$type = factor(names(labels), levels = Annot)
			network <- plot_ly(node.data, x = ~x3D, y = ~y3D, z = ~z3D, opacity = 1, color = ~type, colors = Pal, marker = list(size = ~size, opacity = 1, alpha = 1, line = list(width = 0.1 * node.size, alpha = 0.5, color = ~vCol.border)), text = node.annotations, mode = "markers", hoverinfo = "text", type = "scatter3d")		
			p <- plotly::layout(network, title = title, shapes = edge_shapes, xaxis = axis, yaxis = axis, showlegend = TRUE, legend = list(marker = list(marker.size = 10)))
		}		
	} else {				
		if(is.null(Annot)) {
		    node.data$vCol = vCol
			network <- plot_ly(node.data, x = ~x, y = ~y, opacity = 1, marker = list(color = ~vCol, size = ~size, opacity = 1, alpha = 1, line = list(width = 0.1 * node.size, alpha = 0.5, color = ~vCol.border)), text = node.annotations, mode = "markers", type = "scatter", hoverinfo = "text", showlegend = FALSE)
			p <- plotly::layout(network, title = title, shapes = edge_shapes, xaxis = axis, yaxis = axis)
		} else {
		    node.data$type = factor(names(labels), levels = Annot)
			network <- plot_ly(node.data, x = ~x, y = ~y, opacity = 1, color = ~type, colors = Pal, marker = list(size = ~size, opacity = 1, alpha = 1, line = list(width = 0.1 * node.size, alpha = 0.5, color = ~vCol.border)), text = node.annotations, mode = "markers", type = "scatter", hoverinfo = "text")
			p <- plotly::layout(network, title = title, shapes = edge_shapes, xaxis = axis, yaxis = axis, showlegend = TRUE, legend = list(marker = list(marker.size = 10)))
		}
	}
		
    p
    
}

plot.marker.boxplot <- function(ACTIONet.out, sce, marker.genes, Labels, node.size = 3, CPal = ACTIONet.color.bank, export_path = NA, 
    thread_no = 8, prune = FALSE, scale.factor = 2, p.threshold = NA) {
    require(igraph)
    require(ACTIONet)
    require(ggpubr)
    
    if (!is.factor(Labels)) {
        Labels = factor(Labels, levels = sort(unique(Labels)))
    }
    
    if (!sum(sapply(marker.genes, length) != 1) & is.null(names(marker.genes))) {
        names(marker.genes) = marker.genes
    }
    
    gg = unique(unlist(marker.genes))
    all.marker.genes = sort(intersect(gg, rownames(sce)))
    
    imputed.markers = impute.genes.using.ACTIONet(ACTIONet.out, sce, all.marker.genes, prune = prune)
    
    if (!is.na(p.threshold)) {
        wilc.out = presto::wilcoxauc(t(imputed.markers), Labels)
        
        sig.mask = wilc.out$logFC > 0 & wilc.out$auc > 0.5 & wilc.out$padj < p.threshold
        sig.genes = sort(unique(wilc.out$feature[sig.mask]))
        
        
        updated.names = sapply(colnames(imputed.markers), function(gene) ifelse(gene %in% sig.genes, sprintf("%s*", gene), gene))
        
        colnames(imputed.markers) = updated.names
    }
    
    imputed.markers.df = as.data.frame(log2(imputed.markers * nrow(imputed.markers)))
    imputed.markers.df$Celltype = Labels
    
    
    
    
    if (1 < length(CPal)) {
        Pal = CPal[1:length(levels(Labels))]
    } else {
        Pal = ggpubr::get_palette(CPal, length(levels(Labels)))
    }
    names(Pal) = levels(Labels)
    
    sapply(colnames(imputed.markers), function(gene.name) {
        gp = ggboxplot(imputed.markers.df, x = "Celltype", y = gene.name, fill = "Celltype", palette = Pal) + theme(axis.text.x = element_text(angle = 90, 
            hjust = 1))
        print(gp)
        
        if (!is.na(export_path)) {
            fname = sprintf("%s/%s.pdf", export_path, gene.name)
            pdf(fname)
            print(gp)
            dev.off()
        }
        
    })
    
}


plot.marker.dotplot <- function(ACTIONet.out, sce, marker.genes, Labels, CPal = "YlOrRd", thread_no = 8, prune = FALSE, font.size = 1, 
    p.threshold = NA, keep.order = FALSE, row.rescale = FALSE) {
    require(ACTIONet)
    library(corrplot)
    library(seriation)
    
    if (!is.factor(Labels)) {
        Labels = factor(Labels, levels = sort(unique(Labels)))
    }
    
    if (!sum(sapply(marker.genes, length) != 1) & is.null(names(marker.genes))) {
        names(marker.genes) = marker.genes
    }
    
    gg = unique(unlist(marker.genes))
    all.marker.genes = sort(intersect(gg, rownames(sce)))
    
    imputed.markers = impute.genes.using.ACTIONet(ACTIONet.out, sce, all.marker.genes, prune = prune)
    
    
    X = apply(imputed.markers, 2, function(x) return((x - min(x))/(max(x) - min(x))))
    
    IDX = split(1:nrow(imputed.markers), Labels)
    mean.expr = sapply(IDX, function(idx) as.numeric(Matrix::colMeans(X[idx, ])))
    rownames(mean.expr) = colnames(imputed.markers)
    
    
    if (keep.order == TRUE) {
        perm = order(match(rownames(mean.expr), marker.genes))
    } else {
        set.seed(0)
        perm = seriation::get_order(seriate(mean.expr, "BEA_TSP"))
    }
    
    
    if (!is.na(p.threshold)) {
        wilc.out = presto::wilcoxauc(t(X), Labels)
        
        sig.mask = wilc.out$logFC > 0 & wilc.out$auc > 0.5 & wilc.out$padj < p.threshold
        sig.genes = sort(unique(wilc.out$feature[sig.mask]))
        
        
        updated.names = sapply(rownames(mean.expr), function(gene) ifelse(gene %in% sig.genes, sprintf("%s*", gene), gene))
        
        rownames(mean.expr) = updated.names
    }
    
    
    Pal = ggpubr::get_palette(CPal, 11)
    
    if (row.rescale == TRUE) {
        X = t(scale(t(mean.expr[perm, ])))
    } else {
        X = mean.expr[perm, ]
    }
    corrplot(X, is.corr = FALSE, method = "circle", tl.col = "black", cl.pos = "n", col = Pal, tl.cex = font.size, cl.cex = font.size)
}

plot.ACTIONet.gradient <- function(ACTIONet.out, x, max.update.iter = 3, CPal = "magma", node.size = 3, prune = FALSE, nonparameteric = FALSE, 
    transparency.attr = NA, trans.fact = 1, title = NA, alpha.val = 0.5) {
    require(igraph)
    require(colorspace)
    require(viridis)
    
    add.vertex.shape("fcircle", clip = igraph.shape.noclip, plot = mycircle, parameters = list(vertex.frame.color = 1, vertex.frame.width = 1))
    
    # Normalize z = scale(x) z[z > 3] = 3 x = exp(z) x = x / sum(x) x = page_rank(ACTIONet.out$ACTIONet, personalized = x, damping =
    # 0.15)$vector
    
    
    x[x < 0] = 0
    if (prune == TRUE) {
        x = prune.cell.scores(ACTIONet.out, x, alpha_val = alpha.val, transform = FALSE)
    }
    
    if (CPal %in% c("inferno", "magma", "viridis", "BlGrRd", "RdYlBu", "Spectral")) {
        Pal_grad = switch(CPal, inferno = inferno(500, alpha = 0.8), magma = magma(500, alpha = 0.8), viridis = viridis(500, alpha = 0.8), 
            BlGrRd = colorRampPalette(c("blue", "grey", "red"))(500), Spectral = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, 
                name = "Spectral"))))(100), RdYlBu = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"))))(100))
    } else {
        lg = rgb(0.95, 0.95, 0.95)
        Pal_grad = colorRampPalette(c(lg, CPal))(500)
    }
    
    if (is.igraph(ACTIONet.out)) 
        ACTIONet = ACTIONet.out else ACTIONet = ACTIONet.out$ACTIONet
    
    
    sketch.graph = ACTIONet
    sketch.graph = delete.edges(sketch.graph, E(sketch.graph))
    
    V(sketch.graph)$name = ""
    V(sketch.graph)$shape = "fcircle"
    
    vCol = rgb(0.95, 0.95, 0.95)
    
    if (nonparameteric == TRUE) {
        vCol = colorspace::darken((scales::col_numeric(Pal_grad, domain = NULL))(rank(x)), 0.2)
    } else {
        vCol = colorspace::darken((scales::col_numeric(Pal_grad, domain = NULL))(x), 0.2)
    }
    # vCol[x == 0] = rgb(0.95, 0.95, 0.95)
    
    
    if (is.numeric(transparency.attr)) {
        # beta = 1 / (1 + exp(-trans.fact*(scale(transparency.attr)+1)))
        z = (transparency.attr - median(transparency.attr))/mad(transparency.attr)
        beta = 1/(1 + exp(-trans.fact * (z + 1)))
        beta[z > -1] = 1
        # beta = (transparency.attr - min(transparency.attr)) / (max(transparency.attr) - min(transparency.attr))
        beta = beta^5
        
        # plot(density(1-beta)) vCol = colorspace::lighten(vCol, 1-beta, method = 'relative', space = 'HCL') vCol.border =
        # colorspace::darken(vCol, 0.2*beta)#colorspace::lighten(vCol.border, 1-beta, method = 'relative', space = 'HLS')
        
        vCol.border = scales::alpha(colorspace::darken(vCol, 0.3), beta)
        vCol = scales::alpha(vCol, beta)
        
        # vCol.border = colorspace::darken(vCol, 0.5*beta)#colorspace::lighten(vCol.border, 1-beta, method = 'relative', space = 'HLS')
        
    } else {
        vCol.border = colorspace::darken(vCol, 0.3)
    }
    
    V(sketch.graph)$size = node.size
    V(sketch.graph)$frame.width = 0.33 * node.size
    
    V(sketch.graph)$color = vCol
    V(sketch.graph)$frame.color = vCol.border
    
    
    plot(sketch.graph, main = title)
}

visualize.markers <- function(ACTIONet.out, sce, marker.genes, max.update.iter = 3, CPal = ACTIONet.color.bank, node.size = 3, adjust.node.size = FALSE, 
    alpha_val = 0.9, export_path = NA, prune = TRUE) {
    require(igraph)
    add.vertex.shape("fcircle", clip = igraph.shape.noclip, plot = mycircle, parameters = list(vertex.frame.color = 1, vertex.frame.width = 1))
    
    
    if (!sum(sapply(marker.genes, length) != 1) & is.null(names(marker.genes))) {
        names(marker.genes) = marker.genes
    }
    Pal = ggpubr::get_palette(CPal, length(names(marker.genes)))
    names(Pal) = names(marker.genes)
    
    gg = unique(unlist(marker.genes))
    all.marker.genes = sort(intersect(gg, rownames(sce)))
    
    
    imputed.marker.expression = impute.genes.using.ACTIONet(ACTIONet.out, sce, all.marker.genes, prune = prune, alpha_val = alpha_val)
    
    # return(imputed.marker.expression)
    
    ACTIONet = ACTIONet.out$ACTIONet
    sketch.graph = ACTIONet
    sketch.graph = delete.edges(sketch.graph, E(sketch.graph))
    
    V(sketch.graph)$name = ""
    V(sketch.graph)$shape = "fcircle"
    
    
    eps = 1e-16
    A = as(get.adjacency(ACTIONet, attr = "weight"), "dgTMatrix")
    rs = Matrix::rowSums(A)
    P = sparseMatrix(i = A@i + 1, j = A@j + 1, x = A@x/rs[A@i + 1], dims = dim(A))
    
    gradPal = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"))))(100)
    
    lapply(all.marker.genes, function(gene) {
        print(gene)
        if (!(gene %in% colnames(imputed.marker.expression))) 
            return()
        
        idx = which(sapply(marker.genes, function(gs) gene %in% gs))[1]
        celltype.name = names(marker.genes)[idx]
        
        x = imputed.marker.expression[, gene]
        
        
        # Normalize between [0, 1]
        x = (x - min(x))/(max(x) - min(x))
        
        # Find effective # of nonzeros using participation ratio
        V(sketch.graph)$color = rgb(0.95, 0.95, 0.95)
        V(sketch.graph)$frame.color = rgb(0.9, 0.9, 0.9)
        
        lg = rgb(0.95, 0.95, 0.95)
        Pal_grad = colorRampPalette(c(lg, Pal[[celltype.name]]))(500)
        
        
        V(sketch.graph)$color = colorspace::darken((scales::col_numeric(Pal_grad, domain = NULL))(scale(x)), 0.15)
        V(sketch.graph)$frame.color = colorspace::darken(V(sketch.graph)$color, 0.15)
        
        
        
        if (adjust.node.size == TRUE) 
            V(sketch.graph)$size = 3 * node.size * x else V(sketch.graph)$size = node.size
        
        V(sketch.graph)$frame.width = node.size * 0.1
        
        plot(sketch.graph, main = gene)
        
        
        if (!is.na(export_path)) {
            fname = sprintf("%s/%s.pdf", export_path, ifelse(celltype.name == gene, gene, sprintf("%s_%s", celltype.name, gene)))
            pdf(fname)
            plot(sketch.graph, main = gene)
            dev.off()
        }
    })
}
