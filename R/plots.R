ACTIONet.color.bank = c("#1f78b4", "#33a02c", "#e31a1c", "#6a3d9a", "#d95f02", "#e7298a", "#feaf16", "#a6761d", "#1b9e77", "#808080", "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f", "#006E71", "#000080", "#8C564BFF", "#800000", "#e6194b", "#ffe119", "#AA4488", "#f032e6", "#bcf60c")

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


plot.ACTIONet <- function(ACTIONet.out, labels = NULL, transparency.attr = NULL, trans.z.threshold = -0.5, trans.fact = 3, 
	node.size = 1, CPal = ACTIONet.color.bank, add.text = TRUE, text.halo.width = 0.1, label.text.size = 0.8, 
    suppress.legend = TRUE, legend.pos = "bottomright", add.states = F, title = "", highlight = F) {
    
    node.size = node.size * 0.5
    
    if (is.igraph(ACTIONet.out)) 
        ACTIONet = ACTIONet.out else ACTIONet = ACTIONet.out$ACTIONet
    
    coors = cbind(V(ACTIONet)$x, V(ACTIONet)$y)

	if( (length(labels) == 1) & (is.character(labels)) ) {
		annotation_name = labels
	} else {
		annotation_name = "plot.ACTIONet-tmp"
	}	
    
	labels = preprocess.labels(ACTIONet.out, labels)

	if( !is.null(labels) & (highlight == T) ) {
		if( annotation_name == "plot.ACTIONet-tmp") {
			ACTIONet.out = add.cell.annotations(ACTIONet.out, cell.annotations = labels, annotation.name = annotation_name)
		}
		
		if( is.null(ACTIONet.out$annotations[[annotation_name]]$highlight) ) {
			print("Labels are not highlighted ... generating highlight on the fly")
			ACTIONet.out = highlight.annotations(ACTIONet.out, annotation.name = annotation_name)			
		}
		label.hightlight = ACTIONet.out$annotations[[annotation_name]]$highlight
		
		if( !is.null(transparency.attr) ) {
			print(class(transparency.attr))
			print("highlight=T while transparency.attr is not NULL. Overwriting transparency.attr with highlight values")
		}
		transparency.attr = label.hightlight$connectivity.scores		
	}
	
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
    	require(wordcloud)
        centroids = t(sapply(Annot, function(l) {
            idx = which(names(labels) == l)
            if(length(idx) == 1) {
				return(as.numeric(coors[idx, ]))
			} 

            sub.coors = coors[idx, ]
            anchor.coor = as.numeric(apply(sub.coors, 2, function(x) mean(x, trim = 0.5)))
#             sub.coors.sq = sub.coors^2
# 			norm.sq = Matrix::rowSums(sub.coors.sq)
# 			anchor.idx = which.min(sapply(1:nrow(sub.coors.sq), function(i) { 
# 				dd = norm.sq[i] + norm.sq - 2* sub.coors %*% sub.coors[i, ]
# 				mean.dist.sq = median(dd)
# 				return(mean.dist.sq)
# 			}))
            
                        
            # D = as.matrix(dist(sub.coors))
            # stats = Matrix::rowMeans(D)
            # anchor.idx = which.min(stats)
            
			# anchor.coor = as.numeric(sub.coors[anchor.idx, ])            
            
			return(anchor.coor)
        }))
        layout.labels(x = centroids[, 1], y = centroids[, 2], labels = Annot, col = colorspace::darken(Pal, 0.5), bg = "#eeeeee", r = text.halo.width, cex = label.text.size)
        #wordcloud::textplot(x = centroids[, 1], y = centroids[, 2], new = F, words = Annot, col = colorspace::darken(Pal, 0.5), bg = "#eeeeee") 
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
        
    scatterplot3js(x = coor[, 1], y = coor[, 2], z = coor[, 3], axis.scales = FALSE, size = node.size, axis = F, grid = F, color = as.character(vCol), 
        stroke = as.character(vCol.border), bg = "black")
}



plot.ACTIONet.feature.view <- function(ACTIONet.out, feature.enrichment.table, top.features = 5, CPal = NULL, title = "Feature view", label.text.size = 1) {
	if(ncol(feature.enrichment.table) != nrow(ACTIONet.out$unification.out$H.core)) {
		feature.enrichment.table = Matrix::t(feature.enrichment.table)
	}
	
	if(max(feature.enrichment.table) > 50)
		feature.enrichment.table = log1p(feature.enrichment.table)

	feature.enrichment.table = doubleNorm(feature.enrichment.table)
	
	selected.features = sort(unique(as.character(apply(feature.enrichment.table, 2, function(x) rownames(feature.enrichment.table)[order(x, decreasing = T)[1:top.features]]))))
	
	M = Matrix::t(as(ACTIONet.out$unification.out$H.core, 'sparseMatrix'))
	cs = Matrix::colSums(M)
	M = scale(M, center = FALSE, scale = cs)
	
	core.coors = t(t(ACTIONet.out$vis.out$coordinates) %*% M)
	X = t(feature.enrichment.table[selected.features, ])
	cs = colSums(X)
	cs[cs == 0] = 1
	X = scale(X, center = F, scale = cs)
	feature.coors = t(X) %*% core.coors

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
    feature.color.Lab = t(X) %*% core.Lab
    feature.colors = rgb(grDevices::convertColor(color = feature.color.Lab, from = "Lab", to = "sRGB"))
    names(feature.colors) = selected.features


	x = feature.coors[, 1]
	y = feature.coors[, 2]
	words = selected.features
    plot(x, y, type = "n", col = feature.colors, axes = FALSE, xlab = "", ylab = "", main = title)
    lay <- wordlayout(x, y, words, label.text.size)
    for (i in 1:length(x)) {
        xl <- lay[i, 1]
        yl <- lay[i, 2]
        w <- lay[i, 3]
        h <- lay[i, 4]
        if (x[i] < xl || x[i] > xl + w || y[i] < yl || y[i] > 
            yl + h) {
            points(x[i], y[i], pch = 16, col = "black", cex = 0.75*label.text.size)
            nx <- xl + 0.5 * w
            ny <- yl + 0.5 * h
            lines(c(x[i], nx), c(y[i], ny), col = colorspace::darken(feature.colors[[i]], 0.5))
        }
    }
    text(lay[, 1] + 0.5 * lay[, 3], lay[, 2] + 0.5 * lay[, 4], words, col = feature.colors, cex = label.text.size, xlab = "", ylab = "", main = title)

}

plot.ACTIONet.gene.view <- function(ACTIONet.out, top.genes = 5, CPal = NULL, blacklist.pattern = "\\.|^RPL|^RPS|^MRP|^MT-|^MT|^RP|MALAT1|B2M|GAPDH", title = "Gene view") {
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

	x = gene.coors[, 1]
	y = gene.coors[, 2]
	words = selected.genes
    plot(x, y, type = "n", col = gene.colors, axes = FALSE, xlab = "", ylab = "", main = title)
    lay <- wordlayout(x, y, words, label.text.size)
    for (i in 1:length(x)) {
        xl <- lay[i, 1]
        yl <- lay[i, 2]
        w <- lay[i, 3]
        h <- lay[i, 4]
        if (x[i] < xl || x[i] > xl + w || y[i] < yl || y[i] > 
            yl + h) {
            points(x[i], y[i], pch = 16, col = "black", cex = 0.75*label.text.size)
            nx <- xl + 0.5 * w
            ny <- yl + 0.5 * h
            lines(c(x[i], nx), c(y[i], ny), col = colorspace::darken(gene.colors[[i]], 0.5))
        }
    }
    text(lay[, 1] + 0.5 * lay[, 3], lay[, 2] + 0.5 * lay[, 4], words, col = gene.colors, cex = label.text.size, xlab = "", ylab = "", main = title)
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
			} else {
				cell.scores = NULL
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
    
    if(threeD == TRUE) {
		if(is.null(Annot)) {
		    node.data$vCol = vCol
			node.data$vCol.border = vCol.border
			network <- plot_ly(node.data, x = ~x3D, y = ~y3D, z = ~z3D, opacity = 1, marker = list(color = ~vCol, size = ~size, opacity = 1, alpha = 1, line = list(width = 0.1 * node.size, alpha = 0.5, color = ~vCol.border)), text = node.annotations, mode = "markers", hoverinfo = "text", type = "scatter3d", showlegend = FALSE)		
			p <- plotly::layout(network, title = title, shapes = edge_shapes, xaxis = axis, yaxis = axis)
		} else {
			border.Pal = colorspace::darken(Pal, 0.5)
		    node.data$type = factor(names(labels), levels = Annot)
			network <- plot_ly(node.data, x = ~x3D, y = ~y3D, z = ~z3D, opacity = 1, color = ~type, colors = Pal, marker = list(size = ~size, opacity = 1, alpha = 1, line = list(width = 0.1 * node.size, alpha = 0.5, color = ~type, colors = border.Pal)), text = node.annotations, mode = "markers", hoverinfo = "text", type = "scatter3d")		
			p <- plotly::layout(network, title = title, shapes = edge_shapes, xaxis = axis, yaxis = axis, showlegend = TRUE, legend = list(marker = list(marker.size = 10)))
		}		
	} else {				
		if(is.null(Annot)) {
		    node.data$vCol = vCol
			node.data$vCol.border = vCol.border
			network <- plot_ly(node.data, x = ~x, y = ~y, marker = list(color = ~vCol, size = ~size, opacity = 1, alpha = 1, line = list(width = 0.1 * node.size, alpha = 0.5, color = ~vCol.border)), text = node.annotations, mode = "markers", type = "scatter", hoverinfo = "text", showlegend = FALSE)
			p <- plotly::layout(network, title = title, shapes = edge_shapes, xaxis = axis, yaxis = axis)
		} else {
			border.Pal = colorspace::darken(Pal, 0.5)
		    node.data$type = factor(names(labels), levels = Annot)
			network <- plot_ly(node.data, x = ~x, y = ~y, color = ~type, colors = Pal, marker = list(size = ~size, line = list(width = 0.1 * node.size, color = ~type, colors = border.Pal)), text = node.annotations, mode = "markers", type = "scatter", hoverinfo = "text")
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

plot.ACTIONet.gradient <- function(ACTIONet.out, x, transparency.attr = NULL, trans.z.threshold = -0.5, trans.fact = 3, node.size = 1, CPal = "magma", title = "", prune = F, alpha_val = 0.5, nonparameteric = FALSE) {

    node.size = node.size * 0.5

    if (is.igraph(ACTIONet.out)) 
        ACTIONet = ACTIONet.out else ACTIONet = ACTIONet.out$ACTIONet
    coors = cbind(V(ACTIONet)$x, V(ACTIONet)$y)

    NA.col = "#eeeeee"
        
    ## Create color gradient generator
    if (CPal %in% c("inferno", "magma", "viridis", "BlGrRd", "RdYlBu", "Spectral")) {
        Pal_grad = switch(CPal, inferno = inferno(500, alpha = 0.8), magma = magma(500, alpha = 0.8), viridis = viridis(500, alpha = 0.8), 
            BlGrRd = colorRampPalette(c("blue", "grey", "red"))(500), Spectral = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, 
                name = "Spectral"))))(100), RdYlBu = (grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"))))(100))
    } else {
        Pal_grad = colorRampPalette(c(NA.col, CPal))(500)
    }	

    ## Scale/prune scorees, if needed
    x[x < 0] = 0
	if(max(x) > 50)
		x = log1p(x)
	
    if (prune == TRUE) {
        x = prune.cell.scores(ACTIONet.out, x, alpha_val = alpha.val, transform = FALSE)
    }
    
    if (nonparameteric == TRUE) {
        vCol = scales::col_bin(Pal_grad, domain = NULL, bin = 100, na.color = NA.col)(rank(x))
    } else {
        vCol = scales::col_bin(Pal_grad, domain = NULL, bin = 100, na.color = NA.col)(x)
    }	

    if (!is.null(transparency.attr)) {
        z = (transparency.attr - median(transparency.attr))/mad(transparency.attr)
        beta = 1/(1 + exp(-trans.fact * (z - trans.z.threshold)))
        beta[z > trans.z.threshold] = 1
        beta = beta^trans.fact
        
        vCol.border = scales::alpha(colorspace::darken(vCol, 0.1), beta)
        vCol = scales::alpha(vCol, beta)
    } else {
        vCol.border = colorspace::darken(vCol, 0.1)
    }
	
    idx = order(x, decreasing = F)
    plot(coors[idx, 1], coors[idx, 2], bg = vCol[idx], col = vCol.border[idx], cex = node.size, pch = 21, axes = FALSE, xlab = "", ylab = "", main = title)
	
}

visualize.markers <- function(ACTIONet.out, sce, marker.genes, transparency.attr = NULL, trans.z.threshold = -0.5, trans.fact = 3, node.size = 1, CPal = ACTIONet.color.bank1,  alpha_val = 0.9, export_path = NA, prune = FALSE) {
    require(igraph)
    
    
    if (!sum(sapply(marker.genes, length) != 1) & is.null(names(marker.genes))) {
        names(marker.genes) = marker.genes
    }
    Pal = ggpubr::get_palette(CPal, length(names(marker.genes)))
    names(Pal) = names(marker.genes)
    
    gg = unique(unlist(marker.genes))
    all.marker.genes = sort(intersect(gg, rownames(sce)))
    
    
    imputed.marker.expression = impute.genes.using.ACTIONet(ACTIONet.out, sce, all.marker.genes, prune = FALSE, alpha_val = alpha_val)
    
    lapply(all.marker.genes, function(gene) {
        print(gene)
        if (!(gene %in% colnames(imputed.marker.expression))) 
            return()
        
        idx = which(sapply(marker.genes, function(gs) gene %in% gs))[1]
        celltype.name = names(marker.genes)[idx]
        
        x = imputed.marker.expression[, gene]

		plot.ACTIONet.gradient(ACTIONet.out, x, transparency.attr, trans.z.threshold, trans.fact, node.size, CPal = Pal[[celltype.name]], title = gene, prune = prune, alpha_val = alpha_val)
    })
}

plotOrdered.Heatmap <- function(W, row_title = "Cell states", col_title = "Annotations", measure_name = "Enrichment", scale = T) {
	
	if(scale == T) {
		W = t(scale(t(W)))
	}
	
	require(ComplexHeatmap)
	require(RColorBrewer)
	require(seriation)
	CC = cor(t(W))
	CC[is.na(CC)] = 0
	D = as.dist(1-CC)
	row.perm = get_order(seriate(D, method = "OLO"))


	CC = cor(W)
	CC[is.na(CC)] = 0
	D = as.dist(1-CC)
	col.perm = get_order(seriate(D, method = "OLO"))

	gradPal = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdYlBu")))(100)

	Heatmap(W[row.perm, col.perm], name = measure_name, cluster_rows = F, cluster_columns = F, col = gradPal, row_title = row_title, column_title = col_title, column_names_gp = gpar(fontsize = 8, fontface="bold"), row_names_gp = gpar(fontsize = 8, fontface="bold"), column_title_gp = gpar(fontsize = 10, fontface="bold"), row_title_gp = gpar(fontsize = 10, fontface="bold"))
}
