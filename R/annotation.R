annotate.archetypes.using.labels <- function(ACTIONet.out, Labels, rand_perm_no = 1000, core = F) {
    if (!is.factor(Labels)) {
        Labels = factor(Labels)
    }
    
    if(core == T) {
		profile = ACTIONet.out$unification.out$H.core		
	} else {
		profile = ACTIONet.out$reconstruct.out$H_stacked
	}
    
    Enrichment.Z = sapply(levels(Labels), function(label) {
        mask = Labels == label
        class.profile = profile[, mask]
        null.profile = profile[, !mask]
        
        N.class = sum(mask)
        N.null = sum(!mask)
        
        mu.class = Matrix::rowMeans(class.profile)
        mu.null = Matrix::rowMeans(null.profile)
        
        sigma_sq.class = apply(class.profile, 1, var)
        sigma_sq.null = apply(null.profile, 1, var)
        
        
        delta.mean = mu.class - mu.null
        t.stat = delta.mean/sqrt((sigma_sq.class/N.class) + (sigma_sq.null/N.null))
        return(t.stat)
    })
    
    archetypeLabels = levels(Labels)[apply(Enrichment.Z, 1, which.max)]
    archetypeLabels = factor(archetypeLabels, levels = levels(Labels))
    
    out.list = list(Labels = archetypeLabels, Enrichment = Enrichment.Z)
    
    return(out.list)
}

annotate.clusters.using.labels <- function(clusters, Labels) {
    if (!is.factor(Labels)) {
        Labels = factor(Labels, levels = sort(unique(Labels)))
    }
    clusters = as.numeric(clusters)
    
    pop.size = length(Labels)
    pos.size = table(Labels)
    
    logPvals = sapply(sort(unique(clusters)), function(i) {
        idx = which(clusters == i)
        sample.size = length(idx)
        success.size = table(Labels[idx])
        
        logPval = HGT_tail(pop.size, pos.size, sample.size, success.size)
        
        return(logPval)
    })
    
    clusterLabels = factor(levels(Labels)[apply(logPvals, 2, which.max)])
    cellLabels = factor(clusterLabels[clusters], levels = levels(clusterLabels))
    
    fullLabels = sapply(sort(unique(clusters)), function(i) {
        return(sprintf("Cluster %d (%s)", i, clusterLabels[[i]]))
    })
    fullLabels = factor(fullLabels, levels = fullLabels)
    
    cellFullLabels = factor(fullLabels[clusters], levels = levels(fullLabels))
    
    return(list(Labels = clusterLabels, cellLabels = cellLabels, fullLabels = fullLabels, cellFullLabels = cellFullLabels, Enrichment = logPvals))
}

annotate.archetypes.using.markers <- function(ACTIONet.out, marker.genes, rand.sample.no = 1000, core = F) {
    require(ACTIONet)
    require(igraph)
    require(Matrix)
    require(stringr)
    

	if(core == T) {
		if (("unification.out" %in% names(ACTIONet.out))) {
			print("Using unification.out$DE.core (merged archetypes)")
			archetype.panel = as.matrix(log1p(t(ACTIONet.out$unification.out$DE.core)))
		} else {
			print("unification.out is not in ACTIONet.out. Please run unify.cell.states() first.")
			return()
		}
	} else {
		if (("archetype.differential.signature" %in% names(ACTIONet.out))) {
			print("Using archetype.differential.signature (all archetypes)")
			archetype.panel = as.matrix(log1p(t(ACTIONet.out$archetype.differential.signature)))
		} else {
			print("archetype.differential.signature is not in ACTIONet.out. Please run compute.archetype.gene.specificity() first.")
			return()
		}
	}      

	    
    GS.names = names(marker.genes)
    if (is.null(GS.names)) {
        GS.names = sapply(1:length(GS), function(i) sprintf("Celltype %s", i))
    }
    
    markers.table = do.call(rbind, lapply(names(marker.genes), function(celltype) {
        genes = marker.genes[[celltype]]
        if (length(genes) == 0) 
            return(data.frame())
        
        
        signed.count = sum(sapply(genes, function(gene) grepl("\\+$|-$", gene)))
        is.signed = signed.count > 0
        
        if (!is.signed) {
            df = data.frame(Gene = toupper(genes), Direction = +1, Celltype = celltype)
        } else {
            # pos.genes = toupper(as.character(sapply(genes[sapply(genes, function(gene) grepl('\\+$', gene))]))) neg.genes =
            # toupper(as.character(sapply(genes[sapply(genes, function(gene) grepl('-$', gene))])))
            
            pos.genes = toupper(as.character(sapply(genes[grepl("+", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene, 
                stringr::fixed("+"), ""))))
            neg.genes = toupper(as.character(sapply(genes[grepl("-", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene, 
                stringr::fixed("-"), ""))))
            
            df = data.frame(Gene = c(pos.genes, neg.genes), Direction = c(rep(+1, length(pos.genes)), rep(-1, length(neg.genes))), 
                Celltype = celltype)
        }
    }))
    markers.table = markers.table[markers.table$Gene %in% colnames(archetype.panel), ]
    
    if (dim(markers.table)[1] == 0) {
        print("No markers are left")
        return()
    }    
               
    IDX = split(1:dim(markers.table)[1], markers.table$Celltype)
    
    print("Computing significance scores")
    set.seed(0)
    Z = sapply(IDX, function(idx) {
        markers = toupper(as.character(markers.table$Gene[idx]))
        directions = markers.table$Direction[idx]
        mask = markers %in% colnames(archetype.panel)
        
        A = as.matrix(archetype.panel[, markers[mask]])
        sgn = as.numeric(directions[mask])
        stat = A %*% sgn
        
        rand.stats = sapply(1:rand.sample.no, function(i) {
            rand.samples = sample.int(dim(archetype.panel)[2], sum(mask))
            rand.A = as.matrix(archetype.panel[, rand.samples])
            rand.stat = rand.A %*% sgn
        })
        
        cell.zscores = as.numeric((stat - apply(rand.stats, 1, mean))/apply(rand.stats, 1, sd))
        
        return(cell.zscores)
    })
    
    Z[is.na(Z)] = 0
    Labels = colnames(Z)[apply(Z, 1, which.max)]
    
    L = names(marker.genes)
    L = L[L %in% Labels]
    Labels = factor(Labels, levels = L)
    Labels.conf = apply(Z, 1, max)
    
    names(Labels) = rownames(archetype.panel)
    names(Labels.conf) = rownames(archetype.panel)
    rownames(Z) = rownames(archetype.panel)
    
    out.list = list(Labels = Labels, Labels.confidence = Labels.conf, Enrichment = Z, archetype.panel = archetype.panel)
    
    return(out.list)
}



annotate.archetypes.using.markers.fromMatrix <- function(ACTIONet.out, marker.mat, rand.sample.no = 1000, core = F) {
    require(ACTIONet)
    require(igraph)
    require(Matrix)
    require(stringr)
        
	if(core == T) {
		if (("unification.out" %in% names(ACTIONet.out))) {
			print("Using unification.out$DE.core (merged archetypes)")
			archetype.panel = as.matrix(log1p(t(ACTIONet.out$unification.out$DE.core)))
		} else {
			print("unification.out is not in ACTIONet.out. Please run unify.cell.states() first.")
			return()
		}
	} else {
		if (("archetype.differential.signature" %in% names(ACTIONet.out))) {
			print("Using archetype.differential.signature (all archetypes)")
			archetype.panel = as.matrix(log1p(t(ACTIONet.out$archetype.differential.signature)))
		} else {
			print("archetype.differential.signature is not in ACTIONet.out. Please run compute.archetype.gene.specificity() first.")
			return()
		}
	}        
	   
    #archetype.panel = archetype.panel[, common.genes]
    common.genes = intersect(colnames(archetype.panel), rownames(marker.mat))
    marker.mat = marker.mat[common.genes, ]		

    if (dim(marker.mat)[1] == 0) {
        print("No markers are left")
        return()
    }    
    
    
    print("Computing significance scores")
    set.seed(0)
    Z = apply(marker.mat, 2, function(x) {
		mask = (x != 0)
		if(sum(mask) < 3) {
			return(rep(0, nrow(archetype.panel)))
		}
		markers = rownames(marker.mat)[mask]
        
        A = as.matrix(archetype.panel[, markers])
        sgn = sign(x[mask])
        stat = A %*% sgn
        
        rand.stats = sapply(1:rand.sample.no, function(i) {
            rand.samples = sample.int(dim(archetype.panel)[2], sum(mask))
            rand.A = as.matrix(archetype.panel[, rand.samples])
            rand.stat = rand.A %*% sgn
        })
        
        cell.zscores = as.numeric((stat - apply(rand.stats, 1, mean))/apply(rand.stats, 1, sd))
        
        return(cell.zscores)
    })
    
    Z[is.na(Z)] = 0
    Labels = colnames(Z)[apply(Z, 1, which.max)]
    
    L = colnames(marker.mat)
    L = L[L %in% Labels]
    Labels = factor(Labels, levels = L)
    Labels.conf = apply(Z, 1, max)
    
    names(Labels) = rownames(archetype.panel)
    names(Labels.conf) = rownames(archetype.panel)
    rownames(Z) = rownames(archetype.panel)
    
    out.list = list(Labels = Labels, Labels.confidence = Labels.conf, Enrichment = Z, archetype.panel = archetype.panel)
    
    return(out.list)
}



annotate.core.archetypes.using.markers <- function(ACTIONet.out, sce, marker.genes, rand.sample.no = 1000) {
    require(ACTIONet)
    require(igraph)
    require(Matrix)
    require(stringr)
    
    GS.names = names(marker.genes)
    if (is.null(GS.names)) {
        GS.names = sapply(1:length(GS), function(i) sprintf("Celltype %s", i))
    }
    
    markers.table = do.call(rbind, lapply(names(marker.genes), function(celltype) {
        genes = marker.genes[[celltype]]
        if (length(genes) == 0) 
            return(data.frame())
        
        
        signed.count = sum(sapply(genes, function(gene) grepl("\\+$|-$", gene)))
        is.signed = signed.count > 0
        
        if (!is.signed) {
            df = data.frame(Gene = toupper(genes), Direction = +1, Celltype = celltype)
        } else {
            # pos.genes = toupper(as.character(sapply(genes[sapply(genes, function(gene) grepl('\\+$', gene))]))) neg.genes =
            # toupper(as.character(sapply(genes[sapply(genes, function(gene) grepl('-$', gene))])))
            
            pos.genes = toupper(as.character(sapply(genes[grepl("+", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene, 
                stringr::fixed("+"), ""))))
            neg.genes = toupper(as.character(sapply(genes[grepl("-", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene, 
                stringr::fixed("-"), ""))))
            
            df = data.frame(Gene = c(pos.genes, neg.genes), Direction = c(rep(+1, length(pos.genes)), rep(-1, length(neg.genes))), 
                Celltype = celltype)
        }
    }))
    markers.table = markers.table[markers.table$Gene %in% colnames(archetype.panel), ]
    if (dim(markers.table)[1] == 0) {
        print("No markers are left")
        return()
    }
    
    
    print("Computing signifcance of genes in archetypes")
    if (!("unification.out" %in% names(ACTIONet.out))) {
        print("Computing unification.out")
        ACTIONet.out$unification.out = unify.cell.states(ACTIONet.out, sce)
    }
    archetype.panel = t(ACTIONet.out$unification.out$signature.core)
    colnames(archetype.panel) = rownames(sce)
    
    
    IDX = split(1:dim(markers.table)[1], markers.table$Celltype)
    
    print("Computing significance scores")
    set.seed(0)
    Z = sapply(IDX, function(idx) {
        markers = toupper(as.character(markers.table$Gene[idx]))
        directions = markers.table$Direction[idx]
        mask = markers %in% colnames(archetype.panel)
        
        A = as.matrix(archetype.panel[, markers[mask]])
        sgn = as.numeric(directions[mask])
        stat = A %*% sgn
        
        rand.stats = sapply(1:rand.sample.no, function(i) {
            rand.samples = sample.int(dim(archetype.panel)[2], sum(mask))
            rand.A = as.matrix(archetype.panel[, rand.samples])
            rand.stat = rand.A %*% sgn
        })
        
        cell.zscores = as.numeric((stat - apply(rand.stats, 1, mean))/apply(rand.stats, 1, sd))
        
        return(cell.zscores)
    })
    
    Z[is.na(Z)] = 0
    Labels = colnames(Z)[apply(Z, 1, which.max)]
    
    L = names(marker.genes)
    L = L[L %in% Labels]
    Labels = factor(Labels, levels = L)
    Labels.conf = apply(Z, 1, max)
    
    names(Labels) = rownames(archetype.panel)
    names(Labels.conf) = rownames(archetype.panel)
    rownames(Z) = rownames(archetype.panel)
    
    out.list = list(Labels = Labels, Labels.confidence = Labels.conf, Enrichment = Z, archetype.panel = archetype.panel)
    
    return(out.list)
}

annotate.clusters.using.markers <- function(sce, clusters, marker.genes, rand.sample.no = 1000) {
    require(ACTIONet)
    require(igraph)
    require(Matrix)
    require(stringr)
    
    GS.names = names(marker.genes)
    if (is.null(GS.names)) {
        GS.names = sapply(1:length(GS), function(i) sprintf("Celltype %s", i))
    }
    
    markers.table = do.call(rbind, lapply(names(marker.genes), function(celltype) {
        genes = marker.genes[[celltype]]
        if (length(genes) == 0) 
            return(data.frame())
        
        
        signed.count = sum(sapply(genes, function(gene) grepl("\\+$|-$", gene)))
        is.signed = signed.count > 0
        
        if (!is.signed) {
            df = data.frame(Gene = toupper(genes), Direction = +1, Celltype = celltype)
        } else {
            # pos.genes = toupper(as.character(sapply(genes[sapply(genes, function(gene) grepl('\\+$', gene))]))) neg.genes =
            # toupper(as.character(sapply(genes[sapply(genes, function(gene) grepl('-$', gene))])))
            
            pos.genes = toupper(as.character(sapply(genes[grepl("+", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene, 
                stringr::fixed("+"), ""))))
            neg.genes = toupper(as.character(sapply(genes[grepl("-", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene, 
                stringr::fixed("-"), ""))))
            
            df = data.frame(Gene = c(pos.genes, neg.genes), Direction = c(rep(+1, length(pos.genes)), rep(-1, length(neg.genes))), 
                Celltype = celltype)
        }
    }))
    markers.table = markers.table[markers.table$Gene %in% rownames(sce), ]
    if (dim(markers.table)[1] == 0) {
        print("No markers are left")
        return()
    }
    
    print("Computing signifcance of genes in cluster")
    X = compute.annotations.gene.specificity(sce, clusters)
    if (is.null(X)) {
        print("Cannot compute cluster DE. Returning")
        return()
    }
    expression.panel = Matrix::t(X)
    
    colnames(expression.panel) = rownames(sce)
    
    
    IDX = split(1:dim(markers.table)[1], markers.table$Celltype)
    
    print("Computing significance scores")
    set.seed(0)
    Z = sapply(IDX, function(idx) {
        markers = toupper(as.character(markers.table$Gene[idx]))
        directions = markers.table$Direction[idx]
        mask = markers %in% colnames(expression.panel)
        
        A = as.matrix(expression.panel[, markers[mask]])
        sgn = as.numeric(directions[mask])
        stat = A %*% sgn
        
        rand.stats = sapply(1:rand.sample.no, function(i) {
            rand.samples = sample.int(dim(expression.panel)[2], sum(mask))
            rand.A = as.matrix(expression.panel[, rand.samples])
            rand.stat = rand.A %*% sgn
        })
        
        cell.zscores = as.numeric((stat - apply(rand.stats, 1, mean))/apply(rand.stats, 1, sd))
        
        return(cell.zscores)
    })
    
    clusterLabels = factor(names(marker.genes)[apply(Z, 1, which.max)], levels = names(marker.genes))
    cellLabels = factor(clusterLabels[clusters], levels = levels(clusterLabels))
    
    fullLabels = sapply(sort(unique(clusters)), function(i) {
        return(sprintf("Cluster %d (%s)", i, clusterLabels[[i]]))
    })
    fullLabels = factor(fullLabels, levels = fullLabels)
    
    cellFullLabels = factor(fullLabels[clusters], levels = levels(fullLabels))
    
    res = list(Labels = clusterLabels, cellLabels = cellLabels, fullLabels = fullLabels, cellFullLabels = cellFullLabels, Enrichment = Z, 
        Cluster.DE = X)
    return(res)
}




annotate.cells.using.markers <- function(ACTIONet.out, sce, marker.genes, alpha_val = 0.9, rand.sample.no = 100, thread_no = 8, imputation = "PageRank") {
    require(ACTIONet)
    require(igraph)
    require(Matrix)
    require(stringr)
    
    rownames(sce) = toupper(rownames(sce))
    
    GS.names = names(marker.genes)
    if (is.null(GS.names)) {
        GS.names = sapply(1:length(GS), function(i) sprintf("Celltype %s", i))
    }
    
    markers.table = do.call(rbind, lapply(names(marker.genes), function(celltype) {
        genes = marker.genes[[celltype]]
        if (length(genes) == 0) 
            return(data.frame())
        
        
        signed.count = sum(sapply(genes, function(gene) grepl("\\+$|-$", gene)))
        is.signed = signed.count > 0
        
        if (!is.signed) {
            df = data.frame(Gene = toupper(genes), Direction = +1, Celltype = celltype)
        } else {
            # pos.genes = toupper(as.character(sapply(genes[sapply(genes, function(gene) grepl('\\+$', gene))]))) neg.genes =
            # toupper(as.character(sapply(genes[sapply(genes, function(gene) grepl('-$', gene))])))
            
            pos.genes = toupper(as.character(sapply(genes[grepl("+", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene, 
                stringr::fixed("+"), ""))))
            neg.genes = toupper(as.character(sapply(genes[grepl("-", genes, fixed = TRUE)], function(gene) stringr::str_replace(gene, 
                stringr::fixed("-"), ""))))
            
            df = data.frame(Gene = c(pos.genes, neg.genes), Direction = c(rep(+1, length(pos.genes)), rep(-1, length(neg.genes))), 
                Celltype = celltype)
        }
    }))
    markers.table = markers.table[markers.table$Gene %in% rownames(sce), ]
    if (dim(markers.table)[1] == 0) {
        print("No markers are left")
        return()
    }
    
    rows = match(markers.table$Gene, rownames(sce))
    if (imputation == "PageRank") {
        # PageRank-based imputation
        print("Using PageRank for imptation of marker genes")
        imputed.marker.expression = impute.genes.using.ACTIONet(ACTIONet.out, sce, markers.table$Gene, alpha_val, thread_no, prune = TRUE, 
            rescale = TRUE)
    } else {
        # PCA-based imputation
        print("Using archImpute for imptation of marker genes")
        imputed.marker.expression = t(ACTIONet.out$signature.profile[rows, ACTIONet.out$core.out$core.archs] %*% ACTIONet.out$core.out$H)
    }
    colnames(imputed.marker.expression) = toupper(colnames(imputed.marker.expression))
    
    
    IDX = split(1:dim(markers.table)[1], markers.table$Celltype)
    
    print("Computing significance scores")
    set.seed(0)
    Z = sapply(IDX, function(idx) {
        markers = toupper(as.character(markers.table$Gene[idx]))
        directions = markers.table$Direction[idx]
        mask = markers %in% colnames(imputed.marker.expression)
        
        A = as.matrix(imputed.marker.expression[, markers[mask]])
        sgn = as.numeric(directions[mask])
        stat = A %*% sgn
        
        rand.stats = sapply(1:rand.sample.no, function(i) {
            rand.samples = sample.int(dim(imputed.marker.expression)[2], sum(mask))
            rand.A = as.matrix(imputed.marker.expression[, rand.samples])
            rand.stat = rand.A %*% sgn
        })
        
        cell.zscores = as.numeric((stat - apply(rand.stats, 1, mean))/apply(rand.stats, 1, sd))
        
        return(cell.zscores)
    })
    
    Z[is.na(Z)] = 0
    Labels = colnames(Z)[apply(Z, 1, which.max)]
    
    L = names(marker.genes)
    L = L[L %in% Labels]
    Labels = factor(Labels, levels = L)
    Labels.conf = apply(Z, 1, max)
    
    
    names(Labels) = ACTIONet.out$log$cells   
	
	if(! ('annotations' %in% names(ACTIONet.out)) ) {
		ACTIONet.out$annotations = list()
	}

	time.stamp = as.character(Sys.time())
	if(is.null(annotation.name)) {
		annotation.name = sprintf('%s', time.stamp)
	}
	h = hashid_settings(salt = time.stamp, min_length = 8)
	annotation.hashtag = ENCODE(length(ACTIONet.out$annotations)+1, h)
	
	res = list(Labels = Labels, Labels.confidence = Labels.conf, DE.profile = NULL, highlight = NULL, cells = ACTIONet.out$log$cells, time.stamp = time.stamp, annotation.name = annotation.hashtag, type = "annotate.cells.using.markers", Enrichment = Z)
	
	cmd = sprintf("ACTIONet.out$annotations$\"%s\" = res", annotation.name)	
	eval(parse(text=cmd))
	

	if( !ignore.DE ) {
		cmd = sprintf("ACTIONet.out = compute.annotations.gene.specificity(ACTIONet.out, sce.red, \"%s\")", annotation.name)	
		eval(parse(text=cmd))
	}	

    return(ACTIONet.out)    
}


map.cell.scores.from.archetype.enrichment <- function(ACTIONet.out, Enrichment, core = T) {
    if( (core == T) | (nrow(Enrichment) == nrow(ACTIONet.out$unification.out$H.core)) ) {
		cell.scores.mat = ACTIONet.out$unification.out$H.core		
	} else {
		cell.scores.mat = ACTIONet.out$reconstruct.out$H_stacked
	}    
    cs = Matrix::colSums(cell.scores.mat)
    cs[cs == 0] = 1
    cell.scores.mat = t(scale(cell.scores.mat, center = F, scale = cs))
    

    if (ncol(Enrichment) == ncol(cell.scores.mat)) {
		print("Flipping enrichment matrix")
        Enrichment = t(Enrichment)
    }


	Z = (Enrichment - mean(Enrichment)) / sd(Enrichment)
	Enrichment = 1 / (1 + exp(-Z))
	Enrichment[is.na(Enrichment)] = 0
	
	rs = sqrt(Matrix::rowSums(Enrichment))
	rs[rs == 0] = 1
	D_r = Matrix::Diagonal(nrow(Enrichment), 1/rs)
	
	cs = sqrt(Matrix::colSums(Enrichment))
	cs[cs == 0] = 1
	D_c = Matrix::Diagonal(ncol(Enrichment), 1/cs)
	
	Enrichment.scaled = as.matrix(D_r %*% Enrichment %*% D_c)

    
    cell.Enrichment.mat = cell.scores.mat %*% Enrichment.scaled
    colnames(cell.Enrichment.mat) = colnames(Enrichment)
    rownames(cell.Enrichment.mat) = ACTIONet.out$log$cells
    
    return(cell.Enrichment.mat)
}
	
annotate.cells.from.archetype.enrichment <- function(ACTIONet.out, Enrichment, core = T, annotation.name = NULL, ignore.DE = FALSE) {
    cell.Enrichment.mat = map.cell.scores.from.archetype.enrichment(ACTIONet.out, Enrichment)
    
    cell.Labels = colnames(cell.Enrichment.mat)[apply(cell.Enrichment.mat, 1, which.max)]
    cell.Labels.factor = factor(cell.Labels)
    cell.Labels.conf = apply(cell.Enrichment.mat, 1, max)
    

    names(cell.Labels.factor) = ACTIONet.out$log$cells   
	
	if(! ('annotations' %in% names(ACTIONet.out)) ) {
		ACTIONet.out$annotations = list()
	}

	time.stamp = as.character(Sys.time())
	if(is.null(annotation.name)) {
		annotation.name = sprintf('%s', time.stamp)
	}
	h = hashid_settings(salt = time.stamp, min_length = 8)
	annotation.hashtag = ENCODE(length(ACTIONet.out$annotations)+1, h)
	
	res = list(Labels = cell.Labels.factor, Labels.confidence = cell.Labels.conf, DE.profile = NULL, highlight = NULL, cells = ACTIONet.out$log$cells, time.stamp = time.stamp, annotation.name = annotation.hashtag, type = "annotate.cells.from.archetype.enrichment", Enrichment = cell.Enrichment.mat)
	
	cmd = sprintf("ACTIONet.out$annotations$\"%s\" = res", annotation.name)	
	eval(parse(text=cmd))
	

	if( !ignore.DE ) {
		cmd = sprintf("ACTIONet.out = compute.annotations.gene.specificity(ACTIONet.out, sce.red, \"%s\")", annotation.name)	
		eval(parse(text=cmd))
	}	

    return(ACTIONet.out)
}

annotate.cells.from.archetypes.using.markers <- function(ACTIONet.out, marker.genes, annotation.name = NULL, rand.sample.no = 1000, ignore.DE = FALSE) {
    arch.annot = annotate.archetypes.using.markers(ACTIONet.out, marker.genes, rand.sample.no = rand.sample.no, core = T)
    
    Enrichment = arch.annot$Enrichment
    ACTIONet.out = annotate.cells.from.archetype.enrichment(ACTIONet.out, Enrichment, core = T, annotation.name = annotation.name, ignore.DE = ignore.DE)

    return(ACTIONet.out)
}


prioritize.celltypes <- function(ACTIONet.out, species = "Human", min.score = 3, plot = T) {
	if(tolower(species) == "human") {
		data("CellMarkerDB_human")
		CellMarker.annot.out = annotate.archetypes.using.markers.fromMatrix(ACTIONet.out, CellMarkerDB_human, core = T)
	} else if(tolower(species) == "mouse"){
		data("CellMarkerDB_mouse")
		CellMarker.annot.out = annotate.archetypes.using.markers.fromMatrix(ACTIONet.out, CellMarkerDB_mouse, core = T)
	} else {
		R.utils::printf("Unknown species: %s\n", species)
	}
	CellMarker.Enrichment = t(scale(t(CellMarker.annot.out$Enrichment)))
	
	CC = cor(CellMarker.Enrichment)
	diag(CC) = 0
	
#	max.vals = apply(scale(t(CellMarker.Enrichment)), 1 , max)
	max.vals = apply(CellMarker.Enrichment, 2 , max)
	
	CC.clusters = signed_cluster(as(CC, 'sparseMatrix'), seed = 0)
	IDX = split(1:length(CC.clusters), CC.clusters)
	scores = sapply(IDX, function(idx) {median(max.vals[idx])})
	scores[scores < 0] = 0
	
	nnz = round( (sum(scores^2)^2) / (sum(scores^4)) )
	
	perm = order(scores, decreasing = T)[1:nnz]
	ordered.labels = sapply(IDX[perm], function(idx) {  colnames(CellMarker.Enrichment)[idx[which.max(max.vals[idx])]]  })
	
	df = data.frame(Celltype = ordered.labels, Score = scores[perm])
	
	if(plot == T) {
		require(RColorBrewer)
		require(ComplexHeatmap)
		subEnrichment = t(scale(CellMarker.Enrichment[, as.character(df$Celltype)]))
		gradPal = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
		ht = Heatmap(subEnrichment, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8), name = "Enrichment", row_title = "Putative Cell-types", column_title = "Cell States", col=gradPal)	
		show(ht)
	}
	
	res = list(ranked.celltypes = df, Enrichment = CellMarker.Enrichment)
	return(res)
}

highlight.annotations <- function(ACTIONet.out, annotation.tag = NULL, z.threshold = -1) {
    if (is.null(Labels)) {
        clusters = cluster.ACTIONet.highRes(ACTIONet.out)
    } else {
        clusters = as.numeric(factor(Labels))
    }
    
    IDX = split(1:length(clusters), clusters)
    cluster.cell.connectivity = vector("list", length(IDX))
    cluster.cell.connectivity.smoothed = vector("list", length(IDX))
    cluster.pruned.cells = vector("list", length(IDX))
    
    for (i in 1:length(IDX)) {
        idx = IDX[[i]]
        
        sub.ACTIONet = igraph::induced.subgraph(ACTIONet.out$ACTIONet, V(ACTIONet.out$ACTIONet)[idx])
        sub.cn = coreness(sub.ACTIONet)
        
        if (mad(sub.cn) > 0) {
            z = (sub.cn - median(sub.cn))/mad(sub.cn)
        } else if (sd(sub.cn) > 0) {
            z = as.numeric(scale(sub.cn))
        } else {
            z = (as.numeric(rep(0, length(idx))))
        }
        cluster.cell.connectivity[[i]] = z
        
        cluster.pruned.cells[[i]] = idx[z < z.threshold]
    }
    all.cell.connectivity.scores = as.numeric(sparseVector(unlist(cluster.cell.connectivity), unlist(IDX), length(clusters)))
    all.pruned.cells = sort(unique(unlist(cluster.pruned.cells)))
    
    if (length(all.pruned.cells) > 0) {
        is.pruned = as.numeric(sparseVector(1, all.pruned.cells, length(clusters)))
    } else {
        is.pruned = rep(0, ncol(sce))
    }
    
    
    out = list(connectivity.scores = all.cell.connectivity.scores, pruned.cells = all.pruned.cells, cluster.connectivity.scores = cluster.cell.connectivity, 
        cluster.pruned.cells = cluster.pruned.cells, is.pruned = is.pruned)
    
    return(out)
}

