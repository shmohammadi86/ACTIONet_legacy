mark.doublet.cells <- function(sce) {
    doublet.scores = scran::doubletCells(sce)
    names(doublet.scores) = sce$cell.hashtag
    
    return(doublet.scores)
}

mark.doublet.clusters <- function(sce, clusters) {
    dbl = scran::doubletCluster(sce, clusters)
    doublet.score = dbl$N
    
    is.doublet = scater::isOutlier(doublet.score, log = TRUE, type = "lower")
    
    res = list(doublet.score = doublet.score, is.doublet = is.doublet)
    return(res)
}


highlight.ACTIONet <- function(ACTIONet.out, Labels = NULL, z.threshold = -1) {
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

assess.batch.mixing <- function(ACTIONet.out, sce, batch.vec, thread_no = 8) {
    require(kBET)
    
    knn.size = round(1/4 * mean(table(batch.vec)))
    KNN = ACTIONet::computeNearestDist_edgeList(ACTIONet.out$reconstruct.out$H_stacked, knn.size, thread_no)
    KNN$nn.index = KNN$idx[, -1]
    
    data = t(logcounts(sce))
    batch.estimate <- kBET::kBET(data, batch.vec, knn = KNN, k0 = knn.size, plot = FALSE)
    
    # plot.data <- data.frame(class=rep(c('observed', 'expected'), each=length(batch.estimate$stats$kBET.observed)), data =
    # c(batch.estimate$stats$kBET.observed, batch.estimate$stats$kBET.expected)) g <- ggplot(plot.data, aes(class, data)) +
    # geom_boxplot() + labs(x='Test', y='Rejection rate',title='kBET test results') + theme_bw() + scale_y_continuous(limits=c(0,1))
    # plot(g)
    
    return(batch.estimate)
}
