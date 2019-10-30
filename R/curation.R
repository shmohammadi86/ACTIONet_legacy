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
