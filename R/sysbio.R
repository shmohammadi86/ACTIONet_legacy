run.archetype.SCINET <- function(ACTIONet.out, thread_no = 8, core = T) {
	require(SCINET)
	
	if( !exists('PCNet') ) {
		data("PCNet")
	}	
       
	if(core == T) {
		if (("unification.out" %in% names(ACTIONet.out))) {
			print("Using unification.out$DE.core (merged archetypes)")
			DE.profile = as.matrix(log1p(ACTIONet.out$unification.out$DE.core@assays[["significance"]]))
		} else {
			print("unification.out is not in ACTIONet.out. Running unify.cell.states() first ...")
			ACTIONet.out$unification.out = unify.cell.states(ACTIONet.out, sce, reduction_slot = reduction_slot, sce.data.attr = sce.data.attr)    
			DE.profile = as.matrix(log1p(ACTIONet.out$unification.out$DE.core@assays[["significance"]]))
		}
	} else {
		if (("archetype.differential.signature" %in% names(ACTIONet.out))) {
			print("Using archetype.differential.signature (all archetypes)")
			DE.profile = as.matrix(log1p(ACTIONet.out$archetype.differential.signature))
		} else {
			print("archetype.differential.signature is not in ACTIONet.out. Running compute.archetype.feature.specificity() first ...")
			ACTIONet.out$archetype.differential.signature = compute.archetype.feature.specificity(ACTIONet.out, sce, mode = specificity.mode, sce.data.attr = sce.data.attr)
			DE.profile = as.matrix(log1p(ACTIONet.out$archetype.differential.signature))
		}
	}                  
	
	
	common.genes = intersect(rownames(DE.profile), rownames(PCNet))
	
	A = DE.profile[common.genes, ]
	G = PCNet[common.genes, common.genes]
	
	gene.activity.scores = SCINET::compute_gene_activities_full(A = A, thread_no = thread_no)
	cellstate.nets = SCINET::construct_cell_networks(net = G, gene_activities = gene.activity.scores, thread_no = thread_no)
	cellstate.nets.list = as.list(cellstate.nets[, 1])
	
	if(core == T) {
		ACTIONet.out$unification.out$SCINET.nets = cellstate.nets.list
	} else {
		ACTIONet.out$SCINET.nets = cellstate.nets
	}                  
	
	return(ACTIONet.out)
}

run.cluster.SCINET <- function() {
}
