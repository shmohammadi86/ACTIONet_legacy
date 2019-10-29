#include <actionetcore.h>

namespace ACTIONetcore {
	multilevel_archetypal_decomposition reconstructArchetypes(sp_mat S, vector<mat> C_trace, vector<mat> H_trace, double z_threshold = -1.0) {	
		mat C_stacked;
		mat H_stacked;
		int sample_no = S.n_cols;
		int depth = H_trace.size();
		
		multilevel_archetypal_decomposition results;
		
		// Vector contains an element for k==0, this have to -1
		printf("Joining the trace of C& H matrices (depth = %d) ... ", depth-1); fflush(stdout);
		// Group H and C matrices for different values of k (#archs) into joint matrix
		for(int k = 0; k < depth; k++) {						
			if(H_trace[k].n_rows == 0)
				continue;
			
			if(H_stacked.n_elem == 0) {
				C_stacked = C_trace[k];			
				H_stacked = H_trace[k];			
				
			}
			else {
				C_stacked = join_rows(C_stacked, C_trace[k]);
				H_stacked = join_cols(H_stacked, H_trace[k]);
			}
		}	
		int total_archs = H_stacked.n_rows;
		printf("done\n");



		// Construct backbone
		printf("Constructing backbone graph ... ");
		
		/*
		// Compute distances between cells associated with each archetype (driving cells?)
		uvec all_cells;
		vector <uvec> cells(total_archs);
		for(int i = 0; i < total_archs; i++) {
			cells[i] = find(C_stacked.col(i) > 0);			
			all_cells = join_cols(all_cells, cells[i]);
		}
		all_cells = unique(all_cells);
		
		mat subH = H_stacked.cols(all_cells);		
		mat subD = computeFullDist(subH, 1, 0);	
		
		double cell_pair_no = all_cells.n_elem*(all_cells.n_elem-1)/2;
		double mean_null_d = sum(sum(trimatu(subD))) / cell_pair_no;
		printf("\tNull d = %.2f\n", mean_null_d);
		
		ivec CC;
		uvec idx_i, idx_j, tmp;		
	
		mat backbone = zeros(total_archs, total_archs);
		for(int i = 0; i < total_archs; i++) {
			for(int j = i+1; j < total_archs; j++) {
				uvec cells_i = cells[i];
				uvec cells_j = cells[j];				
				
				intersect(CC, idx_i, tmp, conv_to<ivec>::from(all_cells), conv_to<ivec>::from(cells_i));						
				intersect(CC, idx_j, tmp, conv_to<ivec>::from(all_cells), conv_to<ivec>::from(cells_j));						

				vec d_inter = vectorise(subD(idx_i, idx_j));
				double mean_d = mean(d_inter);
				
				if(mean_d < mean_null_d) 
					backbone(i, j) = 1.0 - mean_d;
			}
		}
		*/

		
		mat backbone = cor(trans(H_stacked));
		backbone.diag().zeros();
		backbone.transform( [](double val) { return (val< 0?0:val); } );
		
		printf("done\n");


		// Barrat weighted transitivity: formulation from "Clustering Coefficients for Weighted Networks" (Kalna)
		printf("Pruning non-specific archetypes ...");
		vec pruned = zeros(total_archs);

		vec transitivity = zeros(total_archs);
		vec s = sum(backbone, 1); // strength of nodes
		vec d = vec(sum(spones(sp_mat(backbone)), 1));
		for(int k = 0; k < total_archs; k++) {
			
			double sum = 0;
			for(int i = 0; i < total_archs; i++) {
				double w_ki = backbone(k, i);
				for(int j = 0; j < total_archs; j++) {
					double w_kj = backbone(k, j);
					
					double mean_weight = (w_ki + w_kj) / 2.0;					
					double triangle_mask = backbone(i, k)*backbone(k, j)*backbone(j, i) > 0?1:0;
					
					sum += mean_weight*triangle_mask;
				}
			}
			transitivity(k) = sum / (s(k)*(d(k)-1));
		}
			
		vec transitivity_z = zscore(transitivity);
		uvec nonspecific_idx = find(transitivity_z < z_threshold);
		pruned(nonspecific_idx).ones(); 
		printf("done (%d archs pruned)\n", nonspecific_idx.n_elem);
		
		// Find landmark cells, i.e., closest cells to each multi-level archetype (its projection on to the cell space, ish)
		printf("Identifying landmark cells ... "); fflush(stdout);
		double epsilon = 1e-6;
		int bad_archs = 0;
		vec landmark_cells = -ones(total_archs);
		for(int i = 0; i < total_archs; i++) {
			vec h = trans(H_stacked.row(i));
			vec c = C_stacked.col(i);
			
			uvec h_landmarks = find( (max(h) - h) < epsilon );
			uvec c_landmarks = find( 0 < c );
			uvec common_landmarks = intersect(h_landmarks, c_landmarks);
			
			if(0 < common_landmarks.n_elem) { // They don't agree on any samples!
				landmark_cells(i) = common_landmarks(index_max(c(common_landmarks)));
			}
			else { // Potentially noisy archetype
				pruned(i) = 1;
				bad_archs++;
			}			
		}				
		printf("done (further %d archs removed)\n", bad_archs); fflush(stdout);
		
		uvec selected_archs = find(pruned == 0);		
		results.backbone = backbone(selected_archs, selected_archs);
		results.selected_archs = selected_archs;
		results.landmark_cells = landmark_cells(selected_archs);
		results.C_stacked = C_stacked.cols(selected_archs);		
		results.H_stacked = H_stacked.rows(selected_archs);		
		
		
		// Original, raw archetype profile (gene x total archs)
		printf("Reconstructing archetypes ... "); fflush(stdout);
		results.archetype_profile = mat(S*sp_mat(results.C_stacked));
		printf("done\n"); fflush(stdout);
		
		return(results);
	}
	
	
	sp_mat mergeArchetypes(mat C_stacked, mat H_stacked) {
		int arch_no = H_stacked.n_rows;
		
		// Merge all archetype-associated cells and compute their inner distances
		uvec all_cells;
		vector <uvec> cells(arch_no);
		for(int i = 0; i < arch_no; i++) {
			cells[i] = find(C_stacked.col(i) > 0);			
			all_cells = join_cols(all_cells, cells[i]);
		}
		all_cells = unique(all_cells);
		
		printf("Computing distances between %d influential cells\n", all_cells.n_elem);
		mat subH = H_stacked.cols(all_cells);		
		mat subD = computeFullDist(subH, 1, 0);		

		double sigma = stddev(vectorise(subD));
		
		ivec CC;
		uvec idx_i, idx_j, tmp;		

		printf("Computing distances between archetypes\n");
		mat interD = -datum::inf*ones(arch_no, arch_no);
		for(int i = 0; i < arch_no; i++) {
			uvec cells_i = cells[i];
			intersect(CC, idx_i, tmp, conv_to<ivec>::from(all_cells), conv_to<ivec>::from(cells_i));						
			
			for(int j = i+1; j < arch_no; j++) {
				uvec cells_j = cells[j];
				intersect(CC, idx_j, tmp, conv_to<ivec>::from(all_cells), conv_to<ivec>::from(cells_j));

				vec d_i = vectorise(subD(idx_i, idx_i));
				vec d_j = vectorise(subD(idx_j, idx_j)); 

				vec d_intra = join_cols(d_i, d_j);
				double mu_intra = mean(d_intra);
				double n_intra = d_intra.n_elem;

				vec d_inter = vectorise(subD(idx_i, idx_j)); 
				double mu_inter = mean(d_inter);
				double n_inter = d_inter.n_elem;

				double delta = (mu_inter - mu_intra);
				if(delta != 0)
					interD(i, j) = delta / (sigma*sqrt( (1.0/n_intra) + (1.0/n_inter) ));
				
			}
		}
		
		interD = max(interD, trans(interD));
		interD.diag() = datum::inf*ones(arch_no);
		
		vec merged = zeros(arch_no);
		sp_mat T(arch_no, arch_no);
		
		printf("Merge archetypes\n");
		for(int k = 0; k < arch_no-1; k++) {
			int idx = interD.index_min();
			uvec sub = ind2sub( size(interD), idx);
			int i = sub(0); int j = sub(1);
			
			int min_idx = min(i, j);
			int max_idx = max(i, j);
			
			T(max_idx, min_idx) = interD(idx);
			
			// Merge min_idx into max_idx
			cells[min_idx] = unique(join_vert(cells[min_idx], cells[max_idx]));
			merged[max_idx] = 1;
			
			interD.col(max_idx) = datum::inf*ones(arch_no);
			interD.row(max_idx) = trans(datum::inf*ones(arch_no));
			
			intersect(CC, idx_i, tmp, conv_to<ivec>::from(all_cells), conv_to<ivec>::from(cells[min_idx]));						
			for(int l = 0; l < arch_no; l++) {
				if(merged[l] || (l == min_idx))
					continue;
				else {
					intersect(CC, idx_j, tmp, conv_to<ivec>::from(all_cells), conv_to<ivec>::from(cells[l]));						

					vec d_i = vectorise(subD(idx_i, idx_i));
					vec d_j = vectorise(subD(idx_j, idx_j)); 

					vec d_intra = join_cols(d_i, d_j);
					double mu_intra = mean(d_intra);
					double n_intra = d_intra.n_elem;

					vec d_inter = vectorise(subD(idx_i, idx_j)); 
					double mu_inter = mean(d_inter);
					double n_inter = d_inter.n_elem;

					if(mu_inter == 0) {
						interD(min_idx, l) = -datum::inf;
						continue;
					}

					double delta = (mu_inter - mu_intra);
					if(delta != 0)
						interD(min_idx, l) = interD(l, min_idx) = delta / (sigma*sqrt( (1.0/n_intra) + (1.0/n_inter) ));
				}
			}
		}

		return(T);
	}
		
}
