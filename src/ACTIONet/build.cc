#include <actionetcore.h>
#include <atria.h>
#include <vptree.h>

namespace ACTIONetcore {
	mat computeFullDist(mat &H_stacked, int thread_no=-1, int verbose = 1) {	
		if(verbose > 0)
			printf("Building full distance matrix (returns full distance matrix)\n");
			
		int sample_no = H_stacked.n_cols;
		int archetype_no = H_stacked.n_rows;;

		mat D(sample_no, sample_no);
		
		H_stacked.transform( [](double val) { return (val < 0?0:val); } );
		mat H_stacked_norm = normalise(H_stacked, 1);		
		
		mat logProb = log(H_stacked_norm);
		logProb.transform( [](double val) { return (__builtin_isinf(val)?0:val); } );
		
		vec entropies = -trans(sum(H_stacked_norm % logProb, 0));
				
		double scaler = 1.0 / (2.0*log(2.0));

		int perc = 1;
		int total_counts = 0;
		#pragma omp parallel for num_threads(thread_no) 			
		for(int v = 0; v < sample_no; v++) {
			total_counts ++;
			if(round(100*(double)total_counts / sample_no) > perc) {
				if(verbose > 0)
					printf("%d %%\n", perc);
				perc++;
			}			
			
			vec p = H_stacked_norm.col(v);					
			mat logM(H_stacked_norm.n_rows, sample_no);
			for(int c = 0; c < sample_no; c++) {
				logM.col(c) = log(0.5*(p + H_stacked_norm.col(c)));
			}		
			logM.transform( [](double val) { return ((__builtin_isinf(val) || __builtin_isnan(val))?0:val); } );
					
			vec DpM = trans(-p.t()*logM - entropies(v));
			vec DQM = trans(-sum(H_stacked_norm % logM, 0)) - entropies;
			
			vec JS_div = scaler*(DpM + DQM);
			JS_div(v) = 0;
		
			D.col(v) = sqrt(JS_div); // Sqrt of JS Div, not JS Div, is a metric
		}		
		D.transform( [](double val) { return (__builtin_isnan(val)?0:val); } );
		D.transform( [](double val) { val = val < 0? 0:val; val = 1 < val? 1:val; return (val); } ); // Make sure that scores are normalized properly
		
		if(verbose > 0)
			printf("done\n");
			
		return D;
	}	
	
	
	sp_mat computeNearestDist(mat &H_stacked, int kNN, int thread_no=-1) {	
		printf("Building distance matrix of the %d nearest neighbors of each node (returns sparse distance matrix)\n", kNN);
		
		double epsilon = 1e-10;
		int sample_no = H_stacked.n_cols;
		
		if(kNN >= sample_no || kNN <= 1)
			kNN = min(30, sample_no-1);
			
		umat subs(2, kNN*sample_no);
		vec vv(kNN*sample_no);

		H_stacked.transform( [](double val) { return (val < 0?0:val); } );

		mat H_stacked_norm = normalise(H_stacked, 1, 0); // make the norm (sum) of each column 1			
		int archetype_no = H_stacked_norm.n_rows;
		
		// build ball tree on set
		VpTree<DataPoint, JSDiv_sqrt_distance>* tree = new VpTree<DataPoint, JSDiv_sqrt_distance>();
		std::vector<DataPoint> samples(sample_no); //, DataPoint(archetype_no, -1, data));
		for (int i = 0; i < sample_no; i++) {
			samples[i] = DataPoint(archetype_no, i, H_stacked_norm.colptr(i));
			//(H_stacked_norm.col(i)).print("col");
		}
		tree -> create(samples, 0);
		
		
		int perc = 1;
		int total_counts = 1;
		#pragma omp parallel num_threads(thread_no) 
		{
			vector< vector<int> > ind_arr(sample_no, std::vector<int>(kNN+1));
			vector< vector<double> > dist_arr(sample_no, std::vector<double>(kNN+1));
			#pragma omp for
			for (int v = 0; v < sample_no; v++) {
				total_counts ++;
				if(round(100*(double)total_counts / sample_no) > perc) {
					printf("%d %%\n", perc);
					perc++;
				}

				
				tree -> search(samples[v], kNN+1, &ind_arr[v], &dist_arr[v]);
				
				int base = v*kNN;
				for(int i = 1; i <= kNN; i++) {
					double d = dist_arr[v][i]; // To match with old scores
					d = d < epsilon?epsilon:d;
						
					subs(0, base + i-1) = ind_arr[v][i]-1;
					subs(1, base + i-1) = v;
					vv(base + i-1) = d;			
				}				
			}
		}
		samples.clear();
		delete tree;
			
		sp_mat D(subs, vv, sample_no, sample_no);	

		return(D);
	}
	

	field<mat> computeNearestDist_edgeList(mat &H_stacked, int kNN, int thread_no=-1) {	
		printf("Building distance matrix of the %d nearest neighbors of each node (returns edge list)\n", kNN);
		
		double epsilon = 1e-10;
		int sample_no = H_stacked.n_cols;
		
		if(kNN >= sample_no || kNN <= 1)
			kNN = min(30, sample_no-1);
			
		umat subs(2, kNN*sample_no);
		vec vv(kNN*sample_no);

		H_stacked.transform( [](double val) { return (val < 0?0:val); } );
		mat H_stacked_norm = normalise(H_stacked, 1, 0); // make the norm (sum) of each column 1			
		
		int archetype_no = H_stacked_norm.n_rows;
		
		// build ball tree on set
		VpTree<DataPoint, JSDiv_sqrt_distance>* tree = new VpTree<DataPoint, JSDiv_sqrt_distance>();
		std::vector<DataPoint> samples(sample_no); //, DataPoint(archetype_no, -1, data));
		for (int i = 0; i < sample_no; i++) {
			samples[i] = DataPoint(archetype_no, i, H_stacked_norm.colptr(i));
			//(H_stacked_norm.col(i)).print("col");
		}
		tree -> create(samples, 0);
		
		
		mat idx = zeros(sample_no, kNN+1);
		mat dist = zeros(sample_no, kNN+1);

		int perc = 1;
		int total_counts = 0;
		#pragma omp parallel num_threads(thread_no) 			
		{
			vector< vector<int> > ind_arr(sample_no, std::vector<int>(kNN+1));
			vector< vector<double> > dist_arr(sample_no, std::vector<double>(kNN+1));
			#pragma omp for
			for (int v = 0; v < sample_no; v++) {
				total_counts ++;
				if(round(100*(double)total_counts / sample_no) > perc) {
					printf("%d %%\n", perc);
					perc++;
				}
				
				tree -> search(samples[v], kNN+1, &ind_arr[v], &dist_arr[v]);
				idx(v, 0) = v+1;
				dist(v, 0) = 0;
							
				int base = v*kNN;
				for(int i = 1; i <= kNN; i++) {
					double d = dist_arr[v][i]; 
					d = d < epsilon?epsilon:d;
					
					#pragma omp atomic write 
					idx(v, i) = ind_arr[v][i];

					#pragma omp atomic write 
					dist(v, i) = d;
				}				
			}
		}			
		
		samples.clear();
		delete tree;
			
		field<mat> output(2);
		output(0) = idx;
		output(1) = dist;
		
		return(output);
	}


	
	field<sp_mat> buildACTIONet(mat &H_stacked, int kNN, int thread_no=-1) {	
		printf("Building ACTIONet\n");
		int sample_no = H_stacked.n_cols;		
		if(kNN >= sample_no || kNN < 1)
			kNN = min(30, sample_no-1);

		sp_mat D = computeNearestDist(H_stacked, kNN, thread_no);

		sp_mat G = D;
		for(sp_mat::iterator it = G.begin(); it != G.end(); ++it) {
		  (*it) = 1.0 - (*it);
		}					
		
		
		field<sp_mat> output(2);
		output(0) = sqrt(G % trans(G));
		output(1) = G;
		
		return(output);	
	}


	// k^{*}-Nearest Neighbors: From Global to Local (NIPS 2016)
	field<sp_mat> buildAdaptiveACTIONet(mat &H_stacked, double LC = 1.0, double epsilon = 0.0, int thread_no=4, bool auto_adjust_LC = false, int sym_method = ACTIONet_AND) {

		printf("Building adaptive ACTIONet (Eps = %.2f, LC = %.2f, Auto adjust LC = %d)\n", epsilon, LC, auto_adjust_LC);

		H_stacked.transform( [](double val) { return (val < 0?0:val); } );
		H_stacked = normalise(H_stacked, 1, 0); // make the norm (sum) of each column 1			

		double kappa = 5.0;
		int sample_no = H_stacked.n_cols;		
		int kNN = min(sample_no-1, (int)(kappa*round(sqrt(sample_no)))); // start with uniform k=sqrt(N) ["Pattern Classification" book by Duda et al.]


		mat idx = zeros(sample_no, kNN+1);
		mat dist = zeros(sample_no, kNN+1);
		
		int total_counts = 0, perc = 0;
		#pragma omp parallel num_threads(thread_no) 
		{			
			Searcher *searcher = new Searcher(H_stacked, "jensen", 0, 64, 0);			
			#pragma omp for
			for (long n = 0; n < sample_no; n++) {
				if(round(100*(double)total_counts / sample_no) > perc) {
					printf("%d %%\n", perc);
					perc++;
				}

				
				vector<neighbor> v;
				//mat::col_iterator it = H_stacked.begin_col(n);
				vec h = H_stacked.col(n);
				searcher->search_k_neighbors(v, kNN+1, h.begin(), -1, -1, epsilon);

									 
				for (long d = 0; d < v.size(); d++) {
					#pragma omp critical
					idx(n, d) = v[d].index() + 1; // Convert back to one-based indexing.
					
					#pragma omp critical
					dist(n, d) = v[d].dist();
				}
				
				total_counts ++;
			}
			
			delete searcher;
		}

		
		dist = clamp(dist, 0.0, 1.0); 
		idx = clamp(idx, 0, sample_no - 1);
		
		printf("\tConstructing adaptive-nearest neighbor graph ... \n");
		mat Delta;
		do {
			mat beta = LC*dist;
			vec beta_sum = zeros(sample_no);
			vec beta_sq_sum = zeros(sample_no);
			
			mat lambda = zeros(size(beta));
			//lambda.col(0) = datum::inf*ones(sample_no);
			//lambda.col(1) = beta.col(1) + 1;			
			register int k;
			for(k = 1; k <= kNN; k++) {
				beta_sum += beta.col(k);
				beta_sq_sum += square(beta.col(k));
				
				lambda.col(k) = (1.0/(double)k) * ( beta_sum + sqrt(k + square(beta_sum) - k*beta_sq_sum) );
			}
			lambda.replace(datum::nan, 0); 
			
			
			lambda = trans(lambda);
			vec node_lambda = zeros(sample_no);
			beta = trans(beta);

			
			Delta = lambda - beta;		
			Delta.shed_row(0);
			
/*
			vec saturation = sum(Delta < 0)
			vec saturation_sgn = saturation.transform( [](double val) { return (val == 0? 1:0; } );			
			double saturated_vertices = sum(saturation_sgn);
			*/
			vec saturation_mask = conv_to<vec>::from(sum(Delta < 0) == 0);
			double saturated_vertices_count = (double)sum(saturation_mask);
			if(auto_adjust_LC && saturated_vertices_count > round(0.01*H_stacked.n_cols)) {
				LC *= 1.1;
				printf("\t\t# saturated vertices = %.1f. Increasing LC to %.2f\n", saturated_vertices_count, LC);
			}
			else {
				break;
			}
		} while(1);
		
		
		sp_mat G(sample_no, sample_no);		
//		# pragma omp parallel for shared(G) num_threads(thread_no)
		for(int v = 0; v < sample_no; v++) {				
			vec delta = Delta.col(v);		
					
			//uvec rows = find(delta > 0, 1, "last");
			uvec rows = find(delta < 0, 1, "first");
			int neighbor_no = rows.n_elem == 0?kNN:(rows(0));
			
			int dst = v;								
			rowvec v_dist = dist.row(v);
			rowvec v_idx = idx.row(v) - 1;
			for (int i = 1; i < neighbor_no; i++) {				
				int src = v_idx(i);
					
				G(src, dst) = 1.0 - v_dist(i);
			} 
		}
		printf("\tdone\n");
		

		
		printf("\tFinalizing network ... ");
		G.replace(datum::nan, 0);  // replace each NaN with 0

		sp_mat Gt = trans(G);	
		
		sp_mat G_sym;
		if(sym_method == ACTIONet_AND) {
			printf("-- AND -- symmetrize\n");
			G_sym = sqrt(G % Gt);
		} else if(sym_method == ACTIONet_OR) {
			printf("-- OR -- symmetrize\n");
			G_sym = (G + Gt);
			G_sym.for_each( [](sp_mat::elem_type& val) { val /= 2.0; } );			
		} else { // Default to MNN
			printf("-- AND -- symmetrize (default)\n");
			sp_mat G_sym = sqrt(G % Gt);
		}
		
		
		sp_mat G_asym = normalise(G, 1, 0);
		
		field<sp_mat> output(2);
		output(0) =	G_sym;
		output(1) = G_asym;
		
		printf("done\n");
		
		return(output);		
	}


	mat constructKstarNN_fromDist(mat &D, double LC = 1.0) {
		printf("Running k*-nearest neighbors algorithm (LC = %f)\n", LC);


		mat G_prime = zeros(size(D));
		
		int nV = D.n_rows;
		for(int i = 0; i < nV; i++) {
			printf("col %d\b", i);
			vec d = D.col(i);
			
			uvec perm = sort_index(d, "ascend");			
			vec beta = LC * d(perm);		
			
			double k = 0, Sum_beta = 0, Sum_beta_square = 0;			
			double lambda = 0;
			
			double last_lambda = 0;
			for(k = 1; k <= nV; k++) {
				last_lambda = lambda;
				
				Sum_beta += beta(k-1);
				Sum_beta_square += std::pow(beta(k-1), 2);
				lambda = (1.0 / k) * ( Sum_beta + sqrt( k  + std::pow(Sum_beta, 2) - k * Sum_beta_square ) ) ; 				
				
				printf("\t%d- %f, %f, %f, %f, %f\n", (int)k, Sum_beta, Sum_beta_square, k * Sum_beta_square, std::pow(Sum_beta, 2) - k * Sum_beta_square, lambda);
				
				//if(lambda < beta(k-1))
					//break;
			}
			/*int knn = (int)(k-1);
						
			vec sub_beta = beta(span(0, knn-1));
			vec w = last_lambda-sub_beta;//last_lambda - beta(span(0, knn));
			
			w /= sum(w);
									
			vec v = zeros(nV);
			v(perm(span(0, knn-1))) = w;
			
			G_prime.col(i) = v;*/
		}

/*	
		printf("LC = %.2f\n", LC);
		
		int sample_no = D.n_cols;
		mat D_sorted = sort(D, "ascend", 1);
		mat beta = LC*D_sorted;
		vec beta_sum = zeros(sample_no);
		vec beta_sq_sum = zeros(sample_no);
		
		mat lambda = zeros(size(beta));
		//lambda.col(0) = datum::inf*ones(sample_no);
		//lambda.col(1) = beta.col(1) + 1;			
		for(int j = 0; j < sample_no; j++) {
			double k = j + 1.0;
			beta_sum += beta.col(j);
			beta_sq_sum += square(beta.col(j));
			
			lambda.col(j) = (1.0/k) * ( beta_sum + sqrt(k + square(beta_sum) - k*beta_sq_sum) );
		}
		return(lambda);
		
		lambda.replace(datum::nan, 0); 		
		
		mat delta = lambda - beta;				
		return(delta);

		delta.transform( [](double val) { return (val < 0?0:val); } );
		delta = normalise(trans(delta), 1);

		
		mat G_prime = zeros(size(delta));
		for(int v = 0; v < sample_no; v++) {				
			vec delta = delta.col(v);		

			vec d = D.col(v);			
			uvec perm = sort_index(d, "ascend");			
			uvec rows = find(delta > 0, 1, "last");
			
			vec x = zeros(sample_no);
			x(perm(rows))  = delta(rows);
			G_prime.col(v) = x;
		}



		
		//mat G_sym = sqrt(G_prime % trans(G_prime));

		//return(G_sym);
		
		return(G_prime);
		*/
	}

}






