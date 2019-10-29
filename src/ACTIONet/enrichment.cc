#include <actionetcore.h>

mat sampleUnif(int l, int m, double a, double b, int seed) {
	std::default_random_engine gen (seed);	
	std::uniform_real_distribution<double> unif(a, b);
	
	
	mat R(l, m);
	for (register int j = 0; j < m; j++) {
		for(register int i = 0; i < l; i++) {
			R(i, j) = unif(gen);
		}
	}
	return R;
}

namespace ACTIONetcore {
	mat assessFeatureSets(sp_mat &S, field<uvec> feature_sets, int rand_perm_no) {	

		int feature_no = S.n_rows;
		int sample_no = S.n_cols;
		int feature_set_no = feature_sets.n_elem;

		std::default_random_engine gen (0);	
		std::uniform_real_distribution<double> sample_next(0, feature_no-1);
	
		printf("Aggregating scores (S: %d x %d)\n", feature_no, sample_no); fflush(stdout);
		
		
		S = S.t(); // Faster to access column-wise
		mat scores = zeros(sample_no, feature_set_no);
		for (int i = 0; i < feature_set_no; i++) {
			uvec indices = feature_sets(i);
			printf("Computing stat for featureset %d/%d (size %d)\n", i+1, feature_set_no, indices.n_elem); fflush(stdout);
			
			vec stat = zeros(sample_no);
			for(int k = 0; k < indices.n_elem; k++) {
				stat += vec(S.col(indices(k)));
			}
			
			mat rand_stat = zeros(sample_no, rand_perm_no);
			for (int j = 0; j < rand_perm_no; j++) {
				
				for(int k = 0; k < indices.n_elem; k++) {
					//int rand_idx = rand() % feature_no;	
					int rand_idx = (int)round(sample_next(gen));
					
					rand_stat.col(j) += vec(S.col(rand_idx));
				}
				
			}						
			
			scores.col(i) = (stat - mean(rand_stat, 1)) / stddev(rand_stat, 0, 1);
		}
		
		return scores;
	}
	
	/*
	mat featureEnrichment(sp_mat &S, mat &feature_scores, int rand_perm_no) {		
		
		int feature_no = S.n_rows;
		int sample_no = S.n_cols;
		int condition_no = feature_scores.n_cols;
		mat feature_scores_norm = normalise(feature_scores, 1, 0);
		
		
		if(rand_perm_no == 0) {
			// Uniform_sum formulation From: Distribution and moments of the weighted sum of uniforms random variables, with applications in reducing monte carlo simulations		
			mat Z = zeros(S.n_rows, feature_scores.n_cols);		
			rowvec sigma = sqrt(sum(square(feature_scores_norm))/12);

			for (int i = 0; i < S.n_rows; i++) {
				rowvec r = conv_to<rowvec>::from(sort_index(S.row(i), "ascend"));
				r = r / r.n_elem;
				rowvec stats = r * feature_scores_norm;
				
				rowvec z = (stats - 0.5) / sigma;
				Z.row(i) = z;
			}
			
			return(Z);
		}
		else {
			// Compute the statistic
			mat stat = S * feature_scores_norm;
			
			mat randstat_sums = zeros(feature_no, condition_no);
			mat randstat_sums_sq = zeros(feature_no, condition_no);
			mat rand_samples = sampleUnif(sample_no, rand_perm_no, 0.0, 1.0, 0);
			
			uvec rand_perm = sort_index(rand_samples.col(0));			
			mat rand_stat = S * feature_scores_norm.rows(rand_perm);			
			mat K = rand_stat;
			
			for (int i = 1; i < rand_perm_no; i++) {
				rand_perm = sort_index(rand_samples.col(i));			
				rand_stat = S * feature_scores_norm.rows(rand_perm);
				
				mat shifted_val = (rand_stat - K);
				randstat_sums += shifted_val;
				randstat_sums_sq += (shifted_val % shifted_val);
			}
					
			mat mu = (randstat_sums + (double)(rand_perm_no-1)*K) / (double)rand_perm_no;	
			mat sigma = sqrt((randstat_sums_sq - square(randstat_sums)/((double)rand_perm_no)) / ((double)rand_perm_no-1));
			
			// arch x phenotype
			mat Z = (stat - mu) / sigma;
			
			return(Z);
		}
	}
	*/
	mat assessFeatureSets_archs(mat &archs, field<uvec> feature_sets, int rand_perm_no) {		
		srand (0);	

		int feature_no = archs.n_rows;
		int sample_no = archs.n_cols;
		int feature_set_no = feature_sets.n_elem;
	
		printf("Aggregating scores (archs: %d x %d)\n", feature_no, sample_no); fflush(stdout);
		
		
		archs = archs.t(); // Faster to access column-wise
		mat scores = zeros(sample_no, feature_set_no);
		for (int i = 0; i < feature_set_no; i++) {
			uvec indices = feature_sets(i);
			printf("Computing stat for featureset %d/%d (size %d)\n", i+1, feature_set_no, indices.n_elem); fflush(stdout);
			
			vec stat = zeros(sample_no);
			for(int k = 0; k < indices.n_elem; k++) {
				stat += vec(archs.col(indices(k)));
			}
			
			mat rand_stat = zeros(sample_no, rand_perm_no);
			for (int j = 0; j < rand_perm_no; j++) {
				
				for(int k = 0; k < indices.n_elem; k++) {
					int rand_idx = rand() % feature_no;					
					rand_stat.col(j) += vec(archs.col(rand_idx));
				}
				
			}						
			
			scores.col(i) = (stat - mean(rand_stat, 1)) / stddev(rand_stat, 0, 1);
		}
		
		return scores;
	}

	mat assessFeatureSets_decoupled(mat &archetype_profile, mat &H_stacked, field<uvec> feature_sets, int rand_perm_no) {
		
		sp_mat archetype_profile_sparse = sp_mat(archetype_profile);
		mat arch_scores = assessFeatureSets(archetype_profile_sparse, feature_sets, rand_perm_no);
		
		mat scores = trans(normalise(H_stacked, 1, 0)) *arch_scores;

		return scores;
	}



	mat phenotypeEnrichment(mat &H_stacked, mat &phenotype_associations, int rand_perm_no) {

		int archetype_no = H_stacked.n_rows;
		int sample_no = H_stacked.n_cols;
		int phenotype_no = phenotype_associations.n_cols;
		mat phenotype_associations_norm = normalise(phenotype_associations, 1, 0);
		
		
		if(rand_perm_no == 0) {
			// Uniform_sum formulation From: Distribution and moments of the weighted sum of uniforms random variables, with applications in reducing monte carlo simulations		
			mat Z = zeros(H_stacked.n_rows, phenotype_associations.n_cols);		
			rowvec sigma = sqrt(sum(square(phenotype_associations_norm))/12);

			for (int i = 0; i < H_stacked.n_rows; i++) {
				rowvec r = conv_to<rowvec>::from(sort_index(H_stacked.row(i), "ascend"));
				r = r / r.n_elem;
				rowvec stats = r * phenotype_associations_norm;
				
				rowvec z = (stats - 0.5) / sigma;
				Z.row(i) = z;
			}
			
			return(Z);
		}
		else {
			// Compute the statistic
			mat stat = H_stacked * phenotype_associations_norm;
			
			mat randstat_sums = zeros(archetype_no, phenotype_no);
			mat randstat_sums_sq = zeros(archetype_no, phenotype_no);
			mat rand_samples = sampleUnif(sample_no, rand_perm_no, 0.0, 1.0, 0);
			
			uvec rand_perm = sort_index(rand_samples.col(0));			
			mat rand_stat = H_stacked * phenotype_associations_norm.rows(rand_perm);			
			mat K = rand_stat;
			
			for (int i = 1; i < rand_perm_no; i++) {
				rand_perm = sort_index(rand_samples.col(i));			
				rand_stat = H_stacked * phenotype_associations_norm.rows(rand_perm);
				
				mat shifted_val = (rand_stat - K);
				randstat_sums += shifted_val;
				randstat_sums_sq += (shifted_val % shifted_val);
			}
					
			mat mu = (randstat_sums + (double)(rand_perm_no-1)*K) / (double)rand_perm_no;	
			mat sigma = sqrt((randstat_sums_sq - square(randstat_sums)/((double)rand_perm_no)) / ((double)rand_perm_no-1));
			
			// arch x phenotype
			mat Z = (stat - mu) / sigma;
			
			return(Z);
		}
	}

	sp_mat shuffleEdges(umat subs, vec vals, int nV, int num_shuffles = 10000, int seed = 0) {
		std::default_random_engine gen (seed);	
		std::uniform_real_distribution<double> unif(0, vals.n_elem);
		
		umat extra_edges(2, num_shuffles);
		vec extra_vals(num_shuffles);
		int idx = 0;
		
		for (register int k = 0; k < num_shuffles; k++) {
			int i = (int)round(unif(gen));
			int j = (int)round(unif(gen));				
			
			if(i == j) 
				continue;
				
			int a = subs(0, i);
			int b = subs(0, j);
			
			int a_prime = subs(1, i);
			int b_prime = subs(1, j);
			
			double w_a = vals(i);
			double w_b = vals(j);
			if(w_a == w_b) {
				continue;
			}
			
			subs(1, i) = b_prime;
			subs(1, j) = a_prime;
			
			if(w_b < w_a) {
				vals(i) = vals(j) = w_b;
				
				extra_vals(idx) = w_a-w_b;
				extra_edges(0, idx) = a;
				extra_edges(1, idx) = a_prime;								
				idx++;
			}
			else if(w_a < w_b) {
				vals(i) = vals(j) = w_a;
				
				extra_vals(idx) = w_b-w_a;
				extra_edges(0, idx) = b;
				extra_edges(1, idx) = b_prime;				
				idx++;
			}			
		}
		
		subs = join_horiz(subs, extra_edges.cols(0, idx-1));
		vals = join_vert(vals, extra_vals(span(0, idx-1)));
		
		sp_mat G_rand(true, subs, vals, nV, nV);
		return(G_rand);
	}
	
	// G is the cell-cell network, scores is a cell x geneset matrix
	field<vec> computeAutocorrelation (sp_mat &G, mat &scores, int rand_perm_no, int num_shuffles = 10000) {
		int nV = G.n_rows;
		int feature_set_no = scores.n_cols;
		
		printf("Computing auto-correlation over network of %d samples for %d scores\n", nV, feature_set_no);
		int nnz = G.n_nonzero;
		int idx = 0;		

		vec vals(nnz);
		umat subs(2, nnz);
		for(sp_mat::iterator it= G.begin(); it != G.end(); ++it) {
			vals(idx) = *it;
			subs(0, idx) = it.row();
			subs(1, idx) = it.col();
			idx++;
		}

		//double scale_factor = sum(sum(G)) / (sample_no-1); // 2W / (N-1)
		double total_weight = sum(vals);

		// Compute graph Laplacian
		vec d = vec(trans(sum(G)));
		
		sp_mat L(nV, nV);
		L.diag() = d;
		L -= G;
		
		// Summary stats
		vec Cstat = zeros(feature_set_no);


		
		printf("Computing autocorrelations ..."); fflush(stdout);
		for(int i = 0; i < feature_set_no; i++) {
			vec x = scores.col(i);
			double stat = dot(x, L*x); // x'Lx = sum_{i,j} (w_{ij}(x_i - x_j)^2)						
			double norm_fact = var(x)*total_weight;
			Cstat(i) = 1 - (stat / norm_fact);
		}
		printf("done\n"); fflush(stdout);
		
		
		vec mu, sigma, Cstat_Z;
		if(rand_perm_no > 0) {
			int perc = 0;
			mat rand_stat = zeros(rand_perm_no, feature_set_no);
			printf("Computing significance of autocorrelations ...\n"); fflush(stdout);		
			for(int k = 0; k < rand_perm_no; k++) {
				if(round(100.0*k / rand_perm_no) > perc) {
					printf("%d %%\n", perc);
					perc++;
				}
				
				
				/*
				sp_mat G_rand = shuffleEdges(subs, vals, nV, num_shuffles, 1000*k);			
								
				sp_mat L_rand(nV, nV);	
				L_rand.diag() = d;
				L_rand -= G_rand;
				*/
				
				for(int i = 0; i < feature_set_no; i++) {
					vec x = scores.col(i);
					uvec rand_perm = sort_index(randu(size(x)));
					vec x_permuted = x(rand_perm);					
					
					double stat = dot(x_permuted, L*x_permuted); // x'Lx = sum_{i,j} (w_{ij}(x_i - x_j)^2)						
					double norm_fact = var(x_permuted)*total_weight;
					rand_stat(k, i) = 1 - (stat / norm_fact);
				}
			}
			printf("done\n"); fflush(stdout);
			
			mu = trans(mean(rand_stat));
			sigma = trans(stddev(rand_stat));
			Cstat_Z = (Cstat - mu) / sigma;
		}
		else {
			mu = zeros(feature_set_no);
			sigma = zeros(feature_set_no);
			Cstat_Z = zeros(feature_set_no);
		}
		
		field<vec> results(4);
		results(0) = Cstat;		
		results(1) = mu;
		results(2) = sigma;
		results(3) = Cstat_Z;
		
		return(results);
	}

}
