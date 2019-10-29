#include <actionetcore.h>

#include <prpack_utils.h>
#include <prpack_result.h>
#include <prpack_solver.h>

vec pagerank(sp_mat &G, vec u_vec, double alpha = 0.85, double tol = 1e-6) {

	prpack::prpack_base_graph g(&G);

	g.normalize_weights(); 
	
	
    prpack::prpack_solver solver(&g, false);

	u_vec = normalise(u_vec, 1);


	double* u = new double[u_vec.size()];
	memcpy(u, u_vec.memptr(), u_vec.size()*sizeof(double));
    double* v = u;
    
    
    const prpack::prpack_result* res = solver.solve(alpha, tol, u, v, "");
    

    vec pr(res->x, res->num_vs);
        
    return pr;
}


vec pagerank_scaled(sp_mat &G, vec u_vec, double alpha = 0.85, double tol = 1e-6) {
	
	prpack::prpack_base_graph g(&G);
	g.normalize_weights(); 
    prpack::prpack_solver solver(&g, false);

	u_vec = normalise(u_vec, 1);
	
	double* u = new double[u_vec.size()];
	memcpy(u, u_vec.memptr(), u_vec.size()*sizeof(double));
    double* v = u;
    
    const prpack::prpack_result* res = solver.solve(alpha, tol, u, v, "");    
    vec pr(res->x, res->num_vs);
    
    
    vec o = ones(size(u_vec)) / (double)u_vec.n_elem;	
	memcpy(u, o.memptr(), o.size()*sizeof(double));
    v = u;
        
    res = solver.solve(alpha, tol, u, v, "");    
    vec pr0(res->x, res->num_vs);
    
    vec log_ratio = log2(pr / pr0);
    
    log_ratio.replace(datum::nan, 0); 

    return log_ratio;
}

namespace ACTIONetcore {
	mat batch_zoned_diffusion(sp_mat &G, uvec &zones, mat &U, double alpha = 0.85, int thread_no = 8, double tol = 1e-6) {
		mat ZonedDiff = zeros(size(U));
		
		if(U.n_rows != G.n_rows) {
			fprintf(stderr, "batch_zoned_diffusion:: Number of rows in U doesn't match with the number of vertices in G\n");
			return(ZonedDiff);
		}
		U.transform( [](double val) { return (val < 0? 0:val); } );
		U = normalise(U, 1);
		
		G.diag().zeros();
		
		uvec uz = unique(zones);					
		printf("Zoned diffusion (%d zones)\n", uz.n_elem);
		
		uvec vertex_id(size(zones));	
		for(int i = 0; i < uz.n_elem; i++) {
			uvec idx = find(zones == uz(i));

			for(int j = 0; j < idx.n_elem; j++) {
				vertex_id(idx(j)) = j;
			}
					
			sp_mat subG(idx.n_elem, idx.n_elem);
			
			sp_mat::iterator it     = G.begin();
			sp_mat::iterator it_end = G.end();
			for(; it != it_end; ++it) {
				if( (zones(it.row()) != uz(i)) || (zones(it.col()) != uz(i)) )
					continue;
				
				int ii = vertex_id(it.row());
				int jj = vertex_id(it.col());

				subG(ii, jj) = subG(jj, ii) = (*it);
			}
			
			prpack::prpack_base_graph g(&subG);		
			g.normalize_weights(); 
			prpack::prpack_solver solver(&g, false);

			
			#pragma omp parallel for num_threads(thread_no) 			
			for(int k = 0; k < U.n_cols; k++) {
				vec init_pr = U.col(k);
				//double total_max = max(init_pr);
				double weight = sum(init_pr(idx));				
				
				
				vec sub_init_pr = normalise(init_pr(idx), 1);
				double* u = new double[sub_init_pr.size()];
				memcpy(u, sub_init_pr.memptr(), sub_init_pr.size()*sizeof(double));
				double* v = u;
			
				const prpack::prpack_result* res = solver.solve(alpha, tol, u, v, "");
											
				vec pr(res->x, res->num_vs);
				pr.replace(datum::nan, 0);

				vec x = ZonedDiff.col(k);
				x(idx) = weight*pr;
				ZonedDiff.col(k) = x;
			}			
		}

		ZonedDiff = normalise(ZonedDiff, 1);
		return(ZonedDiff);	
	}


	
	
	
	
	mat batchPR(sp_mat &G, mat &U, double alpha = 0.85, int thread_no = 8, double tol = 1e-6) {
		
		prpack::prpack_base_graph g(&G);		
		
		g.normalize_weights(); 
		prpack::prpack_solver solver(&g, false);
		
		U = normalise(U, 1);
		
		mat PR = zeros(size(U));
		
		
		#pragma omp parallel for num_threads(thread_no) 			
		for(int i = 0; i < U.n_cols; i++) {			
			vec u_vec = U.col(i);
			double* u = new double[u_vec.size()];
			memcpy(u, u_vec.memptr(), u_vec.size()*sizeof(double));
			double* v = u;
			
			const prpack::prpack_result* res = solver.solve(alpha, tol, u, v, "");
			
			vec pr(res->x, res->num_vs);
			PR.col(i) = pr;
		}	
		
		PR.replace(datum::nan, 0);
		uvec isolated_vertices = find(sum(PR) < 1e-6);
		PR.cols(isolated_vertices) = normalise(ones(PR.n_rows, isolated_vertices.n_elem), 1);
		
		return PR;
	}

	
	mat extractArchetypeAssociatedSamples(sp_mat &ACTIONet, mat H_stacked, double alpha) {
		int nV = ACTIONet.n_rows;
		int arch_no = H_stacked.n_rows;
		int sample_no = H_stacked.n_cols;
		
		double tol = 1e-6;

		double *u, *v;
		u = new double[nV];
		
		prpack::prpack_base_graph g(&ACTIONet);
		g.normalize_weights(); 
		prpack::prpack_solver solver(&g, false);

		vec o = ones(nV) / (double)nV;
		memcpy(u, o.memptr(), o.size()*sizeof(double));
		v = u;			
		prpack::prpack_result* res = solver.solve(alpha, tol, u, v, "");    
		vec pr(res->x, res->num_vs);


		mat scores = zeros(sample_no, arch_no);
		for(int i = 0; i < arch_no; i++) {
			printf("Archetype %d/%d\n", i, arch_no);
			vec h = trans(normalise(H_stacked.row(i), 1));
	
			memcpy(u, h.memptr(), h.size()*sizeof(double));
			double* v = u;
    
			res = solver.solve(alpha, tol, u, v, "");    
			vec ppr(res->x, res->num_vs);
        
    
			vec log_ratio = log2(ppr / pr);    
			log_ratio.replace(datum::nan, 0); 
			scores.col(i) = log_ratio;			
		}
	
		return(scores);
	}	
	
	vec extractArchetypeAssociatedSamples(sp_mat &ACTIONet, vec h, double alpha) {
		int nV = ACTIONet.n_rows;

		printf("Computing personalized PageRank vector ... "); fflush(stdout);
		vec log_ratio = pagerank_scaled(ACTIONet, h, alpha);	
		printf("done\n"); fflush(stdout);
		
		log_ratio.transform( [](double val) { return (val < 0? 0:val); } );


/*
		printf("Computing personalized PageRank vector ... "); fflush(stdout);
		vec pr = pagerank(ACTIONet, h, alpha);		
		printf("done\n"); fflush(stdout);

		vec w = vec(sum(ACTIONet, 1));		
		pr = pr / w; // normalize by degrees		
		uvec indices = sort_index(pr, "descend");
		double W = sum(w);
		
		printf("Computing modularity of prefix sets ... "); fflush(stdout);		
		vec modularity = zeros(nV);
		for (int i = 1; i < indices.n_elem-1; i++) {
			modularity(i) = modularity(i-1);
			
			vec N = vec(ACTIONet.col(indices(i)));
			double w_i = w(indices(i));
			for(int j = 0; j < i; j++) {
				double w_ij = N(indices(j));				
				double w_j = w(indices(j));
				modularity(i) += (w_ij - (w_i*w_j / W));
			}
		}		
		int max_idx = index_max(modularity);
		printf("done (idx = %d, modularity = %.2e)\n", max_idx, modularity(max_idx)); fflush(stdout);
		
		pr(indices(span(max_idx+1, nV-1))) = zeros(nV-max_idx-1);


		
		printf("Computing conductance of prefix sets ... "); fflush(stdout);
		vec cond(nV);

		double cut_C = w(indices(0));
		double links_C = cut_C;		
		cond(0) = cut_C / std::min(links_C, W-links_C);		
		for (int i = 1; i < indices.n_elem-1; i++) {
			vec c = vec(ACTIONet.col(indices(i)));
			
			cut_C += sum(c(indices(span(i+1, nV-1))));			
			links_C += w(indices(i));
			
			cond(i) = cut_C / std::min(links_C, W-links_C);
		}		
		int min_idx = index_min(cond(span(0, nV-2)));
		printf("done (idx = %d, min = %.2e)\n", min_idx, cond(min_idx)); fflush(stdout);
		
		pr(indices(span(min_idx+1, nV-1))) = zeros(nV-min_idx-1);
		
		pr = pr / max(pr);
		return pr;
*/
		log_ratio /= max(log_ratio);
		
		return(log_ratio);
	}
	
	
	/*
	mat regressMultilevelArchetypes(mat S_r, mat C_stacked) {
		mat A = zscore(S_r);
		mat B = zscore(S_r * C_stacked);
				
		mat K = t(B)*B; // Gram matrix -- C
		mat R = t(B) * A;
		
		mat I = eye(size(K));
		mat KpI_inv = inv_sympd(K + I);
		
		mat denom_fact = K*CpI_inv
		
		mat Hs = zeros(C_stacked.n_cols, S_r.n_cols);
		
		for(int i = 0; i < Hs.n_cols; i++) {
			vec r = R.col(i);
			double factor = dot(r, r) / (t(r) * denom_fact * r);
			
		}
		
	}
	*/


	double *l1, *l2, *w;	
	int *match1, *match2, *v1, *v2;	
	int *s, *tt, *deg, *offset, *list;

	/**
	 * n the number of nodes
	 * m the number of nodes
	 * nedges the number of edges
	 * vv1 is the source for each of the nedges 
	 * vv2 is the target for each of the nedges
	 * weight is the weight of each of the nedges
	 * out1 is a vector of length at most min(n,m),
	 * out2 is a vector of length at most min(n,m),
	 * noutedges is the number of out edges
	 */
	double MWM_driver(int n, int m, int nedges, double *vv1, double *vv2, double *weight, double *out1, double *out2, int *noutedges) {	
		double ret, al;
		int i, j, k, p, q, r, t1, t2;
			

		for (i = 0; i < nedges; i++) {
			v1[i] = (int)(vv1[i] + .5);
			v2[i] = (int)(vv2[i] + .5);
		}
		for (i = 0; i < n; i++) {
			offset[i] = 0;
			deg[i] = 1;
		}
		for (i = 0; i < nedges; i++) deg[v1[i]]++;
		for (i = 1; i < n; i++) offset[i] = offset[i-1] + deg[i-1];
		for (i = 0; i < n; i++) deg[i] = 0;
		for (i = 0; i < nedges; i++) {
			list[offset[v1[i]] + deg[v1[i]]] = v2[i];
			w[offset[v1[i]] + deg[v1[i]]] = weight[i];
			deg[(int)v1[i]]++;
		}
		for (i = 0; i < n; i++) {
			list[offset[i] + deg[i]] = m + i;
			w[offset[i] + deg[i]] = 0;
			deg[i]++;
		}
		for (i = 0; i < n; i++) {
			l1[i] = 0;
			for (j = 0; j < deg[i]; j++) {
				if (w[offset[i]+j] > l1[i]) l1[i] = w[offset[i] + j];
			}
		}
		for (i = 0; i < n; i++) {
			match1[i] = -1;
		}
		for (i = 0; i < n + m; i++) {
			l2[i] = 0;
			match2[i] = -1;
		}
		for (i = 0; i < n; i++) {
			for(j = 0; j < n + m; j++) tt[j] = -1;
			s[p = q = 0] = i;
			for(; p <= q; p++) {
				if (match1[i] >= 0) break;
				k = s[p];
				for (r = 0; r < deg[k]; r++) {
					if (match1[i] >= 0) break;
					j = list[offset[k] + r];
					if (w[offset[k] + r] < l1[k] + l2[j] - 1e-8) continue;
					if (tt[j] < 0) {
						s[++q] = match2[j];
						tt[j] = k;
						if (match2[j] < 0) {
							for(; j>=0 ;) {
								k = match2[j] = tt[j];
								p = match1[k];
								match1[k] = j;
								j = p;
							}
						}
					}
				}
			}
			if (match1[i] < 0) {
				al = 1e20;
				for (j = 0; j < p; j++) {
					t1 = s[j];
					for (k = 0; k < deg[t1]; k++) {
						t2 = list[offset[t1] + k];
						if (tt[t2] < 0 && l1[t1] + l2[t2] - w[offset[t1] + k] < al) {
							al = l1[t1] + l2[t2] - w[offset[t1] + k];
						}
					}
				}
				for (j = 0; j < p; j++) l1[s[j]] -= al;
				for (j = 0; j < n + m; j++) if (tt[j] >= 0) l2[j] += al;
				i--;
				continue;
			}
		}
		ret = 0;
		for (i = 0; i < n; i++) {
			for (j = 0; j < deg[i]; j++) {
				if (list[offset[i] + j] == match1[i]) {
					ret += w[offset[i] + j];
				}
			}
		}
		*noutedges = 0;
		for (i = 0; i < n; i++) {
			if (match1[i] < m) (*noutedges)++;
		}
		*noutedges = 0;
		for (i = 0; i < n; i++) {
			if (match1[i] < m) {
				out1[*noutedges] = i;
				out2[*noutedges] = match1[i];
				(*noutedges)++;
			}
		}
		
		return ret;
	}


	mat MWM(mat &G) {
		int n = G.n_rows;
		int m = G.n_cols;
		int N = m+n;
		l1 = new double[N];
		l2 = new double[N];
		v1 = new int[m*n];
		v2 = new int[m*n];
		s = new int[N];
		tt = new int[N];
		match1 = new int[N];
		match2 = new int[N];
		offset = new int[N];
		deg = new int[N];
		list = new int[m*n+N];
		w = new double[m*n+N];		
		
		
		mat G_matched = zeros(size(G));

		
		uvec idx = find(G);
		if(idx.n_elem == 0)
			return G_matched;
		
		int nedges = idx.n_elem;
			
		double *vv1 = new double[nedges];
		double *vv2 = new double[nedges];
		double *weight = new double[nedges];
		
		umat subs = ind2sub( size(G), idx );
		for (int i = 0; i < nedges; i++) {
			weight[i] = G(idx(i));
			vv1[i] = subs(0, i);
			vv2[i] = subs(1, i);		
		}
		
		int match_no = std::min(m, n);
		double *ii = new double[match_no];
		double *jj = new double[match_no];
		
		int matched_edge_no;
		
		MWM_driver(n, m, nedges, vv1, vv2, weight, ii, jj, &matched_edge_no);	
		
		for (int k = 0; k < matched_edge_no; k++) {
			G_matched(ii[k], jj[k]) = G(ii[k], jj[k]);
		}
		
		
		delete [] l1;
		delete [] l2;
		delete [] v1;
		delete [] v2;
		delete [] s;
		delete [] tt;
		delete [] match1;
		delete [] match2;
		delete [] offset;
		delete [] deg;
		delete [] list;
		delete [] w;			
		
		return G_matched;
	}


	umat Rank1_matching(vec u, vec v, double u_threshold, double v_threshold) {
		int pair_no = std::min(u.n_elem, v.n_elem);
		
		// Rank-1 matching for each paired-archetypes
		vec u_sorted = sort( u, "descend"); 
		uvec u_sorted_idx = sort_index( u, "descend"); 
		
		vec v_sorted = sort( v, "descend"); 
		uvec v_sorted_idx = sort_index( v, "descend"); 

	/*	
		// To randomly break ties
		u_sorted = u_sorted + stddev(u)*randn(size(u_sorted))/30;				
		u_sorted.transform( [](double val) { return (val < 0? 0:val); } );
		
		v_sorted = v_sorted + stddev(v)*randn(size(v_sorted))/30;
		v_sorted.transform( [](double val) { return (val < 0? 0:val); } );
	*/
			
		// To reduce false positves by aligning only high-quality pairs
		//vec uv_prod = sqrt(u_sorted(span(0, pair_no-1)) % v_sorted(span(0, pair_no-1)));
		
		
		
		int top_rank = min(sum(u > u_threshold), sum(v > v_threshold)); //sum(uv_prod > threshold);
		umat subs(2, top_rank);
			
		if(top_rank == 0)
			return subs;
			
		// Rank-based matching
		uvec rows = u_sorted_idx(span(0, top_rank-1));				
		uvec cols = v_sorted_idx(span(0, top_rank-1));
		
		subs.row(0) = (rows).t();
		subs.row(1) = (cols).t();	
		
		return subs;
	}		

	vec sweepcut(sp_mat &A, vec s) {
		int top_ignore = 5;
		
		A.diag().zeros();		
		int nV = A.n_rows;
		
		vec w = vec(sum(A, 1));
		double total_vol = sum(w);

		vec conductance = datum::inf*ones(w.n_elem);				
		
		uvec perm = sort_index( s, "descend" );
		vec x = zeros(nV);
		x(perm(span(0, top_ignore-1))).ones();
		double vol = sum(w(perm(span(0, top_ignore-1))));

		double cut_size = vol;
		for(int i = 0; i < top_ignore; i++) {
			for(int j = 0; j < top_ignore; j++) {
				cut_size -= A(i, j);
			}
		}
		
		for(register int i = top_ignore; i < nV-top_ignore-1; i++) {
			int u = perm(i);
			vol += w[u];
			
			x(u) = 1;			

			sp_mat::col_iterator it = A.begin_col( u );							
			for(; it != A.end_col( u ); it++) {
				int v = it.row();
				if(x[v] == 0) { // v is in S_prime (not yet selected)
					cut_size += (*it);
				}	
				else {
					cut_size -= (*it);
				}			
			}
						
			double vol_prime = total_vol - vol;
			conductance(i) = cut_size / min(vol, vol_prime);
		}

		
		return(conductance);
	}
	
	vec nonlinear_diffusion(sp_mat A, vec s, double p = 0.5, double h = 0.001, int k = 10) {
		sp_mat A_norm = normalise(A, 1, 0);
		sp_mat I = speye(size(A));
		sp_mat L = I - A_norm;
		
		printf("h = %e, p = %e\n", h, p);
		
		vec u = s;		
		for(register int i = 0; i < k; i++) {
			vec delta = h * L * pow(u, p);
			u = u - delta;
			u = clamp( u, 2.0e-8, 1 );
		}
		
		return(u);
	}



}
