#include "ACTION.h"

#include <thread>
#include <atomic>

template<class Function>
inline void ParallelFor(size_t start, size_t end, size_t numThreads, Function fn) {
    if (numThreads <= 0) {
        numThreads = std::thread::hardware_concurrency();
    }

    if (numThreads == 1) {
        for (size_t id = start; id < end; id++) {
            fn(id, 0);
        }
    } else {
        std::vector<std::thread> threads;
        std::atomic<size_t> current(start);

        // keep track of exceptions in threads
        // https://stackoverflow.com/a/32428427/1713196
        std::exception_ptr lastException = nullptr;
        std::mutex lastExceptMutex;

        for (size_t threadId = 0; threadId < numThreads; ++threadId) {
            threads.push_back(std::thread([&, threadId] {
                while (true) {
                    size_t id = current.fetch_add(1);

                    if ((id >= end)) {
                        break;
                    }

                    try {
                        fn(id, threadId);
                    } catch (...) {
                        std::unique_lock<std::mutex> lastExcepLock(lastExceptMutex);
                        lastException = std::current_exception();
                        /*
                         * This will work even when current is the largest value that
                         * size_t can fit, because fetch_add returns the previous value
                         * before the increment (what will result in overflow
                         * and produce 0 instead of current + 1).
                         */
                        current = end;
                        break;
                    }
                }
            }));
        }
        for (auto &thread : threads) {
            thread.join();
        }
        if (lastException) {
            std::rethrow_exception(lastException);
        }
    }


}



namespace ACTION {
	void simplexRegression(double *A_ptr, int A_cols, double *B_ptr, int B_rows, int B_cols, double *X_ptr);
	void AA (double *A_ptr, int A_rows, int A_cols, double *W0_ptr, int W0_cols, double *C_ptr, double *H_ptr);
	Projection orthoPCA(sp_mat &profile, int PCA_dim, int iter, int seed);

	Projection reduceGeneExpression(sp_mat &expression, int reduced_dim = DEFAULT_PCA_DIM, int method = ACTIONplusPCA, int iter = 10) {
		printf("Reducing expression matrix\n");
		Projection projection;
		
		switch(method) {
			case ACTIONplusPCA:	// Uses dense formulation
			{
				printf("\tReduce expression matrix using orthogonalization followed by PCA (k = %d) using sparse formulation ... \n", reduced_dim); fflush(stdout);
				projection = reducedKernel(expression, reduced_dim, iter, 1365);				

				printf("done\n"); fflush(stdout);
				
			}
			break;
				
			default:
				fprintf(stderr, "Unknown RNA reduction method code %d\n", method);
		}

		return projection;
	}
	
	SPA_results SPA(mat A, int k) {	
		
		SPA_results res;
		
		int n = A.n_cols;
		uvec K(k); // selected columns from A
			
	 
		rowvec normM = sum(A % A, 0); 
		rowvec normM1 = normM;
		
		mat U(A.n_rows, k);
		
		vec norm_trace = zeros(k);
		double eps = 1e-6;
		for (int i = 0; i < k; i++) {		
			// Find the column with maximum norm. In case of having more than one column with almost very small diff in norm, pick the one that originally had the largest norm
			double a = max(normM); 
			norm_trace(i) = a;
			
			uvec b = find((a*ones(1, n)-normM)/a <= eps); 
			
			if(b.n_elem > 1) {
				uword idx = index_max(normM1(b)); 
				K(i) = b(idx);
			}
			else {
				K(i) = b(0);			
			}			
			
			// Pick column
			U.col(i) = A.col(K(i));

			// Orthogonalize with respect to current basis
			for (int j = 0; j < i-1; j++) {
				U.col(i) = U.col(i) - dot(U.col(j), U.col(i)) * U.col(j);
			}
			U.col(i) = U.col(i)/ norm(U.col(i), 2); 
			
			// Update column norms
			vec u = U.col(i);            
			for (int j = i-1; 0 <= j; j--) {
				u = u - dot(U.col(j), u)*U.col(j); 
			}
			rowvec r = u.t()*A;
			normM = normM - (r % r);
		}
			
		res.selected_columns = K;
		res.column_norms = norm_trace;
		
		return res;
	}

	SPA_results DCSS(sp_mat A, int k, int dim = 5) {	
		
		field<mat> SVD_out = frSVD(A, dim, 5, 0);
		mat V = SVD_out(2);
		vec normV = sum(V % V, 1); 
		
		printf("V:%d x %d\n", V.n_rows, V.n_cols);

		
		uvec indices = stable_sort_index(normV , "descend" );
		
		
		SPA_results res;		
		res.selected_columns = indices(span(0, k-1));		
		res.column_norms = vec(normV(res.selected_columns));
		
		return res;
	}


	field<mat> AA (mat &A, mat &W0) {
		double *A_ptr = A.memptr();
		double *W0_ptr = W0.memptr();
		
		int A_rows = A.n_rows;
		int A_cols = A.n_cols;
		int W0_cols = W0.n_cols;
		
		double *C_ptr = (double *)calloc(A_cols*W0_cols, sizeof(double));
		double *H_ptr = (double *)calloc(A_cols*W0_cols, sizeof(double));

		AA(A_ptr, A_rows, A_cols, W0_ptr, W0_cols, C_ptr, H_ptr);
		
		mat C = mat(C_ptr, A_cols, W0_cols);
		mat H = mat(H_ptr, W0_cols, A_cols);

		C = clamp(C, 0, 1);
		C = normalise(C, 1);
		H = clamp(H, 0, 1);
		H = normalise(H, 1);
		
		field<mat> decomposition(2,1);
		decomposition(0) = C;
		decomposition(1) = H;
		
		return decomposition;
	}

	void simplexRegression(mat &A, mat &B, double *X_ptr) { // min(|| AX - B ||) s.t. simplex constraint
		double *A_ptr = A.memptr();
		double *B_ptr = B.memptr();
		
		int A_cols = A.n_cols;
		int B_rows = B.n_rows;
		int B_cols = B.n_cols;
		
		simplexRegression(A_ptr, A_cols, B_ptr, B_rows, B_cols, X_ptr);
	}

	mat oneHot_encoding(vec batches) {
		vec uniue_batches = sort(unique(batches));
		
		mat encoding = zeros(uniue_batches.n_elem, batches.n_elem);
		for(int i = 0; i < uniue_batches.n_elem; i++) {
			uvec idx = find(batches == uniue_batches[i]);
			vec batch_encoding = zeros(batches.n_elem);
			batch_encoding.elem(idx).ones();

			encoding.row(i) = trans(batch_encoding);
		}
		
		return(encoding);
	}
	
	ACTION_results runACTION_withBatch(mat S_r, vec batch, int k_min = 2, int k_max=20, int max_correction_rounds = 3, double lambda = 1, int numThreads = 4) {
		int feature_no = S_r.n_rows;
		
		printf("Running ACTION (with batch correction)\n");

		mat Phi = oneHot_encoding(batch);		
		mat Phi_moe = join_vert(ones(1, S_r.n_cols), Phi);
		
		k_min = std::max(k_min, 2);
		k_max = std::min(k_max, (int)S_r.n_cols);	
					
		ACTION_results trace; 
		trace.H.resize(k_max + 1);
		trace.C.resize(k_max + 1);
		trace.selected_cols.resize(k_max + 1);
		
		mat S_r_scaled = normalise(S_r, 1);
		 
		printf("Iterating from k=%d ... %d\n", k_min, k_max);
		
		field<mat> AA_res(2,1);
		for(int kk = k_min; kk <= k_max; kk++) {
			printf("\tK = %d\n", kk);
			
			printf("\t\tRunning SPA ... ");
			SPA_results SPA_res = SPA(S_r_scaled, kk);
			trace.selected_cols[kk] = SPA_res.selected_columns;
			printf("Done\n");
			
			printf("\t\tRunning AA with %d (+1) rounds of correction ...\n", max_correction_rounds);
			mat W = S_r_scaled.cols(trace.selected_cols[kk]);			
			
			// First round is just AA using raw input
			mat S_r_corr_scaled = S_r_scaled;
			for(int correction_round = 0; correction_round <= max_correction_rounds; correction_round++) {
				printf("\t\t\tRound %d\n", correction_round);
				AA_res = AA(S_r_corr_scaled, W);
				
				// Correction using mixture of experts -- Adopted from the Harmony method
				printf("\t\t\t\tCorrection ... ");
				mat H =  AA_res(1);				

				mat S_r_corr = S_r;
				for(int k = 0; k < H.n_rows; k++) {
					rowvec h = H.row(k);
					mat Phi_Rk = Phi_moe.each_row() % h;
					
					mat beta = arma::inv(Phi_Rk * Phi_moe.t() + lambda) * Phi_Rk * S_r.t();
					beta.row(0).zeros(); // do not remove the intercept 
					S_r_corr -= beta.t() * Phi_Rk;
				}
				S_r_corr_scaled = normalise(S_r_corr, 1);				
				printf("Done\n");
				
				mat C = AA_res(0);
				//mat C = normalise(trans(AA_res(1)), 1, 0);
				W = S_r_corr_scaled * C; // Start the next round of AA from current state 				
			}
						
			trace.C[kk] = AA_res(0);
			trace.H[kk] = AA_res(1);
		}

		return trace;
	}		
	

	ACTION_results runACTION(mat S_r, int k_min, int k_max, int numThreads) {
		int feature_no = S_r.n_rows;
		
		printf("Running ACTION\n");
		
		k_min = std::max(k_min, 2);
		k_max = std::min(k_max, (int)S_r.n_cols);	
					
		ACTION_results trace; 
		trace.H.resize(k_max + 1);
		trace.C.resize(k_max + 1);
		trace.selected_cols.resize(k_max + 1);
		
		mat X_r = normalise(S_r, 1);
		 
		printf("Iterating from k=%d ... %d\n", k_min, k_max);
		
		int total = 0;
		field<mat> AA_res(2,1);
		//for(int kk = k_min; kk <= k_max; kk++) {
		ParallelFor(k_min, k_max+1, numThreads, [&](size_t kk, size_t threadId) {			
			total++;
			printf("\t%d/%d\n", total, k_max-k_min+1);
			SPA_results SPA_res = SPA(X_r, kk);
			trace.selected_cols[kk] = SPA_res.selected_columns;
			
			mat W = X_r.cols(trace.selected_cols[kk]);
	//		W(span(0, 5), span(0, kk-1)).print("W");
			
			AA_res = AA(X_r, W);
			
			trace.C[kk] = AA_res(0);
			trace.H[kk] = AA_res(1);
		});

		return trace;
	}	
}
