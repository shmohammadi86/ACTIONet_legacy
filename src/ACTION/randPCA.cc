#include "ACTION.h"

namespace ACTION {
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
	
	void gram_schmidt(mat& A) {
		for(uword i = 0; i < A.n_cols; ++i) {
			for(uword j = 0; j < i; ++j) {
				double r = dot(A.col(i), A.col(j));
				A.col(i) -= r * A.col(j);
			}
			
			double col_norm = norm(A.col(i), 2);
			
			if(col_norm < 1E-4) {
				for(uword k = i; k < A.n_cols; ++k)
					A.col(k).zeros();

				return;
			}
			A.col(i) /= col_norm;
		}
	}

	mat applyA(const sp_mat &At, const mat &X) {
		
		mat Xt = X.t();
		mat prod = Xt*At;	
		
		return prod.t();
	}

	mat applyAt(const sp_mat &At, const mat &X) {
				
		mat prod = At*X;
		
		return prod;
	}


	mat applyA_c(const sp_mat &At, const mat &X, const vec &ct) {
	
		mat Xt = X.t();

		mat prod = Xt*At;	
		vec v = vec(Xt*ct);	

		prod.each_col() -= v;
		
		return prod.t();
	}

	mat applyAt_c(const sp_mat &At, const mat &X, const vec &ct) {
				
		mat prod = At*X;
		
		rowvec v = rowvec(sum(X, 0));		
		prod -= ct*v;
		
		return prod;
	}


	field<mat> PCA(sp_mat &A, int dim, int max_it, int seed) {	
		
		// Input matrix is gene/peak x cell. We want to reduce the gene dimension, so we transpose to perform column-wise operations	
		int n = A.n_rows;
		int m = A.n_cols;
		int l = dim + 2;		
					
		vec ct = vec(mean(A, 1));
		
		
		vec s;
		mat R, Q;
		mat U, V, X;
		printf("\t\tRunning randomized PCA. Matrix size: %d x %d\n", A.n_rows, A.n_cols); fflush(stdout);
		
		
		printf("\t\tInitializing random core matrix\n");
		if (m < n) {
	//		R = (2*randu(l, m))-1;
			//R = stats::runif<arma::mat>(l, m, -1.0, 1.0, 0);
			R = sampleUnif(l, m, -1.0, 1.0, 0);
			printf("\t\t\tR: %d x %d\n", R.n_rows, R.n_cols); fflush(stdout);

			/*
			for(int i = 0; i < l; i++) {
				char str[1024];
				sprintf(str, "col %d", i);
				rowvec u = R.col(i).t();
				u(span(l-10, l-1)).print(str);
			}

			rowvec u = R.col(R.n_cols-1).t();
			u(span(l-10, l-1)).print("last col");
			*/
			
			mat Rt = R.t();
			
			Q = applyAt_c(A, Rt, ct); // Apply the adjoint of A to a random matrix, obtaining Q.
		}
		else {
			//R = (2*randu(n,l))-1;
			//R = stats::runif<arma::mat>(n, l, -1.0, 1.0, seed);
			R = sampleUnif(n, l, -1.0, 1.0, 0);

						
			printf("\t\t\tR: %d x %d\n", R.n_rows, R.n_cols); fflush(stdout);
			Q = applyA_c(A, R, ct); // Apply A to a random matrix, obtaining Q.
		}					
		
		// Form a matrix Q whose columns constitute a well-conditioned basis for the columns of the earlier Q.			
		gram_schmidt(Q);
		//Q = orth(Q);
		
		
		if (m < n) {
			// Conduct normalized power iterations.
			for(int it = 1; it <= max_it; it++) {
				printf("\t\tIteration %d\n", it);			

				Q = applyA_c(A, Q, ct); // Apply A to a random matrix, obtaining Q.				
				gram_schmidt(Q); 
				
				Q = applyAt_c(A, Q, ct); // Apply the adjoint of A to a random matrix, obtaining Q.
				gram_schmidt(Q);
				//Q = orth(Q);					

			}

			X = mat(applyA_c(A, Q, ct));
			printf("\t\tReduced SVD ... ");
			
			svd_econ( U, s, V,  X);
			printf("done\n");
			V = Q*V;		
		}				
		else {
			// Conduct normalized power iterations.
			for(int it = 1; it <= max_it; it++) {
				printf("\t\tIteration %d\n", it);
				
				Q = applyAt_c(A, Q, ct); // Apply the adjoint of A to a random matrix, obtaining Q.
				gram_schmidt(Q);
				//Q = orth(Q);
									
				Q = applyA_c(A, Q, ct); // Apply A to a random matrix, obtaining Q.
				gram_schmidt(Q);			
				//Q = orth(Q);
				
			}
			
			// SVD Q' applied to the centered A to obtain approximations to the singular values and right singular vectors of the A;
			
			X = mat(trans(applyAt_c(A, Q, ct)));
			
			printf("\t\tReduced SVD ... ");
			svd_econ( U, s, V,  X);
			printf("done\n");
			U = Q*U;
		}

		U.shed_cols(dim, dim+1);
		s = s(span(0, dim-1));
		V.shed_cols(dim, dim+1);

		vec sigma_sq = square(s);
		vec lambda = sigma_sq / (m-1);
		vec explained_var = cumsum(sigma_sq) / sum(sigma_sq);

		mat scores(n, dim);
		scores = U;
		for (int i = 0; i < dim; i++) {
			scores.col(i) = s(i)*scores.col(i);
		}
		
		field<mat> results(4);
		results(0) = scores.t();;
		results(1) = V;
		results(2) = lambda;
		results(3) = explained_var;
		
		/*
		results(0) = scores.t();
		results(1) = V;
		results(2) = lambda;
		results(3) = explained_var;
		*/
		return results;	
	}



	field<mat> SVD(sp_mat &A, int dim, int max_it, int seed) {	
		field<mat> results(3);

		int m = A.n_rows;
		int n = A.n_cols;
		int l = dim + 2;		
				
		vec s;
		mat R, Q;
		mat U, V, X;
		
		printf("\t\tRunning randomized SVD. Matrix size: %d x %d\n", A.n_rows, A.n_cols); fflush(stdout);
		
		if (m < n) {
			//R = stats::runif<arma::mat>(l, m, -1.0, 1.0, seed);
			R = sampleUnif(l, m, -1.0, 1.0, 0);

			sp_mat At = A.t();
			Q = At*R.t(); 
		}
		else {
			//R = stats::runif<arma::mat>(n, l, -1.0, 1.0, seed);
			R = sampleUnif(n, l, -1.0, 1.0, 0);
			Q = A*R;
		}				


		// Form a matrix Q whose columns constitute a well-conditioned basis for the columns of the earlier Q.			
		gram_schmidt(Q);
		//Q = orth(Q);
		
		if (m < n) {
			// Conduct normalized power iterations.
			for(int it = 1; it <= max_it; it++) {
				printf("\t\tIteration %d\n", it);
				
				
				Q = A*Q; 
				gram_schmidt(Q);
				//Q = orth(Q);

				Q = A.t()*Q; 
				gram_schmidt(Q);								
				//Q = orth(Q);
			}

			X = mat(A*Q);
			printf("\t\tReduced SVD ... ");
			svd_econ( U, s, V,  X);
			printf("done\n");
			V = Q*V;		
		}				
		else {
			// Conduct normalized power iterations.
			for(int it = 1; it <= max_it; it++) {
				printf("\t\tIteration %d\n", it);
				
				Q = A.t()*Q;
				gram_schmidt(Q);
				//Q = orth(Q);
					
				Q = A*Q; // Apply A to a random matrix, obtaining Q.
				gram_schmidt(Q);			
				//Q = orth(Q);
			}
			
			// SVD Q' applied to the centered A to obtain approximations to the singular values and right singular vectors of the A;
			
			X = mat(Q.t()*A);
			svd_econ( U, s, V,  X);
			U = Q*U;		
		}

		U.shed_cols(dim, dim+1);
		s = s(span(0, dim-1));
		V.shed_cols(dim, dim+1);
		
		results(0) = U;
		results(1) = s;
		results(2) = V;
		
		return results;	
	}


	Projection orthoPCA(sp_mat &profile, int PCA_dim, int iter, int seed) {			
		int n = profile.n_rows;
		//profile = normalise(profile, 2);    

		printf("\tRunning ACTION+PCA. Matrix size: %d x %d\n", profile.n_rows, profile.n_cols); fflush(stdout);
		
		// Update 1: Orthogonalize columns w.r.t. background (mean)
		vec mu = vec(mean(profile, 1));
		vec v = mu / norm(mu, 2);
		vec a1 = v;
		vec b1 = -trans(profile)*v;
		
		// Update 2: Center columns of orthogonalized matrix before performing SVD (i.e., PCA)
		vec c = vec(trans(mean(profile, 0)));
		double a1_mean = mean(a1);
		vec a2 = ones(profile.n_rows);
		vec b2 = -(a1_mean*b1 + c);

		mat A = join_rows(a1, a2);
		mat B = join_rows(b1, b2);
		
		printf("\tPerform SVD on the original matrix\n"); fflush(stdout);
		field<mat> SVD_results = SVD(profile, PCA_dim, iter, 1365);	

		mat U = SVD_results(0);
		vec s = SVD_results(1);
		mat V = SVD_results(2);

		printf("\tUpdate SVD ..."); fflush(stdout);
		vec s_prime;
		mat U_prime, V_prime;
			
		mat M = U.t()*A; 
		mat A_ortho_proj = A - U*M;   
		mat P = A_ortho_proj;// = orth(A_ortho_proj);
		gram_schmidt(P);		
		mat R_P = P.t()*A_ortho_proj;
		
		
		mat N = V.t()*B; 
		mat B_ortho_proj = B - V*N; 
		mat Q = B_ortho_proj; //orth(B_ortho_proj); 
		gram_schmidt(Q);
		mat R_Q = Q.t()*B_ortho_proj;	
		
		mat K1 = zeros(s.n_elem+2, s.n_elem+2);
		for(int i = 0; i < s.n_elem; i++) {
			K1(i, i) = s(i);
		}

		mat K2 = join_vert(M, R_P)*trans(join_vert(N, R_Q));

		
		mat K = K1 + K2;

		svd( U_prime, s_prime, V_prime, K );
		
		mat U_updated = join_horiz(U, P)*U_prime;
		mat V_updated = join_horiz(V, Q)*V_prime;
		printf("done.\n"); fflush(stdout);
		
		
		Projection projection;		
		projection.S_r = trans(V_updated.cols(0, PCA_dim-1));
		for(int i = 0; i < PCA_dim; i++) {
			projection.S_r.row(i) *= s_prime(i);
		}
		projection.V = U_updated.cols(0, PCA_dim-1);	
		
		vec sigma_sq = square(s_prime);

		projection.lambda = sigma_sq / (n-1);

		vec explained_var = cumsum(sigma_sq) / sum(sigma_sq);
		projection.exp_var = explained_var(span(0, PCA_dim-1));
		
		return projection;	
	}
}
