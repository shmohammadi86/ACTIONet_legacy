#include <ACTIONet.h>
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

#define ARMA_USE_CXX11_RNG

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List reduceGeneExpression(sp_mat expression, int reduced_dim = 30, int method = 1, int iters = 10) {
	arma_rng::set_seed(1365);	
	
	Projection projection = ACTION::reduceGeneExpression(expression, reduced_dim, method, iters);
	
	List res;	
	res["S_r"] = projection.S_r;		
	res["V"] = projection.V;
	res["lambda"] = projection.lambda;
	res["explained_var"] = projection.exp_var;	
		
	return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat runsimplexRegression(mat A, mat B) {	

	mat X(A.n_cols, B.n_cols);
	ACTION::simplexRegression(A, B, X.memptr());

	X.transform( [](double val) { return (min(1.0, max(0.0, val))); } );
	X = normalise(X, 1);		
	
	return X;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List runSPA(mat A, int k) {	

	SPA_results res = ACTION::SPA(A, k);
	uvec selected_columns = res.selected_columns;
	
	vec cols(k);
	for(int i = 0; i < k; i++) {
		cols[i] = selected_columns[i] + 1;
	}
	

	List out;	
	out["selected_columns"] = cols;		
	out["norms"] = res.column_norms;
		
	return out;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List runDCSS(sp_mat &A, int k, int dim = 10) {	

	SPA_results res = ACTION::DCSS(A, k, dim);
	uvec selected_columns = res.selected_columns;
	
	vec cols(k);
	for(int i = 0; i < k; i++) {
		cols[i] = selected_columns[i] + 1;
	}
			
	List out;	
	out["selected_columns"] = cols;		
	out["norms"] = res.column_norms;
		
	return out;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List runAA(mat A, mat W0) {	

	field<mat> decomposition = ACTION::AA (A, W0);
	mat C = decomposition(0);
	mat H = decomposition(1);
	
	List out_list;		
	out_list["C"] = C;
	out_list["H"] = H;
			
	return out_list;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List runACTION(mat S_r, int k_min = 2, int k_max=20, int thread_no = 4) {	

	ACTION_results trace = ACTION::runACTION(S_r, k_min, k_max, thread_no);

	List res;
	
	List C(k_max);
	for (int i = k_min; i <= k_max; i++) {
		C[i-1] = trace.C[i];
	}
	res["C"] = C;	

	List H(k_max);
	for (int i = k_min; i <= k_max; i++) {
		H[i-1] = trace.H[i];
	}
	res["H"] = H;
	
		
	return res;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List runACTION_withBatch(mat S_r, vec batch, int k_min = 2, int k_max=20, int max_correction_rounds = 3, double lambda = 1, int numThreads = 4) {

	ACTION_results trace = ACTION::runACTION_withBatch(S_r, batch, k_min, k_max, max_correction_rounds, lambda, numThreads);
	
	List res;
	
	List C(k_max);
	for (int i = k_min; i <= k_max; i++) {
		C[i-1] = trace.C[i];
	}
	res["C"] = C;	

	List H(k_max);
	for (int i = k_min; i <= k_max; i++) {
		H[i-1] = trace.H[i];
	}
	res["H"] = H;
	
		
	return res;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List frSVD(sp_mat &A, int dim, int iters = 5, int seed = 0) {	

	field<mat> SVD_out = ACTION::frSVD(A, dim, iters, seed);
	
	List res;
	res["U"] = SVD_out(0);	
	res["sigma"] = SVD_out(1);	
	res["V"] = SVD_out(2);	
		
	return res;
}
