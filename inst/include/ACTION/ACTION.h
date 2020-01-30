#ifndef ACTION_H
#define ACTION_H

#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <errno.h>
#include <sys/stat.h>
#include <getopt.h>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <unordered_map>
#include <ctime>


#define ARMA_DONT_USE_WRAPPER
#undef ARMA_BLAS_CAPITALS
#define ARMA_BLAS_UNDERSCORE

/*
#undef ARMA_BLAS_CAPITALS
#define ARMA_BLAS_UNDERSCORE
#define ARMA_64BIT_WORD
#define ARMA_BLAS_LONG_LONG
*/

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace std;

#define DEFAULT_PCA_DIM 50

// Expression reduction methods
#define PCA_only 0
#define ACTIONplusPCA 1

// Generic structure returned by all reduction methods
struct Projection {
	mat S_r;
	mat V;
	vec lambda;
	vec exp_var;
};


struct SPA_results {
	uvec selected_columns;
	vec column_norms;
};


struct ACTION_results {
	vector<uvec> selected_cols;
	vector<mat> H;
	vector<mat> C;
};

namespace ACTION {
	SPA_results SPA(mat M, int k); // Solve convex-NMF using the Successive Projection Algorithm (SPA)
	SPA_results DCSS(sp_mat A, int k, int dim); // Deterministic column subset selection
	
	void simplexRegression(mat &A, mat &B, double *X_ptr); // min_{X} (|| AX - B ||) s.t. simplex constraint
	field<mat> AA (mat &X, mat &Z0); // Robust archetypal analysis method
	ACTION_results runACTION(mat S_r, int k_min, int k_max, int numThreads); // Main ACTION function	
	ACTION_results runACTION_withBatch(mat S_r, vec batch, int k_min, int k_max, int max_correction_rounds, double lambda, int numThreads);
	
	Projection reduceGeneExpression(sp_mat &expression, int reduced_dim, int method, int iter);	

	Projection reducedKernel(sp_mat &profile, int PCA_dim, int iter, int seed);
	field<mat> frSVD(sp_mat &A, int dim, int iters, int seed);
	
}

#endif
