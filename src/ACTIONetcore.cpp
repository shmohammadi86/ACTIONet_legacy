#define ARMA_USE_CXX11_RNG
#define ARMA_64BIT_WORD

#include <actionetcore.h>
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// set seed
// [[Rcpp::export]]
void set_seed(double seed) {
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(std::floor(std::fabs(seed)));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]

List buildAdaptiveACTIONet(mat &H_stacked, double LC = 1.0, double M = 16, double ef_construction = 200, double ef = 10, int thread_no=4, string sym_method = "OR") {
	int sm = sym_method == "AND"? ACTIONet_AND:ACTIONet_OR;	
	printf("Sym method: %s (%d)\n", sym_method.c_str(), sm);
		
    field<sp_mat> res = ACTIONetcore::buildAdaptiveACTIONet(H_stacked, LC, M, ef_construction, ef, thread_no, sm);

	List out_list;		
	out_list["ACTIONet"] = res(0);
	out_list["ACTIONet_asym"] = res(1);

    return out_list;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List layoutACTIONet(sp_mat G,
	mat S_r,
	int compactness_level= 50,
	unsigned int n_epochs = 100,
	int thread_no = 8) {

	field<mat> res = ACTIONetcore::layoutACTIONet(G, S_r, compactness_level, n_epochs, thread_no);
    
	List out_list;		
	out_list["coordinates"] = res(0);
	out_list["coordinates_3D"] = res(1);
	out_list["colors"] = res(2);

    return out_list;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List reconstructArchetypes(sp_mat &S, const List& C_trace, const List& H_trace, double z_threshold = 1.0) {
	
	int n_list = H_trace.size();
	vector<mat> C_trace_vec(n_list+1);
	vector<mat> H_trace_vec(n_list+1);
	for (int i = 0; i < n_list; i++) {
		if(Rf_isNull(H_trace[i])) {
			continue;
		}
		C_trace_vec[i+1] = (as<mat>(C_trace[i]));
		H_trace_vec[i+1] = (as<mat>(H_trace[i]));
	}

	ACTIONetcore::multilevel_archetypal_decomposition results = ACTIONetcore::reconstructArchetypes(S, C_trace_vec, H_trace_vec, z_threshold = -1.0);

	List out_list;		
	
	
	for(int i = 0; i < results.selected_archs.n_elem; i++) results.selected_archs[i]++;
	out_list["selected_archs"] = results.selected_archs; 

	for(int i = 0; i < results.landmark_cells.n_elem; i++) results.landmark_cells[i]++;
	out_list["landmark_cells"] = results.landmark_cells; 

	out_list["C_stacked"] = results.C_stacked;
	out_list["H_stacked"] = results.H_stacked;
	
	out_list["archetype_profile"] = results.archetype_profile;

	out_list["backbone"] = results.backbone;

    return out_list;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat assessFeatureSets(sp_mat &S, List index_sets, int rand_perm = 100) {
	
	int n_list = index_sets.size();
	field<uvec> feature_sets(n_list);
	
	for (int i = 0; i < n_list; i++) {
		vec v = as<vec>(index_sets[i]);
		uvec u(v.n_elem);
		for(int j = 0; j < v.n_elem; j++) {
			u(j) = v(j) - 1;
		}
		feature_sets(i) = u;
	}
	
	mat scores = ACTIONetcore::assessFeatureSets(S, feature_sets, rand_perm);		

    return scores;
}




// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat phenotypeEnrichment(mat &H_stacked, mat &phenotype_associations, int rand_perm_no) {	
	mat Z = ACTIONetcore::phenotypeEnrichment (H_stacked, phenotype_associations, rand_perm_no);

    return Z;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat update_layout_2D(mat coors,
	int compactness_level= 50,
	unsigned int n_epochs = 100,
	int thread_no = 8) {
	
	if(coors.n_cols == 2)
		coors = trans(coors);
	
		
	mat updated_coors = ACTIONetcore::update_layout_2D(coors, compactness_level, n_epochs, thread_no);

    return updated_coors;
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat PageRank_iter(sp_mat &G, sp_mat &X0, double alpha = 0.85, int max_it = 3, int thread_no = 8) {	
	mat PR = ACTIONetcore::PR_iter (G, X0, alpha, max_it, thread_no);

    return PR;
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat MWM(mat &G) {	
	mat G_matched = ACTIONetcore::MWM(G);

    return G_matched;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
vec sweepcut(sp_mat &A, vec s) {
    vec conductance = ACTIONetcore::sweepcut(A, s);

    return conductance;
}

