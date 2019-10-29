#include <actionetcore.h>


namespace ACTIONetcore {
	vec FastLasso(vec ipy, mat ipx, double lambda, size_t N){
	  double tol = 1e-3; //threshold for stopping
	  size_t p = ipy.n_elem;
	  double dbm;  //maximum of beta difference for while loop
	  vec beta = zeros(p); //initialize output beta vector with length p and filled with 0.0
	  vec gc = zeros(p); //initialize grandient components
	  
	  //update beta vector with coordinate descent, covariance updates
	  do{
		dbm = 0;
		for(size_t j = 0; j < p; j++){
		  double z = (ipy[j] - gc[j])/N + beta[j];
		  double beta_tmp = max(0.0, z - lambda) - max(0.0, -z - lambda);
		  double diffBeta = beta_tmp - beta[j];
		  double diffabs = abs(diffBeta);
		  if (diffabs > 0){
			beta[j] = beta_tmp;
			for (size_t i = 0; i < p; i++){
			  gc[i] = gc[i] + ipx(i,j) * diffBeta;
			}
			dbm = max(dbm, diffabs);
		  }
		}
	  }
	  while(dbm >= tol);
	  
	  return beta;
	}


	field<mat> FastGGM(mat x, double lambda=-1, int thread_no = 8){
	  field<mat> results(6);

	  // Stop threshold and parameter for edge decision
	  double tol = 1e-3;
	  
	  size_t n = x.n_rows;
	  size_t p = x.n_cols;
	  
	  //penalty parameter for L1 norm of betas
	  if(lambda == -1){
		lambda = sqrt(2 * log(p/sqrt(n))/n);
		cout << "Use default lambda = sqrt(2*log(p/sqrt(n))/n)" << endl;
	  }
	  cout << "In this case, lambda = " << lambda << endl;
	  
	  mat precision = zeros(p,p);
	  mat p_precision = zeros(p,p);
	  mat partialCor = zeros(p,p);
	  mat p_partialCor = zeros(p,p);
	  mat CI_low_parCor = zeros(p,p);
	  mat CI_high_parCor = zeros(p,p);
	  
	  // Center matrix
	  cout << "Center each column." << endl;
	  mat x_cent(n,p);
	  for(size_t i = 0; i < p; i++){
		x_cent.col(i) = x.col(i) - mean(x.col(i));
	  }
	  
	  // Standardize matrix
	  cout << "Standardize each column." << endl;
	  vec x_l2norm(p);
	  for(size_t i = 0; i < p; i++){
		//x_l2norm(i) = sqrt(sum(square(x_cent.col(i))));
		x_l2norm(i) = norm(x_cent.col(i), 2);
		if(x_l2norm(i) == 0){
			fprintf(stderr, "some variables have 0 variance, please remove them at first.");
			return(results);
		}
	  }
	  mat x_stan(n,p);
	  for(size_t i = 0; i < p; i++){
		x_stan.col(i) = x_cent.col(i) / x_l2norm[i] * sqrt(n);
	  }
	  
	  // Pre-calculate inner product matrixes
	  cout << "Pre-calculate inner product matrixes." << endl;
	  mat IP_YX(p,p);
	  for(size_t i = 0; i < p; i++){
		for(size_t j = 0; j < p; j++){
		  IP_YX(i,j) = dot(x_cent.col(i), x_stan.col(j));
		}
	  }
	  mat IP_XX(p,p);
	  for(size_t i = 0; i < p; i++){
		IP_XX.row(i) = IP_YX.row(i) / x_l2norm[i] * sqrt(n);
	  }
	  
	  // Pre-calculate scaled Lasso for each variable
	  cout << "Pre-calculate scaled Lasso for each variable." << endl;
	  mat beta_pre = zeros(p,p);
	  mat epsi_pre = zeros(n,p);
	  
	  #pragma omp parallel for num_threads(thread_no) 			
	  for(size_t i = 0; i < p; i++){
		double sigma = 1.0;
		vec beta(p-1);
		vec epsi(n);
		double reg, diffBeta, sigma_tmp, diffSigma;
		
		// Extract IP_XX[-i,-i] and IP_YX[i,-i]
		mat tmpxx(p-1,p-1);
		vec tmpyx(p-1);
		size_t t = 0;
		for (size_t k = 0; k < p; k++){
		  if (k != i){
			tmpyx[t] = IP_YX(i,k);
			size_t t1 = 0;
			for (size_t l = 0; l < p; l++){ // l is index of row, k is index of col
			  if (l != i){
				tmpxx(t1,t) = IP_XX(l,k);
				t1++;
			  }
			}
			t++;
		  }
		}
		// Extract x_stan[,-i]
		mat tmp_x(n, p-1);
		t=0;
		for (size_t k = 0; k < p; k++){
		  if (k != i){
			tmp_x.col(t) = x_stan.col(k);
			t++;
		  }
		}
		
		//Scaled_lasso
		size_t iter = 0;
		do {
		  //Update beta when fix sigma
		  reg = sigma * lambda;
		  vec beta_tmp = FastLasso(tmpyx, tmpxx, reg, n);
		  diffBeta = max(abs(beta_tmp - beta));
		  beta = beta_tmp;
		  //Update sigma when fix beta
		  vec tmp_y(n);
		  for (size_t k = 0; k < n; k++ ){
			tmp_y[k] = dot(tmp_x.row(k), beta);
		  } //Multiplication of x.stan[,-i] and beta
		  epsi = x_cent.col(i) - tmp_y;
		  sigma_tmp = sqrt(sum(square(epsi))/n);
		  diffSigma = abs(sigma_tmp - sigma);
		  sigma = sigma_tmp;
		  iter++;
		}
		while((diffSigma >= tol || diffBeta >= tol) && iter < 10); //Stop iteration
		epsi_pre.col(i) = epsi;
		t = 0;
		for (size_t k = 0; k < p; k++){
		  if (k != i){
			beta_pre(k,i) = beta[t];
			t++;
		  }
		}
	  }
	  
	  // Fast pairwise GGM based on the pre-calculated scaled Lasso
	  cout << "Perform pairwise GGM." << endl;
	  
	  int perc = 1;
	  int total_counts = 0;
	  #pragma omp parallel for num_threads(thread_no) 			
	  for(size_t i = 0; i < (p-1); i++){
		total_counts ++;
		if(round(100*total_counts/(p-1)) > perc) {
			printf("%d %%\n", perc);
			perc++;
		}
		for(size_t j = (i+1); j < p; j++){
		  vec epsi_i(n);
		  vec epsi_j(n);
		  
		  // Solve scaled Lasso i without variable j
		  if(beta_pre(j,i) == 0){
			epsi_i = epsi_pre.col(i);
		  } else{
			double sigma = 1.0;
			vec beta(p-2);
			vec epsi(n);
			double reg, diffBeta, sigma_tmp, diffSigma;
			
			// Extract IP_XX[-c(i,j),-c(i,j)] and IP_YX[i,-c(i,j)]
			mat tmpxx(p-2,p-2);
			vec tmpyx(p-2);
			size_t t = 0;
			for (size_t k = 0; k < p; k++){
			  if (k != i && k != j){
				tmpyx[t] = IP_YX(i,k);
				size_t t1 = 0;
				for (size_t l = 0; l < p; l++){
				  if ( l != i && l != j){
					tmpxx(t1,t) = IP_XX(l,k);
					t1++;
				  }
				}
				t++;
			  }
			}
			// Extract x_stan[,-c(i,j)]
			mat tmp_x(n, p-2);
			t=0;
			for (size_t k = 0; k < p; k++){
			  if (k != i && k != j){
				tmp_x.col(t) = x_stan.col(k);
				t++;
			  }
			}
			
			//Scaled_lasso
			size_t iter = 0;
			do {
			  //Update beta when fix sigma
			  reg = sigma * lambda;
			  vec beta_tmp = FastLasso(tmpyx, tmpxx, reg, n);
			  diffBeta = max(abs(beta_tmp - beta));
			  beta = beta_tmp;
			  //Update sigma when fix beta
			  vec tmp_y(n);
			  for (size_t k = 0; k < n; k++){
				tmp_y[k] = dot(tmp_x.row(k), beta);
			  }
			  epsi = x_cent.col(i) - tmp_y;
			  sigma_tmp = sqrt(sum(square(epsi))/n);
			  diffSigma = abs(sigma_tmp - sigma);
			  sigma = sigma_tmp;
			  iter++;
			}
			while((diffSigma >= tol || diffBeta >= tol) && iter < 10);
			epsi_i = epsi;
		  }
		  
		  // Solve scaled Lasso j without variable i
		  if(beta_pre(i,j) == 0){
			epsi_j = epsi_pre.col(j);
		  } else{
			double sigma = 1.0;
			vec beta(p-2);
			vec epsi(n);
			double reg, diffBeta, sigma_tmp, diffSigma;
			
			mat tmpxx(p-2,p-2);
			vec tmpyx(p-2);
			size_t t = 0;
			for (size_t k = 0; k < p; k++){
			  if ( k != i && k != j){
				tmpyx[t] = IP_YX(j,k);
				size_t t1 = 0;
				for (size_t l = 0; l < p; l++){
				  if ( l != i && l != j){
					tmpxx(t1,t) = IP_XX(l,k);
					t1++;
				  }
				}
				t++;
			  }
			}
			mat tmp_x(n, p-2);
			t = 0;
			for (size_t k = 0; k < p; k++){
			  if (k != i && k != j){
				tmp_x.col(t) = x_stan.col(k);
				t++;
			  }
			}
			
			//Scaled_lasso
			size_t iter = 0;
			do {
			  //Update beta when fix sigma
			  reg = sigma * lambda;
			  vec beta_tmp = FastLasso(tmpyx, tmpxx, reg, n);
			  diffBeta = max(abs(beta_tmp - beta));
			  beta = beta_tmp;
			  //Update sigma when fix beta
			  vec tmp_y(n);
			  for (size_t k = 0; k < n; k++){
				tmp_y[k] = dot(tmp_x.row(k), beta);
			  }
			  epsi = x_cent.col(j) - tmp_y;
			  sigma_tmp = sqrt(sum(square(epsi))/n);
			  diffSigma = abs(sigma_tmp - sigma);
			  sigma = sigma_tmp;
			  iter++;
			}
			while((diffSigma >= tol || diffBeta >= tol) && iter < 10);
			epsi_j = epsi;
		  }
		  
		  // Precision, solve(t(epsi_ij)%*%epsi_ij/n)
		  mat omega_tmp(2,2);
		  mat omega(2,2);
		  omega_tmp(0,0) = dot(epsi_i, epsi_i)/n;
		  omega_tmp(1,1) = dot(epsi_j, epsi_j)/n;
		  omega_tmp(0,1) = dot(epsi_i, epsi_j)/n;
		  omega_tmp(1,0) = omega_tmp(0,1);
		  // Inverse matrix of omega_tmp
		  double tmp = omega_tmp(0,0) * omega_tmp(1,1) - omega_tmp(0,1) * omega_tmp(1,0);
		  omega(0,0) = omega_tmp(1,1)/tmp;
		  omega(0,1) = -omega_tmp(0,1)/tmp;
		  omega(1,0) = -omega_tmp(1,0)/tmp;
		  omega(1,1) = omega_tmp(0,0)/tmp;
		  precision(i,i) = omega(0,0);
		  precision(j,j) = omega(1,1);
		  precision(i,j) = omega(0,1);
		  precision(j,i) = omega(1,0);
		  
		  // P-value of precision
		  double var_new = (pow(omega(0,1),2) + omega(0,0) * omega(1,1))/n;
		  double zscore = abs(omega(0,1)/sqrt(var_new));
		  p_precision(i,j) = 2 * normcdf(-zscore, 0.0, 1.0);
		  p_precision(j,i) = p_precision(i,j);
		  
		  // Partial correlation
		  partialCor(i,j) = -omega(0,1)/sqrt(omega(0,0) * omega(1,1));
		  partialCor(j,i) = partialCor(i,j);
		  
		  // P-value of partial correlation
		  var_new = pow((1-pow(partialCor(i,j), 2)), 2)/n;
		  zscore = abs(partialCor(i,j))/sqrt(var_new);
		  p_partialCor(i,j) = 2 * normcdf(-zscore, 0.0, 1.0);
		  p_partialCor(j,i) = p_partialCor(i,j);
		  
		  // 95% confidence interval of partial correlation
		  double z_95CI = 1.96; // z-score of 95% CI
		  CI_low_parCor(i,j) = partialCor(i,j) - z_95CI * (1 - pow(partialCor(i,j), 2))/sqrt(n);
		  CI_low_parCor(j,i) = CI_low_parCor(i,j);
		  CI_high_parCor(i,j) = partialCor(i,j) + z_95CI * (1 - pow(partialCor(i,j), 2))/sqrt(n);
		  CI_high_parCor(j,i) = CI_high_parCor(i,j);
		}
	  }
	  
	  results(0) = partialCor;
	  results(1) = p_partialCor;
	  results(2) = precision;
	  results(3) = p_precision;
	  results(4) = CI_low_parCor;
	  results(5) = CI_high_parCor;

	  return(results);
	}
}
