#include <actionetcore.h>

#include <hnswlib.h>
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

namespace ACTIONetcore {
	// k^{*}-Nearest Neighbors: From Global to Local (NIPS 2016)
	field<sp_mat> buildAdaptiveACTIONet(mat &H_stacked, double LC = 1.0, double M = 16, double ef_construction = 200, double ef = 10, int thread_no=8, int sym_method = ACTIONet_AND) {

		printf("Building adaptive ACTIONet (LC = %.2f)\n", LC);

		H_stacked = clamp(H_stacked, 0, 1);
		H_stacked = normalise(H_stacked, 1, 0); // make the norm (sum) of each column 1			

		double kappa = 5.0;
		int sample_no = H_stacked.n_cols;		
		int kNN = min(sample_no-1, (int)(kappa*round(sqrt(sample_no)))); // start with uniform k=sqrt(N) ["Pattern Classification" book by Duda et al.]
	
		
		int dim = H_stacked.n_rows;
		int max_elements = H_stacked.n_cols;
		hnswlib::JSDSpace* space = new hnswlib::JSDSpace(dim);		
		hnswlib::HierarchicalNSW<double> *appr_alg = new hnswlib::HierarchicalNSW<double>(space, max_elements, M, ef_construction);
		//std::unique_ptr<hnswlib::JSDSpace> space = std::unique_ptr<hnswlib::JSDSpace>(new hnswlib::JSDSpace(dim));
		//std::unique_ptr<hnswlib::HierarchicalNSW<double>> appr_alg = std::unique_ptr<hnswlib::HierarchicalNSW<double>>(new hnswlib::HierarchicalNSW<double>(space.get(), max_elements, M, ef_construction));
				
		ParallelFor(0, max_elements, thread_no, [&](size_t j, size_t threadId) {
			appr_alg->addPoint(H_stacked.colptr(j), static_cast<size_t>(j));
			
		});


		mat idx = zeros(sample_no, kNN+1);
		mat dist = zeros(sample_no, kNN+1);
//		for(int i = 0; i < sample_no; i++) {
		ParallelFor(0, sample_no, thread_no, [&](size_t i, size_t threadId) {
			
			std::priority_queue<std::pair<double, hnswlib::labeltype>> result = appr_alg->searchKnn(H_stacked.colptr(i), kNN+1);		
			
			if (result.size() != (kNN+1)) {
			  printf("Unable to find %d results. Probably ef (%f) or M (%d) is too small\n", kNN, ef, M);
			}
			
			for (size_t j = 0; j <= kNN; j++) {
				auto &result_tuple = result.top();
				dist(i, kNN-j) = result_tuple.first;
				idx(i, kNN-j) = result_tuple.second;
				
				result.pop();
			}

		});

		
/*
		for(int j = 0; j < max_elements; j++) {
			appr_alg->addPoint(H_stacked.colptr(j), static_cast<size_t>(j));
		}

		mat idx = zeros(sample_no, kNN+1);
		mat dist = zeros(sample_no, kNN+1);
		for(int i = 0; i < sample_no; i++) {
			std::priority_queue<std::pair<double, hnswlib::labeltype>> result = appr_alg->searchKnn(H_stacked.colptr(i), kNN+1);		
			
			if (result.size() != (kNN+1)) {
			  printf("Unable to find %d results. Probably ef (%f) or M (%d) is too small\n", kNN, ef, M);
			}
			
			for (size_t j = 0; j <= kNN; j++) {
				auto &result_tuple = result.top();
				dist(i, kNN-j) = result_tuple.first;
				idx(i, kNN-j) = result_tuple.second;
				
				result.pop();
			}

		}
*/
		
		delete(appr_alg);
		delete(space);

		dist = clamp(dist, 0.0, 1.0); 
		idx = clamp(idx, 0, sample_no - 1);
		
		printf("\tConstructing adaptive-nearest neighbor graph ... \n");
		mat Delta;
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

		
		
		sp_mat G(sample_no, sample_no);		
		//for(int v = 0; v < sample_no; v++) {				
		ParallelFor(0, sample_no, 1, [&](size_t v, size_t threadId) {
			vec delta = Delta.col(v);		
					
			//uvec rows = find(delta > 0, 1, "last");
			uvec rows = find(delta < 0, 1, "first");
			int neighbor_no = rows.n_elem == 0?kNN:(rows(0));
			
			int dst = v;								
			rowvec v_dist = dist.row(v);
			rowvec v_idx = idx.row(v);
			for (int i = 1; i < neighbor_no; i++) {				
				int src = v_idx(i);
				G(src, dst) = 1.0 - v_dist(i);
			}
		});
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

}






