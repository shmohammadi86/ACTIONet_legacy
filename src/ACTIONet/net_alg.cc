#include <actionetcore.h>

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

	mat PR_linsys(sp_mat &G, sp_mat &X, double alpha = 0.85, int thread_no = -1) {		
		X = normalise(X, 1, 0);
		
		/*
		rowvec d = sum(G, 0);
		uvec idx = find(c != 0);
		d[idx] = 1 / d[idx];
		
		sp_mat D;		
		D.diag() = d;
		
		sp_mat I = speye(size(G));
		*/
		sp_mat P = normalise(G, 1, 0);		
		sp_mat I = speye(size(P));
		sp_mat A = I - alpha*P;		
		//mat PR = (1-alpha)*spsolve(A, mat(X), "superlu"); 
		mat PR = (1-alpha)*spsolve(A, mat(X)); 

		
		return(PR);
	}


	mat PR_iter(sp_mat &G, sp_mat &X0, double alpha = 0.85, int max_it = 3, int thread_no = 1) {		
		
		int N = G.n_rows;
		vec z = ones(N);
		vec c = vec(trans(sum(G, 0)));
		uvec idx = find(c);
		z(idx) = ones(idx.n_elem)*(1.0 - alpha);
		z = z / N;
		
		sp_mat P = alpha*normalise(G, 1, 0);
		X0 = normalise(X0, 1, 0);
		mat X = mat(X0);
				
		X0 *= N;			
		rowvec zt = trans(z);
		ParallelFor(0, X.n_cols, thread_no, [&](size_t i, size_t threadId) {			
			X.col(i) = P*X.col(i) + X0.col(i)*(zt*X.col(i));
		});
		//X = normalise(X, 1)
		
		return(X);
	}	
	
	
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
}





