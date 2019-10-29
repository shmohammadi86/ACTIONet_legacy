#include <iostream>
#include <cstdio>
#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/connected_components.hpp>
// #include <boost/unordered_set.hpp>
#include <ctime>
#include <stack>
#include <queue>
#include <fstream>
// #include <boost/range/algorithm/set_algorithm.hpp>
#include <boost/functional/hash.hpp>
#include <random>
#include <unordered_set>
#include <unordered_map>
// using namespace std;


#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace boost;

typedef std::vector < std::vector <int> > vvi;

typedef adjacency_list< listS, vecS, undirectedS, no_property, property<edge_weight_t, float> > Graph;
// typedef adjacency_list< listS,vecS, undirectedS > Graph;
typedef graph_traits < Graph> ::vertex_descriptor vertex_descriptor;
typedef graph_traits < Graph> ::edge_descriptor edge_descriptor;
typedef Graph::vertex_iterator vertex_iterator;
typedef Graph::edge_iterator edge_iterator;
typedef std::pair <int, int> Edge;
typedef graph_traits < Graph> ::adjacency_iterator my_adjacency_iterator;
property_map<Graph, edge_weight_t>::type weight;


std::unordered_set<int> union_(std::unordered_set<int> &set1, std::unordered_set<int> &set2)
{
	// returns the union of two unordered sets of integers
	std::unordered_set<int> new_set = set1;

	for(auto thing : set2)
		new_set.insert(thing);
	return new_set;
}

std::vector< std::vector <int> > sample(std::vector <int> &things, int &k)
{
	// Returns the set of sampled nodes needed for the whole process of spanner construction as vector of integer vectors 
	std::vector< std::vector<int> > master_v;
	std::vector<int> v = things;

	int i = 1, n = v.size(), s;
	float f;

	std::vector<int> counts;

	for(i = 1; i < k; i ++)
	{
		f = i / (float) k;
		s = ceil(pow(n, 1 - f));
		counts.push_back(s);
	}

	master_v.push_back(v);

	std::default_random_engine engine(std::random_device {}());
	std::shuffle(master_v[0].begin(), master_v[0].end(), engine);

	for(i = 1; i < k; i ++)
	{
		std::vector<int> temp;
		for(int j = 0; j < counts[i - 1]; j ++)
			temp.push_back(master_v[i - 1][j]);
		std::shuffle(temp.begin(), temp.end(), engine);
		master_v.push_back(temp);
	}
	return master_v;
}


tuple<Graph, Graph, Graph> read_graphs(arma::sp_mat &G_adj) {
	Graph G, G_comm, G_prime;
	int u, v;
	float w;

	arma::sp_mat::iterator it     = G_adj.begin();
	arma::sp_mat::iterator it_end = G_adj.end();

	for(; it != it_end; ++it) {
		u = it.row();
		v = it.col();
		w = (*it);
		
		add_edge(u, v, w, G);
		add_edge(u, v, w, G_comm);
		add_edge(u, v, w, G_prime);
	}
	
	return make_tuple(G, G_comm, G_prime);
}

namespace ACTIONetcore {
	// "A Simple and Linear Time Randomized Algorithm for Computing Sparse Spanners in Weighted Graphs"
	arma::sp_mat make_spanner(arma::sp_mat &G_adj, int k) {
		
		// constructs the (2k - 1)-spanner using Baswana and Sen's randomized algorithm 
		Graph G, G_comm, G_prime;
		Graph G_spanner;

		tie(G, G_comm, G_prime) = read_graphs(G_adj);

		std::unordered_set <std::pair<int, int>, hash< std::pair<int, int> > > old_Ei;
		std::unordered_set <std::pair<int, int>, hash< std::pair<int, int> > > new_Ei;

		std::unordered_set <int> v_prime_flag;
		std::unordered_map <int, int> membership_oldc;
		std::unordered_map <int, int> membership_newc;
		std::unordered_map <int, int> neighborhood;

		std::unordered_map <int, std::unordered_set <int> > old_c;

		std::pair <vertex_iterator, vertex_iterator> vp;

		std::vector <int> all_nodes;

		for(vp = vertices(G); vp.first != vp.second; ++vp.first) // iterating over all the nodes
		{
			auto node = *vp.first;
			old_c[node].insert(node);
			v_prime_flag.insert(node);
			membership_oldc[node] = node;
			membership_newc[node] = -1;
			all_nodes.push_back(node);
		}

		my_adjacency_iterator start, end;

		std::vector < std::vector<int> > R_vector; //this stores all the randomly sampled cluster heads for all iterations

		R_vector = sample(all_nodes, k); // R_vector stores the set of sampled nodes required at each iteration

		for(int i = 1; i < k; i ++)
		{
			/*****1. Forming a sample of clusters **********************/

			std::unordered_set<int> R_i(R_vector[i].begin(), R_vector[i].end()); // R_vector[i] gives the randomly picked cluster heads for ith iteration
			std::unordered_map <int, std::unordered_set <int> > new_c;

			for(auto item : old_c)
			{
				int v = item.first;

				if(R_i.find(v) != R_i.end())
					new_c[v] = old_c[v];
			}

			std::unordered_set <int> sampled_nodes;

			for(auto v : R_i)
				sampled_nodes = union_(sampled_nodes, old_c[v]);

			std::unordered_set <int> unsampled_nodes;


			for(vp = vertices(G_prime); vp.first != vp.second; ++vp.first)
			{
				int node = *(vp.first);
				if(sampled_nodes.find(node) == sampled_nodes.end())
					unsampled_nodes.insert(node);
			}

			neighborhood.clear();
			for(auto node : unsampled_nodes)
				neighborhood[node] = -1;

			if(i == 1)
				for(auto v : R_i)
					membership_newc[v] = v;

			else
			{
				std::pair <vertex_iterator, vertex_iterator> vp;
				for(vp = vertices(G_prime); vp.first != vp.second; ++ vp.first)
				{
					auto v = *(vp.first);
					if(membership_newc[v] != -1)
						if(R_i.find(membership_newc[v]) == R_i.end())
							membership_newc[v] = -1;
				}
			}

			new_Ei = old_Ei;
			int u, v;

			for(auto edge : old_Ei)
			{
				tie(u, v) = edge;
				if(R_i.find(membership_newc[u]) == R_i.end() || R_i.find(membership_newc[v]) == R_i.end())
				{
					std::pair<int, int> e = std::make_pair(u, v);
					new_Ei.erase(e);
				}
			}

			/***********************end of step 1**************************/

			/********2. finding nearest neighboring sampled cluster*********/

			for(int v : unsampled_nodes)
			{
				if(v_prime_flag.find(v) == v_prime_flag.end()) // if v is not in G_prime, we don't consider it
					continue;

				std::unordered_map <int, std::pair<float, int> > min_table;
				edge_descriptor e;
				vertex_descriptor src, dest;

				for(std::tie(start, end) = adjacent_vertices(v, G_prime); start != end; ++ start)
				{
					int neighbor = *start;
					float wt;
					src = vertex(v, G);
					dest = vertex(neighbor, G);

					std::tie(e, std::ignore) = edge(src, dest, G);

					wt = get(weight, e);

					if(membership_newc[neighbor] != -1) // checking if the neighbor is sampled
					{
						if(min_table.find(membership_newc[neighbor]) == min_table.end())
							min_table[membership_newc[neighbor]] = std::make_pair(wt, neighbor);


						else if(wt < min_table[membership_newc[neighbor]].first)
						{
							min_table[membership_newc[neighbor]].first = wt;
							min_table[membership_newc[neighbor]].second = neighbor;
						}
					}
				}

				float min_ = INT_MAX;
				int nearest_neighbor = -1;

				for(auto item : min_table)
				{
					if(item.second.first  < min_)
					{
						min_ = item.second.first;
						nearest_neighbor = item.second.second;
					}
				}

				neighborhood[v] = nearest_neighbor;
			}

			/************end of step 2 ********************************/


			/* **********3. adding edges to spanner ******************/

			for(int v : unsampled_nodes)
			{
				if(v_prime_flag.find(v) == v_prime_flag.end())
					continue;

				if(neighborhood[v] == -1) //v is not adjacent to any sampled nodes
				{
					std::unordered_map <int, std::pair<float, int> > min_table;

					edge_descriptor e;
					vertex_descriptor src, dest;
					std::vector <std::pair<int, int> > edges_to_be_removed;
					my_adjacency_iterator start, end;


					for(std::tie(start, end) = adjacent_vertices(v, G_prime); start != end; ++ start)
					{
						int neighbor = *start;
						float wt;

						src = vertex(v, G);
						dest = vertex(neighbor, G);

						std::tie(e, std::ignore) = edge(src, dest, G);

						wt = get(weight, e);

						if(membership_oldc[neighbor] != -1)
						{
							if(min_table.find(neighbor) == min_table.end())
								min_table[neighbor] = std::make_pair(wt, neighbor);

							if(wt <= min_table[neighbor].first)
							{
								min_table[neighbor].first = wt;
								min_table[neighbor].second = neighbor;

								edges_to_be_removed.push_back(std::make_pair(v, neighbor));
							}
						}
					}

					for(auto edge : edges_to_be_removed)
						remove_edge(edge.first, edge.second, G_prime);

					for(auto item : min_table)
					{
						add_edge(v, item.second.second, item.second.first, G_spanner);
						remove_edge(v, item.second.second, G_comm);
					}
				}


				else
				{
					my_adjacency_iterator start, end;
					float wt;

					edge_descriptor e;
					std::tie(e, std::ignore) = edge(v, neighborhood[v], G);
					wt = get(weight, e);

					add_edge(v, neighborhood[v], wt, G_spanner);
					remove_edge(v, neighborhood[v], G_comm);

					new_Ei.insert(std::make_pair(v, neighborhood[v]));

					std::vector <int> neighbors;
					for(std::tie(start, end) = adjacent_vertices(v, G_prime); start != end; ++ start)
					{
						int neighbor = *start;
						neighbors.push_back(neighbor);
					}

					for(auto neighbor : neighbors)
					{
						if(membership_newc[neighbor] == membership_newc[neighborhood[v]])
							remove_edge(v, neighbor, G_prime);
					}

					std::unordered_map <int, std::pair<float, int> > min_table;
					edge_descriptor e1;
					vertex_descriptor src, dest;

					for(auto neighbor : neighbors)
					{
						float wt1, wt2;
						src = vertex(v, G);
						dest = vertex(neighbor, G);

						std::tie(e1, std::ignore) = edge(src, dest, G);
						wt1 = get(weight, e1);

						dest = vertex(neighborhood[v], G);
						std::tie(e1, std::ignore) = edge(src, dest, G);
						wt2 = get(weight, e1);

						if(wt1 < wt2)
						{
							if(min_table.find(membership_oldc[neighbor]) == min_table.end())
								min_table[membership_oldc[neighbor]] = std::make_pair(wt1, neighbor);

							else if(wt1 < min_table[membership_oldc[neighbor]].first)
							{
								min_table[membership_oldc[neighbor]].first = wt2;
								min_table[membership_oldc[neighbor]].second = neighbor;
								remove_edge(v, neighbor, G_prime);
							}
						}
					}

					for(auto item : min_table)
					{
						add_edge(v, item.second.second, item.second.first, G_spanner);
						remove_edge(v, item.second.second, G_comm);
						remove_edge(v, item.second.second, G_prime);
					}


				}

			}


			/************end of step 3*****************************************/


			/******4. removing intra cluster edges ***********************/

			std::vector< std::pair <int, int> > intracluster_edges;

			for(auto item : new_Ei)
			{
				int u = item.first;
				int v = item.second;

				if(membership_newc[u] != -1)
				{
					membership_newc[v] = membership_newc[u];
					new_c[membership_newc[u]].insert(v);
				}

				else if(membership_newc[v] != -1)
				{
					membership_newc[u] = membership_newc[v];
					new_c[membership_newc[v]].insert(u);
				}
			}

			edge_iterator e_start, e_end;

			for(std::tie(e_start, e_end) = edges(G_prime); e_start != e_end; ++ e_start)
			{
				int u = source(*e_start, G_prime), v = target(*e_start, G_prime);

				if(membership_newc[u] == membership_newc[v])
					intracluster_edges.push_back(std::make_pair(u, v));
			}



			for(auto edge : intracluster_edges)
			{
				// // std:: cout << "\n(" << edge.first << ", " << edge.second << ")";
				remove_edge(edge.first, edge.second, G_prime);
			}

			/*************end of step 4*************************************************/

			/**** updations needed at the end of iteration******************************/

			v_prime_flag.clear();

			int count = 0;
			// // std:: cout << "E_prime: \n";
			for(std::tie(e_start, e_end) = edges(G_prime); e_start != e_end; ++ e_start)
			{
				int u = source(*e_start, G_prime), v = target(*e_start, G_prime);
				v_prime_flag.insert(u);
				v_prime_flag.insert(v);
				count ++;
			}

			for(auto item : new_Ei)
			{
				v_prime_flag.insert(item.first);
				v_prime_flag.insert(item.second);
			}

			old_c = new_c;
			membership_oldc = membership_newc;

			/***********The end of iteration i *****************************/
		}

		/******* Phase 2 **************************/

		vertex_iterator v_start, v_end;

		for(std::tie(v_start, v_end) = vertices(G_prime); v_start != v_end; ++ v_start)
		{
			int v = *v_start;
			if(v_prime_flag.find(v) == v_prime_flag.end())
				continue;

			std::unordered_map <int, std::pair<float, int> > min_table;
			edge_descriptor e;
			vertex_descriptor src, dest;

			std::vector <std::pair<int, int> > edges_to_be_removed;


			for(std::tie(start, end) = adjacent_vertices(v, G_prime); start != end; ++ start)
			{
				int neighbor = *start;
				float wt;
				src = vertex(v, G);
				dest = vertex(neighbor, G);

				std::tie(e, std::ignore) = edge(src, dest, G);

				wt = get(weight, e);

				if(membership_newc[neighbor] == -1)
					continue;

				if(min_table.find(membership_newc[neighbor]) == min_table.end())
					min_table[membership_newc[neighbor]] = std::make_pair(wt, neighbor);

				if(wt <= min_table[membership_newc[neighbor]].first)
				{
					min_table[membership_newc[neighbor]].first = wt;
					min_table[membership_newc[neighbor]].second = neighbor;

					edges_to_be_removed.push_back(std::make_pair(v, neighbor));
				}

				if(min_table.find(membership_newc[neighbor]) == min_table.end())
					edges_to_be_removed.push_back(std::make_pair(v, neighbor));
			}

			for(auto edge : edges_to_be_removed)
				remove_edge(edge.first, edge.second, G_prime);

			for(auto item : min_table)
			{
				add_edge(v, item.second.second, item.second.first, G_spanner);
				remove_edge(v, item.second.second, G_comm);
			}
		}
		

		arma::sp_mat G_sparse(G_adj.n_rows, G_adj.n_cols);
		
		typedef property_map<Graph, vertex_index_t>::type IndexMap;
		IndexMap index = get(vertex_index, G_comm);
		
		graph_traits<Graph>::edge_iterator ei, ei_end;
		for (tie(ei, ei_end) = edges(G_spanner); ei != ei_end; ++ei) {
			int u = index[source(*ei, G_spanner)];
			int v = index[target(*ei, G_spanner)];
			
			G_sparse(u, v) = G_sparse(v, u) = G_adj(u, v);        
		}

		return(G_sparse);
	}
}


