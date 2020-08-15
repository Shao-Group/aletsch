/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "scallop.h"
#include "constants.h"
#include "essential.h"

#include <cstdio>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <climits>
#include <cfloat>
#include <algorithm>

scallop::scallop(splice_graph &g, hyper_set &h, const parameters &c)
	: gr(g), hs(h), cfg(c)
{
	round = 0;
	//gr.draw(gr.gid + "." + tostring(round++) + ".tex");
	gr.get_edge_indices(i2e, e2i);
	//add_pseudo_hyper_edges();
	//hs.build_confident_nodes(cfg.min_confident_phasing_path_weight);
	hs.build(gr, e2i);
	init_super_edges();
	init_vertex_map();
	init_inner_weights();
	init_nonzeroset();
}

scallop::~scallop()
{
}

int scallop::assemble()
{
	int c = classify();
	if(cfg.verbose >= 2) printf("process splice graph %s type = %d, vertices = %lu, edges = %lu, phasing paths = %lu\n", gr.gid.c_str(), c, gr.num_vertices(), gr.num_edges(), hs.edges.size());

	//resolve_negligible_edges(false, cfg.max_decompose_error_ratio[NEGLIGIBLE_EDGE]);

	while(true)
	{	
		if(gr.num_vertices() > cfg.max_num_exons) break;

		/*
		printf("---------\n");
		gr.print();

		printf("---------\n");
		print_super_edges();
		//printf("---------\n");
		//hs.print_edges();
		//print_phasing_paths(hs);
		//printf("---------\n");
		printf("=========\n");
		*/

		bool b = false;

		b = resolve_broken_vertex();
		if(b == true) continue;

		b = resolve_trivial_vertex(1, true, cfg.max_decompose_error_ratio[TRIVIAL_VERTEX]);
		if(b == true) continue;

		b = resolve_unsplittable_vertex(UNSPLITTABLE_SINGLE, 1, 0.01);
		if(b == true) continue;

		b = resolve_unsplittable_vertex(UNSPLITTABLE_SINGLE, INT_MAX, cfg.max_decompose_error_ratio[UNSPLITTABLE_SINGLE]);
		if(b == true) continue;

		b = resolve_splittable_vertex(SPLITTABLE_PURE, 1, cfg.max_decompose_error_ratio[SPLITTABLE_PURE]);
		if(b == true) continue;

		b = resolve_splittable_vertex(SPLITTABLE_HYPER, 1, cfg.max_decompose_error_ratio[SPLITTABLE_HYPER]);
		if(b == true) continue;

		b = resolve_splittable_vertex(SPLITTABLE_SIMPLE, 1, cfg.max_decompose_error_ratio[SPLITTABLE_SIMPLE]);
		if(b == true) continue;

		b = resolve_smallest_edge(cfg.max_decompose_error_ratio[SMALLEST_EDGE], 0);
		if(b == true) continue;

		b = resolve_smallest_edge(cfg.max_decompose_error_ratio[SMALLEST_EDGE], 1);
		if(b == true) continue;

		b = resolve_smallest_edge(cfg.max_decompose_error_ratio[SMALLEST_EDGE], 2);
		if(b == true) continue;

		b = resolve_splittable_vertex(SPLITTABLE_PURE, INT_MAX, cfg.max_decompose_error_ratio[SPLITTABLE_PURE]);
		if(b == true) continue;

		b = resolve_splittable_vertex(SPLITTABLE_HYPER, INT_MAX, cfg.max_decompose_error_ratio[SPLITTABLE_HYPER]);
		if(b == true) continue;

		b = resolve_splittable_vertex(SPLITTABLE_SIMPLE, INT_MAX, cfg.max_decompose_error_ratio[SPLITTABLE_SIMPLE]);
		if(b == true) continue;

		b = resolve_smallest_edge(cfg.max_decompose_error_ratio[SMALLEST_EDGE], 3);
		if(b == true) continue;

		/*
		b = resolve_unsplittable_vertex(UNSPLITTABLE_MULTIPLE, 1, 0.01);
		if(b == true) continue;

		b = resolve_unsplittable_vertex(UNSPLITTABLE_MULTIPLE, INT_MAX, cfg.max_decompose_error_ratio[UNSPLITTABLE_MULTIPLE]);
		if(b == true) continue;
		*/

		b = resolve_hyper_edge(2);
		if(b == true) continue;

		b = resolve_hyper_edge(1);
		if(b == true) continue;

		b = resolve_trivial_vertex(2, true, cfg.max_decompose_error_ratio[TRIVIAL_VERTEX]);
		if(b == true) continue;

		// process mixed vertices
		break_divided_phases();

		b = resolve_mixed_vertex(MIXED_DIVIDED);
		if(b == true) continue;

		b = resolve_mixed_smallest_edges();
		if(b == true) continue;

		b = resolve_mixed_vertex(MIXED_TRIVIAL);
		if(b == true) continue;

		b = resolve_mixed_vertex(MIXED_BLOCKED);
		if(b == true) continue;

		b = resolve_mixed_vertex(MIXED_SPLITTABLE);
		if(b == true) continue;

		b = resolve_mixed_vertex(MIXED_TANGLED);
		if(b == true) continue;

		break;
	}

	collect_existing_st_paths();

	greedy_decompose();

	build_transcripts();

	if(cfg.verbose >= 2) 
	{
		for(int i = 0; i < paths.size(); i++) paths[i].print(i);
		printf("finish assemble bundle %s\n\n", gr.gid.c_str());
	}

	return 0;
}

bool scallop::resolve_broken_vertex()
{
	vector<int> vv(nonzeroset.begin(), nonzeroset.end());
	int x = -1;
	for(int k = 0; k < vv.size(); k++)
	{
		int i = vv[k];
		if(i == 0) continue;
		if(i == gr.num_vertices() - 1) continue;
		if(gr.in_degree(i) >= 1 && gr.out_degree(i) >= 1) continue;
		x = i;
		break;
	}

	if(x == -1) return false;

	vector<int> ve;
	PEEI pei = gr.in_edges(x);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = (*it);
		assert(e2i.find(e) != e2i.end());
		ve.push_back(e2i[e]);
	}
	pei = gr.out_edges(x);
	for(edge_iterator it = pei.first; it != pei.second; it++)
	{
		edge_descriptor e = (*it);
		assert(e2i.find(e) != e2i.end());
		ve.push_back(e2i[e]);
	}

	assert(ve.size() >= 1);
	for(int k = 0; k < ve.size(); k++)
	{
		int e = ve[k];
		remove_edge(e);
		hs.remove(e);
	}

	assert(gr.degree(x) == 0);
	nonzeroset.erase(x);

	return true;
}

bool scallop::resolve_smallest_edge(double max_ratio, int degree)
{
	int root = -1;
	double best_ratio = max_ratio;
	vector<int> vv(nonzeroset.begin(), nonzeroset.end());
	for(int k = 0; k < vv.size(); k++)
	{
		int i = vv[k];
		double ratio;
		bool b = false;
		if(degree == 0) b = remove_single_smallest_edge(i, 0.01, ratio);
		if(degree == 1) b = thread_single_smallest_edge(i, 0.01, ratio);
		if(degree >= 2) b = thread_single_smallest_edge(i, 0.01, ratio, degree);

		//printf("resolve %d, ratio = %.3lf, b = %c | best = %.3lf, root = %d\n", i, ratio, b ? 'T' : 'F', best_ratio, root);

		if(b == true) return true;

		if(ratio >= 0 && ratio < best_ratio)
		{
			best_ratio = ratio;
			root = i;
		}
	}

	if(root >= 0 && degree == 0) return remove_single_smallest_edge(root, max_ratio, best_ratio);
	if(root >= 0 && degree == 1) return thread_single_smallest_edge(root, max_ratio, best_ratio);
	if(root >= 0 && degree >= 2) return thread_single_smallest_edge(root, max_ratio, best_ratio, degree);
	return false;
}

bool scallop::thread_single_smallest_edge(int x, double max_ratio, double &ratio, int degree)
{
	ratio = -1;
	if(degree <= 1) return false;
	if(gr.out_degree(x) <= 1) return false;
	if(gr.in_degree(x) <= 1) return false;
	if(gr.mixed_strand_vertex(x)) return false;

	double r;
	int e = compute_smallest_edge(x, r);
	if(e == -1) return false;

	edge_descriptor ee = i2e[e];
	int es = ee->source();
	int et = ee->target();
	assert(es == x || et == x);

	MI s;
	if(et == x) s = hs.get_successors(e);
	if(es == x) s = hs.get_predecessors(e);
	if(s.size() != degree) return false;

	if(r > max_ratio)
	{
		ratio = r;
		return false;
	}

	/*
	PEEI pe1, pe2;
	edge_iterator it1, it2, ot1, ot2;
	for(pe1 = gr.in_edges(x), it1 = pe1.first, it2 = pe1.second; it1 != it2; it1++)
	{
		int e1 = e2i[*it1];
		double w1 = gr.get_edge_weight(*it1);
		printf("in edges of %d: %d - %.3lf\n", x, e1, w1);
	}
	for(pe2 = gr.out_edges(x), ot1 = pe2.first, ot2 = pe2.second; ot1 != ot2; ot1++)
	{
		int e2 = e2i[*ot1];
		double w2 = gr.get_edge_weight(*ot1);
		printf("out edges of %d: %d - %.3lf\n", x, e2, w2);
	}
	*/

	double we = gr.get_edge_weight(ee);
	double sum;
	for(auto &a : s)
	{
		int f = a.first;
		edge_descriptor ff = i2e[f];
		double w = gr.get_edge_weight(ff);
		sum += w;
	}

	double ww1 = gr.get_vertex_weight(x) * r;
	double ww2 = gr.get_vertex_weight(x) - ww1;

	// vertex-n => new sink vertex
	// vertex-(n-1) => splitted vertex for xe and ye
	// vertex-x => splitted vertex for xe2 and ye2

	int n = gr.num_vertices();
	assert(v2v.size() == n);
	gr.add_vertex();
	assert(nonzeroset.find(n - 1) == nonzeroset.end());
	nonzeroset.insert(n - 1);
	gr.set_vertex_info(n, gr.get_vertex_info(n - 1));
	gr.set_vertex_info(n - 1, gr.get_vertex_info(x));
	gr.set_vertex_weight(n, gr.get_vertex_weight(n - 1));
	gr.set_vertex_weight(n - 1, ww1);
	gr.set_vertex_weight(x, ww2);

	v2v.push_back(v2v[n - 1]);
	v2v[n - 1] = v2v[x];

	// use vertex-n instead of vertex-(n-1) as sink vertex
	exchange_sink(n - 1, n);
	if(et == n - 1) et = n;

	assert(gr.degree(n - 1) == 0);

	// attach edges in xe and ye to vertex-(n-1)
	if(et == x) gr.move_edge(ee, es, n - 1);
	if(es == x) gr.move_edge(ee, n - 1, et);

	// create new edges
	for(auto &a : s)
	{
		int f = a.first;
		edge_descriptor ff = i2e[f];
		int fs = ff->source();
		int ft = ff->target();

		edge_descriptor pe = null_edge;
		if(et == x) pe = gr.add_edge(n - 1, ft);
		if(es == x) pe = gr.add_edge(fs, n - 1);

		int z = i2e.size();
		i2e.push_back(pe);
		e2i.insert(PEI(pe, z));

		if(et == x) hs.replace(e, f, e, z);
		if(es == x) hs.replace(f, e, z, e);

		double wf = gr.get_edge_weight(ff);
		double wa = wf / sum * we;
		double wb = wf - wa;
		if(wa < cfg.min_guaranteed_edge_weight) wa = cfg.min_guaranteed_edge_weight;
		if(wb < cfg.min_guaranteed_edge_weight) wb = cfg.min_guaranteed_edge_weight;

		gr.set_edge_weight(pe, wa);
		gr.set_edge_weight(ff, wb);
		gr.set_edge_info(pe, gr.get_edge_info(ff));

		// TODO, might be buggy
		vector<int> v0 = mev[ff];
		if(mev.find(pe) != mev.end()) mev[pe] = v0;
		else mev.insert(PEV(pe, v0));

		double dd = med[ff] * wf / sum;
		if(med.find(pe) != med.end()) med[pe] = dd;
		else med.insert(PED(pe, dd));
		med[ff] = med[ff] - dd;

		int l = mei[ff];
		if(mei.find(pe) != mei.end()) mei[pe] = l;
		else mei.insert(PEI(pe, l));

		borrow_edge_strand(z, f);

		//printf("add edge (%d, %d), id = %d, x = %d, f = %d, fs = %d, ft = %d, z = %d, dd = %.3lf, l = %d, ind = %d, outd = %d\n", pe->source(), pe->target(), z, x, f, fs, ft, z, dd, f, gr.in_degree(n - 1), gr.out_degree(n - 1));
	}

	/*
	for(pe1 = gr.in_edges(x), it1 = pe1.first, it2 = pe1.second; it1 != it2; it1++)
	{
		int e1 = e2i[*it1];
		double w1 = gr.get_edge_weight(*it1);
		printf("in edges of %d: %d - %.3lf\n", x, e1, w1);
	}
	for(pe2 = gr.out_edges(x), ot1 = pe2.first, ot2 = pe2.second; ot1 != ot2; ot1++)
	{
		int e2 = e2i[*ot1];
		double w2 = gr.get_edge_weight(*ot1);
		printf("out edges of %d: %d - %.3lf\n", x, e2, w2);
	}

	for(pe1 = gr.in_edges(n-1), it1 = pe1.first, it2 = pe1.second; it1 != it2; it1++)
	{
		int e1 = e2i[*it1];
		double w1 = gr.get_edge_weight(*it1);
		printf("in edges of %d: %d - %.3lf\n", n-1, e1, w1);
	}
	for(pe2 = gr.out_edges(n-1), ot1 = pe2.first, ot2 = pe2.second; ot1 != ot2; ot1++)
	{
		int e2 = e2i[*ot1];
		double w2 = gr.get_edge_weight(*ot1);
		printf("out edges of %d: %d - %.3lf\n", n-1, e2, w2);
	}
	*/

	if(cfg.verbose >= 2) 
	{
		int32_t p1 = gr.get_vertex_info(x).rpos;
		int32_t p2 = gr.get_vertex_info(x).lpos;
		printf("split smallest edge, edge = %d, degree = %d, weight = %.2lf, ratio = %.3lf, vertex = %d, pos = (%d, %d)\n", 
				e, degree, we, r, x, p1, p2);
	}
	return true;
}

bool scallop::thread_single_smallest_edge(int i, double max_ratio, double &ratio)
{
	ratio = -1;
	if(gr.out_degree(i) <= 1) return false;
	if(gr.in_degree(i) <= 1) return false;
	if(gr.mixed_strand_vertex(i)) return false;

	double r;
	int e = compute_smallest_edge(i, r);
	if(e == -1) return false;

	edge_descriptor ee = i2e[e];
	int es = ee->source();
	int et = ee->target();
	assert(es == i || et == i);

	MI s;
	if(et == i) s = hs.get_successors(e);
	if(es == i) s = hs.get_predecessors(e);
	if(s.size() != 1) return false;

	if(r > max_ratio)
	{
		ratio = r;
		return false;
	}

	int f = s.begin()->first;
	edge_descriptor ff = i2e[f];
	int fs = ff->source();
	int ft = ff->target();

	double we = gr.get_edge_weight(ee);
	double wf = gr.get_edge_weight(ff);

	if(wf <= we + cfg.min_guaranteed_edge_weight) wf = we + cfg.min_guaranteed_edge_weight;
	gr.set_edge_weight(ff, wf);

	if(et == i)
	{
		assert(et == fs);
		int g = merge_adjacent_edges(e, f);
		assert(g >= 0);
		hs.replace(e, f, g);
		hs.replace(e, g);
	}
	if(es == i)
	{
		assert(es == ft);
		int g = merge_adjacent_edges(f, e);
		assert(g >= 0);
		hs.replace(f, e, g);
		hs.replace(e, g);
	}

	if(cfg.verbose >= 2) 
	{
		int32_t p1 = gr.get_vertex_info(i).rpos;
		int32_t p2 = gr.get_vertex_info(i).lpos;
		printf("thread smallest edge, edge = %d, weight = %.2lf, ratio = %.3lf, vertex = %d, pos = (%d, %d)\n", 
				e, we, r, i, p1, p2);
	}
	return true;
}

bool scallop::remove_single_smallest_edge(int i, double max_ratio, double &ratio)
{
	double r;
	ratio = -1;
	int e = compute_smallest_edge(i, r);
	if(e == -1) return false;

	int s = i2e[e]->source();
	int t = i2e[e]->target();
	assert(s == i || t == i);

	if(gr.out_degree(s) <= 1) return false;
	if(gr.in_degree(t) <= 1) return false;

	if(hs.right_extend(e) && hs.left_extend(e)) return false;
	if(t == i && hs.right_extend(e)) return false;
	if(s == i && hs.left_extend(e)) return false;

	//if(gr.mixed_strand_vertex(i)) return false;
	// check if there are at least two possible ways 
	// for e to go to the other side
	vector<int> vs = gr.get_strand_degree(i);
	int z = gr.get_edge_info(i2e[e]).strand;
	if(s == i && z >= 1 && vs[0] + vs[z + 0] <= 1) return false;
	if(t == i && z >= 1 && vs[3] + vs[z + 3] <= 1) return false;

	ratio = r;
	if(r > max_ratio) return false;

	double w = gr.get_edge_weight(i2e[e]);
	if(cfg.verbose >= 2) 
	{
		int32_t p1 = gr.get_vertex_info(s).rpos;
		int32_t p2 = gr.get_vertex_info(t).lpos;
		printf("resolve small edge, edge = %d, weight = %.2lf, ratio = %.2lf, vertex = (%d, %d), degree = (%d, %d), pos = (%d, %d)\n", 
				e, w, r, s, t, gr.out_degree(s), gr.in_degree(t), p1, p2);
	}
	remove_edge(e);
	hs.remove(e);
	return true;
}

bool scallop::resolve_splittable_vertex(int type, int degree, double max_ratio)
{
	int root = -1;
	int min_balance = INT_MAX;
	double min_ratio = max_ratio;
	equation eqn;
	//vector<equation> eqns;
	//for(set<int>::iterator it = nonzeroset.begin(); it != nonzeroset.end(); it++)
	//for(int i = 1; i < gr.num_vertices() - 1; i++)
	vector<int> vv(nonzeroset.begin(), nonzeroset.end());
	for(int k = 0; k < vv.size(); k++)
	{
		int i = vv[k];
		if(gr.in_degree(i) <= 1) continue;
		if(gr.out_degree(i) <= 1) continue;

		MPII mpi = hs.get_routes(i, gr, e2i);
		router rt(i, gr, e2i, i2e, mpi, cfg);
		rt.classify();

		if(rt.type != type) continue;
		if(rt.degree > degree) continue;

		rt.build();

		printf("type = %d, degree = %d, rt.degree = %d\n", type, degree, rt.degree);

		assert(rt.eqns.size() == 2);

		int balance = abs((int)(rt.eqns[0].s.size() - rt.eqns[0].t.size())) + abs((int)(gr.in_degree(i) - rt.eqns[0].s.size()) - (int)(gr.out_degree(i) - rt.eqns[0].t.size()));

		if(rt.ratio > max_ratio) continue;
		if(balance > min_balance) continue;
		if(balance == min_balance && rt.ratio > min_ratio) continue;

		root = i;
		min_ratio = rt.ratio;
		min_balance = balance;
		eqn = rt.eqns[0];
	}

	if(root == -1) return false;

	if(cfg.verbose >= 2) printf("resolve splittable vertex, type = %d, degree = %d, vertex = %d, %d-%d, ratio = %.2lf, balance = %d, degree = (%d, %d), eqn = (%lu, %lu)\n",
			type, degree, root, gr.get_vertex_info(root).lpos, gr.get_vertex_info(root).rpos, min_ratio, min_balance, gr.in_degree(root), gr.out_degree(root), eqn.s.size(), eqn.t.size());

	split_vertex(root, eqn.s, eqn.t);

	return true;
}

bool scallop::resolve_unsplittable_vertex(int type, int degree, double max_ratio)
{
	int root = -1;
	MPID pe2w;
	double ratio = max_ratio;
	bool flag = false;
	//for(int i = 1; i < gr.num_vertices() - 1; i++)
	//for(set<int>::iterator it = nonzeroset.begin(); it != nonzeroset.end(); it++)
	vector<int> vv(nonzeroset.begin(), nonzeroset.end());
	for(int k = 0; k < vv.size(); k++)
	{
		int i = vv[k];
		if(gr.in_degree(i) <= 1) continue;
		if(gr.out_degree(i) <= 1) continue;

		MPII mpi = hs.get_routes(i, gr, e2i);
		router rt(i, gr, e2i, i2e, mpi, cfg);
		rt.classify();

		if(rt.type != type) continue;
		if(rt.degree > degree) continue;

		rt.build();
		//rt.print();

		if(rt.ratio < 0.01)
		{
			if(cfg.verbose >= 2) printf("resolve unsplittable vertex, type = %d, degree = %d, vertex = %d, %d-%d, ratio = %.3lf, degree = (%d, %d)\n",
					type, degree, i, gr.get_vertex_info(i).lpos, gr.get_vertex_info(i).rpos, rt.ratio, gr.in_degree(i), gr.out_degree(i));
			decompose_vertex_extend(i, rt.pe2w);
			flag = true;
			continue;
		}

		if(rt.ratio > ratio) continue;

		root = i;
		ratio = rt.ratio;
		pe2w = rt.pe2w;
	}

	if(flag == true) return true;
	if(root == -1) return false;

	if(cfg.verbose >= 2) printf("resolve unsplittable vertex, type = %d, degree = %d, vertex = %d, %d-%d, ratio = %.3lf, degree = (%d, %d)\n",
			type, degree, root, gr.get_vertex_info(root).lpos, gr.get_vertex_info(root).rpos, ratio, gr.in_degree(root), gr.out_degree(root));

	decompose_vertex_extend(root, pe2w);
	return true;
}

bool scallop::resolve_hyper_edge(int fsize)
{
	edge_iterator it1, it2;
	PEEI pei;
	vector<int> v1, v2;
	int root = -1;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int e = e2i[*it1];
		int vs = (*it1)->source();
		int vt = (*it1)->target();

		MI s;
		s = hs.get_successors(e);
		//if(s.size() >= 2 && hs.right_extend(get_keys(s)) == false && (hs.left_extend(e) == false || gr.out_degree(vs) == 1))
		if(gr.mixed_strand_vertex(vt) == false && s.size() >= fsize && (hs.left_extend(e) == false || gr.out_degree(vs) == 1))
		{
			v1.push_back(e);
			v2 = get_keys(s);
			root = vt;
			break;
		}

		s = hs.get_predecessors(e);
		//if(s.size() >= 2 && hs.left_extend(get_keys(s)) == false && (hs.right_extend(e) == false || gr.in_degree(vt) == 1))
		if(gr.mixed_strand_vertex(vs) == false && s.size() >= fsize && (hs.right_extend(e) == false || gr.in_degree(vt) == 1))
		{
			v1 = get_keys(s);
			v2.push_back(e);
			root = vs;
			break;
		}
	}

	if(root < 0) return false;
	if(gr.in_degree(root) <= 0) return false;
	if(gr.out_degree(root) <= 0) return false;
	if(v1.size() == 0 || v2.size() == 0) return false;
	assert(v1.size() == 1 || v2.size() == 1);

	if(cfg.verbose >= 2) printf("resolve hyper edge, fsize = %d, vertex = %d, %d-%d, degree = (%d, %d), hyper edge = (%lu, %lu)\n",
			fsize, root, gr.get_vertex_info(root).lpos, gr.get_vertex_info(root).rpos, gr.in_degree(root), gr.out_degree(root), v1.size(), v2.size());

	balance_vertex(root);

	vector<double> w1, w2;
	double sum1 = 0, sum2 = 0;
	for(int i = 0; i < v1.size(); i++)
	{
		double w = gr.get_edge_weight(i2e[v1[i]]);
		w1.push_back(w);
		sum1 += w;
	}
	for(int i = 0; i < v2.size(); i++)
	{
		double w = gr.get_edge_weight(i2e[v2[i]]);
		w2.push_back(w);
		sum2 += w;
	}

	double r1 = (sum1 < sum2) ? 1.0 : sum2 / sum1;
	double r2 = (sum1 > sum2) ? 1.0 : sum1 / sum2;

	for(int i = 0; i < w1.size(); i++) w1[i] *= r1;
	for(int i = 0; i < w2.size(); i++) w2[i] *= r2;

	for(int i = 0; i < w1.size(); i++)
	{
		for(int j = 0; j < w2.size(); j++)
		{
			if(w1[i] < cfg.min_guaranteed_edge_weight)
			{
				double e1 = gr.get_edge_weight(i2e[v1[i]]) + cfg.min_guaranteed_edge_weight - w1[i];
				double e2 = gr.get_edge_weight(i2e[v2[j]]) + cfg.min_guaranteed_edge_weight - w1[i];
				w2[j] = w2[j] + cfg.min_guaranteed_edge_weight - w1[i];
				w1[i] = cfg.min_guaranteed_edge_weight;
				gr.set_edge_weight(i2e[v1[i]], e1);
				gr.set_edge_weight(i2e[v2[j]], e2);
			}

			if(w2[j] < cfg.min_guaranteed_edge_weight)
			{
				double e1 = gr.get_edge_weight(i2e[v1[i]]) + cfg.min_guaranteed_edge_weight - w2[j];
				double e2 = gr.get_edge_weight(i2e[v2[j]]) + cfg.min_guaranteed_edge_weight - w2[j];
				w1[i] = w1[i] + cfg.min_guaranteed_edge_weight - w2[j];
				w2[j] = cfg.min_guaranteed_edge_weight;
				gr.set_edge_weight(i2e[v1[i]], e1);
				gr.set_edge_weight(i2e[v2[j]], e2);
			}
		}
	}

	set<int> ss;
	bool flag = false;
	for(int i = 0; i < w1.size(); i++)
	{
		for(int j = 0; j < w2.size(); j++)
		{
			double w = (w1[i] < w2[j]) ? w1[i] : w2[j];
			assert(w >= cfg.min_guaranteed_edge_weight - SMIN);

			flag = true;
			int k1 = split_edge(v1[i], w);
			int k2 = split_edge(v2[j], w);
			int x = merge_adjacent_equal_edges(k1, k2);

			//printf(" split (%d, %d), w = %.2lf, weight = (%.2lf, %.2lf), merge (%d, %d) -> %d\n", v1[i], v2[j], w, w1[i], w2[j], k1, k2, x);

			hs.replace(v1[i], v2[j], x);
			if(k1 == v1[i]) hs.remove(v1[i]);
			if(k2 == v2[j]) hs.remove(v2[j]);
			//if(k1 == v1[i]) hs.replace(v1[i], x);
			//if(k2 == v2[j]) hs.replace(v2[j], x);
		}
	}
	return flag;
}

bool scallop::resolve_trivial_vertex(int type, bool fast, double jump_ratio)
{
	int root = -1;
	double ratio = DBL_MAX;
	int se = -1;
	bool flag = false;
	//for(int i = 1; i < gr.num_vertices() - 1; i++)
	//for(set<int>::iterator it = nonzeroset.begin(); it != nonzeroset.end(); it++)
	vector<int> vv(nonzeroset.begin(), nonzeroset.end());
	for(int k = 0; k < vv.size(); k++)
	{
		int i = vv[k];
		if(gr.in_degree(i) <= 0) continue;
		if(gr.out_degree(i) <= 0) continue;
		if(gr.mixed_strand_vertex(i)) continue;
		assert(gr.in_degree(i) >= 1);
		assert(gr.out_degree(i) >= 1);

		if(gr.in_degree(i) >= 2 && gr.out_degree(i) >= 2) continue;
		if(classify_trivial_vertex(i, fast) != type) continue;

		int e;
		double r = compute_balance_ratio(i);

		if(r < 1.02)
		{
			if(cfg.verbose >= 2) printf("resolve trivial vertex %d, %d-%d, type = %d, ratio = %.2lf, degree = (%d, %d), fast = %c\n", 
					i, gr.get_vertex_info(i).lpos, gr.get_vertex_info(i).rpos, type, r, gr.in_degree(i), gr.out_degree(i), fast ? 'T' : 'F');
			decompose_trivial_vertex(i);
			flag = true;
			continue;
		}

		if(ratio < r) continue;

		root = i;
		ratio = r;
		se = e;

		if(ratio < jump_ratio) break;
	}

	if(flag == true) return true;
	if(root == -1) return false;

	if(cfg.verbose >= 2) printf("resolve trivial vertex %d, %d-%d, type = %d, ratio = %.2lf, degree = (%d, %d)\n", root, gr.get_vertex_info(root).lpos, gr.get_vertex_info(root).rpos, type, 
			ratio, gr.in_degree(root), gr.out_degree(root));

	decompose_trivial_vertex(root);

	assert(gr.degree(root) == 0);
	return true;
}

bool scallop::resolve_single_trivial_vertex(int i, double jump_ratio)
{
	if(gr.in_degree(i) <= 0) return false;
	if(gr.out_degree(i) <= 0) return false;
	if(gr.in_degree(i) >= 2 && gr.out_degree(i) >= 2) return false;
	if(gr.mixed_strand_vertex(i)) return false;
	if(classify_trivial_vertex(i, false) != 1) return false;

	double r = compute_balance_ratio(i);
	if(r >= jump_ratio) return false;

	if(cfg.verbose >= 2) printf("resolve trivial vertex fast, vertex = %d, %d-%d, ratio = %.2lf, degree = (%d, %d)\n",
			i, gr.get_vertex_info(i).lpos, gr.get_vertex_info(i).rpos, r, gr.in_degree(i), gr.out_degree(i));

	decompose_trivial_vertex(i);
	assert(gr.degree(i) == 0);

	return true;
}

bool scallop::resolve_mixed_vertex(int type)
{
	int root = -1;
	vector<int> vv(nonzeroset.begin(), nonzeroset.end());
	for(int k = 0; k < vv.size(); k++)
	{
		int i = vv[k];
		if(gr.out_degree(i) <= 0) continue;
		if(gr.in_degree(i) <= 0) continue;
		if(gr.mixed_strand_vertex(i) == false) continue;

		MPII mpi = hs.get_routes(i, gr, e2i);
		router rt(i, gr, e2i, i2e, mpi, cfg);
		rt.classify();
		if(rt.type != type) continue;

		root = i;

		if(cfg.verbose >= 2)
		{
			printf("resolve mixed vertex, type = %d, root = %d, %d-%d, v2v[root] = %d, degree = (%d, %d)\n", type, root, gr.get_vertex_info(root).lpos, gr.get_vertex_info(root).rpos, v2v[root], gr.in_degree(root), gr.out_degree(root));
			gr.print_vertex(root);
		}

		if(type == MIXED_DIVIDED)
		{
			terminate_blocked_edges(root);
			assert(gr.degree(root) == 0);
			return true;
		}

		if(type == MIXED_TRIVIAL)
		{
			decompose_trivial_vertex(root);
			assert(gr.degree(root) == 0);
			return true;
		}

		if(type == MIXED_SPLITTABLE)
		{
			assert(rt.eqns.size() >= 1);
			split_vertex(root, rt.eqns[0].s, rt.eqns[0].t);
			return true;
		}

		if(type == MIXED_BLOCKED)
		{
			terminate_blocked_edges(root);
			return true;
		}

		if(type == MIXED_TANGLED)
		{
			terminate_smallest_edge(root);
			return true;
		}
	}
	return false;
}

bool scallop::resolve_mixed_smallest_edges()
{
	vector<int> vv(nonzeroset.begin(), nonzeroset.end());
	for(int k = 0; k < vv.size(); k++)
	{
		int i = vv[k];
		if(gr.out_degree(i) <= 0) continue;
		if(gr.in_degree(i) <= 0) continue;

		double r;
		int e = compute_smallest_edge(i, r);
		int s = i2e[e]->source();
		int t = i2e[e]->target();
		assert(s == i || t == i);

		if(gr.out_degree(s) <= 1) continue;
		if(gr.in_degree(t) <= 1) continue;

		if(hs.right_extend(e) && hs.left_extend(e)) continue;
		if(t == i && hs.right_extend(e)) continue;
		if(s == i && hs.left_extend(e)) continue;

		vector<int> vs = gr.get_strand_degree(i);
		int z = gr.get_edge_info(i2e[e]).strand;
		if(s == i && z >= 1 && vs[0] + vs[z + 0] >= 1) continue;
		if(t == i && z >= 1 && vs[3] + vs[z + 3] >= 1) continue;

		if(s == i) assert(gr.out_degree(i) >= 2);
		if(t == i) assert(gr.in_degree(i) >= 2);

		double w = gr.get_edge_weight(i2e[e]);
		if(cfg.verbose >= 2) printf("remove blocked smallest edge, edge = %d, weight = %.2lf, ratio = %.2lf, vertex = (%d, %d), degree = (%d, %d)\n", 
				e, w, r, s, t, gr.out_degree(s), gr.in_degree(t));

		remove_edge(e);
		hs.remove(e);

		return true;
	}
	return false;
}

int scallop::terminate_blocked_edges(int root)
{
	PEEI pi = gr.in_edges(root);
	PEEI po = gr.out_edges(root);
	int xi[3] = {0, 0, 0};
	int xo[3] = {0, 0, 0};

	for(edge_iterator it = pi.first; it != pi.second; it++)
	{
		edge_descriptor e = (*it);
		int s = gr.get_edge_info(e).strand;
		xi[s]++;
	}
	for(edge_iterator it = po.first; it != po.second; it++)
	{
		edge_descriptor e = (*it);
		int s = gr.get_edge_info(e).strand;
		xo[s]++;
	}

	vector<edge_descriptor> vei;
	vector<edge_descriptor> veo;
	for(edge_iterator it = pi.first; it != pi.second; it++)
	{
		edge_descriptor e = (*it);
		int s = gr.get_edge_info(e).strand;
		if(s == 0) continue;
		if(xo[0] >= 1) continue;
		if(xo[s] >= 1) continue;
		vei.push_back(e);
	}

	for(edge_iterator it = po.first; it != po.second; it++)
	{
		edge_descriptor e = (*it);
		int s = gr.get_edge_info(e).strand;
		if(s == 0) continue;
		if(xi[0] >= 1) continue;
		if(xi[s] >= 1) continue;
		veo.push_back(e);
	}

	assert(vei.size() + veo.size() >= 1);

	for(int i = 0; i < vei.size(); i++) terminate_edge_sink(vei[i]);
	for(int i = 0; i < veo.size(); i++) terminate_edge_source(veo[i]);

	if(gr.degree(root) == 0) nonzeroset.erase(root);
	else assert(gr.in_degree(root) >= 1 && gr.out_degree(root) >= 1);

	return 0;
}

int scallop::terminate_smallest_edge(int root)
{
	if(gr.in_degree(root) >= 2)
	{
		edge_descriptor e = gr.min_in_edge(root);
		terminate_edge_sink(e);
	}
	else if(gr.out_degree(root) >= 2)
	{
		edge_descriptor e = gr.min_out_edge(root);
		terminate_edge_source(e);
	}
	else
	{
		assert(gr.in_degree(root) == 1);
		assert(gr.out_degree(root) == 1);
		edge_descriptor e1 = gr.min_in_edge(root);
		edge_descriptor e2 = gr.min_out_edge(root);
		terminate_edge_sink(e1);
		terminate_edge_source(e2);
		assert(gr.degree(root) == 0);
		nonzeroset.erase(root);
	}
	return 0;
}

int scallop::terminate_edge_sink(edge_descriptor e)
{
	int t = e->target();
	int n = gr.num_vertices() - 1;
	double vw = gr.get_vertex_weight(t);
	double ew = gr.get_edge_weight(e);
	double ww = gr.get_out_weights(t) + gr.get_in_weights(t);
	double sw = vw * ew / ww;

	hs.right_break(e2i[e]);
	gr.move_edge(e, e->source(), n);
	mev[e].push_back(t);

	vertex_info vi = gr.get_vertex_info(t);
	assert(med.find(e) != med.end());
	assert(mei.find(e) != mei.end());
	med[e] += sw * (vi.rpos - vi.lpos);
	mei[e] += vi.rpos - vi.lpos;
	assert(vw >= sw - SMIN);
	gr.set_vertex_weight(t, vw - sw);

	return 0;
}

int scallop::terminate_edge_source(edge_descriptor e)
{
	int s = e->source();
	double vw = gr.get_vertex_weight(s);
	double ew = gr.get_edge_weight(e);
	double ww = gr.get_out_weights(s) + gr.get_in_weights(s);
	double sw = vw * ew / ww;

	hs.left_break(e2i[e]);
	gr.move_edge(e, 0, e->target());
	mev[e].insert(mev[e].begin(), s);

	vertex_info vi = gr.get_vertex_info(s);
	assert(med.find(e) != med.end());
	assert(mei.find(e) != mei.end());
	med[e] += sw * (vi.rpos - vi.lpos);
	mei[e] += vi.rpos - vi.lpos;
	assert(vw >= sw - SMIN);
	gr.set_vertex_weight(s, vw - sw);

	return 0;
}

int scallop::break_divided_phases()
{
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		vector<int> vs(6, 0);
		PEEI pei = gr.in_edges(i);
		PEEI peo = gr.out_edges(i);
		for(edge_iterator ti = pei.first; ti != pei.second; ti++)
		{
			edge_descriptor e1 = *ti;
			int s1 = gr.get_edge_info(e1).strand;
			if(s1 == 0) continue;
			for(edge_iterator to = peo.first; to != peo.second; to++)
			{
				edge_descriptor e2 = *to;
				int s2 = gr.get_edge_info(e2).strand;
				if(s2 == 0) continue;
				if(s1 == s2) continue;
				hs.insert_between(e2i[e1], e2i[e2], -1);
			}
		}
	}
	return 0;
}

int scallop::summarize_vertices()
{
	for(set<int>::iterator it = nonzeroset.begin(); it != nonzeroset.end(); it++)
	{
		int i = (*it);
		assert(gr.in_degree(i) >= 1);
		assert(gr.out_degree(i) >= 1);
		
		if(gr.in_degree(i) == 1 || gr.out_degree(i) == 1)
		{
			int c = classify_trivial_vertex(i, false);
			printf("summary: trivial vertex %d, type = %d\n", i, c);
		}
		else
		{
			MPII mpi = hs.get_routes(i, gr, e2i);
			MPII mpx = hs.get_routes(i, gr, e2i);
			set<int> s1;
			set<int> s2;
			for(MPII::iterator it = mpi.begin(); it != mpi.end(); it++)
			{
				PI p = it->first;
				s1.insert(p.first);
				s2.insert(p.second);
			}
			router rt(i, gr, e2i, i2e, mpi, cfg);
			rt.classify();
			rt.build();
			printf("summary: nontrivial vertex %d, degree = (%d, %d), hyper edges = %lu, graph degree = (%lu, %lu), type = %d, degree = %d, ratio = %.3lf\n", 
					i, gr.in_degree(i), gr.out_degree(i), mpi.size(), s1.size(), s2.size(), rt.type, rt.degree, rt.ratio);
		}
	}
	return 0;
}

int scallop::classify()
{
	assert(gr.num_vertices() >= 2);
	if(gr.num_vertices() == 2) return TRIVIAL;

	string s;	

	long p0 = gr.compute_num_paths();
	long p1 = gr.num_edges() - gr.num_vertices() + 2;
	for(int i = 0; i < gr.num_vertices(); i++) 
	{
		if(gr.degree(i) == 0) p1++;
	}

	//assert(p0 >= p1);
	bool b = (p0 <= p1) ? true : false;

	if(p0 == p1) return TRIVIAL;
	else return NORMAL;
}

int scallop::add_pseudo_hyper_edges()
{
	for(int k = 1; k < gr.num_vertices() - 1; k++)
	{
		int s = -1, t = -1;
		double w1 = 0, w2 = 0;
		edge_iterator it1, it2;
		PEEI pei;
		for(pei = gr.in_edges(k), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			double w = gr.get_edge_weight(*it1);
			if(w <= w1) continue;
			w1 = w;
			s = (*it1)->source();
		}
		for(pei = gr.out_edges(k), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			double w = gr.get_edge_weight(*it1);
			if(w <= w2) continue;
			w2 = w;
			t = (*it1)->target();
		}
		if(s == -1 || t == -1) continue;
		if(w1 <= 10.0 || w2 <= 10.0) continue;
		if(s == 0) continue;
		if(t == gr.num_vertices() - 1) continue;

		vector<int> v;
		v.push_back(s - 1);
		v.push_back(k - 1);
		v.push_back(t - 1);
		
		hs.add_node_list(v, 1);
	}
	return 0;
}

int scallop::init_super_edges()
{
	mev.clear();
	med.clear();
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		vector<int> v;
		//int s = (*it1)->source();
		//v.push_back(s);
		mev.insert(PEV(*it1, v));
		med.insert(PED(*it1, 0));
		mei.insert(PEI(*it1, 0));
	}
	return 0;
}

int scallop::init_inner_weights()
{
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		double w = gr.get_edge_weight(e);
		edge_info ei = gr.get_edge_info(e);
		ei.weight = w;
		gr.set_edge_info(e, ei);
	}
	return 0;
}

int scallop::init_vertex_map()
{
	v2v.clear();
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		v2v.push_back(i);
	}
	return 0;
}

int scallop::init_nonzeroset()
{
	nonzeroset.clear();
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.degree(i) <= 0) continue;
		nonzeroset.insert(i);
	}
	return 0;
}

int scallop::decompose_vertex_extend(int root, MPID &pe2w)
{
	PEEI pei;
	edge_iterator it1, it2;

	/*
	// print 
	printf(" in-weights: ");
	for(pei = gr.in_edges(root), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		double w = gr.get_edge_weight(e);
		printf("%d:%.2lf, ", e2i[e], w);
	}
	printf("\n");
	printf(" out-weights: ");
	for(pei = gr.out_edges(root), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		double w = gr.get_edge_weight(e);
		printf("%d:%.2lf, ", e2i[e], w);
	}
	printf("\n");
	printf(" decompose weights: ");
	for(MPID::iterator it = pe2w.begin(); it != pe2w.end(); it++)
	{
		printf("%d:%d:%.2lf, ", it->first.first, it->first.second, it->second);
	}
	printf("\n");
	*/
	// end print

	// compute degree of each edge
	map<int, int> mdegree;
	for(MPID::iterator it = pe2w.begin(); it != pe2w.end(); it++)
	{
		PI p = it->first;
		if(mdegree.find(p.first) == mdegree.end()) mdegree.insert(PI(p.first, 1));
		else mdegree[p.first]++;
		if(mdegree.find(p.second) == mdegree.end()) mdegree.insert(PI(p.second, 1));
		else mdegree[p.second]++;
	}

	// compute weight of each edge
	double total_weight = 0;
	map<int, double> mweight;
	for(MPID::iterator it = pe2w.begin(); it != pe2w.end(); it++)
	{
		PI p = it->first;
		double w = it->second;
		assert(w >= cfg.min_guaranteed_edge_weight - SMIN);
		total_weight += w;
		if(mdegree.find(p.first) == mdegree.end()) mdegree.insert(PI(p.first, w));
		else mdegree[p.first] += w;
		if(mdegree.find(p.second) == mdegree.end()) mdegree.insert(PI(p.second, w));
		else mdegree[p.second] += w;
	}

	// distribute
	vertex_info root_info = gr.get_vertex_info(root);
	double vertex_weight = gr.get_vertex_weight(root) * (root_info.rpos - root_info.lpos);
	for(map<int, double>::iterator it = mweight.begin(); it != mweight.end(); it++)
	{
		double w = it->second / total_weight * vertex_weight;
		it->second = w;
	}

	// add edge-vertex for each adjacent edge of root
	// for those with mdegree >= 2
	// map adjacent edges of root to vertices [m, n)
	int m = gr.num_vertices() - 1;
	int n = m;
	map<int, int> ev1, ev2;
	for(pei = gr.in_edges(root), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		int ei = e2i[e];
		int s = e->source();
		int t = e->target();
		assert(t == root);
		assert(mdegree.find(ei) != mdegree.end());
		if(mdegree[ei] >= 2) ev1.insert(PI(ei, n++));
	}
	for(pei = gr.out_edges(root), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		int ei = e2i[e];
		int s = e->source();
		int t = e->target();
		assert(s == root);
		assert(mdegree.find(ei) != mdegree.end());
		if(mdegree[ei] >= 2) ev2.insert(PI(ei, n++));
	}

	// add vertices and exchange sink
	for(int i = m; i < n; i++) 
	{
		gr.add_vertex();
		assert(nonzeroset.find(i) == nonzeroset.end());
		nonzeroset.insert(i);
		v2v.push_back(-1);
	}
	v2v[n] = v2v[m];
	gr.set_vertex_info(n, gr.get_vertex_info(m));
	exchange_sink(m, n);

	// set vertex info for new vertices
	// detach edge from root to new vertex
	for(map<int, int>::iterator it = ev1.begin(); it != ev1.end(); it++)
	{
		edge_descriptor e = i2e[it->first];
		int k = it->second;

		vertex_info vi;
		vi.lpos = gr.get_vertex_info(e->source()).rpos;
		vi.rpos = gr.get_vertex_info(e->source()).rpos;

		gr.move_edge(e, e->source(), k);
		gr.set_vertex_info(k, vi);
		gr.set_vertex_weight(k, 0);
		v2v[k] = -2; //v2v[root];
	}
	for(map<int, int>::iterator it = ev2.begin(); it != ev2.end(); it++)
	{
		edge_descriptor e = i2e[it->first];
		int k = it->second;

		vertex_info vi;
		vi.lpos = gr.get_vertex_info(e->target()).lpos;
		vi.rpos = gr.get_vertex_info(e->target()).lpos;

		gr.move_edge(e, k, e->target());
		gr.set_vertex_info(k, vi);
		gr.set_vertex_weight(k, 0);
		v2v[k] = -2;
	}

	// connecting edges according to pe2w
	for(MPID::iterator it = pe2w.begin(); it != pe2w.end(); it++)
	{
		int e1 = it->first.first;
		int e2 = it->first.second;
		assert(consistent_strands(e1, e2));
		double w = it->second;
		assert(w >= cfg.min_guaranteed_edge_weight - SMIN);

		if(mdegree[e1] == 1)
		{
			assert(mdegree[e2] >= 2);
			assert(ev1.find(e1) == ev1.end());
			assert(ev2.find(e2) != ev2.end());
			edge_descriptor p = i2e[e1];
			borrow_edge_strand(e1, e2);
			int v1 = p->source();
			int v2 = ev2[e2];
			gr.move_edge(p, v1, v2);

			mev[p].push_back(root);

			assert(med.find(p) != med.end());
			assert(mei.find(p) != mei.end());
			med[p] += mweight[e1];
			mei[p] += root_info.rpos - root_info.lpos;
		}
		else if(mdegree[e2] == 1)
		{
			assert(mdegree[e1] >= 2);
			assert(ev1.find(e1) != ev1.end());
			assert(ev2.find(e2) == ev2.end());
			edge_descriptor p = i2e[e2];
			borrow_edge_strand(e2, e1);
			int v1 = ev1[e1];
			int v2 = p->target();
			gr.move_edge(p, v1, v2);

			mev[p].insert(mev[p].begin(), root);

			assert(med.find(p) != med.end());
			assert(mei.find(p) != mei.end());
			med[p] += mweight[e2];
			mei[p] += root_info.rpos - root_info.lpos;
		}
		else
		{
			assert(mdegree[e1] >= 2);
			assert(mdegree[e2] >= 2);

			assert(ev1.find(e1) != ev1.end());
			assert(ev2.find(e2) != ev2.end());
			int v1 = ev1[e1];
			int v2 = ev2[e2];

			edge_descriptor p = gr.add_edge(v1, v2);
		
			int z = i2e.size();
			i2e.push_back(p);
			e2i.insert(PEI(p, z));

			gr.set_edge_weight(p, w);
			gr.set_edge_info(p, edge_info());

			vector<int> v0;
			v0.push_back(root);
			if(mev.find(p) != mev.end()) mev[p] = v0;
			else mev.insert(PEV(p, v0));

			double wd = w / total_weight * vertex_weight;
			if(med.find(p) != med.end()) med[p] = wd;
			else med.insert(PED(p, wd));

			int l = root_info.rpos - root_info.lpos;
			if(mei.find(p) != mei.end()) mei[p] = l;
			else mei.insert(PEI(p, l));

			borrow_edge_strand(z, e1);
			borrow_edge_strand(z, e2);
			hs.insert_between(e1, e2, z);
		}
	}

	assert(gr.degree(root) == 0);
	nonzeroset.erase(root);

	for(map<int, int>::iterator it = ev1.begin(); it != ev1.end(); it++)
	{
		int k = it->second;
		resolve_single_trivial_vertex(k, cfg.max_decompose_error_ratio[TRIVIAL_VERTEX]);
	}
	for(map<int, int>::iterator it = ev2.begin(); it != ev2.end(); it++)
	{
		int k = it->second;
		resolve_single_trivial_vertex(k, cfg.max_decompose_error_ratio[TRIVIAL_VERTEX]);
	}

	return 0;
}

bool scallop::consistent_strands(int e1, int e2)
{
	int s1 = gr.get_edge_info(i2e[e1]).strand;
	int s2 = gr.get_edge_info(i2e[e1]).strand;
	if(s1 == 1 && s2 == 2) return false;
	if(s1 == 2 && s2 == 1) return false;
	return true;
}

int scallop::borrow_edge_strand(int e1, int e2)
{
	// set the strand for e1 using e2
	int s1 = gr.get_edge_info(i2e[e1]).strand;
	int s2 = gr.get_edge_info(i2e[e2]).strand;
	if(s2 == 0) return 0;
	edge_info ei = gr.get_edge_info(i2e[e1]);
	ei.strand = s2;
	gr.set_edge_info(i2e[e1], ei);
	return 0;
}

int scallop::decompose_vertex_replace(int root, MPID &pe2w)
{
	// reassign weights
	MID md;
	for(MPID::iterator it = pe2w.begin(); it != pe2w.end(); it++)
	{
		int e1 = it->first.first;
		int e2 = it->first.second;
		double w = it->second;
		assert(w >= cfg.min_guaranteed_edge_weight - SMIN);
		if(md.find(e1) == md.end()) md.insert(PID(e1, w));
		else md[e1] += w;
		if(md.find(e2) == md.end()) md.insert(PID(e2, w));
		else md[e2] += w;
	}
	for(MID::iterator it = md.begin(); it != md.end(); it++)
	{
		edge_descriptor e = i2e[it->first];
		double w = it->second;
		gr.set_edge_weight(e, w);
	}

	// assert that all edges are covered
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr.in_edges(root), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int e = e2i[*it1];
		assert(md.find(e) != md.end());
	}
	for(pei = gr.out_edges(root), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int e = e2i[*it1];
		assert(md.find(e) != md.end());
	}

	// remove hyper-edges that are not covered
	MPII mpi = hs.get_routes(root, gr, e2i);
	for(MPII::iterator it = mpi.begin(); it != mpi.end(); it++)
	{
		assert(pe2w.find(it->first) != pe2w.end());
		if(pe2w.find(it->first) != pe2w.end()) continue;
		hs.remove_pair(it->first.first, it->first.second);
	}

	// collect degree
	map<int, int> m;
	for(MPID::iterator it = pe2w.begin(); it != pe2w.end(); it++)
	{
		int e1 = it->first.first;
		int e2 = it->first.second;
		if(m.find(e1) == m.end()) m.insert(PI(e1, 1));
		else m[e1]++;
		if(m.find(e2) == m.end()) m.insert(PI(e2, 1));
		else m[e2]++;
	}

	// decompose
	for(MPID::iterator it = pe2w.begin(); it != pe2w.end(); it++)
	{
		int e1 = it->first.first;
		int e2 = it->first.second;
		double w = it->second;
		assert(w >= cfg.min_guaranteed_edge_weight - SMIN);

		int e = merge_adjacent_edges(e1, e2, w);

		hs.replace(e1, e2, e);

		if(m[e1] == 1) hs.replace(e1, e);
		if(m[e2] == 1) hs.replace(e2, e);
	}

	for(MPID::iterator it = pe2w.begin(); it != pe2w.end(); it++)
	{
		int e1 = it->first.first;
		int e2 = it->first.second;
		
		// test
		if(hs.left_extend(e1) == true && hs.right_extend(e1) == true)
		{
			printf("pe2w:\n");
			for(MPID::iterator ix = pe2w.begin(); ix != pe2w.end(); ix++)
			{
				int ee1 = ix->first.first;
				int ee2 = ix->first.second;
				printf(" (%d, %d) -> %.2lf\n", ee1, ee2, ix->second);
			}

			gr.print();
			hs.print_nodes();
			hs.print_edges();
			fflush(stdout);
		}

		// test
		if(hs.left_extend(e2) == true && hs.right_extend(e2) == true)
		{
			printf("pe2w:\n");
			for(MPID::iterator ix = pe2w.begin(); ix != pe2w.end(); ix++)
			{
				int ee1 = ix->first.first;
				int ee2 = ix->first.second;
				printf(" (%d, %d) -> %.2lf\n", ee1, ee2, ix->second);
			}

			gr.print();
			hs.print_nodes();
			hs.print_edges();
			fflush(stdout);
		}

		//assert(hs.left_extend(e1) == false || hs.right_extend(e1) == false);
		//assert(hs.left_extend(e2) == false || hs.right_extend(e2) == false);
		hs.remove(e1);
		hs.remove(e2);
	}

	assert(gr.degree(root) == 0);
	nonzeroset.erase(root);
	return 0;
}

int scallop::decompose_trivial_vertex(int x)
{
	balance_vertex(x);

	MPID pe2w;
	edge_iterator it1, it2;
	edge_iterator ot1, ot2;
	PEEI pe1, pe2;
	for(pe1 = gr.in_edges(x), it1 = pe1.first, it2 = pe1.second; it1 != it2; it1++)
	{
		int e1 = e2i[*it1];
		double w1 = gr.get_edge_weight(*it1);
		for(pe2 = gr.out_edges(x), ot1 = pe2.first, ot2 = pe2.second; ot1 != ot2; ot1++)
		{
			int e2 = e2i[*ot1];
			double w2 = gr.get_edge_weight(*ot1);
			double w = w1 <= w2 ? w1 : w2;

			pe2w.insert(PPID(PI(e1, e2), w));
		}
	}
	decompose_vertex_replace(x, pe2w);
	return 0;
}

int scallop::classify_trivial_vertex(int x, bool fast)
{
	assert(gr.mixed_strand_vertex(x) == false);
	int d1 = gr.in_degree(x);
	int d2 = gr.out_degree(x);
	if(d1 != 1 && d2 != 1) return -1;

	edge_iterator it1 = gr.in_edges(x).first;
	int e1 = e2i[*it1];
	it1 = gr.out_edges(x).first;
	int e2 = e2i[*it1];

	if(d1 == 1)
	{
		int s = i2e[e1]->source();
		if(gr.out_degree(s) == 1) return 1;
		if(fast && hs.right_dominate(e1) == true) return 1;
	}
	
	if(d2 == 1)
	{
		int t = i2e[e2]->target();
		if(gr.in_degree(t) == 1) return 1;
		if(fast && hs.left_dominate(e2) == true) return 1;
	}

	return 2;
}

int scallop::exchange_sink(int old_sink, int new_sink)
{
	VE ve;
	edge_iterator it1, it2;
	PEEI pei;
	for(pei = gr.in_edges(old_sink), it1 = pei.first, it2 = pei.second; it1 != it2; it1++) ve.push_back(*it1);

	for(int i = 0; i < ve.size(); i++)
	{
		edge_descriptor e = ve[i];
		int s = e->source(); 
		int t = e->target();
		assert(t == old_sink);
		gr.move_edge(e, s, new_sink);
	}
	assert(gr.degree(old_sink) == 0);
	return 0;
}

int scallop::split_merge_path(const VE &p, double wx)
{
	vector<int> v;
	for(int i = 0; i < p.size(); i++)
	{
		assert(p[i] != null_edge);
		assert(e2i.find(p[i]) != e2i.end());
		v.push_back(e2i[p[i]]);
		double w = gr.get_edge_weight(p[i]);
	}
	return split_merge_path(v, wx);
}

int scallop::split_merge_path(const vector<int> &p, double ww)
{
	if(p.size() == 0) return -1;
	int ee = split_edge(p[0], ww);
	for(int i = 1; i < p.size(); i++)
	{
		int x = split_edge(p[i], ww);
		ee = merge_adjacent_equal_edges(ee, x);
	}
	return ee;
}

int scallop::merge_adjacent_equal_edges(int x, int y)
{
	if(i2e[x] == null_edge) return -1;
	if(i2e[y] == null_edge) return -1;

	edge_descriptor xx = i2e[x];
	edge_descriptor yy = i2e[y];

	int xs = (xx)->source();
	int xt = (xx)->target();
	int ys = (yy)->source();
	int yt = (yy)->target();

	if(xt != ys && yt != xs) return -1;
	if(yt == xs) return merge_adjacent_equal_edges(y, x);
	
	assert(xt == ys);
	assert(consistent_strands(x, y));

	edge_descriptor p = gr.add_edge(xs, yt);

	int n = i2e.size();
	i2e.push_back(p);
	assert(e2i.find(p) == e2i.end());
	e2i.insert(PEI(p, n));

	double wx0 = gr.get_edge_weight(xx);
	double wy0 = gr.get_edge_weight(yy);
	assert(fabs(wx0 - wy0) <= SMIN);

	int lx1 = gr.get_edge_info(xx).length;
	int ly1 = gr.get_edge_info(yy).length;
	int lxt = gr.get_vertex_info(xt).length;
	int lxy = lx1 + ly1 + lxt;

	gr.set_edge_weight(p, wx0 * 0.5 + wy0 * 0.5);
	gr.set_edge_info(p, edge_info(lxy));

	borrow_edge_strand(n, x);
	borrow_edge_strand(n, y);

	vector<int> v = mev[xx];
	v.push_back(xt);
	v.insert(v.end(), mev[yy].begin(), mev[yy].end());

	if(mev.find(p) != mev.end()) mev[p] = v;
	else mev.insert(PEV(p, v));

	double sum1 = gr.get_in_weights(xt);
	double sum2 = gr.get_out_weights(xt);
	double sum = (sum1 + sum2) * 0.5;
	double r1 = gr.get_vertex_weight(xt) * (wx0 + wy0) * 0.5 / sum;
	double r2 = gr.get_vertex_weight(xt) - r1;
	gr.set_vertex_weight(xt, r2);

	// set up med/mei
	vertex_info root_info = gr.get_vertex_info(xt);
	int mi = root_info.rpos - root_info.lpos + mei[xx] + mei[yy];
	double md = mi * r1 + med[xx] + med[yy];

	if(med.find(p) != med.end()) med[p] = md;
	else med.insert(PED(p, md));

	if(mei.find(p) != mei.end()) mei[p] = mi;
	else mei.insert(PED(p, mi));

	assert(i2e[n] == p);
	assert(e2i.find(p) != e2i.end());
	assert(e2i[p] == n);
	assert(e2i[i2e[n]] == n);

	remove_edge(x);
	remove_edge(y);

	if(gr.in_degree(xt) == 0 && gr.out_degree(xt) == 0)
	{
		assert(gr.degree(xt) == 0);
		nonzeroset.erase(xt);
	}

	return n;
}

int scallop::remove_edge(int e)
{
	edge_descriptor ee = i2e[e];
	assert(ee != null_edge);
	int s = ee->source();
	int t = ee->target();

	e2i.erase(ee);
	i2e[e] = null_edge;
	gr.remove_edge(ee);

	return 0;
}

int scallop::merge_adjacent_edges(int x, int y, double ww)
{
	assert(ww >= cfg.min_guaranteed_edge_weight - SMIN);
	if(i2e[x] == null_edge) return -1;
	if(i2e[y] == null_edge) return -1;

	edge_descriptor xx = i2e[x];
	edge_descriptor yy = i2e[y];

	int xs = xx->source();
	int xt = xx->target();
	int ys = yy->source();
	int yt = yy->target();

	if(xt != ys) return merge_adjacent_edges(y, x, ww);
	assert(xt == ys);

	int x1 = split_edge(x, ww);
	int y1 = split_edge(y, ww);
	int xy = merge_adjacent_equal_edges(x1, y1);

	return xy;
}

int scallop::merge_adjacent_edges(int x, int y)
{
	if(i2e[x] == null_edge) return -1;
	if(i2e[y] == null_edge) return -1;

	edge_descriptor xx = i2e[x];
	edge_descriptor yy = i2e[y];

	double wx = gr.get_edge_weight(xx);
	double wy = gr.get_edge_weight(yy);
	double ww = (wx <= wy) ? wx : wy;

	return merge_adjacent_edges(x, y, ww);
}

int scallop::split_edge(int ei, double w)
{
	assert(w >= cfg.min_guaranteed_edge_weight - SMIN);
	assert(i2e[ei] != null_edge);
	edge_descriptor ee = i2e[ei];

	double ww = gr.get_edge_weight(ee);
	if(fabs(ww - w) <= SMIN) return ei;

	/*
	if(ww <= w + cfg.min_guaranteed_edge_weight)
	{
		gr.set_edge_weight(ee, w);
		return ei;
	}
	*/

	int s = ee->source();
	int t = ee->target();

	edge_descriptor p2 = gr.add_edge(s, t);
	edge_info eif = gr.get_edge_info(ee);

	double www = ww - w;
	if(www <= cfg.min_guaranteed_edge_weight) www = cfg.min_guaranteed_edge_weight;

	gr.set_edge_weight(ee, www);		// old edge
	gr.set_edge_info(ee, eif);			// old edge
	gr.set_edge_weight(p2, w);			// new edge
	gr.set_edge_info(p2, eif);			// new edge

	if(mev.find(p2) != mev.end()) mev[p2] = mev[ee];
	else mev.insert(PEV(p2, mev[ee]));

	// set up med/mdi
	int mi = mei[ee];
	double md = med[ee] * w / ww;

	if(med.find(p2) != med.end()) med[p2] = md;
	else med.insert(PED(p2, md));

	if(mei.find(p2) != mei.end()) mei[p2] = mi;
	else mei.insert(PEI(p2, mi));

	int n = i2e.size();
	i2e.push_back(p2);
	e2i.insert(PEI(p2, n));

	return n;
}

int scallop::balance_vertex(int v, const vector<int> &ve1, const vector<int> &ve2)
{
	if(gr.degree(v) <= 0) return 0;
	if(ve1.size() <= 0) return 0;
	if(ve2.size() <= 0) return 0;

	vector<double> v1;
	vector<double> v2;
	double w1 = 0, w2 = 0;
	for(int i = 0; i < ve1.size(); i++)
	{
		double w = gr.get_edge_weight(i2e[ve1[i]]);
		assert(w >= cfg.min_guaranteed_edge_weight - SMIN);
		v1.push_back(w);
		w1 += w;
	}
	for(int i = 0; i < ve2.size(); i++)
	{
		double w = gr.get_edge_weight(i2e[ve2[i]]);
		assert(w >= cfg.min_guaranteed_edge_weight - SMIN);
		v2.push_back(w);
		w2 += w;
	}

	// use sqrt-meature
	double ww = sqrt(w1 * w2);

	double r1 = ww / w1;
	double r2 = ww / w2;

	double m1 = 0, m2 = 0;
	for(int i = 0; i < v1.size(); i++)
	{
		edge_descriptor e = i2e[ve1[i]];
		double wx = gr.get_edge_weight(e);
		double wy = wx * r1;
		if(wy < cfg.min_guaranteed_edge_weight)
		{
			m1 += cfg.min_guaranteed_edge_weight - wy;
			wy = cfg.min_guaranteed_edge_weight;
		}
		gr.set_edge_weight(e, wy);
	}
	for(int i = 0; i < v2.size(); i++)
	{
		edge_descriptor e = i2e[ve2[i]];
		double wx = gr.get_edge_weight(e);
		double wy = wx * r2;
		if(wy < cfg.min_guaranteed_edge_weight)
		{
			m2 += cfg.min_guaranteed_edge_weight - wy;
			wy = cfg.min_guaranteed_edge_weight;
		}
		gr.set_edge_weight(e, wy);
	}

	if(m1 > m2)
	{
		edge_descriptor e = i2e[ve2.front()];
		double w = gr.get_edge_weight(e);
		gr.set_edge_weight(e, w + m1 - m2);
	}
	else if(m1 < m2)
	{
		edge_descriptor e = i2e[ve1.front()];
		double w = gr.get_edge_weight(e);
		gr.set_edge_weight(e, w + m2 - m1);
	}

	return 0;
}

int scallop::balance_vertex(int v)
{
	if(gr.in_degree(v) <= 0) return 0;
	if(gr.out_degree(v) <= 0) return 0;

	PEEI pei;
	edge_iterator it1, it2;
	vector<int> ve1;
	vector<int> ve2;
	for(pei = gr.in_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		ve1.push_back(e2i[*it1]);
	}
	for(pei = gr.out_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		ve2.push_back(e2i[*it1]);
	}
	return balance_vertex(v, ve1, ve2);
}

double scallop::compute_balance_ratio(int v)
{
	edge_iterator it1, it2;
	PEEI pei;
	double w1 = 0, w2 = 0;
	for(pei = gr.in_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		w1 += w;
	}
	for(pei = gr.out_edges(v), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		w2 += w;
	}

	assert(w1 >= SMIN);
	assert(w2 >= SMIN);

	//double ww = sqrt(0.5 * w1 * w1 + 0.5 * w2 * w2);
	double ww = 0.5 * w1 + 0.5 * w2;

	if(w1 >= w2) return w1 / w2;
	else return w2 / w1;
}

int scallop::split_vertex(int x, const vector<int> &xe, const vector<int> &ye)
{
	assert(x != 0);
	assert(x != gr.num_vertices() - 1);
	if(xe.size() <= 0) return 0;
	if(ye.size() <= 0) return 0;

	double w1 = 0, w2 = 0;
	for(int i = 0; i < xe.size(); i++) w1 += gr.get_edge_weight(i2e[xe[i]]);
	for(int i = 0; i < ye.size(); i++) w2 += gr.get_edge_weight(i2e[ye[i]]);
	double r1 = w1 / gr.get_in_weights(x);
	double r2 = w2 / gr.get_out_weights(x);
	double ww1 = (r1 + r2) * 0.5 * gr.get_vertex_weight(x);
	double ww2 = gr.get_vertex_weight(x) - ww1;

	int n = gr.num_vertices();
	assert(v2v.size() == n);

	// vertex-n => new sink vertex
	// vertex-(n-1) => splitted vertex for xe and ye
	// vertex-x => splitted vertex for xe2 and ye2

	gr.add_vertex();
	assert(nonzeroset.find(n - 1) == nonzeroset.end());
	nonzeroset.insert(n - 1);
	gr.set_vertex_info(n, gr.get_vertex_info(n - 1));
	gr.set_vertex_info(n - 1, gr.get_vertex_info(x));
	gr.set_vertex_weight(n, gr.get_vertex_weight(n - 1));
	gr.set_vertex_weight(n - 1, ww1);
	gr.set_vertex_weight(x, ww2);

	v2v.push_back(v2v[n - 1]);
	v2v[n - 1] = v2v[x];

	// use vertex-n instead of vertex-(n-1) as sink vertex
	exchange_sink(n - 1, n);

	// revise phasing paths
	set<int> sx(xe.begin(), xe.end());
	set<int> sy(ye.begin(), ye.end());
	PEEI pi = gr.in_edges(x);
	PEEI po = gr.out_edges(x);
	edge_iterator it1, it2;
	for(it1 = pi.first; it1 != pi.second; it1++)
	{
		int e1 = e2i[*it1];
		for(it2 = po.first; it2 != po.second; it2++)
		{
			int e2 = e2i[*it2];
			if(sx.find(e1) == sx.end() && sy.find(e2) == sy.end()) continue;
			if(sx.find(e1) != sx.end() && sy.find(e2) != sy.end()) continue;
			hs.remove_pair(e1, e2);
		}
	}

	// attach edges in xe and ye to vertex-(n-1)
	for(int i = 0; i < xe.size(); i++)
	{
		edge_descriptor e = i2e[xe[i]];
		assert(e != null_edge);
		int s = e->source();
		int t = e->target();
		assert(t == x);
		gr.move_edge(e, s, n - 1);
	}
	for(int i = 0; i < ye.size(); i++)
	{
		edge_descriptor e = i2e[ye[i]];
		assert(e != null_edge);
		int s = e->source();
		int t = e->target();
		assert(s == x);
		gr.move_edge(e, n - 1, t);
	}

	return 0;
}

vector<int> scallop::topological_sort()
{
	vector<PI> v;
	for(int i = 0; i < v2v.size(); i++)
	{
		v.push_back(PI(v2v[i], i));
	}
	sort(v.begin(), v.end());

	vector<int> vv;
	for(int i = 0; i < v.size(); i++)
	{
		vv.push_back(v[i].second);
	}

	return vv;
}

int scallop::collect_existing_st_paths()
{
	for(int i = 0; i < i2e.size(); i++)
	{
		if(i2e[i] == null_edge) continue;
		if(i2e[i]->source() != 0) continue;
		if(i2e[i]->target() != gr.num_vertices() - 1) continue;
		collect_path(i);
	}
	return 0;
}

int scallop::collect_phasing_paths()
{
	for(int i = 0; i < i2e.size(); i++)
	{
		if(i2e[i] == null_edge) continue;
		int s = i2e[i]->source();
		int t = i2e[i]->target();
		collect_phasing_path(i, s, t);
	}
	return 0;
}

int scallop::collect_path(int e)
{
	assert(mev.find(i2e[e]) != mev.end());

	vector<int> v0 = mev[i2e[e]];
	vector<int> v;

	int mi = 0;
	for(int i = 0; i < v0.size(); i++) 
	{
		if(v2v[v0[i]] < 0) continue;
		v.push_back(v2v[v0[i]]);
		vertex_info vi = gr.get_vertex_info(v2v[v0[i]]);
		mi += vi.rpos - vi.lpos;
	}

	assert(mei[i2e[e]] == mi);

	sort(v.begin(), v.end());

	int n = v2v[gr.num_vertices() - 1];
	assert(v[0] > 0);
	assert(v[v.size() - 1] < n);
	v.insert(v.begin(), 0);
	v.push_back(n);

	path p;
	p.length = mi;
	p.abd = gr.get_edge_weight(i2e[e]);
	p.reads = med[i2e[e]];
	//p.abd = med[i2e[e]] / mi;
	//p.reads = gr.get_edge_weight(i2e[e]);
	p.v = v;
	if(gr.get_edge_info(i2e[e]).strand == 1) p.strand = '+';
	if(gr.get_edge_info(i2e[e]).strand == 2) p.strand = '-';
	if(p.strand == '.') p.strand = gr.strand;
	paths.push_back(p);

	gr.remove_edge(i2e[e]);
	e2i.erase(i2e[e]);
	i2e[e] = null_edge;

	return 0;
}

int scallop::collect_phasing_path(int e, int s, int t)
{
	assert(mev.find(i2e[e]) != mev.end());

	vector<int> v0 = mev[i2e[e]];
	vector<int> v;
	for(int i = 0; i < v0.size(); i++) 
	{
		if(v2v[v0[i]] < 0) continue;
		v.push_back(v2v[v0[i]]);
	}

	sort(v.begin(), v.end());

	/*
	printf("prepare to collect: ");
	printv(v);
	printf(", s = %d, t = %d, v2v[s] = %d, v2v[t] = %d\n", s, t, v2v[s], v2v[t]);
	*/

	if(v2v[t] >= 0) v.push_back(v2v[t]);

	path p;
	p.abd = gr.get_edge_weight(i2e[e]);
	p.v = v;
	if(gr.get_edge_info(i2e[e]).strand == 1) p.strand = '+';
	if(gr.get_edge_info(i2e[e]).strand == 2) p.strand = '-';
	if(p.strand == '.') p.strand = gr.strand;

	paths.push_back(p);

	gr.remove_edge(i2e[e]);
	e2i.erase(i2e[e]);
	i2e[e] = null_edge;

	return 0;
}

int scallop::greedy_decompose()
{
	if(gr.num_edges() == 0) return 0;

	for(int i = 1; i < gr.num_vertices() - 1; i++) balance_vertex(i);
	for(int i = 1; i < gr.num_vertices() - 1; i++) balance_vertex(i);

	int cnt = 0;
	int n1 = paths.size();
	while(true)
	{
		VE v;
		double w = gr.compute_maximum_path_w(v);
		if(w < 0) break;
		if(w <= cfg.min_transcript_coverage) break;

		int e = split_merge_path(v, w);
		collect_path(e);
		cnt++;
	}
	int n2 = paths.size();
	if(cfg.verbose >= 2) printf("greedy decomposing produces %d / %d paths\n", n2 - n1, n2);
	return 0;
}

int scallop::compute_smallest_edge(int x, double &ratio)
{
	int e = -1;
	ratio = DBL_MAX;
	edge_iterator it1, it2;
	PEEI pei;
	double sum1 = 0;
	double sum2 = 0;
	for(pei = gr.in_edges(x), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		sum1 += w;
	}
	for(pei = gr.out_edges(x), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		sum2 += w;
	}

	assert(sum1 >= SMIN);
	assert(sum2 >= SMIN);
	for(pei = gr.in_edges(x), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		double r = w * 1.0 / sum1;
		if(r >= ratio) continue;
		ratio = r;
		e = e2i[*it1];
	}
	for(pei = gr.out_edges(x), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		double w = gr.get_edge_weight(*it1);
		double r = w * 1.0 / sum2;
		if(r >= ratio) continue;
		ratio = r;
		e = e2i[*it1];
	}
	assert(e >= 0);
	return e;
}

int scallop::print()
{
	/*
	int n0 = 0;
	int n1 = 0;
	for(int i = 1; i < gr.num_vertices() - 1; i++) 
	{
		if(gr.degree(i) == 0) n0++;
		if(gr.degree(i) >= 1) n1++;
	}
	assert(nonzeroset.size() == n1);
	*/

	printf("statistics: %lu edges, %lu vertices, %lu nonzero vertices\n", gr.num_edges(), gr.num_vertices(), nonzeroset.size());

	//int p1 = gr.compute_num_paths();
	//int p2 = gr.compute_decomp_paths();
	//printf("statistics: %lu edges, %d vertices, total %d paths, %d required\n", gr.num_edges(), n, p1, p2);

	//printf("finish round %d\n\n", round);
	//round++;

	return 0;
}

int scallop::print_super_edges()
{
	for(int i = 0; i < i2e.size(); i++)
	{
		if(i2e[i] == null_edge) continue;
		edge_descriptor e = i2e[i];
		assert(mev.find(e) != mev.end());
		vector<int> &v = mev[e];
		double w = gr.get_edge_weight(e);
		double r = med[e];
		int l = mei[e];
		printf("super-edge %d (%d, %d), w = %.2lf, r = %.2lf, len = %d, list = ( ", i, e->source(), e->target(), w, r, l);
		printv(v);
		printf(")\n");
	}
	return 0;
}

int scallop::print_phasing_paths(hyper_set &hh)
{
	for(int k = 0; k < hh.edges.size(); k++)
	{
		vector<int> &v = hh.edges[k];
		int c = hh.ecnts[k];
		printf("phase %d: count = %d, list = ", k, c);
		for(int i = 0; i < v.size(); i++)
		{
			int e = v[i];
			if(e < 0 || i2e[e] == null_edge) printf("%d (-1, -1) -> ", e);
			else printf("%d (%d, %d) -> ", e, i2e[e]->source(), i2e[e]->target());

			if(i + 1 < v.size() && v[i + 1] > 0 && e > 0)
			{
				int e1 = v[i + 1];

				assert(i2e[e] != null_edge);

				//if(i2e[e1] == null_edge) printf("*****killer = %d\n", e1);
				assert(i2e[e1] != null_edge);

				//printf("e = %d (%d, %d), e1 = %d (%d, %d)\n", e, i2e[e]->source(), i2e[e]->target(), e1, i2e[e1]->source(), i2e[e1]->target());

				assert(i2e[e]->target() == i2e[e1]->source());
			}
		}
		printf("end\n");
	}
	return 0;
}


int scallop::draw_splice_graph(const string &file) 
{
	MIS mis;
	char buf[10240];
	for(int i = 0; i < gr.num_vertices(); i++)
	{
		vertex_info vi = gr.get_vertex_info(i);
		double w = gr.get_vertex_weight(i);
		int l = vi.length;
		//string s = gr.get_vertex_string(i);
		//sprintf(buf, "%d:%.0lf:%s", i, w, s.c_str());
		sprintf(buf, "%d:%.0lf:%d", i, w, l);
		mis.insert(PIS(i, buf));
	}

	MES mes;
	for(int i = 0; i < i2e.size(); i++)
	{
		if(i2e[i] == null_edge) continue;
		double w = gr.get_edge_weight(i2e[i]);
		edge_info ei = gr.get_edge_info(i2e[i]);
		int l = ei.length;
		sprintf(buf, "%d:%.0lf", i, w);
		mes.insert(PES(i2e[i], buf));
	}
	
	vector<int> tp = topological_sort();
	gr.draw(file, mis, mes, 4.5, tp);
	return 0;
}

int scallop::build_transcripts()
{
	trsts.clear();
	for(int i = 0; i < paths.size(); i++)
	{
		string tid = gr.gid + "." + tostring(i);
		transcript trst;
		path &p = paths[i];
		build_transcript(gr, trst, p.v, p.strand, p.abd, tid);
		trsts.push_back(trst);
	}
	return 0;
}
