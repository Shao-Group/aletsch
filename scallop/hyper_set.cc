/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#include "essential.h"
#include "constants.h"
#include "hyper_set.h"
#include "parameters.h"
#include <algorithm>
#include <cstdio>

hyper_set::hyper_set()
{}

hyper_set::hyper_set(splice_graph &gr, const phase_set &ps)
{
	for(MVII::const_iterator it = ps.pmap.begin(); it != ps.pmap.end(); it++)
	{
		const vector<int32_t> &v = it->first;
		int c = it->second;
		vector<int> vv;
		bool b = build_path_from_exon_coordinates(gr, v, vv);
		if(b == false) continue;
		for(int k = 0; k < vv.size(); k++) vv[k]--;
		add_node_list(vv, c);
	}
}

int hyper_set::clear()
{
	nodes.clear();
	edges.clear();
	e2s.clear();
	ecnts.clear();
	return 0;
}

int hyper_set::add_node_list(const set<int> &s)
{
	return add_node_list(s, 1);
}

int hyper_set::add_node_list(const set<int> &s, int c)
{
	vector<int> v(s.begin(), s.end());
	return add_node_list(v, c);
}

int hyper_set::add_node_list(const vector<int> &s, int c)
{
	vector<int> v = s;
	sort(v.begin(), v.end());
	for(int i = 0; i < v.size(); i++) v[i]++;
	if(nodes.find(v) == nodes.end()) nodes.insert(PVII(v, c));
	else nodes[v] += c;
	return 0;
}

int hyper_set::merge_node_list(const vector<int> &s, int c)
{
	vector<int> v = s;
	sort(v.begin(), v.end());
	bool useful = false;
	for(MVII::iterator it = nodes.begin(); it != nodes.end(); it++)
	{
		int z = compare_phasing_paths(it->first, v);
		if(z == 0) return 0;
		if(z == 1) continue;
		if(z == 2) continue;
		if(z == 3) return 0;
		if(z == 4) useful = true;
		if(z == 5) useful = true;
		if(z == 6) useful = true;
		if(z == 7) continue;
		if(z == 8) continue;
		if(z == 9) continue;	// TODO
	}

	if(useful == false) return 0;

	printf("add extra phasing path, c = %d, list = ( ", c);
	printv(v);
	printf(")\n");

	assert(nodes.find(v) == nodes.end());
	nodes.insert(PVII(v, c));

	return 0;
}

int hyper_set::compare(const hyper_set &hx)
{	
	for(MVII::const_iterator it = nodes.begin(); it != nodes.end(); it++)
	{
		const vector<int> &x = it->first;
		for(MVII::const_iterator ti = hx.nodes.begin(); ti != hx.nodes.end(); ti++)
		{
			const vector<int> &y = ti->first;
			int z = compare_phasing_paths(x, y);
			printf("compare z = %d, label = %s: ( ", z, position_names[z].c_str());
			printv(x);
			printf(") vs ( ");
			printv(y);
			printf(")\n");
		}
	}
	return 0;
}

int hyper_set::extend(const hyper_set &hx)
{	
	vector<PVII> vx(nodes.begin(), nodes.end());
	vector<PVII> vy(hx.nodes.begin(), hx.nodes.end());

	vector< vector<int> > cxx;
	vector< vector<int> > cxy;

	cxx.resize(vx.size());
	for(int i = 0; i < vx.size(); i++)
	{
		cxx[i].assign(vx.size(), -1);
		for(int j = 0; j < vx.size(); j++)
		{
			cxx[i][j] = compare_phasing_paths(vx[i].first, vx[j].first);
		}
	}

	cxy.resize(vx.size());
	for(int i = 0; i < vx.size(); i++)
	{
		cxy[i].assign(vy.size(), -1);
		for(int j = 0; j < vy.size(); j++)
		{
			cxy[i][j] = compare_phasing_paths(vx[i].first, vy[j].first);
		}
	}

	int added = 0;
	for(int i = 0; i < vx.size(); i++)
	{
		bool extend_left = false;
		for(int j = 0; j < vx.size(); j++)
		{
			if(cxx[i][j] == EXTEND_LEFT) extend_left = true;
			if(extend_left == true) break;
		}
		if(extend_left == true) continue;


		int bestj = -1;
		int besto = -1;
		int bestw = -1;
		for(int j = 0; j < vy.size(); j++)
		{
			if(cxy[i][j] != EXTEND_LEFT && cxy[i][j] != CONTAINING) continue;

			vector<int>::const_iterator r1 = lower_bound(vy[j].first.begin(), vy[j].first.end(), vx[i].first.front());
			if(r1 == vy[j].first.begin()) continue;
			//if(vy[j].first.front() == 0) continue;		// additional constraint

			int overlap = vy[j].first.end() - r1;
			//if(overlap > vx[i].first.size()) continue;
			if(overlap > vx[i].first.size()) overlap = vx[i].first.size();

			if(overlap < besto) continue;
			if(overlap == besto && vy[j].second < bestw) continue;

			bestj = j;
			besto = overlap;
			bestw = vy[j].second;
		}

		if(bestj == -1) continue;
		int j = bestj;

		printf("phasing path %d: c = %d, ( ", i, vx[i].second);
		printv(vx[i].first);
		printf("), ");

		printf(" best-left-extend %d: c = %d, overlap = %d, ( ", j, bestw, besto);
		printv(vy[j].first);
		printf(")\n");

		if(nodes.find(vy[j].first) == nodes.end()) 
		{
			added++;
			nodes.insert(vy[j]);
		}
	}

	for(int i = 0; i < vx.size(); i++)
	{
		bool extend_right = false;
		for(int j = 0; j < vx.size(); j++)
		{
			if(cxx[i][j] == EXTEND_RIGHT) extend_right = true;
			if(extend_right == true) break;
		}
		if(extend_right == true) continue;

		int bestj = -1;
		int besto = -1;
		int bestw = -1;
		for(int j = 0; j < vy.size(); j++)
		{
			if(cxy[i][j] != EXTEND_RIGHT && cxy[i][j] != CONTAINING) continue;

			vector<int>::const_iterator r1 = lower_bound(vy[j].first.begin(), vy[j].first.end(), vx[i].first.back());
			if(r1 == vy[j].first.end() - 1) continue;

			int overlap = r1 - vy[j].first.begin() + 1;
			//if(overlap > vx[i].first.size()) continue;
			if(overlap > vx[i].first.size()) overlap = vx[i].first.size();

			if(overlap < besto) continue;
			if(overlap == besto && vy[j].second < bestw) continue;

			bestj = j;
			besto = overlap;
			bestw = vy[j].second;
		}

		if(bestj == -1) continue;

		int j = bestj;

		printf("phasing path %d: c = %d, ( ", i, vx[i].second);
		printv(vx[i].first);
		printf("), ");

		printf(" best-right-extend %d: c = %d, overlap = %d, ( ", j, bestw, besto);
		printv(vy[j].first);
		printf(")\n");

		if(nodes.find(vy[j].first) == nodes.end()) 
		{
			added++;
			nodes.insert(vy[j]);
		}
	}

	printf("summary of extend: hs = %lu, hx = %lu, added = %d\n", vx.size(), vy.size(), added);

	return 0;
}

int hyper_set::merge(const hyper_set &hx)
{
	for(MVII::const_iterator it = hx.nodes.begin(); it != hx.nodes.end(); it++)
	{
		const vector<int> &x = it->first;
		int c = it->second;
		merge_node_list(x, c);
	}
	return 0;
}

int hyper_set::build(directed_graph &gr, MEI& e2i)
{
	build_edges(gr, e2i);
	build_index();
	return 0;
}

int hyper_set::build_edges(directed_graph &gr, MEI& e2i)
{
	edges.clear();
	for(MVII::iterator it = nodes.begin(); it != nodes.end(); it++)
	{
		int c = it->second;
		//if(c < min_hyper_count) continue;

		const vector<int> &vv = it->first;
		if(vv.size() <= 1) continue;

		vector<int> ve;
		bool b = true;
		for(int k = 0; k < vv.size() - 1; k++)
		{
			PEB p = gr.edge(vv[k], vv[k + 1]);
			if(p.second == false) b = false;
			if(p.second == false) ve.push_back(-1);
			else ve.push_back(e2i[p.first]);
		}

		if(b == true && ve.size() >= 2)
		{
			edges.push_back(ve);
			ecnts.push_back(c);
		}
		continue;
	}
	return 0;
}

int hyper_set::filter_nodes(splice_graph &gr)
{
	MVII mv;
	for(MVII::iterator it = nodes.begin(); it != nodes.end(); it++)
	{
		int c = it->second;
		const vector<int> &vv = it->first;
		if(vv.size() <= 1) continue;

		bool b = check_valid_path(gr, vv);
		if(b == false) continue;
		mv.insert(PVII(it->first, it->second));
	}
	nodes = mv;
	return 0;
}

int hyper_set::filter()
{
	vector<PVII> vv(nodes.begin(), nodes.end());
	vector<bool> bb(vv.size(), false);

	// build index with the first element
	map< int, vector<int> > m;
	for(int i = 0; i < vv.size(); i++)
	{
		if(vv[i].first.size() <= 0) continue;
		int s = vv[i].first.front();

		if(m.find(s) == m.end())
		{
			vector<int> v;
			v.push_back(i);
			m.insert(pair<int, vector<int> >(s, v));
		}
		else
		{
			m[s].push_back(i);
		}
	}

	for(int i = 0; i < vv.size(); i++)
	{
		if(bb[i] == true) continue;
		vector<int> &v = vv[i].first;

		for(int a = 0; a < v.size(); a++)
		{
			int s = v[a];
			if(m.find(s) == m.end()) continue;

			vector<int> &x = m[s];
			for(int k = 0; k < x.size(); k++)
			{
				int j = x[k];
				vector<int> &u = vv[j].first;

				// check whether i covers j
				assert(u.front() == s);
				if(j == i) continue;
				if(v.size() - a < u.size()) continue;
				bool b = identical(v, a, a + u.size() - 1, u, 0, u.size() - 1);
				if(b == false) continue;

				bb[j] = true;
				continue;
			}
		}
	}

	nodes.clear();
	for(int i = 0; i < vv.size(); i++)
	{
		if(bb[i] == true) continue;
		nodes.insert(vv[i]);
	}

	return 0;
}

int hyper_set::build_index()
{
	e2s.clear();
	for(int i = 0; i < edges.size(); i++)
	{
		vector<int> &v = edges[i];
		for(int j = 0; j < v.size(); j++)
		{
			int e = v[j];
			if(e == -1) continue;
			if(e2s.find(e) == e2s.end())
			{
				set<int> s;
				s.insert(i);
				e2s.insert(PISI(e, s));
			}
			else
			{
				e2s[e].insert(i);
			}
		}
	}
	return 0;
}

int hyper_set::update_index()
{
	vector<int> fb1;
	for(MISI::iterator p = e2s.begin(); p != e2s.end(); p++)
	{
		int e = p->first;
		set<int> &ss = p->second;
		vector<int> fb2;
		for(set<int>::iterator it = ss.begin(); it != ss.end(); it++)
		{
			vector<int> &v = edges[*it];
			for(int i = 0; i < v.size(); i++)
			{
				if(v[i] != e) continue;
				bool b1 = false, b2 = false;
				if(i == 0 || v[i - 1] == -1) b1 = true;
				if(i == v.size() - 1 || v[i + 1] == -1) b2 = true;
				if(b1 == true && b2 == true) fb2.push_back(*it);
				break;
			}
		}
		for(int i = 0; i < fb2.size(); i++) ss.erase(fb2[i]);
		if(ss.size() == 0) fb1.push_back(e);
	}
	for(int i = 0; i < fb1.size(); i++) e2s.erase(fb1[i]);
	return 0;
}

set<int> hyper_set::get_intersection(const vector<int> &v)
{
	set<int> ss;
	if(v.size() == 0) return ss;
	assert(v[0] >= 0);
	if(e2s.find(v[0]) == e2s.end()) return ss;
	ss = e2s[v[0]];
	vector<int> vv(ss.size());
	for(int i = 1; i < v.size(); i++)
	{
		assert(v[i] >= 0);
		set<int> s;
		if(e2s.find(v[i]) == e2s.end()) return s;
		s = e2s[v[i]];
		vector<int>::iterator it = set_intersection(ss.begin(), ss.end(), s.begin(), s.end(), vv.begin());
		ss = set<int>(vv.begin(), it);
	}
	return ss;
}

MI hyper_set::get_successors(int e)
{
	MI s;
	if(e2s.find(e) == e2s.end()) return s;
	set<int> &ss = e2s[e];
	for(set<int>::iterator it = ss.begin(); it != ss.end(); it++)
	{
		vector<int> &v = edges[*it];
		int c = ecnts[*it];
		for(int i = 0; i < v.size(); i++)
		{
			if(v[i] != e) continue;
			if(i >= v.size() - 1) continue;
			int k = v[i + 1];
			if(k == -1) continue;
			if(s.find(k) == s.end()) s.insert(PI(k, c));
			else s[k] += c;
		}
	}
	return s;
}

MI hyper_set::get_predecessors(int e)
{
	MI s;
	if(e2s.find(e) == e2s.end()) return s;
	set<int> &ss = e2s[e];
	for(set<int>::iterator it = ss.begin(); it != ss.end(); it++)
	{
		vector<int> &v = edges[*it];
		int c = ecnts[*it];
		for(int i = 0; i < v.size(); i++)
		{
			if(v[i] != e) continue;
			if(i == 0) continue;
			int k = v[i - 1];
			if(k == -1) continue;
			if(s.find(k) == s.end()) s.insert(PI(k, c));
			else s[k] += c;
		}
	}
	return s;
}

MPII hyper_set::get_routes(int x, directed_graph &gr, MEI &e2i)
{
	MPII mpi;
	edge_iterator it1, it2;
	PEEI pei;
	vector<PI> v;
	for(pei = gr.in_edges(x), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		assert(e2i.find(*it1) != e2i.end());
		int e = e2i[*it1];
		MI s = get_successors(e);
		for(MI::iterator it = s.begin(); it != s.end(); it++)
		{
			PI p(e, it->first);
			mpi.insert(PPII(p, it->second));
		}
	}
	return mpi;
}

/*
int hyper_set::get_routes(int x, directed_graph &gr, MEI &e2i, MPII &mpi)
{
	edge_iterator it1, it2;
	mpi.clear();
	int total = 0;
	for(tie(it1, it2) = gr.in_edges(x); it1 != it2; it1++)
	{
		assert(e2i.find(*it1) != e2i.end());
		int e = e2i[*it1];

		if(e2s.find(e) == e2s.end()) continue;
		set<int> &ss = e2s[e];
		for(set<int>::iterator it = ss.begin(); it != ss.end(); it++)
		{
			int k = *it;
			assert(k >= 0 && k < edges.size());
			assert(k >= 0 && k < ecnts.size());
			vector<int> &v = edges[k];
			int cnt = ecnts[k];
			for(int i = 0; i < v.size(); i++)
			{
				if(v[i] != e) continue;
				if(i == v.size() - 1) continue;
				if(v[i + 1] == -1) continue;
				PI p(e, v[i + 1]);
				total += cnt;
				if(mpi.find(p) != mpi.end()) mpi[p] += cnt;
				else mpi.insert(PPII(p, cnt));
			}
		}
	}
	return total;
}
*/

int hyper_set::replace(int x, int e)
{
	vector<int> v;
	v.push_back(x);
	replace(v, e);
	return 0;
}

int hyper_set::replace(int x, int y, int e)
{
	vector<int> v;
	v.push_back(x);
	v.push_back(y);
	replace(v, e);
	return 0;
}

int hyper_set::replace(const vector<int> &v, int e)
{
	if(v.size() == 0) return 0;
	set<int> s = get_intersection(v);

	vector<int> fb;
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		vector<int> bv = consecutive_subset(vv, v);

		if(bv.size() <= 0) continue;

		//assert(bv.size() == 1); 

		sort(bv.begin(), bv.end());
		for(int j = bv.size() - 1; j >= 0; j--)
		{
			int b = bv[j];
			vv[b] = e;

			vv.erase(vv.begin() + b + 1, vv.begin() + b + v.size());
		}

		fb.push_back(k);

		if(e2s.find(e) == e2s.end())
		{
			set<int> ss;
			ss.insert(k);
			e2s.insert(PISI(e, ss));
		}
		else
		{
			e2s[e].insert(k);
		}
	}

	if(v.size() != 1) return 0;

	for(int i = 0; i < v.size(); i++)
	{
		int u = v[i];
		if(e2s.find(u) == e2s.end()) continue;
		for(int k = 0; k < fb.size(); k++) e2s[u].erase(fb[k]);
		if(e2s[u].size() == 0) e2s.erase(u);
	}
	return 0;
}

int hyper_set::replace_strange(const vector<int> &v, int e)
{
	if(v.size() == 0) return 0;
	set<int> s = get_intersection(v);
	
	vector<int> fb;
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		vector<int> bv = consecutive_subset(vv, v);

		if(bv.size() <= 0) continue;
		assert(bv.size() == 1); 

		int b = bv[0];
		vv[b] = e;

		bool b1 = useful(vv, 0, b);
		bool b2 = useful(vv, b + v.size() - 1, vv.size() - 1);

		if(b1 == false && b2 == false)
		{
			fb.push_back(k);
			continue;
		}

		vv.erase(vv.begin() + b + 1, vv.begin() + b + v.size());

		if(e2s.find(e) == e2s.end())
		{
			set<int> ss;
			ss.insert(k);
			e2s.insert(PISI(e, ss));
		}
		else
		{
			e2s[e].insert(k);
		}
	}

	for(int i = 0; i < v.size(); i++)
	{
		int u = v[i];
		if(e2s.find(u) == e2s.end()) continue;
		for(int k = 0; k < fb.size(); k++) e2s[u].erase(fb[k]);
		if(e2s[u].size() == 0) e2s.erase(u);
	}
	return 0;
}

int hyper_set::remove(const set<int> &s)
{
	return remove(vector<int>(s.begin(), s.end()));
}

int hyper_set::remove(const vector<int> &v)
{
	for(int i = 0; i < v.size(); i++) remove(v[i]);
	return 0;
}

int hyper_set::remove(int e)
{
	if(e2s.find(e) == e2s.end()) return 0;
	set<int> &s = e2s[e];
	vector<int> fb;
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		assert(vv.size() >= 1);

		for(int i = 0; i < vv.size(); i++)
		{
			if(vv[i] != e) continue;

			vv[i] = -1;

			/*
			bool b1 = useful(vv, 0, i - 1);
			bool b2 = useful(vv, i + 1, vv.size() - 1);
			if(b1 == false && b2 == false) fb.push_back(k);
			*/
			 
			fb.push_back(k);
			//break;
		}
	}

	for(int i = 0; i < fb.size(); i++) s.erase(fb[i]);
	e2s.erase(e);
	return 0;
}

int hyper_set::remove_pair(int x, int y)
{
	insert_between(x, y, -1);
	return 0;

	if(e2s.find(x) == e2s.end()) return 0;
	set<int> &s = e2s[x];
	vector<int> fb;
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		assert(vv.size() >= 1);

		for(int i = 0; i < vv.size(); i++)
		{
			if(i == vv.size() - 1) continue;
			if(vv[i] != x) continue;
			if(vv[i + 1] != y) continue;

			bool b1 = useful(vv, 0, i);
			bool b2 = (b1 == true) ? true : useful(vv, i + 1, vv.size() - 1);

			if(b1 == false && b2 == false) fb.push_back(k);
			else vv.insert(vv.begin() + i + 1, -1);

			break;
		}
	}

	for(int i = 0; i < fb.size(); i++) s.erase(fb[i]);
	if(s.size() == 0) e2s.erase(x);

	return 0;
}

bool hyper_set::useful(const vector<int> &v, int k1, int k2)
{
	for(int i = k1; i < k2; i++)
	{
		if(v[i] >= 0 && v[i + 1] >= 0) return true;
	}
	return false;
}

int hyper_set::insert_between(int x, int y, int e)
{
	if(e2s.find(x) == e2s.end()) return 0;
	set<int> s = e2s[x];
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		assert(vv.size() >= 1);

		for(int i = 0; i < vv.size(); i++)
		{
			if(i == vv.size() - 1) continue;
			if(vv[i] != x) continue;
			if(vv[i + 1] != y) continue;
			vv.insert(vv.begin() + i + 1, e);

			//if(e == -1) break;
			if(e == -1) continue;

			if(e2s.find(e) == e2s.end())
			{
				set<int> ss;
				ss.insert(k);
				e2s.insert(PISI(e, ss));
			}
			else
			{
				e2s[e].insert(k);
			}

			//printf("line %d: insert %d between (%d, %d) = (%d, %d, %d)\n", k, e, x, y, vv[i], vv[i + 1], vv[i + 2]);

			//break;
		}
	}
	return 0;
}

bool hyper_set::extend(int e)
{
	return (left_extend(e) || right_extend(e));
}

bool hyper_set::left_extend(int e)
{
	if(e2s.find(e) == e2s.end()) return false;
	set<int> s = e2s[e];
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		assert(vv.size() >= 1);

		for(int i = 1; i < vv.size(); i++)
		{
			if(vv[i] == e && vv[i - 1] != -1) return true; 
		}
	}
	return false;
}

bool hyper_set::right_extend(int e)
{
	if(e2s.find(e) == e2s.end()) return false;
	set<int> s = e2s[e];
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		assert(vv.size() >= 1);

		for(int i = 0; i < vv.size() - 1; i++)
		{
			if(vv[i] == e && vv[i + 1] != -1) return true; 
		}
	}
	return false;
}

bool hyper_set::left_extend(const vector<int> &s)
{
	for(int i = 0; i < s.size(); i++)
	{
		if(left_extend(s[i]) == true) return true;
	}
	return false;
}

bool hyper_set::right_extend(const vector<int> &s)
{
	for(int i = 0; i < s.size(); i++)
	{
		if(right_extend(s[i]) == true) return true;
	}
	return false;
}

bool hyper_set::left_dominate(int e)
{
	// for each appearance of e
	// if right is not empty then left is also not empty
	if(e2s.find(e) == e2s.end()) return true;

	set<PI> x1;
	set<PI> x2;
	set<int> s = e2s[e];
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		assert(vv.size() >= 1);

		for(int i = 0; i < vv.size() - 1; i++)
		{
			if(vv[i] != e) continue;
			if(vv[i + 1] == -1) continue;

			if(i == 0 || vv[i - 1] == -1)
			{
				if(i + 2 < vv.size()) x1.insert(PI(vv[i + 1], vv[i + 2]));
				else x1.insert(PI(vv[i + 1], -1));
			}
			else
			{
				x2.insert(PI(vv[i + 1], -1));
				if(i + 2 < vv.size()) x2.insert(PI(vv[i + 1], vv[i + 2]));
			}
		}
	}

	for(set<PI>::iterator it = x1.begin(); it != x1.end(); it++)
	{
		PI p = (*it);
		if(x2.find(p) == x2.end()) return false;
	}
	return true;
}

bool hyper_set::right_dominate(int e)
{
	// for each appearance of e
	// if left is not empty then right is also not empty
	if(e2s.find(e) == e2s.end()) return true;
	set<PI> x1;
	set<PI> x2;
	set<int> s = e2s[e];
	for(set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		assert(vv.size() >= 1);
		for(int i = 1; i < vv.size(); i++)
		{
			if(vv[i] != e) continue;
			if(vv[i - 1] == -1) continue;

			if(i == vv.size() - 1 || vv[i + 1] == -1)
			{
				if(i - 2 >= 0) x1.insert(PI(vv[i - 1], vv[i - 2]));
				else x1.insert(PI(vv[i - 1], -1));
			}
			else
			{
				x2.insert(PI(vv[i - 1], -1));
				if(i - 2 >= 0) x2.insert(PI(vv[i - 1], vv[i - 2]));
			}
		}
	}

	for(set<PI>::iterator it = x1.begin(); it != x1.end(); it++)
	{
		PI p = (*it);
		if(x2.find(p) == x2.end()) return false;
	}

	return true;
}

int hyper_set::print_nodes()
{
	//printf("PRINT HYPER_SET\n");
	for(MVII::iterator it = nodes.begin(); it != nodes.end(); it++)
	{
		const vector<int> &v = it->first;
		int c = it->second;
		printf("hyper-edge (nodes), counts = %d, list = ( ", c); 
		printv(v);
		printf(")\n");
	}
	return 0;
}

int hyper_set::print_edges()
{
	for(int i = 0; i < edges.size(); i++)
	{
		printf("hyper-edge (edges) %d, counts = %d, edge = ( ", i, ecnts[i]);
		printv(edges[i]);
		printf(")\n");
	}
	return 0;
}

int hyper_set::write(ostream &os) const
{
	os<<fixed;
	os.precision(2);
	
	for(MVII::const_iterator it = nodes.begin(); it != nodes.end(); it++)
	{
		const vector<int> &v = it->first;
		int n = v.size();
		if(n <= 2) continue;

		int c = it->second;
		os << "path " << n;
		for(int k = 0; k < n; k++) os << " " << v[k];
		os << " " << c << " 1";
		os << endl;
	}
	//os << endl;
	return 0;
}
