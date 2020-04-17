/*
Part of meta-scallop
(c) 2020 by Mingfu Shao, The Pennsylvania State University
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdint.h>
#include <map>
#include <set>
#include <sstream>
#include <cassert>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>

#include "constants.h"

using namespace std;

// macros: using int64_t for two int32_t
#define pack(x, y) (int64_t)((((int64_t)(x)) << 32) | ((int64_t)(y)))
#define high32(x) (int32_t)((x) >> 32)
#define low32(x) (int32_t)(((x) << 32) >> 32)

// definitions
typedef map<int32_t, int32_t> MI32;
typedef pair<int32_t, int32_t> PI32;
typedef map<int32_t, int> MPI;
typedef pair<int32_t, int> PPI;
typedef pair<int, int> PI;
typedef map<int, int> MI;

// common small functions
template<typename T>
string tostring(T t)
{
	ostringstream s;
	s << t;
	return s.str();
}

template<typename T>
T compute_overlap(const pair<T, T> &x, const pair<T, T> &y)
{
	assert(x.first <= x.second);
	assert(y.first <= y.second);
	if(x.first > y.first) return compute_overlap(y, x);
	assert(x.first <= y.first);
	if(y.first >= x.second) return x.second - y.first;
	if(x.second <= y.second) return x.second - y.first;
	else return y.second - y.first;
}

template<typename T>
int reverse(vector<T> &x)
{
	if(x.size() == 0) return 0;
	int i = 0;
	int j = x.size() - 1;
	while(i < j)
	{
		T t = x[i];
		x[i] = x[j];
		x[j] = t;
		i++;
		j--;
	}
	return 0;
}

template<typename T>
int max_element(const vector<T> &x)
{
	if(x.size() == 0) return -1;
	int k = 0;
	for(int i = 1; i < x.size(); i++)
	{
		if(x[i] <= x[k]) continue;
		k = i;
	}
	return k;
}

template<typename T>
int min_element(const vector<T> &x)
{
	if(x.size() == 0) return -1;
	int k = 0;
	for(int i = 1; i < x.size(); i++)
	{
		if(x[i] >= x[k]) continue;
		k = i;
	}
	return k;
}

template<typename T>
int printv(const vector<T> &x)
{
	for(int i = 0; i < x.size(); i++)
	{
		cout<< x[i] <<" ";
	}
	return 0;
}

template<typename T>
int compute_mean_dev(const vector<T> &v, int si, int ti, double &ave, double &dev)
{
	ave = -1;
	dev = -1;
	if(si >= ti) return 0;

	assert(si >= 0 && si < v.size());
	assert(ti > 0 && ti <= v.size());

	T sum = 0;
	for(int i = si; i < ti; i++)
	{
		sum += v[i];
	}

	ave = sum * 1.0 / (ti - si);

	double var;
	for(int i = si ; i < ti; i++)
	{
		var += (v[i] - ave) * (v[i] - ave);
	}

	dev = sqrt(var / (ti - si));
	return 0;
}

template<typename T>
vector<int> consecutive_subset(const vector<T> &ref, const vector<T> &x)
{
	vector<int> v;
	if(x.size() == 0) return v;
	if(ref.size() == 0) return v;
	if(x.size() > ref.size()) return v;
	for(int i = 0; i <= ref.size() - x.size(); i++)
	{
		if(ref[i] != x[0]) continue;
		int k = i;
		bool b = true;
		for(int j = 0; j < x.size(); j++)
		{
			if(x[j] != ref[j + k]) b = false;
			if(b == false) break;
		}
		if(b == false) continue;
		v.push_back(k);
	}
	return v;
}

template<typename K, typename V>
vector<K> get_keys(const map<K, V> &m)
{
	vector<K> v;
    typedef typename std::map<K,V>::const_iterator MIT;
	for(MIT it = m.begin(); it != m.end(); it++)
	{
		v.push_back(it->first);
	}
	return v;
}

template<typename T>
bool check_increasing_sequence(const vector<T> &x)
{
	if(x.size() <= 1) return true;
	for(int k = 0; k < x.size() - 1; k++)
	{
		if(x[k] >= x[k + 1]) return false;
	}
	return true;
}

bool check_identical(const vector<int> &x, int x1, int x2, const vector<int> &y, int y1, int y2);

// compare phasing paths / intron chains
template<typename T>
int compare_two_sorted_sequences(const vector<T> &ref, const vector<T> &qry)
{
	assert(ref.size() >= 1);
	assert(qry.size() >= 1);
	if(ref.back() < qry.front()) return FALL_RIGHT;
	if(ref.front() > qry.back()) return FALL_LEFT;

	typedef typename vector<T>::const_iterator VCI;
	VCI r1 = lower_bound(ref.begin(), ref.end(), qry.front());
	VCI q1 = lower_bound(qry.begin(), qry.end(), ref.front());
	assert(r1 != ref.end());
	assert(q1 != qry.end());

	int kr1 = r1 - ref.begin();
	int kq1 = q1 - qry.begin();
	assert(kr1 == 0 || kq1 == 0);

	VCI q2 = lower_bound(qry.begin(), qry.end(), ref.back());
	VCI r2 = lower_bound(ref.begin(), ref.end(), qry.back());
	assert(r2 != ref.end() || q2 != qry.end());

	int kr2 = r2 - ref.begin();
	int kq2 = q2 - qry.begin();

	if(*q1 == ref.front() || *r1 == qry.front())
	{
		if(r2 != ref.end() && q2 != qry.end())
		{
			assert(ref.back() == qry.back());
			assert(kr2 == ref.size() - 1);
			assert(kq2 == qry.size() - 1);
			bool b = check_identical(ref, kr1, kr2, qry, kq1, kq2);
			if(b == false) return CONFLICTING;
			if(b == true && kr1 == 0 && kq1 == 0) return IDENTICAL;
			if(b == true && kr1 >= 1 && kq1 == 0) return CONTAINED;
			if(b == true && kr1 == 0 && kq1 >= 1) return CONTAINING;
			assert(false);
		}
		else if(r2 != ref.end() && q2 == qry.end())
		{
			bool b = check_identical(ref, kr1, kr2, qry, kq1, qry.size() - 1);
			if(b == false) return CONFLICTING;
			if(b == true && kq1 == 0) return CONTAINED;
			if(b == true && kq1 >= 1) return EXTEND_LEFT;
			assert(false);
		}
		else if(r2 == ref.end() && q2 != qry.end())
		{
			bool b = check_identical(ref, kr1, ref.size() - 1, qry, kq1, kq2);
			if(b == false) return CONFLICTING;
			if(b == true && kr1 == 0) return CONTAINING;
			if(b == true && kr1 >= 1) return EXTEND_RIGHT;
		}
	}
	else if(*r1 > qry.front() && r2 == r1 && *r2 > qry.back()) return NESTED;
	else if(*q1 > ref.front() && q2 == q1 && *q2 > ref.back()) return NESTING;
	return CONFLICTING;
}

template<typename T>
bool merge_two_sorted_sequences(const vector<T> &ref, const vector<T> &qry, vector<T> &merged)
{
	merged.clear();
	if(ref.size() == 0) merged = qry;
	if(ref.size() == 0) return true;
	if(qry.size() == 0) merged = ref;
	if(qry.size() == 0) return true;

	typedef typename vector<T>::const_iterator VCI;

	merged.clear();
	int t = compare_two_sorted_sequences<int32_t>(ref, qry);

	if(t == CONFLICTING) return false;
	if(t == NESTED) return false;
	if(t == NESTING) return false;

	if(t == IDENTICAL) merged = ref;
	if(t == CONTAINED) merged = ref;
	if(t == CONTAINING) merged = qry;

	if(t == FALL_RIGHT) 
	{
		assert(ref.back() < qry.front());
		merged = ref;
		merged.insert(merged.end(), qry.begin(), qry.end());
	}
	if(t == FALL_LEFT) 
	{
		assert(qry.back() < ref.front());
		merged = qry;
		merged.insert(merged.end(), ref.begin(), ref.end());
	}
	if(t == EXTEND_LEFT)
	{
		VCI q1 = lower_bound(qry.begin(), qry.end(), ref.front());
		assert(*q1 == ref.front());
		merged.insert(merged.end(), qry.begin(), q1);
		merged.insert(merged.end(), ref.begin(), ref.end());
	}
	if(t == EXTEND_RIGHT)
	{
		VCI q2 = lower_bound(qry.begin(), qry.end(), ref.back());
		assert(*q2 == ref.back());
		merged.insert(merged.end(), ref.begin(), ref.end());
		merged.insert(merged.end(), q2 + 1, qry.end());
	}
	return true;
}

template<typename T>
bool overlap_two_sorted_sequences(const vector<T> &ref, const vector<T> &qry, vector<T> &overlap)
{
	overlap.clear();
	if(ref.size() == 0 || qry.size() == 0) return true;

	typedef typename vector<T>::const_iterator VCI;

	overlap.clear();
	int t = compare_two_sorted_sequences<int32_t>(ref, qry);

	if(t == CONFLICTING) return false;
	if(t == NESTED) return false;
	if(t == NESTING) return false;
	if(t == FALL_RIGHT) return false;
	if(t == FALL_LEFT) return false;

	if(t == IDENTICAL) overlap = ref;
	if(t == CONTAINED) overlap = qry;
	if(t == CONTAINING) overlap = ref;

	if(t == EXTEND_LEFT)
	{
		VCI q1 = lower_bound(qry.begin(), qry.end(), ref.front());
		assert(*q1 == ref.front());
		overlap.insert(overlap.end(), q1, qry.end());
	}
	if(t == EXTEND_RIGHT)
	{
		VCI q2 = lower_bound(qry.begin(), qry.end(), ref.back());
		assert(*q2 == ref.back());
		overlap.insert(overlap.end(), qry.begin(), q2 + 1);
	}
	return true;
}

vector<int> get_random_permutation(int n);
size_t string_hash(const std::string& str);
size_t vector_hash(const vector<int32_t> &str);
vector<string> split_string(const string& str, const string& delim);
vector<int> project_vector(const vector<int> &v, const map<int, int> &m);
int transform_vertex_set_map(const set<int> &s, map<int, int> &m);

#endif
