#ifndef __COMBINER_H__
#define __COMBINER_H__

#include "interval_map.h"

typedef pair<int32_t, int32_t> PI32;
typedef pair<double, int> DI;

class combined_graph
{
public:
	combined_graph();
	combined_graph(const string &line);

public:
	split_interval_map imap;
	map<PI32, DI> junctions;
	map<int32_t, DI> sbounds;
	map<int32_t, DI> tbounds;
	vector<int32_t> splices;
	map<vector<int32_t>, DI> phase;
	map<vector<int32_t>, DI> paths;
	int num_combined;
	string chrm;
	char strand;
	string hline;

	vector<combined_graph> children;

public:
	int combine(const combined_graph &gt);
	int combine_regions(const combined_graph &gt);
	int combine_junctions(const combined_graph &gt);
	int combine_start_bounds(const combined_graph &gt);
	int combine_end_bounds(const combined_graph &gt);
	int combine_paths(const combined_graph &gt);
	int combine_phase(const combined_graph &gt);
	int combine_splice_positions(const combined_graph &gt);

	PI32 get_bounds();
	int get_overlapped_splice_positions(const vector<int32_t> &v) const;

	int build(istream &is, const string &chrm, char c);
	int write(ostream &os, int index, bool headers = false);
	int write(ostream &os);

	int print(int index);
	int analyze(int index);
};

#endif
