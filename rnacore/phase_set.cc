#include "essential.h"
#include "constants.h"
#include "phase_set.h"

int phase_set::add_node_list(const set<int> &s)
{
	return add_node_list(s, 1);
}

int phase_set::add_node_list(const set<int> &s, int c)
{
	vector<int> v(s.begin(), s.end());
	return add_node_list(v, c);
}

int phase_set::add_node_list(const vector<int> &s, int c)
{
	// make sure the given list is reference to a splice graph
	// i.e., 0 and n are start and end vertices
	vector<int> v = s;
	sort(v.begin(), v.end());
	//for(int i = 0; i < v.size(); i++) v[i]++;
	if(nodes.find(v) == nodes.end()) nodes.insert(PVII(v, c));
	else nodes[v] += c;
	return 0;
}

int phase_set::clear()
{
	nodes.clear();
	return 0;
}
