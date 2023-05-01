/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __EDGE_INFO__
#define __EDGE_INFO__

#include <set>
#include <string>
using namespace std;

class edge_info
{
public:
	edge_info();
	edge_info(int l);
	edge_info(const edge_info &ei);

public:
	double stddev;
	int length;
	int type;
	int jid;		// junction id
	int count;		// #supporting subgraphs
    int cntsam;     // #supporting sampels
	double weight;	// new weight from hyper-edges
	int strand;		// strandness ./+/- => 0,1,2

    set<string> supports;
    set<int> samples;
};

#endif
