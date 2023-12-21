/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __EDGE_INFO__
#define __EDGE_INFO__

#include <set>
#include <unordered_map>
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
	int count;		// #supporting samples
	double weight;	// new weight from hyper-edges
	int strand;		// strandness ./+/- => 0,1,2
    double confidence; //log of reliability of every choice
	//string feature;	// feature

    set<int> samples;
    unordered_map<int, double> spAbd; //map sample_id to abd
    double abd;
};

#endif
