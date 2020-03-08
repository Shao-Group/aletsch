#ifndef __DECOMPOSER_H__
#define __DECOMPOSER_H__

#include "merged_graph.h"
#include "hyper_set.h"
#include "transcript.h"
#include <iostream>
#include <mutex>

using namespace std;

int assemble();
int assemble(merged_graph cm, vector<merged_graph> children, map< size_t, vector<transcript> > &trsts, mutex &mylock);
int index_transcript(map< size_t, vector<transcript> > &mt, const transcript &t);
bool query_transcript(const map< size_t, vector<transcript> > &mt, const transcript &t);

#endif
