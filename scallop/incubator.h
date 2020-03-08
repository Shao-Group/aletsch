#ifndef __INCUBATOR_H__
#define __INCUBATOR_H__

#include "combined_graph.h"
#include "hyper_set.h"
#include "transcript.h"
#include <iostream>
#include <mutex>

using namespace std;

int assemble();
int assemble(combined_graph cm, vector<combined_graph> children, map< size_t, vector<transcript> > &trsts, mutex &mylock);
int index_transcript(map< size_t, vector<transcript> > &mt, const transcript &t);
bool query_transcript(const map< size_t, vector<transcript> > &mt, const transcript &t);
int test(int k);

#endif
