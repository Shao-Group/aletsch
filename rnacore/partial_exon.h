/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __PARTIAL_EXON_H__
#define __PARTIAL_EXON_H__

#include <stdint.h>
#include <vector>
#include <string>

using namespace std;

class partial_exon
{
public:
	partial_exon(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype);

public:
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	int ltype;						// type of the left boundary
	int rtype;						// type of the right boundary

	double ave;						// average abundance
	double dev;						// standard-deviation of abundance
	double max;						// largest coverage in this partial exon
	double pvalue;					// significance

    double indel_sum_cov; //sum of indels in the partial exon
    double indel_ratio; //ratio of indel coverags over partial exon coverage
    int left_indel; //nearest distance of indel to left boundary
    int right_indel; //nearest distance of indel to right boundary

public:
	string label() const;
	int print(int index) const;
};

#endif
