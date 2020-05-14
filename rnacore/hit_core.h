#ifndef __HIT_CORE_H__
#define __HIT_CORE_H__

#include <string>
#include <vector>

#include "hit.h"

using namespace std;


/*! @typedef
 @abstract Structure for core alignment information.
 @field  tid     chromosome ID, defined by bam_hdr_t
 @field  pos     0-based leftmost coordinate
 @field  bin     bin calculated by bam_reg2bin()
 @field  qual    mapping quality
 @field  l_qname length of the query name
 @field  flag    bitwise flag
 @field  n_cigar number of CIGAR operations
 @field  l_qseq  length of the query sequence (read)
 @field  mtid    chromosome ID of next read in template, defined by bam_hdr_t
 @field  mpos    0-based leftmost coordinate of next read in template

typedef struct {
    int32_t tid;
    int32_t pos;
    uint32_t bin:16, qual:8, l_qname:8;
    uint32_t flag:16, n_cigar:16;
    int32_t l_qseq;
    int32_t mtid;
    int32_t mpos;
    int32_t isize;
} bam1_core_t;
*/

class hit_core: public bam1_core_t
{
public:
	hit_core(const hit &h);
	hit_core(const hit_core &h);

public:
	int32_t rpos;							// right position mapped to reference [pos, rpos)
	int32_t nh;								// NH aux in sam
	int32_t hi;								// HI aux in sam
	char xs;								// XS aux in sam
	string qname;							// query name
};

#endif
