/*
Part of Scallop Transcript Assembler
(c) 2017 by  Mingfu Shao, Carl Kingsford, and Carnegie Mellon University.
See LICENSE for licensing.
*/

#ifndef __HIT_H__
#define __HIT_H__

#include <string>
#include <vector>

#include "htslib/sam.h"

using namespace std;

/*! @typedef
 @abstract Structure for core alignment information.
 @field  tid     chromosome ID, defined by bam_hdr_t
 @field  pos     0-based leftmost coordinate
 @field  bin     bin calculated by bam_reg2bin()
 @field  qual    mapping quality
 @field  l_qname length of the query name
 @field  flag    bitwise flag
 @field  l_extranul length of extra NULs between qname & cigar (for alignment)
 @field  n_cigar number of CIGAR operations
 @field  l_qseq  length of the query sequence (read)
 @field  mtid    chromosome ID of next read in template, defined by bam_hdr_t
 @field  mpos    0-based leftmost coordinate of next read in template

typedef struct {
    int32_t tid;
    int32_t pos;
    uint16_t bin;
    uint8_t qual;
    uint8_t l_qname;
    uint16_t flag;
    uint8_t unused1;
    uint8_t l_extranul;
    uint32_t n_cigar;
    int32_t l_qseq;
    int32_t mtid;
    int32_t mpos;
    int32_t isize;
} bam1_core_t;
*/

/*! @typedef
 @abstract Structure for one alignment.
 @field  core       core information about the alignment
 @field  l_data     current length of bam1_t::data
 @field  m_data     maximum length of bam1_t::data
 @field  data       all variable-length data, concatenated; structure: qname-cigar-seq-qual-aux

 @discussion Notes:

 1. qname is terminated by one to four NULs, so that the following
 cigar data is 32-bit aligned; core.l_qname includes these trailing NULs,
 while core.l_extranul counts the excess NULs (so 0 <= l_extranul <= 3).
 2. l_qseq is calculated from the total length of an alignment block
 on reading or from CIGAR.
 3. cigar data is encoded 4 bytes per CIGAR operation.
 4. seq is nybble-encoded according to bam_nt16_table.
 */

class hit: public bam1_core_t
{
public:
	hit(bam1_t *b, int id);
	hit(const hit &h);
	virtual ~hit();
	virtual bool operator<(const hit &h) const;
	virtual hit& operator=(const hit &h);

public:
	int hid;								// unique id for this hit, < 0 means removed
	int32_t rpos;							// right position mapped to reference [pos, rpos)
	int32_t nh;								// NH aux in sam
	int32_t hi;								// HI aux in sam
	int32_t nm;								// NM aux in sam
	char strand;							// strandness
	char xs;								// XS aux in sam
	char ts;								// ts tag used in minimap2
	string qname;							// query names

public:
	int set_tags(bam1_t *b);
	int set_strand(int lib_type);
	int print() const;
	size_t get_qhash() const;
	bool get_concordance() const;
	vector<int32_t> extract_splices(bam1_t *b) const;
	//int get_aligned_intervals(vector<int64_t> &v) const;
	//size_t get_phash() const;
	//bool equal(const hit &h) const;
};

#endif
