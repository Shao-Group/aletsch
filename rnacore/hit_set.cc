#include "hit_set.h"
#include "util.h"
#include "essential.h"

#include <algorithm>

hit_set::hit_set(vector<hit> &h)
	: hits(h)
{} 

int hit_set::build_paired_reads_clusters()
{
	typedef pair< vector<int32_t>, vector<int> > PVV;
	map<PVV, int> findex;			// index for clusters

	vprc.clear();
	paired.assign(hits.size(), false);

	vector<PI> fs;
	build_paired_reads(hits, fs);

	for(int i = 0; i < fs.size(); i++)
	{
		int h1 = fs[i].first;
		int h2 = fs[i].second;
		vector<int32_t> v1;
		vector<int32_t> v2;
		for(int k = 0; k < hits[h1].spos.size(); k++)
		{
			int32_t p1 = high32(hits[h1].spos[k]);
			int32_t p2 = low32(hits[h1].spos[k]);
			v1.push_back(p1);
			v1.push_back(p2);
		}
		for(int k = 0; k < hits[h2].spos.size(); k++)
		{
			int32_t p1 = high32(hits[h2].spos[k]);
			int32_t p2 = low32(hits[h2].spos[k]);
			v2.push_back(p1);
			v2.push_back(p2);
		}

		paired[h1] = true;
		paired[h2] = true;

		int32_t k1l = hits[h1].pos;
		int32_t k1r = hits[h1].rpos;
		int32_t k2l = hits[h2].pos;
		int32_t k2r = hits[h2].rpos;

		PVV pvv(v1, v2);
		if(findex.find(pvv) == findex.end())
		{
			rcluster r1;
			rcluster r2;
			r1.vv = v1;
			r2.vv = v2;
			r1.vl.push_back(k1l);
			r1.vr.push_back(k1r);
			r2.vl.push_back(k2l);
			r2.vr.push_back(k2r);

			findex.insert(pair<PVV, int>(pvv, vprc.size()));
			vprc.push_back(PRC(r1, r2));
		}
		else
		{
			int k = findex[pvv];
			assert(vprc[k].first.vv == v1);
			assert(vprc[k].second.vv == v2);
			vprc[k].first.vl.push_back(k1l);
			vprc[k].first.vr.push_back(k1r);
			vprc[k].second.vl.push_back(k2l);
			vprc[k].second.vr.push_back(k2r);
		}
	}

	return 0;
}
