#include "combined_graph.h"
#include "essential.h"

combined_graph::combined_graph()
{
	num_combined = 0;
}

combined_graph::combined_graph(const string &line)
{
	hline = line;
	num_combined = 0;
}

int combined_graph::clear()
{
	imap.clear();
	junctions.clear();
	sbounds.clear();
	tbounds.clear();
	splices.clear();
	phase.clear();
	reads.clear();
	num_combined = 0;
	chrm = "";
	strand = '.';
	hline = "";
	for(int i = 0; i < children.size(); i++) children[i].clear();
	return 0;
}

PI32 combined_graph::get_bounds()
{
	if(sbounds.size() == 0 || tbounds.size() == 0) return PI32(-1, -1);
	map<int32_t, DI>::iterator x1 = sbounds.begin();
	map<int32_t, DI>::iterator x2 = tbounds.end();
	x2--;
	return PI32(x1->first, x2->first);
}

int combined_graph::print(int index)
{
	PI32 p = get_bounds();
	printf("combined-graph %d: #combined = %d, children = %lu, chrm = %s, strand = %c, #regions = %lu, #sbounds = %lu, #tbounds = %lu, #junctions = %lu, #phase = %lu, #splices = %lu, boundary = [%d, %d)\n", 
			index, num_combined, children.size(), chrm.c_str(), strand, std::distance(imap.begin(), imap.end()), sbounds.size(), tbounds.size(), junctions.size(), phase.size(), splices.size(), p.first, p.second);
	return 0;
}

int combined_graph::analyze(int index)
{
	return 0;
}
