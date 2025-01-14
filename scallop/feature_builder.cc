#include "feature_builder.h"
#include <cfloat>

feature_builder::feature_builder(const parameters & p)
	: cfg(p)
{}

int feature_builder::build_input_gtf(splice_graph &gr, const vector<transcript> &trsts, const map<int64_t, vector<int>> & tmap)
{
	return 0;
}

int feature_builder::build_transcripts(splice_graph &gr, vector<path> &paths, vector<transcript> &trsts)
{
	trsts.clear();
    //for(int i = 0; i < paths.size(); i++) paths[i].print(i);
    //
    string prefix = "v"+to_string(cfg.min_num_exons)+"-"+to_string(cfg.max_num_exons);
    string path_file = prefix + ".path.csv";
    ofstream outputPath(path_file, ios::app);
    if(!outputPath.is_open()) 
    {
        printf("open file %s error.\n", path_file.c_str());
        return 0;
    }

	for(int i = 0; i < paths.size(); i++)
	{
		string tid = "chr" + gr.chrm + "." + gr.gid + "." + tostring(i);
		transcript trst;
		path &p = paths[i];
        update_trst_features(gr, trst, i, paths);
		build_transcript(gr, trst, p, tid);
		trsts.push_back(trst);

        if(p.junc.size() == 0) continue;
        string chr_gid = "chr" + gr.chrm + "." + gr.gid;
        outputPath << gr.chrm << "," << chr_gid << "," << tid << ",\"";
        for(int j = 1; j < p.v.size()-1; j++)
        {
            outputPath << p.v[j]-1;
            if(j < p.v.size()-2) outputPath << ",";

			// Update how final paths support raw vertices and edges
			int v1 = p.v[j-1];
			int v2 = p.v[j];
			assert(gr.edge(v1, v2).second);
			edge_descriptor e = gr.edge(v1, v2).first;
			edge_info &ei = gr.get_editable_edge_info(e);
			vertex_info &vi2 = gr.get_editable_vertex_info(v2);
			ei.trstSupport = 1;
			vi2.trstSupport = 1;
        }
        outputPath << "\",\"";

		//output splice junction's source vertex
		for(int j = 0; j < p.junc.size(); j++)
        {
            outputPath << p.junc[j].first-1;
            if(j < p.junc.size()-1) outputPath << ",";
        }
		outputPath << "\",\"";

		//output splice junction's target vertex
		for(int j = 0; j < p.junc.size(); j++)
        {
            outputPath << p.junc[j].second-1;
            if(j < p.junc.size()-1) outputPath << ",";
        }
		outputPath << "\",";
		outputPath << p.weight << "\n";
	}

    outputPath.close();
	return 0;
}

int feature_builder::outputPhasingPath(splice_graph &gr, hyper_set &hs)
{
    string prefix = "v"+to_string(cfg.min_num_exons)+"-"+to_string(cfg.max_num_exons);
    string path_file = prefix + ".phasing.csv";
    ofstream outputPath(path_file, ios::app);
    if(!outputPath.is_open()) 
    {
        printf("open file %s error.\n", path_file.c_str());
        return 0;
    }

	// Downsample
    MVII selectedPaths = hs.nodes;
    if (gr.num_vertices() > 100 && hs.nodes.size() > gr.num_vertices()) {
        selectedPaths = downsamplePhasingPaths(hs.nodes, gr.num_vertices());
		printf("Downsample phasing path %ld/%ld", selectedPaths.size(), hs.nodes.size());
    }

	int i = 0;
	for(MVII::iterator it = selectedPaths.begin(); it != selectedPaths.end(); it++)
	{
		const vector<int> &v = it->first;
		if(v.size() <= 2) continue;

		int c = it->second;
        string pid = "chr" + gr.chrm + "." + gr.gid + ".phasing." + tostring(i);
        string chr_gid = "chr" + gr.chrm + "." + gr.gid;

        outputPath << gr.chrm << "," << chr_gid << "," << pid << ",\"";
        for(int j = 0; j < v.size(); j++)
        {
            outputPath << v[j]-1;
            if(j < v.size()-1) outputPath << ",";
        }
        outputPath << "\"," << c << "\n";
		i++;
	}
	outputPath.close();
	return 0;
}

MVII feature_builder::downsamplePhasingPaths(const MVII& paths, size_t targetSize) {

    int maxCount = 0;
    size_t maxLength = 0;
    for (const auto& p : paths) {
        maxCount = max(maxCount, p.second);
        maxLength = max(maxLength, p.first.size());
    }

    vector<pair<double, pair<vector<int>, int>>> scoredPaths;
    scoredPaths.reserve(paths.size());

    // Calculate scores by lengths and counts
    for (const auto& p : paths) {
        if (p.first.size() <= 2) continue;
        double score = static_cast<double>(p.first.size()) / maxLength + static_cast<double>(p.second) / maxCount;
        scoredPaths.push_back({score, {p.first, p.second}});
    }

    // Sort by score
    sort(scoredPaths.begin(), scoredPaths.end(),
         [](const pair<double, pair<vector<int>, int>>& a, const pair<double, pair<vector<int>, int>>& b) { return a.first > b.first; });

    MVII result;
    for (size_t i = 0; i < targetSize && i < scoredPaths.size(); i++) {
        result[scoredPaths[i].second.first] = scoredPaths[i].second.second;
    }

    return result;
}

int feature_builder::update_trst_features(splice_graph &gr, transcript &trst, int pid, vector<path> &paths)
{
    path &p = paths[pid];

    int n = p.v.size();
    //printf("Path size:%d\n", n);
    assert(n >= 3);
    trst.features.num_vertices = n-2;
    trst.features.num_edges = n-3;
    trst.features.gr_vertices = gr.num_vertices();
    trst.features.gr_edges = gr.num_edges();
    trst.features.gr_reads = gr.reads;
    trst.features.gr_subgraph = gr.subgraph;
    trst.features.max_mid_exon_len = 0;

    int junc = p.junc.size();
    if(junc == 0) return 0;//ignore single exon

    int start_splicing_v = p.junc.front().first;
    int end_splicing_v = p.junc.back().second;
    auto it_s = lower_bound(p.v.begin(), p.v.end(), start_splicing_v);
    auto it_t = lower_bound(p.v.begin(), p.v.end(), end_splicing_v);
    if(it_s == p.v.end() || *it_s != start_splicing_v || it_t == p.v.end() || *it_t != end_splicing_v) assert(false);
    trst.features.junc_ratio = 1.0*junc/(it_t-it_s);

    //printf("path%d, #junc=%d\n", pid, junc);
    for(int i = 1; i < junc; i++)
    {
        int exon_len = gr.get_vertex_info(p.junc[i].first).rpos-gr.get_vertex_info(p.junc[i-1].second).lpos;
        trst.features.max_mid_exon_len = max(trst.features.max_mid_exon_len, exon_len);
    }
    
    const vertex_info & svi =  gr.get_vertex_info(p.v[1]);
    const vertex_info & evi = gr.get_vertex_info(p.v[n-2]);
    trst.features.start_loss1 = svi.boundary_loss1;
    trst.features.start_loss2 = svi.boundary_loss2;
    trst.features.start_loss3 = svi.boundary_loss3;
    trst.features.end_loss1 = evi.boundary_loss1;
    trst.features.end_loss2 = evi.boundary_loss2;
    trst.features.end_loss3 = evi.boundary_loss3;
    trst.features.start_merged_loss = svi.boundary_merged_loss;
    trst.features.end_merged_loss = evi.boundary_merged_loss;

    trst.features.uni_junc = unique_junc(paths, pid);
    trst.features.introns = 0;
    trst.features.start_introns = 0;
    trst.features.end_introns = 0;
    trst.features.intron_ratio = 0.0;
    trst.features.start_intron_ratio = 0.0;
    trst.features.end_intron_ratio = 0.0;
    if(junc == 0) return 0;
    for(int i = 0; i < paths.size(); i++)
    {
        if(i == pid) continue;

        vector<pair<int, int>>& junc1 = p.junc, junc2 = paths[i].junc;
        if(junc1.size()<2 || junc2.size()<1) continue;
        int intron_cnt = 0;
        int start_intron = 0;
        int end_intron = 0;
        for (size_t i = 0; i < junc1.size(); ++i) 
        {
            for (size_t j = 0; j < junc2.size(); ++j) 
            {
                if(i == 0)
                {
                    if(junc2[j].first >= p.v[1] && junc2[j].second <= junc1[0].first)
                    {
                        start_intron++;
                        int v1 = junc2[j].first;
                        int v2 = junc2[j].second;
                        assert(gr.edge(v1, v2).second);
                        edge_descriptor e = gr.edge(v1, v2).first;
                        assert(gr.edge(v1, v1+1).second);
                        assert(gr.edge(v2-1, v2).second);
                        edge_descriptor e1 = gr.edge(v1, v1+1).first;
                        edge_descriptor e2 = gr.edge(v2-1, v2).first;
                        trst.features.start_intron_ratio = max(trst.features.start_intron_ratio,gr.get_edge_weight(e)/min(gr.get_edge_weight(e1), gr.get_edge_weight(e2)));
                    }
                }
                else if(junc2[j].second <= junc1[i].first && junc2[j].first >= junc1[i-1].second)
                {
                    //intron_cnt += junc2[j].second - junc2[j].first - 1;
                    intron_cnt++;

                    int v1 = junc2[j].first;
                    int v2 = junc2[j].second;
                    assert(gr.edge(v1, v2).second);
                    edge_descriptor e = gr.edge(v1, v2).first;
                    assert(gr.edge(v1, v1+1).second);
                    assert(gr.edge(v2-1, v2).second);
                    edge_descriptor e1 = gr.edge(v1, v1+1).first;
                    edge_descriptor e2 = gr.edge(v2-1, v2).first;
                    trst.features.intron_ratio = max(trst.features.intron_ratio, gr.get_edge_weight(e)/min(gr.get_edge_weight(e1), gr.get_edge_weight(e2)));
                }

                if(i == junc1.size()-1)
                {
                    if(junc2[j].first >= junc1[i].second && junc2[j].second <= p.v[n-2])
                    {
                        end_intron++;
                        int v1 = junc2[j].first;
                        int v2 = junc2[j].second;
                        assert(gr.edge(v1, v2).second);
                        edge_descriptor e = gr.edge(v1, v2).first;
                        assert(gr.edge(v1, v1+1).second);
                        assert(gr.edge(v2-1, v2).second);
                        edge_descriptor e1 = gr.edge(v1, v1+1).first;
                        edge_descriptor e2 = gr.edge(v2-1, v2).first;
                        trst.features.end_intron_ratio = max(trst.features.end_intron_ratio,gr.get_edge_weight(e)/min(gr.get_edge_weight(e1), gr.get_edge_weight(e2)));
                    }
                }

            }
        }
        trst.features.introns = max(trst.features.introns, intron_cnt);
        trst.features.start_introns = max(trst.features.start_introns, start_intron);
        trst.features.end_introns = max(trst.features.end_introns, end_intron);


    }

    trst.features.seq_min_wt = DBL_MAX;
    trst.features.seq_min_cnt = INT_MAX;
    trst.features.seq_min_abd = DBL_MAX;
    trst.features.seq_min_ratio = 1.0;
    trst.features.seq_max_wt = 0;
    trst.features.seq_max_cnt = 0;
    trst.features.seq_max_abd = 0;
    trst.features.seq_max_ratio = 0;


    trst.features.unbridge_start_coming_count = 0;
    trst.features.unbridge_start_coming_ratio = 0;
    trst.features.unbridge_end_leaving_count = 0;
    trst.features.unbridge_end_leaving_ratio = 0;
    
    //gr.print();
    for(int i = 1; i < n; i++)
    {
        int v1 = p.v[i-1];
        int v2 = p.v[i];
        assert(gr.edge(v1, v2).second);
        edge_descriptor e = gr.edge(v1, v2).first;
        edge_info ei = gr.get_edge_info(e);
        vertex_info vi1 = gr.get_vertex_info(v1);
        vertex_info vi2 = gr.get_vertex_info(v2);

        //if(v1 >= p.junc[0].first && v2 <= p.junc[junc-1].second)
        
        trst.features.seq_min_wt = min(trst.features.seq_min_wt, gr.get_edge_weight(e));
        trst.features.seq_min_cnt = min(trst.features.seq_min_cnt, ei.count);
        trst.features.seq_min_abd = min(trst.features.seq_min_abd, ei.abd);
        trst.features.seq_min_ratio = min(trst.features.seq_min_ratio, gr.get_edge_weight(e)/max(gr.get_in_weights(v2), gr.get_out_weights(v1)));

        trst.features.seq_max_wt = max(trst.features.seq_max_wt, gr.get_edge_weight(e));
        trst.features.seq_max_cnt = max(trst.features.seq_max_cnt, ei.count);
        trst.features.seq_max_abd = max(trst.features.seq_max_abd, ei.abd);
        trst.features.seq_max_ratio = max(trst.features.seq_max_ratio, gr.get_edge_weight(e)/max(gr.get_in_weights(v2), gr.get_out_weights(v1)));

        
        if(i == 1)
        {
            trst.features.unbridge_start_coming_count = vi2.unbridge_coming_count;
            trst.features.unbridge_start_coming_ratio = vi2.unbridge_coming_ratio;
            trst.features.start_cnt = ei.count;
            trst.features.start_weight = gr.get_edge_weight(e);
            trst.features.start_abd = ei.abd;
        }
        else if(i == n-2)
        {
            trst.features.unbridge_end_leaving_count = vi2.unbridge_leaving_count;
            trst.features.unbridge_end_leaving_ratio = vi2.unbridge_leaving_ratio;
        }
        else if(i == n-1)
        {
            trst.features.end_cnt = ei.count;
            trst.features.end_weight = gr.get_edge_weight(e);
            trst.features.end_abd = ei.abd;
        }
    }

    return 0;
}

int feature_builder::infer_introns(const vector<pair<int, int>>& junc1, const vector<pair<int, int>>& junc2) 
{
    if(junc1.size()<2 || junc2.size()<1) return 0;

    int sum = 0;
    for (size_t i = 1; i < junc1.size(); ++i) 
    {
        for (size_t j = 0; j < junc2.size(); ++j) 
        {
            if(junc2[j].second <= junc1[i].first && junc2[j].first >= junc1[i-1].second)
            {
                sum += junc2[j].second - junc2[j].first - 1;
            }

        }
    }
    return sum;
}

int feature_builder::unique_junc(const vector<path>& paths, int i) 
{
    map<pair<int, int>, int> juncUni;

    for (size_t idx = 0; idx < paths.size(); ++idx) {
        for (const auto& pair : paths[idx].junc) {
            // If the pair is not in the map, add it
            if (juncUni.find(pair) == juncUni.end()) {
                juncUni[pair] = idx;
            } else if (juncUni[pair] != idx && juncUni[pair] != -1) {
                // If the pair is already in the map but from a different path, mark it as -1
                juncUni[pair] = -1;
            }
        }
    }

    // Count the unique pairs from paths[i]
    int uniqueCount = 0;
    for (const auto& pair1 : paths[i].junc) {
        if (juncUni.find(pair1) != juncUni.end() && juncUni[pair1] == i) {
            uniqueCount++;
        }
    }

    return uniqueCount;
}

