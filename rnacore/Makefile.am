noinst_LIBRARIES = librnacore.a

UTIL_INCLUDE = $(top_srcdir)/util
GRAPH_INCLUDE = $(top_srcdir)/graph
GTF_INCLUDE = $(top_srcdir)/gtf

librnacore_a_CPPFLAGS = -std=c++11 -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE) -I$(GTF_INCLUDE)

librnacore_a_SOURCES = splice_graph.h splice_graph.cc \
					   vertex_info.h vertex_info.cc \
					   edge_info.h edge_info.cc \
					   interval_map.h interval_map.cc \
					   binomial.h binomial.cc \
					   hit.h hit.cc \
					   hit_core.h hit_core.cc \
					   partial_exon.h partial_exon.cc \
					   region.h region.cc \
					   junction.h junction.cc \
					   essential.h essential.cc \
					   chain_set.h chain_set.cc \
					   pereads_cluster.h pereads_cluster.cc \
					   phase_set.h phase_set.cc \
					   transcript_set.h transcript_set.cc \
					   filter.h filter.cc \
					   sample_profile.h sample_profile.cc \
					   bundle_base.h bundle_base.cc \
					   disjoint_set.h disjoint_set.cc \
					   graph_builder.h graph_builder.cc \
					   graph_cluster.h graph_cluster.cc \
					   graph_reviser.h graph_reviser.cc
