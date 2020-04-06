noinst_LIBRARIES = librnacore.a

UTIL_INCLUDE = $(top_srcdir)/util
GRAPH_INCLUDE = $(top_srcdir)/graph
GTF_INCLUDE = $(top_srcdir)/gtf

librnacore_a_CPPFLAGS = -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE) -I$(GTF_INCLUDE)

librnacore_a_SOURCES = splice_graph.h splice_graph.cc \
					   vertex_info.h vertex_info.cc \
					   edge_info.h edge_info.cc \
					   interval_map.h interval_map.cc \
					   hit.h hit.cc \
					   partial_exon.h partial_exon.cc \
					   region.h region.cc \
					   junction.h junction.cc \
					   essential.h essential.cc \
					   pereads_cluster.h pereads_cluster.cc \
					   phase_set.h phase_set.cc \
					   bundle.h bundle.cc \
					   graph_builder.h graph_builder.cc \
					   graph_cluster.h graph_cluster.cc
