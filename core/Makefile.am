noinst_LIBRARIES = libcore.a

GTF_INCLUDE = $(top_srcdir)/lib/gtf
UTIL_INCLUDE = $(top_srcdir)/lib/util
GRAPH_INCLUDE = $(top_srcdir)/lib/graph

#core_CPPFLAGS = -std=c++11 -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE)
libcore_a_CPPFLAGS = -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE)

libcore_a_SOURCES = splice_graph.h splice_graph.cc \
					vertex_info.h vertex_info.cc \
					edge_info.h edge_info.cc \
					interval_map.h interval_map.cc \
					hit.h hit.cc \
					partial_exon.h partial_exon.cc \
					hyper_set.h hyper_set.cc \
					region.h region.cc \
					junction.h junction.cc \
					bundle_base.h bundle_base.cc \
					path.h path.cc \
					constants.h constants.cc
