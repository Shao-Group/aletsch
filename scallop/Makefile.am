noinst_LIBRARIES = libscallop.a

GTF_INCLUDE = $(top_srcdir)/lib/gtf
UTIL_INCLUDE = $(top_srcdir)/lib/util
GRAPH_INCLUDE = $(top_srcdir)/lib/graph

#scallop_CPPFLAGS = -std=c++11 -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE)
libscallop_a_CPPFLAGS = -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE)

libscallop_a_SOURCES = splice_graph.h splice_graph.cc \
					   super_graph.h super_graph.cc \
					   vertex_info.h vertex_info.cc \
					   edge_info.h edge_info.cc \
					   interval_map.h interval_map.cc \
					   config.h config.cc \
					   hit.h hit.cc \
					   partial_exon.h partial_exon.cc \
					   hyper_set.h hyper_set.cc \
					   subsetsum.h subsetsum.cc \
					   router.h router.cc \
					   region.h region.cc \
					   junction.h junction.cc \
					   bundle_base.h bundle_base.cc \
					   bundle.h bundle.cc \
					   path.h path.cc \
					   equation.h equation.cc \
					   gtf.h gtf.cc \
					   scallop.h scallop.cc \
					   previewer.h previewer.cc \
					   filter.h filter.cc \
					   cluster.h cluster.cc \
					   graph_revise.h graph_revise.cc \
					   hyper_graph.h hyper_graph.cc
