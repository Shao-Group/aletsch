noinst_LIBRARIES = libmerge.a

GTF_INCLUDE = $(top_srcdir)/lib/gtf
UTIL_INCLUDE = $(top_srcdir)/lib/util
GRAPH_INCLUDE = $(top_srcdir)/lib/graph
SCALLOP_INCLUDE = $(top_srcdir)/scallop

#merge_CPPFLAGS = -std=c++11 -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE)
libmerge_a_CPPFLAGS = -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE) -I$(SCALLOP_INCLUDE)

libmerge_a_SOURCES = combined_graph.h combined_graph.cc \
					 combined_group.h combined_group.cc 
