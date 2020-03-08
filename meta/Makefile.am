bin_PROGRAMS = meta-scallop

GTF_INCLUDE = $(top_srcdir)/lib/gtf
UTIL_INCLUDE = $(top_srcdir)/lib/util
GRAPH_INCLUDE = $(top_srcdir)/lib/graph
SCALLOP_INCLUDE = $(top_srcdir)/scallop
MERGE_INCLUDE = $(top_srcdir)/merge

GTF_LIB = $(top_builddir)/lib/gtf
UTIL_LIB = $(top_builddir)/lib/util
GRAPH_LIB = $(top_builddir)/lib/graph
SCALLOP_LIB = $(top_builddir)/scallop
MERGE_LIB = $(top_builddir)/merge

#scallop_CPPFLAGS = -std=c++11 -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE)
meta_scallop_CPPFLAGS = -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE) -I$(SCALLOP_INCLUDE) -I$(MERGE_INCLUDE)
meta_scallop_LDFLAGS = -pthread -L$(GTF_LIB) -L$(GRAPH_LIB) -L$(UTIL_LIB) -L$(SCALLOP_LIB) -L$(MERGE_LIB)
meta_scallop_LDADD = -lmerge -lscallop -lgtf -lgraph -lutil
meta_scallop_SOURCES = decomposer.h decomposer.cc \
					   generator.h generator.cc \
					   merged_graph.h merged_graph.cc \
					   merger.h merger.cc \
					   config.h config.cc \
					   main.cc
