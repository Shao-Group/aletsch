noinst_LIBRARIES = libmeta.a

GTF_INCLUDE = $(top_srcdir)/gtf
UTIL_INCLUDE = $(top_srcdir)/util
GRAPH_INCLUDE = $(top_srcdir)/graph
RNACORE_INCLUDE = $(top_srcdir)/rnacore
SCALLOP_INCLUDE = $(top_srcdir)/scallop
BRIDGE_INCLUDE = $(top_srcdir)/bridge

libmeta_a_CPPFLAGS = -std=c++11 -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE) -I$(SCALLOP_INCLUDE) -I$(RNACORE_INCLUDE) -I$(BRIDGE_INCLUDE)
libmeta_a_SOURCES = combined_graph.h combined_graph.cc \
					graph_group.h graph_group.cc \
					bundle.h bundle.cc \
					bundle_group.h bundle_group.cc \
					generator.h generator.cc \
					assembler.h assembler.cc \
					previewer.h previewer.cc \
					incubator.h incubator.cc
