bin_PROGRAMS = meta-scallop

GTF_INCLUDE = $(top_srcdir)/gtf
UTIL_INCLUDE = $(top_srcdir)/util
GRAPH_INCLUDE = $(top_srcdir)/graph
RNACORE_INCLUDE = $(top_srcdir)/rnacore
SCALLOP_INCLUDE = $(top_srcdir)/scallop
BRIDGE_INCLUDE = $(top_srcdir)/bridge

GTF_LIB = $(top_builddir)/gtf
UTIL_LIB = $(top_builddir)/util
GRAPH_LIB = $(top_builddir)/graph
RNACORE_LIB = $(top_builddir)/rnacore
SCALLOP_LIB = $(top_builddir)/scallop
BRIDGE_LIB = $(top_builddir)/bridge

meta_scallop_CPPFLAGS = -std=c++11 -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE) -I$(SCALLOP_INCLUDE) -I$(RNACORE_INCLUDE) -I$(BRIDGE_INCLUDE)
meta_scallop_LDFLAGS = -pthread -L$(GTF_LIB) -L$(GRAPH_LIB) -L$(UTIL_LIB) -L$(SCALLOP_LIB) -L$(RNACORE_LIB) -L$(BRIDGE_LIB)
meta_scallop_LDADD = -lscallop -lgtf -lgraph -lutil -lrnacore -lbridge
meta_scallop_SOURCES = combined_graph.h combined_graph.cc \
					   combined_group.h combined_group.cc \
					   generator.h generator.cc \
					   previewer.h previewer.cc \
					   sample_profile.h sample_profile.cc \
					   incubator.h incubator.cc \
					   main.cc
