noinst_LIBRARIES = libscallop.a

GTF_INCLUDE = $(top_srcdir)/lib/gtf
UTIL_INCLUDE = $(top_srcdir)/lib/util
GRAPH_INCLUDE = $(top_srcdir)/lib/graph
RNACORE_INCLUDE = $(top_srcdir)/rnacore

#scallop_CPPFLAGS = -std=c++11 -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE)
libscallop_a_CPPFLAGS = -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE) -I$(RNACORE_INCLUDE)

libscallop_a_SOURCES = super_graph.h super_graph.cc \
					   config.h config.cc \
					   subsetsum.h subsetsum.cc \
					   router.h router.cc \
					   bundle.h bundle.cc \
					   equation.h equation.cc \
					   gtf.h gtf.cc \
					   scallop.h scallop.cc \
					   previewer.h previewer.cc \
					   filter.h filter.cc \
					   cluster.h cluster.cc \
					   graph_revise.h graph_revise.cc
