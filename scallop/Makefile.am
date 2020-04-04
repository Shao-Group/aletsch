noinst_LIBRARIES = libscallop.a

GTF_INCLUDE = $(top_srcdir)/lib/gtf
UTIL_INCLUDE = $(top_srcdir)/lib/util
GRAPH_INCLUDE = $(top_srcdir)/lib/graph
RNACORE_INCLUDE = $(top_srcdir)/rnacore

#scallop_CPPFLAGS = -std=c++11 -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE)
libscallop_a_CPPFLAGS = -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE) -I$(RNACORE_INCLUDE)

libscallop_a_SOURCES = config.h config.cc \
					   subsetsum.h subsetsum.cc \
					   router.h router.cc \
					   equation.h equation.cc \
					   scallop.h scallop.cc \
					   filter.h filter.cc \
					   cluster.h cluster.cc \
					   path.h path.cc \
					   hyper_set.h hyper_set.cc \
					   graph_reviser.h graph_reviser.cc
