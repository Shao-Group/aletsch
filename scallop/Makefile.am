noinst_LIBRARIES = libscallop.a

GTF_INCLUDE = $(top_srcdir)/gtf
UTIL_INCLUDE = $(top_srcdir)/util
GRAPH_INCLUDE = $(top_srcdir)/graph
RNACORE_INCLUDE = $(top_srcdir)/rnacore

libscallop_a_CPPFLAGS = -std=c++11 -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE) -I$(RNACORE_INCLUDE)

libscallop_a_SOURCES = subsetsum.h subsetsum.cc \
					   router.h router.cc \
					   equation.h equation.cc \
					   scallop.h scallop.cc \
					   cluster.h cluster.cc \
					   path.h path.cc \
					   hyper_set.h hyper_set.cc
