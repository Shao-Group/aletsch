noinst_LIBRARIES = libbridge.a

GTF_INCLUDE = $(top_srcdir)/lib/gtf
UTIL_INCLUDE = $(top_srcdir)/lib/util
GRAPH_INCLUDE = $(top_srcdir)/lib/graph
RNACORE_INCLUDE = $(top_srcdir)/rnacore

libbridge_a_CPPFLAGS = -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE) -I$(RNACORE_INCLUDE)

libbridge_a_SOURCES = fragment.h fragment.cc \
					  fcluster.h fcluster.cc \
					  path.h path.cc \
					  config.h config.cc \
					  bridger.h bridger.cc \
					  pier.h pier.cc 
