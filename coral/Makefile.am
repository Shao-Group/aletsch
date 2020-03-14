noinst_LIBRARIES = libcoral.a

GTF_INCLUDE = $(top_srcdir)/lib/gtf
UTIL_INCLUDE = $(top_srcdir)/lib/util
GRAPH_INCLUDE = $(top_srcdir)/lib/graph
RNACORE_INCLUDE = $(top_srcdir)/rnacore

libcoral_a_CPPFLAGS = -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE) -I$(RNACORE_INCLUDE)

libcoral_a_SOURCES = config.h config.cc \
					 hit.h hit.cc \
					 fragment.h fragment.cc \
					 fcluster.h fcluster.cc \
					 bundle_base.h bundle_base.cc \
					 bundle2.h bundle2.cc \
					 path.h path.cc \
					 previewer.h previewer.cc \
					 bridger.h bridger.cc \
					 reference.h reference.cc \
					 coral.h coral.cc \
					 hyper_node.h hyper_node.cc
