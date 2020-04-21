noinst_LIBRARIES = libbridge.a

GTF_INCLUDE = $(top_srcdir)/gtf
UTIL_INCLUDE = $(top_srcdir)/util
GRAPH_INCLUDE = $(top_srcdir)/graph
RNACORE_INCLUDE = $(top_srcdir)/rnacore

libbridge_a_CPPFLAGS = -std=c++11 -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE) -I$(RNACORE_INCLUDE) -I$(GTF_INCLUDE)

libbridge_a_SOURCES = bridge_path.h bridge_path.cc \
					  pier.h pier.cc \
					  bridge_solver.h bridge_solver.cc
