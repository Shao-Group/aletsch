noinst_LIBRARIES = libgraph.a

UTILDIR = $(top_srcdir)/util

libgraph_a_CPPFLAGS = -I$(UTILDIR) -std=c++11 -O2
libgraph_a_SOURCES = graph_base.cc graph_base.h \
					 directed_graph.cc directed_graph.h \
					 edge_base.cc edge_base.h \
					 undirected_graph.cc undirected_graph.h \
					 vertex_base.cc vertex_base.h \
					 draw.h draw.cc
