# Makefile for G0 Analysis Code
# (modified version of $(ROOTSYS)/test/Makefile)
#
# March 3, 2000  Author: W. Korsch
# April 25, 2000 	 E. Beise, modify to work with Alpha DUX4.0

ObjSuf        = o
SrcSuf        = C
ExeSuf        =
DllSuf        = so
EVENTLIB      = $(EVENTO)
OutPutOpt     = -o  

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

CXX           = g++
CXXFLAGS      = -Wall -fPIC -O3
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) -lMinuit
GLIBS         = $(ROOTGLIBS)
LIBS         += -lstdc++ -lz


######################################################################
#
#  Generic compilation and linking step to make an executable from
#  a single *.C file
#
%: %.$(SrcSuf)
	$(CXX) $(CXXFLAGS) -o $@ $< $(LIBS) $(GLIBS)

krTotal: krTotal.C
	$(CXX) $(CXXFLAGS) -o krTotal krTotal.C $(LIBS) $(GLIBS)

clean:
		@rm -f $(OBJS) *~  core

cleanall:
		@rm -f $(OBJS) $(EXE) *~  core








