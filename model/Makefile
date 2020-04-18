#!/bin/bash

# include src/snap/Makefile.config
# include src/snap/examples/node2vec/Makefile.ex
# include src/snap/examples/Makefile.exmain


#
#	Global configuration for SNAP makefiles
#

GLIB = src/glib-core
SNAP = src/snap-core
GLIBADV = src/glib-adv
SNAPADV = src/snap-adv
SNAPEXP = src/snap-exp

CGLIB = ./$(GLIB)
CSNAP = ./$(SNAP)

EXGLIB = ./$(GLIB)
EXSNAP = ./$(SNAP)
EXGLIBADV = ./$(GLIBADV)
EXSNAPADV = ./$(SNAPADV)
EXSNAPEXP = ./$(SNAPEXP)

UNAME := $(shell uname)


# ifeq ($(UNAME), Linux)
#   # Linux flags
#   CC = g++
#   CXXFLAGS += -std=c++11 -Wall
#   CXXFLAGS += -O3 -DNDEBUG -fopenmp
#   # turn on for crash debugging, get symbols with <prog> 2>&1 | c++filt
#   #CXXFLAGS += -g -rdynamic
#   #CXXFLAGS += -ggdb
#   # turn on for OpenMP
#   CXXOPENMP = 
#   LDFLAGS +=
#   LIBS += -lrt

# else ifeq ($(UNAME), Darwin)
#   # OS X flags
#   CC = g++
#   CXXFLAGS += -std=c++98 -Wall -Wno-unknown-pragmas
#   CXXFLAGS += -O3 -DNDEBUG
#   CLANG := $(shell g++ -v 2>&1 | grep clang | cut -d " " -f 2)
#   ifneq ($(CLANG), LLVM)
#     CXXFLAGS += -fopenmp
#     #CXXOPENMP += -fopenmp
#   else
#     CXXFLAGS += -DNOMP
#     CXXOPENMP =
#   endif
#   LDFLAGS +=
#   LIBS +=

# else ifeq ($(shell uname -o), Cygwin)
#   # Cygwin flags
#   CC = g++
#   CXXFLAGS += -Wall -D__STDC_LIMIT_MACROS
#   CXXFLAGS += -O3 -DNDEBUG
#   CXXOPENMP = -fopenmp
#   LDFLAGS +=
#   LIBS +=

# endif


CXX=g++
CXXFLAGS= -ansi -O5 -ggdb3 -Wall
CXXFLAGS += -std=c++11 -Wall
CXXFLAGS += -O3 -DNDEBUG -fopenmp
CXXOPENMP = 
LDFLAGS +=
LIBS += -lrt
DIRSRC= ./src/
EXEC=louvain convert hierarchy matrix
OBJ1= $(DIRSRC)graph_binary.o $(DIRSRC)louvain.o $(DIRSRC)quality.o $(DIRSRC)modularity.o $(DIRSRC)zahn.o $(DIRSRC)owzad.o $(DIRSRC)goldberg.o $(DIRSRC)condora.o $(DIRSRC)devind.o $(DIRSRC)devuni.o $(DIRSRC)dp.o $(DIRSRC)shimalik.o $(DIRSRC)balmod.o
OBJ2= $(DIRSRC)graph.o
DEPH = $(EXSNAPADV)/n2v.h $(EXSNAPADV)/word2vec.h $(EXSNAPADV)/biasedrandomwalk.h
DEPCPP = $(EXSNAPADV)/n2v.cpp $(EXSNAPADV)/word2vec.cpp $(EXSNAPADV)/biasedrandomwalk.cpp
CXXFLAGS += $(CXXOPENMP)

all: $(EXEC)

louvain : $(OBJ1) $(DIRSRC)main_louvain.o $(DEPCPP) $(EXSNAP)/Snap.o 
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ1) $(DIRSRC)main_louvain.o $(DEPCPP) $(EXSNAP)/Snap.o -I$(EXSNAP) -I$(EXSNAPADV) -I$(EXGLIB) -I$(EXSNAPEXP) $(LDFLAGS) $(LIBS)


convert : $(OBJ2) $(DIRSRC)main_convert.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

hierarchy : $(DIRSRC)main_hierarchy.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

matrix : $(DIRSRC)main_matrix.o
	$(CXX) -o $@ $^ $(CXXFLAGS)

##########################################
# Generic rules
##########################################


$(EXSNAP)/Snap.o: 
	make -C $(EXSNAP)


# $(DIRSRC)main_louvain.o: $(DIRSRC)main_louvain.cpp $(DIRSRC)main_louvain.h
# 	$(CXX) -o  $@ -c $< $(CXXFLAGS)

$(DIRSRC)main_louvain.o: $(DIRSRC)main_louvain.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS)



# $(DIRSRC)main_convert.o: $(DIRSRC)main_convert.cpp $(DIRSRC)main_convert.h
# 	$(CXX) -o  $@ -c $< $(CXXFLAGS)

$(DIRSRC)main_convert.o: $(DIRSRC)main_convert.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS)



# $(DIRSRC)main_hierarchy.o: $(DIRSRC)main_hierarchy.cpp $(DIRSRC)main_hierarchy.h
# 	$(CXX) -o  $@ -c $< $(CXXFLAGS)

$(DIRSRC)main_hierarchy.o: $(DIRSRC)main_hierarchy.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS)


# $(DIRSRC)main_matrix.o: $(DIRSRC)main_matrix.cpp $(DIRSRC)main_matrix.h
# 	$(CXX) -o  $@ -c $< $(CXXFLAGS)

$(DIRSRC)main_matrix.o: $(DIRSRC)main_matrix.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS)


%.o: %.cpp %.h
	$(CXX) -o  $@ -c $< $(CXXFLAGS)

%.o: %.cpp
	$(CXX) -o $@ -c $< $(CXXFLAGS)

clean:
	rm -f $(DIRSRC)*.o

mrproper: clean
	rm -f *~ $(EXEC)
