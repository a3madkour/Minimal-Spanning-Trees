CXX = g++
CXXFLAGS = -std=c++11 
BINARIES = NCuts test DPCuts 
LAPACKFLAGS = -llapacke -lblas

all: $(BINARIES)

NCuts: NCuts.o Graph.o
	$(CXX)  $(CXXFLAGS) $(LAPACKFLAGS) -o NCuts NCuts.o Graph.o 

NCuts.o: NCuts.cpp Graph.h
	$(CXX) $(CXXFLAGS) -c NCuts.cpp

Graph.o: Graph.h Graph.cpp
	$(CXX) $(CXXFLAGS) -c Graph.cpp

clean:
	$(RM) $(BINARIES) *.o

