CXXFLAGS=-Wall # -pedantic -ansi
ROOT_LIB:=`root-config --libs --glibs`
ROOT_FLAGS:=`root-config --cflags --ldflags`
ROOT_INCLUDE:=`root-config --incdir`

DEPS= interface/setOutputTree.h interface/setOutputTreeSynch.h interface/METzCalculator.h interface/analysisUtils.h interface/setInputTree.h
DEPS_OBJ= lib/setOutputTree.o lib/setOutputTreeSynch.o lib/METzCalculator.o lib/analysisUtils.o lib/setInputTree.o

CC = g++
CFLAGS = -Wall

lib/%.o: src/%.cc $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $< $(ROOT_LIB) $(ROOT_FLAGS)

all: produceWWNtuples.exe produceWWSynchNtuples.exe

produceWWNtuples.exe: bin/produceWWNtuples.cpp $(DEPS_OBJ)
	g++ $(CFLAGS) -o $@ $^ $(ROOT_LIB) $(ROOT_FLAGS)

produceWWSynchNtuples.exe: bin/produceWWSynchNtuples.cpp $(DEPS_OBJ)
	g++ $(CFLAGS) -o $@ $^ $(ROOT_LIB) $(ROOT_FLAGS)

clean:
	rm -f lib/*.o
	rm -f lib/*.d
	rm -f *.exe