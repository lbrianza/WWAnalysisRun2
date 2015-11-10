CXXFLAGS=-Wall # -pedantic -ansi
ROOT_LIB:=`root-config --libs --glibs`
ROOT_FLAGS:=`root-config --cflags --ldflags`
ROOT_INCLUDE:=`root-config --incdir`

DEPS= interface/setOutputTree.h interface/setOutputTreeSynch.h interface/METzCalculator.h interface/analysisUtils.h interface/setInputTree.h interface/METzCalculator_Run2.h interface/pseudodataNtuple.h interface/setInputPseudodata.h interface/PUWeight.h
DEPS_OBJ= lib/setOutputTree.o lib/setOutputTreeSynch.o lib/METzCalculator.o lib/analysisUtils.o lib/setInputTree.o lib/METzCalculator_Run2.o lib/pseudodataNtuple.o lib/setInputPseudodata.o lib/PUWeight.o

CC = g++
CFLAGS = -Wall

lib/%.o: src/%.cc $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $< $(ROOT_LIB) $(ROOT_FLAGS)

all: produceWWNtuples.exe produceWWSynchNtuples.exe ntupleConverter.exe produceWWNtuples_finalSynch.exe produceWWpseudodata.exe PUWeightCalculator.exe synchNtupleConverter.exe

synchNtupleConverter.exe: bin/synchNtupleConverter.cpp $(DEPS_OBJ)
	g++ $(CFLAGS) -o $@ $^ $(ROOT_LIB) $(ROOT_FLAGS)

PUWeightCalculator.exe: bin/PUWeightCalculator.cpp $(DEPS_OBJ)
	g++ $(CFLAGS) -o $@ $^ $(ROOT_LIB) $(ROOT_FLAGS)

produceWWNtuples.exe: bin/produceWWNtuples.cpp $(DEPS_OBJ)
	g++ $(CFLAGS) -o $@ $^ $(ROOT_LIB) $(ROOT_FLAGS)

produceWWNtuples_finalSynch.exe: bin/produceWWNtuples_finalSynch.cpp $(DEPS_OBJ)
	g++ $(CFLAGS) -o $@ $^ $(ROOT_LIB) $(ROOT_FLAGS)

produceWWSynchNtuples.exe: bin/produceWWSynchNtuples.cpp $(DEPS_OBJ)
	g++ $(CFLAGS) -o $@ $^ $(ROOT_LIB) $(ROOT_FLAGS)

ntupleConverter.exe: bin/ntupleConverter.cpp $(DEPS_OBJ)
	g++ $(CFLAGS) -o $@ $^ $(ROOT_LIB) $(ROOT_FLAGS)

produceWWpseudodata.exe: bin/produceWWpseudodata.cpp $(DEPS_OBJ)
	g++ $(CFLAGS) -o $@ $^ $(ROOT_LIB) $(ROOT_FLAGS)

clean:
	rm -f lib/*.o
	rm -f lib/*.d
	rm -f *.exe