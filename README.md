# WWAnalysisRun2

The package contains a code to produce ntuple for WW semileptonic final state.
It takes in input ntuples produced from miniAOD with the TreeMaker at this link: https://github.com/lbrianza/RA2_2014


Instructions:

git clone https://github.com/lbrianza/WWAnalysisRun2;

cd WWAnalysisRun2/;

make;

python python/produceWWNtuples.py -n ReducedSelection_RSGraviton4000.root -o RSGraviton4000.root -mc True