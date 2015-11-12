#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import subprocess

currentDir = os.getcwd();
CMSSWDir = currentDir+"/../";
ReducedTreeDir = "";
ntupleDir = currentDir+"/output/output_"

name = ["RSGraviton600", "RSGraviton800", "RSGraviton1000", "RSGraviton1200", "RSGraviton1400", "RSGraviton1600", "RSGraviton1800", "RSGraviton2000",
        "RSGraviton2500", "RSGraviton3000", "RSGraviton3500", "RSGraviton4000", "RSGraviton4500",
        "WJets","WJets100", "WJets200", "WJets400", "WJets600", "WJets600bis", "WJets800", "WJets1200", "WJets2500", "TTbar_amcatnlo", "sch", "tch", "tch_bar", "tWch", "tWch_bar", "WW", "ZZ", "WZ",
        "WW_excl", "WZ_excl", "ZZ_excl", "TTbar_powheg", "TTbar_madgraph",
#        "WJets_50ns", "WW_50ns", "WZ_50ns", "ZZ_50ns", "TTbar_amcatnlo_50ns", "TTbar_powheg_50ns", "TTbar_madgraph_50ns",
#        "tch_50ns", "tch_bar_50ns", "tWch_bar_50ns", "tWch_50ns",
        "BulkGraviton800", "BulkGraviton1000", "BulkGraviton1200", "BulkGraviton1400", "BulkGraviton1600", "BulkGraviton1800", 
        "BulkGraviton2000", "BulkGraviton2500", "BulkGraviton3000", "BulkGraviton4000", "BulkGraviton4500", "data", "data_silver"];

category = ["mu","el"];

#MC
for a in range(len(category)):
    for i in range(len(name)):
        fn = "Job_synch/Job_"+name[i]+"_"+category[a];
        outScript = open(fn+".sh","w");
        command = "python python/synchNtupleConverter.py -n "+name[i]+" -l "+category[a]+" -o "+name[i]+" -i "+ntupleDir;
        print command;
        outScript.write('#!/bin/bash');
        outScript.write("\n"+'cd '+CMSSWDir);
        outScript.write("\n"+'eval `scram runtime -sh`');
        outScript.write("\n"+'cd '+currentDir);
        outScript.write("\n"+command);
        outScript.close();
        os.system("chmod 777 "+currentDir+"/"+fn+".sh");
        command2 = "bsub -q cmscaf1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
        os.system(command2);
        print command2
                
