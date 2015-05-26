#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import subprocess

currentDir = os.getcwd();
CMSSWDir = currentDir+"../";
ReducedTreeDir = "";

name = ["WJets100", "WJets200", "WJets400", "WJets600", "TTbar","RSGraviton1000", "RSGraviton2000", "RSGraviton3000", "RSGraviton4000", "ZJets", "sch", "sch_bar", "tch", "tch_bar", "tWch", "tWch_bar","WW"];
category = ["mu","el"];
xSecWeight = ["1817.", "471.6", "55.61", "18.81","831.76", "2.37","0.04797","0.00292","0.0002739", "1.000", "7.2", "4.16", "136.05", "80.97", "35.6", "35.6","0.29"];
N = ["5262265.","4936077.","4640594.","4581841.","25446993.","29205.","29854","28983","28687.","1.","500000.","250000.","3991000.","1999800.","986100.","971800.","19800"];

for a in range(len(category)):
    for i in range(len(name)):
        fn = "Job/Job_"+name[i]+"_"+category[a];
        outScript = open(fn+".sh","w");
        command = "python python/produceWWNtuples.py -i /gwteray/users/brianza/WWNtupleRun2/ReducedTree/ -n ReducedSelection_"+name[i]+".root -o WWTree_"+name[i]+".root -l "+category[a]+" -w "+xSecWeight[i]+" -no "+N[i];
        print command;
        outScript.write('#!/bin/bash');
        outScript.write("\n"+'cd '+CMSSWDir);
        outScript.write("\n"+'eval `scram runtime -sh`');
        outScript.write("\n"+'cd '+currentDir);
        outScript.write("\n"+"unbuffer "+command+" > "+fn+"_output.txt");
        outScript.close();
        os.system("chmod 777 "+currentDir+"/"+fn+".sh");
        os.system("qsub -V -d "+currentDir+" -q shortcms "+currentDir+"/"+fn+".sh");
