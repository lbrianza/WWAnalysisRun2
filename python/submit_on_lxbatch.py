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

name = ["WJets","TTbar","RSGraviton1000"];
category = ["mu","el"];

for a in range(len(category)):
    for i in range(len(name)):
        fn = "Job/Job_"+name[i]+"_"+category[a];
        outScript = open(fn+".sh","w");
        command = "python python/produceWWNtuples.py -n ReducedSelection_"+name[i]+".root -o WWTree_"+name[i]+"_"+category[a]+".root -l "+category[a];
        outScript.write('#!/bin/bash');
        outScript.write("\n"+'cd '+CMSSWDir);
        outScript.write("\n"+'eval `scram runtime -sh`');
        outScript.write("\n"+'cd '+currentDir);
        outScript.write("\n"+"unbuffer "+command+" > "+fn+"_output.txt");
        outScript.close();
        os.system("chmod 777 "+currentDir+"/"+fn+".sh");
    	os.system("bsub -q 8nh -cwd "+currentDir+" "+fn+".sh");

