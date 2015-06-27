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

name = ["RSGraviton600", "RSGraviton800", "RSGraviton1000", "RSGraviton1200", "RSGraviton1400", "RSGraviton1600", "RSGraviton1800", "RSGraviton2000",
        "RSGraviton2500", "RSGraviton3000", "RSGraviton3500", "RSGraviton4000", "RSGraviton4500",
        "WJets100", "WJets200", "WJets400", "WJets600", "TTbar", "sch", "tch", "tch_bar", "tWch", "tWch_bar", "WW", "ZZ", "WZ",
        "BulkGraviton800", "BulkGraviton1000", "BulkGraviton1200", "BulkGraviton1400", "BulkGraviton1600", "BulkGraviton1800", 
        "BulkGraviton2000", "BulkGraviton2500", "BulkGraviton3000", "BulkGraviton4000", "BulkGraviton4500"];
category = ["mu","el"];
xSecWeight = ["4.76735", "1.16691", "0.377865", "0.144482", "0.0616708", "0.0288651", "0.0141334", "0.00751431",
              "0.00167726", "0.000443483", "0.000133915", "0.0000424117", "0.0000130705",
              "1292.", "385.9", "47.9", "19.9", "831.76", "10.11", "216.99", "80.95", "38.09", "38.09", "63.21", "10.32", "22.82",
              "0.001332687", "0.000359194", "0.000119842", "0.000045798", "0.", "0.", 
              "0.000004197", "0.000000786", "0.000000172", "0.", "0."];
N = ["32354.", "31906.", "32448.", "32252.", "32275.", "31971.", "32021.", "31295.",
     "32032.", "31374.", "32194.", "32207.", "31551.",
     "10142187.", "5231856.", "1901705.", "1036108.", "42730273.", "984400.", "2966200.", "1695400.", "995600.", "1000000.", "994416.", "996168.", "991232.",
     "50000", "50000", "50000", "50000", "49200", "48400", 
     "50000", "48400", "49800", "50000", "50000"];

for a in range(len(category)):
    for i in range(len(name)):
        fn = "Job/Job_"+name[i]+"_"+category[a];
        outScript = open(fn+".sh","w");
        command = "python python/produceWWNtuples_finalSynch.py -i root://eoscms//eos/cms/store/caf/user/lbrianza/WWReducedTree_MC74x/ -n ReducedSelection_"+name[i]+".root -o WWTree_"+name[i]+"_synch.root -l "+category[a]+" -w "+xSecWeight[i]+" -no "+N[i];
        print command;
        outScript.write('#!/bin/bash');
        outScript.write("\n"+'cd '+CMSSWDir);
        outScript.write("\n"+'eval `scram runtime -sh`');
        outScript.write("\n"+'cd '+currentDir);
        outScript.write("\n"+"unbuffer "+command);
        outScript.close();
        os.system("chmod 777 "+currentDir+"/"+fn+".sh");
        command2 = "bsub -q 8nh -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
        os.system(command2);
        print command2
