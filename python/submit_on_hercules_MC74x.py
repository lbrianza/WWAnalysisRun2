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

name = ["RSGraviton600", "RSGraviton800", "RSGraviton1000", "RSGraviton1200", "RSGraviton1400", "RSGraviton1600", "RSGraviton1800", "RSGraviton2000",
        "RSGraviton2500", "RSGraviton3000", "RSGraviton3500", "RSGraviton4000", "RSGraviton4500", 
        "WJets100", "WJets200", "WJets400", "WJets600", "TTbar", "sch", "tch", "tWch", "tWch_bar", "WW", "ZZ"];
category = ["mu","el"];
xSecWeight = ["4.76735", "1.16691", "0.377865", "0.144482", "0.0616708", "0.0288651", "0.0141334", "0.00751431",
              "0.00167726", "0.000443483", "0.000133915", "0.0000424117", "0.0000130705",
              "1292.", "385.9", "47.9", "19.9", "831.76", "10.11", "216.99", "38.09", "38.09", "63.21", "10.32"];
N = ["32354.", "31906.", "32448.", "32252.", "32275.", "31971.", "32021.", "31295.",
     "32032.", "31374.", "32194.", "32207.", "31551.",
     "10142187.", "5231856.", "1901705.", "1036108.", "42730273.", "984400.", "2966200.", "995600.", "1000000.", "994416.", "996168."];
mass = ["600", "800", "1000", "1200", "1400", "1600", "1800", "2000",
        "2500", "3000", "3500", "4000", "4500", 
        "0", "0", "0", "0", "0", "0", "0", "0", "0", "0", "0"];

for a in range(len(category)):
    for i in range(len(name)):
        fn = "Job/Job_"+name[i]+"_"+category[a];
        outScript = open(fn+".sh","w");
        command = "python python/produceWWNtuples.py -i /gwteray/users/brianza/WWNtupleRun2/ReducedTree_MC74x_puppi/ -n ReducedSelection_"+name[i]+".root -o WWTree_"+name[i]+".root -l "+category[a]+" -w "+xSecWeight[i]+" -no "+N[i]+" -mass "+mass[i];
        print command;
        outScript.write('#!/bin/bash');
        outScript.write("\n"+'cd '+CMSSWDir);
        outScript.write("\n"+'eval `scram runtime -sh`');
        outScript.write("\n"+'cd '+currentDir);
        outScript.write("\n"+"unbuffer "+command+" > "+fn+"_output.txt");
        outScript.close();
        os.system("chmod 777 "+currentDir+"/"+fn+".sh");
        os.system("qsub -V -d "+currentDir+" -q shortcms "+currentDir+"/"+fn+".sh");
