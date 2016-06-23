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

category = ["mu","el"];
xSecWeight = [
    #for the signal - CX(GF) from here: https://github.com/acarvalh/Cross_sections_CMS/blob/master/WED/bulk_KKgrav_decay.txt
    # and BR(WW) from here: https://github.com/acarvalh/Cross_sections_CMS/blob/master/WED/bulk_KKgrav_LHC13.txt
    # you compute the numbers as: CX(GF)*BR(WW)*BR(lvqq)*(k'/0.1)^2, where k'=0.5 here BR(lvqq)=4.38048000000000048e-01
    "1.77400363645440040e-01","3.31548423242400067e-02", "8.99252925667200220e-03", "2.99686333939200040e-03", "1.15022191555440019e-03", "4.83220517270400096e-04", "2.18919227679360021e-04", 
    "2.23940299210069489e-05", "4.19924014292041532e-06", "9.19110485489025090e-07", "2.25464182856213949e-07", "5.82648753152116075e-08","0.",
    "7.22813367744000179e-02","4.79197963636363647e-02","1.67515749198720032e-02",
              ];

name = [
    "BulkGraviton600","BulkGraviton800", "BulkGraviton1000", "BulkGraviton1200", "BulkGraviton1400", "BulkGraviton1600", "BulkGraviton1800", 
    "BulkGraviton2000", "BulkGraviton2500", "BulkGraviton3000", "BulkGraviton3500", "BulkGraviton4000", "BulkGraviton4500",
    "BulkGraviton700","BulkGraviton750","BulkGraviton900",
        ];

N = [
    "50000.", "48000.", "50000.", "50000.", "50000.", "50000.", "50000.", 
    "50000.", "50000.", "50000.", "50000.", "49800.", "50000.", #RICONTROLLA...
    "48000.","46000","49000",
    ];

mass = [
    "600", "800", "1000", "1200", "1400", "1600", "1800",
    "2000", "2500", "3000", "3500", "4000", "4500",
    "700","750","900",
    ];

#nameDataMu = ['data_mu_prompt_v4_25ns_runD_4']


nameDataMu = [
    "data_mu_2016_runB_v2_1",
    "data_mu_2016_runB_v2_2",
    "data_mu_2016_runB_v2_3",
    "data_mu_2016_runB_v2_4",
    "data_mu_2016_runB_v2_5",
    "data_mu_2016_runB_v2_6",
    "data_mu_2016_runB_v2_7",
    "data_mu_2016_runB_v2_8",
    "data_mu_2016_runB_v2_9",
    "data_mu_2016_runB_v2_10",
    "data_mu_2016_runB_v2_11",
    "data_mu_2016_runB_v2_12",
    "data_mu_2016_runB_v2_13",
    "data_mu_2016_runB_v2_14",
    "data_mu_2016_runB_v2_15",
    "data_mu_2016_runB_v2_16",
    "data_mu_2016_runB_v2_17",
    "data_mu_2016_runB_v2_18",
    "data_mu_2016_runB_v2_19",
    "data_mu_2016_runB_v2_20",
    "data_mu_2016_runB_v2_21",
    "data_mu_2016_runB_v2_22",
    "data_mu_2016_runB_v2_23",
    "data_mu_2016_runB_v2_24",
    "data_mu_2016_runB_v2_25",
    "data_mu_2016_runB_v2_26"
];

nameDataEl = [
    "data_el_2016_runB_v2_1",
    "data_el_2016_runB_v2_2",
    "data_el_2016_runB_v2_3",
    "data_el_2016_runB_v2_4",
    "data_el_2016_runB_v2_5",
    "data_el_2016_runB_v2_6",
    "data_el_2016_runB_v2_7",
    "data_el_2016_runB_v2_8",
    "data_el_2016_runB_v2_9",
    "data_el_2016_runB_v2_10",
    "data_el_2016_runB_v2_11",
    "data_el_2016_runB_v2_12",
    "data_el_2016_runB_v2_13",
    "data_el_2016_runB_v2_14",
    "data_el_2016_runB_v2_15",
    "data_el_2016_runB_v2_16",
    "data_el_2016_runB_v2_17",
    "data_el_2016_runB_v2_18",
    "data_el_2016_runB_v2_19",
    "data_el_2016_runB_v2_20",
    "data_el_2016_runB_v2_21",
    "data_el_2016_runB_v2_22",
    "data_el_2016_runB_v2_23",
    "data_el_2016_runB_v2_24",
    "data_el_2016_runB_v2_25",
    "data_el_2016_runB_v2_26"
];


#MC
for a in range(len(category)):
    for i in range(len(name)):
        fn = "Job/Job_"+name[i]+"_"+category[a];
        outScript = open(fn+".sh","w");
        command = "python python/produceWWNtuples.py -n "+name[i]+" -o WWTree_"+name[i]+" -l "+category[a]+" -w "+xSecWeight[i]+" -no "+N[i]+" -mass "+mass[i]+" --ismc 1 -trig 0";
        print command;
        outScript.write('#!/bin/bash');
        outScript.write("\n"+'cd '+CMSSWDir);
        outScript.write("\n"+'eval `scram runtime -sh`');
        outScript.write("\n"+'cd '+currentDir);
        outScript.write("\n"+command);
        outScript.close();
        os.system("chmod 777 "+currentDir+"/"+fn+".sh");
#        command2 = "bsub -q cmscaf1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
        command2 = "bsub -q cmscaf1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
        os.system(command2);
        print command2


'''
#data mu
for i in range(len(nameDataMu)):
        fn = "Job/Job_"+nameDataMu[i]+"_mu";
        outScript = open(fn+".sh","w");
        command = "python python/produceWWNtuples.py -n "+nameDataMu[i]+" -o WWTree_"+nameDataMu[i]+" -l mu"+" -w 1. -no 1. -mass 0 --ismc 0 -trig 0";
        print command;
        outScript.write('#!/bin/bash');
        outScript.write("\n"+'cd '+CMSSWDir);
        outScript.write("\n"+'eval `scram runtime -sh`');
        outScript.write("\n"+'cd '+currentDir);
        outScript.write("\n"+command);
        outScript.close();
        os.system("chmod 777 "+currentDir+"/"+fn+".sh");
#        command2 = "bsub -q cmscaf1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
        command2 = "bsub -q cmscaf1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
        os.system(command2);
        print command2


#data el
for i in range(len(nameDataEl)):
        fn = "Job/Job_"+nameDataEl[i]+"_el";
        outScript = open(fn+".sh","w");
        command = "python python/produceWWNtuples.py -n "+nameDataEl[i]+" -o WWTree_"+nameDataEl[i]+" -l el"+" -w 1. -no 1. -mass 0 --ismc 0 -trig 0";
        print command;
        outScript.write('#!/bin/bash');
        outScript.write("\n"+'cd '+CMSSWDir);
        outScript.write("\n"+'eval `scram runtime -sh`');
        outScript.write("\n"+'cd '+currentDir);
        outScript.write("\n"+command);
        outScript.close();
        os.system("chmod 777 "+currentDir+"/"+fn+".sh");
#        command2 = "bsub -q cmscaf1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
        command2 = "bsub -q cmscaf1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
        os.system(command2);
        print command2
'''
