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

inputFolder = "/eos/cms/store/caf/user/lbrianza/WWReducedTree_run2";
outputFolder = currentDir+"/output/";
exeName = currentDir+"/produceWWNtuples.exe"

dryRun = False;
doMC = True;
doData = True;

category = ["mu","el"];
#category = ["el"];
#category = ["mu"];



samples = [
    #(  118.70000000000,               "WW",        0,    0.),
    #(   16.50000000000,               "ZZ",        0,    0.),
    #(   47.13000000000,               "WZ",        0,    0.),
    (   49.99700000000,          "WW_excl",  1951000,    0.),
    #(   10.71000000000,        "WZ_excl_1", 25996157,    0.), #before neg. event subtraction
    #(   10.71000000000,        "WZ_excl_2", 25996157,    0.), #before neg. event subtraction
    #(   10.71000000000,        "WZ_excl_3", 25996157,    0.), #before neg. event subtraction
    (   10.71000000000,        "WZ_excl_1", 20781321,    0.),
    (   10.71000000000,        "WZ_excl_2", 20781321,    0.),
    (   10.71000000000,        "WZ_excl_3", 20781321,    0.),
    #(    3.22000000000,        "ZZ_excl_1", 15498581,    0.), #before neg. event subtraction
    #(    3.22000000000,        "ZZ_excl_2", 15498581,    0.), #before neg. event subtraction
    #(    3.22000000000,        "ZZ_excl_3", 15498581,    0.), #before neg. event subtraction
    (    3.22000000000,        "ZZ_excl_1", 12641708,    0.),
    (    3.22000000000,        "ZZ_excl_2", 12641708,    0.),
    (    3.22000000000,        "ZZ_excl_3", 12641708,    0.),
    #(   10.11009000000,              "sch",  1000000,    0.), #before neg. event subtraction
    (   10.11009000000,              "sch",  811495,    0.),
    (   43.80000000000,              "tch",        0,    0.),
    (   26.07000000000,          "tch_bar",  1682400,    0.),
    (   35.60000000000,             "tWch",   998400,    0.),
    (   35.60000000000,         "tWch_bar",   985000,    0.),
    #(  831.76000000000, "TTbar_amcatnlo_1", 38475776,    0.),
    #(  831.76000000000, "TTbar_amcatnlo_2", 38475776,    0.),
    #(  831.76000000000, "TTbar_amcatnlo_3", 38475776,    0.),
    #(  831.76000000000, "TTbar_amcatnlo_4", 38475776,    0.),
    #(  831.76000000000, "TTbar_amcatnlo_5", 38475776,    0.),
    (  831.76000000000,   "TTbar_powheg_1", 182123200,   0.),
    (  831.76000000000,   "TTbar_powheg_2", 182123200,   0.),
    (  831.76000000000,   "TTbar_powheg_3", 182123200,   0.),
    (  831.76000000000,   "TTbar_powheg_4", 182123200,   0.),
    (  831.76000000000,   "TTbar_powheg_5", 182123200,   0.),
    (  831.76000000000,   "TTbar_powheg_6", 182123200,   0.),
    (  831.76000000000,   "TTbar_powheg_7", 182123200,   0.),
    (  831.76000000000,   "TTbar_powheg_8", 182123200,   0.),
    (  831.76000000000,   "TTbar_powheg_9", 182123200,   0.),
    (  831.76000000000,  "TTbar_powheg_10", 182123200,   0.),
    (  831.76000000000,  "TTbar_powheg_11", 182123200,   0.),
    (  831.76000000000,  "TTbar_powheg_12", 182123200,   0.),
    (  831.76000000000,  "TTbar_powheg_13", 182123200,   0.),
    (  831.76000000000,  "TTbar_powheg_14", 182123200,   0.),
    (  831.76000000000,  "TTbar_powheg_15", 182123200,   0.),
    (  831.76000000000,  "TTbar_powheg_16", 182123200,   0.),
    (  831.76000000000,  "TTbar_powheg_17", 182123200,   0.),
    (  831.76000000000,  "TTbar_powheg_18", 182123200,   0.),
    (  831.76000000000,  "TTbar_powheg_19", 182123200,   0.),
    (  831.76000000000,  "TTbar_powheg_20", 182123200,   0.),
    (  831.76000000000,  "TTbar_powheg_21", 182123200,   0.),
    (  831.76000000000,  "TTbar_powheg_22", 182123200,   0.),
    (  831.76000000000,  "TTbar_powheg_23", 182123200,   0.),
    #(  831.76000000000, "TTbar_madgraph_1", 10215131,    0.),
    #(  831.76000000000, "TTbar_madgraph_2", 10215131,    0.),
    #(  831.76000000000, "TTbar_madgraph_3", 10215131,    0.),
    #(61526.70000000000,            "WJets",        0,    0.),
    ( 1629.87000000000,       "WJets100_1", 27546978,    0.),
    ( 1629.87000000000,       "WJets100_2", 27546978,    0.),
    ( 1629.87000000000,       "WJets100_3", 27546978,    0.),
    (  435.60000000000,         "WJets200",  4963240,    0.),
    (   59.16900000000,         "WJets400",  1963464,    0.),
    #(   22.71170000000,         "WJets600",        0,    0.),
    (   12.80000000000,      "WJets600bis",  3722395,    0.),
    (    5.26000000000,         "WJets800",  6314257,    0.),
    (    1.33000000000,        "WJets1200",   246737,    0.),
    (    0.03080000009,        "WJets2500",   253561,    0.),
    #(61526.70000000000, "WJets_madgraph_1",        0,    0.),
    #(61526.70000000000, "WJets_madgraph_2",        0,    0.),
    #(61526.70000000000, "WJets_madgraph_3",        0,    0.),
    #(    1.35200000000,         "Higgs650",   399600,  650.),
    #(    0.63980000000,         "Higgs750",    96200,  750.),
    #(    0.12330000000,        "Higgs1000",   400000, 1000.),
    #(    0.27190000000,      "VBFHiggs650",   400000,  650.),
    #(    0.19150000000,      "VBFHiggs750",   100000,  750.),
    #(    0.08732000000,     "VBFHiggs1000",   399634, 1000.),
    (    1.77400363645440040e-01,  "BulkGraviton600",    50000,  600.),
    (    7.22813367744000179e-02, "BulkGraviton700",     48000,  700.),
    (    4.79197963636363647e-02, "BulkGraviton750",     46000,  750.),
    (    3.31548423242400067e-02,  "BulkGraviton800",    48000,  800.),
    (    1.67515749198720032e-02, "BulkGraviton900",     49000,  900.),
    (    8.99252925667200220e-03, "BulkGraviton1000",    50000, 1000.),
    (    2.99686333939200040e-03, "BulkGraviton1200",    50000, 1200.),
    (    1.15022191555440019e-03, "BulkGraviton1400",    50000, 1400.),
    (    4.83220517270400096e-04, "BulkGraviton1600",    50000, 1600.),
    (    2.18919227679360021e-04, "BulkGraviton1800",    50000, 1800.),
    (    2.23940299210069489e-05, "BulkGraviton2000",    50000, 2000.),
    (    4.19924014292041532e-06, "BulkGraviton2500",    50000, 2500.),
    (    9.19110485489025090e-07, "BulkGraviton3000",    50000, 3000.),
    (    2.25464182856213949e-07, "BulkGraviton3500",    50000, 3500.),
    (    5.82648753152116075e-08, "BulkGraviton4000",    49800, 4000.),
    (    0.000000000000000000000, "BulkGraviton4500",    50000, 4500.)
    ]

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
    "data_mu_2016_runB_v2_26",
    "data_mu_2016_runB_v2_27",
    "data_mu_2016_runB_v2_28",
    "data_mu_2016_runB_v2_29"
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
    "data_el_2016_runB_v2_26",
    "data_el_2016_runB_v2_27",
    "data_el_2016_runB_v2_28",
    "data_el_2016_runB_v2_29"
];

nameData = {"el": nameDataEl, "mu":nameDataMu};



for a in range(len(category)):
    
    #MC
    if( doMC ):
        for i in range(len(samples)):
            fn = "jobs/job_"+samples[i][1]+"_"+category[a];
            outScript = open(fn+".sh","w");
            command = "python "+currentDir+"/python/produceWWNtuples.py --exe "+exeName+" -i "+inputFolder+" -n "+str(samples[i][1])+" -o WWTree_"+str(samples[i][1])+" -l "+category[a]+" -w "+str(samples[i][0])+" -no "+str(samples[i][2])+" -mass "+str(samples[i][3])+" --ismc 1 -trig 0";
            print command;
            outScript.write('#!/bin/bash');
            outScript.write("\n"+'cd '+CMSSWDir);
            outScript.write("\n"+'eval `scram runtime -sh`');
            #outScript.write("\n"+'cd '+currentDir);
            #outScript.write("\n"+'source scripts/setup.sh');
            outScript.write("\n"+"cd -");
            outScript.write("\n"+"cp "+currentDir+"/*.txt ./");
            outScript.write("\n"+"cp "+currentDir+"/*.root ./");
            outScript.write("\n"+command);
            outScript.write("\n"+"mv WWTree_"+(samples[i][1])+".root "+outputFolder+"/output_"+category[a]);
            outScript.write("\n");
            outScript.close();
            os.system("chmod 777 "+currentDir+"/"+fn+".sh");
            command2 = "bsub -q cmscaf1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
            print command2
            if( dryRun != True ):
                os.system(command2);
                time.sleep(1)
    
    #data
    if( doData ):
        for i in range(len(nameData[category[a]])):
            fn = "jobs/job_"+(nameData[category[a]])[i]+"_mu";
            outScript = open(fn+".sh","w");
            command = "python "+currentDir+"/python/produceWWNtuples.py --exe "+exeName+" -i "+inputFolder+" -n "+(nameData[category[a]])[i]+" -o WWTree_"+(nameData[category[a]])[i]+" -l mu"+" -w 1. -no 1. -mass 0 --ismc 0 -trig 0";
            print command;
            outScript.write('#!/bin/bash');
            outScript.write("\n"+'cd '+CMSSWDir);
            outScript.write("\n"+'eval `scram runtime -sh`');
            #outScript.write("\n"+'cd '+currentDir);
            #outScript.write("\n"+'source scripts/setup.sh');
            outScript.write("\n"+"cd -");
            outScript.write("\n"+"cp "+currentDir+"/*.txt ./");
            outScript.write("\n"+"cp "+currentDir+"/*.root ./");
            outScript.write("\n"+command);
            outScript.write("\n"+"mv WWTree_"+(nameData[category[a]])[i]+".root "+outputFolder+"/output_"+category[a]);
            outScript.write("\n");
            outScript.close();
            os.system("chmod 777 "+currentDir+"/"+fn+".sh");
            command2 = "bsub -q cmscaf1nd -cwd "+currentDir+" "+currentDir+"/"+fn+".sh";
            print command2
            if( dryRun != True ):
                os.system(command2);
                time.sleep(1)
