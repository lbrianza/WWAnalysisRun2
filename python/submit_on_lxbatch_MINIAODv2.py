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
#    "1.352", "0.6398", "0.1233", "0.2719", "0.1915", "0.08732", #https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBSMAt13TeV
    "61526.7", "1629.87","1629.87","1629.87", "435.6", "59.169", "22.7117", 
#    "831.76", "831.76", "831.76", "831.76", "831.76", 
    "10.11", "43.8", "26.07", "35.6", "35.6", "118.7", "16.5", "47.13",
    "49.997", "10.71", "10.71", "10.71", "3.22", "3.22", "3.22", 
    "831.76", "831.76", "831.76", "831.76", "831.76", 
    "831.76", "831.76", "831.76", "831.76", "831.76", 
    "831.76", "831.76", "831.76", "831.76", "831.76", 
    "831.76", "831.76", "831.76", "831.76", "831.76", 
    "831.76", "831.76", "831.76",
#    "831.76", "831.76", "831.76",
    "12.8", "5.26", "1.33", "0.03089", "61526.7", "61526.7", "61526.7", 
    #for the signal - CX(GF) from here: https://github.com/acarvalh/Cross_sections_CMS/blob/master/WED/bulk_KKgrav_decay.txt
    # and BR(WW) from here: https://github.com/acarvalh/Cross_sections_CMS/blob/master/WED/bulk_KKgrav_LHC13.txt
    # you compute the numbers as: CX(GF)*BR(WW)*BR(lvqq)*(k'/0.1)^2, where k'=0.5 here BR(lvqq)=4.38048000000000048e-01
    "1.77400363645440040e-01","3.31548423242400067e-02", "8.99252925667200220e-03", "2.99686333939200040e-03", "1.15022191555440019e-03", "4.83220517270400096e-04", "2.18919227679360021e-04", 
    "2.23940299210069489e-05", "4.19924014292041532e-06", "9.19110485489025090e-07", "2.25464182856213949e-07", "5.82648753152116075e-08","0.",
    "7.22813367744000179e-02","4.79197963636363647e-02","1.67515749198720032e-02",
#    "0.","0.153250","0.102236","0.056357","0.031283","0.017844","0.010458",
 #    "0.006273","0.0018842","0.0006039","0.00019991","0.","0.",
 #    "0.00052133977","0.00010210343", "0.000030476512", "0.", "0.", "0.", "0.", 
 #    "0.", "0.", "0.", "0.", "0.", "0."
              ];

name = [
#    "Higgs650", "Higgs750", "Higgs1000", "VBFHiggs650", "VBFHiggs750", "VBFHiggs1000",
     "WJets","WJets100_1","WJets100_2","WJets100_3", "WJets200", "WJets400", "WJets600", 
 #"TTbar_amcatnlo_1","TTbar_amcatnlo_2","TTbar_amcatnlo_3","TTbar_amcatnlo_4","TTbar_amcatnlo_5", 
     "sch", "tch", "tch_bar", "tWch", "tWch_bar", "WW", "ZZ", "WZ",
     "WW_excl", "WZ_excl_1", "WZ_excl_2", "WZ_excl_3", "ZZ_excl_1", "ZZ_excl_2", "ZZ_excl_3", 
     "TTbar_powheg_1", "TTbar_powheg_2", "TTbar_powheg_3" , "TTbar_powheg_4" , "TTbar_powheg_5" ,
     "TTbar_powheg_6", "TTbar_powheg_7", "TTbar_powheg_8" , "TTbar_powheg_9" , "TTbar_powheg_10" ,
     "TTbar_powheg_11", "TTbar_powheg_12", "TTbar_powheg_13" , "TTbar_powheg_14" , "TTbar_powheg_15" ,
     "TTbar_powheg_16", "TTbar_powheg_17", "TTbar_powheg_18" , "TTbar_powheg_19" , "TTbar_powheg_20" ,
     "TTbar_powheg_21", "TTbar_powheg_22", "TTbar_powheg_23" ,
 #    "TTbar_madgraph_1", "TTbar_madgraph_2", "TTbar_madgraph_3",
     "WJets600bis", "WJets800", "WJets1200", "WJets2500", "WJets_madgraph_1", "WJets_madgraph_2", "WJets_madgraph_3",
    "BulkGraviton600","BulkGraviton800", "BulkGraviton1000", "BulkGraviton1200", "BulkGraviton1400", "BulkGraviton1600", "BulkGraviton1800", 
    "BulkGraviton2000", "BulkGraviton2500", "BulkGraviton3000", "BulkGraviton3500", "BulkGraviton4000", "BulkGraviton4500",
    "BulkGraviton700","BulkGraviton750","BulkGraviton900",
#    "WprimeToWZ600","WprimeToWZ800","WprimeToWZ1000","WprimeToWZ1200","WprimeToWZ1400","WprimeToWZ1600", "WprimeToWZ1800",
 #    "WprimeToWZ2000","WprimeToWZ2500","WprimeToWZ3000","WprimeToWZ3500","WprimeToWZ4000","WprimeToWZ4500",
 #    "VBFBulkGraviton600","VBFBulkGraviton800", "VBFBulkGraviton1000", "VBFBulkGraviton1200", "VBFBulkGraviton1400", "VBFBulkGraviton1600", "VBFBulkGraviton1800", 
 #    "VBFBulkGraviton2000", "VBFBulkGraviton2500", "VBFBulkGraviton3000", "VBFBulkGraviton3500", "VBFBulkGraviton4000", "VBFBulkGraviton4500"
        ];

N = [
#    "399600.", "96200.", "400000.", "400000.", "100000.", "399634.",
     "0.","27546978.","27546978.","27546978.", "4963240.", "1963464.", "0.", 
 #"38475776.", "38475776.", "38475776.", "38475776.", "38475776.", 
     "1000000.", "0.", "1682400.", "998400.", "985000.", "0.", "0.", "0.",
     "1951000.", "25996157.", "25996157.", "25996157.", "15498581.", "15498581.", "15498581.", 
 #    "998400.", "3299200.", "1630900.", "1000000.", "999400.", "988418.", "985600.", "1000000.", #MCATNLO
 #    "1924400.", "25704656.", "25704656.", "25704656.", "15301695.", "15301695.", "15301695.", #MCATNLO
     "182123200.", "182123200.", "182123200.", "182123200.", "182123200.", 
     "182123200.", "182123200.", "182123200.", "182123200.", "182123200.", 
     "182123200.", "182123200.", "182123200.", "182123200.", "182123200.", 
     "182123200.", "182123200.", "182123200.", "182123200.", "182123200.", 
     "182123200.", "182123200.", "182123200.", 
 #    "10215131.", "10215131.", "10215131.",
     "3722395.", "6314257.", "246737.", "253561.", "0.", "0.", "0.",
    "50000.", "48000.", "50000.", "50000.", "50000.", "50000.", "50000.", 
    "50000.", "50000.", "50000.", "50000.", "49800.", "50000.", #RICONTROLLA...
    "48000.","46000","49000",
#    "49600.","50000.","50000.","50000.","50000.","49200.","50000.",
 #    "50000.","49800.","50000.","50000.","50000.","50000.",
 #    "48400.","50000.", "50000.", "50000.", "50000.", "50000.", "50000.", 
 #    "50000.", "50000.", "50000.", "49200.", "50000.", "50000."
    ];

mass = [
#    "650", "750", "1000", "650", "750", "1000",
     "0","0","0","0", "0", "0", "0", 
 #"0", "0", "0", "0", "0", 
     "0", "0", "0", "0", "0", "0", "0", "0",
     "0", "0", "0", "0", "0", "0", "0", 
     "0", "0", "0", "0", "0", 
     "0", "0", "0", "0", "0", 
     "0", "0", "0", "0", "0", 
     "0", "0", "0", "0", "0", 
     "0", "0", "0",
 #    "0", "0", "0",
     "0", "0", "0", "0", "0", "0", "0", 
    "600", "800", "1000", "1200", "1400", "1600", "1800",
    "2000", "2500", "3000", "3500", "4000", "4500",
    "700","750","900",
#    "600", "800", "1000", "1200", "1400", "1600", "1800",
 #    "2000","2500", "3000", "3500", "4000", "4500",
 #    "600", "800", "1000", "1200", "1400", "1600", "1800",
 #    "2000", "2500", "3000", "3500", "4000", "4500"
    ];

#nameDataMu = ['data_mu_prompt_v4_25ns_runD_4']


nameDataMu = [
    "data_mu_2016_runB_v2_1", "data_mu_2016_runB_v2_2", "data_mu_2016_runB_v2_3", "data_mu_2016_runB_v2_4", "data_mu_2016_runB_v2_5",
    "data_mu_2016_runB_v2_6", "data_mu_2016_runB_v2_7", "data_mu_2016_runB_v2_8", "data_mu_2016_runB_v2_9", "data_mu_2016_runB_v2_10",
    "data_mu_2016_runB_v2_11", "data_mu_2016_runB_v2_12", "data_mu_2016_runB_v2_13", "data_mu_2016_runB_v2_14", "data_mu_2016_runB_v2_15",
    "data_mu_2016_runB_v2_16", "data_mu_2016_runB_v2_17", "data_mu_2016_runB_v2_18", "data_mu_2016_runB_v2_19", "data_mu_2016_runB_v2_20",
    "data_mu_2016_runB_v2_21", "data_mu_2016_runB_v2_22", "data_mu_2016_runB_v2_23", "data_mu_2016_runB_v2_24", "data_mu_2016_runB_v2_25",
    "data_mu_2016_runB_v2_26", "data_mu_2016_runB_v2_27", "data_mu_2016_runB_v2_28", "data_mu_2016_runB_v2_29"
];

nameDataEl = [
    "data_el_2016_runB_v2_1", "data_el_2016_runB_v2_2", "data_el_2016_runB_v2_3", "data_el_2016_runB_v2_4", "data_el_2016_runB_v2_5",
    "data_el_2016_runB_v2_6", "data_el_2016_runB_v2_7", "data_el_2016_runB_v2_8", "data_el_2016_runB_v2_9", "data_el_2016_runB_v2_10",
    "data_el_2016_runB_v2_11", "data_el_2016_runB_v2_12", "data_el_2016_runB_v2_13", "data_el_2016_runB_v2_14", "data_el_2016_runB_v2_15",
    "data_el_2016_runB_v2_16", "data_el_2016_runB_v2_17", "data_el_2016_runB_v2_18", "data_el_2016_runB_v2_19", "data_el_2016_runB_v2_20",
    "data_el_2016_runB_v2_21", "data_el_2016_runB_v2_22", "data_el_2016_runB_v2_23", "data_el_2016_runB_v2_24", "data_el_2016_runB_v2_25",
    "data_el_2016_runB_v2_26", "data_el_2016_runB_v2_27", "data_el_2016_runB_v2_28", "data_el_2016_runB_v2_29"
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
        time.sleep(3)

#data mu
for i in range(len(nameDataMu)):
        fn = "Job/Job_"+nameDataMu[i]+"_mu";
        outScript = open(fn+".sh","w");
        command = "python python/produceWWNtuples.py -n "+nameDataMu[i]+" -o WWTree_"+nameDataMu[i]+" -l mu"+" -w 1. -no 1. -mass 0 --ismc 0 -trig 1";
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
        time.sleep(3)

#data el
for i in range(len(nameDataEl)):
        fn = "Job/Job_"+nameDataEl[i]+"_el";
        outScript = open(fn+".sh","w");
        command = "python python/produceWWNtuples.py -n "+nameDataEl[i]+" -o WWTree_"+nameDataEl[i]+" -l el"+" -w 1. -no 1. -mass 0 --ismc 0 -trig 1";
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
        time.sleep(3)

