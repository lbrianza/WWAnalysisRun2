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
    "1.0149", "0.11877", "0.320", "0.125",
    "4.76735", "1.16691", "0.377865", "0.144482", "0.0616708", "0.0141334", "0.00751431",
              "0.00167726", "0.000443483", "0.000133915", "0.0000424117", "0.0000130705",
              "61526.7", "1629.87", "435.6", "59.169", "22.7117", "831.76", "831.76", "831.76", "831.76", "831.76", 
    "10.11", "43.8", "26.07", "35.6", "35.6", "118.7", "16.5", "47.13",
              "49.997", "10.71", "10.71", "10.71", "3.22", "3.22", "3.22", 
    "831.76", "831.76", "831.76", "831.76", "831.76", "831.76", "831.76", "831.76",
    "12.8", "5.26", "1.33", "0.03089", "61526.7", "61526.7", "61526.7", 
    "0.", "0.", "0.", "0.",
#              "61526.7", "118.7", "15.4", "66.1", "831.76", "831.76", "831.76",
#              "43.8", "26.07", "35.6", "35.6",
              "0.406830851","0.076058293", "0.020499693", "0.000119842", "0.000045798", "0.", "0.",
              "0.000000786", "0.", "0.000000172", "0.", "0.","0.",
              "0.","0.153250","0.102236","0.056357","0.031283","0.017844","0.010458",
              "0.006273","0.0018842","0.0006039","0.00019991","0.","0.",
    "0.165094354","0.110500282","0.0382161"
              ];

name = [
    "Higgs650", "Higgs1000", "VBFHiggs650", "VBFHiggs1000",
    "RSGraviton600", "RSGraviton800", "RSGraviton1000", "RSGraviton1200", "RSGraviton1400", "RSGraviton1800", "RSGraviton2000",
        "RSGraviton2500", "RSGraviton3000", "RSGraviton3500", "RSGraviton4000", "RSGraviton4500",
        "WJets","WJets100", "WJets200", "WJets400", "WJets600", "TTbar_amcatnlo_1","TTbar_amcatnlo_2","TTbar_amcatnlo_3","TTbar_amcatnlo_4","TTbar_amcatnlo_5", 
    "sch", "tch", "tch_bar", "tWch", "tWch_bar", "WW", "ZZ", "WZ",
        "WW_excl", "WZ_excl_1", "WZ_excl_2", "WZ_excl_3", "ZZ_excl_1", "ZZ_excl_2", "ZZ_excl_3", 
    "TTbar_powheg_1", "TTbar_powheg_2", "TTbar_powheg_3" , "TTbar_powheg_4" , "TTbar_powheg_5" ,"TTbar_madgraph_1", "TTbar_madgraph_2", "TTbar_madgraph_3",
    "WJets600bis", "WJets800", "WJets1200", "WJets2500", "WJets_madgraph_1", "WJets_madgraph_2", "WJets_madgraph_3", 
    "TTbar600", "TTbar800", "TTbar1200", "TTbar2500",
#        "WJets_50ns", "WW_50ns", "WZ_50ns", "ZZ_50ns", "TTbar_amcatnlo_50ns", "TTbar_powheg_50ns", "TTbar_madgraph_50ns",
#        "tch_50ns", "tch_bar_50ns", "tWch_bar_50ns", "tWch_50ns",
        "BulkGraviton600","BulkGraviton800", "BulkGraviton1000", "BulkGraviton1200", "BulkGraviton1400", "BulkGraviton1600", "BulkGraviton1800", 
        "BulkGraviton2000", "BulkGraviton2500", "BulkGraviton3000", "BulkGraviton3500", "BulkGraviton4000", "BulkGraviton4500",
        "WprimeToWZ600","WprimeToWZ800","WprimeToWZ1000","WprimeToWZ1200","WprimeToWZ1400","WprimeToWZ1600", "WprimeToWZ1800",
        "WprimeToWZ2000","WprimeToWZ2500","WprimeToWZ3000","WprimeToWZ3500","WprimeToWZ4000","WprimeToWZ4500",
    "BulkGraviton700","BulkGraviton750","BulkGraviton900"
        ];

N = [
    "399600.", "400000.", "398400.", "400000.",
    "32354.", "31906.", "32448.", "32252.", "32275.", "32021.", "31295.",
     "32032.", "31374.", "32194.", "32207.", "32351.",
     "24184766.","10152718.", "5221599.", "1745914.", "1039152.", "42784971.", "42784971.", "42784971.", "42784971.", "42784971.", 
    "613384.", "3299800.", "1680200.", "995600.", "988500.", "993640.", "996944.", "978512.", #    "984400.", "3299800.", "1680200.", "995600.", "988500.", "993640.", "996944.", "978512.", #MCATNLO samples
     "1951600.", "14346866.", "14346866.", "14346866.", "11863244.", "11863244.", "11863244.", #"1951600.", "24714550.", "24714550.", "24714550.", "18790122.", "18790122.", "18790122.", #MCATNLO samples 
    "19757190.", "19757190.", "19757190.", "19757190.", "19757190.", "11344206.", "11344206.", "11344206.",
    "4041997.", "1574633.", "255637.", "253036.", "72207128.", "72207128.", "72207128.", 
    "5119009.", "3510897.", "1014678.", "507842.",
#     "24089991.", "989608.", "998848.", "996920.", "4994250.", "19665194.", "4992231.",
#     "1273800.", "681900.", "1000000.", "998400.",
     "49600.", "50000.", "50000.", "50000.", "50000.", "49200.", "48400.", 
     "14800.", "48400.", "49800.", "49700.", "50000.", "47600.",
              "0.","50000.","50000.","50000.","50000.","49200.","49800.",
              "50000.","50000.","49600.","50000.","50000.","50000.",
    "49800.","50000.","50000."
     ];

mass = [
    "650", "1000", "650", "1000",
    "600", "800", "1000", "1200", "1400", "1800", "2000",
        "2500", "3000", "3500", "4000", "4500",
     "0","0", "0", "0", "0", "0", "0", "0", "0", "0", 
    "0", "0", "0", "0", "0", "0", "0", "0",
     "0", "0", "0", "0", "0", "0", "0", 
    "0", "0", "0", "0", "0", "0", "0", "0",
    "0", "0", "0", "0", "0", "0", "0", 
    "0", "0", "0", "0",
#        "0", "0", "0", "0", "0", "0", "0",
#        "0", "0", "0", "0",
        "600", "800", "1000", "1200", "1400", "1600", "1800",
        "2000", "2500", "3000", "3500", "4000", "4500",
        "600", "800", "1000", "1200", "1400", "1600", "1800",
        "2000","2500", "3000", "3500", "4000", "4500",
    "700","750","900"
        ];

#nameDataMu = ['data_mu_prompt_v4_25ns_runD_4']


nameDataMu = [
    "data_mu_05Oct_25ns_runD_1",
    "data_mu_05Oct_25ns_runD_2",
    "data_mu_05Oct_25ns_runD_3",
    "data_mu_05Oct_25ns_runD_4",
    "data_mu_05Oct_25ns_runD_5",
    "data_mu_05Oct_25ns_runD_6",
    "data_mu_05Oct_25ns_runD_7",
    "data_mu_05Oct_25ns_runD_8",
    "data_mu_05Oct_25ns_runD_9",
    "data_mu_05Oct_25ns_runD_10",
    "data_mu_05Oct_25ns_runD_11",
    "data_mu_05Oct_25ns_runD_12",
    "data_mu_05Oct_25ns_runD_13",
    "data_mu_05Oct_25ns_runD_14",
    "data_mu_05Oct_25ns_runD_15",
    "data_mu_05Oct_25ns_runD_16",
    "data_mu_05Oct_25ns_runD_17",
    "data_mu_prompt_v4_25ns_runD_1",
    "data_mu_prompt_v4_25ns_runD_2",
    "data_mu_prompt_v4_25ns_runD_3",
    "data_mu_prompt_v4_25ns_runD_4",
    "data_mu_prompt_v4_25ns_runD_5",
    "data_mu_prompt_v4_25ns_runD_6",
    "data_mu_prompt_v4_25ns_runD_7",
    "data_mu_prompt_v4_25ns_runD_8",
    "data_mu_prompt_v4_25ns_runD_9",
    "data_mu_prompt_v4_25ns_runD_10",
    "data_mu_prompt_v4_25ns_runD_11",
    "data_mu_prompt_v4_25ns_runD_12",
    "data_mu_prompt_v4_25ns_runD_13",
    "data_mu_prompt_v4_25ns_runD_14",
    "data_mu_prompt_v4_25ns_runD_15",
    "data_mu_prompt_v4_25ns_runD_16",
    "data_mu_prompt_v4_25ns_runD_17",
    "data_mu_prompt_v4_25ns_runD_18",
    "data_mu_prompt_v4_25ns_runD_19", 
    "data_mu_prompt_v4_25ns_runD_20",
    "data_mu_prompt_v4_25ns_runD_21",
    "data_mu_prompt_25ns_runC"];

nameDataEl = [
    "data_el_05Oct_25ns_runD_1",
    "data_el_05Oct_25ns_runD_2",
    "data_el_05Oct_25ns_runD_3",
    "data_el_05Oct_25ns_runD_4",
    "data_el_05Oct_25ns_runD_5",
    "data_el_05Oct_25ns_runD_6",
    "data_el_05Oct_25ns_runD_7",
    "data_el_05Oct_25ns_runD_8",
    "data_el_05Oct_25ns_runD_9",
    "data_el_05Oct_25ns_runD_10",
    "data_el_05Oct_25ns_runD_11",
    "data_el_05Oct_25ns_runD_12",
    "data_el_05Oct_25ns_runD_13",
    "data_el_05Oct_25ns_runD_14",
    "data_el_05Oct_25ns_runD_15",
    "data_el_05Oct_25ns_runD_16",
    "data_el_05Oct_25ns_runD_17",
    "data_el_prompt_v4_25ns_runD_1",
    "data_el_prompt_v4_25ns_runD_2",
    "data_el_prompt_v4_25ns_runD_3",
    "data_el_prompt_v4_25ns_runD_4",
    "data_el_prompt_v4_25ns_runD_5",
    "data_el_prompt_v4_25ns_runD_6",
    "data_el_prompt_v4_25ns_runD_7",
    "data_el_prompt_v4_25ns_runD_8",
    "data_el_prompt_v4_25ns_runD_9",
    "data_el_prompt_v4_25ns_runD_10",
    "data_el_prompt_v4_25ns_runD_11",
    "data_el_prompt_v4_25ns_runD_12",
    "data_el_prompt_v4_25ns_runD_13",
    "data_el_prompt_v4_25ns_runD_14",
    "data_el_prompt_v4_25ns_runD_15",
    "data_el_prompt_v4_25ns_runD_16",
    "data_el_prompt_v4_25ns_runD_17",
    "data_el_prompt_v4_25ns_runD_18",
    "data_el_prompt_v4_25ns_runD_19", 
    "data_el_prompt_v4_25ns_runD_20",
    "data_el_prompt_v4_25ns_runD_21",
    "data_el_prompt_25ns_runC"];

#MC

for a in range(len(category)):
    for i in range(len(name)):
        fn = "Job/Job_"+name[i]+"_"+category[a];
        outScript = open(fn+".sh","w");
        command = "python python/produceWWNtuples.py -n "+name[i]+" -o WWTree_"+name[i]+" -l "+category[a]+" -w "+xSecWeight[i]+" -no "+N[i]+" -mass "+mass[i]+" --ismc 1 -trig 1";
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

