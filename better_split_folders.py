#! /usr/bin/env python                                           
import os
import glob
import math
from array import array
import sys
import time
import subprocess

channel=["mu","el"]

print "qui"

for l in range(len(channel)):
    for i in range(1,39):
        command = "mkdir data_"+channel[l]+"_2016_runB_v2_"+str(i)
        print command
        os.system(command);

    c3=1
    c2=1
    for c1 in range(1,7464):
        command2 = "mv data_"+channel[l]+"_2016_runB_v2/ReducedSelection_"+str(c1)+".root data_"+channel[l]+"_2016_runB_v2_"+str(c3)
        print command2
        os.system(command2)
        c2=c2+1    
        if (c2==200): 
            c3=c3+1
            c2=1


for l in range(len(channel)):
    for i in range(1,14):
        command = "mkdir data_"+channel[l]+"_2016_runC_v2_"+str(i)
        print command
        os.system(command);

    c3=1
    c2=1
    for c1 in range(1,2484):
        command2 = "mv data_"+channel[l]+"_2016_runC_v2/ReducedSelection_"+str(c1)+".root data_"+channel[l]+"_2016_runC_v2_"+str(c3)
        print command2
        os.system(command2)
        c2=c2+1    
        if (c2==200): 
            c3=c3+1
            c2=1



for l in range(len(channel)):
    for i in range(1,22):
        command = "mkdir data_"+channel[l]+"_2016_runD_v2_"+str(i)
        print command
        os.system(command);

    c3=1
    c2=1
    for c1 in range(1,4167):
        command2 = "mv data_"+channel[l]+"_2016_runD_v2/ReducedSelection_"+str(c1)+".root data_"+channel[l]+"_2016_runD_v2_"+str(c3)
        print command2
        os.system(command2)
        c2=c2+1    
        if (c2==200): 
            c3=c3+1
            c2=1

#mv data_el_2016_runB_v2/ReducedSelection_74* data_el_2016_runB_v2_37/

#mv data_el_2016_runB_v2/ReducedSelection_* data_el_2016_runB_v2_1/


