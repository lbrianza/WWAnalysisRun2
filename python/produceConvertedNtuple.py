#!/usr/bin/python
import sys
import os
import commands
from commands import getstatusoutput
import datetime
import argparse
import string
if __name__ == '__main__':
    parser = argparse.ArgumentParser (description = 'produce ntuples with WW semileptonic final state')
    parser.add_argument ('-i', '--inputFolder' , default = 'pseudoData/' , help='input folder with the reduced trees')
    parser.add_argument ('-o', '--output' , default = 'Pseudodata.root', help='output file')
    parser.add_argument ('-l', '--lepton' , default = 'mu', help='lepton category (mu or el)')
    parser.add_argument ('-t', '--tree' , default = 'PKUTree', help='name of the input tree')
    parser.add_argument ('-n', '--name' , default = 'mu_PKUTree_pdata.root' , help='input file')
    args = parser.parse_args ()

    print 'ntupleConverter '+args.inputFolder+' '+args.output+' '+args.lepton+' '+args.tree+' '+args.name
    os.system('./ntupleConverter.exe '+args.inputFolder+' '+args.output+' '+args.lepton+' '+args.tree+' '+args.name)
