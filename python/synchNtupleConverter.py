#!/usr/bin/python
import sys
import os
import commands
from commands import getstatusoutput
import datetime
import argparse
import string
if __name__ == '__main__':
    parser = argparse.ArgumentParser (description = 'produce synch ntuples with WW semileptonic final state')
    parser.add_argument ('-i', '--inputFolder' , default = './output/output_' , help='input folder with the reduced trees')
    parser.add_argument ('-o', '--output' , default = 'RSGraviton1000', help='output file')
    parser.add_argument ('-l', '--lepton' , default = 'mu', help='lepton category (mu or el)')
    parser.add_argument ('-n', '--name' , default = 'RSGraviton1000' , help='input file')
    args = parser.parse_args ()

    print 'synchNtupleConverter.exe '+args.inputFolder+' '+args.output+' '+args.lepton+' '+args.name
    os.system('./synchNtupleConverter.exe '+args.inputFolder+' '+args.output+' '+args.lepton+' '+args.name)
