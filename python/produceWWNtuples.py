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
    parser.add_argument ('-i', '--input' , default = '../ReducedSelection.root' , help='input file')
    parser.add_argument ('-o', '--output' , default = 'output', help='output file')
    parser.add_argument ('-mc', '--ismc' , default = 'False', help='is MC or not')
    args = parser.parse_args ()

    print 'produceWWNtuples '+args.input+' '+args.output+' '+args.ismc
    os.system('./produceWWNtuples.exe '+args.input+' '+args.output+' '+args.ismc)
