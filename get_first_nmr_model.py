#!/usr/bin/env python
# -*- coding: utf-8 -*-
#A program to get the first NMR model (Model 1) from NMR structure
#Input: inputpdb, outputpdb

import sys, os, argparse

def get_first_nmr_model(inputpdb,outputpdb):
    pdb_read = open(inputpdb, 'r').read()
    
    #check if the structure is an NMR structure
    if 'MODEL        1' in pdb_read:
        print (inputpdb," - NMR structure, processing..")
        pdb_lines = pdb_read.split('MODEL        2')[0].split('\n')
        with open(outputpdb, 'w') as pdbout:
            for l in pdb_lines:
                pdbout.write(l+'\n')
    else:
        print (inputpdb, " - Not an NMR structure.")

# Create the parser
my_parser = argparse.ArgumentParser(description='Get first model from NMR structure')

# Add the arguments
my_parser.add_argument('inputpdb',type=str,help='Input PDB')
my_parser.add_argument('outputpdb',type=str,help='Output PDB')

# Execute the parse_args() method
args = my_parser.parse_args()
print ("Get first model from NMR structure")

get_first_nmr_model(args.inputpdb,args.outputpdb)

#e.g. python3 get_first_nmr_model.py 1a0s.pdb 1a0s.pdb


