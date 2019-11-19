#!/usr/bin/env python
"""
Script to select only teh first model in nmr structure
!! Require Biopython.
PDB formatted file only
"""
import sys
import os
from Bio.PDB import *

def get_first_model(infile,outfile):
	pdbid = infile.split("/")[-1][:4]
	parser = PDBParser()
	structure = parser.get_structure(pdbid, infile)
	model = structure[0]
	#header = parser.get_header()
	io = PDBIO()
	io.set_structure(model)
	io.save(outfile)

"""
Command line arguments
"""
if len(sys.argv)<2: 
	print "not enough arguments.\n e.g. select_first_model_nmr.py infile outfile"
else:
	get_first_model(sys.argv[1],sys.argv[2])


#try:
#open file
#infile = "1g03.pdb"
#outfile = "1g03_nmr.pdb"






