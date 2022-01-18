#!/usr/bin/env python
# -*- coding: utf-8 -*-
#A program to transform a PDB-formatted structure (pattern or whole PDB), given the matrix (rotation, translation)
#Requires Bio module (Biopython)
#Input: inputpdb, matrixlist, outputpdb

import itertools, sys, os, re, operator, argparse, csv
from Bio import PDB
import numpy

def transform_pdb(inputpdb,matrixlist,outputpdb)
    #open PDB structure
    parser = PDB.PDBParser()
    pdbname = inputpdb.split('.')[0]
    structure = parser.get_structure(pdbname, inputpdb)
    
    #get matrix (rotation, translation) from the list
    matrix_dict = {
        'rotation': [[r[0]]+[r[1]]+[r[2]], [r[3]]+[r[4]]+[r[5]], [r[6]]+[r[7]]+[r[8]]],
        'translation': [r[9]]+[r[10]]+[r[11]]
    }
    rotation = numpy.array([[float(s) for s in k] for k in matrix_dict['rotation']])
    translation = [float(k) for k in matrix_dict['translation']]
    
    #perform transformation on each atom
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom.transform(rotation, translation)
                
    #save the transformed structure
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(outputpdb)

# Create the parser
my_parser = argparse.ArgumentParser(description='Transform a PDB-formatted structure')

# Add the arguments
my_parser.add_argument('inputpdb',type=str,help='Input PDB')
my_parser.add_argument('matrixlist',type=str,help='Matrix in a list (3x3 rotation and translation) (joined by comma)')
my_parser.add_argument('outputpdb',type=str,help='Output PDB')

# Execute the parse_args() method
args = my_parser.parse_args()
print ("Transform a PDB-formatted structure")

transform_pdb(args.inputpdb,args.matrixlist,args.outputpdb)

#e.g. transform_pdb("1a0s.pdb","0.3786,-0.7283,-0.5712,0.8695,0.4913,-0.0500,0.3170,-0.4777,0.8193,93.254,19.092,37.797","1a0s_transformed.pdb")


