#!/usr/bin/env python
# -*- coding: utf-8 -*-
#Author: Nur Syatila Ab Ghani
#A program to create ligand binding sites
#Input: inputpdb, ligand, distance_cutoff

import os
import sys
import glob
import re
import collections
import itertools
import argparse

trans = {'ALA':'A','CYS':'C','CYH':'C','CSS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G','HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N','PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V','TRP':'W','TYR':'Y','UNK':'X','MSE':'M'}

atom_list = {
            "ALA":('N  ','C  ','O  ','CA ','CB '), "CYS":('N  ','C  ','O  ','CA ','CB ','SG '),"ASP":('N  ','C  ','O  ','CA ','CB ','CG ','OD1','OD2'),"GLU":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','OE1','OE2'),"PHE":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2','CE1','CE2','CZ '),"GLY":('N  ','C  ','O  ','CA '),"HIS":('N  ','C  ','O  ','CA ','CB ','CG ','ND1','CD2','CE1','NE2'),
            "ILE":('N  ','C  ','O  ','CA ','CB ','CG1','CG2','CD1'),"LYS":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','CE ','NZ '),"LEU":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2'),"MET":('N  ','C  ','O  ','CA ','CB ','CG ','SD ','CE '),"ASN":('N  ','C  ','O  ','CA ','CB ','CG ','OD1','ND2'),"PRO":('N  ','C  ','O  ','CA ','CB ','CG ','CD '),"GLN":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','OE1','NE2'),"ARG":('N  ','C  ','O  ','CA ','CB ','CG ','CD ','NE ','CZ ','NH1','NH2'),
            "SER":('N  ','C  ','O  ','CA ','CB ','OG '),"THR":('N  ','C  ','O  ','CA ','CB ','OG1','CG2'),"VAL":('N  ','C  ','O  ','CA ','CB ','CG1','CG2'),"TRP":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2'),"TYR":('N  ','C  ','O  ','CA ','CB ','CG ','CD1','CD2','CE1','CE2','CZ ','OH ')}
sc_list = {
            "ALA":('CB '), "CYS":('CB ','SG '),"ASP":('CB ','CG ','OD1','OD2'),"GLU":('CB ','CG ','CD ','OE1','OE2'),"PHE":('CB ','CG ','CD1','CD2','CE1','CE2','CZ '),"GLY":('CA '),"HIS":('CB ','CG ','ND1','CD2','CE1','NE2'),
            "ILE":('CB ','CG1','CG2','CD1'),"LYS":('CB ','CG ','CD ','CE ','NZ '),"LEU":('CB ','CG ','CD1','CD2'),"MET":('CB ','CG ','SD ','CE '),"ASN":('CB ','CG ','OD1','ND2'),"PRO":('CB ','CG ','CD '),"GLN":('CB ','CG ','CD ','OE1','NE2'),"ARG":('CB ','CG ','CD ','NE ','CZ ','NH1','NH2'),
            "SER":('CB ','OG '),"THR":('CB ','OG1','CG2'),"VAL":('CB ','CG1','CG2'),"TRP":('CB ','CG ','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2'),"TYR":('CB ','CG ','CD1','CD2','CE1','CE2','CZ ','OH ')}

#calculation of distance between coordinates of x,y,z
def get_distance(a,b):
    dist = ((a[1]-b[1])**2 + (a[2]-b[2])**2 + (a[3]-b[3])**2)**0.5
    return dist
    
#calculation of distance between coordinates of x,y,z
def get_distance2(a,b):
    dist = ((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)**0.5
    return dist
    
#get the center of coordinates of x,y,z
def get_centroid_coord(all_coords):
    x_coord = sum([x[0] for x in all_coords])/len(all_coords)
    y_coord = sum([y[1] for y in all_coords])/len(all_coords)
    z_coord = sum([z[2] for z in all_coords])/len(all_coords)
    centroid = [x_coord]+[y_coord]+[z_coord]
    return centroid

#get binding sites for given ligand
def get_binding_sites_pdb(pdbfile,ligid,distance_cutoff):
    print ("creating directory..")
    if not os.path.exists("binding_sites"): os.mkdir("binding_sites")
    
    print ("processing PDB structure..")
    pdb_read = open(pdbfile, 'r').read()
    if 'MODEL        1' in pdb_read:
        print ("processing nmr models - ",pdbfile)
        pdb_lines = pdb_read.split('MODEL        2')[0].split('\n')
        f_out = open(pdbfile, 'w')
        for l in pdb_lines: f_out.write(l+'\n')
    pdb = pdbfile.split("/")[-1].replace(".pdb","")
    
    print ("processing list of ligands..")
    ligandid = " "*(3-len(ligid))+ligid
    matched_line = sorted(list(set([line.replace("\n","")[17:20] for line in open(pdbfile,"r").readlines() if line[:6] == "HETATM" if line[17:20] ==ligandid])))
    if matched_line !=[]:
        print ("generate binding sites..")
        atom_coord = []
        het_coord = []
        atom_lines = [line.replace("\n","") for line in open(pdbfile,"r").readlines() if line[:4] == "ATOM"]
        for line in atom_lines:
            for sc in sc_list.keys():
                if line[17:20] == sc:
                    atom_coord.append([line[13:26], float(line[30:38]), float(line[38:46]),float(line[46:54])])
        het_lines = [line.replace("\n","") for line in open(pdbfile,"r").readlines() if line[:6] == "HETATM" if line[17:20]==ligandid]
        for line in het_lines:
            het_coord.append([line[13:26], float(line[30:38]), float(line[38:46]),float(line[46:54])])
        E = collections.defaultdict(list)
        for h in het_coord: E[h[0][4:]].append(h)
        het_records = [list(i) for i in E.items()]
        bs_residues = []
        for het in het_records:
            combinations_aa_het = [i for i in itertools.product(atom_coord,het[1])]
            ligand_binding_residues = []
            for l in combinations_aa_het:
                if get_distance(l[0],l[1])<float(distance_cutoff):
                    ligand_binding_residues.append(l[0][0][4:])
            ligand_binding_residues2 = sorted(list(collections.OrderedDict.fromkeys(list(set(ligand_binding_residues)))),key=lambda x:int(x[5:]))
            if len(ligand_binding_residues2)>2: bs_residues.append([het[0]]+ligand_binding_residues2)
        for per_set in bs_residues:
            sc_coord2 = []
            Ksc2 = []
            scfound = set()
            bs_file_output = "binding_sites/"+pdb+"_"+per_set[0].replace(" ","")+"_"+str(len(per_set[1:]))+".pdb"
            with open(bs_file_output,"w") as output:
                atoms = [atom_list[u[:3]] for u in per_set[1:]]
                V = [i for k in atoms for i in k]
                with open(pdbfile,"r") as pdb_input2:
                    pdb_lines = pdb_input2.readlines()
                    t = [line.replace("\n","") for line in pdb_lines if line[:4] =="ATOM" and any(item in line[17:26] for item in per_set[1:]) and any(thing in line[13:16] for thing in V)]
                    for i in t: sc_coord2.append(i)
                    s = [line.replace("\n","") for line in pdb_lines if line[:6] == "HETATM" and line[17:26] == per_set[0]]
                    for r in s: sc_coord2.append(r)
                Ksc1 = [i[:16]+i[16:17].replace(i[16:17]," ")+i[17:] for i in sc_coord2]
                for item in Ksc1:
                    n = item[13:26]
                    if n not in scfound:
                        Ksc2.append(item)
                        scfound.add(n)
                for i in Ksc2:
                    output.write(i.replace("\n","")+"\n")
        print (pdb,"- completed")
    else: print (pdb, "No ligand (from the list) found.")

# Create the parser
my_parser = argparse.ArgumentParser(description='Transform a PDB-formatted structure')

# Add the arguments
my_parser.add_argument('pdbfile',type=str,help='Input PDB')
my_parser.add_argument('ligid',type=str,help='Ligand ID (2 or 3-letter code')
my_parser.add_argument('distance_cutoff',type=str,help='Distance cut-off for the distance between residue and ligand')

# Execute the parse_args() method
args = my_parser.parse_args()
print ("Transform a PDB-formatted structure")

get_binding_sites_pdb(args.pdbfile,args.ligid,float(args.distance_cutoff))

#distance_cutoff = 4.0
#pdbfile = "1a4g.pdb"
#ligid = "ZMR"
#get_binding_sites_pdb(pdbfile,ligid,distance_cutoff)

