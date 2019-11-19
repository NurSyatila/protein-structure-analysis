#!/usr/bin/env python
"""
Script to perform dalilite search and obtain the matrix file using Python.
User have to obtain licence for DaliLite program and install it.
Make sure the DAT directory is available (read the documentation).
Linux system available.
"""
import sys
import os
import re
import glob
import shutil
import collections
import subprocess
import itertools

"""
Part 1: Convert PDB file to DAT file for DaliLite search.
Dalilite search DAT file and compare against other DAT files.
To make it easier for a large scale search, create a directory containing all DAT files matching to the database of PDB.
DAT file will be stored in DAT directory.
"""
def dalilite_readdat(pdbdir):
	list_pdbs = [fi for fi in glob.glob("{}/*.pdb".format(pdbdir))]
	for pdbfile in list_pdbs:
		proc = subprocess.Popen('DaliLite -readbrk {} {}'.format(pdbfile, pdbfile[:4]),shell=True,stdout=subprocess.PIPE) #create DAt file for each pdbchain
		out=proc.communicate()
		if "dali.default" in [fi for fi in glob.glob("*")]: os.remove("dali.default") #if error, delete error files and continue for next session
		if "dali.lock" in [fi for fi in glob.glob("*")]: os.remove("dali.lock")
		#os.remove("{}.dssp".format(pdbfile[:4])) #dssp file will be generated (optional to remove file)

"""
Part 2: Pairwise comparison between two pdbchain.
Make sure the PDB file had been converted into DAT file.
"""
def dalilite_pairwise(pdbone,pdbtwo):
	proc = subprocess.Popen('DaliLite -align {} {}'.format(pdbone,pdbtwo),shell=True,stdout=subprocess.PIPE) #perform pairwise fold comparison
	out=proc.communicate()
	if "dali.default" in [fi for fi in glob.glob("*")]: os.remove("dali.default")
	if "dali.lock" in [fi for fi in glob.glob("*")]: os.remove("dali.lock")
	combinations = ["{}\t{}".format(pdbone,pdbtwo)]
	parse_dccp(combinations)
	
	
"""
Part 3: One-against-all comparison
The list of pdbchains for query should be stored in a text document.
"""
def dalilite_one_all(pdbone,textfile):
	proc = subprocess.Popen('DaliLite -list {} {}'.format(pdbone,textfile),shell=True,stdout=subprocess.PIPE) #perform pairwise fold comparison
	out=proc.communicate()
	if "dali.default" in [fi for fi in glob.glob("*")]: os.remove("dali.default")
	if "dali.lock" in [fi for fi in glob.glob("*")]: os.remove("dali.lock")
	chains = [chain.replace("\n","") for chain in open(textfile,"r").readlines()]
	combinations = [pdbone+"\t"+chain for chain in chains]
	parse_dccp(combinations)

"""
Part 4: One-against-all comparison
The list of pdbchains for query should be stored in a text document.
"""

def dalilite_all_all(textfile):
	proc = subprocess.Popen('DaliLite -AllAll {}'.format(textfile),shell=True,stdout=subprocess.PIPE) #perform pairwise fold comparison
	out=proc.communicate()
	if "dali.default" in [fi for fi in glob.glob("*")]: os.remove("dali.default")
	if "dali.lock" in [fi for fi in glob.glob("*")]: os.remove("dali.lock")
	chains = [chain.replace("\n","") for chain in open(textfile,"r").readlines()]
	combinations=["\t".join(k) for k in itertools.combinations(chains,2)]

"""
Part 4: Parse DCCP file (output file from DaliLite search) to generate matrix of fold
"""
def parse_dccp(combinations):
	if [files for files in glob.glob("*.dccp")] !=[]:
		N = []
		output = open("fold_comparison.txt","w") #get list of pairwise fold similarity value in Z-score
		for dccp_file in glob.glob("*.dccp"):
			input = open(dccp_file,"r")
			dccp_1 = [line[69:80].split(" ")+[line[1:5]]+[int(line[5:9])]+[float(line[9:18])]+[float(line[18:22])]+[int(line[22:26])]+[float(line[26:34])]+[int(line[43:46])]+[int(line[46:53])] for line in input.readlines() if line.startswith(" DCCP")]
			dccp_2 = list(set(['\t'.join(i[:2]+[str(i[7])]+[str(i[8])]) for i in dccp_1]))
			for i in dccp_2: N.append(i)
		dc = sorted(sorted(sorted(list(set(N)),key=lambda x: x.split("\t")[3]),key=lambda x: x.split("\t")[2]),key=lambda x: sorted(x.split("\t")[:2]))
		for i in dc: N.append(i)	
		N2 = [n.split("\t") for n in sorted(list(set(N)))]
		dc_dict = collections.defaultdict(list)
		for n in N2: dc_dict["\t".join(n[:2])].append("\t".join(n[2:]))
		dc_dict2 = {k:sorted(sorted(v,key=lambda x:float(x.split("\t")[0]),reverse=True),key=lambda x:float(x.split("\t")[1]),reverse=True)[0] for k,v in dc_dict.iteritems()}
		dc_dict3 = sorted(["\t".join(list(i)) for i in dc_dict2.items()])
		print "\nDone! Fold similarity is stored in fold_comparison.txt\n"
		print>>output, "pdb1\tpdb2\tZ-score\tseqid\n"+"\n".join(dc_dict3)
	for dssp_file in glob.glob("*.dssp"): os.remove(dssp_file)
	if [files for files in glob.glob("*.dccp")] ==[]:
		print "No output generated. Search failed or no fold similarity.."
"""
Command line arguments
"""
if len(sys.argv)==1: 
	print "not enough arguments.\n e.g. convert pdb to dat:: readdat pdbdir\n e.g. pairwise comparison:: pairwise pdbone pdbtwo\n e.g. one against all comparison:: one-all pdbone textfile\n e.g. all-against-all comparison:: all-all textfile\n"
else:
	if sys.argv[1]=="readdat":
		print "Converting pdb file to dat file.."
		dalilite_readdat(sys.argv[2])
	if sys.argv[1]=="pairwise":
		print "Pairwise comparison.."
		dalilite_pairwise(sys.argv[2],sys.argv[3])
	if sys.argv[1]=="one-all":
		print "One against all comparison.."
		dalilite_one_all(sys.argv[2],sys.argv[3])	
	if sys.argv[1]=="all-all":
		print "All-against-all comparison.."
		dalilite_all_all(sys.argv[2])
