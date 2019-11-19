#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A program to rename chains in biological assemblies
Each model (subunit) will have different chain names
Only PDB-formatted biological assemblies
"""
import itertools, sys, os, re, collections, subprocess, glob, optparse

chains=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','0','1','2','3','4','5','6','7','8','9']

def create_modified_biological_assemblies_pdbid(pdbid,outputfile):
	#pdbid = "4ex1"
	#outputfile = "4ex1_modified.pdb"
	records = "https://files.rcsb.org/download/{}.pdb1".format(pdbid)
	proc = subprocess.Popen("wget {}".format(records), shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
	out = proc.communicate()
	if pdbid+".pdb1" in glob.glob("*"):
		ba_file = pdbid+".pdb1"
		pdb_lines=open(ba_file,"r").read().split("MODEL       ")
		print "This structure contain {} models..".format(len(pdb_lines[1:]))
		calculated_chains=[k for n in [sorted(list(set([model.split("\n")[0].strip()+line[21:22] for line in model.split("\n") if line[:4]=="ATOM"]))) for model in pdb_lines[1:]] for k in n]
		if len(list(set(calculated_chains)))<63:
			atom_lines = []
			output=open(outputfile,"w")
			#header
			header_lines = pdb_lines[0][:-1]
			#enumerate all models
			for index, model in enumerate(pdb_lines[1:]):
				model_chains=sorted(list(set([line[21:22] for line in model.split("\n") if line[:4]=="ATOM"]+[line[21:22] for line in model.split("\n") if line[:6]=="HETATM"])))
				new_chains = chains[:len(model_chains)]
				chain_dictionary = dict(zip(model_chains, new_chains))
				model_lines_atom = ["MODEL       {}".format(str(index+1))]+[line[:21]+chain_dictionary[line[21:22]][0]+line[22:] for line in model.split("\n") if line[:4]=="ATOM"]+[line[:21]+chain_dictionary[line[21:22]][0]+line[22:] for line in model.split("\n") if line[:6]=="HETATM"]+["ENDMDL"]
				for a in model_lines_atom: atom_lines.append(a)
				for s in new_chains: chains.remove(s)
			finalized_lines = header_lines+"\n"+"\n".join(atom_lines)
			print>>output, finalized_lines
		else:
			print "Large number of chains detected. Cannot rename chains.."
	else:
		print "Cannot download pdb file containing biological assemblies. Please check manually.."

def create_modified_biological_assemblies_file(inputfile,outputfile):
	pdb_lines=open(inputfile,"r").read().split("MODEL       ")
	print "This structure contain {} models..".format(len(pdb_lines[1:]))
	calculated_chains=[k for n in [sorted(list(set([model.split("\n")[0].strip()+line[21:22] for line in model.split("\n") if line[:4]=="ATOM"]))) for model in pdb_lines[1:]] for k in n]
	if len(list(set(calculated_chains)))<63:
		atom_lines = []
		output=open(outputfile,"w")
		#header
		header_lines = pdb_lines[0][:-1]
		#enumerate all models
		for index, model in enumerate(pdb_lines[1:]):
			model_chains=sorted(list(set([line[21:22] for line in model.split("\n") if line[:4]=="ATOM"]+[line[21:22] for line in model.split("\n") if line[:6]=="HETATM"])))
			new_chains = chains[:len(model_chains)]
			chain_dictionary = dict(zip(model_chains, new_chains))
			model_lines_atom = ["MODEL       {}".format(str(index+1))]+[line[:21]+chain_dictionary[line[21:22]][0]+line[22:] for line in model.split("\n") if line[:4]=="ATOM"]+[line[:21]+chain_dictionary[line[21:22]][0]+line[22:] for line in model.split("\n") if line[:6]=="HETATM"]+["ENDMDL"]
			for a in model_lines_atom: atom_lines.append(a)
			for s in new_chains: chains.remove(s)
		finalized_lines = header_lines+"\n"+"\n".join(atom_lines)
		print>>output, finalized_lines
	else:
		print "Large number of chains detected. Cannot rename chains.."


#13-command line arguments 
if __name__ == '__main__': # 13- PARSE COMMAND LINE ARGUMENTS
	p1="\n USAGE: \n\n python rename_chains.py -i pdbid or inputfile -o outputfile\n"
	p2="\n e.g. python rename_chains.py -i 4ex1 -o 4ex1_modified.pdb\n"
	p3="\n e.g. python rename_chains.py -i 4ex1.pdb -o 4ex1_modified.pdb\n"
	usage=p1+p2+p3
	#usage = "\n USAGE: \n\n python auto_imaaagine.py -q search -r residuenumber -n residuelist \n OR \n python auto_imaaagine.py -q search -r residuenumber -p exactpattern \n OR \n python auto_imaaagine.py -q analyse -o output_dir -f filename \n OR \n python auto_imaaagine.py -q analyse -o output_dir -f directory \n"
	parser = optparse.OptionParser(usage=usage)
	parser.add_option("-i", "--input", dest="input", type=str, help="pdbid or inputfile")
	parser.add_option("-o", "--output", dest="output", type=str, help="outputfile")
	(options, args) = parser.parse_args()

	if (options.input and options.output):
		if ".pdb" in options.input:
			create_modified_biological_assemblies_file(options.input,options.output)
		else:
			create_modified_biological_assemblies_pdbid(options.input,options.output)
	else:
		print "Not enough input arguments supplied"
		print usage
		quit()
