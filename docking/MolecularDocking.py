#!/usr/bin/env python
import sys, os, chimera, argparse
from chimera import runCommand
from AddCharge import initiateAddCharges
initiateAddCharges(method='am1-bcc', nogui=True)
from chimera.tkgui import saveReplyLog

# Perform Docking Against the Whole Protein Structure (Blind Docking)
# Example 1 (id): /Applications/Chimera.app/Contents/MacOS/chimera --script "MolecularDocking.py --dock blind --t 4ey7 --c A --l 44159808 --v /usr/local/bin/vina"
# Example 2 (file): /Applications/Chimera.app/Contents/MacOS/chimera --script "MolecularDocking.py --dock blind --t target.cif --c A --l ligand.sdf --v /usr/local/bin/vina"

def BlindDocking(target,chainselection,ligand,vinapath):
    #target
    if chainselection=='all': runCommand('open {}; del #0.2-100; sel protein; sel invert sel; del sel; del solvent; del ligand; del ions; wait'.format(target))
    else: runCommand('open {}; del #0.2-100; sel #0:.{} & protein; sel invert sel; del sel; del solvent; del ligand; del ions; wait'.format(target,chainselection))
    #runCommand('~sel; swapaa same #0 preserve true ignoreOtherModels true; wait')
    runCommand('addh spec #0; addcharge all spec #0 method am1; wait')
    #ligand
    if os.path.exists(ligand): runCommand('open {}; wait'.format(ligand))
    else: runCommand('open pubchem:{}; wait'.format(ligand))
    #runCommand('swapaa same #1 preserve true ignoreOtherModels true; wait')
    runCommand('addh spec #1; addcharge all spec #1 method am1; wait')
    #docking
    runCommand('vina docking receptor #0 ligand #1 output dock prep false backend local location {} wait true; wait'.format(vinapath))
    #save session
    runCommand('save DockSession.py; wait')
    runCommand('combine #0,2.1 refSpec #0 newchainids true; wait; write format pdb 3 complex.pdb; wait')
    saveReplyLog('DockSessionLog.txt')
    runCommand('stop')

#Perform Molecular Docking against the Binding Site (use binding site as the reference / grid box)
# Example 1 (id,reference ligand for selection of bs): /Applications/Chimera.app/Contents/MacOS/chimera --script "MolecularDocking.py --dock bsite --t 4ey7 --c A --r E20 --l 44159808 --v /usr/local/bin/vina"
# Example 2 (file, reference residues for selection of bs): /Applications/Chimera.app/Contents/MacOS/chimera --script "MolecularDocking.py --dock bsite --t target.cif --c A --r 124.a,297.a,202.a,86.a,286.a,337.a,341.a --l ligand.sdf --v /usr/local/bin/vina"

def BSDocking(target,chainselection,bsres,ligand,vinapath):
    #target, bsres
    if chainselection=='all':
        runCommand('open {}; measure center :{}; del #0.2-100; sel protein; sel invert sel; del sel; del solvent; del ligand; del ions; wait'.format(target,bsres))
    else:
        runCommand('open {}; measure center :{}.{}; del #0.2-100; sel #0:.{} & protein; sel invert sel; del sel; del solvent; del ligand; del ions; wait'.format(target,bsres,chainselection,chainselection))
    saveReplyLog('COM.txt')
    runCommand('wait;')
    #get gridbox
    if os.path.exists('COM.txt'):
        com = [r for r in open('COM.txt','r').readlines() if 'coordinate system' in r][0]
        com2 = ','.join(com.replace('\n','').split('(')[-1].replace(')','').split(', '))
        print (com2)
    #runCommand('~sel; swapaa same #0 preserve true ignoreOtherModels true; wait')
    runCommand('addh spec #0; addcharge all spec #0 method am1; wait')
    #ligand
    if os.path.exists(ligand): runCommand('open {}; wait'.format(ligand))
    else: runCommand('open pubchem:{}; wait'.format(ligand))
    #runCommand('swapaa same #1 preserve true ignoreOtherModels true; wait')
    runCommand('addh spec #1; addcharge all spec #1 method am1; wait')
    #docking
    #runCommand('vina docking receptor #0 ligand #1 output dock prep false backend local location {} wait true; wait'.format(vinapath))
    runCommand('vina docking receptor #0 ligand #1 output dock prep false search_center {} search_size 30,30,30 backend local location {} wait true; wait'.format(com2, vinapath))
    #save session
    runCommand('save DockSession.py; wait')
    runCommand('combine #0,2.1 refSpec #0 newchainids true; wait; write format pdb 3 complex.pdb; wait')
    saveReplyLog('DockSessionLog.txt')
    runCommand('stop')
    
# Create the parser
my_parser = argparse.ArgumentParser(description='Run a molecular docking given a protein target (file or pdbid) and a ligand molecule (file or pubchem id')

# Add the arguments
my_parser.add_argument('--dock',type=str,help='Type of docking (blind/bsite)')
my_parser.add_argument('--t',type=str,help='Input file (.cif, .pdb) or PDB ID of protein target')
my_parser.add_argument('--c',type=str,help='Chain of protein target, e.g. "A" for chain A only, or "all" for all chain')
my_parser.add_argument('--l',type=str,help='Input file (.pdb, .sdf) or PubChem ID of ligand molecule')
my_parser.add_argument('--v',type=str,help='Path to executable vina')
my_parser.add_argument('--r',type=str,help='Residue for the selection of binding site')


# Execute the parse_args() method
args = my_parser.parse_args()
print ("Run a molecular docking for protein-ligand interaction")
docktype = args.dock
target = args.t
chainselection = args.c
ligand = args.l
vinapath = args.v
bres = args.r

if docktype == 'blind':
    BlindDocking(target, chainselection, ligand, vinapath)
if docktype == 'bsite':
    BSDocking(target, chainselection, bres, ligand, vinapath)

#target = 'target.cif'
#chainselection = 'A'
#ligand = 'ligand.sdf'
#vinapath = '/usr/local/bin/vina'
#BlindDocking(target,chainselection,ligand,vinapath)
