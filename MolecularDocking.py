#!/usr/bin/env python
import sys, os, chimera, argparse
from chimera import runCommand
from AddCharge import initiateAddCharges
initiateAddCharges(method='am1-bcc', nogui=True)
from chimera.tkgui import saveReplyLog

# Example 1: /Applications/Chimera.app/Contents/MacOS/chimera --script "MolecularDocking.py 4ey7 A 44159808 /usr/local/bin/vina"
# Example 2: /Applications/Chimera.app/Contents/MacOS/chimera --script "MolecularDocking.py target.cif A ligand.sdf /usr/local/bin/vina"

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
    

# Create the parser
my_parser = argparse.ArgumentParser(description='Run a blind molecular docking given a protein target (file or pdbid) and a ligand molecule (file or pubchem id')

# Add the arguments
my_parser.add_argument('target',type=str,help='Input file (.cif, .pdb) or PDB ID of protein target')
my_parser.add_argument('chainselection',type=str,help='Chain of protein target, e.g. "A" for chain A only, or "all" for all chain')
my_parser.add_argument('ligand',type=str,help='Input file (.pdb, .sdf) or PubChem ID of ligand molecule')
my_parser.add_argument('vinapath',type=str,help='Path to executable vina')

# Execute the parse_args() method
args = my_parser.parse_args()
print ("Run a blind molecular docking for protein-ligand interaction")
target = args.target
chainselection = args.chainselection
ligand = args.ligand
vinapath = args.vinapath

BlindDocking(target,chainselection,ligand,vinapath)


#target = 'target.cif'
#chainselection = 'A'
#ligand = 'ligand.sdf'
#vinapath = '/usr/local/bin/vina'
#BlindDocking(target,chainselection,ligand,vinapath)
