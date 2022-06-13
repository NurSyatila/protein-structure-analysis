# protein-structure-analysis
Scripts to parse PDB file and analyse protein structures and sequences

### 1) transform_pdb.py: 
  - A program to transform a PDB-formatted structure (pattern or whole PDB), given the matrix (rotation, translation)
  - Requires Bio module (Biopython)
  - Arguments: inputpdb, matrixlist, outputpdb
  - Usage: transform_pdb.py inputpdb matrixlist outputpdb
  - E.g. transform_pdb.py 1a0s.pdb 0.3786,-0.7283,-0.5712,0.8695,0.4913,-0.0500,0.3170,-0.4777,0.8193,93.254,19.092,37.797 1a0s_transformed.pdb
  
### 2) get_first_nmr_model.py
  - A program to get the first NMR model (Model 1) from NMR structure
  - Arguments: inputpdb, outputpdb
  - Usage: get_first_nmr_model.py inputpdb outputpdb
  - E.g. get_first_nmr_model.py 1g03.pdb 1g03_nmr.pdb
  
### 3) generate_binding_site.py
  - A program to create ligand binding sites
  - Arguments: pdbfile, ligid, distance_cutoff
  - Usage: generate_binding_site.py pdbfile ligid distance_cutoff
  - E.g. generate_binding_site.py 1a4g.pdb ZMR 4.0
  
### 4) CreateNRPDB.py (requires BatchDownload.sh)
  - A program to create a data set of non-redundant PDB from a list of representative PDB entities (in PDB format) 
  - Based on results of the weekly clustering of protein sequences in the PDB by MMseqs2 at 30%, 40%, 50%, 70%, 90%, 95%, or 100% sequence identity
  - Arguments: SeqId (Sequence identity : 30, 40, 50, 70, 90, 95 or 100)
  - Usage: CreateNRPDB.py 30
  - The following directories will be created: 
      pdb (original pdb files)
      pdbCAonly (pdb with CA atoms only)
      pdb<SeqID> (pdb files containing representative chains)
  
### 5) MolecularDocking.py
  - Run a molecular docking for protein-ligand interaction, given a protein target (file or pdbid) and a ligand molecule (file or pubchem id)
  - Note: 
    a) The script must be used through UCSF Chimera
    b) vina must be installed   
  - Option 1: Blind Docking:  
    - Arguments: target (file or id), chainselection (selected chain, e.g. 'A' or 'all'), ligand (file or id), vinapath. 
  - Option 2: Binding Site Docking:  
    - Arguments: target (file or id), chainselection (selected chain, e.g. 'A' or 'all'), ligand (file or id), vinapath. 
  - Usage:  
    ##### Option 1: Blind Docking. 
    ~ Perform Docking Against the Whole Protein Structure (Blind Docking).   
    ~ Example 1 (id): /Applications/Chimera.app/Contents/MacOS/chimera --script "MolecularDocking.py --dock blind --t 4ey7 --c A --l 44159808 --v /usr/local/bin/vina".  
    ~ Example 2 (file): /Applications/Chimera.app/Contents/MacOS/chimera --script "MolecularDocking.py --dock blind --t target.cif --c A --l ligand.sdf --v /usr/local/bin/vina".  
  
    ##### Option 2: Binding Site Docking.  
    ~ Perform Molecular Docking against the Binding Site (use binding site as the reference / grid box).  
    ~ Example 1 (id,reference ligand for selection of bs): /Applications/Chimera.app/Contents/MacOS/chimera --script "MolecularDocking.py --dock bsite --t 4ey7 --c A --r E20 --l 44159808 --v /usr/local/bin/vina".  
    ~ Example 2 (file, reference residues for selection of bs): /Applications/Chimera.app/Contents/MacOS/chimera --script "MolecularDocking.py --dock bsite --t target.cif --c A --r 124.a,297.a,202.a,86.a,286.a,337.a,341.a --l ligand.sdf --v /usr/local/bin/vina".  
