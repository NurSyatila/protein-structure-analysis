# protein-structure-analysis
Scripts to parse PDB file and analyse protein structures and sequences

1) transform_pdb.py: 
  - A program to transform a PDB-formatted structure (pattern or whole PDB), given the matrix (rotation, translation)
  - Requires Bio module (Biopython)
  - Input: inputpdb, matrixlist, outputpdb
  - Usage: transform_pdb.py inputpdb matrixlist outputpdb
  - E.g. transform_pdb.py 1a0s.pdb 0.3786,-0.7283,-0.5712,0.8695,0.4913,-0.0500,0.3170,-0.4777,0.8193,93.254,19.092,37.797 1a0s_transformed.pdb
  
2) get_first_nmr_model.py
  - A program to get the first NMR model (Model 1) from NMR structure
  - Input: inputpdb, outputpdb
  - Usage: get_first_nmr_model.py inputpdb outputpdb
  - E.g. get_first_nmr_model.py 1g03.pdb 1g03_nmr.pdb
  
3) generate_binding_site.py
  - A program to create ligand binding sites
  - Input: pdbfile, ligid, distance_cutoff
  - Usage: generate_binding_site.py pdbfile ligid distance_cutoff
  - E.g. generate_binding_site.py 1a4g.pdb ZMR 4.0
  
4)
