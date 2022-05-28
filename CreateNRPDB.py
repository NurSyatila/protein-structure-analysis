#!/usr/bin/env python
# Author: NSAG
# Aim: To create a data set of non-redundant PDB from a list of representative PDB entities
# Based on results of the weekly clustering of protein sequences in the PDB by MMseqs2 at 30%, 40%, 50%, 70%, 90%, 95%, and 100% sequence identity

import sys,os,shutil,argparse,glob,collections,requests,re,subprocess
import urllib.request as request
from contextlib import closing

#   1) Download the list of sequence clusters from RCSB PDB
def DownloadClusterFile(SeqId):
    URL = 'https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-{}.txt'.format(SeqId)
    response = requests.get(URL)
    open('clusters-by-entity-{}.txt'.format(SeqId), 'wb').write(response.content)
    
#   2) Remove PDB bundle / large structures (large mmcif files / structures with > 62 chains)
def RemoveLargeStructures(SeqId):
    ClusterFile = 'clusters-by-entity-{}.txt'.format(SeqId)
    if os.path.isfile('clusters-by-entity-{}.txt'.format(SeqId)):
        # Download PDB bundle IDs
        URL = 'ftp://ftp.wwpdb.org/pub/pdb/compatible/pdb_bundle/pdb_bundle_index.txt'
        with closing(request.urlopen(URL)) as response:
            with open('pdb_bundle_index.txt', 'wb') as file:
                shutil.copyfileobj(response, file)
        # Create a new list without large structures (structures > 62 chains)
        if os.path.isfile('pdb_bundle_index.txt'):
            PDBBundles = [p.replace('\n','').upper() for p in open('pdb_bundle_index.txt','r').readlines()]
            NewClusterList = [[c for c in cl.replace('\n','').split() if c[:4] not in PDBBundles] for cl in open(ClusterFile,'r').readlines()]
            NewClusterList = [cl for cl in NewClusterList if cl !=[]]
        else:
            print ('Error: MMCIF Bundle File cannot be downloaded.')
            NewClusterList = []
    else:
        print ('Error: Sequence Cluster file cannot be not found.')
        NewClusterList = []
    return (NewClusterList)

#   3) For individual clusters, sort by xray resolution and get representative entities and IDs
def SortByResolution(NewClusterList,SeqId):
    if NewClusterList !=[]:
        # Download PDB bundle IDs
        URL = 'https://ftp.wwpdb.org/pub/pdb/derived_data/index/resolu.idx'
        response = requests.get(URL)
        open('resolution.txt', 'wb').write(response.content)
        if os.path.isfile('resolution.txt'):
            Resolution = [[r.split('\t')[0]]+[r.split('\t')[-1]] for r in re.findall(r'[0-9][a-zA-Z0-9][a-zA-Z0-9][a-zA-Z0-9][\t][;][\t].+',open('resolution.txt').read())]
            ResoDict = collections.defaultdict(list)
            for r in Resolution: ResoDict[r[0]].append(r[1])
            ResoClusterList = [[c+'|'+ResoDict[c[:4]][0] if c[:4] in ResoDict.keys() else c+'|'+'-2.0' for c in cl] for cl in NewClusterList]
            SortedResoClusterList = [sorted(s,key=lambda x:float(x.split('|')[1]),reverse=True) for s in ResoClusterList]
            RepresentativeEntities = [s[0].split('|')[0] for s in SortedResoClusterList]
            with open('{}-NRPDBEntities.txt'.format(SeqId),'w') as EntityListOut: EntityListOut.write('\n'.join(sorted(list(set(RepresentativeEntities))))) #comment this line for testing purpose
            #with open('{}-NRPDBEntities.txt'.format(SeqId),'w') as EntityListOut: EntityListOut.write('\n'.join(sorted(list(set(RepresentativeEntities[:10]))))) #uncomment this line for testing purpose
            print ('Representative PDB Entities: ',len(sorted(list(set(RepresentativeEntities)))))
            with open('{}-NRPDBIDs.txt'.format(SeqId),'w') as PDBListOut: PDBListOut.write(','.join(sorted(list(set([s[:4] for s in RepresentativeEntities]))))) #comment this line for testing purpose
            #with open('{}-NRPDBIDs.txt'.format(SeqId),'w') as PDBListOut: PDBListOut.write(','.join(sorted(list(set([s[:4] for s in RepresentativeEntities][:10]))))) #uncomment this line for testing purpose
            print ('Representative PDB IDs: ',len((sorted(list(set([s[:4] for s in RepresentativeEntities]))))))
        else:
            print ('Error: Resolution File cannot be downloaded.')
            SortedResoClusterList = []
    else:
        print ('Error: Filtered Sequence Clusters cannot be generated.')
        SortedResoClusterList = []
    return SortedResoClusterList

#   4) Download the PDBs (PDB format)
def DownloadPDBFiles(SeqId):
    if not os.path.exists('pdb'): os.mkdir('pdb')
    DownloadList = '{}-NRPDBIDs.txt'.format(SeqId)
    p = subprocess.Popen('./BatchDownload.sh -f {} -p -o pdb'.format(DownloadList), shell=True)
    sStdout, sStdErr = p.communicate(input=None,timeout=None)
    os.chdir('pdb')
    for pdb in glob.glob('*.pdb'): os.rename(pdb, pdb.lower())
    os.chdir('../')

#   5) Get the first model from each NMR structures
def NMRCheck():
    for PDBfile in glob.glob('pdb/*.pdb'):
        PDBread = open(PDBfile, 'r').read()
        if 'MODEL        1' in PDBread:
            PDBLines = PDBread.split('MODEL        2')[0].split('\n')
            NewPDB = open(PDBfile, 'w')
            for l in PDBLines: NewPDB.write(l+'\n')

#    6) Remove PDB that contains only CA atoms (Only select PDB with side chains)
def PDBCheck(SortedResoClusterList, SeqId):
    ClusterDict = collections.defaultdict(list)
    PDBDict = collections.defaultdict(list)
    for index,s in enumerate(SortedResoClusterList):
        ClusterDict['cluster_'+str(index)].append(s)
        for x in s:
            PDBDict[x.split('|')[0]].append('cluster_'+str(index))
    if not os.path.exists('pdbCAonly'): os.mkdir('pdbCAonly')
    for p in glob.glob('pdb/*.pdb'):
        plines = [line for line in open(p,'r').readlines() if line[:6]=='MDLTYP']
        if plines !=[]:
            if 'CA ATOMS ONLY' in plines[0]:
                shutil.move(p,p.replace('pdb/','pdbCAonly/'))

#    7) Generate a directory with representative PDB chains based on entity IDs
def CreateNRPDBDirectory(SeqId):
    # get list of representative entities
    NRDirectory = 'pdb'+SeqId
    if not os.path.exists(NRDirectory): os.mkdir(NRDirectory)
    RepresentativeEntities = [rep.replace('\n','') for rep in open('{}-NRPDBEntities.txt'.format(SeqId),'r').readlines()]
    PDBDict = collections.defaultdict(list)
    for rep in RepresentativeEntities: PDBDict[rep[:4].lower()].append(rep[5:])
    #for each pdb, get only representative chains based on entity IDs (ATOM and HETATM)
    SummaryFile = 'NRPDBSummary.txt'
    with open(SummaryFile,'a') as SumOut:
        for r in sorted(PDBDict.keys()):
            PDBOri = 'pdb/'+r+'.pdb'
            PDBNew = '{}/'.format(NRDirectory)+r+'.pdb'
            if os.path.isfile(PDBOri):
                NREntities =  ['MOL_ID: '+e+';' for e in PDBDict[r]]
                EntityLines = ''.join([line for line in open(PDBOri,'r').readlines() if line[:6]=='COMPND'])
                NRChains = []
                for e in NREntities:
                    Chain = re.findall(r'{}.*\n.*MOLECULE.*\n.*CHAIN:.*'.format(e),EntityLines,re.M)[0].split('CHAIN: ')[1].replace(';','').split(', ')[0]
                    NRChains.append(Chain)
                NRChains = sorted(list(set(NRChains)))
                print (r, 'NR PDB chain/s: ', ','.join(NRChains))
                #get chains to be removed
                PDBLines = open(PDBOri,'r').readlines()
                ChainsToRemove = sorted(list(set([line[21:22] for line in PDBLines if line[:4]=='ATOM' or line[:6]=='HETATM'])-set(NRChains)))
                if ChainsToRemove == []: shutil.copy(PDBOri,PDBNew)
                if ChainsToRemove !=[]:
                    LinesToRemove = [line for line in PDBLines if (line[:4]=='ATOM' or line[:6]=='HETATM') and line[21:22] in ChainsToRemove]
                    for l in LinesToRemove: PDBLines.remove(l)
                    with open(PDBNew,'w') as PDBOut:
                        for line in PDBLines: PDBOut.write(line)
                SumOut.write(r+'\t'+','.join(NRChains)+'\n')
            else:
                if os.path.isfile('pdbCAonly/'+r+'.pdb'):
                    SumOut.write(r+'\t'+'CA ATOMS ONLY.\n')
                else:
                    SumOut.write(r+'\t'+'PDB NOT FOUND.\n')
                
def GenerateNRPDB(SeqId):
    print ('Get NR for {}% sequence identity'.format(SeqId))
    print ('Download cluster file from RCSB PDB')
    DownloadClusterFile(SeqId)
    print ('Remove large structures (select only PDB with < 62 chains')
    NewClusterList = RemoveLargeStructures(SeqId)
    print ('Sort and select representatives based on resolution value')
    SortedResoClusterList = SortByResolution(NewClusterList,SeqId)
    print ('Download PDB files. Download will takes time.')
    DownloadPDBFiles(SeqId)
    print ('Select only the first model for NMR structures')
    NMRCheck()
    print ('Select PDB with side chain atoms, discard PDB with only CA')
    PDBCheck(SortedResoClusterList, SeqId)
    print ('Generate NR directory')
    CreateNRPDBDirectory(SeqId)
    print ('Done!')

# Create the parser
my_parser = argparse.ArgumentParser(description='Create a Non-redundant PDB Directory (Side Chains) based on the clustering of protein sequences in the PDB by MMseqs2.')

# Add the arguments
my_parser.add_argument('SeqId',type=str,help='Sequence Identity: 30, 40, 50, 70, 90, 95, or 100')

# Execute the parse_args() method
args = my_parser.parse_args()

GenerateNRPDB(args.SeqId)

#e.g. Python CreateNRPDB.py 30
