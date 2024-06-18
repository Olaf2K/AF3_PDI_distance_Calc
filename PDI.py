
import os
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from scipy.spatial.distance import cdist
from ase.io import read, write

workdir = r'' # set workdir
os.chdir(workdir)
dirsOfInterst = [r''] # dir with all your cif files

# ----------- Fucntions ------------ #
def ConvertCIFtoPDB(inputfile):
    inputfile =  inputfile.split('.cif')[0]
    parser = MMCIFParser()
    structure = parser.get_structure("alphafold_structure", inputfile+".cif")
    io = PDBIO()
    io.set_structure(structure)
    io.save(inputfile+".pdb")
    return inputfile

def readPDB(inputfile):
    parser = PDBParser()
    structure = parser.get_structure("protein_dna_complex", inputfile+".pdb")
    return structure

def findSequenceIndexEffector(lowerlimit,upperlimit,chain_id):
    ProteinMotifResidues = list(range(lowerlimit,upperlimit))
    return ProteinMotifResidues 


def findSequenceIndexDNA(DNAbindingMotif,chain_id):
    chain = structure[0][chain_id]
   
    residue_ids = []
    residueString = ''
    for residue in chain:
        residueString+=re.split('',residue.get_resname())[-2]
    dna_motif_residuesB = residueString.find(DNAbindingMotif)
    if dna_motif_residuesB == -1:
        
        DNAbindingMotif = reverse_complement(DNAbindingMotif)
        residue_ids = []
        residueString = ''
        for residue in chain:
            residueString+=re.split('',residue.get_resname())[-2]
        dna_motif_residuesB = residueString.find(DNAbindingMotif) 
    else:
        pass 
    dna_motif_residuesB = range(dna_motif_residuesB,dna_motif_residuesB+len(DNAbindingMotif))
    return dna_motif_residuesB 

def ExtractAtomCoords(structure, chainOfInterest,RangeOfInterest):
    atoms = []
    for model in structure:
        for chain in model:
            if chain.id == chainOfInterest: 
                for residue in chain:
                    if residue.id[1] in RangeOfInterest:
                        for atom in residue:
                            atoms.append(atom.get_coord())
    return atoms

def calcDistances(atoms_dna_motif,atoms_other_part):
    distances = calculate_distances(atoms_dna_motif, atoms_other_part)
    min_distance = distances.min()
    max_distance = distances.max()
    avg_distance = distances.mean()
    return avg_distance,max_distance,min_distance

def reportLowestOnly(min_distanceA,min_distanceB):
    avg_distanceA = min_distanceA[0]
    max_distanceA = min_distanceA[1]
    min_distanceA = min_distanceA[2]
    avg_distanceB = min_distanceB[0]
    max_distanceB = min_distanceB[1]
    min_distanceB = min_distanceB[2]
    if min_distanceA < min_distanceB:
        return avg_distanceA,max_distanceA,min_distanceA    
    else:
        return avg_distanceB,max_distanceB,min_distanceB  
        
from natsort import natsorted
def plotandsave(infile):
    import seaborn as sns
    import pandas as pd
    import matplotlib.pyplot as plt
    df = pd.read_csv(infile)

    long_df = df.melt(id_vars=['inputfileBasename', 'model'], 
                    value_vars=['AverageDist', 'MaxDist', 'MinDist'], 
                    var_name='DistanceType', 
                    value_name='DistanceValue')
    long_df['inputfileBasename'] = pd.Categorical(
        long_df['inputfileBasename'],
        categories=natsorted(long_df['inputfileBasename'].unique()),
        ordered=True
    )


    fig = sns.catplot(
        data=long_df, 
        x='inputfileBasename', 
        y='DistanceValue',
        hue='inputfileBasename', 
        col='DistanceType', 
        kind='boxen',
        col_wrap=1, 
        sharey=False,  
        height=5,   
        aspect=0.5   
    )


    for ax in fig.axes.flat:
        for label in ax.get_xticklabels():
            label.set_rotation(90)


    plt.tight_layout()

    plt.show()
    fig.savefig(infile+"seq.svg") 
    plt.close()

    fig = plt.figure(figsize=(3, 5))
    ax = sns.boxplot(
        data=long_df, 
        x='DistanceType', 
        y='DistanceValue',
        hue='inputfileBasename'
    )


    plt.tight_layout()




    fig.figure.savefig(infile+"AAseqCombined.svg") 
    plt.close()
    fig = sns.violinplot(
        data=long_df, x='DistanceType', y='DistanceValue',
        hue='DistanceType')
    fig.figure.savefig(infile+"CombinedViolin.svg") 
    plt.close()

    fig = sns.displot(
    data=long_df[long_df.DistanceType == 'MinDist'], 
    x="DistanceValue", 
    hue='inputfileBasename', 
    kde=True, 
    col='inputfileBasename', 
    col_wrap=1
)
    fig.figure.savefig(infile+"KDE.svg") 



# ----------- Main ------------ #

import re
from Bio.Seq import Seq


for i in dirsOfInterst:
    os.chdir(i)

    directory = os.getcwd()
    files = os.listdir(directory)
    cif_files = [file for file in files if file.endswith(".cif")]
    infileName = 'DistancesBindingDomainFromMotif10AAseq.csv'
    with open(infileName,'w') as f:
        f.write(','.join(['inputfileBasename','model', 'AverageDist','MaxDist','MinDist'])+'\n')
        for CIFfile in cif_files:
            inputfile = ConvertCIFtoPDB(CIFfile)
            structure = readPDB(inputfile)
            print(inputfile)
            searchString = (inputfile.split('_model')[0]) 
            zip_fileName = [file for file in files if file.endswith(".zip") and searchString in file]
            TargetSequenceDNA =zip_fileName[0].split('.zip')[0].split('_')[-1]

            dna_motif_residuesA = findSequenceIndexEffector(0,10,'A')
            dna_motif_residuesB = findSequenceIndexDNA(TargetSequenceDNA,'B')
            dna_motif_residuesC = findSequenceIndexDNA(TargetSequenceDNA,'C')

            AtomLocA = ExtractAtomCoords(structure, 'A' ,dna_motif_residuesA)
            AtomLocB = ExtractAtomCoords(structure, 'B' ,dna_motif_residuesB)
            AtomLocC = ExtractAtomCoords(structure, 'C' ,dna_motif_residuesC)

            DistancesB = calcDistances(AtomLocB,AtomLocA)
            DistancesC = calcDistances(AtomLocC,AtomLocA)

            res = reportLowestOnly(DistancesB,DistancesC)

            basename = inputfile.split('_v1_')[0] 
            modelNO = inputfile.split('model_')[-1]
            f.write(','.join([basename,modelNO, str(res[0]),str(res[1]),str(res[2])])+'\n')
    plotandsave(infileName)








