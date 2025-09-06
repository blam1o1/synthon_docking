from openeye import oechem
from openeye.oechem import *
import pandas as pd
import argparse
import os
import pandas as pd

synthon_data = '/lustre/fs6/lyu_lab/scratch/blam/work/synthon/rescorePoses/synthonDict/m_22bcb_RescoreSynthon.csv'

def grabPose(id, mol2):
    
    ''' Grab mol from mol2 file based on title of mol '''

    ifs = oechem.oemolistream()
    if not ifs.open(mol2):
        oechem.OEThrow.Fatal('Unable to parse mol2 for grabbing Pose')

    #Iterate through mols and check if title matches id, if it does return it
    for mol in ifs.GetOEMols():
        title = mol.GetTitle()
        if title == id:
            return mol
    
    return None


def main(fullMoleculeSynthonMol2, synthonMol2):

    ''' Iterate through fullMoleculeSynthonMol2 and lookup in synthonMol2 and calculate Mol2 '''

    #Initiate mol2 iterator
    ifs = oechem.oemolistream()
    if not ifs.open(fullMoleculeSynthonMol2):
        oechem.OEThrow.Fatal('Unable to parse fullMol2')
    
    rmsdEntries = []

    for mol in ifs.GetOEMols():

        #Grab mol title and extract synthonID
        
        fullMcule, synthonTitle = mol.GetTitle().split('_')[1:]

        print(fullMcule)
        print(synthonTitle)
        
        synthon_mol = grabPose(synthonTitle,synthonMol2)

        #If synthon can be found in synthonMol2, calculate RMSD

        if synthon_mol:
            rmsd = OERMSD(synthon_mol, mol)
            entry = {'fullMolecule': fullMcule, 'synthon': synthonTitle, 'RMSD': rmsd}

        #Else, return None

        else:
            entry = {'fullMolecule': fullMcule, 'synthon': synthonTitle, 'RMSD': None}  

        #Append entry to entriesList 

        rmsdEntries.append(entry)


    #Create DataFrame and save
    
    save_name = fullMoleculeSynthonMol2.split('/')[-1].split('.')[0]

    df = pd.DataFrame(rmsdEntries)
    df.to_csv(f'./{save_name}_rmsd_out.csv',index=False)

    return None


if __name__ == '__main__':

    #Initiate parser and add arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('fullMol2', type=str, help='Fullmolecule mol2')
    parser.add_argument('synthonMol2', type=str, help='Synthon mol2')

    args = parser.parse_args()

    synthonMol2 = os.path.join(os.getcwd(), args.synthonMol2)
    fullMol2 = os.path.join(os.getcwd(), args.fullMol2)

    main(fullMol2, synthonMol2)

    