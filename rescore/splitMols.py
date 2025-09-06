from openeye import oechem
from openeye.oechem import *
import pandas as pd
import argparse
import os
import pandas as pd
import numpy as np

def main(mol2, chunk, prefix):

    ifs = oechem.oemolistream()
    if not ifs.open(mol2):
        oechem.OEThrow.Fatal('Unable to open mol2')

    if chunk == 1:

        for mol in ifs.GetOEMols():
            save_title = mol.GetTitle().split()[0]
            save_mol = oechem.OEGraphMol(mol)

            save_name = f'{save_title}.mol2'
            with oemolostream(save_name) as ofs:
                OEWriteMolecule(ofs, save_mol)

    else:

        listy = []

        for mol in ifs.GetOEMols():
             listy.append(oechem.OEGraphMol(mol))

        chunk_size = len(listy) // chunk + (len(listy) % chunk > 0)
        chunks = [listy[i:i + chunk_size] for i in range(0, len(listy), chunk_size)]

        for i,chunk in enumerate(chunks):

            filename = f'{prefix}__{i+1}.mol2'
            with oechem.oemolostream(filename) as ofs:
                for mol in chunk:
                    oechem.OEWriteMolecule(ofs,mol)
            
    return None

if __name__ == '__main__':

    #Initiate parser and add arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('mol2File', type=str, help='mol2 file')
    parser.add_argument('splitNum', type=int, help='# to chunk mol into', default=1, nargs='?')
    parser.add_argument('prefix', type=str, help='prefix to name  mol chunks', default='chunk', nargs='?')


    args = parser.parse_args()

    main(args.mol2File, args.splitNum, args.prefix)