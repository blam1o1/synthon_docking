from openeye import oechem
from openeye.oechem import *
import pandas as pd
import argparse
import os

synthon_data = '/path/to/synthon_data.txt'

def prepareSynthonData(path):
    df = pd.read_csv(path)
    return dict(zip(df['id'],zip(df['smiles'],df['component'])))

def patternMatch(mol, synthonSmiles):

    ''' Match mol onto synthonSmiles and then return match '''

    #patternMol = oechem.OEGraphMol(mol)

    ss = oechem.OESubSearch(synthonSmiles)
    if not ss.Init(synthonSmiles):

        oechem.OEThrow.Fatal("Unable to parse SMILES: %s" % synthonSmiles)

    oechem.OEAddExplicitHydrogens(mol)
    oechem.OEPrepareSearch(mol, ss)
    matches = list(ss.Match(mol))

    if matches:
        return matches[0]
    else:
        for atom in mol.GetAtoms():
            atom.SetAromatic(False)
    
        for bond in mol.GetBonds():
            bond.SetAromatic(False)

        matches = list(ss.Match(mol))
        
        if matches:
            return matches[0]

        else:
            return None
                       

def generateMol(match, mol, synthonData, synthonOneID, moleculeName):

    ''' Create mol from match atoms and original mol  and return new_mol'''
    
    match_atoms = [ma.target for ma in match.GetAtoms()]
    atmIdx = [atom.GetIdx() for atom in match_atoms]
    
    original_atoms = [molat for molat in mol.GetAtoms() if molat.GetIdx() in atmIdx]
    hydrogens = [hat for hat in mol.GetAtoms() if hat.GetAtomicNum() == 1 and any(neighbor.GetIdx() in atmIdx for neighbor in hat.GetAtoms())]
    
    original_atoms = original_atoms + hydrogens

    new_mol = OEGraphMol()
    atom_mapping = {}
    
    # Copy atoms and their 3D coordinates
    for atom in original_atoms:

        new_atom = new_mol.NewAtom(atom.GetAtomicNum())
        new_atom.SetFormalCharge(atom.GetFormalCharge())
        new_atom.SetIsotope(atom.GetIsotope())
        new_atom.SetMapIdx(atom.GetMapIdx())
        new_atom.SetName(atom.GetName())
        new_atom.SetType(atom.GetType())
        atom_mapping[atom] = new_atom

        # Copy 3D coordinates
        coords = mol.GetCoords(atom)
        new_mol.SetCoords(new_atom, coords) 

    # Copy bonds only between matched atoms from mol
    for bond in mol.GetBonds():
        a1 = bond.GetBgn()
        a2 = bond.GetEnd()
        if a1 in atom_mapping and a2 in atom_mapping:
            new_bond = new_mol.NewBond(atom_mapping[a1], atom_mapping[a2], bond.GetOrder())
            new_bond.SetType(bond.GetType())
    
    new_mol.SetTitle(f'extract_{moleculeName}_{synthonOneID}')

    OEAssignMDLHydrogens(new_mol)
    
    #print(f"New molecule has {new_mol.NumAtoms()} atoms and {new_mol.NumBonds()} bonds")    
    hydrogen_count = sum(1 for atom in new_mol.GetAtoms() if atom.GetAtomicNum() == 1)
    #print(f"Number of hydrogen atoms in target molecule: {hydrogen_count}")
    oechem.OEAddExplicitHydrogens(new_mol)
    oechem.OESet3DHydrogenGeom(new_mol)
    hydrogen_count = sum(1 for atom in new_mol.GetAtoms() if atom.GetAtomicNum() == 1)
    #print(f"Number of hydrogen atoms in target molecule: {hydrogen_count}")

    OEAssignAromaticFlags(new_mol)

    return new_mol

def main(in_mol2):

    synthonData = prepareSynthonData(synthon_data)

    ''' Instantiate mol2 iterator '''
    ifs = oechem.oemolistream()
    if not ifs.open(in_mol2):
        oechem.OEThrow.Fatal("Unable to open")

    synthonOneMol2 = []
    synthonTwoMol2 = []

    ''' Iterate through mol2 and pattern match both SYN1 and SYN2 '''

    successful = 0
    failure = 0
    fail_mols = []

    cwd = os.getcwd()
    #Create synthon mol dir
    os.mkdir(os.path.join(cwd,'extractedMols'))




    for mol in ifs.GetOEMols():

        #Grab molecule name
        moleculeName = mol.GetTitle().split('.')[0]
        print(moleculeName)

        #Grab synthon info
        synthonOneID = moleculeName[2:5]
        synthonTwoID = moleculeName[5:8]
        
        ids = [synthonOneID,synthonTwoID]
    


        for id in ids:

            listdirs = os.listdir(os.path.join(cwd,'extractedMols'))
            if id not in listdirs:
                os.mkdir(os.path.join(cwd,'extractedMols',id))

            print(id)
            copy_mol = oechem.OEGraphMol(mol)

            print(id in synthonData)
            if id not in synthonData:
                #print('here')
                continue
            
            match = patternMatch(copy_mol, synthonData[id][0])
            print(match)


            if match:
                successful += 1
                extractMol = generateMol(match, copy_mol,  synthonData, id, moleculeName)
            else: 
                failure += 1
                fail_mols.append(oechem.OEGraphMol(mol))
                print(f'Failed_match -> Molecule: {moleculeName} | SynthonID: {id} | Component: {synthonData[id][1]}')
                continue

            filename = f"{os.path.join(os.getcwd(),'extractedMols',id)}/{extractMol.GetTitle()}_{synthonData[id][1]}.mol2"
            with oemolostream(filename) as ofs:
                OEWriteMolecule(ofs, extractMol)

    #print(fail_mols)
    ofs = oechem.oemolostream("./failure.mol2")
    for fail_pose in fail_mols:
        oechem.OEWriteMolecule(ofs, fail_pose)
    ofs.close()
    print(f' Success: {successful} | Fails {failure}')


if __name__ == '__main__':

    #Initiate parser and add_arguments
    parser = argparse.ArgumentParser()
    
    parser.add_argument('input_mol', type=str, help='input mol2 file')
    #parser.add_argument('output_mol', type=str, help='output mol2 file')

    args = parser.parse_args()

    in_mol2 = args.input_mol
    
    main(in_mol2)


