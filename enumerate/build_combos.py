
import json 
import rdkit
from rdkit import Chem
import os
from rdkit.Chem import AllChem
import pandas as pd
import sys
import csv

#FILE OF PAIRS, LIB FOR DICT, RXN

class Reaction:
    import csv 
    def grab_rxn(self,path):
        with open(path,'r') as rfile:
            reader = csv.reader(rfile,delimiter='\t')
            for row in reader:
                return row[2]
    def grab_name(self,path):
        with open(path,'r') as rfile:
            reader = csv.reader(rfile,delimiter='\t')
            for row in reader:
                return row[0]

    def __init__(self,path) -> None:
        self.scheme = self.grab_rxn(path)
        self.name = self.grab_name(path)


#
def combo_2_mol(pair,three_comp):
    if three_comp == False:
        s1 = Chem.MolFromSmiles(pair[0])
        s2 = Chem.MolFromSmiles(pair[1])
        return(s1,s2)
    else:
        s1 = Chem.MolFromSmiles(pair[0])
        s2 = Chem.MolFromSmiles(pair[1])
        s3 = Chem.MolFromSmiles(pair[2])
        return(s1,s2,s3)

def run_rxn(pair,three_comp,reaction=Reaction):
    rxn_input = AllChem.ReactionFromSmarts(reaction.scheme)
    return Chem.MolToSmiles(rxn_input.RunReactants(combo_2_mol(pair,three_comp))[0][0])

def gen_code(pair,dict,three_comp,reaction=Reaction):
    if three_comp == False:
        return f'{reaction.name}_{dict[pair[0]]}_{dict[pair[1]]}'
    else:
        return f'{reaction.name}_{dict[pair[0]]}_{dict[pair[1]]}_{dict[pair[2]]}'

def main():
    #json of id combos to enumerate
    file = sys.argv[1]
    #Trimmed lib of just reaction specific synthons
    trim_lib = sys.argv[2]
    #Reaction scheme
    rxn = sys.argv[3]
    #Save_dir
    directory_save=sys.argv[4]


    #Create DataFrame of trimmed library
    df = pd.read_csv(os.path.join(os.getcwd(),trim_lib),header=None,sep='\t')
    #Grab key_value_pairs of synthon,smile from df
    key_value_pairs = dict(zip(df.iloc[:,0], df.iloc[:,1]))

    #Instantiate reaction class
    rxn = Reaction(os.path.join(os.getcwd(),rxn))


    #Read through json file and grab all id combos
    with open(os.path.join(os.getcwd(),file),'r') as rfile:
        combos = json.load(rfile)



    #Where we will save enumerated smiles

    #Path to save smiles file
    
    base_name = os.path.basename(file)

    save_path = os.path.join(os.getcwd(),directory_save)
    save_file = os.path.join(save_path,f'{base_name}_enumerated.smi')
    print(save_path)
    print(save_file)


    #Open save path and write out smile and supplier code
    with open(save_file,'w') as wfile:
        for pair in combos:
            if len(pair)==3:
                three_comp=True
            else:
                three_comp=False
            line = f'{run_rxn(pair,three_comp,rxn)}\t{gen_code(pair,key_value_pairs,three_comp,rxn)}\n'
            wfile.write(line)

if __name__ == "__main__":
    main()





    




