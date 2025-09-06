import csv
import sys
import pandas as pd
import os
from itertools import product
import json

def chunk_list(lst, chunk_size):
    """Split the list into chunks."""
    for i in range(0, len(lst), int(chunk_size)):
        yield lst[i:i + int(chunk_size)]


def save_chunk(chunks,lib,save_path):
    for i,chunk in enumerate(chunks):
        print(f'Saving chunk {i}')
        save = f'{lib.split(".")[0]}_chunk-{i}.json'
        savefile_path = os.path.join(save_path,save)
        with open(savefile_path,'w')as file:
            json.dump(chunk,file)

def main():
    lib = sys.argv[1]
    chunk_size = sys.argv[2]
    directory_save = sys.argv[3]

    save_path = os.path.join(os.getcwd(),directory_save)
    os.makedirs(save_path,exist_ok=True)

    path = os.path.join(os.getcwd(),lib)

    #Generate two sets. One for each component of the reaction type.
    df = pd.read_csv(path,header=None,sep='\t')

    s1 = set(df[df.iloc[:,2] == 1].iloc[:,0])
    s2 = set(df[df.iloc[:,2] == 2].iloc[:,0])
    df3 = df[df.iloc[:,2] == 3]
    if  df3.empty:
        #Dealing with two component product
        all_combos = list(product(s1,s2))
        save_chunk(chunk_list(all_combos,chunk_size),lib,save_path)

    else:
        #Dealing with three component product
        s3 = set(df[df.iloc[:,2] == 3].iloc[:,0])
    
        all_combos = list(product(s1,s2,s3))
        save_chunk(chunk_list(all_combos,chunk_size),lib,save_path)



if __name__ == "__main__":
    main()










    