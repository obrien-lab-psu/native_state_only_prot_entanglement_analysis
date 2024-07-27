from Bio import PDB
import argparse
import os
import numpy as np
import pandas as pd
import sys


def scan_disulfide_bonds(pdb_file: str) -> None:

    """
    determines all contacts of side chain heavy atoms that are 4.5A apart or less and are atleast 4 residues apart along the primary structure

    """

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    
    contacts = []

    # Iterate over each model
    for model in structure:
        # Iterate over each chain
        for chain in model:
            residues = list(chain)
            num_residues = len(residues)

            for i in range(num_residues):
                res1 = residues[i]
                res_id1 = res1.get_id()[1]
                res_name1 = res1.get_resname()
                atoms1 = [atom for atom in res1.get_atoms() if atom.get_id()[0] not in 'H']  # Exclude hydrogen atoms

                for j in range(i + 1, num_residues):
                    res2 = residues[j]
                    res_id2 = res2.get_id()[1]
                    res_name2 = res2.get_resname()
                    atoms2 = [atom for atom in res2.get_atoms() if atom.get_id()[0] not in 'H']  # Exclude hydrogen atoms

                    min_distance = float('inf')

                    for atom1 in atoms1:
                        coords1 = atom1.get_coord()
                        for atom2 in atoms2:
                            coords2 = atom2.get_coord()
                            distance = calculate_distance(coords1, coords2)
                            if distance < min_distance:
                                min_distance = distance

                    if min_distance <= 4.5:
                        if abs(res_id2 - res_id1) > 4:
                            contacts.append({
                                'chainID': chain.get_id(),
                                'resid 1': res_id1,
                                'resname 1': res_name1,
                                'resid 2': res_id2,
                                'resname 2': res_name2,
                                'min_dist': min_distance
                            })

                            ## Check if there is a CYS - CYS bond
                            if res_name1 == 'CYS' and res_name2 == 'CYS':
                                print(f'FOUND CYS-CYS @ {res_id1}-{res_id2} with a min dist of {min_distance} in chain {chain.get_id()} in {pdb_file}')

    return pd.DataFrame(contacts)
        

def calculate_distance(atom1, atom2):
    """
    Calculates the Euclidean distance between two 3D points.
    """
    return np.sqrt((atom1[0] - atom2[0])**2 + (atom1[1] - atom2[1])**2 + (atom1[2] - atom2[2])**2)



##########################################################################
if __name__ == "__main__":

    # parse user arguments
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--PDB", type=str, required=True, help="Path to PDB file you want to generate raw entanglments for")
    args = parser.parse_args()

    # parse some of the default parameters
    global pdb_dir
    pdb_dir = args.PDB
    
    if os.path.isdir(pdb_dir):
        input_data = [f'{pdb_dir}{f}' for f in os.listdir(f"{pdb_dir}")]
    elif os.path.isfile(pdb_dir):
        input_data = [pdb_dir]

    for f in input_data:
        print(f)
        results = scan_disulfide_bonds(f)
        print(results)
        out = f.replace('.pdb', '_SSheavy_contacts.pdb')
        results.to_csv(out, index=False)
        print(f'SAVED: {out}')

    print(f'NORMAL TERMINATION')

  
