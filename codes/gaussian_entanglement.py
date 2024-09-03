import itertools
import os
from operator import itemgetter
from warnings import filterwarnings
import numpy as np
import MDAnalysis as mda
import pandas as pd
from numba import njit
from scipy.spatial.distance import pdist, squareform
from topoly import lasso_type  # used pip
import re
import sys

filterwarnings("ignore")

def pre_processing_pdb(pdb_file: str) -> None:

    """
    Pre-processing the PDB files (.pdb, .pdb1 (model 1) and alphafold pdbs) by removing everything after the last TER. 

    Function removes all model instances except the first one in the case of .pdb1 files. 
    
    Why remove the last TER?
        Residues after the TER are non-protein

    """

    pdb_data = np.loadtxt(f"PDBs/{pdb_file}", dtype=str, delimiter="\n")

    if pdb_file in os.listdir("Before_TER_PDBs/"):
        print(f"\033[4m{pdb_file}\033[0m is already processed. Please look in Before_TER_PDBs/")
        return

    if pdb_data[1].split()[1] == "ALPHAFOLD":
        print(f"ALPHAFOLD PDB \033[4m{pdb_file}\033[0m DOES NOT REQUIRE PROCESSING; COPYING PDB FILE TO Before_TER_PDBs/")
        os.system(f"cp PDBs/{pdb_file} Before_TER_PDBs/{pdb_file}")
        return
    
    else:
        print(f"PROCESSING \033[4m{pdb_file}\033[0m")

    last_TER_index = 0
    Model_indices = []

    for i in range(len(pdb_data)):

        if pdb_data[i].startswith("MODEL"):
            Model_indices.append(i)

    if len(Model_indices) > 1:

        for j in range(Model_indices[1]):

            if pdb_data[j].startswith("TER"):

                if j > last_TER_index:
                    last_TER_index = j

    elif len(Model_indices) == 1 or len(Model_indices) == 0:

        for j in range(len(pdb_data)):

            if pdb_data[j].startswith("TER"):

                if j > last_TER_index:
                    last_TER_index = j

    if last_TER_index == 0 and len(Model_indices) == 1:

        last_TER_index = len(pdb_data)
    
    elif last_TER_index == 0 and len(Model_indices) > 1:

        last_TER_index = Model_indices[1]

    with open(f"Before_TER_PDBs/{pdb_file}", "w") as f:

        for k in range(last_TER_index):
        
            f.write(f"{pdb_data[k]}\n")

@njit(fastmath=True)
def helper_dot(Runit: np.ndarray, dR_cross: np.ndarray) -> list:

    """
    Numba function to speed up dot product calculation. Ability 
    to use current GPU (if available) and CPU
    
    """

    return [np.dot(x,y) for x,y in zip(Runit,dR_cross)]

def point_rounding(num: float) -> float:

    """
    This function perform rounding to the nearest 0.6. 

    Ex: 0.61 rounds up to 1 but 0.59 rounds down to 0

    """
    if len(str(num).split("e")) == 2: 
        num = 0.0

    if num % 1 >= g_threshold:
        rounded_num = round(num)
    else:
        rounded_num = int(str(num).split(".")[0])
    
    return rounded_num

def get_entanglements(coor: np.ndarray, l: int, termini_threshold: list, pdb_file: str, resids: np.ndarray, 
        resnames: np.ndarray,resid_index_to_ref_allatoms_idx: dict, ca_coor: np.ndarray, resid_index_to_resid: dict) -> dict:

    """
    Find proteins containing non-covalent lasso entanglements.

    Entanglements are composed of loops (defined by native contacts) and crossing residue(s).

    """
    Nterm_thresh = termini_threshold[0]
    Cterm_thresh = termini_threshold[1]

    # make native contact contact map
    dist_matrix = squareform(pdist(coor))
    if Calpha == False:
        native_cmap = np.where(dist_matrix <= 4.5, 1, 0) # if true then 1 will appear otherwise zero
    elif Calpha == True:
        native_cmap = np.where(dist_matrix <= 8, 1, 0) # if true then 1 will appear otherwise zero
    native_cmap = np.triu(native_cmap, k=4) # element below the 4th diagonal starting from middle are all zeros; # protein contact map

    num_res = len(resid_index_to_ref_allatoms_idx.keys())

    assert num_res == len(resids), f"Something's wrong with {pdb_file}"

    res_ncmap = np.zeros((num_res, num_res))
    resid_pairs = list(itertools.product(np.arange(num_res), np.arange(num_res)))

    for pair in resid_pairs:
        # check that the resid are greater than 4 apart
        if abs(resid_index_to_resid[pair[1]] - resid_index_to_resid[pair[0]]) > 4:
            if pair[0] in resid_index_to_ref_allatoms_idx and pair[1] in resid_index_to_ref_allatoms_idx:
                res1_atoms = resid_index_to_ref_allatoms_idx[pair[0]]
                res2_atoms = resid_index_to_ref_allatoms_idx[pair[1]]
                res1_atoms_start = min(res1_atoms)
                res1_atoms_end = max(res1_atoms)
                res2_atoms_start = min(res2_atoms)
                res2_atoms_end = max(res2_atoms)
                contact = np.sum(native_cmap[res1_atoms_start:res1_atoms_end + 1, res2_atoms_start:res2_atoms_end + 1])

                if contact > 0:
                    res_ncmap[pair[0], pair[1]] = 1
                    #print(f'Found contact: {resid_index_to_resid[pair[0]]} & {resid_index_to_resid[pair[1]]}')        
    del native_cmap
    native_cmap = res_ncmap 

    nc_indexs = np.stack(np.nonzero(native_cmap)).T # stack indices based on rows

    # make R coordinate and gradient of R length N-1
    range_l = np.arange(0, l-1)
    range_next_l = np.arange(1,l)

    ca_coor = ca_coor.astype(np.float32)
    R = 0.5*(ca_coor[range_l] + ca_coor[range_next_l])
    dR = ca_coor[range_next_l] - ca_coor[range_l]

    #make dRcross matrix
    pair_array = np.asarray(list(itertools.product(dR,dR))) # combination of elements within array

    x = pair_array[:,0,:]
    y = pair_array[:,1,:]

    dR_cross = np.cross(x, y)

    #make Rnorm matrix
    pair_array = np.asarray(list(itertools.product(R,R)))
    diff = pair_array[:,0,:] - pair_array[:,1,:]
    diff = diff.astype(np.float32)

    Runit = diff / np.linalg.norm(diff, axis=1)[:,None]**3 
    Runit = Runit.astype(np.float32)

    #make final dot matrix
    dot_matrix = helper_dot(Runit, dR_cross)
    dot_matrix = np.asarray(dot_matrix)
    dot_matrix = dot_matrix.reshape((l-1,l-1))

    nc_gdict = {} 

    for i,j in nc_indexs:

        loop_range = np.arange(i,j)
        nterm_range = np.arange(Nterm_thresh,i-5)
        cterm_range = np.arange(j+6,l-(Cterm_thresh + 1))

        gn_pairs_array = np.fromiter(itertools.chain(*itertools.product(nterm_range, loop_range)), int).reshape(-1, 2)
        gc_pairs_array = np.fromiter(itertools.chain(*itertools.product(loop_range, cterm_range)), int).reshape(-1, 2)

        if gn_pairs_array.size != 0:
            
            gn_vals = dot_matrix[gn_pairs_array[:,0],gn_pairs_array[:,1]]
            gn_vals = gn_vals[~np.isnan(gn_vals)] 
            gn_val = np.sum(gn_vals) / (4.0 * np.pi)
        
        else:
            gn_val = 0

        if gc_pairs_array.size != 0:
            
            gc_vals = dot_matrix[gc_pairs_array[:,0],gc_pairs_array[:,1]]
            gc_vals = gc_vals[~np.isnan(gc_vals)] 
            gc_val = np.sum(gc_vals) / (4.0 * np.pi)
        
        else:
            gc_val = 0
        
        rounded_gc_val = point_rounding(np.float64(abs(gc_val)))
        rounded_gn_val = point_rounding(np.float64(abs(gn_val)))

        if np.abs(rounded_gn_val) >= 1 or np.abs(rounded_gc_val) >= 1:
            #print(f'({i}, {j}) with gn: {gn_val} and gc: {gc_val}')
            nc_gdict[ (int(i), int(j)) ] = (gn_val, gc_val, rounded_gn_val, rounded_gc_val)

    missing_residues = find_missing_residues(resids)
    #print(f'missing_residues: {missing_residues}')

    filtered_nc_gdict = loop_filter(nc_gdict, resids, missing_residues)
    #print(f'size filtered_nc_gdict after accounting for missing residues: {len(filtered_nc_gdict)}')

    entangled_res = find_crossing(ca_coor.tolist(), filtered_nc_gdict, resids)
    #print(f'entangled_res: {entangled_res}')

    filtered_entangled_res = crossing_filter(entangled_res, missing_residues)
    #print(f'filtered_entangled_res: {filtered_entangled_res}')

    disulfide_keys = []
    for ent_idx, ent in enumerate(filtered_entangled_res):

        native_i = ent[0]
        native_j = ent[1]

        nc_i_resname = resnames[np.where(resids == native_i)][0]

        nc_j_resname = resnames[np.where(resids == native_j)][0]

        if nc_i_resname == "CYS" and nc_j_resname == "CYS":
            disulfide_keys += [ent]

    print(f'disulfide_keys: {disulfide_keys} to remove')
    # delet any keys with disulfide bonds
    for dk in disulfide_keys:
        del filtered_entangled_res[dk]

    return filtered_entangled_res, missing_residues


def find_missing_residues(resids:np.ndarray) -> np.ndarray:

    """
    Find missing residues in pdb file

    """

    check_all_resids = np.arange(resids[0], resids[-1] + 1)

    missing_residues = np.setdiff1d(check_all_resids, resids)

    return missing_residues

def loop_filter(native_contacts: dict, resids: np.ndarray, missing_res: np.ndarray) -> dict:

    """
    Remove loops if there are three or more consecutive missing residues
    or the amount of any missing residues exceed 5% of the loop length 

    """

    for ij, values in native_contacts.items():

        native_i = resids[ij[0]]

        native_j = resids[ij[1]]

        rounded_gn = values[-2]

        rounded_gc = values[-1]

        check_loop = np.arange(native_i , native_j + 1) 

        loop_length = check_loop.size

        missing_res_loop = np.intersect1d(check_loop, missing_res)

        for index, diff_resid_index in itertools.groupby(enumerate(missing_res_loop), lambda ix : ix[0] - ix[1]):

            conseuctive_missing_residues = list(map(itemgetter(1), diff_resid_index))

            if len(conseuctive_missing_residues) >= 3 or len(missing_res_loop) > 0.05 * loop_length:

                native_contacts[ij] = None

    return native_contacts

def find_crossing(coor: np.ndarray, nc_data: dict, resids: np.ndarray) -> dict:

    """
    Use Topoly to find crossing(s) based on partial linking number

    """

    entangled_res = {}

    native_contacts = [[ij[0], ij[1]] for ij, values in nc_data.items() if values is not None]

    # reduction:
        # 1. each crossing must be 10 residues apart [default]
        # 2. first crossing should be at least 6 residues from the loop
        # 3. first crossing should be at least 5 residues from the closest termini

    # note this functionality doesnt work whats 
    data = lasso_type(coor, loop_indices=native_contacts, more_info=True, density=density, min_dist=[10, 6, 5])
    # high precision, low denisty

    for native_contact in native_contacts:

        crossings = []

        native_contact = tuple(native_contact)

        if abs(nc_data[native_contact][-2]) >= 1: # if rounded_gn >= 1

            crossingN = [f"{cr[0]}{resids[int(cr[1:])]}" for cr in data[native_contact]["crossingsN"]]

            crossings += crossingN

        if abs(nc_data[native_contact][-1]) >= 1: # if rounded_gc >= 1
            
            crossingC = [f"{cr[0]}{resids[int(cr[1:])]}" for cr in data[native_contact]["crossingsC"]]
            
            crossings += crossingC

        gn = nc_data[native_contact][0]
        
        gc = nc_data[native_contact][1]

        ij_gN_gC = (resids[native_contact[0]], resids[native_contact[1]]) + (gn, gc) 

        entangled_res[ij_gN_gC] = np.unique(crossings)
        
    return entangled_res

def crossing_filter(entanglements: dict, missing_res: np.ndarray) -> dict:

    """
    Remove entanglements if there are any missing residues plus-and-minus 10 of the crossing(s)

    """
    
    for ij_gN_gC, crossings in entanglements.items():

        if crossings.size:

            check_crossing = []

            for crossing in crossings:

                reg_exp = re.split("\\+|-", crossing, maxsplit=1) # split the chirality

                check_crossing.append(np.arange(int(reg_exp[1]) - 10 , int(reg_exp[1]) + 11))

            check_crossing = np.concatenate(check_crossing)

            missing_res_cr = np.intersect1d(check_crossing, missing_res)

            if missing_res_cr.size:

                entanglements[ij_gN_gC] = None

    filtered_entanglements = {nc: re_cr for nc, re_cr in entanglements.items() if re_cr is not None and len(re_cr) > 0} 

    return filtered_entanglements

def calculate_native_entanglements(pdb_file: str) -> None:

    """
    Driver function that outputs native lasso-like self entanglements and missing residues for pdb and all of its chains if any

    """

    pdb = pdb_file.split('/')[-1].split(".")[0]

    if Calpha == False:
        if f"{pdb}_GE.txt" in os.listdir("unmapped_GE/"):
            print(f"\033[4m{pdb}\033[0m has non-covalent lasso entanglements. Please look in unmapped_GE/")
            return
    elif Calpha == True:
        if f"{pdb}_GE.txt" in os.listdir("unmapped_GE_CA/"):
            print(f"\033[4m{pdb}\033[0m has non-covalent lasso entanglements. Please look in unmapped_GE/")
            return


    ref_univ = mda.Universe(f"{pdb_file}", format="PDB")

    print(f"COMPUTING ENTANGLEMENTS FOR \033[4m{pdb}\033[0m")

    if Calpha == False:
        ref_allatoms_dups = ref_univ.select_atoms("not name H*")
    else:
        ref_allatoms_dups = ref_univ.select_atoms("name CA")

    termini_threshold = [5, 5]

    chains_to_analyze = set(ref_univ.segments.segids)

    for chain in chains_to_analyze:

        atom_data = []
        check = set()

        for atom in ref_allatoms_dups.select_atoms(f"segid {chain}").atoms:

            atom_data.append((atom.resid, atom.name))
            check.add((atom.resid, atom.name))

        temp_df = pd.DataFrame(atom_data, columns=["resid", "name"])

        unique_rows = temp_df.drop_duplicates()
        unique_indices = unique_rows.index.tolist()

        assert len(check) == len(unique_indices), "You did not remove dup atoms!"

        ref_allatoms_unique = ref_allatoms_dups.select_atoms(f"segid {chain}")[unique_indices]

        ref_ca_unique = ref_allatoms_unique.select_atoms("name CA")

        resid_index_to_ref_allatoms_idx = {}
        resid_index_to_resid = {}
        ref_allatoms_idx_to_resid = {}
        atom_ix = 0
        res_ix = 0
        PDB_resids = ref_ca_unique.resids
        new_atm_idx = []

        if len(PDB_resids) == 0 or len(PDB_resids) == 1:
            print(f"Skipping over chain {chain} for \033[4m{pdb}\033[0m since chain has only one alpha carbon or none")
            continue

        ref_allatoms_unique_list = []
        for resid in PDB_resids:
            sel = ref_allatoms_unique.select_atoms(f"resid {resid} or resid {resid}A or resid {resid}C")
            #print(f'\n{resid} {sel} {len(sel)}')
            if len(sel) == 0:
                raise ValueError(f'There was an empty atom group at resid: {resid} that is most likely due to resid names have letters next to them')
            ref_allatoms_unique_list += [sel]
        ref_allatoms_unique = mda.Merge(*ref_allatoms_unique_list)
        #ref_allatoms_unique = mda.Merge(*[ref_allatoms_unique.select_atoms(f"resid {resid} or resid {resid}A or resid {resid}C") for resid in PDB_resids])

        for atom in ref_allatoms_unique.atoms:
            new_atm_idx.append(atom_ix)
            ref_allatoms_idx_to_resid[atom_ix] = [atom.resid]

            if atom_ix == new_atm_idx[0]:
                resid = atom.resid
                resid_index_to_ref_allatoms_idx[res_ix] = [atom_ix]
                resid_index_to_resid[res_ix] = resid
                atom_ix += 1
            else:
                if atom.resid == resid:
                    resid_index_to_ref_allatoms_idx[res_ix] += [atom_ix]
                    resid_index_to_resid[res_ix] = resid
                    resid = atom.resid
                    atom_ix += 1
                else:
                    res_ix += 1
                    resid_index_to_ref_allatoms_idx[res_ix] = [atom_ix]
                    resid_index_to_resid[res_ix] = resid
                    atom_ix += 1
                    resid = atom.resid
        
        # Quality check that if Calpha is True there is 1-to-1 mapping of resid index to allatom indexs
        #print(f'\nresid_index_to_ref_allatoms_idx:\n{resid_index_to_ref_allatoms_idx}')
        if Calpha == True:
            for k,v in resid_index_to_ref_allatoms_idx.items():
                if len(v) != 1:
                    raise ValueError(f'When Calpha is specified there should only be one resid index for each all atom index: resid index {k} has {v}')
        #print(f'\nref_allatoms_idx_to_resid:\n{ref_allatoms_idx_to_resid}')
        #print(f'\nresid_index_to_resid:\n{resid_index_to_resid}')

        assert len(new_atm_idx) == np.concatenate(list(resid_index_to_ref_allatoms_idx.values())).size, f"Not enough atom indicies! {pdb_file}"

        # x y z cooridnates of chain
        coor = ref_allatoms_unique.atoms.positions[new_atm_idx]

        for resid_idx, all_atom_idx in resid_index_to_ref_allatoms_idx.items():

            resid = PDB_resids[resid_idx]

            check_coor = coor[all_atom_idx]

            structure_coor = ref_allatoms_unique.select_atoms(f"resid {resid} or resid {resid}A").positions
            
            #assert np.all(check_coor == structure_coor), f"Coordinates do not match up! Resid {resid} {pdb_file}"
            try:
                np.all(check_coor == structure_coor)
            except:
                raise ValueError(f'Error in checking residue coordinates: most likely caused by resides with letters after them. Check resid: {resid} in PDB')

            if not np.all(check_coor == structure_coor):
                raise ValueError(f"Coordinates do not match up! Resid {resid} {pdb_file}")

        ca_coor = ref_ca_unique.positions

        resnames = ref_ca_unique.resnames

        chain_res = PDB_resids.size

        if PDB_resids.size:

            ent_result, missing_residues = get_entanglements(coor, chain_res, termini_threshold, pdb_file, PDB_resids, resnames, resid_index_to_ref_allatoms_idx, ca_coor, resid_index_to_resid)
            
            if ent_result: 
                if Calpha == True:
                    outfile = f'unmapped_GE_CA/{pdb}_GE.txt'
                else:
                    outfile = f'unmapped_GE/{pdb}_GE.txt'
                print(f'WRITING: {outfile}')
                for ij_gN_gC, crossings in ent_result.items():
                    if crossings.size:
                        with open(outfile, "a") as f:
                            f.write(f"Chain {chain} | ({ij_gN_gC[0]}, {ij_gN_gC[1]}, {crossings}) | {ij_gN_gC[2]} | {ij_gN_gC[3]}\n")
            else:
                print(f'NO ENT DETECTED for {pdb}')
            
            if len(missing_residues):
                print(f'WRITING: {pdb}_M.txt')
                with open(f"unmapped_missing_residues/{pdb}_M.txt", "a") as f:
                    f.write(f"Chain {chain}: ")
                    for m_res in missing_residues:
                        f.write(f"{m_res} ")
                    f.write("\n")
        

if __name__ == "__main__":

    import multiprocessing as mp
    import sys,os
    import argparse

    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--PDB", type=str, required=True, help="Path to PDB file you want to generate raw entanglments for")
    parser.add_argument("--GLN_threshold", type=float, required=False, help="Threshold applied to the absoluate value of the GLN to determine if an entanglement is present")
    parser.add_argument("--Calpha", type=str, required=False, help="use CA 8A cutoff instead of defualt 4.5A heavy atom cutoff for native contacts")
    parser.add_argument("--topoly_density", type=int, required=False, help="Density of the triangulation of minimal loop surface for determining pericing. Default=0 to speed up calculations but might cause unrealistic crossings in AF structures with large disorderd loops. Increase to 1 if that is the case")
    args = parser.parse_args()
    print(args)

    # parse some of the default parameters
    global pdb_dir, density, g_threshold, Calpha
    pdb_dir = args.PDB

    g_threshold = args.GLN_threshold
    if g_threshold == None:
        g_threshold = 0.6

    density = args.topoly_density
    if density == None:
        density = 0
    
    Calpha = args.Calpha
    if Calpha == None:
        Calpha = False
    else:
        if Calpha == 'True':
            Calpha = True
        else:
            Calpha = False
    
    print(f'pdb_dir: {pdb_dir}')
    print(f'g_threshold: {g_threshold}')
    print(f'density: {density}')
    print(f'Calpha: {Calpha}')

    cores = len(os.sched_getaffinity(0))

    result_obj = set()

    if Calpha == True:
        directories = {"unmapped_GE_CA", "unmapped_missing_residues", "Before_TER_PDBs", "clustered_unmapped_GE"}
    else:
        directories = {"unmapped_GE", "unmapped_missing_residues", "Before_TER_PDBs", "clustered_unmapped_GE"}

    folder_exists = set(os.listdir("."))

    for folder in directories:
        if folder not in folder_exists:
            os.mkdir(f"{os.getcwd()}/{folder}") 

    if os.path.isdir(pdb_dir):
        input_data = [f'{pdb_dir}{f}' for f in os.listdir(f"{pdb_dir}")]
    elif os.path.isfile(pdb_dir):
        input_data = [pdb_dir]

    for f in input_data:
        results = calculate_native_entanglements(f)


    print(f'NORMAL TERMINATION')

  
