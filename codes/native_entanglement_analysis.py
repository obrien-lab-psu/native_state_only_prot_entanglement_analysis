#!/usr/bin/env python3
################################################################################################################
script_title='native_entanglement_analysis'
script_author='Ian Sitarik'
script_updatelog=f"""Update log for {script_title}

                   Date: 08.20.2021
                   Note: Started covertion of topology_anal codes

                   Date: 09.01.2021
                   Note: v1.5 reproduced paper. This version will include a separate feature to iterate through all loops

                   Date: 09.11.2021
                   note: implimented joblib for multiprocessing support ~O(|native contacts|)

                   Date: 10.06.2021
                   note: added ability to track resid from PDB for entanglement

                   Date: 10.24.2021
                   note:- calculates two change in entanglement methods:
                        1. if there is a change in the un rounded partial linking number of any tail for a given
                        NC greater than the threshold provided by the user
                        2. if the total linking number changes

                        - added the ability to catch missing residues in two ways.
                        1. if there are any missing residues between the first and last residue in the PDB
                        2. from the PDB file its self as they should be reported unless the user removed the info
                            as happens sometimes

                        NOTE::: make sure reference and traj have same resid labeling

                   Date: 10.28.2021
                   note: added automatic restart functionaily

                   Date: 10.30.2021
                   note: added control file read in

                   Date: 10.30.2021
                   note: added fraction of native contacts analysis (automatic)

                   Date: 11.17.2021
                   note: added minimal loop finder

                   Date: 11.22.2021
                   note: added mask to select certain part of PDB

                   Date: 11.24.2021
                   note: added correction to output loop_analysis for ref state as well
                   note: added unique entanglement finder

                   Date: 02.08.2022
                   note: REBIRTH of script from OP_v8.0.py
                   note: overhaul to simplify inputs and calculations

                   Date: 03.05.2022
                   note: simplified output to two files
                   1. the time series for fraction of native contacts (Q) and the fraction of native contacts with a change in entanglement (G)
                   2. a pickle binary file containing a single dictionary with one entry for the reference state and a single entry for each frame analyzed in the trajectory
                        this file can be examined by launching an interactive python session and using the following commands

                        import pickle
                        with open('./test_outfiles/entanglement_analysis_1.4/output/6u32.pkl', 'rb') as fh:
                            data = pickle.load(fh)

                        for k,v in data.items():
                            print(k,v)


                        the top level of keys will be integers ranging from -1 to the number of frames analyzed minus 1. For eample if you
                        analyzed a trajectory with 10 frames the dictionary would have a total of 11 entries with the following keys

                        -1 = reference state results
                        0 = first frame results
                        1 = second frame results
                        ...
                        9 = tenth frame results


                        in each of these entires the value is another dictionary containing one entry for each native contact that was detected to have a change in entanglement

                        for the refernce state key = -1 the inner dictionary will be structured like this
                        key = (residues invovled in native contact)
                        value = [array containing partial linking values for the N and C termini for this native contact,
                                 array containing the N and C terminal crossings for this native contact,
                                 residues within 8A of the crossing residues]

                        so for example if the key value pair in the reference state returned this
                        (4, 50) [array([0.        , 0.84160559]), [[], [61]], [[], [24, 25, 26, 27, 28, 29, 34, 35, 36, 37, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 62, 63]]]
                        the native contact formed by the residues 4 and 50 has an entanglement present.
                        the partial linking value of the N terminus is 0 indicating no entanglement
                        the partial linking value of the C terminus is 0.84160559 indicating a entanglement is present
                        Therefor the N terminus should not have a crossing while the C terminus should and indeed does at residue 61
                        The residues who have alpha carbons within 8A of the crossing residues are reported last


                        for the frames in the traj key >= 0 the inner dictionary will be structured like this
                        key = (residues invovled in native contact)
                        value = [array containing partial linking values for the N and C termini for this native contact,
                                 array containing [change type, refernce state partial linking value, frame partial linking value] for the N and C termini for this native contact,
                                 array containing the N and C terminal crossings for this native contact,
                                 residues within 8A of the crossing residues]

                                change_type = 0 if gain of entanglement and no change in chirality
                                change_type = 1 if gain of entanglement and change in chirality
                                change_type = 2 if loss of entanglement and no change in chirality
                                change_type = 3 if loss of entanglement and change in chirality
                                change_type = 4 if only chirality

                        so for example if the key value pair in a frame returned this
                        (108, 135) [array([0.96495744, 0.00888425]), [[0, 0.06142746756886508, 0.9649574387268217], []], [[12, 106]], [[10, 11, 13, 14, 104, 105, 107, 108, 109, 110, 111, 112, 113, 114, 131, 132, 133, 134, 135, 163, 164, 167]]]
                        the native contact formed by the residues 108 and 135 has an entanglement present.
                        the partial linking value of the N terminus is 0.96495744 indicating a entanglement
                        the partial linking value of the C terminus is 0.00888425 indicating no entanglement is present
                        This signifies a gain of entanglement with no change in chirality so the change_type = 0
                        and the patial linking value in the reference state for the N termini is 0.06142746756886508 while in the misfolded state is 0.9649574387268217
                        Therefor the C terminus should not have a crossing while the N terminus should and indeed does at residue 12 and 106
                        The residues who have alpha carbons within 8A of the crossing residues are reported last


                   Date: 06.15.2022
                   note: added ability to detect missing residues in reference state files and correct the entanglement output to remove entanglements if
                         1. there are any missing residues within +/- 10 of any crossing
                         2. more than three consecutive missing residues in the loop
                         3. more than 5% of the loop has missing residues

                         This should not affect the trajectory caclulations as those structure should be complete and a complete reference state should be used when analyzing them.
                         This only affects the case when you are analyzing a single structure for entanglements and not comparing to a treatment structure or trajectory.

                    Date: 08.16.2022
                    note: changed from finding surr residues using mdanalysis to manual calculation

                    Date: 08.19.2022
                    note: fixed error in fraction of native contact analysis that did not include 1.2 cutoff for thermal flux

                    Date: 11.21.2022
                    note: added ability to skip topoly if you do not need crossings or surrounding residues and only need GQ output

                    Date: 4.19.2023
                    note: converted to be only native state analysis and no traj analysis. Also simplified output

                    Date: 10.7.2023
                    note: added csv output for easy user reading
                  """

################################################################################################################

import os
import sys
import numpy as np
import time
import itertools
from MDAnalysis import *
from scipy.spatial.distance import pdist, squareform
from itertools import product, combinations
from itertools import chain as iterchain
from joblib import Parallel, delayed
import configparser
import pickle
from topoly import *
import more_itertools as mit
import pandas as pd

##################################################################################################################
### START argument parse ###
##################################################################################################################
usage = """
[1] = num processors
[2] = ref coor file
[3] = ref atom selection mask
[4] = out_path"""

if len(sys.argv) != 5:
    print(f'{usage}')
    quit()


global nproc
nproc = int(sys.argv[1])
print(f'nproc: {nproc}')

ref_path= sys.argv[2]
print(f'ref_path: {ref_path}')
if ref_path == None: print(f'\n{script_updatelog}\n'); sys.exit()

ref_mask = sys.argv[3]
print(f'ref_mask: {ref_mask}')
if ref_mask == None: ref_mask = 'all'

out_path=sys.argv[4]
print(f'out_path: {out_path}')
if out_path == None: print(f'\n{script_updatelog}\n{usage}'); sys.exit()


global S
S = 5
print(f'S: {S}')

global TS
TS = 5
print(f'TS: {TS}')

global ent_gval_threshold
ent_gval_threshold = 0.6
print(f'ent_gval_threshold: {ent_gval_threshold}')

outfile_basename=ref_path.split('.pdb')[0]
print(f'outfile_basename: {outfile_basename}')
if outfile_basename == None or outfile_basename == '': print(f'\n{script_updatelog}\n'); sys.exit()



##################################################################################################################
### START initial loading of structure files and qualtiy control ###
##################################################################################################################

### START dir declaration ###

if os.path.exists(f'{out_path}/'):
    print(f'{out_path}/ does exists and will be used')
    pass
else:
    os.makedirs(f'{out_path}/')

if os.path.exists(f'{out_path}{script_title}/'):
    print(f'{out_path}{script_title}/ does exists and will be used')
    pass
else:
    os.makedirs(f'{out_path}{script_title}/')

if os.path.exists(f'{out_path}{script_title}/output/'):
    print(f'{out_path}{script_title}/output/ does exists and will be used')
    pass
else:
    os.makedirs(f'{out_path}{script_title}/output/')

### END dir declaration ###

### START preference setting ###

start_time=time.time() #time since epoch
print('time since epoch = '+str(start_time))

np.set_printoptions(precision=4, suppress=True, linewidth=sys.maxsize, threshold=sys.maxsize)
np.seterr(divide='ignore', invalid='ignore')

### END preference setting ###

######################################################################################################################
# USER DEFINED FUNCTIONS                                                                                             #
######################################################################################################################

def gen_nc_gdict(coor, coor_cmap, **kwargs):
    dom_nc_gdict = {}
    dom_gn_dict = {}
    dom_contact_ent = {}
    global dot_matrix, l

    nc_indexs = np.stack(np.nonzero(coor_cmap)).transpose()

    l = len(coor)
    #print(f'l: {l}')

    #make R and dR waves of length N-1
    range_l = np.arange(0, l-1)
    range_next_l = np.arange(1,l)

    coor = coor.astype(np.float32)
    R = 0.5*(coor[range_l] + coor[range_next_l])
    dR = coor[range_next_l] - coor[range_l]

    #make dRcross matrix
    pair_array = np.asarray(list(itertools.product(dR,dR)))

    x = pair_array[:,0,:]
    y = pair_array[:,1,:]

    dR_cross = np.cross(x,y)

    #make Rnorm matrix
    pair_array = np.asarray(list(itertools.product(R,R)))

    diff = pair_array[:,0,:] - pair_array[:,1,:]
    diff = diff.astype(np.float32)
    Runit = diff / np.linalg.norm(diff, axis=1)[:,None]**3
    Runit = Runit.astype(np.float32)

    #make final dot matrix
    dot_matrix = [np.dot(x,y) for x,y in zip(Runit,dR_cross)]
    dot_matrix = np.asarray(dot_matrix)
    dot_matrix = dot_matrix.reshape((l-1,l-1))

    coor = [list(x) for x in coor]
    contact_ent = list(Parallel(n_jobs=nproc)(delayed(g_calc)(i, j, coor) for i,j in nc_indexs if j >= i + 10))
    contact_ent = {k: v for d in contact_ent for k, v in d.items()}

    return (contact_ent)


def g_calc(i,j, coor):
    loop_range = np.arange(i,j)
    crossings = [list(), list()]

    nterm_range = np.arange(TS,i-S)
    gn_parital_link, gn, gn_j1, gn_j2 = sub_lists(nterm_range,loop_range)
    if abs(gn_parital_link) > ent_gval_threshold:
        Ncrossings = lasso_type(coor, [(i,j)], more_info=True)
        crossings[0] = list(np.abs(np.unique(np.asarray([x.replace('*','') for x in Ncrossings[(i,j)]['crossingsN']], dtype=int))))

    cterm_range = np.arange(j+S,l-TS-1)
    gc_parital_link, gc, gc_j1, gc_j2 = sub_lists(cterm_range,loop_range)
    if abs(gc_parital_link) > ent_gval_threshold:

        Ccrossings = lasso_type(coor, [(i,j)], more_info=True)
        crossings[1] = list(np.abs(np.unique(np.asarray([x.replace('*','') for x in Ccrossings[(i,j)]['crossingsC']], dtype=int))))

    out = {(i, j):np.asarray([[gn_parital_link,gc_parital_link],crossings], dtype='O')}

    return out



#@nb.njit(fastmath=True)
def helper_func(g_vals: np.ndarray):

    return abs(g_vals.sum()/(4.0*np.pi))


def sub_lists(thread, loop):

    if len(thread) > 10:
        pairs = np.fromiter(itertools.chain(*itertools.product(thread, loop)), int).reshape(-1, 2)
        parital_link = dot_matrix[pairs[:,0], pairs[:,1]].sum()/(4.0*np.pi)

        return parital_link, parital_link, thread[0], thread[-1]

    else:
        return 0, 0, 0, 0


def ent_cmap(cor, ref = True, restricted = True, cut_off = 8.0, bb_buffer = 4, **kwargs):

    distance_map=squareform(pdist(cor,'euclidean'))
    distance_map=np.triu(distance_map,k=bb_buffer)

    contact_map = np.where((distance_map < cut_off) & (distance_map > 0), 1, 0)

    contact_num=contact_map.sum()

    return contact_map, contact_num


######################################################################################################################
# MAIN                                                                                                               #
######################################################################################################################

outdata = {}
#load ref structs and get entanglement
print('\n######################################## START analysis ########################################\n')

#load data into a universe and use the mask the suer provided
print(f'Loading: {ref_path}')
ref = Universe(ref_path)
ref_calphas = ref.select_atoms(f'{ref_mask}')
print('REF ATOMS: ', ref_calphas, len(ref_calphas))

#get mapping for coor_idx to PDB resid
global ref_cooridx2pdbresid
ref_cooridx2pdbresid = {i:res.resid for i,res in enumerate(ref_calphas.residues)}
ref_pdbresid2cooridx = {v: k for k, v in ref_cooridx2pdbresid.items()}
#print(ref_cooridx2pdbresid)
#for corid, resid in ref_cooridx2pdbresid.items():
#    print(corid, resid)

#get coordinate positions
ref_coor = ref_calphas.positions
ref_distance_map=squareform(pdist(ref_coor,'euclidean'))
#print(ref_coor[:10])


#get cmap for G and restricted cmap for Q if nec_elems file was specified
ref_cmap, ref_num_contacts = ent_cmap(ref_coor)
ref_distance_map=squareform(pdist(ref_coor,'euclidean'))

#GLN analysis
ref_cont_ent_data = gen_nc_gdict(ref_coor, ref_cmap)
for nc,(gvals, crossings) in ref_cont_ent_data.copy().items():

    nc = tuple(ref_cooridx2pdbresid[n] for n in nc)

    #get residues surrounding the crossings in the N terminus
    Ncrossings = crossings[0]
    if len(Ncrossings) != 0:
        distance_map = np.where(ref_distance_map < 8.0, 1, 0)[Ncrossings, :]
        Nsurr = list(np.where(distance_map == 1)[1])
        Nsurr = [ref_cooridx2pdbresid[s] for s in Nsurr]
    else:
        Nsurr = []

    #get residues surrounding the crossings in the C terminus
    Ccrossings = crossings[1]
    if len(Ccrossings) != 0:
        distance_map = np.where(ref_distance_map < 8.0, 1, 0)[Ccrossings, :]
        Csurr = list(np.where(distance_map == 1)[1])
        Csurr = [ref_cooridx2pdbresid[s] for s in Csurr]
    else:
        Csurr = []

    #map residue idx to resnum from PDB
    Ncrossings = [ref_cooridx2pdbresid[c] for c in Ncrossings]
    Ccrossings = [ref_cooridx2pdbresid[c] for c in Ccrossings]
    #print(nc,gvals, crossings)
    outdata[nc] = {'gval_N': gvals[0], 'gval_C': gvals[1], 'Ncrossings': Ncrossings, 'Ccrossings': Ccrossings, 'Nsurr': Nsurr, 'Csurr': Csurr}


for nc, ent in outdata.items():
    print(nc, ent)

print('\n######################################## END analysis ########################################\n')
### END loading of analysis universe ###


### output
outfilename = f'{out_path}{script_title}/output/{outfile_basename}.pkl'
with open(outfilename, "wb") as fh:
    pickle.dump(outdata, fh)

print(f'Saved: {outfilename}')


### make a dataframe to output a simple csv file as well for user readin
df = {'nc': [], 'gval_N': [], 'gval_C': [], 'Ncrossings': [], 'Ccrossings': [], 'Nsurr': [], 'Csurr': []}
for nc, ent_info in outdata.items():
    df['nc'] += [nc]
    df['gval_N'] += [ent_info['gval_N']]
    df['gval_C'] += [ent_info['gval_C']]
    df['Ncrossings'] += [ent_info['Ncrossings']]
    df['Ccrossings'] += [ent_info['Ccrossings']]
    df['Nsurr'] += [ent_info['Nsurr']]
    df['Csurr'] += [ent_info['Csurr']]
df = pd.DataFrame(df)
print(df)
outfilename = f'{out_path}{script_title}/output/{outfile_basename}.csv'
df.to_csv(outfilename, index=False)
print(f'SAVED: {outfilename}')

######################################################################################################################
comp_time=time.time()-start_time
print(f'computation time: {comp_time}')
print(f'NORMAL TERMINATION')
