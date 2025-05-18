#!/usr/bin/env python3
from collections import defaultdict
import numpy as np
import itertools
from geom_median.numpy import compute_geometric_median
from scipy.spatial.distance import cdist, squareform
from functools import cache
import re
import random
import pandas as pd
import pickle
import sys, getopt, math, os, time, traceback, glob, copy
from scipy.cluster.hierarchy import fcluster, linkage, cophenet
import parmed as pmd
import mdtraj as mdt
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.colors as mcolors
from scipy.stats import mode
import pathlib
from dataclasses import dataclass

matplotlib.use('Agg')
pd.set_option('display.max_rows', 500)

class ClusterNativeEntanglements:
    """
    Class to calculate native entanglements given either a file path to an entanglement file or an entanglement object
    """

    ##########################################################################################################################################################
    def __init__(self, organism: str = 'Ecoli', cut_off: int = 57) -> None:
        """
        Constructor for GaussianEntanglement class.

        Parameters
        ----------
        """

        if organism == 'Human':
            self.cut_off = 52
        elif organism == 'Ecoli':
            self.cut_off = 57
        elif organism == 'Yeast':
            self.cut_off = 49
        else:
            self.cut_off = cut_off
        self.organism = organism
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def loop_distance(self, entangled_A: tuple, entangled_B: tuple):

        # remove chiralites then perform euclidean distance
        new_cr_A = [int(cr_A[1:]) for cr_A in entangled_A[3:-1]]
        new_entangled_A = (entangled_A[0], entangled_A[1], entangled_A[2], *new_cr_A)

        new_cr_B = [int(cr_B[1:]) for cr_B in entangled_B[3:-1]]
        new_entangled_B = (entangled_B[0], entangled_B[1], entangled_B[2], *new_cr_B)

        return math.dist(new_entangled_A[1:], new_entangled_B[1:])
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    def check_step_ij_kl_range(self, ent1: tuple, ent2: tuple):

        # check if i or j of (i,j) reside within the range (inclusive) of (k,l), and vice versa

        nc_pair_1 = ent1[1:3]
        nc_pair_1_range = np.arange(ent1[1:3][0], ent1[1:3][1] + 1)

        nc_pair_2 = ent2[1:3]
        nc_pair_2_range = np.arange(ent2[1:3][0], ent2[1:3][1] + 1)

        #return True if (nc_pair_1[0] in nc_pair_2_range or nc_pair_1[1] in nc_pair_2_range or 

        if nc_pair_1[0] in nc_pair_2_range or nc_pair_1[1] in nc_pair_2_range:
            return True
        elif nc_pair_2[0] in nc_pair_1_range or nc_pair_2[1] in nc_pair_1_range:
            return True
        else:
            return False
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    #@cache
    def Cluster_NativeEntanglements(self, GE_filepath: str, outdir: str='./', outfile: str='Cluster_NativeEntanglements.txt'):
        
        """
        PARAMS:
            GE_file: str 
            cut_off: int

        1. Identify all unique "residue crossing set and chiralites"

            1b. sort the residues along with the chiralities
        
        2. Find the minimal loop encompassing a given "residue crossing set and chiralites"
            
            2b
            
            i. Identify entanglements that have any crossing residues between them that are 
                less than or equal to 3 residues apart and have the same chirality.
        
            ii. Then check if i or j of (i,j) reside within the range (inclusive) of (k,l), or vice versa;
            
            iii. If yes, then check if any crossing residues are in the range of min(i,j,k,l) to max(i,j,k,l); 
            if yes, skip rest of 2
            
            iv. If no, then check if the number of crossing residues, in each residue set, are different;
            
            v. All crossing residue(s) in the entanglement with the fewer crossings need to have a distance <= 20 
            with the crossings in the other entanglement. Do this by the "brute force" approach and 
            the true distance formula. This means, calculate the distances and take the minimal distance 
            as the distance you check that is less than or equal to 20.
            
            If yes, then keep the {i,j} {r} that have the greatest number of crossing residues;
            If not, then keep the two entanglements separate. 
        
        3. For at least two entanglements each with 1 or more crossings. 
            Loop over the entanglments two at time (avoid double counting)
            Check if i or j of (i,j) reside within the range (inclusive) of (k,l), or  vice versa;
            If yes, check if number of crossing residues is the same (and it is 1 or more)
            If yes, calculate the distances between all crossing residues
                and if both have the same chiralities. 
                (Do NOT use brute force, just compare 1-to-1 index of crossing residues).
            If all the distances are less than or equal to 20, then determine which 
                entanglement has the smaller loop, remove the entanglement with the larger loop
                
        4. Spatially cluster those outputs that have (i) the same number of crossings and (ii) the same chiralities

        """
        print(f'Clustering {self.organism} Native Entanglements with dist_cutoff: {self.cut_off}')
        GE_file = GE_filepath.split('/')[-1]
        print(GE_file, self.cut_off, outdir)

        full_entanglement_data = defaultdict(list)

        ent_data = defaultdict(list)

        rep_ID_ent = defaultdict(list) 

        grouped_entanglement_data = defaultdict(list)

        Before_cr_dist = defaultdict(list)
        
        After_cr_dist = defaultdict(list)

        entanglement_partial_g_data = {}
        
        ## Check if the clustering file is already made and if so use it
        outfilepath = os.path.join(f'{outdir}', f'{outfile}')
        if os.path.exists(outfilepath):
            print(f'{outfilepath} ALREADY EXISTS AND WILL BE LOADED')
            outdf = pd.read_csv(outfilepath, sep='|')
            return {'outfile':outfilepath, 'ent_result':outdf}

        print(f'Loading {GE_filepath}')
        column_names = ['ID', 'ijr', 'gn', 'gc', 'CCbond']
        GE_data = pd.read_csv(GE_filepath, sep='|', dtype={'c': str}, names=column_names)
        print(GE_data)
        #GE_data = GE_data[GE_data['ENT'] == True].reset_index(drop=True)
        # if Quality is a column name then only get the High Quality raw entanglements
        #if 'Quality' in GE_data.keys():
        #    GE_data = GE_data[GE_data['Quality'] == 'High'].reset_index(drop=True)
        #print(GE_data)

        ### STEP 1 INITAL LOADING AND MERGING ################################################################################################################
        ############################################################################################
        ## parse the entanglement file and
        ## get those native contacts that are disulfide bonds
        print(f'# Step 1')
        CCBonds = []
        num_raw_ents = {}
        for rowi, row in GE_data.iterrows():
            ID = row['ID']
            #print(row)
            ijr = row['ijr']
            #print(ijr)
            native_contact_i, native_contact_j, crossing_res = ijr.split(', ')
            native_contact_i = native_contact_i.replace('(', '')
            crossing_res = crossing_res.replace(']) ', '') 
            crossing_res = crossing_res.replace('[', '')
            crossing_res = crossing_res.replace("'", '')

            native_contact_i, native_contact_j = int(native_contact_i.strip(' ')), int(native_contact_j.strip(' '))
            #print(native_contact_i, native_contact_j, crossing_res)
     
            #native_contact_i, native_contact_j, crossing_res = row['i'], row['j'], row['c']
            gn, gc = row['gn'], row['gc']
            CCbond = row['CCbond']

            # keep track of number of raw ents for QC purposes
            if ID not in num_raw_ents:
                num_raw_ents[ID] = 1
            else:
                num_raw_ents[ID] += 1

            #native_contact_i, native_contact_j, crossing_res = line[1], line[2], line[3]
            #native_contact_i = int(native_contact_i)
            #native_contact_j = int(native_contact_j)

            reformat_cr = crossing_res.split(' ')
            #print(reformat_cr)

            reformat_cr = sorted(reformat_cr, key = lambda x: int(re.split("\\+|-", x, maxsplit= 1)[1]))
            #print(native_contact_i, native_contact_j, reformat_cr)


            # Step 1 and 1b
            grouped_entanglement_data[(ID, *reformat_cr)].append((native_contact_i, native_contact_j))

            entanglement_partial_g_data[(native_contact_i, native_contact_j, *reformat_cr)] = (gn, gc)

            #print(f'CCbond: {CCbond}')
            if CCbond == True:
                CCBonds += [(native_contact_i, native_contact_j)]

        #print(f'num_raw_ents: {num_raw_ents}')        

        #print(f'Step 1 results')
        Step1_QC_counter = 0
        for k,v in grouped_entanglement_data.items():
            print(k,v)
            Step1_QC_counter += len(v)
        print(f'CCBonds: {CCBonds}')
    
        ### STEP 2 ################################################################################################################
        ############################################################################################
        # Step 2 Get the minimal loop encompassing each set of unique crossings
        print(f'\n# Step 2a')
        for ID_cr, loops in grouped_entanglement_data.items():

            ID = ID_cr[0]

            crossings = np.asarray(list(ID_cr[1:]))

            loop_lengths = [nc[1] - nc[0] for nc in loops]

            minimum_loop_length = min(loop_lengths)

            minimum_loop_length_index = loop_lengths.index(minimum_loop_length)

            minimum_loop_nc_i, minimum_loop_nc_j = loops[minimum_loop_length_index]

            ent_data[ID].append((len(loops), minimum_loop_nc_i, minimum_loop_nc_j, *crossings, ';'.join(['-'.join([str(loop[0]), str(loop[1])]) for loop in loops])))

        #print(f'Step 2a results')
        for ID, ents in ent_data.items():
            Step2a_QC_counter = 0
            for ent_i, ent in enumerate(ents):
                #print(ID, ent_i, ent)
                Step2a_QC_counter += ent[0]

            ## QC that the number of tracked entanglements after step 2a is still valid
            #print(f'Step2a_QC_counter: {Step2a_QC_counter} should = {num_raw_ents[ID]}')
            if Step2a_QC_counter != num_raw_ents[ID]:
                raise ValueError(f'The number of tracked entaglements after Step 2a {Step2a_QC_counter} != {num_raw_ents[ID]}')

        ############################################################################################
        # Step 2b: 
        print(f'\n# Step 2b')
        merged_ents = []
        for ID, ents in ent_data.items():
            orig_ents = ents.copy()
            comb_ents = itertools.combinations(ents, 2)

            # for each pair of ents
            for each_ent_pair in comb_ents:
                #print(f'\nAnalyzing pair: {each_ent_pair}')
                if each_ent_pair[0] == each_ent_pair[1]:
                    print(f'Ents are the same: {each_ent_pair}')
                    continue

                distance_thresholds = []

                ent1 = each_ent_pair[0]
                ent2 = each_ent_pair[1]
                #print(f'\n{"#"*100}\nent1: {ent1} | ent2: {ent2}')

                # get crossings from ent pair without chiralities
                cr1 = set([int(ent_cr_1[1:]) for ent_cr_1 in list(ent1[3:-1])])
                cr2 = set([int(ent_cr_2[1:]) for ent_cr_2 in list(ent2[3:-1])])
                #print(cr1, cr2)

                # get all possible pairs of the crossings
                all_cr_pairs = itertools.product(ent1[3:-1], ent2[3:-1])

                # get the distances between all pairs of crossings
                cr_dist_same_chiral = np.abs([int(pr[0][1:]) - int(pr[1][1:]) for pr in all_cr_pairs if pr[0][0] == pr[1][0]])
                #print(cr_dist_same_chiral)
                
                
                # if any of those distances is less than 3 and the number of crossings is not the same and ij in range of kl
                dist_check = np.any(cr_dist_same_chiral <= 3)
                loop_check = self.check_step_ij_kl_range(ent1, ent2)
                diff_cross_size_check = len(cr1) != len(cr2)
                #print(dist_check, loop_check, diff_cross_size_check)

                if np.any(cr_dist_same_chiral <= 3) and len(cr1) != len(cr2) and self.check_step_ij_kl_range(ent1, ent2):
                    #print(f'\nAnalyzing pair: {each_ent_pair}')
                    #print(f'step 2b conditions met')

                    minumum_loop_base = min(ent1[1], ent1[2], ent2[1], ent2[2])

                    maximum_loop_base = max(ent1[1], ent1[2], ent2[1], ent2[2])

                    all_crossings = cr1.union(cr2)

                    min_max_loop_base_range = set(range(minumum_loop_base, maximum_loop_base + 1))

                    # if the crossings are not within the min max loop range covering both entanglements
                    if not min_max_loop_base_range.intersection(all_crossings):

                        fewer_cr = min(cr1, cr2, key = len)
                        more_cr = max(cr1, cr2, key = len)

                        distributive_product = list(itertools.product(fewer_cr, more_cr))

                        slices = itertools.islice(distributive_product, 0, None, len(more_cr))
                    
                        groupings = []

                        for end_point in slices:

                            first_index = distributive_product.index(end_point)

                            groupings.append(distributive_product[first_index:len(more_cr) + first_index])

                        if len(groupings) != 1:

                            all_pair_products = itertools.product(*groupings)
                            all_pair_groupings = set()

                            for pairs in all_pair_products:

                                flag = True

                                # check common elements column wise
                                stacked_pairs = np.stack(pairs)

                                for col in range(stacked_pairs.shape[1]):

                                    if stacked_pairs[:, col].size != len(set(stacked_pairs[:, col])):

                                        flag = False
                                        break
                                
                                if flag:

                                    all_pair_groupings.add(pairs)

                        else: 

                            all_pair_groupings = groupings[0]

                        for condensed_pair in all_pair_groupings:

                            if isinstance(condensed_pair[0], int):

                                # when dealing with ent with one crossing 

                                condensed_pair = [condensed_pair]
                            
                            dist = np.sqrt(sum([(each_ele[0] - each_ele[1]) ** 2 for each_ele in condensed_pair]))

                            distance_thresholds.append(dist)

                        # all_pair_groupings and distance thresholds have the same size
                        if min(distance_thresholds) <= 20:
                            #print(f'ent1: {ent1} | ent2: {ent2}')
                            
                            min_ent = min(ent1, ent2, key = len)
                            max_ent = max(ent1, ent2, key = len)
                            if min_ent == max_ent:
                                print(f'WARNING: Ents are the same. Setting min_ent = ent1 and max_ent = ent2')
                                min_ent = ent1
                                max_ent = ent2
         
                            #print(f'ent1: {ent1} | ent2: {ent2}')
                            #print(f'min_ent: {min_ent}')
                            #print(f'max_ent: {max_ent}')
                         

                            if min_ent in ents and len(ents) > 1:
                            #if min_ent in ents and max_ent in ents and len(ents) > 1:
                                #min_ent_num_ncs = min_ent[0]
                    
                                #print(f'Removing: min_ent {min_ent} at index {ents.index(min_ent)}')
                                del ents[ents.index(min_ent)]
                                #del ents[ents.index(max_ent)]

                                
                                if max_ent == min_ent:
                                    raise ValueError(f'WARNING: Ents are the same\n{min_ent} == {max_ent}')
                                else:
                                    merged_ents += [(max_ent, min_ent)]

   
        print(f'\nStep 2b results')
        # results foor the end of step 2
        for ID, ents in ent_data.items():
            ent_dict = {ent_idx:[ent] for ent_idx, ent in enumerate(ents)}
            for ent_idx, ent in ent_dict.items():
                print(ent_idx, ent)


            ### Update entanglement list with those that got merged
            print(f'\n### Update entanglement dict with those that got merged\nNumber of merged pairs: {len(merged_ents)}')
            while len(merged_ents) != 0:
                pre_num_merged = len(merged_ents)
                #print(f'# merged_ents: {pre_num_merged}')
                for m_ent in merged_ents:
                    #print(f'\n{m_ent}')
                    for ent_idx, ent in ent_dict.copy().items():
                        #print(ent_idx, ent)
                        if m_ent[0] in ent:
                            #print(f'FOUND MATCH for kept ent {ent_idx}')
                            ent_dict[ent_idx] += [m_ent[1]]
                            merged_ents.remove(m_ent)

                #print(f'# merged_ents: {len(merged_ents)}')

                # QC to ensure you dont enter an infinite loop
                if len(merged_ents) == pre_num_merged:
                    raise ValueError('Failed to find a match this cycle and entering infi loop')
                    
            print(f'\n### merge ents to final and reform ent_data')
            updated_ents = []
            for ent_idx, ent in ent_dict.items():
                #print(ent_idx, ent)
                if len(ent) > 1:
                    num_loops = np.sum([e[0] for e in ent])
                    NCs = ';'.join([e[-1] for e in ent])
                    #print(ent, num_loops, NCs)
                    ent = (num_loops, *ent[0][1:-1], NCs)
                    updated_ents += [ent]
                else:
                    updated_ents += ent

            #print(f'Results after adding those that got merged to each representative entanglement')
            Step2b_QC_counter = 0
            for uent in updated_ents:
                #print(uent)
                Step2b_QC_counter += uent[0]

            ## QC that the number of tracked entanglements after step 2a is still valid
            print(f'Step2b_QC_counter: {Step2b_QC_counter} should = {num_raw_ents[ID]}')
            if Step2b_QC_counter != num_raw_ents[ID]:
                raise ValueError(f'The number of tracked entaglements after Step 2b {Step2b_QC_counter} != {num_raw_ents[ID]}')        

            ent_data[ID] = updated_ents
        for ID, ents in ent_data.items():
            print(ID)
            for ent_i, ent in enumerate(ents):
                print(ent_i, ent)
        
        ### STEP 3 ################################################################################################################
        # Step 3
        print(f'\n# Step 3')
        for ID, processed_ents in ent_data.items():

            comb_processed_ents = itertools.combinations(processed_ents, 2)

            keep_track_of_larger_proc_ent = []
            keep_track_of_shorter_proc_ent = []

            for each_processed_ent_pair in comb_processed_ents:
                #print(f'\npair of ents: {each_processed_ent_pair}')

                proc_ent1 = each_processed_ent_pair[0]
                proc_ent2 = each_processed_ent_pair[1]

                proc_ent1_ijr = proc_ent1[1:-1]
                proc_ent2_ijr = proc_ent2[1:-1]
                #print(proc_ent1_ijr, proc_ent2_ijr)

                if proc_ent1_ijr not in keep_track_of_larger_proc_ent and proc_ent2_ijr not in keep_track_of_larger_proc_ent:

                    # without chiralites
                    proc_cr1 = np.asarray([int(ent_cr_1[1:]) for ent_cr_1 in list(proc_ent1[3:-1])])
                    proc_cr2 = np.asarray([int(ent_cr_2[1:]) for ent_cr_2 in list(proc_ent2[3:-1])])

                    if len(proc_ent1[3:-1]) == len(proc_ent2[3:-1]):

                        chirality1 = [chir1[0] for chir1 in proc_ent1[3:-1]]
                        chirality2 = [chir2[0] for chir2 in proc_ent2[3:-1]]

                        if chirality1 == chirality2 and self.check_step_ij_kl_range(proc_ent1, proc_ent2) and np.all(np.abs(proc_cr1 - proc_cr2) <= 20):
                            #print(proc_ent1, proc_ent2)

                            loop1_length = proc_ent1[2] - proc_ent1[1]
                            loop2_length = proc_ent2[2] - proc_ent2[1]
                            #print(loop1_length, loop2_length)

                            check = [loop1_length, loop2_length]

                            maximum_loop_length = max(loop1_length, loop2_length)
                            minimum_loop_length = min(loop1_length, loop2_length)

                            if maximum_loop_length == minimum_loop_length:
                                longer_loop_ent = proc_ent1
                                shorter_loop_ent = proc_ent2
                            else:
                                longer_loop_ent = each_processed_ent_pair[check.index(maximum_loop_length)]
                                shorter_loop_ent = each_processed_ent_pair[check.index(minimum_loop_length)]
                            longer_loop_ent_ijr = longer_loop_ent[1:-1]
                            shorter_loop_ent_ijr = shorter_loop_ent[1:-1]

                            if len(processed_ents) > 1:

                                for long_proc_ent_index, long_proc_ent in enumerate(processed_ents):
                                    long_proc_ent_ijr = long_proc_ent[1:-1]
                                    #print(long_proc_ent_index, long_proc_ent, long_proc_ent_ijr)
                                    if long_proc_ent_ijr == longer_loop_ent_ijr:
                                        break
                                del processed_ents[long_proc_ent_index]
                                # remove the one with larger loop

                                # find the shorter loop and remove it
                                for short_proc_ent_index, short_proc_ent in enumerate(processed_ents):
                                    short_proc_ent_ijr = short_proc_ent[1:-1]
                                    if short_proc_ent_ijr == shorter_loop_ent_ijr:
                                        break
                                del processed_ents[short_proc_ent_index]

                                updated_ent = (short_proc_ent[0] + long_proc_ent[0],  *short_proc_ent[1:-1], short_proc_ent[-1] + ';' + long_proc_ent[-1])
                                processed_ents += [updated_ent]

                            keep_track_of_larger_proc_ent.append(longer_loop_ent_ijr)
                            keep_track_of_shorter_proc_ent.append(shorter_loop_ent_ijr)
        
        #print(f'Step 3 results')
        for ID, ents in ent_data.items():
            Step3_QC_counter = 0
            for ent in ents:
                #print(ent)
                Step3_QC_counter += ent[0]
            
            # QC to ensure number of raw ents was preserved after step 3
            print(f'Step3_QC_counter: {Step3_QC_counter} should = {num_raw_ents[ID]}')
            if Step3_QC_counter != num_raw_ents[ID]:
                raise ValueError(f'The number of tracked entaglements after Step 3 {Step3_QC_counter} != {num_raw_ents[ID]}')    

        ### STEP 4 SPATIAL CLUSTERING ################################################################################################################
        # Step 4 prep
        print(f'\n# Step 4 prep')
        for ID, new_ents in ent_data.items():

            for ent in new_ents:

                number_of_crossings = len(ent[3:-1])

                chiralites = [each_cr[0] for each_cr in ent[3:-1]]

                ID_num_chirality_key = f"{ID}_{number_of_crossings}_{chiralites}"

                full_entanglement_data[ID_num_chirality_key].append(ent)

        reset_counter = []

        # Step 4
        print(f'\n# Step 4')
        for ID_num_chiral in full_entanglement_data.keys():
            #print(ID_num_chiral)
            #ID = ID_num_chiral.split("_")[0]

            if ID not in reset_counter:

                reset_counter.append(ID)

                split_cluster_counter = 0

            length_key = defaultdict(list)
            loop_dist = defaultdict(list)
            dups = []
            clusters = {} 
            cluster_count = 0

            pairwise_entanglements = list(itertools.combinations(full_entanglement_data[ID_num_chiral], 2))

            if pairwise_entanglements:

                for i, pairwise_ent in enumerate(pairwise_entanglements):

                    dist = self.loop_distance(pairwise_ent[0], pairwise_ent[1])
                    
                    if dist <= self.cut_off and pairwise_ent[0] not in dups and pairwise_ent[1] not in dups:
                        # 1. pair must be <= self.cut_off
                        # 2. the neighbor cannot be the next key and it cannot be captured by another key

                        loop_dist[pairwise_ent[0]].append(pairwise_ent[1])
                        dups.append(pairwise_ent[1])
                
                key_list = list(loop_dist.keys())

                for key in key_list:

                    length_key[len(loop_dist[key])].append(key)

                # create clusters

                while len(length_key.values()) > 0:

                    max_neighbor = max(length_key.keys())

                    selected_ent = random.choice(length_key[max_neighbor])

                    cluster = copy.deepcopy(loop_dist[selected_ent])
                    cluster.append(selected_ent)

                    clusters[cluster_count] = cluster
                    cluster_count += 1

                    length_key[max_neighbor].remove(selected_ent)

                    if len(length_key[max_neighbor]) == 0:
                        length_key.pop(max_neighbor)
            
            # create single clusters
            if clusters:
                clusters_ijr_values = list(itertools.chain.from_iterable(list(clusters.values())))
            else:
                clusters_ijr_values = []

            full_ent_values = np.asarray(full_entanglement_data[ID_num_chiral], dtype=object)

            difference_ent = np.zeros(len(full_ent_values), dtype=bool)

            for k, ijr in enumerate(full_ent_values):

                if tuple(ijr) in clusters_ijr_values:
                    difference_ent[k] = True
                else:
                    difference_ent[k] = False

            i = np.unique(np.where(difference_ent == False)[0])

            next_cluster_count = cluster_count

            for single_cluster in full_ent_values[i]:
                
                single_cluster_list = []
                single_cluster_list.append(tuple(single_cluster))

                clusters[next_cluster_count] = single_cluster_list

                next_cluster_count += 1

            # pick representative entanglement per cluster
            #print(f'\nPick rep entanglements')
            for counter, ijr_values in clusters.items():
                #print(f'\nCluster {counter} {ijr_values}')

                # clusters contain many entanglements
                if len(ijr_values) > 1:

                    ijr = np.asarray(ijr_values)
                    #print(f'cluster ijr:\n{ijr}')

                    cr_values = np.asarray([[int(r_value[0][1:])] for r_value in ijr[:, 3:-1]])
                    #print(f'cr_values: {cr_values}')

                    median_cr = compute_geometric_median(cr_values).median
                    #print(f'median_cr: {median_cr}')

                    distances = cdist(cr_values, [median_cr])

                    minimum_distances_i = np.where(distances == min(distances))[0]
                    #print(f'minimum_distances_i: {minimum_distances_i}')

                    possible_cand = ijr[minimum_distances_i]
                    #print(f'possible_cand:\n{possible_cand}')

                    loop_lengths = np.abs(possible_cand[:, 1].astype(int) - possible_cand[:, 2].astype(int))
                    #print(f'loop_lengths: {loop_lengths}')

                    smallest_loop_length = min(loop_lengths)
                    #print(f'smallest_loop_length: {smallest_loop_length}')

                    num_nc_in_cluster = np.sum([int(n[0]) for n in ijr_values])
                    #print(f'num_nc_in_cluster: {num_nc_in_cluster}')

                    loops = ';'.join([n[-1] for n in ijr_values])
                    #print(f'loops: {loops}')
                    
                    #rep_entanglement = possible_cand[random.choice(np.where(smallest_loop_length == loop_lengths)[0])]
                    rep_entanglement = possible_cand[random.choice(np.where(smallest_loop_length == loop_lengths))[0]]
                    rep_entanglement = [str(num_nc_in_cluster), *rep_entanglement[1:-1], loops]
                    #rep_ID_ent[f"{ID}_{split_cluster_counter}"].append(rep_entanglement)
                    rep_ID_ent[(ID, split_cluster_counter)].append(rep_entanglement)

                # clusters with a single entnalgement
                else:
                    #rep_ID_ent[f"{ID}_{split_cluster_counter}"].append(ijr_values[0])
                    rep_ID_ent[(ID, split_cluster_counter)].append(ijr_values[0])
                
                split_cluster_counter += 1
        
        ## QC Step 4 results
        print(f'Step 4 results')
        num_raw_ents_FINAL = {}
        for ID_counter, ijrs in rep_ID_ent.items():
            #print(ID_counter, ijrs)
            ID, counter = ID_counter
            #print(ID_counter, ID, counter, ijrs)

            if ID not in num_raw_ents_FINAL:
                num_raw_ents_FINAL[ID] = 0

            for ijr in ijrs:
                num_nc = int(ijr[0])
                num_raw_ents_FINAL[ID] += num_nc

        ## check the final tracking of raw ents
        for ID, count in num_raw_ents.items():
            if count != num_raw_ents_FINAL[ID]:
                raise ValueError(f'The FINAL # of raw ents {num_raw_ents_FINAL[ID]} != the starting {count} for ID {ID}')
            
        ### STEP 5 OUTPUT FILE ################################################################################################################
        # Step 5
        print(f'\nStep 5: Output file')

        ## set up the outdir for this calculation
        #outdir = f"{os.getcwd()}/{outdir}"
        if not os.path.isdir(outdir):
            os.mkdir(f"{outdir}") 
            print(f"Creating directory: {outdir}")

        outfilepath = os.path.join(f'{outdir}', f'{outfile}')
        #print(f'WRITING: {outfilepath}')

        with open(outfilepath, "w") as f:

            f.write(f'ID|i|j|c|gn|gc|num_contacts|contacts|CCBond\n')
            for ID_counter, ijrs in rep_ID_ent.items():

                ID, counter = ID_counter

                for ijr in ijrs:

                    new_ijr = (int(ijr[1]), int(ijr[2]), *list(ijr[3:-1]))

                    num_nc = int(ijr[0])

                    gn, gc = entanglement_partial_g_data[new_ijr]
                    gn = float(gn)
                    gc = float(gc)

                    ## check for disulfide bonds
                    CCBond_flag = False
                    for CCBond in CCBonds:
                        check1 = f'{CCBond[0]}-{CCBond[1]}' 
                        check2 = f'{CCBond[1]}-{CCBond[0]}'
                        if check1 in ijr[-1] or check2 in ijr[-1]:
                            CCBond_flag = True

                    #line = f"{ID}|{(int(ijr[1]), int(ijr[2]), *ijr[3:-1])}|{gn:.5f}|{gc:.5f}|{num_nc}|{ijr[-1]}|{CCBond_flag}"
                    line = f"{ID}|{int(ijr[1])}|{int(ijr[2])}|{','.join(ijr[3:-1])}|{gn:.5f}|{gc:.5f}|{num_nc}|{ijr[-1]}|{CCBond_flag}"
                    #print(line)
                    f.write(f"{line}\n")
        print(f'SAVED: {outfilepath}')
        outdf = pd.read_csv(outfilepath, sep='|')
        return {'outfile':outfilepath, 'ent_result':outdf}
    ##########################################################################################################################################################

    ##########################################################################################################################################################
    #def run(self,):

    ##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################


if __name__ == "__main__":

    import multiprocessing as mp
    import sys,os
    import argparse

    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--rawENT", type=str, required=True, help="Path to raw entanglement file")
    parser.add_argument("--name", type=str, required=True, help="An id for the PDB to be analyzed")
    parser.add_argument("--outdir", type=str, default='./Clustered_GE', help="Output directory for the clustering results")
    parser.add_argument("--organism", type=str, default='Ecoli', help="Organism to be used for the clustering, Ecoli or Human or Yeast")
    args = parser.parse_args()
    print(args)

    # parse some of the default parameters
    rawENT = args.rawENT
    name = args.name
    organism = args.organism
    outdir = args.outdir
    
    print(f'rawENT: {rawENT}')
    print(f'name: {name}')
    print(f'organism: {organism}')
    print(f'outdir: {outdir}')
    
    
    clustering = ClusterNativeEntanglements(organism=organism)
    print(clustering)

    ## Cluster the native entanglements
    nativeClusteredEnt = clustering.Cluster_NativeEntanglements(rawENT, outdir=outdir, outfile=f'{name}_Clustered_GE.txt')
    print(nativeClusteredEnt)

    print(f'NORMAL TERMINATION')