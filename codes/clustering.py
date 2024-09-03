from collections import defaultdict
import numpy as np
import itertools
from geom_median.numpy import compute_geometric_median
from scipy.spatial.distance import cdist
from functools import cache
import re
import math
import random
import copy
import os

def loop_distance(entangled_A: tuple, entangled_B: tuple):

    # remove chiralites then perform euclidean distance
    new_cr_A = [int(cr_A[1:]) for cr_A in entangled_A[3:-1]]
    new_entangled_A = (entangled_A[0], entangled_A[1], entangled_A[2], *new_cr_A)

    new_cr_B = [int(cr_B[1:]) for cr_B in entangled_B[3:-1]]
    new_entangled_B = (entangled_B[0], entangled_B[1], entangled_B[2], *new_cr_B)

    return math.dist(new_entangled_A[1:], new_entangled_B[1:])


def check_step_ij_kl_range(ent1: tuple, ent2: tuple):

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


@cache
def cluster_entanglements(GE_file_w_cutoff: tuple):

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

    GE_file = GE_file_w_cutoff[0].split('/')[-1]
    cut_off = GE_file_w_cutoff[1]
    filename_split = GE_file.split("_")
    
    #Add filename warning
    if len(filename_split) > 2:
        print('Warning: The unmapped_GE filenames should only contain a single "_" character: ',GE_file) 

    protein = filename_split[0]

    if f"{protein}_clustered_GE.txt" not in os.listdir(f"{outpath}"):
        print(f"CLUSTERING ENTANGLEMENTS FOR \033[4m{protein}\033[0m")
    
    else:
        print(f"\033[4m{protein}\033[0m IS CLUSTERED. PLEASE LOOK IN {outpath}")
        return

    full_entanglement_data = defaultdict(list)

    ent_data = defaultdict(list)

    rep_chain_ent = defaultdict(list) 

    grouped_entanglement_data = defaultdict(list)

    Before_cr_dist = defaultdict(list)
    
    After_cr_dist = defaultdict(list)

    entanglement_partial_g_data = {}
    
    print(f'Loading {GE_file_w_cutoff[0]}')
    GE_data = open(f"{GE_file_w_cutoff[0]}", "r").readlines()
    GE_data = np.asarray(GE_data)
    #if GE_data.size == 1:

    for line in GE_data:
        line = line.split("|")
        print(line)

        if len(line) == 4:

            chain = line[0].strip()

            native_contact_i, native_contact_j, crossing_res = line[1].split(",")

            native_contact_i = int(native_contact_i.replace("(", "").strip())

            native_contact_j = int(native_contact_j.strip())

            reformat_cr = crossing_res.replace(")", "").strip()[1:-1].replace("'", "").split()

            reformat_cr = sorted(reformat_cr, key = lambda x: int(re.split("\\+|-", x, maxsplit= 1)[1]))

            # Step 1 and 1b
            grouped_entanglement_data[(chain, *reformat_cr)].append((native_contact_i, native_contact_j))

            entanglement_partial_g_data[(native_contact_i, native_contact_j, *reformat_cr)] = (line[2].strip(), line[3].strip())
    print(f'Step 1 results')
    for k,v in grouped_entanglement_data.items():
        print(k,v)

    # Step 2
    print(f'\n# Step 2')
    for chain_cr, loops in grouped_entanglement_data.items():

        chain = chain_cr[0]

        crossings = np.asarray(list(chain_cr[1:]))

        loop_lengths = [nc[1] - nc[0] for nc in loops]

        minimum_loop_length = min(loop_lengths)

        minimum_loop_length_index = loop_lengths.index(minimum_loop_length)

        minimum_loop_nc_i, minimum_loop_nc_j = loops[minimum_loop_length_index]

        ent_data[chain].append((len(loops), minimum_loop_nc_i, minimum_loop_nc_j, *crossings, ';'.join(['-'.join([str(loop[0]), str(loop[1])]) for loop in loops])))
    print(f'Step 2 results')
    for chain, ents in ent_data.items():
        for ent_i, ent in enumerate(ents):
            print(ent_i, ent)

    # Step 2b
    print(f'\n# Step 2b')
    merged_ents = []
    for chain, ents in ent_data.items():
        comb_ents = itertools.combinations(ents, 2)

        # for each pair of ents
        for each_ent_pair in comb_ents:
            #print(f'\nAnalyzing pair: {each_ent_pair}')

            distance_thresholds = []

            ent1 = each_ent_pair[0]
            ent2 = each_ent_pair[1]

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
            if np.any(cr_dist_same_chiral <= 3) and len(cr1) != len(cr2) and check_step_ij_kl_range(ent1, ent2):
                #print(f'\nAnalyzing pair: {each_ent_pair}')
                #print(f'step 2b conditions met')

                minumum_loop_base = min(ent1[1], ent1[2], ent2[1], ent2[2])

                maximum_loop_base = max(ent1[1], ent1[2], ent2[1], ent2[2])

                all_crossings = cr1.union(cr2)

                min_max_loop_base_range = set(range(minumum_loop_base, maximum_loop_base + 1))

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

                        min_ent = min(ent1, ent2, key = len)
                        max_ent = max(ent1, ent2, key = len)
                        #print(f'min_ent: {min_ent}')
                        #print(f'max_ent: {max_ent}')

                        if min_ent in ents and len(ents) > 1:
                        #if min_ent in ents and max_ent in ents and len(ents) > 1:
                            #min_ent_num_ncs = min_ent[0]
                
                            #print(f'Removing: min_ent {min_ent} at index {ents.index(min_ent)}')
                            del ents[ents.index(min_ent)]
                            #del ents[ents.index(max_ent)]

                            merged_ents += [(max_ent, min_ent)]


    
    print(f'\nStep 2b results')
    # results foor the end of step 2
    for chain, ents in ent_data.items():
        ent_dict = {ent_idx:[ent] for ent_idx, ent in enumerate(ents)}
        for ent_idx, ent in ent_dict.items():
            print(ent_idx, ent)

        ### Update entanglement list with those that got merged
        print(f'\n### Update entanglement dict with those that got merged\nNumber of merged pairs: {len(merged_ents)}')
        while len(merged_ents) != 0:
            #print(f'# merged_ents: {len(merged_ents)}')
            for m_ent in merged_ents:
                #print(f'\n{m_ent}')
                for ent_idx, ent in ent_dict.copy().items():
                    #print(ent_idx, ent)
                    if m_ent[0] in ent:
                        #print(f'FOUND MATCH for kept ent {ent_idx}')
                        ent_dict[ent_idx] += [m_ent[1]]
                        merged_ents.remove(m_ent)
            #print(f'# merged_ents: {len(merged_ents)}')

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

        print(f'Results after adding those that got merged to each representative entanglement')
        for uent in updated_ents:
            print(uent)

        ent_data[chain] = updated_ents

    # Step 3
    print(f'\n# Step 3')
    for chain, processed_ents in ent_data.items():

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

                    if chirality1 == chirality2 and check_step_ij_kl_range(proc_ent1, proc_ent2) and np.all(np.abs(proc_cr1 - proc_cr2) <= 20):
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
    
    print(f'Step 3 results')
    for chain, ents in ent_data.items():
        for ent in ents:
            print(ent)

    # Step 4 prep
    print(f'\n# Step 4 prep')
    for chain, new_ents in ent_data.items():

        for ent in new_ents:

            number_of_crossings = len(ent[3:-1])

            chiralites = [each_cr[0] for each_cr in ent[3:-1]]

            chain_num_chirality_key = f"{chain}_{number_of_crossings}_{chiralites}"

            full_entanglement_data[chain_num_chirality_key].append(ent)

    reset_counter = []

    print(f'# Step 4 prep results')
    for chain, ents in full_entanglement_data.items():
        for ent in ents:
            print(ent)

    # Step 4
    print(f'\n# Step 4')
    for chain_num_chiral in full_entanglement_data.keys():

        chain = chain_num_chiral.split("_")[0]

        if chain not in reset_counter:

            reset_counter.append(chain)

            split_cluster_counter = 0

        length_key = defaultdict(list)
        loop_dist = defaultdict(list)
        dups = []
        clusters = {} 
        cluster_count = 0

        pairwise_entanglements = list(itertools.combinations(full_entanglement_data[chain_num_chiral], 2))

        if pairwise_entanglements:

            for i, pairwise_ent in enumerate(pairwise_entanglements):

                dist = loop_distance(pairwise_ent[0], pairwise_ent[1])
                
                if dist <= cut_off and pairwise_ent[0] not in dups and pairwise_ent[1] not in dups:
                    # 1. pair must be <= cut_off
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

        full_ent_values = np.asarray(full_entanglement_data[chain_num_chiral], dtype=object)

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
                rep_chain_ent[f"{chain}_{split_cluster_counter}"].append(rep_entanglement)

            # clusters with a single entnalgement
            else:
                rep_chain_ent[f"{chain}_{split_cluster_counter}"].append(ijr_values[0])
            
            split_cluster_counter += 1

    # Step 5
    print(f'\nStep 5: Output file')
    with open(outfile, "w") as f:

        f.write(f'chain|ijr|gn|gc|num_contacts|contacts\n')
        for chain_counter, ijrs in rep_chain_ent.items():

            chain, counter = chain_counter.split("_")

            for ijr in ijrs:

                new_ijr = (int(ijr[1]), int(ijr[2]), *list(ijr[3:-1]))

                num_nc = int(ijr[0])

                gn, gc = entanglement_partial_g_data[new_ijr]
                gn = float(gn)
                gc = float(gc)

                line = f"{chain}|{(int(ijr[1]), int(ijr[2]), *ijr[3:-1])}|{gn:.5f}|{gc:.5f}|{num_nc}|{ijr[-1]}"
                print(line)
                f.write(f"{line}\n")
    print(f'SAVED: {outfile}')
    return rep_chain_ent


if __name__ == "__main__":

    import multiprocessing as mp
    import sys,os
    import argparse

    cores = len(os.sched_getaffinity(0))

    result_obj = set()

    # cutoffs for different organisms
    # 57 for E.coli
    # 49 for Yeast 
    # 52 for Humans

    #input_data = [(GE_file, 57) for GE_file in os.listdir("unmapped_GE/")]
    #unmapped_GE/P00490_1QM5_A_GE.txt
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("--prot_unmapped_GE_file", type=str, required=True, help="Path to unmapped_GE file created by gaussian_entanglement.py for the protein you want to cluster")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="Path to output directory")
    parser.add_argument("--organism", type=str, required=True, help="Human, Ecoli, Yeast")
    args = parser.parse_args()

    input_data = args.prot_unmapped_GE_file
    outpath = args.outpath
    protein = input_data.split('/')[-1].split('_')[0]
    organism = args.organism

    if organism == 'Human':
        dist_cutoff = 52
    elif organism == 'Ecoli':
        dist_cutoff = 57
    elif organism == 'Yeast':
        dist_cutoff = 49
    else:
        raise ValueError("Must specify Human, Yeast, or Ecoli")

    if not os.path.exists(outpath):
        os.mkdir(outpath)
        print(f'Made: {outpath}')

    global outfile
    outfile = f'{outpath}/{protein}_clustered_GE.txt'

    if os.path.exists(outfile):
        print(f'{outfile} EXISTS. Delete if you wish to cluster')
    else:
        result = cluster_entanglements((input_data,dist_cutoff))
print(f'NORMAL TERMINATION')
