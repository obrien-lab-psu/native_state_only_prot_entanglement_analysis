import os
import math
import re
import logging
import argparse
import numpy as np
import pandas as pd
from glob import glob
import mdtraj as md
from Bio.PDB import PDBParser, is_aa
from Bio import PDB
from scipy.spatial.distance import pdist, squareform
np.set_printoptions(linewidth=np.inf, precision=4)
pd.set_option('display.max_rows', None)

class BioDataProcessor:
    """
    Processes biological data including PDB files, sequence data, and interaction potentials.
    """
    #############################################################################################################
    def __init__(self, PDBfile, chain, outpath, tag, log_file):
        self.PDBfile = PDBfile
        self.chain = chain
        self.outpath = outpath
        self.setup_directories(self.outpath)
        self.tag = tag

        # Initiate a Logging file
        logdir = os.path.join(self.outpath, 'logs/')
        logging.basicConfig(filename=f'{logdir}{log_file}', level=logging.INFO, format='%(asctime)s %(message)s')
        logging.info(f'Intializing entanglement feature analysis')
        print(f'Opened log file: {logdir}{log_file}')

        # Load clustered entanglment file



    #############################################################################################################
    def setup_directories(self, path_to_outdir):
        """
        Setup output directories for storing results.
        """
        #self.outpath = path_to_outdir
        required_dirs = ['uent_features_lib', 'logs']
        for directory in required_dirs:
            full_path = os.path.join(path_to_outdir, directory)
            print(full_path)
            if not os.path.exists(full_path):
                os.makedirs(full_path)
                logging.info(f'Created directory: {full_path}')
                print(f'Created directory: {full_path}')

    #############################################################################################################
    def get_AA(self, pdb_file):

        """
        Get the PDB resid to AA mapping for the provided PDB
        """
        three_to_one_letter = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'MSE': 'M', 'PHE': 'F', 
        'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 
        'VAL': 'V'}

        resid2AA = {}
        # Define the path to your PDB file

        # Create a PDB parser
        parser = PDBParser(QUIET=True)

        # Parse the PDB file
        structure = parser.get_structure("protein", pdb_file)

        # Initialize an empty list to store amino acid codes
        amino_acid_codes = []

        # Iterate through the structure and extract amino acid codes
        for model in structure:
            for chain in model:
                for residue in chain:
                    if is_aa(residue):
                        resname = residue.get_resname()
                        resid = residue.get_id()[1]
                        if resname in three_to_one_letter:
                            AA = three_to_one_letter[resname]
                        else:
                            AA = 'NC'
                        #print(resname, resid, AA)
                        resid2AA[resid] = AA
        self.resid2AA = resid2AA
        self.prot_size = len(resid2AA)


#############################################################################################################
class Analyzer:
    """
    Handles sequence and structural analysis for proteins.
    """
    #############################################################################################################
    def __init__(self, PDBfile, outpath, cluster_file, prot_size):
        self.traj = md.load(PDBfile) 
        print(f'traj: {self.traj}')
        self.outpath = outpath
        self.prot_size = prot_size

        ## parse lines to get native contacts, crossings,
        if os.path.exists(cluster_file):
            self.ent_data = np.asarray([x.strip('\n') for x in open(cluster_file, 'r').readlines()])
            print(self.ent_data)
        else:
            raise ValueError(f"{self.cluster_file} does not exits")


    #############################################################################################################
    def get_uent_features(self, chain, tag):
        """
        Get the features for each unique entanglement provided in the clustered_unampped_GE file
        """

        uent_df = {'chain':[], 
                    'ENT-ID':[],
                    'Gn':[],
                    'N_term_thread':[],
                    'Gc':[],
                    'C_term_thread':[],
                    'unmapped-NC':[],
                    'unmapped-NC_wbuff':[],
                    'unmapped-crossings':[], 
                    'unmapped-crossings_wbuff':[], 
                    'loopsize': [], 
                    'num_zipper_nc':[], 
                    'perc_bb_loop':[],
                    'num_loop_contacting_res':[],
                    'num_cross_nearest_neighbors':[],
                    'ent_coverage':[],
                    'min_N_prot_depth_left':[],
                    'min_N_thread_depth_left':[],
                    'min_N_thread_slippage_left':[],
                    'min_C_prot_depth_right':[],
                    'min_C_thread_depth_right':[],
                    'min_C_thread_slippage_right':[], 
                    'prot_size':[], 
                    'ACO':[],
                    'RCO':[]}

        #############################################################################################################################################################################
        ### Load entanglement information if present
        topology = self.traj.topology

        # get mapping of chain letters to chain index 
        chain_ids = {chain.chain_id: chain.index for chain in topology.chains}
        print(chain_ids)
        if chain not in chain_ids:
            raise ValueError(f'chain {chain} not in PDB file')

        # Get alpha carbon indices from PDB so I can find contacting residues
        #haystack_alpha_carbon_indices = self.traj.top.select(f'name CA')
        haystack_alpha_carbon_indices = self.traj.top.select(f'chainid {chain_ids[chain]} and name CA')
        print(f'haystack_alpha_carbon_indices: {haystack_alpha_carbon_indices} {len(haystack_alpha_carbon_indices)}')


        ## parse lines to get native contacts, crossings,
        ent_present = True
        rbuffer = 3
        pdb_NC_list = [] # list of PDB native contact residues +/- rbuffer
        pdb_NC_core_list = [] # list of PDB natvie contact residues
        pdb_crossing_list = [] # list of PDB crossing residues +/- rbuffer
        pdb_crossing_core_list = [] # list of PDB crossing residues
        total_ent_res = set()
        resid_dict = {}

        for index, line in enumerate(self.ent_data):
            if index == 0:
                continue
            logging.info(f'#######: ENT-ID: {index}')
            ent_core = []

            #['Chain C ', " (23, 315, '+7') ", ' 0.7708439147210744 ', ' -0.18254831113041398', 2]
            #['Chain C', "(23, 315, '+7')", '0.77084', '-0.18255', '10', '23-315;26-322;26-325;26-326;26-329;27-325;28-325;22-322;23-321;23-325']
            line = line.split("|")
            print(line)
            logging.info(line)

            ## check that the entanglement isnt in a non-mapped area. if so skip it
            #line = line[1].split(',')
            pdb_NCi_core = line[1].split(',')[0]
            pdb_NCj_core = line[1].split(',')[1]
            pdb_NCi_core = int(float(pdb_NCi_core.replace("(", "").strip()))
            pdb_NCj_core = int(float(pdb_NCj_core.strip()))
            pdb_crossing_res_core = [abs(int(float(re.findall(r'\d+', cross)[0]))) for cross in line[1].split(',')[2:]]
            total_core_ent_res = [pdb_NCi_core, pdb_NCj_core] + pdb_crossing_res_core
            print(pdb_NCi_core, pdb_NCj_core, pdb_crossing_res_core, total_core_ent_res)

            uent_df['chain'] += [chain]
            uent_df['ENT-ID'] += [index]


            #########################################################################
            ## get Gn and Gc and if it is present the cluster size
            if len(line) == 6:
                num_zipper_nc = int(line[-2])
            else:
                num_zipper_nc = np.nan
            Gn = float(line[2])
            Gc = float(line[3])

            # Calcualte the absolute and relative contact orders
            range_strings = line[-1].split(';')
            loops = [(int(x[0]), int(x[1])) for x in [l.split('-') for l in range_strings]]
            #loops = self.parse_ranges(range_strings)
            #print(loops)
            #loops = [(int(i), int(j)) for i,j in [l.rsplit('-', 1) for l in line[-1].split(';')]]
            loop_sizes = [j-i for i,j in loops]
            print(f'loop_sizes: {loop_sizes}')
            ACO = np.sum(loop_sizes)/len(loop_sizes)
            RCO = ACO/self.prot_size
            print(f'Gn: {Gn} | Gc: {Gc} | num_zipper_nc: {num_zipper_nc} | ACO: {ACO} | RCO: {RCO}')
            logging.info(f'Gn: {Gn} | Gc: {Gc} | num_zipper_nc: {num_zipper_nc} | ACO: {ACO} | RCO: {RCO}')
            uent_df['Gn'] += [Gn]
            uent_df['Gc'] += [Gc]
            uent_df['num_zipper_nc'] += [num_zipper_nc]
            uent_df['ACO'] += [ACO]
            uent_df['RCO'] += [RCO]


            #########################################################################
            #get PDB native contact and those +/- rbuffer along the primary structure
            line = line[1].split(',')
            pdb_NCi_core = line[0]
            pdb_NCj_core = line[1]
            pdb_NCi_core = int(float(pdb_NCi_core.replace("(", "").strip()))
            pdb_NCj_core = int(float(pdb_NCj_core.strip()))
            pdb_NC_core = [pdb_NCi_core, pdb_NCj_core]
            pdb_NC_core_list += pdb_NC_core

            pdb_NCi = np.arange(pdb_NCi_core - rbuffer, pdb_NCi_core + rbuffer + 1)
            pdb_NCj = np.arange(pdb_NCj_core - rbuffer, pdb_NCj_core + rbuffer + 1)
            pdb_NC = np.hstack([pdb_NCi, pdb_NCj]).tolist()
            pdb_NC_list += pdb_NC

            logging.info(f'pdb_NC: {pdb_NC}')
            logging.info(f'pdb_NC_core: {pdb_NC_core}')
            uent_df['unmapped-NC'] += [",".join([str(r) for r in pdb_NC_core])]
            uent_df['unmapped-NC_wbuff'] += [",".join([str(r) for r in pdb_NC])]

            loopsize = pdb_NCj_core - pdb_NCi_core
            loop_resids = np.arange(pdb_NCi_core, pdb_NCj_core + 1)
            loop_alpha_carbon_indices = [atom.index for atom in self.traj.top.atoms if (atom.name == 'CA' and atom.segment_id == chain) if atom.residue.resSeq in loop_resids]
            loop_nearest_neighbors_indices = md.compute_neighbors(self.traj, 0.8, query_indices=loop_alpha_carbon_indices, haystack_indices=haystack_alpha_carbon_indices)[0]
            num_loop_contacting_res = len(loop_nearest_neighbors_indices)

            logging.info(f'num_loop_contacting_res: {num_loop_contacting_res}')
            uent_df['loopsize'] += [loopsize]
            uent_df['perc_bb_loop'] += [loopsize/self.prot_size]
            uent_df['num_loop_contacting_res'] += [num_loop_contacting_res]
            #########################################################################


            #########################################################################
            #get PDB crossings and those +/- rbuffer along the primary structure
            pdb_crossing_res_core = [abs(int(float(re.findall(r'\d+', cross)[0]))) for cross in line[2:]]
            pdb_crossing_res = np.hstack([np.arange(int(x) - rbuffer, int(x) + rbuffer + 1) for x in pdb_crossing_res_core]).tolist()
            logging.info(f'pdb_crossing_res: {pdb_crossing_res}')
            logging.info(f'pdb_crossing_res_core: {pdb_crossing_res_core}')

            pdb_crossing_list += pdb_crossing_res
            pdb_crossing_core_list += pdb_crossing_res_core
            uent_df['unmapped-crossings'] += [",".join([str(c) for c in pdb_crossing_res_core])]
            uent_df['unmapped-crossings_wbuff'] += [",".join([str(c) for c in pdb_crossing_res])]

            ### Get residues in contact with crossing residues +/- cbuff
            #cross_alpha_carbon_indices = [atom.index for atom in self.traj.top.atoms if atom.name == 'CA' if atom.residue.resSeq in pdb_crossing_res]
            cross_alpha_carbon_indices = [atom.index for atom in self.traj.top.atoms if (atom.name == 'CA' and atom.segment_id == chain) if atom.residue.resSeq in pdb_crossing_res]
            cross_nearest_neighbors_indices = md.compute_neighbors(self.traj, 0.8, query_indices=cross_alpha_carbon_indices, haystack_indices=haystack_alpha_carbon_indices)[0]
            num_cross_nearest_neighbors = len(cross_nearest_neighbors_indices)

            logging.info(f'num_cross_nearest_neighbors: {num_cross_nearest_neighbors}')
            uent_df['num_cross_nearest_neighbors'] += [num_cross_nearest_neighbors]
            #########################################################################


            #########################################################################
            ## Get number of threads in each termini and depth
            N_term_thread = [c for c in pdb_crossing_res_core if c < pdb_NCi_core]            
            num_N_term_thread = len(N_term_thread)
            C_term_thread = [c for c in pdb_crossing_res_core if c > pdb_NCj_core]            
            num_C_term_thread = len(C_term_thread)
            logging.info(f'N_term_thread: {N_term_thread}')

            logging.info(f'C_term_thread: {C_term_thread}')
            uent_df['N_term_thread'] += [num_N_term_thread]
            uent_df['C_term_thread'] += [num_C_term_thread]

            if num_N_term_thread != 0:
                min_N_thread_slippage_left = min(N_term_thread)
                min_N_thread_depth_left = min_N_thread_slippage_left / pdb_NCi_core
                min_N_prot_depth_left = min_N_thread_slippage_left / self.prot_size
            else:
                min_N_thread_slippage_left = np.nan
                min_N_thread_depth_left = np.nan
                min_N_prot_depth_left = np.nan
            uent_df['min_N_thread_slippage_left'] += [min_N_thread_slippage_left]
            uent_df['min_N_thread_depth_left'] += [min_N_thread_depth_left]
            uent_df['min_N_prot_depth_left'] += [min_N_prot_depth_left]

            if num_C_term_thread != 0:
                min_C_thread_slippage_right = self.prot_size - max(C_term_thread)
                min_C_thread_depth_right = min_C_thread_slippage_right / (self.prot_size - pdb_NCj_core)
                min_C_prot_depth_right = min_C_thread_slippage_right / self.prot_size
            else:
                min_C_thread_slippage_right = np.nan
                min_C_thread_depth_right = np.nan
                min_C_prot_depth_right = np.nan
            uent_df['min_C_thread_slippage_right'] += [min_C_thread_slippage_right]
            uent_df['min_C_thread_depth_right'] += [min_C_thread_depth_right]
            uent_df['min_C_prot_depth_right'] += [min_C_prot_depth_right]
            #########################################################################
            

            #########################################################################
            ### Get entangled residues. Those that are within 8A of the core residues
            print('Get entangled residues. Those that are within 8A of the core residues')
            ent_core = set(pdb_NC).union(set(pdb_crossing_res))
            ent_core = list(ent_core)
            logging.info(f'ent_core: {ent_core}')
            
            ## Get alpha carbon indices
            haystack_alpha_carbon_indices = self.traj.top.select('name CA')
            #print(haystack_alpha_carbon_indices)

            ent_res = []
            for ent_core_res in ent_core:
                alpha_carbon_indices = [atom.index for atom in self.traj.top.atoms if (atom.name == 'CA' and atom.segment_id == chain) and atom.residue.resSeq in [ent_core_res]]

                # Calculate nearest neighbors within 8 angstroms based on alpha carbon coordinates
                nearest_neighbors_indices = md.compute_neighbors(self.traj, 0.8, query_indices=alpha_carbon_indices, haystack_indices=haystack_alpha_carbon_indices, periodic=False)[0]
                for atom_index in nearest_neighbors_indices:
                    dist = math.sqrt(sum((x - y) ** 2 for x, y in zip(self.traj.xyz[0][alpha_carbon_indices[0]], self.traj.xyz[0][atom_index])))

                    print(ent_core_res, alpha_carbon_indices, self.traj.xyz[0][alpha_carbon_indices[0]], self.traj.top.atom(atom_index).residue.resSeq, atom_index, self.traj.xyz[0][atom_index], dist)
                    logging.info(f'DIST_check: {ent_core_res}, {alpha_carbon_indices}, {self.traj.xyz[0][alpha_carbon_indices[0]]}, {self.traj.top.atom(atom_index).residue.resSeq}, {atom_index}, {self.traj.xyz[0][atom_index]}, {dist}')
                    if float(f'{dist:.4f}') > 0.8:
                        print(f'ERROR distance greater than 8A.. See logs')
                        quit()

                nearest_neighbors_resid = [self.traj.top.atom(atom_index).residue.resSeq for atom_index in nearest_neighbors_indices]
                nearest_neighbors_resid = list(set(nearest_neighbors_resid))
                nearest_neighbors_res_idx = [self.traj.top.atom(atom_index).residue.index for atom_index in nearest_neighbors_indices]
                num_nearest_neighbors = len(nearest_neighbors_indices)

                ent_res += nearest_neighbors_resid

    
            ent_res = set(ent_res).union(set(ent_core))
            #print(f'ent_res: {ent_res} {len(ent_res)}')
            logging.info(f'ent_res: {ent_res} {len(ent_res)}')
            uent_df['ent_coverage'] += [len(ent_res)/self.prot_size]
            uent_df['prot_size'] += [self.prot_size]
            total_ent_res = set(ent_res).union(total_ent_res)
            for res in ent_res:

                if res in resid_dict:
                    if isinstance(resid_dict[res], list):
                        resid_dict[res] += [str(index)]

                    else:
                        #print(res, resid_dict[res]) 
                        resid_dict[res] = [str(index)]
                else:
                    resid_dict[res] = [str(index)]
            
            nonent_res = [self.traj.top.atom(atom_index).residue.resSeq for atom_index in haystack_alpha_carbon_indices if atom_index not in nearest_neighbors_indices]
            #print(f'nonent_res: {nonent_res}')
            for res in nonent_res:
                if res in resid_dict:
                    continue
                else:
                    resid_dict[res] = np.nan
            #########################################################################

        ### save file for unique entanglement features
        uent_df = pd.DataFrame(uent_df)
        print(f'uent_df:\n{uent_df}')
        uent_outfile = f'{self.outpath}uent_features_lib/{tag}_native_uents.csv'
        uent_df.to_csv(uent_outfile, sep='|', index=False)
        print(f'SAVED: {uent_outfile}')
        logging.info(f'SAVED: {uent_outfile}')
        ########################################################################################################################




#############################################################################################################
def main():

    # Parse user provided arguments
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-o", "--outpath", type=str, required=True, help=f"path to outdir")
    parser.add_argument("-p", "--PDB", type=str, required=True, help="PDB to process")
    parser.add_argument("-c", "--chain", type=str, required=True, help="Chain of the PDB to use")
    parser.add_argument("-l", "--log_file", type=str, required=True, help="logging file name")
    parser.add_argument("-t", "--tag", type=str, required=True, help="tag for output file name")
    parser.add_argument("--cluster_file", type=str, required=True, help="path to clustered entanglement file")
    args = parser.parse_args()

    # Make BioDataProcessor object and get contact potential data, stride data, LiPMS data, canonical sequence, mapping
    data_processor = BioDataProcessor(
            PDBfile = args.PDB,
            chain = args.chain,
            outpath = args.outpath,
            tag = args.tag, 
            log_file= args.log_file)
    print(data_processor)

    # Get PDB resid to AA mapping
    print(f'\nGet PDB resid to AA mapping')
    data_processor.get_AA(data_processor.PDBfile)
    print(data_processor.resid2AA)

    # Initiate the Analyzer object
    print(f'\nInitiate the Analyzer object')
    analyzer = Analyzer(data_processor.PDBfile, data_processor.outpath, args.cluster_file, data_processor.prot_size)

    # Get unique entanglement features
    print(f'\nGet unique entanglement features')
    analyzer.get_uent_features(args.chain, args.tag)

    logging.info("Finished processing")
    print("Processing completed successfully.")

if __name__ == "__main__":
    main()

