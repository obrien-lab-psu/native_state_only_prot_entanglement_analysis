# Protein Entanglement Detection (gaussian_entanglement.py)

This repository contains a suite of scripts for detecting non-covalent lasso entanglements in protein structures using mass spectrometry data. The primary functionality includes preprocessing PDB files, calculating native entanglements, and identifying crossing residues.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Functions](#functions)
- [Contributing](#contributing)
- [License](#license)

## Installation

### Prerequisites

Ensure you have the following Python packages installed:

- `numpy` (requires 1.2X)
- `pandas`
- `MDAnalysis`
- `scipy`
- `numba`
- `topoly` (can be installed via pip)

You can install the required packages using:

```sh
pip install numpy pandas MDAnalysis scipy numba topoly
```

### Cloning the Repository

Clone this repository to your local machine:

```sh
git clone https://github.com/yourusername/protein-entanglement-detection.git
cd protein-entanglement-detection
```

## Usage

### Command Line Interface

You can also use the provided command line interface:

```sh
python your_script.py --PDB path/to/pdb --GLN_threshold 0.6 --topoly_density 1 --Calpha False
```
--PDB is the required path to the PDB file to analyze  
--GLN_threshold is an optional parameter to vary the threshold we find real entanglements. Any entanglement with a absolute value of the GLN above this threshold will be used. Default = 0.6  
--topoly_density is an optional parameter to vary the density of triangles when drawing the minimal loop surface to determine peircing events   
--Calpha is a True or False optional parameter to use either 8A between alpha carbons for native contacts or 4.5A between heavy atoms (False, default)  
  
## Functions

### `pre_processing_pdb(pdb_file: str) -> None`

Pre-processes the PDB files by removing everything after the last `TER`. Handles `.pdb`, `.pdb1`, and Alphafold PDBs.

### `helper_dot(Runit: np.ndarray, dR_cross: np.ndarray) -> list`

A Numba-accelerated function to speed up dot product calculations.

### `point_rounding(num: float) -> float`

Rounds numbers to the nearest 0.6.

### `get_entanglements(...) -> dict`

Identifies proteins containing non-covalent lasso entanglements.

### `find_missing_residues(resids: np.ndarray) -> np.ndarray`

Finds missing residues in the PDB file.

### `loop_filter(native_contacts: dict, resids: np.ndarray, missing_res: np.ndarray) -> dict`

Filters loops based on missing residues.

### `find_crossing(coor: np.ndarray, nc_data: dict, resids: np.ndarray) -> dict`

Uses Topoly to find crossings based on partial linking numbers.

### `crossing_filter(entanglements: dict, missing_res: np.ndarray) -> dict`

Filters entanglements based on missing residues near crossings.

## Output

### Raw entanglement file

File containing 1 entanglement per line separated by "|".  
[0] ChainID  
[1] (loop native contact residue i, loop native contact residue j, crossings)  
[2] Gn  
[3] Gc  
[4] Whether a disulfide bond was identified at the native contact (two SG atoms within 2.2 A)  

## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any improvements or bug fixes.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

 
 
# Protein Entanglement Clustering (clustering.py)

This repository contains scripts for clustering non-covalent lasso entanglements in protein structures. The primary functionality includes calculating distances between entanglements, identifying minimal loops, and spatially clustering entanglements based on residue crossings and chiralities.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Functions](#functions)
- [Contributing](#contributing)
- [License](#license)

## Installation

### Prerequisites

Ensure you have the following Python packages installed:

- `numpy`
- `scipy`
- `geom_median`
- `collections`
- `itertools`
- `functools`

You can install the required packages using:

```sh
pip install numpy scipy geom_median
```

### Cloning the Repository

Clone this repository to your local machine:

```sh
git clone https://github.com/yourusername/protein-entanglement-clustering.git
cd protein-entanglement-clustering
```

## Usage

### Clustering Entanglements

Cluster non-covalent lasso entanglements in protein structures.

```python
from your_script import cluster_entanglements

# Example usage
result = cluster_entanglements(('path/to/unmapped_GE_file.txt', 57))
```

### Command Line Interface

You can also use the provided command line interface:

```sh
python codes/clustering.py -o clustered_unmapped_GE_HQ/Human_no_pure_slipknots/ --prot_unmapped_GE_file unmapped_GE_HQ/Human_no_pure_slipknots/Q9NR16_AF_A_GE.txt --organism Human
```
-o is the output directory  
--prot_unmapped_GE_file is the path to the raw entanglement file  
--organism is the organism you want to use which is either Human, Ecoli, or Yeast  

## Functions

### `loop_distance(entangled_A: tuple, entangled_B: tuple)`

Calculates the Euclidean distance between two entanglements.

### `check_step_ij_kl_range(ent1: tuple, ent2: tuple)`

Checks if the indices of one entanglement pair fall within the range of another.

### `cluster_entanglements(GE_file_w_cutoff: tuple) -> dict`

Clusters entanglements based on residue crossings and chiralities, and outputs the clustered results.

### Other Utility Functions

- `find_crossing_residues()`: Identifies crossing residues in the protein structure.
- `filter_loops()`: Filters loops based on certain criteria.
- `calculate_distances()`: Calculates distances between residue sets.

## Output

### Raw entanglement file

File containing 1 entanglement per line separated by "|".  
[0] ChainID  
[1] (loop native contact residue i, loop native contact residue j, crossings) representaive entanglement for the cluster  
[2] Gn  
[3] Gc  
[4] number of loop closing contacts  
[4] ";" separated list of the residueIDs for the loop closing contacts  
[4] Whether a disulfide bond was identified in one of the loop forming native contact (two SG atoms within 2.2 A)  
 
## Contributing

Contributions are welcome! Please open an issue or submit a pull request for any improvements or bug fixes.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


# Scan PDB for dislufide bonds 

To scane a PDB for disulfide bonds use the following script.
Residues must have side chain heavy atoms within 4.5A and be atleast 4 residues apart along the primary structure

### Prerequisites

Ensure you have the following Python packages installed:

- `numpy`
- `argparse`
- `os`
- `sys`
- `pandas`
- `Bio`


```python
python codes/scan_disulfide_bonds.py --PDB [path/to/pdb-file or directory containing pdbs]
```
It will print to screen the found CYS-CYS bonds in each chain along with their min distance between sidechain heavy atoms
It will also save a csv of all the contacts found


# Calculate the entanglement features 
Calculate various entanglement parameters

### Prerequisites

Ensure you have the following Python packages installed:

- `numpy`
- `argparse`
- `os`
- `sys`
- `pandas`
- `mdtraj`
- `Bio`

```python
usage: EntFeatures.py [-h] -o OUTPATH -p PDB -c CHAIN -l LOG_FILE -t TAG --cluster_file CLUSTER_FILE

Process user specified arguments

options:
  -h, --help            show this help message and exit
  -o OUTPATH, --outpath OUTPATH
                        path to outdir
  -p PDB, --PDB PDB     PDB to process
  -c CHAIN, --chain CHAIN
                        Chain of the PDB to use
  -l LOG_FILE, --log_file LOG_FILE
                        logging file name
  -t TAG, --tag TAG     tag for output file name
  --cluster_file CLUSTER_FILE
                        path to clustered entanglement file

python codes/EntFeatures.py --PDB test_pdbs/P0AD61_4YNG_C.pdb -c C -l EntFeatures_P0AD61_4YNG_C.log -t P0AD61 -o EntFeatures/ --cluster_file clustered_unmapped_GE/P0AD61_clustered_GE.txt
```
![Entanglement Complexity Definitions](entanglement_complexity_metrics.jpg)


# LASSO regression of entanglement complexity features
Takes a directory of the entanglement complexity features and two lists of uniprot accession IDs and permutation tests for differences in the means and medians of each group and LASSO regression for feature selection to discriminate between each group. 

### Prerequisites

Ensure you have the following Python packages installed:

- `numpy`
- `argparse`
- `os`
- `sys`
- `pandas`
- `sklearn`
- `glob`
- `scipy`

```python
usage: Compare_ent_complexity_metrics_generic.py [-h] -g1 GROUP1_GENE_LIST -g2 GROUP2_GENE_LIST -l LOG_FILE -e UENT_FILES -o OUTPATH -p NUM_PERMUTE

Process user specified arguments

options:
  -h, --help            show this help message and exit
  -g1 GROUP1_GENE_LIST, --Group1_gene_list GROUP1_GENE_LIST
                        path to Group1 gene list used for mask. one uniprot name per line. 
  -g2 GROUP2_GENE_LIST, --Group2_gene_list GROUP2_GENE_LIST
                        path to Group2 gene list used for mask. one uniprot name per line.
  -l LOG_FILE, --log_file LOG_FILE
                        Path to logging file
  -e UENT_FILES, --uent_files UENT_FILES
                        path to unique entanglement files
  -o OUTPATH, --outpath OUTPATH
                        path to output directory. will be made if doesnt exist
  -p NUM_PERMUTE, --num_permute NUM_PERMUTE
                        Number of permutations
```

