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

### Preprocessing PDB Files

Preprocess PDB files by removing everything after the last `TER` and handling multiple model instances.

```python
from your_script import pre_processing_pdb

pre_processing_pdb('example.pdb')
```

### Calculate Native Entanglements

Calculate and output native lasso-like self-entanglements and missing residues for PDB files and their chains.

```python
from your_script import calculate_native_entanglements

calculate_native_entanglements('path/to/pdb/file.pdb')
```

### Command Line Interface

You can also use the provided command line interface:

```sh
python your_script.py --PDB path/to/pdb --GLN_threshold 0.6 --topoly_density 1
```
--PDB is the required path to the PDB file to analyze
--GLN_threshold is an optional parameter to vary the threshold we find real entanglements. Any entanglement with a absolute value of the GLN above this threshold will be used. Default = 0.6
--topoly_density is an optional parameter to vary the density of triangles when drawing the minimal loop surface to determine peircing events 

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
