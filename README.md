# Protein Entanglement Detection

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

- `numpy`
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
python your_script.py --PDB path/to/pdb --GLN_threshold 0.5
```

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

