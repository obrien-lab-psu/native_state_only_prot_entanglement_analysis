# Protein entanglement analysis (Native state only)

This repository provides the means to analyze a PDB for native entanglements.

## Usage of code:  
python codes/native_entanglement_analysis_v1.1.py [1] [2] [3] [4]  
[1] = num processors  
[2] = ref PDB or coor file (MDAnalysis supported input coordinate files)  
[3] = ref atom selection mask (in double quotes following MDAnalysis atom selection commands)  
[4] = out_path  

##### example cmd line execution:  
python codes/native_entanglement_analysis.py 4 1p7l_native.pdb "name CA" ./  

## THINGS TO NOTE:  
(1) Make sure there are no missing residues in your PDB. rebuilt them if they are missing. Large sections of missing residues (more than 3 in a row) can cause phantom (not real) entanglements to be found.  
(2) Make sure there is only one alpha carbon per residue. Many PDBs have multiple conformations of a residue labled locA or locB ect. choose one and delete the others.  

## OUTPUT file(s)  
There will be a single pickle file created native_entanglement_analysis_1.1/output/1p7l_native.pkl

This is a binary file containing a nested dictionary. Each top level key is a pair of residue IDs from the PDB constituting a native contact and each value is a dictionary with information about the entanglement.  
  
##### for an example entry:  
  
~~~
{(256, 296): {'gval_N': 0.671, 'gval_C': -0.074, 'Ncrossings': [125], 'Ccrossings': [], 'Nsurr': [123, 124, 125, 126, 127, 256, 276, 277, 278, 279, 280, 281, 282, 283, 296, 297, 298], 'Csurr': []}}  
~~~  
each element of output dictionary explained  
~~~
(256, 296) native contact 
~~~
~~~
gval_N = parital linking value of N terminal thread with loop formed by native contact
~~~
~~~
gval_C = parital linking value of C terminal thread with loop formed by native contact
~~~
~~~
Ncrossings = residue ID of those residues in the N terminus determined to be piercing the plane of the loop
~~~
~~~
Ccrossings = residue ID of those residues in the C terminus determined to be piercing the plane of the loop
~~~
~~~
Nsurr = residue ID of those residues in the N terminus determined to be within 8A of the crossings residues in the N terminus
~~~
~~~
Csurr = residue ID of those residues in the C terminus determined to be within 8A of the crossings residues in the C terminus
~~~
So for the example entry above we have the following interpretation.  
1. Residues 256 and 296 are in contact forming a loop along the backbone.  
2. The N terminal thread has a partial linking value of (gval_N) 0.671 and the C terminal thread has a partial linking value of (gval_C) -0.074 indicating a N terminal entanglement.  
3. The point at which the N terminal thread pierces the loop formed by the native contacts (Ncrossings) is residue 125 and there is no C terminal crossings.  
4. The residues that have alpha carbons within 8A of the crossings residues in the N terminus (Nsurr) are [123, 124, 125, 126, 127, 256, 276, 277, 278, 279, 280, 281, 282, 283, 296, 297, 298].  
