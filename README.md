Protein entanglement analysis (Native state only)

This repository provides the means to analyze a PDB for native entanglements.

Usage of code:  
python codes/native_entanglement_analysis_v1.1.py [1] [2] [3] [4]  
[1] = num processors  
[2] = ref PDB or coor file (MDAnalysis supported input coordinate files)  
[3] = ref atom selection mask (in double quotes following MDAnalysis atom selection commands)  
[4] = out_path  

example cmd line execution:  
python codes/native_entanglement_analysis_v1.1.py 4 1p7l_native.pdb "name CA" ./  

THINGS TO NOTE:  
(1) Make sure there are no missing residues in your PDB. rebuilt them if they are missing. Large sections of missing residues (more than 3 in a row) can cause phantom (not real) entanglements to be found.  
(2) Make sure there is only one alpha carbon per residue. Many PDBs have multiple conformations of a residue labled locA or locB ect. choose one and delete the others.  
