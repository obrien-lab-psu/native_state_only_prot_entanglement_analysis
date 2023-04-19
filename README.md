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
