# Xq2EFT
A gird-based QM poteintial energy scriong fuction
Get enegy, force and torque of pair of water or other fragment

-`grid.dat` is 400 k QM grid database
- `grid.py` to organize mesh grid  
- `eft_calculator.py` to calculator enegy, force and torque. (EFT)   
- `tools.py` to convert mol. information
- `gen_coors.py` to write input of GAMESS input .inp file
- `grid.save.py` to save grid in a text file  
- `mol2mol.py` to handle diff mol. format, like .inp .pdb  
- `Q.py` to qualify performance of interplation
- `test_eft_calculator_QM.py` to test the QM grid and compare to MM

## compare RMSE(RMSD) of datasets with different data points.   
QM_grid |Energy |Force  |Torque
---     |---    |---    |---   
400k    |0.2322 |1.0039 |1.4302
900k    |0.3656 |0.8293 |0.7501
1700k   |0.1829 |0.6258 |0.4237
the test set have 2000 dimer water conformation.

