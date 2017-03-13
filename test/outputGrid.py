#!/usr/bin/env python
import os,sys
from eft_calculator import EFT_calculator

calculator = EFT_calculator()
calculator.setup()
# Please change the following code to whatever needed to generate the input 
# coordinates files
# Please make sure to carry the id number along with the results
calculator.grid.save(
for id, coors in calculator.gen_atomic_coors(): 
    #print id, coors
    if id%size == 0:
        folder = "EFT_%04d"%(id)
        os.mkdir(folder)
    mol = mol2mol_init(ele)
    for i in range(len(coors)):
        for j in range(3):
            mol[i][j+1]=coors[i][j]
    inf = open("%s/eft.%08d.inp"%(folder,id),"w")
    WriteINP(inf, mol)
    inf.close()

