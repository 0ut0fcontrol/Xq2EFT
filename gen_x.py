#!/usr/bin/env python2
import os,sys
from eft_calculatorRQ import EFT_calculator
from mol2mol import * #GAMESS_Settings, WriteINP, WritePDB
#  The coordinate structure of intermediate data is:
#      [[elem, x, y, z],
#       [elem, x, y, z],
#       ...]
ele = "OHHOHH"
GAMESS_Settings="""! Gamess input file generated by gen_coors.py for HuaHe.
 $CONTRL SCFTYP=RHF RUNTYP=GRADIENT MPLEVL=2 $END
 $SYSTEM MWORDS=500 MEMDDI=500 $END
-$CONTRL EXETYP=CHECK $END
 $CONTRL MAXIT=200 $END
-$STATPT OPTTOL=1.0E-5 $END
-$STATPT NSTEP=300 $END
-$DFT NTHE=36 NPHI=72 $END
 $SCF DIRSCF=.T. $END
 $BASIS GBASIS=N311 NGAUSS=6 NDFUNC=2 NPFUNC=2 DIFFSP=.T. DIFFS=.T. $END
 $GUESS GUESS=HUCKEL $END
 $DATA
--Cartesian coordinates with C1 symmetry as follows:
C1
"""

calculator = EFT_calculator()
#calculator.setup()
# Please change the following code to whatever needed to generate the input 
# coordinates files
# Please make sure to carry the id number along with the results
root = 'pdbRQ.dat'
if not os.path.exists(root):os.mkdir(root)
def mol2mol_init(ele):
    mol = [[i,0.0,0.0,0.0] for i in ele]
    return mol
size = 200
folder_id = 0
file_count = 0
for idx, coors in calculator.gen_PDB(): 
#for id, coors in calculator.gen_atomic_coors(0,10): 
    #print(idx, coors)

    if  file_count%size == 0:
        folder = os.path.join(root,"EFT_%04d"%(folder_id))
        if not os.path.exists(folder):os.mkdir(folder)
        folder_id += 1
    pdb = open("%s/eft.%s.pdb"%(folder,idx),"w")
    pdb.write(coors)
    pdb.close()
    file_count += 1

