#!/usr/bin/env python2
import os,sys
import tools
from eft_calculator import EFT_calculator
calculator = EFT_calculator('wtr','wtr')
#calculator = EFT_calculator('alc','wtr')
#calculator.setup()
# Please change the following code to whatever needed to generate the input 
# coordinates files
# Please make sure to carry the id number along with the results
root = calculator.com.frg + '_' + calculator.probe.frg + '.random.dat'
if not os.path.exists(root):os.mkdir(root)
def mol2mol_init(ele):
    mol = [[i,0.0,0.0,0.0] for i in ele]
    return mol
size = 2500
folder_id = 0
file_count = 0

def yiledconfs():
    for i in range(1,2001):
        name = '/home/yangjc/git/Xq2EFT/test.dat/random/test%04d.inp.log' % i
        eft, coors = calculator._parseQMlog(name)
        X0, q0 = calculator.com.atomic2Xq(coors[:calculator.com.n])
        X1, q1 = calculator.probe.atomic2Xq(coors[calculator.probe.n:])
        X, q = calculator.Xqfrom2Xq(X0, q0, X1, q1)
        pdb = calculator._Xq2PDB(X,q, bfactor=eft[0])
        yield i, pdb
#confs =  calculator.grid.gen_grid_x()
#for idx, coors in calculator.gen_PDB(confs): 
for idx, coors in yiledconfs(): 
    if  file_count%size == 0:
        folder = os.path.join(root,"EFT_%04d"%(folder_id))
        if not os.path.exists(folder):os.mkdir(folder)
        folder_id += 1
    pdb = open("%s/eft.%04d.pdb"%(folder,idx),"w")
    pdb.write(coors)
    pdb.close()
    file_count += 1

