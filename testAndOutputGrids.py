#!/usr/bin/env python2
import numpy as np
import pdb
from random import sample
from time import time
import heapq
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import sys, os
from eft_calculator import EFT_calculator, Water
import tools


def load_coordinates(name):
    lines = open('test.dat/random/'+name).readlines()[-7:-1]
    coors = [[float(item) for item in line.split()[2:5]] for line in lines]
    return np.array(coors)

class Classical_calculator:
    def __init__(self):
        self.eps = [0.12, 0.046, 0.046]
        self.sigma = [1.7, 0.2245, 0.2245]
        self.charge = [-0.834, 0.417, 0.417]

    def eval(self, coors):
        mol = Water()
        coor0 = coors[:3]
        coor1 = coors[3:]
        e = 0.
        f = np.zeros(3)
        t = np.zeros(3)
        com1 = mol.getCOM(coor1)
        eps, sigma, charge = self.eps, self.sigma, self.charge
        for i in range(3):
            for j in range(3):
                ener, force = self.atomicEF(coor0[i], eps[i], sigma[i], charge[i], coor1[j], eps[j], sigma[j], charge[j])
                e += ener
                f += force
                t += np.cross(coor1[j]-com1, force)
       #if e>100.0:
       #    e = 100.0
       #    f = f/np.linalg.norm(f) * 100.0
       #    t = t/np.linalg.norm(t) * 100.0
        return np.array([e, f[0], f[1], f[2], t[0], t[1], t[2]])

    def atomicEF(self, x0, e0, s0, q0, x1, e1, s1, q1):
        k = 138.935456 
        e = np.sqrt(e0 * e1)
        s = s0 + s1
        r = np.linalg.norm(x0 - x1)
        if r <0.1 : return 100.0, np.array([100., 100.,100.,])
        sor6 = (s/r) ** 6
        evdw = e * (sor6**2 - 2 * sor6)
        fvdw = e / r**2 * sor6 * (sor6 - 1) * (x1 - x0)
        eelec = k * q0 * q1 / r
        felec = k * q0 * q1 / r**3 * (x1 - x0)
        ener = evdw + eelec
        force = fvdw + felec
        return ener, force

def test_random_set():
    e0 = []
    e1 = []
    fce0 = []
    fce1 = []
    trq0 = []
    trq1 = []
    all = []
    t1 = time()
    for i in range(1, 2000):
        # load atomic coor 
        name = 'test.dat/random/test%04d.inp.log' % i
        #if i == 1693: pdb.set_trace()
        eft, coors = calculator._parseQMlog(name)
        # evaluate with analytical function
        eft = cc.eval(coors)
        e0.append(eft[0])
        fce0 += list(eft[1:4])
        trq0 += list(eft[4:7])
        # convert atomic coor to r, phi, theta... 
        X0, q0 = calculator.mol.atomic2Xq(coors[:3])
        X1, q1 = calculator.mol.atomic2Xq(coors[3:])
        # evaluate with calculator
        eft = calculator.eval(X0, q0, X1, q1)
        e1.append(eft[0])
        #if eft[0] > 15:
        #    print(coors, name)
        #    print(np.dtype(q1[0]))
        fce1 += list(eft[1:4])
        trq1 += list(eft[4:7])
        #all.append((-np.abs(e0[-1]-e1[-1]), name))
        all.append((-np.linalg.norm(np.array(fce0) - np.array(fce1)), name))
    t2 = time()
    print('took %.1f s to evaluate the random set' % (t2 - t1))
    heapq.heapify(all)
    #for i in range(3):
    #    de, name = heapq.heappop(all)
    #    print -de, name
    """
    for i in range(len(e0)):
        if e1[i]> 100.0:
            e0[i] = e1[i] = 0.0
            for j in range(3):
                fce0[i*3 +j ] = fce1[i*3+j] = trq0[i*3+j] = trq1[i*3+j] = 0.0
   """
   # make a plot
    _, axarr = plt.subplots(1, 3)
    p = np.corrcoef(e0, e1)[0, 1]
    print("Energy: p =", p)
    axarr[0].scatter(e0, e1)
    axarr[0].text(0, 0, 'p=%.4f'%p)
    p = np.corrcoef(fce0, fce1)[0, 1]
    print("Force: p =", p)
    axarr[1].scatter(fce0, fce1)
    axarr[1].text(0, 0, 'p=%.4f'%p)
    p = np.corrcoef(trq0, trq1)[0, 1]
    print("Torque: p =", p)
    axarr[2].scatter(trq0, trq1)
    axarr[2].text(0, 0, 'p=%.4f'%p)
    plt.savefig(figname)

def randomSample():
    root = 'golden.dat'
    if not os.path.exists(root):os.mkdir(root)
    def mol2mol_init(ele):
        mol = [[i,0.0,0.0,0.0] for i in ele]
        return mol
    size = 200
    folder_id = 0
    file_count = 0
    confs = calculator.grid._iter_conf()
    confs = list(confs)
    if len(confs) > 2000:
        confs = sample(list(confs), 2000) 
    for idx, coors in calculator.gen_PDB(confs):
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

def grids_conf():
    root = figname[:-4] + '.grids.dat' 
    if not os.path.exists(root):os.mkdir(root)
    def mol2mol_init(ele):
        mol = [[i,0.0,0.0,0.0] for i in ele]
        return mol
    size = 200
    folder_id = 0
    file_count = 0
    confs = calculator.grid._grid_conf()
    for idx, coors in calculator.gen_PDB(confs):
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

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("\n       Usage:#0 figname.png [datfilename.dat err_cutoff]\n")
        sys.exit()
    figname = sys.argv[1] # a output fig name
    databaseName = sys.argv[2]

    t0 = time()
    cc = Classical_calculator()
    if os.path.exists(databaseName):
        print("loaded a old database")
        calculator = EFT_calculator(databaseName)
    else:
        print("created a new mesh")
        calculator = EFT_calculator()
    if len(sys.argv) == 4:
        error_cutoff = float(sys.argv[3])
        print("set cutoff as %f"%(error_cutoff))
        calculator.fill_grid(cc, databaseName, error_cutoff)
    t1 = time()
    print('took %.1f s to fill the grid' % (t1 - t0))
    test_random_set()
    #randomSample()
    grids_conf()

