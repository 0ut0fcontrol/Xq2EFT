#!/usr/bin/env python2
import sys
import numpy as np
from time import time
import heapq
from matplotlib import pyplot as plt
from Q import uqi, rmse
from eft_calculator import EFT_calculator, Water
import tools


def load_coordinates(name):
    lines = open('/Users/he/Work/rigid_water/qm/new_python/random/'+name).readlines()[-7:-1]
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
        return np.array([e, f[0], f[1], f[2], t[0], t[1], t[2]])

    def atomicEF(self, x0, e0, s0, q0, x1, e1, s1, q1):
        k = 138.935456 
        e = np.sqrt(e0 * e1)
        s = s0 + s1
        r = np.linalg.norm(x0 - x1)
        sor6 = (s/r) ** 6
        evdw = e * (sor6**2 - 2 * sor6)
        fvdw = e / r**2 * sor6 * (sor6 - 1) * (x1 - x0)
        eelec = k * q0 * q1 / r
        felec = k * q0 * q1 / r**3 * (x1 - x0)
        ener = evdw + eelec
        force = fvdw + felec
        return ener, force

def test_random_set(log_list):
    e0 = []
    e1 = []
    e2 = []
    fce0 = []
    fce1 = []
    fce2 = []
    trq0 = []
    trq1 = []
    trq2 = []
    all  =  []
    t1 = time()
    for name in open(log_list, "r").readlines():
        # load atomic coor 
        name = name.rstrip()
        # QM value
        eft, coors = calculator._parseQMlog(name)
        e0.append(eft[0]) 
        fce0 += list(eft[1:4])
        trq0 += list(eft[4:7])
        # convert atomic coor to r, phi, theta... 
        X0, q0 = calculator.mol.atomic2Xq(coors[3:])
        X1, q1 = calculator.mol.atomic2Xq(coors[:3])
        # evaluate with calculator
        eft = calculator.eval(X0, q0, X1, q1)
        e1.append(eft[0])
        fce1 += list(eft[1:4])
        trq1 += list(eft[4:7])
        # evaluate with analytical function
        eft = cc.eval(coors)
        e2.append(eft[0])
        fce2 += list(eft[1:4])
        trq2 += list(eft[4:7])
        #all.append((-np.abs(e0[-1]-e1[-1]), name))
        all.append((-np.linalg.norm(np.array(fce0) - np.array(fce1)), name))
    t2 = time()
    print 'took %.1f s to evaluate the random set' % (t2 - t1)
    heapq.heapify(all)
    #for i in range(3):
    #    de, name = heapq.heappop(all)
    #    print -de, name
    # make a plot
    with open('result_energy.txt',"w") as f:
        f.write("#QM    interpInQMgrid400k    MM\n")
        for i in zip(e0,e1,e2):
            f.write("%f %f %f\n"%i)
    with open('result_force.txt',"w") as f:
        f.write("#QM    interpInQMgrid400k    MM\n")
        for i in zip(fce0,fce1,fce2):
            f.write("%f %f %f\n"%i)
    with open('result_torque.txt',"w") as f:
        f.write("#QM    interpInQMgrid400k    MM\n")
        for i in zip(trq0,trq1,trq2):
            f.write("%f %f %f\n"%i)
    _, axarr = plt.subplots(2, 2,figsize=(8, 8))
    err = rmse(e0,e1)
    q = uqi(e0,e1)
    axarr[0,0].scatter(e0, e1,c='blue')
    axarr[0,0].set_xlabel('qE:rmse=%.4f,uqi=%.4f'%(err,q))
    axarr[0,0].xaxis.set_label_position('top') 

    err = rmse(e0,e2)
    q = uqi(e0,e2)
    axarr[0,1].scatter(e0, e1,c='blue')
    axarr[0,1].scatter(e0, e2,c='red')
    axarr[0,1].set_xlabel('mmE(red):rmse=%.4f, uqi=%.4f'%(err,q))
    axarr[0,1].xaxis.set_label_position('top') 
    
    err1 = rmse(fce0,fce1)
    q1 = uqi(fce0,fce1)
    err2 = rmse(fce0,fce2)
    q2 = uqi(fce0,fce2)
    axarr[1,0].scatter(fce0, fce1,c='blue')
    axarr[1,0].scatter(fce0, fce2,c='red')
    axarr[1,0].set_xlabel('qF:rmse=%.4f, uqi=%.4f\nmmF:rmse=%.4f, uqi=%.4f'%(err1,q1,err2,q2))
    
    err1 = rmse(trq0,trq1)
    q1 = uqi(trq0,trq1)
    err2 = rmse(trq0,trq2)
    q2 = uqi(trq0,trq2)
    axarr[1,1].scatter(trq0, trq1,c='blue')
    axarr[1,1].scatter(trq0, trq2,c='red')
    axarr[1,1].set_xlabel('qT:rmse=%.4f, uqi=%.4f\nmmT:rmse=%.4f, uqi=%.4f'%(err1,q1,err2,q2))
    plt.savefig('corr.png')


if __name__ == '__main__':
    order = 3
    calculator = EFT_calculator(order)
    t0 = time()
    cc = Classical_calculator()
    #acalculator.setup('grid_data.txt')
    calculator.setup('grid.dat')
    #calculator.fill_grid(cc)
    t1 = time()
    print 'took %.1f s to fill the grid' % (t1 - t0)
    test_random_set(sys.argv[1]) # a file with a list of location of QM .log file

