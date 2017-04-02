#!/usr/bin/env python
from eft_calculator842 import EFT_calculator as EFT
import numpy as np
import os, sys
e = EFT()
X = np.array([0, 0, 4.0])
p = 1/2.0
n = -1/2.0
qs = np.array([[p,p,p,p],
               [n,p,p,p],
               [p,n,p,p],
               [p,p,n,p],
               [p,p,p,n],
               [n,n,p,p],
               [p,n,n,p],
               [n,p,n,p]
              ])
qs = np.array([[1.,1.,0.,0.],
               [1.,-1.,0.,0.],
               [1,0.,0.,0.],
              ])
ang = 0.0
with open(sys.argv[1], "w") as f:
    for i in range(len(qs)):
        #if i< 8:
        #    q = qs[i]
        #if i> 8:
        #    q = -qs[i-8]
        q = qs[i]/np.linalg.norm(qs[i])
        m = e._Xq2PDB(X,q)   
        f.write("MODEL  %6d\n"%(i))
        f.write(m)
        f.write("ENDMDL\n")
