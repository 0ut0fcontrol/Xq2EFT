#!/usr/bin/env python
from eft_calculator842 import EFT_calculator as EFT
import numpy as np
import os, sys
e = EFT()
p = np.array([0, 0, 4.0])
a = np.array([1,0,0])
ang = 0.0
with open(sys.argv[1], "w") as f:
    for i in range(17):
        m = e._hopf2PDB(p,a,ang)   
        ang += i*2*np.pi/16
        f.write("MODEL  %6d\n"%(i))
        f.write(m)
        f.write("ENDMDL\n")
