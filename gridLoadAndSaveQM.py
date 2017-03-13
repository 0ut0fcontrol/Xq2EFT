#!/usr/bin/env python2
import numpy as np
from eft_calculator import EFT_calculator, Water
import tools
import sys
from time import time

qmLogList = sys.argv[1]
calculator = EFT_calculator()
t1 = time()
calculator.setup()
t2 = time()
calculator.fill_with_QM(qmLogList)
t3 = time()
calculator.grid.save("grid.dat")
t4 = time()
print('took %.1fs, %.1fs, %.1fs s to setup, fill and save'%(t2-t1, t3-t2, t4-t3))
