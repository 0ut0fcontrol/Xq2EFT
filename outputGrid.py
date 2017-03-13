#!/usr/bin/env python
import os,sys
from eft_calculator import EFT_calculator

calculator = EFT_calculator()
calculator.setup()
# Please change the following code to whatever needed to generate the input 
# coordinates files
# Please make sure to carry the id number along with the results
calculator.grid.save("grid.dat")

