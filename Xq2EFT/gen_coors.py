#!/usr/bin/env python

from eft_calculator import EFT_calculator


calculator = EFT_calculator()
calculator.setup()
# Please change the following code to whatever needed to generate the input 
# coordinates files
# Please make sure to carry the id number along with the results
for id, coors in calculator.gen_atomic_coors(0, 10): 
    print id, coors

