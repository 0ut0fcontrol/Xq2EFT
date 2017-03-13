#!/usr/bin/env python
import os,sys
from eft_calculator import EFT_calculator

calculator = EFT_calculator()
calculator.setup()
calculator.grid.save("grid.dat")

