#!/usr/bin/env python
import numpy as np
from tools import *
vecter = (1,1,1)
angle = 0
print("before hopf2q:")
print(vecter,angle)
print('\n')
a = hopf2q(vecter,angle)
print("the q after hopf2q")
print(a)
print('\n')
b = q2hopf(a)
print("after q2hopf back")
print(b)
print('\n')
c = hopf2q(b[0],b[1])
print("covert back to q")
print(c)
print('\n')

