#! /usr/bin/env python

# import mican class
from pymican import mican

# setup mican object
mican = mican()

# calculate mican alignment
output = mican.align('test1.1.pdb', 'test1.2.pdb')

# print result
print(output)