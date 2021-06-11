#! /usr/bin/env python
from pymican import mican

mican = mican()

output = mican.align('test1.1.pdb', 'test1.2.pdb')

print(output)