#! /usr/bin/env python
import sys
from pymican import mican

#mican = mican(binary='/Users/sminami/MyGithub/mican/mican')
mican = mican()

output = mican.align('test1.1.pdb', 'test1.2.pdb')

print(output)