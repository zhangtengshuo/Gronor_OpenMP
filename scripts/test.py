#!/usr/bin/env python3

import os
import sys

num = len(sys.argv)

print('PYTHON SCRIPT num is ',num,' ',sys.argv[1])

name = sys.argv[1]

comm = 'ls -al ' + name + '.inp'

print('Comm is ',comm)

os.system(comm)
os.system('ls -al '+name+'.inp')

os.system('mpirun -n 3 gronor '+name)

os.system('diff '+name+'.tst '+name+'.ok > '+name+'.dif')

len = os.stat(name+'.dif').st_size

print('Diff size is ',len)
