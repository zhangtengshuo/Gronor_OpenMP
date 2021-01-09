#!/usr/bin/env python3

import os
import sys

name = sys.argv[1]

runfile = name+".run"

runcommand = open(runfile).read().rstrip("\n")

os.system('rm '+name+'.dif')

os.system(runcommand+' '+name)

os.system('diff '+name+'.tst '+name+'.ok > '+name+'.dif')

len = os.stat(name+'.dif').st_size

if len > 0 : sys.exit(1)

sys.exit(0)

