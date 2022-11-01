#!/usr/bin/env python3

import os
import sys

name = sys.argv[1]

runfile = name+".run"

runcommand = open(runfile).read().rstrip("\n")

os.system('rm -f '+name+'.dif')
os.system('rm -f '+name+'.tst')
os.system('rm -f '+name+'.cpr')

os.system(runcommand+' gronor '+name)

os.system('diff '+name+'.tst '+name+'.ok > '+name+'.dif')

len = os.stat(name+'.dif').st_size

if len > 0 : sys.exit(1)

sys.exit(0)

