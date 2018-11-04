#!/usr/bin/env python

import sys,os

for iext in ['xml','p','pdf']:
    os.rename(sys.argv[1]+iext,sys.argv[2]+iext)
