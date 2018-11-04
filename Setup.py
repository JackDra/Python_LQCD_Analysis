#!/usr/bin/env python

import os
import sys
if sys.version_info[0] < 3:
    raise Exception("This branch (Py3 or offshoot) is built with python3, please use older branch if python2 is needed")

command_list =  'conda install python ipython jupyter numpy pandas scipy matplotlib traitsui traits dill'
os.system(command_list)

command_list =  'conda install -c conda-forge openpyxl uncertainties'
os.system(command_list)

command_list =  'conda update --all'
os.system(command_list)
