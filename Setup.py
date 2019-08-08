#!/usr/bin/env python

import os
import sys
if sys.version_info[0] < 3:
    raise Exception("This branch (Py3 or offshoot) is built with python3, please use older branch if python2 is needed")

command_list =  'conda install python ipython jupyter numpy pandas scipy matplotlib traitsui traits dill xarray'
os.system(command_list)

command_list =  'conda install -c conda-forge openpyxl uncertainties'
os.system(command_list)

command_list =  'conda update --all'
os.system(command_list)

try:
    from FitFunctions import WriteAllFuns
    WriteAllFuns()
except Exception as err:
    from Params import this_dir
    print('failed writing pickled fit functions, usually file structure is not set up!')
    raise Exception(err)
