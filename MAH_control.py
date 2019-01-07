#!/usr/bin/env python

import sys
from subprocess import Popen
import time
master_ens_list = ['mpi411','mpi570','mpi701','L16','L20','L28','qu_L24']

nblock_list = []
if len(sys.argv)> 1:
    ens_list = []
    for iarg in sys.argv[1:]:
        if 'mpi' in iarg and 'ens' in iarg:
            ens_list.append(['mpi411','mpi570','mpi701'])
        elif 'latspace' in iarg and 'ens' in iarg:
            ens_list.append(['L16','L20','L28'])
        elif 'quenched' in iarg and 'ens' in iarg:
            ens_list.append(['qu_L24'])
        elif 'nblock' in iarg:
            nblock_list.append(iarg)
        elif iarg not in master_ens_list:
            print('Warning, ensemble name not found: ',iarg)
        else:
            ens_list.append(iarg)
else:
    ens_list = master_ens_list

if len(nblock_list) == 0:
    nblock_list = ['nblock1','nblock2','nblock3']

print()
print('Running over ensembles')
print(', '.join(ens_list))
print()

cprocs = []
for iens in ens_list:
    for iblock in nblock_list:
        print('         starting job ',iens,iblock)
        cprocs.append(Popen(['./Mat_Alpha_Hack.py',iens,iblock]))
        time.sleep(5)

for cp in cprocs:
    cp.wait()
