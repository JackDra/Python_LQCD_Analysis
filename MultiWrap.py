#!/usr/bin/env python

# from pathos.multiprocessing import Pool
from multiprocessing import Pool
from Params import defnProc,defMaxTasks



## runs function over each element of input parameters of listofvars.
## listofvars = [(a1,b1,c1,....) , (a2,b2,c2,....) ...]
## function(a1,b1,c1,...)
## function(a2,b2,c2,...)
## function(a3,b3,c3,...)
## etc...




def DoMulticore(thisfun,thislist,nprocmax = defnProc):
    # print thisfun
    # for ilist in thislist:
    #     print ilist
    if len(thislist) > 1 and nprocmax > 1:
        print('Multicore run, running over ', min(len(thislist),nprocmax) , ' cores')
        thisPool = Pool(min(len(thislist),nprocmax),maxtasksperchild=defMaxTasks)
        output = thisPool.map(thisfun,thislist)
        # thisPool.join()
        thisPool.close()
        return output
    else:
        return [thisfun(iparams) for iparams in thislist]
