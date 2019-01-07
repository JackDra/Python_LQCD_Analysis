#!/usr/bin/env python

import numpy as np
# import PlotData as jpl
import matplotlib.pyplot as pl
import pandas as pa
from BootStrapping_self import BootStrap



def CreateAndreaTest(nt=64,boot=True):
    # energy_states = np.array([0.2,0.4,0.8,1.6])
    # coeff_list = np.array([[1,0.5,0.25,0.125],
    #                        [0.5,0.5,0.25,0.25],
    #                        [0.1,0.2,0.4,0.8],
    #                        [0.1,1,10,100]]).T

    energy_states = np.array([0.5,0.63,0.8])
    coeff_list = np.array([1,0.1,0.01])
    err_per_start = 0.01
    err_per_end = 1
    def sig_fun(t):
        return err_per_start + (err_per_end - err_per_start)*t/32
    this_nboot = 200
    this_ncfg = 200
    booted = boot
    def exp_fun(val,this_size):
        this_sum = np.sum( coeff_list* np.exp(-energy_states*val))
        rand_list = np.exp(np.random.normal(loc=1,scale=sig_fun(val+1),size=this_size))
        return this_sum * rand_list
        # return np.random.lognormal(mean=this_sum,sigma=sig_fun(val+1)*this_sum,size=this_size)
    values = []
    for it in range(nt):
        if booted:
            this_vals = exp_fun(it,this_ncfg)
            values.append(BootStrap(thisnboot = this_nboot,cfgvals=this_vals))
        else:
            values.append(exp_fun(it,this_ncfg))
    return np.array(values)

if __name__ == '__main__':
    nt = 32
    boot = False
    basedir = '/home/jackdra/LQCD/Scripts/Python_Analysis/TestGraphs/'
    raw_data = CreateAndreaTest(nt=nt,boot=boot)
    if boot:
        pl.errorbar(list(range(nt)),[ival.Avg for ival in raw_data],[ival.Std for ival in raw_data],label='test 2pt',fmt='o')
    else:
        pl.errorbar(list(range(nt)),[np.mean(ival) for ival in raw_data],[np.std(ival) for ival in raw_data],label='test 2pt',fmt='o')
    pl.semilogy()
    pl.savefig(basedir+'/TestCorr.pdf')
    pl.clf()
    if boot:
        eff_mass = [(ival/ivalp1).Log() for ival,ivalp1 in zip(raw_data[:-1],raw_data[1:])]
        [ival.Stats() for ival in eff_mass]
    else:
        eff_mass = [np.log(raw_data[:-1]/raw_data[1:])]
    if boot:
        pl.errorbar(list(range(nt-1)),[ival.Avg for ival in eff_mass],[ival.Std for ival in eff_mass],label = 'effmass',fmt='o')
    else:
        pl.errorbar(list(range(nt)),[np.mean(ival) for ival in eff_mass],[np.std(ival) for ival in eff_mass],label = 'effmass',fmt='o')
    pl.savefig(basedir+'/TestEffMass.pdf')
    pl.clf()
