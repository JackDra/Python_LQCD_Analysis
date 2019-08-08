#!/usr/bin/env python
import pystan
import pandas as pd
import numpy as np
import matplotlib.pyplot as pl
from Params import this_dir

# C [ism, jsm, t, b]
C_list = []
for ism in range(3):
    C_list.append([])
    for jsm in range(3):
        C_list[-1].append([])
        for ib in range(200):
            this_file = '...'
            t_list = read_file(this_file)
            C_list[-1][-1].append(t_list)
# C_list [ism, jsm, b, t]
C_list = np.array(C_list)
# C_list [ism, jsm, t, b]
C_list = np.swapaxis(C_list,2,3)

Meff = np.log(C_list[:,:,:-1,:]/C_list[:,:,1:,:])
# Meff [ism, jsm, t, b]
for ism in range(3):
    # this_Meff [t, b] [[...],[...]...]
    this_Meff = Meff[ism,ism]
    avg_t = np.avg(Meff[ism,ism],axis=1)
    std_t = np.std(Meff[ism,ism],axis=1)
    pl.errorbar(range(len(avg_t)),avg_t,std_t)



str(34).zfill(5)

for line in file:
    print(line.strip().replace('<','>').split('>'))
    this_val = float(line.strip().replace('<','>').split('>')[2])


C[:,0,1]

import os
for subdir,dirs,files in os.walk(os.path.abspath('./')):
    for ifile in files:
        if ifile[-3:] == '.py':
            print(os.path.join(subdir, ifile))
# The Stan model as a string.
model_string = """
data {
  # Number of trials
  int nA;
  int nB;
  # Number of successes
  int sA;
  int sB;
}

parameters {
  real<lower=0, upper=1> rateA;
  real<lower=0, upper=1> rateB;
}

model {
  rateA ~ uniform(0, 1);
  rateB ~ uniform(0, 1);
  sA ~ binomial(nA, rateA);
  sB ~ binomial(nB, rateB);
}

generated quantities {
  real rate_diff;
  rate_diff <- rateB - rateA;
}
"""

data_list = dict(nA = 16, nB = 16, sA = 6, sB = 10)

# Compiling and producing posterior samples from the model.
stan_samples = pystan.stan(model_code = model_string, data = data_list)
print(stan_samples)
import matplotlib.pyplot as pl
print(stan_samples.plot())
pl.show()


range(10).next()
import numpy as np
this_nboot = 10
from warnings import warn
warn('this is a warning')
val = 6.00

[f't_f{ival:.3}' for ival in np.arange(0,10.00,0.1)]
f'{val:.2f}'
np.arange(0.01,10.00,0.01)
import numpy as np
data = np.arange(25).reshape(5,5)
(data.T + data)/2

import xarray as xr
el = [2,1,3,4]
lev = [1,2,3,4]
rev = [4,3,2,1]
list(zip(*[[iel,ilev,irev] for iel,ilev,irev in sorted(zip(el,lev,rev))]))
vals = {'vals':1}
list(vals.values())
val = 12
f'{int(val)}'
f'{val:04}'
data_empty = xr.DataArray([])
len(data_empty)
this_coords = {}
this_coords['third_dim'] = ['a','b','c']
this_coords['first_dim'] = ['first','second']
this_coords['second_dim'] = ['one','two','four','three']
class TestClass(object):
    pass
import numpy as np
np.mean(np.random.normal(scale=.5,size=100))
this_data = np.array([[[1,2,3,4],[5,6,7,8]],[[1,2,3,4],[5,6,7,8]],[[1,2,3,4],[5,6,7,8]]])
this_shape = this_data.shape
this_data = np.array([TestClass() for i in this_data.flatten()]).reshape(this_shape)
this_data

def myfun(**kargs):
    print(kargs)
    return kargs
len({})
myfun(**myfun(this_input=5))

this_class_data.dtype

import pandas as pa
from copy import copy
this_series = pa.Series([1,2,3])
this_series[::]
print()
print(this_series)
print(this_series.append())
data = xr.DataArray(this_data,coords=this_coords,dims=list(this_coords.keys()),name='name')
data.to_series()
series_data = data.to_series()
series_data_2 = series_data.copy()
series_data.index.names = ['one','two','three']
series_data_2.index = series_data.index
series_data_2
df_data = series_data.copy().to_frame('pie')
series_data.index = series_data_2.index
df_data.loc[:,'pie2'] = series_data
df_data
data = data.stack(first_and_second=['first_dim','second_dim'])
data.unstack()
data = data.stack(all=data.dims)
data
def this_fun(x):
    print()
    print(x)
    print()
    return x+1
data
data.to_pandas()
data.to_pandas()['first']['one'].index
'a' in data.coords[data.coords.dims[-3]].values
data.dtype == float
data.unstack()
def myfun(val):
    print(val.shape)
    print(type(val))
    return [val[0],val[1]]
data
new_data = xr.apply_ufunc(myfun,data,
                          input_core_dims=[['first_dim']],
                          output_dtypes=[float,float],vectorize=True)
this_iter = enumerate(range(10))
ic,ival = next(this_iter)
ic,ival
new_data
data.isel(first_and_second=2)
list(data.coords)[:-2]

print()
data_ds = data.to_dataset()
data_ds.to_array()
pandas_data = data.sel(first_dim='second',second_dim='two').to_dataframe()
pandas_data
pandas_data[pandas_data.columns[-1]].values
len(data)

from BootStrapping import BlockCfgs
import numpy as np
import pandas as pa
data = np.random.normal(loc=0.5,scale=0.25,size=10000)
n_boot = 3
data_blocked = BlockCfgs(data,n_boot)
this_df = pa.DataFrame()
this_df['data'] = data
this_df['blocked_data'] = data_blocked
this_df
this_df.describe()
np.mean(data)
np.mean(data_blocked)


class MyDataArray(xr.DataArray):
    def __init__(self,*args,**kargs):
        super().__init__(*args,**kargs)




data_empty = MyDataArray([])
len(data_empty)
this_coords = {}
this_coords['first_dim'] = ['first','second']
this_coords['second_dim'] = ['one','two','four','three']
data = MyDataArray([[1,2,3,4],[5,6,7,8]],coords=this_coords,dims=list(this_coords.keys()),name='name')
data
len(data)
t0 = 1
f'{t0:.3}'
import pandas as pa
data.attrs['lat_space'] = 0.0980
data.attrs['cfg_list'] = pa.Series([1,2,3,4],index=['one','two','four','three'])
data
test_file = this_dir+'/TestGraphs/test_netcdf.nc'
data.to_netcdf(test_file)

data_read = xr.open_dataarray(test_file)
data_read

import pandas as pa
data = pa.DataFrame({'pie':[1,2,3,4]})
data = pa.Series([1,2,3,4],name='pie')
data.columns[0]

import numpy as np

data = np.arange(24).reshape(2,3,4)
data = xr.DataArray(data)
data.rename({data.dims[0]:'configs'})
data.sel(0=0)
list(map(int,np.zeros(data.ndim)))
data[tuple(map(int,np.zeros(data.ndim)))]
data[(0,0,0)]
len(data.shape)
data.ndim
data = []
data.append([1,2,3])
data.append([4,5])
data = np.array(data)
type(data[tuple(map(int,np.zeros(data.ndim)))])
data.ndim


data = {'one':1}
data.values()

import numpy as np
import dask.dataframe as dd
import pandas as pa
data = pa.read_csv('./Configs/t2_energy_16X32.dat',delimiter=' ')
data
flow_list = np.arange(3.01,9.02,1.0)
out_list,err_list = [],[]
a20rat = (0.0980/0.1215)**2
a28rat = (0.0685/0.1215)**2
for t16 in flow_list:
    for t20 in flow_list:
        for t28 in flow_list:
            out_list.append([t16,t20,t28])


out_list = np.swapaxes(np.array(out_list),0,1)
# err_list = np.swapaxes(np.array(err_list),0,1)
data = pa.DataFrame()
data['t16'] = out_list[0]
data['t20'] = out_list[1]
data['t28'] = out_list[2]
data['t16_rad'] = np.sqrt(8*out_list[0])*0.1215
data['t20_rad'] = np.sqrt(8*out_list[1])*0.0980
data['t28_rad'] = np.sqrt(8*out_list[2])*0.0685
data
# data = dd.from_pandas(data,npartitions=3)
data
data['err1'] = data.apply(lambda x : abs(x['t16'] - a20rat*x['t20']),axis=1)
data['err2'] = data.apply(lambda x : abs(x['t16'] - a28rat*x['t28']),axis=1)
data['comb'] = data['err1']**2 + data['err2']**2
data = data.compute()
es['x'].compute()
data = data.sort_values(['comb'])
data['err1'].values

data['err1'] = err_list[0]
data['err2'] = err_list[1]
this_err_1 = t16 - a20rat*t20
this_err_2 = t16 - a28rat*t28
err_list.append([this_err_1,this_err_2])

data.mean()


import dask
list = [0,1,2,10,3,4,5,6]
vals = [2,4,1,5,6,8,6,7]

list,vals = zip(*[[y,x] for y,x in sorted(zip(list,vals))])

list
vals
import numpy as np
str.zfill()
zfill(10)
n=5
print(f'{n:03}')


from MiscFuns import ODNested

print('abc',end='')

alpha = ODNested()
alpha['one']['two'] = 5
alpha.values()[0].values()[0]

from uncertainties import ufloat

this_val = ufloat(0.0005,0.0025)

print r'${:.1eSL}$'.format(this_val)
print 'a\rightarrow'.replace('\rightarrow','RIGHTARROW')

1.0 == 1.0000000000000001
xvalssize = (5,5)
xdata_1 = np.array([range(xvalssize[1]) for _ in range(xvalssize[0])])
xdata_2 = xdata_1.swapaxes(0,1).flatten()
xdata_1 = xdata_1.flatten()
import pandas as pa
import numpy as np
len({})
len(pa.Series())
len(pa.DataFrame())
type(slice(None))

from XmlFormatting import MakeValAndErr
import numpy as np
np.exp.__name__
this_array = np.ma.array([1,2,3,-1,-2])
type(this_array)
type(np.ma.filled(np.log(this_array),0))
print '\\\\hello world'
print '\\\\hello world'.replace('\\\\','\\')

print MakeValAndErr(21235,1232)
print MakeValAndErr(0.0002123*10000,0.00000123141*10000)+'e{}'.format(np.log10(1/10000.))

[1,2,3,4][0:1]

't10'.index('1')

import re

re.sub('t[1234]','t=\.','t1')

''.split('_')
from XmlFormatting import LegendFmt

LegendFmt('Vector_t10_t_f1.23_ts12_tss13 frac{F_{3}}{2m_{N}}')

leg_string
a = [1,2,3,4,5]
for ic,ia in enumerate(a):
    a[ic] = ia + 1
    a[ic] = a[ic] + 2
a

'Chi1'.lower()

10 <= 5 <= 10
import operator as op

from MiscFuns import MultiSplit
MultiSplit('one_two@three',['@','_'])
a = 'one_two@three'.replace('@','_').split('_')

op.truediv.__name__

import pandas as pa
import pandas
import numpy as np

[1]*.8

s1 = pa.Series()

s1 = s1.set_value('ikey',10)
np.array([1,2]) == [3,4]
s1
n = 5000
##n = 10
values = np.random.uniform(size=(n, 5))
# values = np.arange(1,10)
# values.shape
# for ival in itertools.combinations(values,1):
#     print ival[0]
# construct a pandas.data
# columns = list(zip(*[np.array(['a','a','b','b','b']),np.array(['a','b','c','d','e']),np.array([1,2,1,2,1])]))
# indicies = pandas.multiindex.from_tuples(columns,names=['upper','lower','number'])

list1 = np.array(['C','C','C','C','A','A','A','A','B','B','B','B'])
list2 = np.array(['a','a','b','b','a','a','b','b','a','a','b','b'])
list3 = np.array(['1','2','1','2','1','2','1','2','1','2','1','2'])
columns = list(zip(list1,list2,list3))
indicies = pandas.MultiIndex.from_tuples(columns,names=['upper','lower','number'])
pa.isnull(this_series['pie']).nonzero()[0]
pandas.isnull(this_series['pie']).nonzero()[0]
this_series = pandas.Series(list(values[:12,0]),index=indicies).to_frame(name='pie')
for icol in this_series:
    print(icol)
this_series['pie'] = this_series.apply(lambda x : x['pie'] + 1,axis=1)
this_series.sort_values(by=['pie'],ascending=False)

this_series['pie']
this_series['pie'].index.shape
arr = this_series.values.reshape(*map(len,this_series.index.levels))
arr.transpose(2,0,1).shape
this_series['pie'].iloc[5] = float('NaN')
for this_row,irow in this_series.iterrows():
    print()
    print(this_row)
    print(irow)
    print()
this_series.isnull()
this_series['pie'].items()
any(this_series.isnull().values)
this_series
a=[1,2,3,4]
a[slice(0,-1)]
list(this_series.index).index(('C','a','2'))
this_series['pie'].reset_index()
this_data = this_series['pie'].values.flatten()
this_data = this_data - np.roll(this_data,-1)
this_data.index
this_series.index.names[0]
this_data
this_data[0]-this_data[1]
thisl = []
for (ilower,iupper),iseries in this_series['pie'].groupby(level=('upper','number')):
    print()
    print(ilower,iupper)
    print(iseries)
    print()

this_series
0.788895-0.466798
def this_fun(this_data):
    print 'this_data',this_data

this_series['pie2'] = list(this_series['pie'].values.flatten()+5)
this_series
this_series.index[0][0]
lower_name = this_series.index.names[1]
this_series['pie'].loc[(slice(None),['a','b'],slice(None))]
this_series.reset_index(inplace=True)
this_series['lower'] = this_series['lower'].apply(lambda x : x+'_add').values
this_series.set_index(['upper','lower','number'],inplace=True)

this_series
this_series['pie']['C',:,'1']
this_series.index.levels[2][0]
this_series.loc[('A','a','3')] = pa.Series([5,10],index=['pie2','pie'])
this_series
this_series.loc[:,'pie2'] = 100

pa.Series(list(zip(this_series.loc[:,'pie'],this_series.loc[:,'pie2'])),index=this_series.index)

for ikey,ival in this_series.groupby('lower'):
    print ikey
    print ival
    print

%timeit this_series.at[('C','a','3'),'pie'] = 5
%timeit this_series.at[('C','a','3'),'pie'] = 5
%timeit this_series.loc[:,'pie2'] = this_series.loc[:,'pie']*5
%timeit this_series.at[:,'pie3'] = this_series.at[:,'pie']*5
%timeit this_val = this_series.at[('C','a','3'),'pie']
%timeit this_val = this_series.loc[('C','a','3'),'pie']
this_series

for irow,ival in this_series.iterrows():
    hold_val = ival['pie']*3
    ival['pie2'] = hold_val
this_series

type(this_series.index)

this_series2 = pa.Series([1,2],index=[3,4])
type(this_series2.index)

list(this_series)

this_series['pie'].loc[(slice(None),slice(None),'2')]
('1','c') in this_series.index.swaplevel(0,2)
this_series[['pie']]
list(this_series.index.swaplevel(0,2))

for ikey,ival in this_series['pie'].loc[(slice(None),slice(None),'2')].iteritems():
    print ikey
    print ival
    print

for ikey,ival in this_series['pie'].groupby(['upper','number']):
    print ikey
    print ival
    print
this_series.drop(('C','a','1'),inplace=True)
this_series
this_series**2
this_series['pie'].values * this_series['pie']

this_series.loc[(slice(None),'a',slice(None)),'pie']

(slice(None),'a',slice(None)) in this_series.index
this_series.index.isin(this_series)
this_series.index.get_level_values(('lower','number')


# hash((slice(None),1,2))
#
#
# qcEDM_dat = '/home/jackdra/PHD/CHROMA/TestVar/scratch/dataqcEDM/t8/RC16x32_B1830Kud013825Ks013710C1761-1-002150_t8.dat'
#
# qcEDM_txt = '/home/jackdra/PHD/CHROMA/TestVar/scratch/dataqcEDM/t8/RC16x32_B1830Kud013825Ks013710C1761-1-002150_t8.txt'
#
# my_dat = '/home/jackdra/PHD/CHROMA/TestVar/scratch/cfunsg5/RC32x64_Kud01375400Ks01364000/threept/sm64PoF0D1GMA4tsink13p000sing/RC32x64_B1900Kud01375400Ks01364000C1715-a-004410_xsrc104_k1375400_tsrc0sm64_PoF0D1GMA4tsink13p000singNDer0.3cf'
#
# from ReadBinaryCfuns import holder
# import numpy as np
# holder(qcEDM_dat,cmplx_dtype=np.complex64)
# holder(my_dat,cmplx_dtype=np.complex128)
#
# from ReadBinaryCfuns import ReadFSCfunPickCHROMA
#
# data = ReadFSCfunPickCHROMA([my_dat],['0 0 0'],['g4'],forcent=64)
# import glob
# data_list = []
# file_list = []
# for ifile in glob.glob('/home/jackdra/PHD/CHROMA/TestVar/scratch/dataqcEDM/t8/*'):
#     file_list.append(ifile)
#     data_list.append(ReadFSCfunPickCHROMA([qcEDM_dat],['0 0 0'],['I'],forcent=32,force_cmplx_type=np.complex64))
# file_list[1],data_list[1].data[0][0]
# np.array(data_list).swapaxis(0,1)
# pl.errorbar(indexl,avgl,stdl)
# pl.title(r'$\frac{t_f<W(t_f)W(0.21)>}{<Q(5.01) W(0.21)>}$')
# pl.xlabel(r'$t_f$')
# pl.ylabel(r'Ratio')
# pl.savefig(this_dir+'/TestGraphs/ForMatAndrea.pdf')

# W_021 = W_data['Op_Stats'].loc['t_f0.21','boot']
# Q_5 = Q_data['Op_Stats'].loc['t_f5.01','boot']
# W_t = W_data['Op_Stats'].loc[:,'boot']
# t_list = [untflowstr(ikey) for ikey in W_data['Op_Stats'].index]
#
# result =

# import dask.dataframe

# import torch
# # from __future__ import print_function
#
#
# print(torch.version.cuda)
# x = torch.rand(5, 3)
# # print x
#
#
#
#
#
# y = torch.rand(5, 3)
# # print y
#
#
#
# torch.cuda.is_available()
# # if torch.cuda.is_available():
# x = x.cuda()
# y = y.cuda()
# x + y
# a = {}
#
# a[]
#
#
#
# from  scipy.optimize import curve_fit
#
# np.arange(24).reshape((2,3,4))
#
#
# import pandas
# import pandas as pa
# import numpy as np
# # import matplotlib.pyplot as pl
# from copy import deepcopy
# import glob
# from itertools import product
# import re,os
# import itertools
# # from glob import glob
# import sys
# import dill as pickle
# import pickle
# import cPickle as pickle
# from collections import ordereddict
# from xmlformatting import rm_tsumfitr
# from miscfuns import removeallboots
#
# thisfile = '/home/jackdra/PHD/CHROMA/TestVar/scratch/resultsSCfgs/Cfgs/RC16x32Kud01382500Ks01371000/G2/Baryon_-12-_CPOdd_tsrc0_ism64_jsm64_Full.cfgs.msg'
#
# this_a_file = '/home/jackdra/PHD/CHROMA/TestVar/scratch/resultsSCfgs/Cfgs/RC16x32Kud01382500Ks01371000/G2/Baryon_-12-_CPOdd_tsrc0_ism64_jsm64_Full.cfgs.msg'
#
# files = glob.glob('/home/jackdra/PHD/CHROMA/TestVar/scratch//resultsSCfgs//Cfgs/RC16x32Kud01382500Ks01371000/G2//Baryon_-*-_CPEven_tsrc0_ism64_jsm64_Full.cfgs.msg')
# import pandas as pa
#
# fo_file = '/home/jackdra/PHD/CHROMA/TestVar/scratch/resultsSCfgs/Cfgs/FlowOps/RC16x32_kud1382500_ks1371000_-12-_TopCharge_Full.cfgs.msg'
# cfg_data = pa.read_msgpack(thisfile)
# fo_data = pa.read_msgpack(fo_file)
# cfg_data
#
# cfg_data['C2'].head().apply(lambda x : np.array(x).shape)
# cfg_data.tsink_list = np.arange(32)
# cfg_data.tsink_list
#
# test_file = this_dir+'/TestGraphs/Testmsg.msg'
# import cPickle as pickle
#
# with open(test_file,'wb') as f:
#     pickle.dump(np.arange(32),f)
#     cfg_data.to_msgpack(f)
#
# with open(test_file,'rb') as f:
#     tsink_list = pickle.load(f)
#     data_new = pa.read_msgpack(f)
#
# list_a = np.array([1,2,3,4])
# list_b = np.array([2,4])
#
# [np.where(list_a == ib)[0][0] for ib in list_b]
# list_a.indices(list_b)
# [list_a.index(ib) for ib in list_b]
#
# set(list_b).issubset(list_a)
#
# tsink_list
# data_new
# data_new
# np.array(fo_data['Op'].values[0]).shape
# files
# print files
# this_file = '/home/jackdra/phd/chroma/testvar/scratch/beta62_flow00999.dat'
# from autocorr import autocorrelate
# import numpy as np
# import pandas as pa
#
# sdata = pa.read_csv(this_file).iloc[1:]
# sdata = sdata.astype(float)
# sparam = np.arange(1,2,0.1)
#
# rho_list = []
# w_list = []
# tau_list = []
# tauerr_list = []
# ilist = []
# itau_list = []
# w_ilist = []
# auto_list = []
#
# for iss in sparam:
#     hold_auto = autocorrelate(fun=(lambda x : x,lambda x : [x/x]),data=sdata,sparam=iss)
#     auto_list.append(hold_auto)
#     skey = 's{:.1f}'.format(iss)
#     for iw,wdata in enumerate(hold_auto.gft/hold_auto.gft[0]):
#         rho_list.append(wdata)
#         ilist.append((skey,'w'+str(iw)))
#     for iw,taudata in enumerate(hold_auto.tauw):
#         tau_list.append(taudata)
#         tauerr_list.append(hold_auto.tauerrw[iw])
#         itau_list.append((skey,'w'+str(iw)))
#     w_list.append(hold_auto.wopt)
#     w_ilist.append(skey.format(iss))
#
# indicies = pa.multiindex.from_tuples(ilist,names=['s_parameter','w_parameter'])
# rho_list = pa.series(rho_list,index=indicies)
# indicies = pa.multiindex.from_tuples(itau_list,names=['s_parameter','w_parameter'])
# tau_list = pa.series(tau_list,index=indicies)
# tauerr_list = pa.series(tauerr_list,index=indicies)
# w_data = pa.series(w_list,index=w_ilist)
# auto_list = pa.series(auto_list,index=w_ilist)
#
#
# fmt_data = auto_list.apply(lambda x : x.avg).to_frame('avg')
# fmt_data['std'] = auto_list.apply(lambda x : x.std)
# fmt_data['formatted'] = auto_list.apply(lambda x : makevalanderr(x.avg,x.std))
# fmt_data
# fmt_data.to_csv(this_dir+'/TestGraphs/mathias_data.csv')
#
#
# rho_list.loc[('s1.5',slice(none))]
#
#
# this_info = pa.series()
# this_info['save_file'] = this_dir+'/TestGraphs/mathias_top_w.pdf'
# this_info['title'] = 'test w auto graph'
# # this_info['xlims'] = [0,10]
# # this_info['ylims'] = [0,15]
# import plotdata as jpl
# data_plot = jpl.plotting(plot_info=this_info)
# plot_series = pa.series()
# plot_series['type'] = 'scatter_vary'
# plot_series['x_data'] = 'from_keys'
# plot_series['y_data'] = rho_list
# plot_series['color'] = 'blue'
# plot_series['shift'] = 0.0
# plot_series['key_select'] = ('s1.5',slice(none))
# data_plot.appenddata(plot_series)
# plot_series = pa.series()
# plot_series['type'] = 'vline_vary'
# plot_series['x_data'] = w_data
# plot_series['key_select'] = 's1.5'
# plot_series['shift'] = 0.0
# plot_series['color'] = 'previous'
# data_plot.appenddata(plot_series)
# data_plot.printdata()
# data_plot.plotall()
#
# tauerr_list
#
# this_info = pa.series()
# this_info['save_file'] = this_dir+'/TestGraphs/mathias_top_tauint.pdf'
# this_info['title'] = 'test tau int auto graph'
# # this_info['xlims'] = [0,10]
# # this_info['ylims'] = [0,15]
# data_plot = jpl.plotting(plot_info=this_info)
# plot_series = pa.series()
# plot_series['type'] = 'error_bar_vary'
# plot_series['x_data'] = 'from_keys'
# plot_series['y_data'] = tau_list
# plot_series['yerr_data'] = tauerr_list
# plot_series['color'] = 'red'
# plot_series['shift'] = 0.0
# plot_series['key_select'] = ('s1.5',slice(none))
# data_plot.appenddata(plot_series)
# plot_series = pa.series()
# plot_series['type'] = 'vline_vary'
# plot_series['x_data'] = w_data
# plot_series['key_select'] = 's1.5'
# plot_series['color'] = 'previous'
# plot_series['shift'] = 0.0
# data_plot.appenddata(plot_series)
# # data_plot.loadpickle(defwipe=false)
# data_plot.printdata()
# data_plot.plotall()
#
#
# from xmlformatting import makevalanderr
# # print makevalanderr(test_data.avg,test_data.std)
#
# data = none
# def test_read(this_file):
#     global data
#     file_read = open(this_file,'rb')
#     data = pickle.load(file_read)
#     file_read.close()
#     return data
#
# def test_write(this_file):
#     global data
#     file_write = open(this_file,'wb')
#     pickle.dump(data,file_write)
#     file_write.close()
#
# def test_read_write(this_file):
#     test_read(this_file)
#     test_write(this_file)
#
# this_file = '/home/jackdra/phd/chroma/testvar/scratch/resultschidistdebug/rc32x64kud01375400ks01364000/alpha/pickle/baryon_-a-_cpodd_tsrc0_ism64_jsm64_topcharge_sr_su_fo2_sy.py3p'
#
#
# this_file = '/home/jackdra/phd/chroma/testvar/scratch/resultschidistdebug/rc32x64kud01375400ks01364000/alpha/pickle/baryon_-a-_g5_tsrc0_ism64_jsm64_topcharge_sr_mt_nnq2sub_g5.py3p'
#
# import dill as pickle
# pickle.settings['byref'] = false
# pickle.settings['recurse'] = false
# # # %timeit test_read(this_file)
# # # %timeit test_write(this_file)
# # # %timeit test_read_write(this_file)
#
#
# import dill as pickle
# pickle.settings['byref'] = false
# pickle.settings['recurse'] = true
# # %timeit test_read(this_file)
# # %timeit test_write(this_file)
# # %timeit test_read_write(this_file)
#
# import dill as pickle
# pickle.settings['byref'] = true
# pickle.settings['recurse'] = false
# # %timeit test_read(this_file)
# # %timeit test_write(this_file)
# # %timeit test_read_write(this_file)
#
# import dill as pickle
# pickle.settings['byref'] = true
# pickle.settings['recurse'] = true
# # %timeit test_read(this_file)
# # %timeit test_write(this_file)
# # %timeit test_read_write(this_file)
# #
# # import pickle
# # # %timeit test_read(this_file)
# # # %timeit test_write(this_file)
# # # %timeit test_read_write(this_file)
#
# #
# # data['nnqfull_stats'].head()
# # data['nnqfull_stats'].loc[('p000','t10','t_f6.01','ts32'),'boot'].removeboots()
# # data['nnqfull_stats'].loc[('p000','t10','t_f6.01','ts32'),'boot'].bootvals
# nnqdata = data['nnqfull_stats']
# nnqdata = nnqdata.applymap(removeallboots)
#
# nnqdata.loc[('p000','t10','t_f6.01','ts32'),'boot'].bootvals
# nnqdata_new = removeallboots(nnqdata.loc[('p000','t10','t_f6.01','ts32'),'boot'])
# # print nnqdata_new.bootvals
#
#
# def constant_fit_fun(x,p):
#    return p
#
# value,covar = curve_fit(constant_fit_fun, xdata = np.array([1,2,3,4]),ydata = np.array([3,3,3,3]))
#
# fit_min,fit_max = 2,11
# min_fitr = 4
# fit_list = []
# for ifitmin in range(fit_min,fit_max+1):
#     for ifitmax in range(fit_min,fit_max+1):
#         if ifitmax - ifitmin >= min_fitr:
#             fit_list.append('fitr'+str(ifitmin)+'-'+str(ifitmax))
# print 'debug: all fit ranges:'
# print '\n'.join(fit_list)
#
#
# # %matplotlib inline
# data = np.arange(0,32,0.01)
# a = 0.20
# b_1 = 1
# e_1 = 0.75
# b_2 = -1
# e_2 = 5
# state1 = b_1*np.exp(-e_1*data)
# state2 = b_2*np.exp(-e_2*data)
# this_fun = a + state1 + state2
# pl.plot(data,this_fun )
# pl.xlim(-5,35)
# pl.ylim(0.1,0.22)
# pl.savefig(this_dir+'/TestGraphs/fit_test.pdf')
# range(10)[slice(1,10,3)]
#
# def pie(a):
#     return a
#
# print pie.__name__
#
# data = pa.read_csv('./debug2pt_cfgs.csv')
# test_dict = {}
# test_dict['a'] = 5
# test_dict['b'] = 'cow'
# test_series = pa.series(test_dict)
# test_series
# info = pa.read_csv('./debug2pt_info.txt')
# data
#
# this_file = '/home/jackdra/phd/chroma/testvar/scratch/resultschipickdebug/rc32x64kud01375400ks01364000/alpha/pickle/baryon_-a-_cpodd_tsrc0_ism64_jsm64_topcharge_ce_su_sy.py3p'
#
#
#
# with open(this_file,'rb') as f:
#     data = pickle.load(f)
# data.keys()
# data['nnqfull_cfgs']
#
#
#
#
# thisdir = '/home/jackdra/phd/chroma/testvar/scratch/debugdata/flow/rc32x64_kud01375400ks01364000-a-/'
#
# glob_dir = glob.glob(thisdir+'cfg*')
#
# for iglob in glob_dir:
#     if os.path.isfile(iglob+'/_meson'):
#         os.rename(iglob+'/_meson',(iglob+'_meson'))
#     if os.path.isfile(iglob+'/_adjoint_meson'):
#         os.rename(iglob+'/_adjoint_meson',(iglob+'_adjoint_meson'))
#     if os.path.isfile(iglob+'/_e'):
#         os.rename(iglob+'/_e',(iglob+'_e'))
#
# r'//a/b/c'.replace(r'////',r'//')
# r'//a/b/c'.replace(r'//',r'/') == r'/a//b/c'.replace(r'//',r'/')
#
# # %matplotlib inline
# data = np.arange(0,2.5,0.01)
# state1 = (np.exp(-(5-data)) - np.exp(-data))
# state2 = 10*(np.exp(-10*(5-data)) - np.exp(-10*data))
# pl.plot(data,10 + state1 + state2 )
#
#
#
# try:
#     raise notimplementederror("no error")
# except exception as e:
#     exc_type, exc_obj, exc_tb = sys.exc_info()
#     fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
#     print exc_type, fname, exc_tb.tb_lineno
#     print e
#     sys.exit()
#
# print rm_tsumfitr('rfitr4-8_tsumfitr2-15')
#
# print np.argmax([10,2,3,4])
#
#
#
#
# # import matplotlib.pyplot as pl
# # # %matplotlib inline
# print '\n'.join(pl.style.available)
#
# ## move somewhere?
# params = {'legend.fontsize': 7,
#           'legend.numpoints': 1,
#           'axes.labelsize' : 10,
#           'figure.autolayout': true,
#           'axes.grid': true,
#           'errorbar.capsize': 5,
#           'axes.xmargin':0.01,
#           'axes.titlesize' : 20,
#           'axes.ymargin':0.01}
# pl.rcparams.update(params)
# pl.errorbar(range(10),range(5,15),[1]*10)
#
#
# a = (1,2)
# b = (1,2)
# a == b
#
# a = [1,2,3]
# b = np.array(a)
# all(a == b)
# a
#
# pl.clf()
#
# this_file = '/home/jackdra/phd/chroma/testvar/scratch/topcharge/kud01370000ks01364000/pergf/q_flow_b1.90_ng4000.out'
#
# data = np.loadtxt(this_file)
#
# pl.plot(data[20:,0],data[20:,1]**2)
#
# this_file = '/home/jackdra/phd/chroma/testvar/scratch/topcharge/kud01370000ks01364000/pergf/q_flow_b1.90_ng3500.out'
#
# data = np.loadtxt(this_file)
#
# pl.plot(data[20:,0],data[20:,1]**2)
#
# this_file = '/home/jackdra/phd/chroma/testvar/scratch/topcharge/kud01370000ks01364000/pergf/q_flow_b1.90_ng3700.out'
#
# data = np.loadtxt(this_file)
#
# pl.plot(data[20:,0],data[20:,1]**2)
#
# pl.savefig('/home/jackdra/phd/chroma/testvar/scripts/python_analysis/testgraphs/flowtest.pdf')
#
# [1,2,3,4][:3]
# glob('/home/jackdra/phd/chroma/testvar/scripts/reworkanalysis/*')
#
#
#
# this_list =  np.array([1,2,3,3,2,5,10])
# this_list.index(np.max(this_list))
# np.power(2,3)
# np.roll(np.array([1,2,3,4,5,6,7,8,9,10]),-7)
#
#
# def testfun(a,b):
#     print a,b
#
# igen  = iter([1,2,3,4])
# testfun(igen.next(),igen.next())
#
# class test(object):
#     def __int__(self,a):
#         self.a
#
# test_class = test()
# test_class.a = 5
# print getattr(test_class,'a')
#
# range(1,4)
#
# def fun(*args):
#     print ','.join(args)
# fun('yes',1)
# np.lib.pad([1,2,3,4],(0,6),'constant',constant_values=(0,5))
#
# data = np.arange(32)
# data.size**(1/5.)
# data.flatten().reshape(*[data.size**(1/1.)]*1)
#
# 63/2 == 63/2.
#
# list(product([1,2,3],[4,5,6],[7,8,0]))
# np.roll(np.array([1,2,3,4]),-2)
# np_list = np.arange(24,dtype=np.float64)
# np_list[[1,-2]] = 0,0
# np_list[none]
# min(1,2)
# np.append(np.append('nan',np_list[::-1][:-1]),'nan')
# np.rollaxis(np.arange(24).reshape(*[2,3,4]),1,3).shape
#
# os.path.isdir('./flowopps.py')
#
# re.findall('-.*-','/home/jackdra/phd/chroma/testvar/scratch/topcharge/kud01375400ks01364000-a-/pergf/')
#
# glob.glob('/home/jackdra/phd/chroma/testvar/scratch/topcharge/kud01375400ks01364000-*-/pergf/')
# data = np.loadtxt('/home/jackdra/phd/chroma/testvar/scripts/reworkanalysis/testtorch/testdata.txt',dtype=str)
#
# import pandas as pa
# data[0].loc[['a','b']]
# data = pa.DataFrame([1,2,3],index=['a','b','c'])
#
# len(data)
# data.loc(slice([1,2,3]))
#
# data
# slice([1,2,3])
# data
# import numpy as np
# np.random.choice(10,size=5,replace=False)
# values = np.random.randint(low=0,high= 10,size=10)
# values
# testdata = [bs.bootstrap(2000,name='test_bootstrap'+str(icv),cfgvals=ival,thisdelval=false) for icv,ival in enumerate(values)]
# def make1(testdata):
#     this_series = pa.dataframe(np.array(testdata),columns=['this_series'])
#     this_series['avg'] = pa.series(np.array([idata.avg for idata in this_series['this_series']]),name='avg')
#     this_series['std'] = pa.series(np.array([idata.std for idata in this_series['this_series']]),name='std')
#
# def make2(testdata):
#     thislist = []
#     thisliststd = []
#     for ikey,ival in this_series.iterrows():
#         thisval = ival['this_series']
#         thislist.append(thisval.avg)
#         thisliststd.append(thisval.std)
#     this_series['avg_v2'] = pa.series(np.array(thislist),name='avg_v2')
#     this_series['std_v2'] = pa.series(np.array(thislist),name='avg_v2')
#
#
# thisseries = pa.series(values[0],index=np.arange(len(values[0]))[::-1])
# isinstance(thisseries.index,pa.index)
# for iseries in thisseries:
#     print iseries
# thisdf = pa.dataframe(index=np.arange(len(values[0]))[::-1])
#
# thisdf['pie'] = thisseries
# thisdf
# thisdf['str'] = map(str,thisseries)
# thisdf
# thisdf.sort_values('str').reset_index(drop=true)
#
# a = pa.series()
# if a is not none:
#     print a
# thislist = [thisdf,deepcopy(thisdf)]
# for ilist in thislist:
#     del ilist['thisseries']
# thislist[1]
# newseries = pa.series(values[1],index=thisdf.index)
# thisdf['thisseries'][4]
# thisdf['test'] = [np.nan]*5
# 6 in thisdf['thisseries'].index
#
# thisseries = pa.series(np.nan,index=[1,2,3,4])
# thisseries.index
# thisdf2 = pa.dataframe()
#
# # # %timeit make1(testdata)
# #
# # # %timeit make2(testdata)
#
#
#
# # generate some data
# n = 5000
# ##n = 10
# values = np.random.uniform(size=(n, 5))
# # values = np.arange(1,10)
# # values.shape
# # for ival in itertools.combinations(values,1):
# #     print ival[0]
# # construct a pandas.data
# # columns = list(zip(*[np.array(['a','a','b','b','b']),np.array(['a','b','c','d','e']),np.array([1,2,1,2,1])]))
# # indicies = pandas.multiindex.from_tuples(columns,names=['upper','lower','number'])
#
# columns = list(zip(*[   np.array(['C','C','C','C','A','A','A','A','B','B','B','B']),
#                         np.array(['a','a','b','b','a','a','b','b','a','a','b','b']),
#                         np.array(map(str,[1,2,1,2,1,2,1,2,1,2,1,2]))]))
# indicies = pandas.MultiIndex.from_tuples(columns,names=['upper','lower','number'])
# indicies.get_level_values('upper')
# indicies.get_loc_level((slice(None),'a',slice(None)))[1]
#
# for ival in indicies[(slice(None),'a',slice(None))]:
#     print ival
#
# indicies.get_level_values()
#
# indicies
#
# this_series = pandas.series(list(values[:12,0]),index=indicies).to_frame(name='pie')
#
# this_series.sort_values('pie')
# this_series
# for ikey,ival in this_series.groupby(level='upper'):
#     print ikey
#     print ival
# this_series['a']
#
# this_series
#
# for icol,ival in this_series.iteritems():
#     print icol
#     for ikey,keyval in ival.iteritems():
#         print ikey, ival
#
# thisdf.loc[('a','a','1'),'pie'] = this_series.loc[('a','a','1'),'pie']
# thisdf.loc[:,'pie'] = this_series.loc[:,'pie']
# thisdf
#
# this_series = this_series.assign(**{'pie2': this_series['pie'].values})
#
# this_series['pie2'] = this_series['pie']
# this_series['pie'][('a','a','4')] = 5
#
# this_series.iteritems
#
# this_series.reorder_levels([1,0,2])
# ('a','c') in this_series.reorder_levels([1,0,2])['pie'].index
# for iindex in this_series.index:
#     print iindex
#
# test_series = pandas.dataframe(columns=['a','b'],index=[0])
# test_series['a'][0] = [1,2,3]
# test_series
# list(this_series.ix[:,0]['a','a'].index)
#
# len(this_series.index)
#
# np.array(['a','b','c']) == 'c'
#
# this_value = this_series['pie']['a','a','1']
#
# this_where = np.where(this_series['pie']['a'].values == 5)
#
# this_where[0]
# this_series.size
# sum([1,2])
# list(this_series.iterrows())[0][1]['pie']
# pa.concat([this_series['pie']['a'],this_series['pie'][:,'b']])
#
# this_series.drop(('a'))
#
# this_series['pie'][:,:,'1']
# this_series['pie']['a','v','1'] = 1
# this_series['pie']
# ('a','3') in this_series['pie'].swaplevel(1,2)
#
# list(this_series.iterrows())[0][1]
# list(this_series['pie']['a'].groupby(level=('number','number')))[0][1]
# this_series + 1
# pa.series(list(product(this_series.values,this_series.values)),list(product(this_series.index,this_series.index)))
# ('c') in this_series.index
# this_series.iiloc[5]
# this_series.index.get_level_values('upper')
#
# true in [4,2,3]
#
# this_series[0]['c','a']
#
# test_series = pa.series(index=indicies)
# test_series
# test_series['a','a','1'] = 5
#
# test_series[('c','a')] = 5
# np.where(np.array([1,2,3,4]) == 6)[0]
# this_lambda = lambda x : [ix for ix in x]
#
# pa.series([','.join(iindex) for iindex in this_series.index],index=this_series.index)
# map(tuple,this_series.index.levels)
# this_series.values
# this_series.index.levels
# this_series.values[10]
# this_series.iloc[np.where(this_series.values == this_series.values[10])[0][0]]
# for ival in zip(this_series.iteritems(),this_series.values):
#     print ival
# this_series.reshape((3,2,2))
# np.pad(this_series.values,((0,1,2),(0,1,2)), 'constant', constant_values=np.nan)
# this_series.tolist()
# this_series[this_series < 0.5]
# this_series.to_frame().values.shape
# dim_len = this_series.index
# len(dim_len[0])
# this_series.values.reshape(dim_len,this_series.size/dim_len).shape
# this_series = this_series.to_frame(name='values')
# this_series.size
# isinstance(this_series.index,pa.multiindex)
# this_series = this_series.sort_values('values')
# this_series.values
# this_series.sort_index(level='upper',sort_remaining=false)
# this_series['values']['a','a',:].astype(str).values
# this_series['b']
#
# # # %timeit this_series.sample(10, replace=true,random_state=10000).mean()
# # # %timeit this_series.sample(10, replace=true).mean()
#
# this_gen = itertools.product(enumerate([1,2,3]),[4,5,6])
# list(this_gen)
# for igen in this_gen:
#     print igen
# this_df = this_series.to_frame(name='values')
# this_df['values'].size
# this_df['backwards'] = this_series.iloc[::-1]
# this_df.apply(lambda x: pa.series([x['values']-1],index=['res']),axis=1)
# this_df[['values','backwards']]
# this_df.index.levels[0]
# this_df.values.shape
# this_df.unstack().unstack()
#
#
#
# this_np = len(np.arange(24).reshape(2,3,4))
# this_np
# this_np = np.arange(450,0,-2)
# for inp in this_np:
#     inp += 1
# np.ones(100)
# this_np
# np.zeros(5)
#
# list(list(this_series.groupby(level='upper'))[0][1].groupby(level='lower'))[0][1]
# this_series['a','a']-np.roll(this_series['a','a'],1)
# this_series['a','a']
# list(this_series['a'].iteritems())
# list(this_series.iteritems())
# ('a','a',1) in this_series
#
# up_size = this_series.index.levels[0].size
# down_size = this_series.index.levels[1].size
# numb_size = this_series.index.levels[2].size
# for icount,((iup,idown,inumb),ival) in enumerate(this_series.iteritems()):
#     icup = list(this_series.index.levels[0]).index(iup)
#     icdown = list(this_series.index.levels[1]).index(idown)
#     icnumb = list(this_series.index.levels[2]).index(inumb)
#     print icup,iup,icdown,idown,icnumb,inumb,ival
# this_series.index.get_level_values('upper')
# this_series.iloc[1:3]
# this_series + np.array([1,2,3,4,5])
#
# np.array([1,2,3]).tolist()
# [1,2,3] == [1,2,4]
# this_df = this_series.to_frame(name='first')
# this_df['backwards'] = this_series.iloc[::-1]
# this_df.ndim
#
# for ikey,ival in this_df.iterrows():
#     print ikey
#     print ival
#     print
#
#
# this_df_numpy = this_df.reindex(list(product(*map(tuple, this_df.index.levels)))).values\
#      .reshape(map(len, this_df.index.levels))
#
# len(this_df.index.get_level_values('upper'))
# any(np.array([1,2,3,4]) > 0)
# any(this_df > 0)
# [this_df[icol].mean() for icol in this_df]
# 'first' in this_df
# str((1,2))
# list(this_series.swaplevel(1,2)['a','1'].iteritems())
# list(this_series.swapaxis(1,2)['a','1'].iteritems())
# np.cos(this_df)
#
# this_df.xf
# this_df['first'].index
#
# columns_s = list(zip(*[ np.array(['a','a','b','b','c','c']),np.array(['a','b','a','b','a','b'])]))
# indicies_s = pandas.multiindex.from_tuples(columns_s,names=['upper','lower'])
# this_series_s = pandas.series(np.arange(6),index=indicies_s)
# this_df_s = this_series_s.to_frame(name='first')
# this_df_s
# list(this_df.groupby(level='upper'))
# this_df['first'][:,:,1] + this_df_s['first']
# this_df_s + this_df.loc[:,:,1]
#
# this_df = this_df.append(this_df,ignore_index=true)
# list(this_df.iterrows())[0][1]['values']
# this_df.xs('a',level='upper')
# def thisfun(val):
#     return val*2
# this_df.apply(thisfun)
# list(this_series.iteritems())
# np.log(this_series)
# this_df[1] = this_series
#
# for ival in this_df.itertuples():
#     print ival[0],ival['upper']
#
# # printlist(this_series.groupby('a'))
# this_series = this_series.to_frame(name='hi')
# this_series = pa.series([1,2,3]).to_frame(name='hi')
# this_series
# this_series['bye'] = pa.series([4,5,6])
# () in this_series.swaplevel(1,2)
# this_series['a','avg'] = 5
# this_series.keys()
#
# thisdict = {'a':[1,2,3],'b':[2,3,4,5]}
# thislist = map(str,[1,2,3,4,5])
# series_from_dict = pa.series(thislist,index=['a','b','c','d','e'])
# series_from_dict.str.replace('[1,2]','10').to_dict
# ordereddict(list(series_from_dict.iteritems()))['a']
# df_from_dict = pa.dataframe.from_items([ (k,pa.series(v)) for k,v in thisdict.iteritems() ])
# df_new = pa.dataframe(thisdict['a'],columns=['a'])
# df_new['b'] = pa.series(thisdict['b'])
# df_new.iloc[:,1]
# df_new.iloc[:,1][df_new.isin(df_new.iloc[:,0].values).iloc[:,1].values]
# df_from_dict['b'][df_from_dict.isin(df_from_dict['a'].values)]
# pa.series(df_from_dict.keys()).str.replace('[1,2]','')
# series_from_dict
# reversed(list(this_series.iteritems()))
# zip(this_series.values,this_series.values)
# for key,ilist in this_series.iteritems():
#     print isinstance(key,tuple)
#     print key
#     print ilist
#     # print np.mean(ilist)
# series2 = pa.series(np.array([1,2,3,4]),index=['a','b','c','d'])
# for ikey,ival in series2.iteritems():
#     print isinstance(ikey,tuple)
#     print ikey,ival
# this_series[('a',)]
#
# df = pa.dataframe()
# df['col1']
# df['col1'][('a','b')] = 5
# df['col1'][('b','b')] = 5
# df
#
# df['a']
# # # %matplotlib inline
# from pandas.tools.plotting import lag_plot
# lag_plot(df['a'])
# df['a'].autocorr(1)
#
# list(df['a'].iteritems())
# df
# df.index
# import sqlalchemy as sql
# # import sqlite3
# db = sql.create_engine('sqlite:///./testgraphs/test_sql.db')
# df.to_sql('test_data',db)
#
# # # %timeit testsql = pandas.read_sql('select * from "test_data" where "b" < "0.2"',db)
# testsql = pandas.read_sql('select * from "test_data" where "b" < "0.2"',db)
# testsql.head()
#
# newseries = pandas.series(df['a'])
# newseries.index
# newseries.iteritems()
# newseries.to_csv(this_dir+'/TestGraphs/test.cvs')
#
# # bootstrap
# # # %timeit df.iloc[np.random.randint(n, size=n)]
#
# thislist = []
# for ival in xrange(df.a.size):
#     thislist.append(df.a.sample(100,random_state=ival).mean())
#
# df['aboot'] = thislist
# df.aboot.head()
#
# abs(df.a.std() - df.aboot.std())
#
# # # %timeit df.merge(pandas.dataframe(index=np.random.randint(n, size=n)), left_index=true, right_index=true, how='right')
#
# # # %matplotlib inline
#
# # generate some data
# n = 5000
# values = np.random.uniform(size=(n, 5))
# # values = np.random.uniform(size=200)
#
# # construct a pandas.dataframe
# columns = ['a', 'b', 'c', 'd', 'e']
# df = pandas.dataframe(values, columns=columns)
#
# df.index
# df_sample = df.merge(pandas.dataframe(index=np.random.randint(n, size=n)), left_index=true, right_index=true, how='right')
#
# df_series = df['a']
# df_series
#
# df_row_series = df.loc[0]
# df_sum = df_series.cumsum()
# df_sum.describe()
# df_sum.max()
#
# df_sum.hist(histtype='stepfilled', label='myname')
# df_sum.plot.kde()
# pl.plot()
# pl.legend()
# pl.savefig(this_dir+'/TestGraphs/pandas.pdf')
#
#
# np.random.random_integers
#
# type(np.float64(np.nan))
#
#
# values = np.arange(10)[::-1].reshape(10)
# mydf = pandas.dataframe(values, columns=['thisval'])
# # %timeit mydf['cfgvals'] = pandas.series(values, name='cfgvals')
# mydf.cfgvals
#
# mydf = pandas.series(np.nan, index=range(len(values[:, 0])))
# mydf.fillna(values)
# np.nan == np.nan
#
#
# # p000 0  bootval= -3.1909569724e-23  bootstd 5.32055487489e-23 diff= -0.127405228631
# # p000 1  bootval= -1.2831697362e-23  bootstd 1.02615667316e-23 diff= -0.0682895920572
# # p000 2  bootval= -3.89243152583e-24 bootstd 4.80872997196e-24 diff= -0.138210463085
# # p000 3  bootval= -1.39928015966e-24 bootstd 2.28868513982e-24 diff= -0.197884607429
# # p000 4  bootval= -1.04647056197e-24 bootstd 1.11936073721e-24 diff= -0.136566188174
# # p000 5  bootval= -5.47401170543e-25 bootstd 5.69872791818e-25 diff= -0.11681324904
# # p000 6  bootval= -2.75258626007e-25 bootstd 2.92908031078e-25 diff= -0.117093987757
# # p000 7  bootval= -2.16764463102e-25 bootstd 1.55171319006e-25 diff= -0.0712928149594
# # p000 8  bootval= -9.94056986228e-26 bootstd 8.23185798879e-26 diff= -0.0934999845112
# # p000 9  bootval= -6.34141465689e-26 bootstd 4.33752300737e-26 diff= -0.0795283228745
# # p000 10 bootval= -1.21953764187e-26 bootstd 2.32561229394e-26 diff= -0.080582596618
# # p000 11 bootval= -9.22955613094e-27 bootstd 1.29396602323e-26 diff= -0.0476167691032
# # p000 12 bootval= -4.55825513177e-28 bootstd 7.57203496659e-27 diff= -0.483434538767
# # p000 13 bootval= 2.96781939208e-27  bootstd 4.35477240732e-27 diff= 0.0283018486612
# # p000 14 bootval= 1.97435254148e-27  bootstd 2.49044461318e-27 diff= 0.00343515937686
# # p000 15 bootval= 1.9354300106e-28   bootstd 1.53456401081e-27 diff= 0.41966933376
# # p000 16 bootval= -2.11840689647e-28 bootstd 9.22734381653e-28 diff= -0.0481331247429
# # p000 17 bootval= -5.20383232088e-28 bootstd 5.82827825936e-28 diff= -0.0192494430265
# # p000 18 bootval= -4.23383041485e-28 bootstd 3.95077171372e-28 diff= -0.0418847865625
# # p000 19 bootval= -2.53927855819e-28 bootstd 2.62972282883e-28 diff= -0.0106836860991
# # p000 20 bootval= -1.92398838753e-28 bootstd 2.42474449148e-28 diff= -0.0056207841866
# # p000 21 bootval= -1.0474394045e-28  bootstd 2.60497023556e-28 diff= -0.153244643657
# # p000 22 bootval= -1.0448864923e-28  bootstd 3.3973685549e-28  diff= -0.352459465871
# # p000 23 bootval= -7.86329989956e-29 bootstd 4.66722912738e-28 diff= -0.365537931543
# # p000 24 bootval= -1.90796563207e-28 bootstd 6.74454907891e-28 diff= -0.0186459443755
# # p000 25 bootval= 3.33695633239e-28  bootstd 1.14835827035e-27 diff= 0.0696299038905
# #
# #
# #
# #
# #
# #
# #
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# print 4.06544602621e-24/-3.1909569724e-23
# print 8.76271378252e-25/-1.2831697362e-23
# print 5.37974763713e-25/-3.89243152583e-249
# print 2.76896005078e-25/-1.39928015966e-24
# print 1.42912495685e-25/-1.04647056197e-24
# print 6.39437092595e-26/-5.47401170543e-25
# print 3.22311301837e-26/-2.75258626007e-25
# print 1.54537487577e-26/-2.16764463102e-25
# print 9.29443128156e-27/-9.94056986228e-26
# print 5.04322072314e-27/-6.34141465689e-26
# print 9.82735098553e-28/-1.21953764187e-26
# print 4.39481643212e-28/-9.22955613094e-27
# print 2.20361796721e-28/-4.55825513177e-28
# print 8.39947752883e-29/2.96781939208e-27
# print 6.78221564609e-30/1.97435254148e-27
# print 8.12240623088e-29/1.9354300106e-28
# print 1.01965543404e-29/-2.11840689647e-28
# print 1.0017087378e-29/ -5.20383232088e-28
# print 1.77333083268e-29/-4.23383041485e-28
# print 2.7128855034e-30/ -2.53927855819e-28
# print 1.0814323504e-30/ -1.92398838753e-28
# print 1.60514478295e-29/-1.0474394045e-28
# print 3.68280134972e-29/-1.0448864923e-28
# print 2.87433438039e-29/-7.86329989956e-29
# print 3.55758210459e-30/-1.90796563207e-28
# print 2.32351948711e-29/3.33695633239e-28
