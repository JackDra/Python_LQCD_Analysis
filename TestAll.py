#!/usr/bin/env python

from PredefFitFuns import C2OSFAntiper,C2OSFAntiperDer
from Params import this_dir

temp = 'pie'
this_file = this_dir+'/TestGraphs/flow_plot_test.py3p'
with open(this_file,'rb') as f:
    with open(this_file,'wb') as g:
        pickle.dump(g,temp)
        data = pickle.load(f)

from PlotData import Test_Papers_FF,Plot_E_Jangho,Plot_one_on_x,Plot_one_on_x2
plot_one_on_x = Plot_one_on_x()
plot_one_on_x2 = Plot_one_on_x2()
plot_jangho = Plot_E_Jangho()
plot_FF = Test_Papers_FF()


from PlotData import Test_Papers_FF,Test_Paper_Single,Test_Plot,Plotting
this_plot = Plotting(plot_info = {'save_file':this_dir+'/TestGraphs/Paper_FF.pdf'})
this_plot.LoadPickle()
this_plot.
plot_FF = Test_Plot()
aspect = plot_FF.GetAspectRatio()
plot_FF.SetAspectRatio(1)
plot_FF.PlotAll()
plot_FF.ShowFigure()
print(aspect)
plot_FF = Test_Paper_Single()
plot_FF = mpi701_plot()

# %load_ext autoreload
# %autoreload 2
import numpy as np
import pandas as pa
from collections import OrderedDict
import pickle
# import dill as pickle



from FitFunctions import TestFit,TestComb
data = TestFit()
data2 = TestFit(low=-10,high=20)
datacomb = TestComb()

from FitFunctions_TF import TestFit_TF
import numpy as np
data = TestFit_TF()
data.data_bs
data.data
data.data_bs.sum(axis=0)/len(data.data_bs.index)
data2 = TestFit_TF(low=-10,high=20)
data2.data_bs
data2.data
data2.data_bs.sum(axis=0)
print()
print(data2.loss_values[-1])
# datacomb = TestComb()



import MomParams as mp
new_df,old_df,calc_df,phys_df = mp.ConvertQsqrd()
phys_df
calc_df
new_df
new_df-old_df
calc_df+new_df
old_df
from Params import defInfo
defInfo['kud'] = 1375400
data = mp.TestMom(defInfo)

data
from XmlFormatting import MakeValAndErr
val = data.hbarcdivlat*0.18903
err = data.hbarcdivlat*0.00079
print(MakeValAndErr(val,err))
import SetsOfFits as sff
data,fit_data = sff.TestSetOfFits()
com_data,rand_list = fit_data.GetChiVariation()


import SetsOfFits as sff
data,fit_data = sff.TestSetOfFits_2D()
fit_data.Fit_Stats_fmt
# print com_data['Combine_Fits']['Par2']

import FFSolve as ff
internal_dir = '/home/jackdra/LQCD/Results/DebugResults/'
ff.FixQsqrd_All(this_dir)

dataFFP = ff.TestFFSolve(DicvefWipe=False,q2list=[1,4],DS='Proton')
dataFFN = ff.TestFFSolve(DefWipe=False,q2list=[1,4],DS='Neutron')
datacomb = dataFFP / dataFFN

datacomb.FitFFs()
datacomb.FF_Fit_Stats
dataFFFlow = ff.TestFFFlowSolve(DefWipe=False,q2list=[3])
this_file = '/home/jackdra/PHD/CHROMA/TestVar/scratch//resultsSCfgs//RC32x64Kud01375400Ks01364000/FormFactors/doub_qmax0_ppmax4//Pickle/GeGm_RFmin3_BestChi.py3p'

this_file = this_dir+'/TestGraphs/testFF.py3p'
with open(this_file,'rb') as f:
    test_data = pickle.load(f)

for ikey,ival in test_data.items():
    print()
    print(type(ival), ikey)
    print(ival)
    with open(this_dir+'/TestGraphs/test.py3p','wb') as f:
        pickle.dump(ival,f)

import RatioCorrelators as rat
dataRat2 = rat.TestRat2(DefWipe=False)
dataRatFO2 = rat.TestRatFO2(DefWipe=False)

dataRat = rat.TestRat(DefWipe=False)
dataRatFO = rat.TestRatFO(DefWipe=False)
# dataRatFOFull = rat.TestRatFOFull(DefWipe=False)

datacomb = dataRat2 + (dataRat*1.23)

datacomb.Stats()


datacomb.Stats()

datacomb.Rat_Stats['boot']
dataRat2.Rat_Stats['boot']
datacomb.ReduceIndicies()

import ThreePtCorrelators as t3pc
dataNJNQFull = t3pc.TestNJNQFull(DefWipe=True)
dataNJNQFull.NJNQFull_Stats


dataNJNQ = t3pc.TestNJNQ(DefWipe=True)
dataNJNQ.FO.FlowFun
dataNJNQ.Write()

data3pt = t3pc.Test3pt(DefWipe=True)
data3pt.Read()
data3pt.C3_cfgs
data3pt.C3_Stats


import FormFactors as FF
dataFF = FF.TestFF(DefWipe=False,this_curr='Scalar',qsqrdlist=[0,1,2,3,4])
dataFF.GetRedFFParms()
for ikey in dataFF.itemsParListRed():
    print(ikey)

dataFF.paramRed_Stats
this_file = '/home/jackdra/PHD/CHROMA/TestVar/scratch/resultsSCfgs/RC32x64Kud01375400Ks01364000/Rat/TopCharge/P3_Neutron_pp000/Pickle/FFRats_P3_g4_Top_q1-11_TopCharge.py3p'

import pickle as pickle

rat_data = pickle.load(open(this_file,'rb'))

rat_data['name']

for icoeff,ikey in zip(dataFF.combine_keys_forRF,dataFF.coeffs_forRF):
    print('')
    for jcoeff,jkey in zip(icoeff,ikey):
        print(jcoeff,jkey)





import TwoPtCorrelators as tpc
import os
pre_folder = '/home/jackdra/LQCD/Results/DebugResults/Configs/'
post_folder = '/home/jackdra/LQCD/Results/DebugResults/Configs_Nucleon/'
try:
    os.mkdir(post_folder)
except Exception as err:
    pass
cfg_list = ['RC16x32Kud01382500Ks01371000/G2/Baryon_-12-_CPEven_tsrc0_ism64_jsm64.cfgs.msg',
            'RC20x40Kud01370000Ks01364000/G2/Baryon_-1234-_CPEven_tsrc0_ism64_jsm64.cfgs.msg',
            'RC28x56Kud01356000Ks01351000/G2/Baryon_-12-_CPEven_tsrc0_ism64_jsm64.cfgs.msg',
            'RC32x64Kud01370000Ks01364000/G2/Baryon_-b-_CPEven_tsrc0_ism64_jsm64.cfgs.msg',
            'RC32x64Kud01372700Ks01364000/G2/Baryon_-b-_CPEven_tsrc0_ism64_jsm64.cfgs.msg',
            'RC32x64Kud01375400Ks01364000/G2/Baryon_-ab-_CPEven_tsrc0_ism64_jsm64.cfgs.msg']

cfg_list_out = ['RC16x32Kud01382500Ks01371000/MOMNUMB/Nucleon_CFGNUMB.txt',
                'RC20x40Kud01370000Ks01364000/MOMNUMB/Nucleon_CFGNUMB.txt',
                'RC28x56Kud01356000Ks01351000/MOMNUMB/Nucleon_CFGNUMB.txt',
                'RC32x64Kud01370000Ks01364000/MOMNUMB/Nucleon_CFGNUMB.txt',
                'RC32x64Kud01372700Ks01364000/MOMNUMB/Nucleon_CFGNUMB.txt',
                'RC32x64Kud01375400Ks01364000/MOMNUMB/Nucleon_CFGNUMB.txt']
in_files = [pre_folder + icfg for icfg in cfg_list]
out_files = [post_folder + icfg for icfg in cfg_list_out]
for iin,iout in zip(in_files,out_files):
    try:
        os.mkdir(iout.replace('MOMNUMB/Nucleon_CFGNUMB.txt',''))
    except Exception as err:
        pass
    tpc.ReformatC2(iin,iout)


import FlowOpps as fo
import os
pre_folder = '/home/jackdra/LQCD/Results/DebugResults/Configs/FlowOps_NewForm/'
post_folder = '/home/jackdra/LQCD/Results/DebugResults/Configs_Q/'
post_folder_picktf = '/home/jackdra/LQCD/Results/DebugResults/Configs_Q_picktf/'
try:
    os.mkdir(post_folder)
except Exception as err:
    pass
try:
    os.mkdir(post_folder_picktf)
except Exception as err:
    pass
cfg_list = ['RC16x32_kud1382500_ks1371000_-12-_TopCharge.cfgs.msg',
            'RC20x40_kud1370000_ks1364000_-1432-_TopCharge.cfgs.msg',
            'RC28x56_kud1356000_ks1351000_-12-_TopCharge.cfgs.msg',
            'RC32x64_kud1370000_ks1364000_-b-_TopCharge.cfgs.msg',
            'RC32x64_kud1372700_ks1364000_-b-_TopCharge.cfgs.msg',
            'RC32x64_kud1375400_ks1364000_-ab-_TopCharge.cfgs.msg']
in_files = [pre_folder + icfg for icfg in cfg_list]

cfg_list_out = ['RC16x32Kud01382500Ks01371000.txt',
                'RC20x40Kud01370000Ks01364000.txt',
                'RC28x56Kud01356000Ks01351000.txt',
                'RC32x64Kud01370000Ks01364000.txt',
                'RC32x64Kud01372700Ks01364000.txt',
                'RC32x64Kud01375400Ks01364000.txt']
out_files = [post_folder + icfg for icfg in cfg_list_out]

# for iin,iout in zip(in_files,out_files):
#     try:
#         os.mkdir(iout.replace('Q_CFGNUMB.txt',''))
#     except Exception as err:
#         pass
#     fo.ReformatFO(iin,iout)

picked_tf = ['t_f3.01',
             't_f5.01',
             't_f9.01',
             't_f6.01',
             't_f6.01',
             't_f6.01']
out_files = [post_folder_picktf + icfg for icfg in cfg_list_out]
for iin,iout,itf in zip(in_files,out_files,picked_tf):
    this_out = iout.replace('.txt','_'+itf+'.txt')
    fo.ReformatFO_picktf(iin,this_out,itf)

data2pt = tpc.TestFlowCorr(DefWipe=True)

data2pt = tpc.TestNNQFullCorr(DefWipe=True)



data2pt = tpc.TestTPCorr(DefWipe=True)
data2pt = tpc.TestNNQCorr(DefWipe=False)

import operator as op
import RatioCorrelators as rat
dataRatFO2 = rat.TestRatFO2(DefWipe=False)
import TwoPtCorrelators as tpc
data2pt = tpc.TestNNQCorr(DefWipe=False)
dataRatFO2.CombineWithAlpha(data2pt,op.truediv)

print(dataRatFO2.Rat_Stats['alpha_rat_div'].loc['g4','q01-1','t_f6.01',:].apply(lambda x : x.Avg))
print(dataRatFO2.Rat_Stats['alpha_rat_div'].loc['g4','q01-1','t_f6.01',:].apply(lambda x : x.Std))
# data2pt.RemoveFuns()
# data2pt.GetFuns()
# for ikey,ival in data2pt.__dict__.iteritems():
#     print ikey,type(ival)
#     pickle.dump(ival,open(this_dir+'/TestGraphs/testpickle.py3p','wb'))
#
#
# for ikey,ival in data2pt.NNCPEven.__dict__.iteritems():
#     print ikey,type(ival)
#     pickle.dump(ival,open(this_dir+'/TestGraphs/testpickle.py3p','wb'))
#
# data2pt.NNCPEven
# for ikey,ival in data2pt.NNCPEven.C2_Stats['Auto'].iloc[0].__dict__.iteritems():
#     print ikey,type(ival)
#     pickle.dump(ival,open(this_dir+'/TestGraphs/testpickle.tstr(self.C3FOclass.NJN.tsink)p','wb'))



import FlowOpps as fo
from BootStrapping import FlattenBootstrapDF


data,data_comb = fo.TestFO(DefWipe=False)
data.Op_Stats
this_data = data.Op_Stats[['boot','Avg','Std']]
flat_data  = FlattenBootstrapDF(this_data)
dataFull = fo.TestFOFull(DefWipe=False)
# dataFull.RemoveFuns()
# dataFull.GetFuns()


# import numpy as np
# from XmlFormatting import untflowstr
# from Params import TflowToPhys
# tflowphys = np.array(map(untflowstr,data.tflowlist))
# tflowphys = TflowToPhys(np.array(tflowphys),0.0908)
# tflowphys = map(str,tflowphys)
# formatted = [itflow + ' '+str(itstr) for itflow,itstr in zip(tflowphys,data.tflowlist)]
# print '\n'.join(formatted)
# data.FlowGetCfgList()
# data.FlowRead()
# data.FlowWrite)(

import BootStrapping as bs
bdata,bdata2,bdata3 = bs.TestBoot(n_block=3)
bdata.Select_Block(3)
bdata.block_sel
bdata.MakeValAndErr()
comb_data = bdata*bdata2
comb_data.blocked_cfgvals

import Autocorr as ac
import numpy as np
from XmlFormatting import MakeValAndErr
adatab,adatab2,adatab3,bdatab,bdatab2,bdatab3,nbb,nbb2,nbb3 = ac.TestBlocking(tau_cutoff=0.5,block_cutoff=10)
print(adatab2.n_block)
print()
print(MakeValAndErr(adatab2.auto_object.tau,adatab2.auto_object.tauerr))
print()
print(bdatab2.MakeValAndErr())
print()
print(MakeValAndErr(adatab2.auto_object.Avg,adatab2.auto_object.Std))

print()
print()
print(MakeValAndErr(adatab2.noblock_auto.tau,adatab2.noblock_auto.tauerr))
print()
print(nbb2.MakeValAndErr())
print()
print(MakeValAndErr(adatab2.noblock_auto.Avg,adatab2.noblock_auto.Std))

bdatab.MakeValAndErr()
abs(bdatab+nbb2).MakeValAndErr()
bdatab2.name
MakeValAndErr(adatab2.auto_object.tau,adatab2.auto_object.tauerr)
print(diff2.MakeValAndErr())
print(bdatab2.MakeValAndErr())
print(nbb2.MakeValAndErr())
adata,adata2,adata3,adatarat = ac.TestAuto()

bdatarat = 100*bdata*bdata2
bdatarat.Stats()

print()
print('combination of data')
print('boot Avg: ',bdatarat.Avg)
print('auto Avg: ',adatarat.Avg)
print('boot Std: ',bdatarat.Std)
print('auto Std_uncorr: ',adatarat.Std_corr)
print('err diff uncorr: ',abs(bdatarat.Std-adatarat.Std_corr))
print('auto Std: ',adatarat.Std)
print('err diff: ',abs(bdatarat.Std-adatarat.Std))
print('auto Wopt: ',adatarat.Wopt)
print('auto taurat: ',adatarat.tauerr/adatarat.tau)


print()
print('regular data')
print('boot Avg: ',bdata.Avg)
print('auto Avg: ',adata.Avg)
print('boot Std: ',bdata.Std)
print('auto Std_uncorr: ',adata.Std_corr)
print('err diff uncorr: ',abs(bdata.Std-adata.Std_corr))
print('auto Std: ',adata.Std)
print('err diff: ',abs(bdata.Std-adata.Std))
print('auto Wopt: ',adata.Wopt)
print('auto tau: ',adata.tauerr/adata.tau)


median = data.bootvals.median()
low = data.bootvals.quantile(0.32)
high = data.bootvals.quantile(1-0.32)
medstd = (high - low)/2
print(median)
print(low)
print(high)
print(medstd)
print(median +medstd)

import pandas as pa
datacomb = pa.Series([data]).to_frame('pie').loc[:,'pie']/float(23)
datacomb
datacomb.apply(lambda x : x.Stats())
print(data)
print(datacomb[0])
import numpy as np


datanp = np.array(data)
datanp

import MomParams as mp
import numpy as np
from Params import defInfo
defInfo['a'] = 0.0907
defInfo['nxyzt'] = [32,32,32,64]
momdata = mp.LatticeParameters()
momdata.LoadPickle()
print(momdata.latspace*np.sqrt(8*4.01))
print(momdata.latspace*np.sqrt(8*5.01))
print(momdata.latspace*np.sqrt(8*6.01))

import TimeStuff as ts

ntime = 100000
test_list = [list(range(1000)),list(range(100))]
thistimer = ts.Timer(test_list,name='Test Timer:')
import itertools
for it in itertools.product(*test_list):
        thistimer.Lap(it)

# momdata.GetMomList(list_type='form')
# del data.Op_cfgs
# data.FlowGetCfgList()
# data.FlowLoadPickle()
# test_data = pickle.load(open(data.flowPickleFile,'rb'))
# test_data['Op_cfgs']
# data.Op_cfgs
# data.Write_Cfgs()__name__

# data.Op_cfgs
# for ikey,ival in data.flowFit_Stats['boot'].iloc[0].__dict__.iteritems():
#     print ikey,type(ival)
#     pickle.dump(ival,open(this_dir+'/TestGraphs/testpickle.py3p','wb'))
#
# %timeit pickle.load(open(this_dir+'/TestGraphs/testpickle.py3p','rb'))
#
# %timeit data.FlowWrite()
# %timeit data.FlowRead(show_timer=False)
#
# for ikey,val in data.__dict__.iteritems():
#     print
#     print ikey
#     print val
# data.Op_cfgs.memory_usage()
# data.Get_Flattened_Cfgs()[0].memory_usage()
#
# %timeit data.FlowLoadPickle()
# %timeit data.FlowLoadPickle()
#
#
# %timeit data.Write_Cfgs(show_timer=False)
# %timeit data.Read_Cfgs(show_timer=False)
#
#
# %timeit pickle.dump(data.Get_Flattened_Cfgs()[0],open(this_dir+'/TestGraphs/testpickle.py3p','wb'))
# %timeit data.Write_Cfgs()
# %timeit -n 20 data.Read_Cfgs(file_type='read_pickle',show_timer=False)
#
#
# %timeit data.Write_Cfgs(file_type='to_feather',show_timer=False)
# %timeit data.Read_Cfgs(file_type='read_feather',show_timer=False)
#
#
# %timeit data.Write_Cfgs(file_type='to_parquet',show_timer=False)
# %timeit -n 20 data.Read_Cfgs(file_type='read_parquet',show_timer=False)
#
#
#
# %timeit data.Write_Cfgs(file_type='to_msgpack',show_timer=False)
# %timeit -n 20 data.Read_Cfgs(file_type='read_msgpack',show_timer=False)
#
#
# %timeit data.Write_Cfgs(file_type='to_stata',show_timer=False)
# %timeit -n 20 data.Read_Cfgs(file_type='read_stata',show_timer=False)
#
#
#
#
# prev_cfglist = data.Op_cfgs
# data.Op_cfgs.head()
# prev_cfglist.head()
#
# all(data.Op_cfgs['stream_cfgs'] == data.Op_cfgs['stream_cfgs'])
#
# cfg_list,dump = data.Get_Flattened_Cfgs()
#
# cfg_list['flow_time'].dtype
# cfg_list.head()
# cfg_list['stream'] = map(np.where,cfg_list['stream'].values)
#
# data.stream_list
#
# cfg_list[cfg_list['flow_time'] == 9]
#
# new_df = pa.DataFrame(columns = data.Op_cfgs.columns[data.Op_cfgs.columns != 'Op'].tolist() + ['flow_time'],index = data.Op_cfgs.index)
# new_df
#
# this_dict = OrderedDict()
# this_dict['b'] = 4,5,6
# this_dict['a'] = 1,2,3
# pa.DataFrame(this_dict,index=['one','two','three'])


# from MiscFuns import CombineListCfgs
from FileIOMat import ReadWithMeta_NP

file2pt = '/home/jackdra/PHD/CHROMA/TestVar/Scripts/Python_Analysis/Configs/RC32x64Kud01375400Ks01364000/G2/TopCharge/Baryon_-ab-_CPEven_tsrc0_ism64_jsm64.cfgs.msg'

fileFO = '/home/jackdra/LQCD/Data/Qconfs/RC32x64_kud1375400_ks1364000_-ab-_TopCharge_Full.cfgs.msg'

meta2,cfglist_2pt = ReadWithMeta(file2pt)
metafo,cfglist_FO = ReadWithMeta(fileFO)
metafo,cfglist_FO = ReadWithMeta_NP(fileFO)
cfglist_FO.shape
np.asarray(metafo).shape
np.array(cfglist_FO['Op'].apply(lambda x : np.array(x)).values)


np.asarray(cfglist_FO['Op'].head(5).apply(lambda x : list(x)).values)
outlist = []
for ival in cfglist_FO['Op'].values:
    outlist.append(list(ival))
np.array(outlist).shape

import numpy as np
for ival in cfglist_FO['Op'].values:
    print(np.array(ival).shape)

#
# class twopt(object):
#     def __init__(self,cfglist):
#         self.test = cfglist
#
# class2pt = twopt(cfglist_2pt)
#
# classfo = twopt(cfglist_FO)
# comb_cfglist = CombineListCfgs([class2pt,classfo],'test')
# for icfg in cfglist_2pt['stream_cfgs'].values:
#     if icfg not in cfglist_FO['stream_cfgs'].values:
#         print icfg
# cfglist_2pt = cfglist_FO.drop(('-a-',6))
# cfglist_2pt['stream_cfgs']
# cfglist_FO['stream_cfgs']
#
# comb_cfglist['stream_cfgs']
