#!/usr/bin/env python

# data = [0.000450407533726283, -9.56943763353717e-22, 0.00125881827032569, 1.45686539151048e-21, 0.00144307889177642, 7.33567031666753e-22, 0.00110937528111652, -3.34872438764729e-22, 0.00087627896599212, -1.05972499343127e-21, 0.000684362650375569, 2.96072498157164e-22, 0.000476396547438195, -2.38728848426779e-22, 0.000311467405236602, 7.3551219857596e-23, 0.000203458813921517, -5.15695413138528e-24, 0.000128290916647265, 1.58230603893444e-23, 7.89872580503493e-05, -5.75633695805945e-24, 5.24843401587603e-05, 1.02084912022841e-23, 3.6917769494426e-05, 0.000631482945704443, -5.74599720040254e-23, 0.000929791000522095, 1.92607066928412e-22, 0.00117727057735764, 2.91128801527452e-23, 0.000940934205012497, 2.73349265118964e-22, 0.000328144847069594, 3.51550144747903e-22]
#
# data = [idata for idata in data if idata > 10e-10]
# import matplotlib.pyplot as pl
# pl.plot(data)
# pl.savefig('/home/jackdra/LQCD/Scripts/Python_Analysis/TestGraphs/flow_plot_test.pdf')
# pl.clf()

import pandas as pa
import numpy as np
from NullPlotData import null_series
from PlotData import Plotting
from MomParams import hbarc
from BootStrapping import BootStrap
from SetsOfFits import SetOfFitFuns
# from PredefFitFuns import c3FitFun,C2OneStateFitFunNoExp,C2TwoStateFitFunNoExp
from PredefFitFuns import c3ContFitFun,c3FitFun
from PredefFitFuns import c3FitFun_NoA
from PredefFitFuns import c3PolyFun,c3PolyFun2,c3PolyFun3
from PredefFitFuns import c3PolyFun_skip1,c3PolyFun2_skip1,c3PolyFun3_skip1
from PredefFitFuns import c3ExpFitFun,c3ExpFitFun2,c3ExpFitFun_nosqrt
from PredefFitFuns import c3ExpFitFun2_nosqrt,c3ExpFitFun3_nosqrt
from PredefFitFuns import C2OSFAntiper,C2OSFAntiperDer
import FitFunctions as ff
from ReadBinaryCfuns import RC2Full

import os
import pickle as pik

this_folder = '/home/jackdra/LQCD/Scripts/Jose_PoS/DataGraphs/FitRes/'
temp_folder = '/home/jackdra/temp/'

Fit_ens_files = {}
Fit_ens_files['M'] = ['DataFrame_C_32X64_1.csv','DataFrame_C_32X64_2.csv','DataFrame_C_32X64_3.csv']
Fit_ens_files['A'] = ['DataFrame_C_16X32.csv','DataFrame_C_20X40.csv','DataFrame_C_28X56.csv']

do_exraps = True
tree_antit = True
do_nucleon = False
double_tree_fix = True
# Mass_ens_files = {}
# Mass_ens_files['M'] = ['Ave_E_32X64_1.csv','Ave_E_32X64_2.csv','Ave_E_32X64_3.csv']
# Mass_ens_files['A'] = ['Ave_E_16X32.csv','Ave_E_20X40.csv','Ave_E_28X56.csv']
#
# AA0_ens_files = {}
# AA0_ens_files['M'] = ['AA0_ave_v2_32X64_1','AA0_ave_v2_32X64_2','AA0_ave_v2_32X64_3']
# AA0_ens_files['A'] = ['AA0_ave_v2_16X32','AA0_ave_v2_20X40','AA0_ave_v2_28X56']
#
# AA0_err_ens_files = {}
# AA0_err_ens_files['M'] = ['AA0_err_v2_32X64_1','AA0_err_v2_32X64_2','AA0_err_v2_32X64_3']
# AA0_err_ens_files['A'] = ['AA0_err_v2_16X32','AA0_err_v2_20X40','AA0_err_v2_28X56']



M_labels = ['mpi410','mpi570','mpi700']
A_labels = ['L16','L20','L28']
# M_labels = ['mpi410']
# A_labels = []
tsink_ens_dict = {}
tsink_ens_dict['mpi410'] =  64
tsink_ens_dict['mpi570'] =  64
tsink_ens_dict['mpi700'] =  64
tsink_ens_dict['L16'] =  32
tsink_ens_dict['L20'] =  40
tsink_ens_dict['L28'] =  56

a_dict = {}
a_dict['mpi410'] =  0.0907
a_dict['mpi570'] =  0.0907
a_dict['mpi700'] =  0.0907
a_dict['L16'] =  0.1095
a_dict['L20'] =  0.0936
a_dict['L28'] =  0.0684
t0_dict = {}
t0_dict['mpi410'] =     2.5377162290
t0_dict['mpi570'] =     2.3999592122
t0_dict['mpi700'] =     2.2590866059
t0_dict['L16'] =        1.3628954894
t0_dict['L20'] =        2.2386627909
t0_dict['L28'] =        4.9885618908

t0_dict_err = {}
t0_dict_err['mpi410'] =  0.0015529386
t0_dict_err['mpi570'] =  0.0011452810
t0_dict_err['mpi700'] =  0.0011742296
t0_dict_err['L16'] =     0.0015817442
t0_dict_err['L20'] =     0.0021913620
t0_dict_err['L28'] =     0.0064980950

# from XmlFormatting import MakeValAndErr
# for ival_err,ival in zip(t0_dict_err.values(),t0_dict.values()):
#     print(MakeValAndErr(ival,ival_err,Dec=2))
use_t0rat = True
if use_t0rat:
    sqrt8t_key = 'sqrt8t_t0rat'
else:
    sqrt8t_key = 'sqrt8t_fm'

def MakeFileList(this_dict):
    output_dict = {}
    for ifile,ilab in zip(this_dict['M'],M_labels):
        output_dict[ilab] = this_folder + ifile
    for ifile,ilab in zip(this_dict['A'],A_labels):
        output_dict[ilab] = this_folder + ifile
    return output_dict

Fit_filelist = MakeFileList(Fit_ens_files)

flow_list = [f'{val:.1f}' for val in np.arange(0,1,0.1)]
flow_list += [f'{val:.1f}' for val in np.arange(1,10,1)]
def ReadTLFileList(this_dict):
    output_dict = {}
    for ikey,ifile in this_dict.items():
        tf_data = []
        t_list,tf_list = [],[]
        for iflow in flow_list:
            this_file = ifile.replace('FLOW_TIME',iflow)
            if do_nucleon:
                this_data = RC2Full([this_file],['0 0 0'],
                                    InterpNumb='9',
                                    MesOrBar='Baryon').data
            else:
                this_data = RC2Full([this_file],['0 0 0'],
                                    InterpNumb='15',
                                    MesOrBar='Meson').data
            this_data = list(np.array(this_data)[0,:,0].flatten())
            tf_data += this_data
            t_list += [f'{it}' for it in range(len(this_data))]
            tf_list += [iflow for it in this_data]
        output_dict[ikey] = pa.DataFrame()
        output_dict[ikey]['t_sink'] = pa.Series(t_list)
        output_dict[ikey]['t_flow'] = pa.Series(tf_list)
        output_dict[ikey]['Corr'] = pa.Series(tf_data)
        if double_tree_fix:
            output_dict[ikey]['Corr'] = output_dict[ikey]['Corr']*2
    return output_dict

def ReadFileList(this_dict):
    output_dict = {}
    for ikey,ifile in this_dict.items():
        with open(ifile,'r') as f:
            output_dict[ikey] = pa.read_csv(f)
    return output_dict

Fit_data = ReadFileList(Fit_filelist)
# this_kappa = 'k0.1250'
# this_kappa = 'k0.1260'
# this_kappa = 'k0.1270'
# tree_kappa_list = ['k0.1250','k0.1260','k0.1270']
# tree_kappa_list = ['k0.1250']
import numpy as np
tree_kappa_list = ['k0.1250','k0.1237','k0.1225','k0.1212','k0.1200']
tree_m0_list = [float(ival.replace('k','')) for ival in tree_kappa_list]
tree_m0_list = [1/ival - 8 for ival in tree_m0_list]
tree_mp_list = [np.log(1+ival) for ival in tree_m0_list]
tree_mpi_list = [2*ival for ival in tree_mp_list]
tree_df_info = pa.DataFrame()
tree_df_info['kappa'] = pa.Series(tree_kappa_list)
tree_df_info['m0'] = pa.Series(tree_m0_list)
tree_df_info['mp'] = pa.Series(tree_mp_list)
tree_df_info['mpi'] = pa.Series(tree_mpi_list)
tree_df_info = tree_df_info.set_index('kappa')

if tree_antit:
    tree_base = '/home/jackdra/LQCD/Scripts/Jose_PoS/DataGraphs/Treelev/'
else:
    tree_base = '/home/jackdra/LQCD/Scripts/Jose_PoS/DataGraphs/Treelev_Tper/'

tree_files = {}
for ikappa in tree_kappa_list:
    for ikey in Fit_filelist.keys():
        if 'mpi' in ikey:
            tree_files[ikey+'_'+ikappa] = tree_base+'/RC32x64_tree/RC32x64_tree_'+ikappa+'_FLOW_TIME.xml'
        else:
            this_RC = int(ikey[1:])
            this_RC = 'RC'+str(this_RC)+'x'+str(this_RC*2)
            tree_files[ikey+'_'+ikappa] =        tree_base+this_RC+'_tree/'+this_RC+'_tree_'+ikappa+'_FLOW_TIME.xml'


tree_data = ReadTLFileList(tree_files)
tree_data['mpi410_k0.1200']
effm_tree_data = {}
from pynverse import inversefunc
for iens,ens_data in tree_data.items():
    this_tsink = tsink_ens_dict[iens.split('_')[0]]
    val_list,AA0_list = [],[]
    flow_0_data = ens_data[ens_data['t_flow'] == '0.0']['Corr'].values
    for itflow in flow_list:
        flow_data = ens_data[ens_data['t_flow'] == itflow]
        flow_data = flow_data['Corr'].values
        for it,(tdata,tdatap1) in enumerate(zip(flow_data[:-1],flow_data[1:])):
            AA0_list.append(tdata/flow_0_data[it])
            if do_nucleon:
                val_list.append(np.log(np.abs(tdata/tdatap1)))
            else:
                Tshift = this_tsink/2-it
                def ThisEffMFun(mass):
                    return np.cosh(-mass*Tshift)/np.cosh(-mass*(Tshift-1))
                val_list.append(np.abs(inversefunc(ThisEffMFun,y_values=tdata/tdatap1)))
        val_list.append(float('NaN'))
        AA0_list.append(float('NaN'))
    if len(val_list) > 0:
        tree_data[iens]['EffM'] = pa.Series(val_list,index=tree_data[iens].index)
        tree_data[iens]['AA0'] = pa.Series(AA0_list,index=tree_data[iens].index)

this_data = pa.DataFrame()
for iens,ens_data in tree_data.items():
    plot_data = ens_data.set_index(['t_sink','t_flow'])['EffM']
    hold_series = null_series
    hold_series['type'] = 'plot_vary'
    hold_series['key_select'] = (slice(None),'0.0')
    hold_series['x_data'] = 'from_keys'
    hold_series['y_data'] = plot_data
    hold_series['label'] = iens
    this_data[iens] = hold_series
for ikappa,mpi_val in tree_df_info['mpi'].items():
    hold_series = null_series
    hold_series['type'] = 'hline'
    hold_series['y_data'] = [mpi_val]
    hold_series['label'] = 'from kappa ' + ikappa
    this_data['from kappa ' + ikappa] = hold_series

this_info = pa.Series()
this_info['save_file'] = '/home/jackdra/LQCD/Scripts/Jose_PoS/DataGraphs/Jacks/Test_EffM.pdf'
if tree_antit:
    this_info['save_file'] = this_info['save_file'].replace('Test_EffM','Test_EffM_antit')
if do_nucleon:
    this_info['save_file'] = this_info['save_file'].replace('Test_EffM','Test_EffM_nucleon')
this_info['title'] = r'$M_{eff}$'
this_info['ylabel'] = r'$M_{eff}$'
this_info['xlabel'] = r'$t_{sink}$'
# this_info['xlims'] = [0,10]
# this_info['ylims'] = [0,15]
data_plot = Plotting(plot_data=this_data,plot_info=this_info)
data_plot.PlotAll()

this_data = pa.DataFrame()
for iens,ens_data in tree_data.items():
    # plot_data = ens_data.set_index(['t_sink','t_flow'])['AA0']
    plot_data = ens_data.set_index(['t_sink','t_flow'])['Corr']
    hold_series = null_series
    hold_series['type'] = 'plot_vary'
    hold_series['key_select'] = (slice(None),'0.0')
    hold_series['x_data'] = 'from_keys'
    hold_series['y_data'] = plot_data
    hold_series['label'] = iens
    this_data[iens] = hold_series
this_info = pa.Series()
this_info['save_file'] = '/home/jackdra/LQCD/Scripts/Jose_PoS/DataGraphs/Jacks/Test_Corr.pdf'
if tree_antit:
    this_info['save_file'] = this_info['save_file'].replace('Test_Corr','Test_Corr_antit')
if do_nucleon:
    this_info['save_file'] = this_info['save_file'].replace('Test_Corr','Test_Corr_nucleon')
this_info['title'] = r'$Correlator$'
this_info['ylabel'] = r'$C$'
this_info['xlabel'] = r'$t_{sink}$'
# this_info['xlims'] = [0,10]
# this_info['ylims'] = [0,15]
data_plot = Plotting(plot_data=this_data,plot_info=this_info)
data_plot.PlotAll()
# tree_data['L16_k0.1250'].set_index(['t_sink','t_flow'])['Corr'][:,'5.0']
# import matplotlib.pyplot as pl
# flow_test = tree_data['mpi410_k0.1270']
# flow_test = flow_test[flow_test['t_flow'] == '0.0']['Corr']
# flow_test = flow_test
# flow_test
# flow_test.plot(color='blue')
# # pl.savefig('/home/jackdra/LQCD/Scripts/Python_Analysis/TestGraphs/flow_plot_test2.pdf')
# # pl.clf()
# flow_test = tree_data['mpi410_k0.1250']
# flow_test = flow_test[flow_test['t_flow'] == '4.0']
# flow_test['Corr'].plot(color='red')
# pl.savefig('/home/jackdra/LQCD/Scripts/Python_Analysis/TestGraphs/flow_plot_test.pdf')
# pl.clf()

def DoLSFit(x_data,y_data,this_name,this_fun):
    testdf = pa.DataFrame()
    testdf.loc[:,'xdata'] = pa.Series([[x_data]])
    testdf.loc[:,'ydata'] = pa.Series([y_data])
    testfit = ff.Fitting(data = testdf, Funs=this_fun,name=this_name)
    return testfit.FitMean()

tree_fit_data = {}
DefWipe = False
for iens,ens_data in tree_data.items():
    tree_fit_file = temp_folder+'/'+iens+'_tree_fit_data.p'
    if os.path.isfile(tree_fit_file) and not DefWipe:
        tree_fit_data[iens] = pa.read_pickle(tree_fit_file)
    else:
        this_tsink = tsink_ens_dict[iens.split('_')[0]]
        def C2atT(t,p):
            return C2OSFAntiper(t,p,this_tsink+1)
        C2atT.file_name = 'C2atT'+str(this_tsink+1)
        def C2DeratT(t,p):
            return C2OSFAntiperDer(t,p,this_tsink+1)
        C2DeratT.file_name = 'C2DeratT'+str(this_tsink+1)
        this_fun = [C2atT,2,C2DeratT]
        all_flow_list,tcut_list = [],[]
        E_list,A_list,chi_list = [],[],[]
        print('\r','Begining fit for '+iens,' '*20)
        for itflow in flow_list:
            print('\r','     tflow:', itflow, ' '*20)
            flow_data = ens_data[ens_data['t_flow'] == itflow]
            for itcut in range(int(this_tsink/2)):
                print('\r','         itcut:', itcut, end=' '*20)
                fit_data = flow_data[flow_data['t_sink'].apply(lambda x : int(x)) > itcut]
                fit_data = fit_data[fit_data['t_sink'].apply(lambda x : int(x)) < this_tsink-itcut+1]
                xdata,ydata = fit_data['t_sink'].apply(lambda x : float(x)),fit_data['Corr']
                test_fit_res = DoLSFit(xdata.values,ydata.values,'test_fit',this_fun)
                all_flow_list.append(itflow)
                tcut_list.append(itcut)
                E_list.append(test_fit_res['Params']['Energy'])
                A_list.append(test_fit_res['Params']['A_Ep'])
                chi_list.append(test_fit_res['Chi2DoF']['Energy'])
        print('fitting tree complete ')
        tree_fit_data[iens] = pa.DataFrame()
        if len(tcut_list) > 0:
            tree_fit_data[iens]['flow t'] = pa.Series(all_flow_list)
            tree_fit_data[iens]['tmin'] = pa.Series(tcut_list)
            tree_fit_data[iens]['A'] = pa.Series(E_list)
            tree_fit_data[iens]['E'] = pa.Series(A_list)
            tree_fit_data[iens]['chisq'] = pa.Series(chi_list)
        tree_fit_data[iens].to_pickle(tree_fit_file)
# tree_fit_data['mpi410_k0.1250']
# # print(tree_fit_data['L16_k0.1250'].to_string())
# tree_fit_data['L16_k0.1250']
# tree_fit_data['L16_k0.1250'][tree_fit_data['L16_k0.1250']['flow t']=='4.0']
# test_fit = test_fit[test_fit['t_sink'].apply(lambda x : int(x)) > 15]
# test_fit = test_fit[test_fit['t_sink'].apply(lambda x : int(x)) < 64-15]
# xdata,ydata = test_fit['t_sink'].apply(lambda x : float(x)),test_fit['Corr']
# this_fun = [C2atT,2,C2DeratT]
# test_fit_res = DoLSFit(xdata.values,ydata.values,'test_fit',this_fun)
# test_fit_res
# fit_eval = np.array([C2atT([ival],[test_fit_res['Params']['A_Ep'],test_fit_res['Params']['Energy']]) for ival in range(15,64-15-1)])
# ydata.values
# np.abs(fit_eval-ydata.values)

def BootStrap_Vals(this_DF):
    A_list,E_list,chi_list,ilist = [],[],[],[]
    for ikey,iDF in this_DF.groupby(('flow t','tmin')):
        A_list.append(BootStrap(thisnboot=len(iDF),bootvals=iDF['A'].values))
        E_list.append(BootStrap(thisnboot=len(iDF),bootvals=iDF['E'].values))
        chi_list.append(BootStrap(thisnboot=len(iDF),bootvals=iDF['chisq'].values))
        ilist.append(ikey)
    if len(ilist) > 0:
        ilist = pa.MultiIndex.from_tuples(ilist,names=('flow t','tmin'))
        this_DF = pa.DataFrame()
        this_DF['A'] = pa.Series(A_list,index = ilist)
        this_DF['E'] = pa.Series(E_list,index = ilist)
        this_DF['chisq'] = pa.Series(chi_list,index = ilist)
    return this_DF.reset_index()

def MakeAdiv(this_DF,flow_str = False):
    out_list,ilist = [],[]
    for it,ival in this_DF.groupby('tmin'):
        if flow_str:
            DF_tf0 = ival['A'][ival['flow t'] == '0.0'].values[0]
        else:
            DF_tf0 = ival['A'][ival['flow t'] == 0].values[0]
        for jkey,jval in ival['A'].items():
            if hasattr(jval,'__div__'):
                outval = jval.__div__(DF_tf0)
                outval.Stats()
            else:
                outval = jval/DF_tf0
            out_list.append(outval)
            ilist.append(jkey)
    if len(ilist) > 0:
        return pa.Series(out_list,index=ilist)
    else:
        return pa.Series()

first_key,second_key = {},{}
# DefWipe = DefWipe or False
DefWipe = True
for ikey,ival in Fit_data.items():
    data_file = temp_folder+'/'+ikey+'_wtree_data.p'
    if os.path.isfile(data_file) and not DefWipe:
        with open(data_file,'rb') as f:
            Fit_data[ikey],first_key[ikey],second_key[ikey] = pik.load(f)
    else:
        Fit_data[ikey] = BootStrap_Vals(ival)
        Fit_data[ikey]['A(t)/A(0)'] = MakeAdiv(Fit_data[ikey])
        Fit_data[ikey]['sqrt8t_fm'] = Fit_data[ikey]['flow t'].apply(lambda x : str(np.sqrt(8*x)*a_dict[ikey]))
        print('REMOVE THIS, volume scaling')
        # Fit_data[ikey]['sqrt8t_t0rat'] = Fit_data[ikey]['flow t'].apply(lambda x : str(np.sqrt(x/t0_dict[ikey])))
        Fit_data[ikey]['sqrt8t_t0rat'] = Fit_data[ikey]['flow t'].apply(lambda x : str(2*np.sqrt(8*x)/tsink_ens_dict[ikey]))
        Fit_data[ikey]['tmin_fm'] = Fit_data[ikey]['tmin'].apply(lambda x : str(x*a_dict[ikey]))
        Fit_data[ikey]['E_GeV'] = Fit_data[ikey]['E'].apply(lambda x : x*hbarc/a_dict[ikey])
        Fit_data[ikey]['E_GeV'].apply(lambda x : x.Stats())
        for ikappa in tree_kappa_list:
            Fit_data[ikey]['A(t)/A(0)_tree_'+ikappa] = MakeAdiv(tree_fit_data[ikey+'_'+ikappa],flow_str = True)
            # Fit_data[ikey]['A(t)/A(0)_div_tree_'+ikappa] = Fit_data[ikey]['A(t)/A(0)']/Fit_data[ikey]['A(t)/A(0)_tree_'+ikappa]
            print('THIS FORCES IT TO tree itself, change comment to divide by tree level')
            Fit_data[ikey]['A(t)/A(0)_div_tree_'+ikappa] = Fit_data[ikey]['A(t)/A(0)']*Fit_data[ikey]['A(t)/A(0)_tree_'+ikappa]/Fit_data[ikey]['A(t)/A(0)']
            Fit_data[ikey]['A(t)/A(0)_div_tree_'+ikappa].apply(lambda x : x.Stats())
            Fit_data[ikey]['E_tree_'+ikappa] = tree_fit_data[ikey+'_'+ikappa]['E']
        first_key[ikey] = Fit_data[ikey][sqrt8t_key].iloc[0]
        second_key[ikey] = Fit_data[ikey]['tmin_fm'].iloc[0]
        Fit_data[ikey] = Fit_data[ikey].set_index([sqrt8t_key,'tmin_fm'])
        with open(data_file,'wb') as f:
            pik.dump((Fit_data[ikey],first_key[ikey],second_key[ikey]),f)
# Fit_data['L20']['A(t)/A(0)_tree_k0.1250'][:,'0.1872'].plot()
# print(Fit_data['mpi410']['A(t)/A(0)'][:,'2.1768'].apply(lambda x : x.MakeValAndErr()).to_string())
# Fit_data['mpi410']['A(t)/A(0)_div_tree_'+ikappa].apply(lambda x : x.Stats())
# Fit_data['mpi410']['A(t)/A(0)_div_tree_'+ikappa].apply(lambda x : x.MakeValAndErr())
fit_info = {}
# fit_info['Funs'] = [c3FitFun_NoA,2]
fit_info['Funs'] = [c3PolyFun,3]
# fit_info['Funs'] = [c3PolyFun,2]
# fit_info['Funs'] = [c3FitFun,3]
# fit_info['Funs'] = [c3FitFun_V2,2]
# fit_info['Funs'] = [C2OneStateFitFunNoExp,2]

fit_info_list = []
fit_info = {}
fit_info['Funs'] = [c3FitFun,3]
fit_info_list.append(fit_info)
# fit_info = {}
# fit_info['Funs'] = [c3PolyFun,1]
# fit_info_list.append(fit_info)
# fit_info = {}
# fit_info['Funs'] = [c3PolyFun2,2]
# fit_info_list.append(fit_info)
# fit_info = {}
# fit_info['Funs'] = [c3PolyFun3,3]
# fit_info_list.append(fit_info)
# fit_info = {}
# fit_info['Funs'] = [c3PolyFun_skip1,1]
# fit_info_list.append(fit_info)
# fit_info = {}
# fit_info['Funs'] = [c3PolyFun2_skip1,2]
# fit_info_list.append(fit_info)
# fit_info = {}
# fit_info['Funs'] = [c3PolyFun3_skip1,3]
# fit_info_list.append(fit_info)
# fit_info['Funs'] = [c3PolyFun,2]
# fit_info = {}
# fit_info['Funs'] = [c3FitFun,3]
# fit_info_list.append(fit_info)
# fit_info = {}
# fit_info['Funs'] = [c3ExpFitFun,1]
# fit_info_list.append(fit_info)
# fit_info = {}
# fit_info['Funs'] = [c3ExpFitFun2,2]
# fit_info_list.append(fit_info)
# fit_info = {}
# fit_info['Funs'] = [c3ExpFitFun_nosqrt,1]
# fit_info_list.append(fit_info)
# fit_info = {}
# fit_info['Funs'] = [c3ExpFitFun2_nosqrt,2]
# fit_info_list.append(fit_info)
# fit_info = {}
# fit_info['Funs'] = [c3ExpFitFun3_nosqrt,3]
# fit_info_list.append(fit_info)
# fit_info['Funs'] = [C2OneStateFitFunNoExp,2]
# fit_info['iGuess'] = [1,1]
# fit_info = {}
# fit_info['Funs'] = [c3ExpFitFun,1]
# fit_info_list = [fit_info]

fit_tmin_min = 1.5
fit_tmin_max = 2.0
tflow_fit_min = 0.1
tflow_fit_max = 0.4
min_par_tflow = 0
min_par_tmin = 0
best_fit_picked = [(c3FitFun.__name__,'fitr4-11','fittwor17-22') for _ in range(3)]
best_fit_picked.append((c3FitFun.__name__,'fitr4-10','fittwor14-14'))
best_fit_picked.append((c3FitFun.__name__,'fitr4-11','fittwor17-18'))
best_fit_picked.append((c3FitFun.__name__,'fitr5-11','fittwor22-26'))
min_hold = tflow_fit_min
max_hold = tflow_fit_max
tflow_fit_min = {}
tflow_fit_max = {}
if use_t0rat:
    print(min_hold,' to ',max_hold, ' in fm corresponds to')
    for ikey in a_dict.keys():
        print(ikey)
        # tflow_fit_min[ikey] = min_hold/(a_dict[ikey]*np.sqrt(8*t0_dict[ikey]))
        # tflow_fit_max[ikey] = max_hold/(a_dict[ikey]*np.sqrt(8*t0_dict[ikey]))
        print('Volume Scaling fix this too!')
        tflow_fit_min[ikey] = 2*min_hold/(a_dict[ikey]*tsink_ens_dict[ikey])
        tflow_fit_max[ikey] = 2*max_hold/(a_dict[ikey]*tsink_ens_dict[ikey])
        print(tflow_fit_min[ikey],' to ',tflow_fit_max[ikey], ' in units of sqrt(8t_{0})')
        print()
else:
    for ikey in a_dict.keys():
        tflow_fit_min[ikey] = min_hold
        tflow_fit_max[ikey] = max_hold

cont_fit_info = {}
cont_fit_info['Funs'] = [c3ContFitFun,3]
tflow_picked = 0.15
tmin_picked = fit_tmin_min

Master_Fit_Data = pa.DataFrame()


A0Fit_data = pa.DataFrame()
fit_ilist,fit_data = [],[]

def GetMass_minfit(this_df,this_tflow,this_t):
    this_df = this_df.reset_index()
    this_df = this_df[this_df[sqrt8t_key].apply(lambda x : float(x)) > this_tflow]
    this_df = this_df[this_df['tmin_fm'].apply(lambda x : float(x)) > this_t]
    return this_df['E_GeV'].iloc[0]

DoTree = True
# picked_kappa = tree_kappa_list[4]
picked_kappa = tree_kappa_list[1]
# picked_kappa = tree_kappa_list[2]
if DoTree:
    # this_arat_key = 'A(t)/A(0)_tree_'+picked_kappa
    # file_kappa = 'AA0_tree_'+picked_kappa.replace('k0.','k')
    this_arat_key = 'A(t)/A(0)_div_tree_'+picked_kappa
    file_kappa = 'AA0_dtree_'+picked_kappa.replace('k0.','k')
else:
    this_arat_key = 'A(t)/A(0)'
    file_kappa = 'AA0'


mpi_dict = {}
DefWipe = DefWipe or False
for ikey,ival in Fit_data.items():
    fit_file = temp_folder+'/fit_'+file_kappa+'_'+ikey+'fitr'+'-'.join([str(int(fit_tmin_min*10)),
                                                         str(int(fit_tmin_max*10))])+'_'+fit_info['Funs'][0].__name__+'.p'
    if os.path.isfile(fit_file) and not DefWipe:
        with open(fit_file,'rb') as f:
            this_fit,this_df,mass_val = pik.load(f)
        this_fit.GetFuns()
        mpi_dict[ikey] = mass_val.Avg
        mpi_dict[ikey+'_err'] = mass_val.Std
    else:
        # ival = ival.reset_index()
        # ival = ival[ival['tmin_fm'].apply(lambda x : float(x)) >fit_tmin_min]
        # ival = ival[fit_tmin_max>ival['tmin_fm'].apply(lambda x : float(x))]
        # ival = ival.set_index([sqrt8t_key,'tmin_fm'])
        # ival = pa.Series(ival.values,index=ival.index)
        this_fit = SetOfFitFuns(data = ival[this_arat_key],name=ikey)
        # fit_info['Funs'] = [ConstantFitFun,1]
        # fit_info['iGuess'] = [0,0]
        # this_fit.ScanBox_FixEnd(tflow_fit_min,tflow_fit_max,
        #                         fit_tmin_min,fit_tmin_max,
        #                         fit_info,min_fit_len_1=min_par_tflow,min_fit_len_2=min_par_tmin,
        #                         from_xdata=True)
        this_fit.ScanBox_maxdim2(tflow_fit_min[ikey],tflow_fit_max[ikey],
                                fit_tmin_min,fit_tmin_max,
                                fit_info_list,min_fit_len_1=min_par_tflow,
                                min_fit_len_2=min_par_tmin,
                                from_xdata=True)
        this_fit.DoFits()
        this_fit.SortChi()
        this_df = ival.reset_index()[[sqrt8t_key,'tmin_fm',this_arat_key]]
        this_df['ens'] = ikey
        this_df['t0_vals'] = 1/t0_dict[ikey]
        mass_val = GetMass_minfit(Fit_data[ikey]['E_GeV'].reset_index(),tflow_fit_min[ikey],fit_tmin_min)
        mpi_dict[ikey] = mass_val.Avg
        mpi_dict[ikey+'_err'] = mass_val.Std
        # this_df = this_df.set_index([sqrt8t_key,'tmin_fm'])
        this_fit.RemoveFuns()
        with open(fit_file,'wb') as f:
            pik.dump((this_fit,this_df,mass_val),f)
        this_fit.GetFuns()
    fit_ilist.append((mpi_dict[ikey],1/t0_dict[ikey]))
    fit_data.append(this_fit)
    Master_Fit_Data = Master_Fit_Data.append(this_df)
if len(fit_data) > 0:
    ilist = pa.MultiIndex.from_tuples(fit_ilist,names=[r'$m_{\pi}$ GeV',r'$a^{2}/t_{0}$ fm'])
    A0Fit_data['Fit'] = pa.Series(fit_data,index=ilist)
    A0Fit_data['ens'] = pa.Series(list(a_dict.keys()),index=ilist)
    best_fit = [ifit.Fit_Stats_fmt.loc[ibest,'Fit'] for ibest,ifit in zip(best_fit_picked,fit_data)]
    eval_list = [ifit.Eval_Function(tflow_picked,tmin_picked) for ifit in best_fit]
    chi_list = [ifit.Get_Chi2DoF() for ifit in best_fit]
    A0Fit_data['Fit_Best_Eval'] = pa.Series(eval_list,index=ilist)
    A0Fit_data['Fit_Best_chi2pdf'] = pa.Series(chi_list,index=ilist)
    if 'a' in best_fit[0].fit_data['Params'].keys():
        par_list_a = [ifit.fit_data['Params']['a'] for ifit in best_fit]
        A0Fit_data['Par_a'] = pa.Series(par_list_a,index=ilist)
    if 'b' in best_fit[0].fit_data['Params'].keys():
        par_list_b = [ifit.fit_data['Params']['b'] for ifit in best_fit]
        A0Fit_data['Par_b'] = pa.Series(par_list_b,index=ilist)
    if 'c' in best_fit[0].fit_data['Params'].keys():
        par_list_c = [ifit.fit_data['Params']['c'] for ifit in best_fit]
        A0Fit_data['Par_c'] = pa.Series(par_list_c,index=ilist)
    if 'd' in best_fit[0].fit_data['Params'].keys():
        par_list_d = [ifit.fit_data['Params']['d'] for ifit in best_fit]
        A0Fit_data['Par_d'] = pa.Series(par_list_d,index=ilist)

Master_Fit_Data = Master_Fit_Data.set_index([sqrt8t_key,'tmin_fm']).sort_index()
# A0Fit_data['Fit'].iloc[5].Fit_Stats_fmt['Fit'].apply(lambda x : x.fit_data['Chi2DoF']).head()
# A0Fit_data['Fit'].iloc[0].Fit_Stats_fmt['Fit'].apply(lambda x : x.fit_data['Params'])
# if sqrt8t_key in Master_Fit_Data.columns:
#     del Master_Fit_Data[sqrt8t_key]
# if 'tmin_fm' in Master_Fit_Data.columns:
#     del Master_Fit_Data['tmin_fm']
fit_file = temp_folder+'/fit_combine'+file_kappa+'_fitr'+'-'.join([str(int(fit_tmin_min*10)),
                                                     str(int(fit_tmin_max*10))])+'_'+fit_info['Funs'][0].__name__+'.p'

# Master_Fit_Data
# Master_Fit_Data = Master_Fit_Data.reset_index()
# Master_Fit_Data['tmin_fm'].apply(lambda x : float(x)) >fit_tmin_min
# Master_Fit_Data = Master_Fit_Data[Master_Fit_Data['tmin_fm'].apply(lambda x : float(x)) >fit_tmin_min]
# Master_Fit_Data = Master_Fit_Data[fit_tmin_max>Master_Fit_Data['tmin_fm'].apply(lambda x : float(x))]
DefWipe = DefWipe or False
if do_exraps:
    par_fit_list = {}
    if os.path.isfile(fit_file) and not DefWipe:
        with open(fit_file,'rb') as f:
            if len(M_labels) > 1:
                combine_fit,cont_fit,par_fit_list = pik.load(f)
                cont_fit.GetFuns()
            else:
                combine_fit,par_fit_list = pik.load(f)
    else:
        combine_fit = SetOfFitFuns(data = Master_Fit_Data[this_arat_key],name='combined')
        combine_fit.ScanBox_FixEnd_maxdim2( tflow_fit_min['mpi410'],tflow_fit_max['mpi410'],
                                    fit_tmin_min,fit_tmin_max,
                                    fit_info_list,min_fit_len_1=min_par_tflow,min_fit_len_2=min_par_tmin,
                                    from_xdata=True)
        combine_fit.DoFits()
        combine_fit.SortChi()
        combine_fit.RemoveFuns()

        if len(M_labels) > 1:
            cont_data = A0Fit_data['Fit_Best_Eval']
            cont_fit = SetOfFitFuns(data = cont_data,name='cont_extrap')
            cont_fit.ScanBox( 0,1,0,10,cont_fit_info,min_fit_len_1=1
                                ,min_fit_len_2=1,from_xdata=True)
            cont_fit.DoFits()
            cont_fit.SortChi()
            cont_fit.RemoveFuns()

        if 'Par_a' in A0Fit_data.keys():
            Par_a_data = A0Fit_data['Par_a']
            Par_a_fit = SetOfFitFuns(data = Par_a_data,name='Par_a')
            Par_a_fit.ScanBox( 0,1,0,10,cont_fit_info,min_fit_len_1=1
                                ,min_fit_len_2=1,from_xdata=True)
            Par_a_fit.DoFits()
            Par_a_fit.SortChi()
            Par_a_fit.RemoveFuns()
            par_fit_list['Par_a'] = Par_a_fit

        if 'Par_b' in A0Fit_data.keys():
            Par_b_data = A0Fit_data['Par_b']
            Par_b_fit = SetOfFitFuns(data = Par_b_data,name='Par_b')
            Par_b_fit.ScanBox( 0,1,0,10,cont_fit_info,min_fit_len_1=1
                                ,min_fit_len_2=1,from_xdata=True)
            Par_b_fit.DoFits()
            Par_b_fit.SortChi()
            Par_b_fit.RemoveFuns()
            par_fit_list['Par_b'] = Par_b_fit

        if 'Par_c' in A0Fit_data.keys():
            Par_c_data = A0Fit_data['Par_c']
            Par_c_fit = SetOfFitFuns(data = Par_c_data,name='Par_c')
            Par_c_fit.ScanBox( 0,1,0,10,cont_fit_info,min_fit_len_1=1
                                ,min_fit_len_2=1,from_xdata=True)
            Par_c_fit.DoFits()
            Par_c_fit.SortChi()
            Par_c_fit.RemoveFuns()
            par_fit_list['Par_c'] = Par_c_fit

        if 'Par_d' in A0Fit_data.keys():
            Par_d_data = A0Fit_data['Par_d']
            Par_d_fit = SetOfFitFuns(data = Par_d_data,name='Par_d')
            Par_d_fit.ScanBox( 0,1,0,10,cont_fit_info,min_fit_len_1=1
                                ,min_fit_len_2=1,from_xdata=True)
            Par_d_fit.DoFits()
            Par_d_fit.SortChi()
            Par_d_fit.RemoveFuns()
            par_fit_list['Par_d'] = Par_d_fit


        with open(fit_file,'wb') as f:
            if len(M_labels) > 1:
                pik.dump((combine_fit,cont_fit,par_fit_list),f)
                cont_fit.GetFuns()
            else:
                pik.dump((combine_fit,par_fit_list),f)
    for ifit in par_fit_list.values():
        ifit.GetFuns()
    combine_fit.GetFuns()


# c = (1-0.66/0.43)/(np.log((0.22**2)/8)-((0.66/0.43)*np.log((0.11**2)/8)))
# b = (-0.43*(-0.11 + np.log((0.11**2)/8)))
# c
this_data = pa.DataFrame()
for ikey,ival in Fit_data.items():
    hold_series = null_series
    hold_series['type'] = 'error_bar_vary'
    hold_series['key_select'] = (slice(None),second_key[ikey])
    hold_series['x_data'] = 'from_keys'
    hold_series['y_data'] = ival[this_arat_key].apply(lambda x : x.Avg)
    hold_series['yerr_data'] = ival[this_arat_key].apply(lambda x : x.Std)
    hold_series['label'] = ikey
    this_data[ikey] = hold_series

    hold_series = null_series
    hold_series['type'] = 'fit_vary'
    # plot_fit = A0Fit_data[ikey].ReduceViaChi(nchi_cutoff=15)[0]
    plot_fit = A0Fit_data[A0Fit_data['ens'] == ikey].values[0][0]
    hold_series['key_select'] = 'First'
    hold_series['fit_class'] = plot_fit
    hold_series['label'] = ikey+' Fit'
    hold_series['xaxis'] = 0
    # hold_series['xdatarange'] = 'Data'
    hold_series['otherXvals'] = [1]
    this_data[ikey+'_Fit'] = hold_series

if do_exraps:
    hold_series = null_series
    hold_series['type'] = 'fit_vary'
    # plot_comb = combine_fit.ReduceViaChi(nchi_cutoff=20)[0]
    plot_comb = combine_fit
    # hold_series['key_select'] = plot_comb.Fit_Stats.index[0]
    hold_series['key_select'] = 'First'
    hold_series['fit_class'] = plot_comb
    hold_series['label'] = 'Combine Fit'
    hold_series['xaxis'] = 0
    hold_series['otherXvals'] = [1]
    # hold_series['xdatarange'] = 'Data'
    this_data['combine_Fit'] = hold_series

this_info = pa.Series()
this_info['save_file'] = '/home/jackdra/LQCD/Scripts/Jose_PoS/DataGraphs/Jacks/'+file_kappa+'_master_'+fit_info['Funs'][0].__name__+'.pdf'
if DoTree:
    this_info['title'] = r'$A(t_{f})/A(0)$ Div Tree ' + picked_kappa
else:
    this_info['title'] = r'$A(t_{f})/A(0)$'
if use_t0rat:
    this_info['xlabel'] = r'$\frac{\sqrt{8t_{f}}}{\sqrt{8t_{0}}}$'
else:
    this_info['xlabel'] = r'$\sqrt{8t_{f}}[fm]$'
this_info['ylabel'] = r'$A(t_{f})/A(0)$'
# this_info['xlims'] = [0,10]
# this_info['ylims'] = [0,15]
data_plot = Plotting(plot_data=this_data,plot_info=this_info)
data_plot.PlotAll()
data_plot.ShowFigure()
data_plot.PrintData()

if do_exraps:
    this_data = pa.DataFrame()
    hold_series = null_series
    hold_series['type'] = 'error_bar'
    # hold_series['key_select'] = (f'{mpi_dict["mpi410"]:.5f}',slice(None))
    hold_series['x_data'] = 'from_keys'
    plot_data = A0Fit_data['Fit_Best_Eval'].reset_index()
    plot_data[r'$m_{\pi}$ GeV'] = plot_data[r'$m_{\pi}$ GeV'].apply(lambda x : f'{x:.5f}')
    plot_data[r'$a^{2}/t_{0}$ fm'] = plot_data[r'$a^{2}/t_{0}$ fm'].apply(lambda x : f'{x:.3f}')
    a_plot_data = plot_data.set_index([r'$a^{2}/t_{0}$ fm'])['Fit_Best_Eval']
    hold_series['y_data'] = a_plot_data.apply(lambda x : x.Avg)
    hold_series['yerr_data'] = a_plot_data.apply(lambda x : x.Std)
    if use_t0rat:
        hold_series['label'] = 'Ens Results vs a' + r'$\frac{\sqrt{8t_{f}}}{\sqrt{8t_{0}}}='+f'{tflow_picked:.2f}'+r'$'
    else:
        hold_series['label'] = 'Ens Results vs a' + r'$\sqrt{8t_{f}}='+f'{tflow_picked:.2f}'+r'$'
    this_data['ens results vs 1_div_t0'] = hold_series

    hold_series = null_series
    hold_series['type'] = 'error_bar'
    # hold_series['key_select'] = (slice(None),)
    mpi_plot_data = plot_data.set_index([r'$m_{\pi}$ GeV'])['Fit_Best_Eval']
    hold_series['x_data'] = 'from_keys'
    hold_series['y_data'] = mpi_plot_data.apply(lambda x : x.Avg)
    hold_series['yerr_data'] = mpi_plot_data.apply(lambda x : x.Std)
    if use_t0rat:
        hold_series['label'] = 'Ens Results vs mpi' + r'$\frac{\sqrt{8t_{f}}}{\sqrt{8t_{0}}}='+f'{tflow_picked:.2f}'+r'$'
    else:
        hold_series['label'] = 'Ens Results vs mpi' + r'$\sqrt{8t_{f}}='+f'{tflow_picked:.2f}'+r'$'
    this_data['ens results vs mpi'] = hold_series

    if len(M_labels) > 1:
        hold_series = null_series
        hold_series['type'] = 'fit_vary'
        # plot_cont = combine_fit.ReduceViaChi(nchi_cutoff=20)[0]
        plot_cont = cont_fit
        # hold_series['key_select'] = plot_comb.Fit_Stats.index[0]
        hold_series['key_select'] = 'First'
        hold_series['fit_class'] = plot_cont
        hold_series['label'] = 'Continum Extrapolation'
        hold_series['xaxis'] = 0
        hold_series['otherXvals'] = [1]
        # hold_series['xdatarange'] = 'Data'
        this_data['combine_Fit'] = hold_series

    this_info = pa.Series()
    this_info['save_file'] = '/home/jackdra/LQCD/Scripts/Jose_PoS/DataGraphs/Jacks/'+file_kappa+'_cont_'+fit_info['Funs'][0].__name__+'.pdf'
    if DoTree:
        this_info['title'] = r'cont. extraps. Div Tree ' + picked_kappa
    else:
        this_info['title'] = r'cont. extraps.'
    this_info['xlabel'] = r'$m_{\pi} [GeV]$'
    this_info['ylabel'] = r'$val(m_{\pi},a^{2}/t_{0})$'
    # this_info['xlims'] = [0,10]
    # this_info['ylims'] = [0,15]
    data_plot = Plotting(plot_data=this_data,plot_info=this_info)
    data_plot.PlotAll()
    data_plot.ShowFigure()
    data_plot.PrintData()

    this_data = pa.DataFrame()
    if 'Par_a' in par_fit_list.keys():
        hold_series = null_series
        hold_series['type'] = 'fit_vary'
        # plot_cont = combine_fit.ReduceViaChi(nchi_cutoff=20)[0]
        # hold_series['key_select'] = plot_comb.Fit_Stats.index[0]
        hold_series['key_select'] = 'First'
        hold_series['fit_class'] = par_fit_list['Par_a']
        hold_series['label'] = 'param a'
        hold_series['xaxis'] = 0
        hold_series['otherXvals'] = [1]
        # hold_series['xdatarange'] = 'Data'
        this_data['fitpar_a'] = hold_series

        hold_series = null_series
        hold_series['type'] = 'error_bar'
        # hold_series['key_select'] = (f'{mpi_dict["mpi410"]:.5f}',slice(None))
        hold_series['x_data'] = 'from_keys'
        plot_data = A0Fit_data['Par_a'].reset_index()
        plot_data[r'$m_{\pi}$ GeV'] = plot_data[r'$m_{\pi}$ GeV'].apply(lambda x : f'{x:.5f}')
        plot_data[r'$a^{2}/t_{0}$ fm'] = plot_data[r'$a^{2}/t_{0}$ fm'].apply(lambda x : f'{x:.3f}')
        a_plot_data = plot_data.set_index([r'$a^{2}/t_{0}$ fm'])['Par_a']
        hold_series['y_data'] = a_plot_data.apply(lambda x : x.Avg)
        hold_series['yerr_data'] = a_plot_data.apply(lambda x : x.Std)
        hold_series['label'] = r'Ens Results vs $a^{2}/t_{0}$, par a'
        this_data['ens results vs 1_div_t0 par a'] = hold_series

        hold_series = null_series
        hold_series['type'] = 'error_bar'
        # hold_series['key_select'] = (slice(None),)
        mpi_plot_data = plot_data.set_index([r'$m_{\pi}$ GeV'])['Par_a']
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = mpi_plot_data.apply(lambda x : x.Avg)
        hold_series['yerr_data'] = mpi_plot_data.apply(lambda x : x.Std)
        hold_series['label'] = 'Ens Results vs mpi, par a'
        this_data['ens results vs mpi par a'] = hold_series


    if 'Par_b' in par_fit_list.keys():
        hold_series = null_series
        hold_series['type'] = 'fit_vary'
        # plot_cont = combine_fit.ReduceViaChi(nchi_cutoff=20)[0]
        # hold_series['key_select'] = plot_comb.Fit_Stats.index[0]
        hold_series['key_select'] = 'First'
        hold_series['fit_class'] = par_fit_list['Par_b']
        hold_series['label'] = 'param b'
        hold_series['xaxis'] = 0
        hold_series['otherXvals'] = [1]
        # hold_series['xdatarange'] = 'Data'
        this_data['fitpar_b'] = hold_series

        hold_series = null_series
        hold_series['type'] = 'error_bar'
        # hold_series['key_select'] = (f'{mpi_dict["mpi410"]:.5f}',slice(None))
        hold_series['x_data'] = 'from_keys'
        plot_data = A0Fit_data['Par_b'].reset_index()
        plot_data[r'$m_{\pi}$ GeV'] = plot_data[r'$m_{\pi}$ GeV'].apply(lambda x : f'{x:.5f}')
        plot_data[r'$a^{2}/t_{0}$ fm'] = plot_data[r'$a^{2}/t_{0}$ fm'].apply(lambda x : f'{x:.3f}')
        a_plot_data = plot_data.set_index([r'$a^{2}/t_{0}$ fm'])['Par_b']
        hold_series['y_data'] = a_plot_data.apply(lambda x : x.Avg)
        hold_series['yerr_data'] = a_plot_data.apply(lambda x : x.Std)
        hold_series['label'] = r'Ens Results vs $a^{2}/t_{0}$, par b'
        this_data['ens results vs 1_div_t0 par b'] = hold_series

        hold_series = null_series
        hold_series['type'] = 'error_bar'
        # hold_series['key_select'] = (slice(None),)
        mpi_plot_data = plot_data.set_index([r'$m_{\pi}$ GeV'])['Par_b']
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = mpi_plot_data.apply(lambda x : x.Avg)
        hold_series['yerr_data'] = mpi_plot_data.apply(lambda x : x.Std)
        hold_series['label'] = 'Ens Results vs mpi, par b'
        this_data['ens results vs mpi par b'] = hold_series


    if 'Par_c' in par_fit_list.keys():
        hold_series = null_series
        hold_series['type'] = 'fit_vary'
        # plot_cont = combine_fit.ReduceViaChi(nchi_cutoff=20)[0]
        # hold_series['key_select'] = plot_comb.Fit_Stats.index[0]
        hold_series['key_select'] = 'First'
        hold_series['fit_class'] = par_fit_list['Par_c']
        hold_series['label'] = 'param c'
        hold_series['xaxis'] = 0
        hold_series['otherXvals'] = [1]
        # hold_series['xdatarange'] = 'Data'
        this_data['fitpar_c'] = hold_series

        hold_series = null_series
        hold_series['type'] = 'error_bar'
        # hold_series['key_select'] = (f'{mpi_dict["mpi410"]:.5f}',slice(None))
        hold_series['x_data'] = 'from_keys'
        plot_data = A0Fit_data['Par_c'].reset_index()
        plot_data[r'$m_{\pi}$ GeV'] = plot_data[r'$m_{\pi}$ GeV'].apply(lambda x : f'{x:.5f}')
        plot_data[r'$a^{2}/t_{0}$ fm'] = plot_data[r'$a^{2}/t_{0}$ fm'].apply(lambda x : f'{x:.3f}')
        a_plot_data = plot_data.set_index([r'$a^{2}/t_{0}$ fm'])['Par_c']
        hold_series['y_data'] = a_plot_data.apply(lambda x : x.Avg)
        hold_series['yerr_data'] = a_plot_data.apply(lambda x : x.Std)
        hold_series['label'] = r'Ens Results vs $a^{2}/t_{0}$, par c'
        this_data['ens results vs 1_div_t0 par c'] = hold_series

        hold_series = null_series
        hold_series['type'] = 'error_bar'
        # hold_series['key_select'] = (slice(None),)
        mpi_plot_data = plot_data.set_index([r'$m_{\pi}$ GeV'])['Par_c']
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = mpi_plot_data.apply(lambda x : x.Avg)
        hold_series['yerr_data'] = mpi_plot_data.apply(lambda x : x.Std)
        hold_series['label'] = 'Ens Results vs mpi, par c'
        this_data['ens results vs mpi par c'] = hold_series

    if 'Par_d' in par_fit_list.keys():
        hold_series = null_series
        hold_series['type'] = 'fit_vary'
        # plot_cont = combine_fit.ReduceViaChi(nchi_cutoff=20)[0]
        # hold_series['key_select'] = plot_comb.Fit_Stats.index[0]
        hold_series['key_select'] = 'First'
        hold_series['fit_class'] = par_fit_list['Par_d']
        hold_series['label'] = 'param d'
        hold_series['xaxis'] = 0
        hold_series['otherXvals'] = [1]
        # hold_series['xdatarange'] = 'Data'
        this_data['fitpar_d'] = hold_series

        hold_series = null_series
        hold_series['type'] = 'error_bar'
        # hold_series['key_select'] = (f'{mpi_dict["mpi410"]:.5f}',slice(None))
        hold_series['x_data'] = 'from_keys'
        plot_data = A0Fit_data['Par_d'].reset_index()
        plot_data[r'$m_{\pi}$ GeV'] = plot_data[r'$m_{\pi}$ GeV'].apply(lambda x : f'{x:.5f}')
        plot_data[r'$a^{2}/t_{0}$ fm'] = plot_data[r'$a^{2}/t_{0}$ fm'].apply(lambda x : f'{x:.3f}')
        a_plot_data = plot_data.set_index([r'$a^{2}/t_{0}$ fm'])['Par_d']
        hold_series['y_data'] = a_plot_data.apply(lambda x : x.Avg)
        hold_series['yerr_data'] = a_plot_data.apply(lambda x : x.Std)
        hold_series['label'] = r'Ens Results vs $a^{2}/t_{0}$, par d'
        this_data[r'ens results vs t_0 par d'] = hold_series

        hold_series = null_series
        hold_series['type'] = 'error_bar'
        # hold_series['key_select'] = (slice(None),)
        mpi_plot_data = plot_data.set_index([r'$m_{\pi}$ GeV'])['Par_d']
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = mpi_plot_data.apply(lambda x : x.Avg)
        hold_series['yerr_data'] = mpi_plot_data.apply(lambda x : x.Std)
        hold_series['label'] = 'Ens Results vs mpi, par d'
        this_data['ens results vs mpi par d'] = hold_series


    this_info = pa.Series()
    this_info['save_file'] = '/home/jackdra/LQCD/Scripts/Jose_PoS/DataGraphs/Jacks/'+file_kappa+'_params_'+fit_info['Funs'][0].__name__+'.pdf'
    if DoTree:
        this_info['title'] = r'cont. extraps. Div Tree ' + picked_kappa
    else:
        this_info['title'] = r'cont. extraps.'
    this_info['xlabel'] = r'$m_{\pi} [GeV]$'
    this_info['ylabel'] = r'$val(m_{\pi},a^{2}/t_{0})$'
    # this_info['xlims'] = [0,10]
    # this_info['ylims'] = [0,15]
    data_plot = Plotting(plot_data=this_data,plot_info=this_info)
    data_plot.PlotAll()
    data_plot.ShowFigure()
    data_plot.PrintData()

this_data = pa.DataFrame()
for ikey,ival in Fit_data.items():
    hold_series = null_series
    hold_series['type'] = 'error_bar_vary'
    hold_series['key_select'] = (first_key[ikey],slice(None))
    hold_series['x_data'] = 'from_keys'
    hold_series['y_data'] = ival['E_GeV'].apply(lambda x : x.Avg)
    hold_series['yerr_data'] = ival['E_GeV'].apply(lambda x : x.Std)
    hold_series['label'] = ikey
    this_data[ikey] = hold_series
this_info = pa.Series()
this_info['save_file'] = '/home/jackdra/LQCD/Scripts/Jose_PoS/DataGraphs/Jacks/Mass_master.pdf'
this_info['title'] = r'Effective Mass'
this_info['xlabel'] = r'$t_{cut} [fm]$'
this_info['ylabel'] = r'$EffM[GeV]$'
# this_info['xlims'] = [0,10]
# this_info['ylims'] = [0,15]
data_plot = Plotting(plot_data=this_data,plot_info=this_info)
data_plot.PlotAll()
data_plot.ShowFigure()
data_plot.PrintData()

this_data = pa.DataFrame()
for ikey,ival in Fit_data.items():
    hold_series = null_series
    hold_series['type'] = 'error_bar_vary'
    hold_series['key_select'] = (first_key[ikey],slice(None))
    hold_series['x_data'] = 'from_keys'
    ivals = []
    for (itf,itc),this_ival in ival['E_GeV'].items():
        ivals.append(this_ival/ival['E_GeV'][('0.0',itc)])
        ivals[-1].Stats()
    this_series = pa.Series(ivals,index=ival['E_GeV'].index)
    hold_series['y_data'] = this_series.apply(lambda x : x.Avg)
    hold_series['yerr_data'] = this_series.apply(lambda x : x.Std)
    hold_series['label'] = ikey
    this_data[ikey] = hold_series
this_info = pa.Series()
this_info['save_file'] = '/home/jackdra/LQCD/Scripts/Jose_PoS/DataGraphs/Jacks/MM0_master.pdf'
this_info['title'] = r'Effective Mass Flowed Ratio'
this_info['xlabel'] = r'$t_{cut} [fm]$'
this_info['ylabel'] = r'$M(t_{f})/M(0)$'
# this_info['xlims'] = [0,10]
# this_info['ylims'] = [0,15]
data_plot = Plotting(plot_data=this_data,plot_info=this_info)
data_plot.PlotAll()
data_plot.ShowFigure()
data_plot.PrintData()

# import matplotlib.pyplot as pl
# import numpy as np
#
# xvals = np.mgrid[0.1:1:100j]
# # yvals = 1 +2.8/(np.log((xvals**2)/8)-0.11)
# yvals1 = c3FitFun_NoA([xvals],[1.08,3.65])
# yvals2 = c3FitFun_NoA([xvals],[1.28,3.65])
# yvals3 = c3FitFun_NoA([xvals],[1.08,2.65])
# yvals4 = c3FitFun_NoA([xvals],[1.28,2.65])
# pl.clf()
# pl.plot(xvals,yvals1,label='1')
# pl.plot(xvals,yvals2,label='2')
# pl.plot(xvals,yvals3,label='3')
# pl.plot(xvals,yvals4,label='4')
# pl.legend()
