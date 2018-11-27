#!/usr/bin/env python

import pandas as pa
import numpy as np
from PlotData import Plotting,null_series
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
from ReadBinaryCfuns import RC2Full

import os
import pickle as pik

this_folder = '/home/jackdra/LQCD/Scripts/Jose_PoS/DataGraphs/FitRes/'
temp_folder = '/home/jackdra/temp/'

Fit_ens_files = {}
Fit_ens_files['M'] = ['DataFrame_C_32X64_1.csv','DataFrame_C_32X64_2.csv','DataFrame_C_32X64_3.csv']
Fit_ens_files['A'] = ['DataFrame_C_16X32.csv','DataFrame_C_20X40.csv','DataFrame_C_28X56.csv']

do_exraps = True
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

from XmlFormatting import MakeValAndErr
for ival_err,ival in zip(t0_dict_err.values(),t0_dict.values()):
    print(MakeValAndErr(ival,ival_err,Dec=2))
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
flow_list += [f'{val:.1f}' for val in np.arange(0,10,1)]
def ReadTLFileList(this_dict):
    output_dict = {}
    for ikey,ifile in this_dict.items():
        tf_data = []
        t_list,tf_list = [],[]
        for iflow in flow_list:
            this_file = ifile.replace('FLOW_TIME',iflow)
            this_data = RC2Full([this_file],['0 0 0'],
                                InterpNumb='15',
                                MesOrBar='Meson').data
            this_data = list(np.array(this_data)[0,:,0].flatten())
            tf_data += this_data
            t_list += [f'{it}' for it in this_data]
            tf_list += [iflow for it in this_data]
        output_dict[ikey] = pa.DataFrame()
        output_dict[ikey]['t_sink'] = pa.Series(t_list)
        output_dict[ikey]['t_flow'] = pa.Series(tf_list)
        output_dict[ikey]['Corr'] = pa.Series(tf_data)
    return output_dict

def ReadFileList(this_dict):
    output_dict = {}
    for ikey,ifile in this_dict.items():
        with open(ifile,'r') as f:
            output_dict[ikey] = pa.read_csv(f)
    return output_dict

Fit_data = ReadFileList(Fit_filelist)
this_kappa = 'k0.1250'
tree_files = {}
for ikey in Fit_filelist.keys():
    if 'mpi' in ikey:
        tree_files[ikey] =        '/home/jackdra/LQCD/Scripts/Jose_PoS/DataGraphs/Treelev/RC32x64_tree/RC32x64_tree_'+this_kappa+'_FLOW_TIME.xml'
    else:
        this_RC = int(ikey[1:])
        this_RC = 'RC'+str(this_RC)+'x'+str(this_RC*2)
        tree_files[ikey] =        '/home/jackdra/LQCD/Scripts/Jose_PoS/DataGraphs/Treelev/'+this_RC+'_tree/'+this_RC+'_tree_'+this_kappa+'_FLOW_TIME.xml'


# tree_data = ReadTLFileList(tree_files)


def BootStrap_Vals(this_DF,this_tree= None):
    A_list,E_list,chi_list,ilist = [],[],[],[]
    for ikey,iDF in this_DF.groupby(('flow t','tmin')):
        if this_tree is None:
            A_list.append(BootStrap(thisnboot=len(iDF),bootvals=iDF['A'].values))
            E_list.append(BootStrap(thisnboot=len(iDF),bootvals=iDF['E'].values))
            chi_list.append(BootStrap(thisnboot=len(iDF),bootvals=iDF['chisq'].values))
        else:
            A_list.append(BootStrap(thisnboot=len(iDF)+1,
                                    bootvals=np.append(iDF['A'].values,this_tree['A'])))
            E_list.append(BootStrap(thisnboot=len(iDF)+1,
                                    bootvals=np.append(iDF['E'].values,this_tree['E'])))
            chi_list.append(BootStrap(thisnboot=len(iDF),bootvals=iDF['chisq'].values))
        ilist.append(ikey)
    if len(ilist) > 0:
        ilist = pa.MultiIndex.from_tuples(ilist,names=('flow t','tmin'))
        this_DF = pa.DataFrame()
        this_DF['A'] = pa.Series(A_list,index = ilist)
        this_DF['E'] = pa.Series(E_list,index = ilist)
        this_DF['chisq'] = pa.Series(chi_list,index = ilist)
    return this_DF.reset_index()

def MakeAdiv(this_DF):
    out_list,ilist = [],[]
    for it,ival in this_DF.groupby('tmin'):
        DF_tf0 = ival['A'][ival['flow t'] == 0].values[0]
        for jkey,jval in ival['A'].items():
            outval = jval.__div__(DF_tf0)
            outval.Stats()
            out_list.append(outval)
            ilist.append(jkey)
    if len(ilist) > 0:
        return pa.Series(out_list,index=ilist)
    else:
        return pa.Series()

first_key,second_key = {},{}
for ikey,ival in Fit_data.items():
    data_file = temp_folder+'/'+ikey+'_data.p'
    if os.path.isfile(data_file):
        with open(data_file,'rb') as f:
            Fit_data[ikey],first_key[ikey],second_key[ikey] = pik.load(f)
    else:
        Fit_data[ikey] = BootStrap_Vals(ival)
        Fit_data[ikey]['A(t)/A(0)'] = MakeAdiv(Fit_data[ikey])
        Fit_data[ikey]['sqrt8t_fm'] = Fit_data[ikey]['flow t'].apply(lambda x : str(np.sqrt(8*x)*a_dict[ikey]))
        Fit_data[ikey]['sqrt8t_t0rat'] = Fit_data[ikey]['flow t'].apply(lambda x : str(np.sqrt(x/t0_dict[ikey])))
        Fit_data[ikey]['tmin_fm'] = Fit_data[ikey]['tmin'].apply(lambda x : str(x*a_dict[ikey]))
        Fit_data[ikey]['E_GeV'] = Fit_data[ikey]['E'].apply(lambda x : x*hbarc/a_dict[ikey])
        Fit_data[ikey]['E_GeV'].apply(lambda x : x.Stats())
        first_key[ikey] = Fit_data[ikey][sqrt8t_key].iloc[0]
        second_key[ikey] = Fit_data[ikey]['tmin_fm'].iloc[0]
        Fit_data[ikey] = Fit_data[ikey].set_index([sqrt8t_key,'tmin_fm'])
        with open(data_file,'wb') as f:
            pik.dump((Fit_data[ikey],first_key[ikey],second_key[ikey]),f)
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
        tflow_fit_min[ikey] = min_hold/(a_dict[ikey]*np.sqrt(8*t0_dict[ikey]))
        tflow_fit_max[ikey] = max_hold/(a_dict[ikey]*np.sqrt(8*t0_dict[ikey]))
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

mpi_dict = {}
for ikey,ival in Fit_data.items():
    fit_file = temp_folder+'/fit_'+ikey+'fitr'+'-'.join([str(int(fit_tmin_min*10)),
                                                         str(int(fit_tmin_max*10))])+'_'+fit_info['Funs'][0].__name__+'.p'
    if os.path.isfile(fit_file):
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
        this_fit = SetOfFitFuns(data = ival['A(t)/A(0)'],name=ikey)
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
        this_df = ival.reset_index()[[sqrt8t_key,'tmin_fm','A(t)/A(0)']]
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
fit_file = temp_folder+'/fit_combine_fitr'+'-'.join([str(int(fit_tmin_min*10)),
                                                     str(int(fit_tmin_max*10))])+'_'+fit_info['Funs'][0].__name__+'.p'

# Master_Fit_Data
# Master_Fit_Data = Master_Fit_Data.reset_index()
# Master_Fit_Data['tmin_fm'].apply(lambda x : float(x)) >fit_tmin_min
# Master_Fit_Data = Master_Fit_Data[Master_Fit_Data['tmin_fm'].apply(lambda x : float(x)) >fit_tmin_min]
# Master_Fit_Data = Master_Fit_Data[fit_tmin_max>Master_Fit_Data['tmin_fm'].apply(lambda x : float(x))]
if do_exraps:
    par_fit_list = {}
    if os.path.isfile(fit_file):
        with open(fit_file,'rb') as f:
            if len(M_labels) > 1:
                combine_fit,cont_fit,par_fit_list = pik.load(f)
                cont_fit.GetFuns()
            else:
                combine_fit,par_fit_list = pik.load(f)
    else:
        combine_fit = SetOfFitFuns(data = Master_Fit_Data['A(t)/A(0)'],name='combined')
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
    hold_series['y_data'] = ival['A(t)/A(0)'].apply(lambda x : x.Avg)
    hold_series['yerr_data'] = ival['A(t)/A(0)'].apply(lambda x : x.Std)
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
this_info['save_file'] = '/home/jackdra/LQCD/Scripts/Jose_PoS/DataGraphs/Jacks/AA0_master_'+fit_info['Funs'][0].__name__+'.pdf'
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
    this_info['save_file'] = '/home/jackdra/LQCD/Scripts/Jose_PoS/DataGraphs/Jacks/AA0_cont_'+fit_info['Funs'][0].__name__+'.pdf'
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
    this_info['save_file'] = '/home/jackdra/LQCD/Scripts/Jose_PoS/DataGraphs/Jacks/AA0_params_'+fit_info['Funs'][0].__name__+'.pdf'
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
    for (itf,itc),this_ival in ival['E_GeV'].iteritems():
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
