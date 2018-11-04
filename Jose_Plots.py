#!/usr/bin/env python

import pandas as pa
import numpy as np
from PlotData import Plotting,null_series
from MomParams import hbarc
from BootStrapping import BootStrap
from SetsOfFits import SetOfFitFuns
from PredefFitFuns import c3FitFun,C2OneStateFitFunNoExp,C2TwoStateFitFunNoExp,c3FitFun_V2
import os
import pickle as pik

this_folder = '/home/jackdra/LQCD/Scripts/Jose_PoS/DataGraphs/FitRes/'
temp_folder = '/home/jackdra/temp/'

Fit_ens_files = {}
Fit_ens_files['M'] = ['DataFrame_C_32X64_1.csv','DataFrame_C_32X64_2.csv','DataFrame_C_32X64_3.csv']
Fit_ens_files['A'] = ['DataFrame_C_16X32.csv','DataFrame_C_20X40.csv','DataFrame_C_28X56.csv']

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

def MakeFileList(this_dict):
    output_dict = {}
    for ifile,ilab in zip(this_dict['M'],M_labels):
        output_dict[ilab] = this_folder + ifile
    for ifile,ilab in zip(this_dict['A'],A_labels):
        output_dict[ilab] = this_folder + ifile
    return output_dict

Fit_filelist = MakeFileList(Fit_ens_files)

def ReadFileList(this_dict):
    output_dict = {}
    for ikey,ifile in this_dict.items():
        with open(ifile,'r') as f:
            output_dict[ikey] = pa.read_csv(f)
    return output_dict

Fit_data = ReadFileList(Fit_filelist)

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
        Fit_data[ikey]['tmin_fm'] = Fit_data[ikey]['tmin'].apply(lambda x : str(x*a_dict[ikey]))
        Fit_data[ikey]['E_GeV'] = Fit_data[ikey]['E'].apply(lambda x : x*hbarc/a_dict[ikey])
        Fit_data[ikey]['E_GeV'].apply(lambda x : x.Stats())
        first_key[ikey] = Fit_data[ikey]['sqrt8t_fm'].iloc[0]
        second_key[ikey] = Fit_data[ikey]['tmin_fm'].iloc[0]
        Fit_data[ikey] = Fit_data[ikey].set_index(['sqrt8t_fm','tmin_fm'])
        with open(data_file,'wb') as f:
            pik.dump((Fit_data[ikey],first_key[ikey],second_key[ikey]),f)
A0Fit_data = {}
fit_info = {}
# fit_info['Funs'] = [c3FitFun,3]
fit_info['Funs'] = [c3FitFun_V2,3]
# fit_info['Funs'] = [C2OneStateFitFunNoExp,2]
# fit_info['iGuess'] = [1,1]
fit_tmin_min = 1
fit_tmin_max = 1.5
tflow_fit_min = 0.2
tflow_fit_max = 1.5
min_par_tflow = 5
min_par_tmin = 0
Master_Fit_Data = pa.DataFrame()
for ikey,ival in Fit_data.items():
    fit_file = temp_folder+'/fit_'+ikey+'fitr'+'-'.join([str(int(fit_tmin_min*10)),
                                                         str(int(fit_tmin_max*10))])+'_'+fit_info['Funs'][0].__name__+'.p'
    if os.path.isfile(fit_file):
        with open(fit_file,'rb') as f:
            this_fit,this_df = pik.load(f)
            this_fit.GetFuns()
    else:
        # ival = ival.reset_index()
        # ival = ival[ival['tmin_fm'].apply(lambda x : float(x)) >fit_tmin_min]
        # ival = ival[fit_tmin_max>ival['tmin_fm'].apply(lambda x : float(x))]
        # ival = ival.set_index(['sqrt8t_fm','tmin_fm'])
        # ival = pa.Series(ival.values,index=ival.index)
        this_fit = SetOfFitFuns(data = ival['A(t)/A(0)'],name=ikey)
        # fit_info['Funs'] = [ConstantFitFun,1]
        # fit_info['iGuess'] = [0,0]
        this_fit.ScanBox_FixEnd(tflow_fit_min,tflow_fit_max,
                                fit_tmin_min,fit_tmin_max,
                                fit_info,min_fit_len_1=min_par_tflow,min_fit_len_2=min_par_tmin,
                                from_xdata=True)
        this_fit.DoFits()
        this_fit.SortChi()
        this_df = ival.reset_index()[['sqrt8t_fm','tmin_fm','A(t)/A(0)']]
        this_df['ens'] = ikey
        # this_df = this_df.set_index(['sqrt8t_fm','tmin_fm'])
        with open(fit_file,'wb') as f:
            this_fit.RemoveFuns()
            pik.dump((this_fit,this_df),f)
            this_fit.GetFuns()
    A0Fit_data[ikey] = this_fit
    Master_Fit_Data = Master_Fit_Data.append(this_df)
Master_Fit_Data = Master_Fit_Data.set_index(['sqrt8t_fm','tmin_fm']).sort_index()
# if 'sqrt8t_fm' in Master_Fit_Data.columns:
#     del Master_Fit_Data['sqrt8t_fm']
# if 'tmin_fm' in Master_Fit_Data.columns:
#     del Master_Fit_Data['tmin_fm']
fit_file = temp_folder+'/fit_combine_fitr'+'-'.join([str(int(fit_tmin_min*10)),
                                                     str(int(fit_tmin_max*10))])+'_'+fit_info['Funs'][0].__name__+'.p'

# Master_Fit_Data
# Master_Fit_Data = Master_Fit_Data.reset_index()
# Master_Fit_Data['tmin_fm'].apply(lambda x : float(x)) >fit_tmin_min
# Master_Fit_Data = Master_Fit_Data[Master_Fit_Data['tmin_fm'].apply(lambda x : float(x)) >fit_tmin_min]
# Master_Fit_Data = Master_Fit_Data[fit_tmin_max>Master_Fit_Data['tmin_fm'].apply(lambda x : float(x))]

if os.path.isfile(fit_file):
    with open(fit_file,'rb') as f:
        combine_fit = pik.load(f)
        combine_fit.GetFuns()
else:
    # Master_Fit_Data = Master_Fit_Data.reset_index()
    # Master_Fit_Data = Master_Fit_Data[Master_Fit_Data['tmin_fm'].apply(lambda x : float(x)) >fit_tmin_min]
    # Master_Fit_Data = Master_Fit_Data[fit_tmin_max>Master_Fit_Data['tmin_fm'].apply(lambda x : float(x))]
    # Master_Fit_Data = Master_Fit_Data.set_index(['sqrt8t_fm','tmin_fm'])
    # Master_Fit_Data = pa.Series(Master_Fit_Data['A(t)/A(0)'].values,index=Master_Fit_Data.index)
    combine_fit = SetOfFitFuns(data = Master_Fit_Data['A(t)/A(0)'],name='combined')
    # fit_info['Funs'] = [ConstantFitFun,1]
    # fit_info['iGuess'] = [0,0]
    combine_fit.ScanBox_FixEnd( tflow_fit_min,tflow_fit_max,
                                fit_tmin_min,fit_tmin_max,
                                fit_info,min_fit_len_1=min_par_tflow+5,min_fit_len_2=min_par_tmin,
                                from_xdata=True)
    # combine_fit.ScanBox(tflow_fit_min,nt_flow,0,ntmin,fit_info,min_fit_len_1=5,min_fit_len_2=5)
    combine_fit.DoFits()
    combine_fit.SortChi()
    with open(fit_file,'wb') as f:
        combine_fit.RemoveFuns()
        pik.dump(combine_fit,f)
        combine_fit.GetFuns()

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
    plot_fit = A0Fit_data[ikey]
    hold_series['key_select'] = 'First'
    hold_series['fit_class'] = plot_fit
    hold_series['label'] = ikey+' Fit'
    hold_series['xaxis'] = 0
    # hold_series['xdatarange'] = 'Data'
    hold_series['otherXvals'] = [1]
    this_data[ikey+'_Fit'] = hold_series

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
this_info['save_file'] = '/home/jackdra/LQCD/Scripts/Jose_PoS/DataGraphs/Jacks/AA0_master.pdf'
this_info['title'] = r'$A(t_{f})/A(0)$'
this_info['xlabel'] = r'$\sqrt{8t_{f}}[fm]$'
this_info['ylabel'] = r'$A(t_{f})/A(0)$'
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
