#!/usr/bin/env python

import pickle as pik
import matplotlib.pyplot as pl
import numpy as np
import pandas as pa
import os
import PlotData as jpl
from XmlFormatting import MakeValAndErr
from Params import scratchdir,graphdir
from MiscFuns import mkdir_p

master_ens_list = ['mpi411','mpi570','mpi701','L16','L20','L28','qu_L24']
this_filelist = [scratchdir + '/fit_WQ_'+iens+'.py3p' for iens in master_ens_list]
block_flags = ['','_nb2','_nb3']
# this_file = '/home/jackdra/LQCD/Results/DebugResults/graphs/flowops_test_NoE_t_s0.36_quL24.tex.py3p'

this_info = pa.Series()
mat_graph_folder = graphdir+'MatHack/'
mkdir_p(mat_graph_folder)
this_info['save_file'] = mat_graph_folder+'FitrComp.pdf'
this_info['title'] = r'fitr comp'
this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
this_info['ylabel'] = r'ratio'
data_plot = jpl.Plotting(plot_info=this_info)
for iens,ifile in zip(master_ens_list,this_filelist):
    for iblock in block_flags:
        this_file = ifile.replace('.py3p',iblock+'.py3p')
        if os.path.isfile(this_file):
            print('Reading: ',this_file)
            with open(this_file,'rb') as f:
                fit_data,dump = pik.load(f)
            data_plot = fit_data.PlotVaryFitr(data_plot)
data_plot.PlotAll()
# data_plot.ShowFigure()
# data_plot.PrintData()



# this_fit.Fitr_Len_Sort()
# this_fit.Cut_max(30)
# # Cut_min_len(this_fit,20)
# # Cut_min_len(this_fit,40)
# fmt_fit = this_fit.Get_Formatted_Table()
# fmt_fit
# fmt_fit.to_latex(escape=False)
# this_fit.Fit_Stats_fmt.reset_index()[['Chi2pdf']]
# this_data = pa.DataFrame()
# hold_series = null_series
# hold_series['type'] = 'error_bar'
# hold_series['x_data'] = 'from_keys'
# hold_series['y_data'] = fit_data.apply(lambda x : x.Avg)
# hold_series['yerr_data'] = fit_data.apply(lambda x : x.Std)
# hold_series['label'] = r'$data$'
# this_data['data'] = hold_series
#
# hold_series = null_series
# hold_series['type'] = 'fit_vary'
# hold_series['key_select'] = 'First'
# hold_series['fit_class'] = this_fit
# hold_series['xaxis'] = 0
# hold_series['label'] = r'$Fit$'
# this_data['Fit'] = hold_series
#
# this_info = pa.Series()
# this_info['save_file'] = '/home/jackdra/LQCD/Scripts/Python_Analysis/TestGraphs/L28_Rat.pdf'
# this_info['title'] = r'$\frac{t<Q(t_{s})W(t)>}{<Q(t_{s})^{2}>}$'
# this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
# this_info['ylabel'] = r'R'
# # this_info['xlims'] = [0,10]
# # this_info['ylims'] = [0,15]
# data_plot = jpl.Plotting(plot_data=this_data,plot_info=this_info)
# data_plot.PlotAll()
# data_plot.ShowFigure()
# data_plot.PrintData()
