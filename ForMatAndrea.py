#!/usr/bin/env python


weinberg_file = '/home/jackdra/PHD/CHROMA/TestVar/Scripts/Python_Analysis/Configs/FlowOps/RC32x64_kud1370000_ks1364000_-b-_Weinberg.cfgs.msg'

top_file = '/home/jackdra/PHD/CHROMA/TestVar/Scripts/Python_Analysis/Configs/FlowOps/RC32x64_kud1370000_ks1364000_-b-_TopCharge.cfgs.msg'

from XmlFormatting import untflowstr
from FileIO import ReadWithMeta
from BootStrapping import BootStrap
import numpy as np
import pandas as pa
from PlotData import Plotting

W_header,W_cfgs = ReadWithMeta(weinberg_file)
Q_header,Q_cfgs = ReadWithMeta(top_file)

tflow_list = np.array([untflowstr(itf) for itf in W_header])
t_f501_index = list(Q_header).index('t_f5.01')
t_f021_index = list(W_header).index('t_f0.21')
W_021 = W_cfgs['Op'].apply(lambda x : x[t_f021_index])
Q_501 = Q_cfgs['Op'].apply(lambda x : x[t_f501_index])

denom = BootStrap(666,cfgvals=W_021*Q_501)


avgl,stdl,indexl = [],[],[]
for ictflow,itflow in enumerate(tflow_list):
    W_t = W_cfgs['Op'].apply(lambda x : x[ictflow])
    W_data = BootStrap(666,cfgvals=W_t*W_021)
    # comb_data = tflow_list*W_data/denom
    comb_data = W_data/denom
    comb_data = itflow*comb_data
    comb_data.Stats()
    stdl.append(comb_data.Std)
    avgl.append(comb_data.Avg)
    indexl.append(itflow)

xdata = indexl
yseries = pa.Series(avgl)
yerrseries = pa.Series(stdl)

this_data = pa.DataFrame()
hold_series = pa.Series()
hold_series['type'] = 'error_bar'
hold_series['x_data'] = xdata
hold_series['y_data'] = yseries
hold_series['yerr_data'] = yerrseries
hold_series['label'] = 'Data'
hold_series['color'] = 'blue'
hold_series['shift'] = 0.0
hold_series['symbol'] = 'o'
this_data['test_1'] = hold_series

this_info = pa.Series()
this_info['save_file'] = './TestGraphs/ForMatAndrea.pdf'
this_info['title'] = r'$\frac{t_f<W(t_f)W(0.21)>}{<Q(5.01) W(0.21)>}$'
this_info['title_pad'] = 20
this_info['xlabel'] = r'$t_f$'
this_info['ylabel'] = r'Ratio'

data_plot = Plotting(plot_data=this_data,plot_info=this_info)
data_plot.PlotAll()
data_plot.ShowFigure()
data_plot.PrintData()
