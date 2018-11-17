#!/usr/bin/env python


# import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
import matplotlib.pyplot as pl
import numpy as np
from Params import plot_style,fillalpha
import pandas as pa
from BootStrapping import BootStrap
from collections import OrderedDict
import SetsOfFits as sff
from MiscFuns import Series_TO_ODict,ODNested,get_val_float_0,fmt_Qt5_col,mkdir_p,DEBUGprint
from XmlFormatting import AvgStdToFormat,minfmt,maxfmt,MakeValAndErr,MakeVal
from FileIO import WriteXml,WritePickle,ReadPickleWrap
from copy import copy
from XmlFormatting import LegendFmt

Py2Read = False


import os
# import seaborn as sns
#
# sns.set_style("darkgrid",{'legend.frameon':True,'errorbar.capsize':5})

def tuple_wrap(val):
    if isinstance(val,(list,np.ndarray)):
        if len(val) == 1:
            return val
        else:
            return tuple(val)
    elif isinstance(val,tuple):
        if len(val) == 1:
            return val[0]
        else:
            return val
    else:
        return val

tpad,tsize,xsize,ysize = 20,40,60,60
legcol,legsize,legalpha = 1,25,0.7
linewidth,dotsize,errorbarlen = 2,10,6
def_tick_size = 30
def_grid = True
# linewidth,dotsize,errorbarlen = 4,20,12
# def_tick_size = 55
# def_grid = False

## Move Somewhere?
params = {'legend.fontsize': 25,
          'xtick.labelsize':def_tick_size,
          'ytick.labelsize':def_tick_size,
          'legend.numpoints': 1,
          'axes.labelsize' : 35,
          'axes.titlesize' : 35,
          'figure.autolayout': True,
          'axes.grid': def_grid,
          'errorbar.capsize': errorbarlen,
          'lines.markersize': dotsize,
          'lines.linewidth': linewidth,
          'font.size': 30,
          # 'lines.markeredgewidth':2,
          'axes.xmargin':0.01,
          'axes.ymargin':0.01}
# params = {'errorbar.capsize': 5}
pl.style.use(plot_style)
pl.rcParams.update(params)


colourset8 = [ '#000080','#B22222','#00B800','#8B008B', '#0000FF','#FF0000','#000000','#32CD32','#FF0066']
markerset = ['o','s','^','v','d']
shiftmax = 3
shiftper = 0.005 ##R function shift
shiftset = [0]
for ish in np.arange(1,shiftmax+1): shiftset += [-ish*shiftper,ish*shiftper]
del ish, shiftmax,shiftper

null_series = pa.Series()
null_series['boot_data'] = None
null_series['x_data'] = None
null_series['x_fun'] = None
null_series['y_data'] = None
null_series['yerr_data'] = None
null_series['xerr_data'] = None
null_series['x_data_median'] = None
null_series['y_data_median'] = None
null_series['yerr_data_up'] = None
null_series['xerr_data_up'] = None
null_series['yerr_data_down'] = None
null_series['xerr_data_down'] = None
null_series['key_select'] = None
null_series['xerr_data'] = None
null_series['type'] = None
null_series['fit_class'] = None
null_series['label'] = ' '
null_series['symbol'] = 'o'
null_series['color'] = None
null_series['shift'] = 0
null_series['shift_overlap'] = 0
null_series['xdatarange'] = 'Data'
null_series['scale'] = None
null_series['Phys'] = None
null_series['ShowPar'] = None
null_series['ShowEval'] = None
null_series['Median'] = None
null_series['suppress_key_index'] = False
null_series['fmt_class'] = None
null_series['supres_legend'] = None
null_series['x_scale'] = 1
null_series['x_range_min'] = None
null_series['x_range_max'] = None
null_series['x_increment'] = None
null_series['hairline'] = False
null_series['extrap_fade'] = 0
null_series['hair_alpha'] = False
null_series['phys_axies'] = False
null_series['FPName'] = ''
null_series['FPUnits'] = ''
null_series['plot_err'] = ''
null_series['arrow_text_size'] = 30
null_series['arrow_color'] = 'black'
null_series['arrow_style'] = 'fancy'


def Fix_Se_Keys(this_key,this_Se):
    new_i = []
    for iindex in this_Se.index:
        temp = []
        for ikey,df_key in zip(this_key,iindex):
            if isinstance(df_key,np.bytes_):
                df_key = df_key.decode('UTF-8')
            temp.append(df_key)
        new_i.append(tuple(temp))
    if len(new_i) > 0:
        name_l = [str(iname) for iname in this_Se.index.names]
        this_Se.index = pa.MultiIndex.from_tuples(new_i,names=name_l)
    out_key = []
    for ikey in this_key:
        if isinstance(ikey,np.bytes_):
            ikey = ikey.decode('UTF-8')
        out_key.append(ikey.replace('b\'','').replace('\'',''))
    return tuple(out_key),this_Se
#
# def Fix_Se_Keys(this_key,this_Se):
#     return this_Se

## this only works if they are in order, Add wrapper to this where it is sorted and unsorted!
def Shift_Duplicates(x_data,this_shiftper=0.01):
    if len(x_data) == 0:
        return x_data
    list_len = max(x_data) - min(x_data)
    this_shiftper = this_shiftper*list_len ##R function shift
    this_shiftset = [0]
    for ish in np.arange(1,50):
        this_shiftset += [-ish*this_shiftper,ish*this_shiftper]
    prev_shifted,shift_index = False,0
    x_out = []
    for ic,ix in enumerate(x_data):
        if ix in list(x_data[:ic]) + list(x_data[(ic + 1):]):
            if ix != prev_shifted:
                prev_shifted = ix
                shift_index = 0
            else:
                shift_index += 1
            x_out.append(ix+this_shiftset[shift_index])
        else:
            x_out.append(ix)
    if isinstance(x_data,np.ndarray):
        return np.array(x_out)
    else:
        return type(x_data)(x_out)

if Py2Read:
    def Un_Fix_Key_Select(this_key_select):
        this_key_select = list(this_key_select)
        for ic,ikey in enumerate(this_key_select):
            if not isinstance(ikey,slice) and 't_f' in ikey:
                # this_key_select[ic] = bytes(ikey.replace('b\'','').replace('\'',''),'ascii').decode('ascii')
                this_key_select[ic] = bytes(ikey.replace('b\'','').replace('\'',''),'ascii')
        return tuple(this_key_select)

    def Fix_Key_Select(this_key_select,this_data,add_b=True):
        first_index = list(this_data.index)[0]
        key_types = map(type,first_index)
        this_key_select = list(this_key_select)
        for ic,(ikey,ikey_sel) in enumerate(zip(key_types,this_key_select)):
            if not isinstance(ikey_sel,slice):
                if not isinstance(ikey_sel,ikey):
                    try:
                        this_key_select[ic] = ikey(this_key_select[ic])
                    except:
                        if 'b\'' in str(list(this_data.index)[0][ic]) and add_b:
                            this_key_select[ic] = ikey('b\''+this_key_select[ic].replace('b\'','').replace('\'','')+'\'','ascii')
                        else:
                            this_key_select[ic] = ikey(this_key_select[ic],'ascii')
        return tuple(this_key_select)
else:
    def Fix_Key_Select(this_key_select,this_data,add_b=True):
        return this_key_select

    def Un_Fix_Key_Select(this_key_select):
        return this_key_select



class Plotting(object):

    """
    Class to organise and duplicate the plotting.

    """


    def __init__(self, plot_data = pa.DataFrame(), plot_info = pa.Series()):

        ## currenlty only implementing single plot per plane
        self.MakePlane()


        """
        plot_data must be pandas DataFrame or list of DataFrames
        if dictionary passed, dictionary will be converted into DataFrame.

        Contains:
        columns:
            - number of data points
        keys:
            x_data, y_data, (xerr_data, yerr_data,)
            type, (errorbar, plot, hist etc..),
            label, color, symbol

        plotting a fit is done through:
        type = 'fit'
        fit_class = instance of Fitting from FitFunctions.py

        """
        self.ImportData(plot_data)

        """
        plot_info must be pandas Series, with elements:

        title, x_label, y_label, etc..

        """

        self.ImportInfo(plot_info)
        self.AddNulls()

        self.colcyc = iter(colourset8)
        self.symcyc = iter(markerset)
        self.shiftcyc = iter(shiftset)

    def ImportInfo(self,plot_info):
        if isinstance(plot_info,pa.Series):
            self.plot_info = plot_info
        elif isinstance(plot_info,(dict,OrderedDict)):
            self.plot_info = pa.Series()
            for ikey,ival in plot_info.items():
                self.plot_info[ikey] = ival
        else:
            raise IOError('pandas Series for plot_info not passed into plotting class')
        if 'style' in plot_info:
            pl.style.use(plot_info['style'])
        if 'leg_loc' not in plot_info:
            self.plot_info['leg_loc'] = 'best'
        if 'leg_ncol' not in plot_info:
            self.plot_info['leg_ncol'] = legcol
        if 'leg_fontsize' not in plot_info:
            self.plot_info['leg_fontsize'] = legsize
        if 'leg_alpha' not in plot_info:
            self.plot_info['leg_alpha'] = legalpha
        if 'y_plot_scale' not in plot_info:
            self.plot_info['y_plot_scale'] = 'None'
        if 'x_plot_scale' not in plot_info:
            self.plot_info['x_plot_scale'] = 'None'
        if 'save_file' in self.plot_info:
            self.HumanFile = self.plot_info['save_file'].replace('.pdf','')+'.xml'
            self.PickleFile = self.plot_info['save_file'].replace('.pdf','')+'.py3p'
            if '.pdf' not in self.plot_info['save_file']:
                self.plot_info['save_file'] = self.plot_info['save_file']+'.pdf'

    def UpdateInfo(self,plot_info):
        self.plot_info.update(plot_info)
        self.ImportInfo(self.plot_info)

    def GetFuns(self):
        if not hasattr(self,"plot_data"):
            raise EnvironmentError('import data before fixing functions')
        for icol,col_data in self.plot_data.items():
            if 'fit_class' in list(col_data.keys()):
                if hasattr(col_data['fit_class'],'GetFuns'):
                    self.plot_data[icol]['fit_class'].GetFuns()
                elif isinstance(col_data['fit_class'],pa.Series):
                    for ikey,ival in col_data['fit_class'].items():
                        if hasattr(ival,'GetFuns'):
                            ival.GetFuns()
                            # self.plot_data[icol]['fit_class'][ikey].GetFuns()
                            # if isinstance(ival.Fun,str) and 'FNF' in ival.Fun:
                            if isinstance(ival.Fun,str):
                                raise EnvironmentError('fit '+ival.Fun)
                            # if ival.FunDer == 'file_not_found':
                            #     raise EnvironmentError('FNF for fitting derivative: \n' + ival.FunDer_name)

    def RemoveFuns(self):
        if not hasattr(self,"plot_data"):
            raise EnvironmentError('import data before fixing functions')
        for icol,col_data in self.plot_data.items():
            if 'fit_class' in list(col_data.keys()):
                if hasattr(col_data['fit_class'],'RemoveFuns'):
                    self.plot_data[icol]['fit_class'].RemoveFuns()
                elif isinstance(col_data['fit_class'],pa.Series):
                    for ikey,ival in col_data['fit_class'].items():
                        if hasattr(ival,'RemoveFuns'):
                            # self.plot_data[icol]['fit_class'][ikey].RemoveFuns()
                            ival.RemoveFuns()

    def GetAspectRatio(self):
        return self.plot_plane.get_aspect()

    def SetAspectRatio(self,ratio='auto'):
        self.plot_plane.set_aspect(ratio)


    def ImportData(self,plot_data):
        if isinstance(plot_data,pa.DataFrame):
            self.plot_data = plot_data
            # print 'DEBUG'
            # for ikey,iplot in self.plot_data.iteritems():
            #     print ikey, iplot['type']
            #     print '\n   '.join(iplot.keys())
            #     if 'y_data' in iplot.keys():
            #         print iplot['y_data']
        elif isinstance(plot_data,pa.Series):
            self.plot_data = plot_data.to_frame(name='data')
        elif isinstance(plot_data,(dict,OrderedDict)):
            self.plot_data = pa.DataFrame()
            for ikey,ival in plot_data.items():
                self.plot_data[ikey] = ival
        elif isinstance(plot_data,(tuple,list,np.ndarray)):
            self.plot_data = pa.DataFrame()
            for ic,iplot in enumerate(plot_data):
                if isinstance(iplot,(tuple,list,np.ndarray)):
                    hold_series = pa.Series()
                    hold_series['x_data'] = iplot[0]
                    hold_series['y_data'] = iplot[1]
                    if len(iplot) == 4:
                        hold_series['xerr_data'] = iplot[2]
                        hold_series['yerr_data'] = iplot[3]
                    elif len(iplot) == 3:
                        hold_series['yerr_data'] = iplot[2]
                    self.plot_data['data'+str(ic+1)] = hold_series
                else:
                    self.plot_data['data'+str(ic+1)] = iplot
        else:
            raise IOError('invalid data type for plot_data passed into plotting class')
        self.FormatSetsOfFits()
        self.GetFuns()
        # self.AddNulls()

    def AppendData(self,plot_data):
        if isinstance(plot_data,pa.Series):
            if 'label' in list(plot_data.keys()):
                plot_data.name = plot_data['label']
            this_name = str(plot_data.name)
            if this_name in list(self.plot_data.keys()):
                key_next = 2
                while str(this_name)+'_'+str(key_next) in list(self.plot_data.keys()):
                    key_next += 1
                plot_data.name = this_name+'_'+str(key_next)
            this_name = str(plot_data.name)
            self.plot_data = self.plot_data.assign(**{this_name : plot_data})
        elif isinstance(plot_data,(dict,OrderedDict)):
            for ikey,ival in plot_data.items():
                this_name = str(ikey)
                if 'label' in list(ival.keys()):
                    this_name = ival['label']
                if this_name in list(self.plot_data.keys()):
                    key_next = 2
                    while str(this_name)+'_'+str(key_next) in list(self.plot_data.keys()):
                        key_next += 1
                    this_name = this_name+'_'+str(key_next)
                self.plot_data = self.plot_data.assign(**{this_name : ival})
        elif isinstance(plot_data,(tuple,list,np.ndarray)):
            for ic,iplot in enumerate(plot_data):
                if isinstance(iplot,(tuple,list,np.ndarray)):
                    hold_series = pa.Series()
                    hold_series['x_data'] = iplot[0]
                    hold_series['y_data'] = iplot[1]
                    if len(iplot) == 4:
                        hold_series['xerr_data'] = iplot[2]
                        hold_series['yerr_data'] = iplot[3]
                    elif len(iplot) == 3:
                        hold_series['yerr_data'] = iplot[2]
                    self.plot_data['data'+str(ic+1)] = hold_series
                else:
                    self.plot_data['data'+str(ic+1)] = iplot
        else:
            raise IOError('invalid data type for plot_data passed into plotting class')
        self.FormatSetsOfFits()
        self.GetFuns()
        self.AddNulls()

    def AddNulls(self):
        self.plot_data = self.plot_data.fillna(False)
        for icol,col_data in self.plot_data.items():
            if 'plot_this' not in col_data:
                self.plot_data[icol]['plot_this'] = True
            if 'type' not in col_data or col_data['type'] is False:
                if 'yerr_data' in col_data:
                    self.plot_data[icol]['type'] = 'error_bar'
                else:
                    self.plot_data[icol]['type'] = 'plot'
            if '_vary' in col_data['type'] and 'key_select' not in list(col_data.keys()):
                if 'y_data' in list(col_data.keys()) and isinstance( col_data['y_data'],pa.Series):
                    self.plot_data[icol]['key_select'] = list(self.plot_data[icol]['y_data'].keys())[0]
                elif 'boot_data' in list(col_data.keys()) and isinstance( col_data['boot_data'],pa.Series):
                    self.plot_data[icol]['key_select'] = list(self.plot_data[icol]['boot_data'].keys())[0]
                elif 'fit_class' in list(col_data.keys()) and isinstance( col_data['fit_class'],pa.Series):
                    self.plot_data[icol]['key_select'] = list(self.plot_data[icol]['fit_class'].keys())[0]
                elif 'x_data' in list(col_data.keys()) and isinstance( col_data['x_data'],pa.Series):
                    self.plot_data[icol]['key_select'] = list(self.plot_data[icol]['x_data'].keys())[0]
                if 'fit' in col_data['type']:
                    if isinstance(self.plot_data[icol]['key_select'],(list,tuple,np.ndarray)):
                        self.plot_data[icol]['key_select'] = tuple(self.plot_data[icol]['key_select'])
                    else:
                        self.plot_data[icol]['key_select'] = self.plot_data[icol]['key_select']
                else:
                    if isinstance(self.plot_data[icol]['key_select'],(list,tuple,np.ndarray)):
                        self.plot_data[icol]['key_select'] = (slice(None),)+tuple(self.plot_data[icol]['key_select'][1:])
                    else:
                        self.plot_data[icol]['key_select'] = slice(None)
            if 'Median' not in col_data:
                self.plot_data[icol]['Median'] = False
            if 'fmt_class' not in col_data:
                self.plot_data[icol]['fmt_class'] = None
            if 'phys_axies' not in col_data:
                self.plot_data[icol]['phys_axies'] = False
            if 'arrow_text_size' not in col_data:
                self.plot_data[icol]['arrow_text_size'] = 30
            if 'arrow_color' not in col_data:
                self.plot_data[icol]['arrow_color'] = 'black'
            if 'arrow_style' not in col_data:
                self.plot_data[icol]['arrow_style'] = 'fancy'
            if 'x_fun' not in col_data:
                self.plot_data[icol]['x_fun'] = None
            if 'label' not in col_data or col_data['label'] is False:
                self.plot_data[icol]['label'] = None
            if 'color' not in col_data or col_data['color'] is False:
                self.plot_data[icol]['color'] = None
            if 'symbol' not in col_data or col_data['symbol'] is False:
                self.plot_data[icol]['symbol'] = 'o'
            if 'shift' not in col_data or col_data['shift'] is False:
                self.plot_data[icol]['shift'] = None
            if 'shift_overlap' not in col_data or col_data['shift_overlap'] is False:
                self.plot_data[icol]['shift_overlap'] = 0.005
            if 'supress_legend' not in col_data:
                self.plot_data[icol]['supress_legend'] = False
            if 'hairline' not in col_data:
                self.plot_data[icol]['hairline'] = False
            if 'extrap_fade' not in col_data or col_data['extrap_fade'] == 0:
                self.plot_data[icol]['extrap_fade'] = 1.5
            if 'hair_alpha' not in col_data:
                self.plot_data[icol]['hair_alpha'] = (20,1)
            if 'FPName' not in col_data:
                self.plot_data[icol]['FPName'] = ''
            if 'FPUnits' not in col_data:
                self.plot_data[icol]['FPUnits'] = ''
            if 'x_range_min' not in col_data:
                self.plot_data[icol]['x_range_min'] = 0
            if 'x_scale' not in col_data:
                self.plot_data[icol]['x_scale'] = 1
            if 'x_range_max' not in col_data:
                self.plot_data[icol]['x_range_max'] = -1
            if 'x_increment' not in col_data:
                self.plot_data[icol]['x_increment'] = 1
            if 'alpha' not in col_data or col_data['alpha'] is False:
                self.plot_data[icol]['alpha'] = fillalpha
            if 'histtype' not in col_data or col_data['histtype'] is False:
                self.plot_data[icol]['histtype'] = 'stepfilled'
            if 'suppress_key_index' not in col_data or not isinstance(col_data['suppress_key_index'],bool):
                self.plot_data[icol]['suppress_key_index'] = False
            if 'scale' not in col_data or col_data['scale'] is False:
                self.plot_data[icol]['scale'] = 1.0
            if 'plot_err' not in col_data or not isinstance(col_data['plot_err'],str):
                self.plot_data[icol]['plot_err'] = ''
        if 'xlims' in self.plot_info:
            for icx,ix in enumerate(self.plot_info['xlims']):
                if ix is False:
                    self.plot_info['xlims'][icx] = None
        if 'ylims' in self.plot_info:
            for icy,iy in enumerate(self.plot_info['ylims']):
                if iy is False:
                    self.plot_info['ylims'][icy] = None
            # if 'fit' in col_data['type']:
            #     if 'fit_class' not in col_data or \
            #     not isinstance(col_data['fit_class'],ff.Fitting):
            #         raise IOError('plotting a fit requires fit class to be in fit_class')

    # def RemoveNulls(self):
    #     for icol,col_data in self.plot_data.iteritems():
    #         if 'label' not in col_data or col_data['label'] is None:
    #             self.plot_data[icol]['label'] = False
    #         if 'color' not in col_data or col_data['color'] is None:
    #             self.plot_data[icol]['color'] = False
    #         if 'symbol' not in col_data or col_data['symbol'] is None:
    #             self.plot_data[icol]['symbol'] = False
    #         if 'shift' not in col_data or col_data['shift'] is None:
    #             self.plot_data[icol]['shift'] = False
    #     if 'xlims' in self.plot_info:
    #         for icx,ix in enumerate(self.plot_info['xlims']):
    #             if ix is None:
    #                 self.plot_info['xlims'][icx] = False
    #     if 'ylims' in self.plot_info:
    #         for icy,iy in enumerate(self.plot_info['ylims']):
    #             if iy is None:
    #                 self.plot_info['ylims'][icy] = False

    def PrintData(self):
        print('INFO:')
        print(self.plot_info)
        print()
        print('DATA:')
        print(self.plot_data)

    def FormatSetsOfFits(self):
        for icol,col_data in self.plot_data.items():
            if 'fit' in col_data['type']:
                if 'fit_class' not in col_data:
                    raise EnvironmentError('fit_class not present even though fit was selected!')
                elif isinstance( col_data['fit_class'],sff.SetOfFitFuns):
                    col_data['fit_class'] = col_data['fit_class'].Fit_Stats_fmt['Fit']
                elif isinstance( col_data['fit_class'].iloc[0],sff.SetOfFitFuns):
                    col_data['fit_class'] = sff.PullSOFSeries(col_data['fit_class'])
        return True
    # def ScaleY(self,this_plot_data):
    #     if 'scale' in this_plot_data.keys():
    #         if this_plot_data['scale'] == False:
    #             this_plot_data['scale'] = 1.0
    #         else:
    #             if 'y_data' in this_plot_data.keys():
    #                 this_plot_data['y_data'] = np.array(this_plot_data['y_data'])*this_plot_data['scale']
    #             if 'yerr_data' in this_plot_data.keys():
    #                 this_yerrdata = np.array(this_yerrdata)*this_plot_data['scale']
    #     else:
    #         this_plot_data['scale'] = 1.0
    #     return this_plot_data

    def XRange_Chop(self,this_plot_data,xdim_list,*data):
        # if 'x_range_max' not in this_plot_data.keys():
        #     this_plot_data['x_range_max'] = -1
        # if 'x_range_min' not in this_plot_data.keys():
        #     this_plot_data['x_range_min'] = 0
        # if 'x_increment' not in this_plot_data.keys():
        #     this_plot_data['x_increment'] = 1
        x_max = this_plot_data['x_range_max']
        if x_max == -1:
            x_max = None
        x_min = this_plot_data['x_range_min']
        x_inc = this_plot_data['x_increment']
        if x_max is False:
            x_max = -1
        if x_min is False:
            x_min = 0
        # elif x_max <= x_min:
        #     x_max = x_min + 1
        if x_inc is False:
            x_inc = 1
        elif x_inc == 0:
            x_inc = 1
        data_out = []
        for idim,idata in enumerate(data):
            if idim in xdim_list and this_plot_data['x_scale'] != 0:
                this_scale = this_plot_data['x_scale']
            else:
                this_scale = 1
            if isinstance(idata,(list,tuple,np.ndarray)):
                data_out.append(np.array(idata)[x_min:x_max:x_inc]*this_scale)
            elif idata is not None:
                data_out.append(idata*this_scale)
            else:
                data_out.append(None)
        return data_out

    def PlotElement(self,this_plot_data,no_key_formatting=False):
        if not this_plot_data['plot_this']:
            return -1
        # this_plot_data = self.ScaleY(this_plot_data)
        is_vary = '_vary' in this_plot_data['type']
        if is_vary and 'key_select' not in list(this_plot_data.keys()):
            raise EnvironmentError('key_select needs to be set when calling _vary')
        med_leg = ''
        tupley_err,tuplex_err = False,False
        this_ydata,this_yerrdata = None,None
        this_xdata,this_xerrdata = None,None
        if 'Median' in this_plot_data and this_plot_data['Median']:
            if 'y_data_median' in this_plot_data:
                med_leg = 'Medy'
                this_ydata = this_plot_data['y_data_median']
            elif 'y_data' in this_plot_data:
                this_ydata = this_plot_data['y_data']
            if 'yerr_data_up' in this_plot_data and 'yerr_data_down' in this_plot_data:
                yerr_bool = this_plot_data['yerr_data_up'] is not None and \
                            this_plot_data['yerr_data_up'] is not False
                yerr_bool = yerr_bool and this_plot_data['yerr_data_down'] is not None and \
                            this_plot_data['yerr_data_down'] is not False
                if yerr_bool:
                    tupley_err = True
                    this_yerrdata = pa.Series(list(zip(this_plot_data['yerr_data_up'],this_plot_data['yerr_data_down'])),index=this_plot_data['yerr_data_up'].index)
                elif 'yerr_data' in this_plot_data:
                    this_yerrdata = this_plot_data['yerr_data']
            elif 'yerr_data' in this_plot_data:
                this_yerrdata = this_plot_data['yerr_data']
            if 'x_data_median' in this_plot_data:
                med_leg = 'Med'+med_leg.replace('Med','') + 'x'
                this_xdata = this_plot_data['x_data_median']
            elif 'x_data' in this_plot_data:
                this_xdata = this_plot_data['x_data']
            if 'xerr_data_up' in this_plot_data and 'xerr_data_down' in this_plot_data:
                xerr_bool = this_plot_data['xerr_data_up'] is not None and \
                            this_plot_data['xerr_data_up'] is not False
                xerr_bool = xerr_bool and this_plot_data['xerr_data_down'] is not None and \
                            this_plot_data['xerr_data_down'] is not False
                if xerr_bool:
                    tuplex_err = True
                    this_xerrdata = pa.Series(list(zip(this_plot_data['xerr_data_up'],this_plot_data['xerr_data_down'])),index=this_plot_data['xerr_data_up'].index)
            elif 'xerr_data' in this_plot_data:
                this_xerrdata = this_plot_data['xerr_data']
                xerr_bool = this_xerrdata is not None and this_xerrdata is not False
            else:
                xerr_bool = False

        else:
            if 'y_data' in this_plot_data:
                this_ydata = this_plot_data['y_data']
            if 'yerr_data' in this_plot_data:
                this_yerrdata = this_plot_data['yerr_data']
            if 'x_data' in this_plot_data:
                this_xdata = this_plot_data['x_data']
            if 'xerr_data' in this_plot_data:
                this_xerrdata = this_plot_data['xerr_data']
                xerr_bool = this_xerrdata is not None and this_xerrdata is not False
            else:
                xerr_bool = False
        if 'fit_class' in this_plot_data:
            this_fit = this_plot_data['fit_class']
        else:
            this_fit = False
        if 'boot_data' in this_plot_data:
            this_boot = this_plot_data['boot_data']
        else:
            this_boot = False


        if this_plot_data['fmt_class'] is not False and this_plot_data['fmt_class'] is not None:
            if 'key_select' in this_plot_data and this_plot_data['key_select'] is not False:
                if no_key_formatting:
                    this_key_select = tuple_wrap(this_plot_data['key_select'])
                else:
                    this_key_select = this_plot_data['fmt_class'].FormatKeySelect(tuple_wrap(this_plot_data['key_select']),
                                                                this_plot_data['phys_axies'])
            if this_boot is not None:
                this_boot = this_plot_data['fmt_class'].FormatSeriesKeys(this_boot,this_plot_data['phys_axies'])
            if this_fit is not None:
                this_fit = this_plot_data['fmt_class'].FormatSeriesKeys(this_fit,this_plot_data['phys_axies'])
            if this_ydata is not None:
                this_ydata = this_plot_data['fmt_class'].FormatSeriesKeys(this_ydata,this_plot_data['phys_axies'])
            if this_xdata is not None:
                this_xdata = this_plot_data['fmt_class'].FormatSeriesKeys(this_xdata,this_plot_data['phys_axies'])
            if this_xerrdata is not None:
                this_xerrdata = this_plot_data['fmt_class'].FormatSeriesKeys(this_xerrdata,this_plot_data['phys_axies'])
            if this_yerrdata is not None:
                this_yerrdata = this_plot_data['fmt_class'].FormatSeriesKeys(this_yerrdata,this_plot_data['phys_axies'])
        else:
            if 'key_select' in this_plot_data and this_plot_data['key_select'] is not False:
                this_key_select = tuple_wrap(this_plot_data['key_select'])
        if len(med_leg) > 0:
            med_leg += ' '
        output_dict = None
        if 'error_band' in this_plot_data['type']:
            # if len(this_xdata) == 0: return -1
            #TODO implement median for xdata too

            if is_vary:
                if not isinstance(this_ydata,pa.Series):
                    raise EnvironmentError(str(type(this_ydata))+'wrong data type for _vary')
                if isinstance(this_plot_data['key_select'],(list,tuple,np.ndarray)):
                    this_leg = []
                    for ikey in this_plot_data['key_select']:
                        if isinstance(ikey,str):
                            this_leg.append(ikey)
                    this_leg = ' '.join(this_leg)
                else:
                    this_leg = ' '
                output_dict = this_ydata.index
                if isinstance(this_xdata,pa.Series):
                    this_xdata = this_xdata.loc[this_key_select].values
                    hold_x = []
                    for ix in this_xdata:
                        if get_val_float_0(ix) in list(map(get_val_float_0,list(this_ydata.loc[this_key_select].index))):
                            hold_x.append(get_val_float_0(ix))
                    this_xdata = hold_x
                elif isinstance(this_xdata,str) and this_xdata == 'from_keys':
                    this_key_select = Fix_Key_Select(this_key_select,this_ydata)
                    try:
                        this_xdata = list(this_ydata.loc[this_key_select].index)
                    except Exception as err:
                        out_str = 'key_select is incorrect\n'
                        out_str += str(this_key_select) + '\n'
                        print(this_ydata.to_string())
                        raise Exception(out_str + str(err))
                    if isinstance(this_xdata[0],(list,tuple,np.ndarray)) and len(this_xdata[0]) > 1:
                        slice_loc = this_key_select.index(slice(None))
                        if isinstance(this_xdata[0],(list,tuple,np.ndarray)):
                            this_xdata = [ix[slice_loc] for ix in this_xdata]
                        this_xdata = list(map(get_val_float_0,this_xdata))
                    elif isinstance(this_xdata[0],(list,tuple,np.ndarray)):
                        this_xdata = [get_val_float_0(ix[0]) for ix in this_xdata]
                    elif isinstance(this_xdata,(list,tuple,np.ndarray)):
                        this_xdata = [get_val_float_0(ix) for ix in this_xdata]
                this_ydata = this_ydata.loc[this_key_select]
                this_yerrdata = this_yerrdata.loc[this_key_select]
            else:
                this_leg = ''
            if len(this_xdata) != len(this_ydata):
                print('Warning, xdata and ydata are different length, manually cutting')
                if len(this_xdata) > len(this_ydata):
                    this_xdata = this_xdata[:len(this_ydata)]
                elif len(this_xdata) < len(this_ydata):
                    this_ydata = this_ydata[:len(this_xdata)]
            if not xerr_bool:
                if len(this_xdata) == 0: return -1
                if isinstance(this_plot_data['shift'],str):
                    this_plot_data['shift'] = 0
                if tupley_err:
                    this_yerrdata = np.array([list(ival) for ival in this_yerrdata.values]).swapaxes(0,1)
                if not this_plot_data['suppress_key_index']:
                    this_leg = med_leg+str(LegendFmt(this_plot_data['label']+' '+this_leg))
                else:
                    this_leg = med_leg+str(LegendFmt(this_plot_data['label']))
                if this_plot_data['x_fun'] is not None and this_plot_data['x_fun'] is not False:
                    this_xdata = [this_plot_data['x_fun'](ix) for ix in this_xdata]
                try:
                    np.array(this_xdata)+float(this_plot_data['shift'])
                except Exception as err:
                    out_str = 'Shift cannot be added to this_xdata\n'
                    out_str += str(type(this_plot_data['shift'])) + ': ' + str(float(this_plot_data['shift'])) + '\n'
                    out_str += str(type(this_xdata )) + ': \n'
                    for ix in this_xdata:
                        out_str += '    '+str(ix) + '\n'
                    raise Exception(out_str +str(err))
                if this_plot_data['scale'] is not False and this_plot_data['scale'] is not None:
                    this_yerrdata = np.array(this_yerrdata)*this_plot_data['scale']
                    this_ydata = np.array(this_ydata)*this_plot_data['scale']
                if this_plot_data['supress_legend']: this_leg = None
                this_boot,this_ydata,this_xdata,this_xerrdata,this_yerrdata = self.XRange_Chop(this_plot_data,[2,3],this_boot,this_ydata,this_xdata,this_xerrdata,this_yerrdata)
                this_xdata = Shift_Duplicates(this_xdata,this_shiftper=float(this_plot_data['shift_overlap']))
                yup = np.array(this_ydata)+np.array(this_yerrdata)
                ydown = np.array(this_ydata)-np.array(this_yerrdata)
                self.plot_plane.plot(   np.array(this_xdata)+float(this_plot_data['shift']),
                                            this_ydata,
                                            label=this_leg,
                                            color=this_plot_data['color'])
                self.plot_plane.fill_between(   np.array(this_xdata)+float(this_plot_data['shift']),
                                            yup,
                                            ydown,
                                            edgecolor='none',
                                            alpha=this_plot_data['alpha'],
                                            color=this_plot_data['color'])


        if 'arrow' in this_plot_data['type']:
            location = (this_xdata[0],this_ydata[0])
            text_loc = (this_xdata[1],this_ydata[1])
            if this_plot_data['arrow_color'] == False:
                this_plot_data['arrow_color'] = 'black'
            if this_plot_data['arrow_style'] == False:
                this_plot_data['arrow_style'] = 'fancy'
            if this_plot_data['arrow_text_size'] == False:
                this_plot_data['arrow_text_size'] = 30
            if this_plot_data['supress_legend']:
                self.plot_plane.annotate('', xy=location, xytext=text_loc,
                    arrowprops={'arrowstyle': this_plot_data['arrow_style'],'color':this_plot_data['arrow_color']}, va='center',size=this_plot_data['arrow_text_size'])
            else:
                self.plot_plane.annotate(LegendFmt(this_plot_data['label']), xy=location, xytext=text_loc,
                    arrowprops={'arrowstyle': this_plot_data['arrow_style'],'color':this_plot_data['arrow_color']}, va='center',size=this_plot_data['arrow_text_size'])
        elif 'error_bar' in this_plot_data['type']:
            # if len(this_xdata) == 0: return -1
            #TODO implement median for xdata too

            if is_vary:
                if not isinstance(this_ydata,pa.Series):
                    raise EnvironmentError(str(type(this_ydata))+'wrong data type for _vary')
                if isinstance(this_plot_data['key_select'],(list,tuple,np.ndarray)):
                    this_leg = []
                    for ikey in this_plot_data['key_select']:
                        if isinstance(ikey,str):
                            this_leg.append(ikey)
                    this_leg = ' '.join(this_leg)
                else:
                    this_leg = ' '
                output_dict = this_ydata.index
                if isinstance(this_xdata,pa.Series):
                    this_xdata = this_xdata.loc[this_key_select].values
                    hold_x = []
                    for ix in this_xdata:
                        if get_val_float_0(ix) in list(map(get_val_float_0,list(this_ydata.loc[this_key_select].index))):
                            hold_x.append(get_val_float_0(ix))
                    this_xdata = hold_x
                elif isinstance(this_xdata,str) and this_xdata == 'from_keys':
                    this_key_select = Fix_Key_Select(this_key_select,this_ydata)
                    try:
                        this_xdata = list(this_ydata.loc[this_key_select].index)
                    except Exception as err:
                        this_key_select = list(this_ydata.index)[0]
                        this_key_select = [slice(None)] + list(this_key_select)[1:]
                        print('key not found, setting to first key',this_key_select)
                        this_key_select = tuple(this_key_select)
                        this_xdata = list(this_ydata.loc[this_key_select].index)
                        # out_str = 'key_select is incorrect\n'
                        # out_str += str(this_key_select) + '\n'
                        # print(this_ydata.to_string())
                        # raise Exception(out_str + str(err))
                    if isinstance(this_xdata[0],(list,tuple,np.ndarray)) and len(this_xdata[0]) > 1:
                        slice_loc = this_key_select.index(slice(None))
                        if isinstance(this_xdata[0],(list,tuple,np.ndarray)):
                            this_xdata = [ix[slice_loc] for ix in this_xdata]
                        this_xdata = list(map(get_val_float_0,this_xdata))
                    elif isinstance(this_xdata[0],(list,tuple,np.ndarray)):
                        this_xdata = [get_val_float_0(ix[0]) for ix in this_xdata]
                    elif isinstance(this_xdata,(list,tuple,np.ndarray)):
                        this_xdata = [get_val_float_0(ix) for ix in this_xdata]
                this_ydata = this_ydata.loc[this_key_select]
                this_yerrdata = this_yerrdata.loc[this_key_select]
            else:
                this_leg = ''
            if (not isinstance(this_xdata,str)) and len(this_xdata) != len(this_ydata):
                print('Warning, xdata and ydata are different length, manually cutting')
                if len(this_xdata) > len(this_ydata):
                    this_xdata = this_xdata[:len(this_ydata)]
                elif len(this_xdata) < len(this_ydata):
                    this_ydata = this_ydata[:len(this_xdata)]
            if xerr_bool:
                if '_vary' in this_plot_data['type']:
                    if isinstance(this_xdata,str) and this_xdata == 'from_keys':
                        raise EnvironmentError('from_keys for x_data cannot be used with xerr_data included')
                    if isinstance(this_xerrdata,pa.Series):
                        this_xerrdata = this_xerrdata[this_key_select].values
                # print this_plot_data['label']
                # print this_xdata
                # print float(this_plot_data['shift'])
                # print np.array(this_xdata)+float(this_plot_data['shift'])
                # print this_xerrdata
                if tuplex_err:
                    this_xerrdata = np.array([list(ival) for ival in this_xerrdata.values]).swapaxes(0,1)
                if tupley_err:
                    this_yerrdata = np.array([list(ival) for ival in this_yerrdata.values]).swapaxes(0,1)
                if not this_plot_data['suppress_key_index']:
                    this_leg = med_leg+str(LegendFmt(this_plot_data['label']+' '+this_leg))
                else:
                    this_leg = med_leg+str(LegendFmt(this_plot_data['label']))
                if this_plot_data['x_fun'] is not None and this_plot_data['x_fun'] is not False:
                    this_xdata = [this_plot_data['x_fun'](ix) for ix in this_xdata]
                if this_plot_data['scale'] is not False and this_plot_data['scale'] is not None:
                    this_yerrdata = np.array(this_yerrdata)*this_plot_data['scale']
                    this_ydata = np.array(this_ydata)*this_plot_data['scale']
                if this_plot_data['supress_legend']: this_leg = None
                this_boot,this_ydata,this_xdata,this_xerrdata,this_yerrdata = self.XRange_Chop(this_plot_data,[2,3],this_boot,this_ydata,this_xdata,this_xerrdata,this_yerrdata)
                if len(this_xdata) == 0: return -1
                if 'abs' in this_plot_data['plot_err'].lower():
                    self.plot_plane.scatter(np.array(this_xdata)+float(this_plot_data['shift']),
                                            this_yerrdata,
                                            label=this_leg,
                                            marker=this_plot_data['symbol'],
                                            color=this_plot_data['color'])
                elif 'rel' in this_plot_data['plot_err'].lower():
                    self.plot_plane.scatter(np.array(this_xdata)+float(this_plot_data['shift']),
                                            np.abs(this_yerrdata/this_ydata),
                                            label=this_leg,
                                            marker=this_plot_data['symbol'],
                                            color=this_plot_data['color'])
                else:
                    this_xdata = Shift_Duplicates(this_xdata,this_shiftper=float(this_plot_data['shift_overlap']))
                    self.plot_plane.errorbar(   np.array(this_xdata)+float(this_plot_data['shift']),
                                                this_ydata,
                                                yerr=this_yerrdata,
                                                xerr=this_xerrdata,
                                                label=this_leg,
                                                fmt=this_plot_data['symbol'],
                                                color=this_plot_data['color'])
            else:
                if len(this_xdata) == 0: return -1
                if isinstance(this_plot_data['shift'],str):
                    this_plot_data['shift'] = 0
                if tupley_err:
                    this_yerrdata = np.array([list(ival) for ival in this_yerrdata.values]).swapaxes(0,1)
                if not this_plot_data['suppress_key_index']:
                    this_leg = med_leg+str(LegendFmt(this_plot_data['label']+' '+this_leg))
                else:
                    this_leg = med_leg+str(LegendFmt(this_plot_data['label']))
                if isinstance(this_xdata,str) and this_xdata == 'from_keys':
                    this_xdata = np.array(list(map(get_val_float_0,this_ydata.index)))
                    if isinstance(this_xdata[0],(list,tuple,np.ndarray)) and len(this_xdata[0]) > 1:
                        this_xdata = this_xdata[:,0]
                try:
                    np.array(this_xdata)+float(this_plot_data['shift'])
                except Exception as err:
                    out_str = 'Shift cannot be added to this_xdata\n'
                    out_str += str(type(this_plot_data['shift'])) + ': ' + str(float(this_plot_data['shift'])) + '\n'
                    out_str += str(type(this_xdata )) + ': \n'
                    for ix in this_xdata:
                        out_str += '    '+str(ix) + '\n'
                    raise Exception(out_str +str(err))
                if this_plot_data['x_fun'] is not None and this_plot_data['x_fun'] is not False:
                    this_xdata = [this_plot_data['x_fun'](ix) for ix in this_xdata]
                if this_plot_data['scale'] is not False and this_plot_data['scale'] is not None:
                    this_yerrdata = np.array(this_yerrdata)*this_plot_data['scale']
                    this_ydata = np.array(this_ydata)*this_plot_data['scale']
                if this_plot_data['supress_legend']: this_leg = None
                this_boot,this_ydata,this_xdata,this_xerrdata,this_yerrdata = self.XRange_Chop(this_plot_data,[2,3],this_boot,this_ydata,this_xdata,this_xerrdata,this_yerrdata)
                if 'abs' in this_plot_data['plot_err'].lower():
                    self.plot_plane.scatter(np.array(this_xdata)+float(this_plot_data['shift']),
                                            this_yerrdata,
                                            label=this_leg,
                                            marker=this_plot_data['symbol'],
                                            color=this_plot_data['color'])
                elif 'rel' in this_plot_data['plot_err'].lower():
                    self.plot_plane.scatter(np.array(this_xdata)+float(this_plot_data['shift']),
                                            np.abs(this_yerrdata/this_ydata),
                                            label=this_leg,
                                            marker=this_plot_data['symbol'],
                                            color=this_plot_data['color'])
                else:
                    this_xdata = Shift_Duplicates(this_xdata,this_shiftper=float(this_plot_data['shift_overlap']))
                    if len(np.array(this_xdata)) != len(this_ydata): return -1
                    self.plot_plane.errorbar(   np.array(this_xdata)+float(this_plot_data['shift']),
                                                this_ydata,
                                                yerr=this_yerrdata,
                                                label=this_leg,
                                                fmt=this_plot_data['symbol'],
                                                color=this_plot_data['color'])
        elif 'plot' in this_plot_data['type']:
            if len(this_xdata) == 0: return -1
            if is_vary:
                if not isinstance(this_ydata,pa.Series):
                    raise EnvironmentError('wrong data type for _vary: ' + str(type(this_ydata)))
                if isinstance(this_plot_data['key_select'],(list,tuple,np.ndarray)):
                    this_leg = []
                    for ikey in this_plot_data['key_select']:
                        if isinstance(ikey,str):
                            this_leg.append(ikey)
                    this_leg = ' '.join(this_leg)
                else:
                    this_leg = ' '
                if isinstance(this_xdata,pa.Series):
                    this_xdata = this_xdata[this_key_select].values
                elif isinstance(this_xdata,str) and this_xdata == 'from_keys':
                    this_key_select = Fix_Key_Select(this_key_select,this_ydata)
                    this_xdata = list(this_ydata.loc[this_key_select].index)
                    slice_loc = this_key_select.index(slice(None))
                    if isinstance(this_xdata[0],(list,tuple,np.ndarray)):
                        this_xdata = [ix[slice_loc] for ix in this_xdata]
                    this_xdata = list(map(get_val_float_0,this_xdata))
                this_ydata = this_ydata.loc[this_key_select].values
            else:
                if isinstance(this_xdata,pa.Series):
                    this_xdata = this_xdata[this_key_select].values
                elif isinstance(this_xdata,str) and this_xdata == 'from_keys':
                    this_xdata = list(this_ydata.index)
                    if isinstance(this_xdata[0],(tuple,list,np.ndarray)):
                        this_xdata = [ix[0] for ix in this_xdata]
                    this_xdata = list(map(get_val_float_0,this_xdata))
                this_leg = ''
            if not this_plot_data['suppress_key_index']:
                this_leg = med_leg+str(LegendFmt(this_plot_data['label']+' '+this_leg))
            else:
                this_leg = med_leg+str(LegendFmt(this_plot_data['label']))
            if this_plot_data['scale'] is not False and this_plot_data['scale'] is not None:
                this_ydata = np.array(this_ydata)*this_plot_data['scale']
            if this_plot_data['supress_legend']: this_leg = None
            this_boot,this_ydata,this_xdata,this_xerrdata,this_yerrdata = self.XRange_Chop(this_plot_data,[2,3],this_boot,this_ydata,this_xdata,this_xerrdata,this_yerrdata)
            self.plot_plane.plot(   np.array(this_xdata)+float(this_plot_data['shift']),
                                        this_ydata,
                                        label=this_leg,
                                        color=this_plot_data['color'])
        elif 'histogram' in this_plot_data['type']:
            if 'boot_data' in list(this_plot_data.keys()):
                if is_vary:
                    if not isinstance(this_boot,pa.Series):
                        raise EnvironmentError('wrong data type for _vary')
                    if isinstance(this_plot_data['key_select'],(list,tuple,np.ndarray)):
                        this_leg = []
                        for ikey in this_plot_data['key_select']:
                            if isinstance(ikey,str):
                                this_leg.append(ikey)
                        this_leg = ' '.join(this_leg)
                    else:
                        this_leg = ' '
                    # if hasattr(this_plot_data['fmt_class'],'FormatSeriesKeys'):
                    #     this_boot = this_plot_data['fmt_class'].FormatSeriesKeys(this_boot,
                    #                                                                 this_plot_data['phys_axies'])
                    this_boot = this_boot[this_key_select]
                else:
                    this_leg = ''
                # if isinstance(this_plot_data,pa.DataFrame):
                #     this_plot_data = this_plot_data.iloc[0,0]
                # if isinstance(this_plot_data,pa.Series):
                #     this_plot_data = this_plot_data.iloc[0]
                # if isinstance(this_plot_data,pa.Series):
                #     this_plot_data = this_plot_data.iloc[0]
                if not isinstance(this_boot,BootStrap):
                    print(type(this_plot_data), 'is not bootstrap class')
                    print(this_boot)
                    return -1
                if not this_plot_data['suppress_key_index']:
                    this_leg = med_leg+str(LegendFmt(this_plot_data['label']+' '+this_leg))
                else:
                    this_leg = med_leg+str(LegendFmt(this_plot_data['label']))
                if this_plot_data['scale'] is not False and this_plot_data['scale'] is not None:
                    this_boot = np.array(this_boot)*this_plot_data['scale']
                if this_plot_data['supress_legend']: this_leg = None
                this_boot,this_ydata,this_xdata,this_xerrdata,this_yerrdata = self.XRange_Chop(this_plot_data,[2,3],this_boot,this_ydata,this_xdata,this_xerrdata,this_yerrdata)
                self.plot_plane = this_boot.BootHistogram(self.plot_plane,
                                    histtype=this_plot_data['histtype'],
                                    label=this_leg,
                                    color=this_plot_data['color'],
                                    alpha=this_plot_data['alpha'],
                                    Median=this_plot_data['Median'])
            elif len(this_xdata) == 0: return -1
            else:
                if is_vary:
                    if not isinstance(this_xdata,pa.Series):
                        raise EnvironmentError('wrong data type for _vary')
                    if isinstance(this_plot_data['key_select'],(list,tuple,np.ndarray)):
                        this_leg = []
                        for ikey in this_plot_data['key_select']:
                            if isinstance(ikey,str):
                                this_leg.append(ikey)
                        this_leg = ' '.join(this_leg)
                    else:
                        this_leg = ' '
                    this_xdata = this_xdata[this_key_select].values
                else:
                    this_leg = ''
                if not this_plot_data['suppress_key_index']:
                    this_leg = med_leg+str(LegendFmt(this_plot_data['label']+' '+this_leg))
                else:
                    this_leg = med_leg+str(LegendFmt(this_plot_data['label']))
                if this_plot_data['x_fun'] is not None and this_plot_data['x_fun'] is not False:
                    this_xdata = [this_plot_data['x_fun'](ix) for ix in this_xdata]
                if this_plot_data['scale'] is not False and this_plot_data['scale'] is not None:
                    this_xdata = np.array(this_xdata)*this_plot_data['scale']
                if this_plot_data['supress_legend']: this_leg = None
                this_boot,this_ydata,this_xdata,this_xerrdata,this_yerrdata = self.XRange_Chop(this_plot_data,[2,3],this_boot,this_ydata,this_xdata,this_xerrdata,this_yerrdata)
                self.plot_plane.hist(   np.hstack(this_xdata),
                                        histtype=this_plot_data['histtype'],
                                        label=this_leg,
                                        color=this_plot_data['color'],
                                        alpha=this_plot_data['alpha'])
        elif 'scatter' in this_plot_data['type']:
            if len(this_xdata) == 0: return -1
            if is_vary:
                if not isinstance(this_ydata,pa.Series):
                    raise EnvironmentError('wrong data type for _vary')
                if isinstance(this_plot_data['key_select'],(list,tuple,np.ndarray)):
                    this_leg = []
                    for ikey in this_plot_data['key_select']:
                        if isinstance(ikey,str):
                            this_leg.append(ikey)
                    this_leg = ' '.join(this_leg)
                else:
                    this_leg = ' '

                if isinstance(this_xdata,pa.Series):
                    this_xdata = this_xdata[this_key_select].values
                elif isinstance(this_xdata,str) and this_xdata == 'from_keys':
                    this_key_select = Fix_Key_Select(this_key_select,this_ydata)
                    this_xdata = list(this_ydata.loc[this_key_select].index)
                    slice_loc = this_key_select.index(slice(None))
                    if isinstance(this_xdata[0],(list,tuple,np.ndarray)):
                        this_xdata = [ix[slice_loc] for ix in this_xdata]
                    this_xdata = list(map(get_val_float_0,this_xdata))
                this_ydata = this_ydata.loc[this_key_select].values
            else:
                this_leg = ''
            if not this_plot_data['suppress_key_index']:
                this_leg = med_leg+str(LegendFmt(this_plot_data['label']+' '+this_leg))
            else:
                this_leg = med_leg+str(LegendFmt(this_plot_data['label']))
            if this_plot_data['x_fun'] is not None and this_plot_data['x_fun'] is not False:
                this_xdata = [this_plot_data['x_fun'](ix) for ix in this_xdata]
            if this_plot_data['scale'] is not False and this_plot_data['scale'] is not None:
                this_ydata = np.array(this_ydata)*this_plot_data['scale']
            if this_plot_data['supress_legend']: this_leg = None
            this_boot,this_ydata,this_xdata,this_xerrdata,this_yerrdata = self.XRange_Chop(this_plot_data,[2,3],this_boot,this_ydata,this_xdata,this_xerrdata,this_yerrdata)
            self.plot_plane.scatter(   np.array(this_xdata)+float(this_plot_data['shift']),
                                        this_ydata,
                                        label=this_leg,
                                        marker=this_plot_data['symbol'],
                                        color=this_plot_data['color'])
        elif 'vline' in this_plot_data['type']:
            if is_vary:
                if not isinstance(this_xdata,pa.Series):
                    raise EnvironmentError('wrong data type for _vary')
                if isinstance(this_plot_data['key_select'],(list,tuple,np.ndarray)):
                    this_leg = []
                    for ikey in this_plot_data['key_select']:
                        if isinstance(ikey,str):
                            this_leg.append(ikey)
                    this_leg = ' '.join(this_leg)
                else:
                    this_leg = ' '
                this_xdata = this_xdata.loc[this_key_select]
            elif isinstance(this_xdata,(list,tuple,np.ndarray)):
                this_xdata = this_xdata[0]
                this_leg = ''
            else:
                this_leg = ''
            if isinstance(  this_xdata,(list,tuple,np.ndarray)) and \
                            len(this_xdata) == 0: return -1
            if 'ACTUAL_VALUE_' in this_plot_data['label']:
                prefix,value = this_plot_data['label'].split('ACTUAL_VALUE_')
                value += '='
                if 'xerr_data' in this_plot_data:
                    if is_vary:
                        err = this_xerrdata[this_key_select]
                    elif isinstance(this_xerrdata,(list,tuple,np.ndarray)):
                        err = this_xerrdata[0]
                    else:
                        err = this_xerrdata
                    if isinstance(err,(list,tuple,np.ndarray)):
                        err = np.mean(err)
                    this_leg += ' '+prefix+r' $'+value+MakeValAndErr(this_xdata,err,latex=True)+r'$'
                else:
                    this_leg += ' '+prefix+r' $'+value+MakeVal(this_xdata,Dec=3)+r'$'
            else:
                this_leg = this_plot_data['label']+' '+this_leg
            if not this_plot_data['suppress_key_index']:
                this_leg = med_leg+str(LegendFmt(this_plot_data['label']+' '+this_leg))
            else:
                this_leg = med_leg+str(LegendFmt(this_plot_data['label']))
            if this_plot_data['scale'] is not False and this_plot_data['scale'] is not None:
                this_xdata = np.array(this_xdata)*this_plot_data['scale']
            if this_plot_data['supress_legend']: this_leg = None
            this_xdata = self.XRange_Chop(this_plot_data,[0],this_xdata)
            self.plot_plane.axvline(    np.array(this_xdata)+float(this_plot_data['shift']),
                                        label=this_leg,
                                        color=this_plot_data['color'])
            if 'xerr_data' in this_plot_data:
                if is_vary:
                    mid = this_xdata[this_key_select].values,
                    err = this_xerrdata[this_key_select]
                elif isinstance(this_xerrdata,(list,tuple,np.ndarray)):
                    mid,err = this_xdata[0],this_xerrdata[0]
                else:
                    mid,err = this_xdata,this_xerrdata
                if isinstance(err,(list,tuple,np.ndarray)):
                    up,down = mid+err[0],mid-err[1]
                else:
                    up,down = mid+err,mid-err
                if this_plot_data['scale'] is not False and this_plot_data['scale'] is not None:
                    down = np.array(down)*this_plot_data['scale']
                    up = np.array(up)*this_plot_data['scale']
                up,down = self.XRange_Chop(this_plot_data,[0,1],up,down)
                self.plot_plane.axvspan(    down+float(this_plot_data['shift']),
                                            up+float(this_plot_data['shift']),
                                            alpha=this_plot_data['alpha'],
                                            color=this_plot_data['color'])

        elif 'hline' in this_plot_data['type']:
            if is_vary:
                if not isinstance(this_ydata,pa.Series):
                    raise EnvironmentError(str(type(this_ydata))+' wrong data type for _vary')
                if isinstance(this_plot_data['key_select'],(list,tuple,np.ndarray)):
                    this_leg = []
                    for ikey in this_plot_data['key_select']:
                        if isinstance(ikey,str):
                            this_leg.append(ikey)
                    this_leg = ' '.join(this_leg)
                else:
                    this_leg = ' '
                this_key_select = Fix_Key_Select(this_key_select,this_ydata)
                this_ydata = this_ydata.loc[this_key_select]
                if isinstance(this_ydata,pa.Series):
                    this_ydata = this_ydata.values[0]
            elif isinstance(this_ydata,(list,tuple,np.ndarray)):
                this_ydata = this_ydata[0]
                this_leg = ''
            else:
                this_leg = ''
            if 'ACTUAL_VALUE_' in this_plot_data['label']:
                prefix,value = this_plot_data['label'].split('ACTUAL_VALUE_')
                value += '='
                if 'yerr_data' in this_plot_data:
                    if is_vary:
                        try:
                            err = this_yerrdata[this_key_select]
                        except:
                            err = this_yerrdata[this_key_select[0]]
                    elif isinstance(this_yerrdata,(list,tuple,np.ndarray)):
                        err = this_yerrdata[0]
                    else:
                        err = this_yerrdata
                    if isinstance(err,(list,tuple,np.ndarray)):
                        err = np.mean(err)
                    this_leg += ' '+prefix+r' $'+value+MakeValAndErr(this_ydata,err,latex=True)+r'$'
                else:
                    this_leg += ' '+prefix+r' $'+value+MakeVal(this_ydata,Dec=3)+r'$'
            else:
                this_leg = this_plot_data['label']+' '+this_leg
            if not this_plot_data['suppress_key_index']:
                this_leg = med_leg+str(LegendFmt(this_plot_data['label']+' '+this_leg))
            else:
                this_leg = med_leg+str(LegendFmt(this_plot_data['label']))
            if this_plot_data['scale'] is not False and this_plot_data['scale'] is not None:
                this_ydata = np.array(this_ydata)*this_plot_data['scale']
            if this_plot_data['supress_legend']: this_leg = None
            # this_ydata = self.XRange_Chop(this_plot_data,this_ydata)
            self.plot_plane.axhline(    this_ydata,
                                        label=this_leg,
                                        color=this_plot_data['color'])
            if 'yerr_data' in this_plot_data:
                if is_vary:
                    try:
                        err = this_yerrdata[this_key_select]
                    except:
                        err = this_yerrdata[this_key_select[0]]
                elif isinstance(this_yerrdata,(list,tuple,np.ndarray)):
                    err = this_yerrdata[0]
                else:
                    err = this_yerrdata
                if isinstance(err,(list,tuple,np.ndarray)):
                    up,down = this_ydata+err[0],this_ydata-err[1]
                elif isinstance(this_ydata,(list,tuple,np.ndarray)):
                    up,down = this_ydata[0]+err,this_ydata[0]-err
                else:
                    up,down = this_ydata+err,this_ydata-err
                if this_plot_data['scale'] is not False and this_plot_data['scale'] is not None:
                    down = np.array(down)*this_plot_data['scale']
                    up = np.array(up)*this_plot_data['scale']
                # up,down = self.XRange_Chop(this_plot_data,up,down)
                self.plot_plane.axhspan(    down,
                                            up,
                                            alpha=this_plot_data['alpha'],
                                            color=this_plot_data['color'])
        elif 'effm_fit' in this_plot_data['type']:
            if 'xdatarange' not in this_plot_data: this_plot_data['xdatarange'] = 'Data'
            if 'thislineres' not in this_plot_data: this_plot_data['thislineres'] = 100
            if 'Phys' not in this_plot_data: this_plot_data['Phys'] = True
            if this_plot_data['label'] is not None:
                this_fit.name = this_plot_data['label']
            if is_vary:
                if not isinstance(this_fit,pa.Series):
                    raise EnvironmentError('wrong data type for _vary')
                # if this_key_select not in list(this_fit.index):
                #     this_key_select = this_fit.index[0]
                # this_key_select = Fix_Key_Select(this_key_select,this_fit)

                try:
                    this_fit = this_fit[this_key_select]
                except:
                    this_key_select,this_fit = Fix_Se_Keys(this_key_select,this_fit)
                    this_fit = this_fit[this_key_select]
            if this_plot_data['fmt_class'] is not False and this_plot_data['fmt_class'] is not None:
                def this_scale(val):
                    return this_plot_data['fmt_class'].FmtTime(val,this_plot_data['phys_axies'])
            else:
                this_scale = lambda x : x

            self.plot_plane = this_fit.PlotEffMass(
                                xdatarange=this_plot_data['xdatarange'],
                                thislineres=this_plot_data['thislineres'],
                                color=this_plot_data['color'],
                                shift=this_plot_data['shift'],
                                Phys=this_plot_data['Phys'],
                                y_scale=this_plot_data['scale'],
                                x_scale=this_plot_data['x_scale'],
                                plot_plane = self.plot_plane)
        elif 'log_fit' in this_plot_data['type']:
            if 'xdatarange' not in this_plot_data: this_plot_data['xdatarange'] = 'Data'
            if 'thislineres' not in this_plot_data: this_plot_data['thislineres'] = 100
            if this_plot_data['label'] is not None:
                this_fit.name = this_plot_data['label']
            if is_vary:
                if not isinstance(this_fit,pa.Series):
                    raise EnvironmentError('wrong data type for _vary')
                this_fit = this_fit[this_key_select]
            self.plot_plane = this_fit.PlotLog(
                                xdatarange=this_plot_data['xdatarange'],
                                thislineres=this_plot_data['thislineres'],
                                scale=this_plot_data['scale'],
                                color=this_plot_data['color'],
                                shift=this_plot_data['shift'],
                                plot_plane = self.plot_plane)
        elif 'fit' in this_plot_data['type']:
            if 'otherXvals' not in this_plot_data: this_plot_data['otherXvals'] = []
            if 'xaxis' not in this_plot_data: this_plot_data['xaxis'] = 0
            if 'xdatarange' not in this_plot_data: this_plot_data['xdatarange'] = 'Data'
            if 'thislineres' not in this_plot_data: this_plot_data['thislineres'] = 100
            if 'flowphys' not in this_plot_data: this_plot_data['flowphys'] = False
            if 'ShowPar' not in this_plot_data: this_plot_data['ShowPar'] = 'None'
            if 'ShowEval' not in this_plot_data: this_plot_data['ShowEval'] = None
            if this_plot_data['label'] is not None:
                this_fit.name = this_plot_data['label']
            if is_vary:
                if not isinstance(this_fit,pa.Series):
                    raise EnvironmentError('wrong data type for _vary')
                this_key = this_key_select
                if any([type(jkey) == slice for jkey in this_key]):
                    this_key = this_fit.index[0]
                if this_key in this_fit.index:
                    this_fit = this_fit[this_key]
                    if isinstance(this_fit,pa.Series):
                        this_fit = this_fit.iloc[0]
                else:
                    if 'First' != this_key:
                        print('Warning, ',this_key, 'not in fit plot keys')
                    this_fit = this_fit.iloc[0]
            if this_plot_data['fmt_class'] is not False and this_plot_data['fmt_class'] is not None:
                def this_scale(val):
                    return this_plot_data['fmt_class'].FmtTime(val,this_plot_data['phys_axies'])
            else:
                this_scale = lambda x : x
            this_fit.name = med_leg+str(LegendFmt(this_plot_data['label']))
            if this_plot_data['ShowEval'] is not None and this_plot_data['ShowEval'] is not False:
                if not isinstance(this_plot_data['xaxis'],int):
                    this_plot_data['xaxis'] = 0
                if isinstance(this_plot_data['ShowEval'],(list,tuple,np.ndarray)):
                    this_plot_data['otherXvals'] = [x for ix,x in enumerate(this_plot_data['ShowEval']) if ix != this_plot_data['xaxis']]
            if isinstance(this_plot_data['hair_alpha'],(tuple,list,np.ndarray)):
                this_fit.hairline_dropoff = this_plot_data['hair_alpha'][0]
                this_fit.hairline_alpha = this_plot_data['hair_alpha'][1]
            else:
                this_fit.hairline_dropoff = this_plot_data['hair_alpha']
                this_fit.hairline_alpha = 1.0

            if isinstance(this_plot_data['hair_alpha'],(tuple,list,np.ndarray)):
                this_fit.hairline_dropoff = this_plot_data['hair_alpha'][0]
                this_fit.hairline_alpha = this_plot_data['hair_alpha'][1]
            else:
                this_fit.hairline_dropoff = this_plot_data['hair_alpha']
                this_fit.hairline_alpha = 1.0
            this_fit.FPName = this_plot_data['FPName']
            this_fit.FPUnits = this_plot_data['FPUnits']
            self.plot_plane = this_fit.PlotFunction(
                                otherXvals=this_plot_data['otherXvals'],
                                xaxis=this_plot_data['xaxis'],
                                xdatarange=this_plot_data['xdatarange'],
                                thislineres=this_plot_data['thislineres'],
                                color=this_plot_data['color'],
                                shift=this_plot_data['shift'],
                                flowphys=this_plot_data['flowphys'],
                                ShowPar=this_plot_data['ShowPar'],
                                ShowEval=this_plot_data['ShowEval'],
                                y_scale=this_plot_data['scale'],
                                x_scale=this_plot_data['x_scale'],
                                suppress_key=this_plot_data['suppress_key_index'],
                                x_fun=this_plot_data['x_fun'],
                                supress_legend=this_plot_data['supress_legend'],
                                hairline=this_plot_data['hairline'],
                                extrap_fade=this_plot_data['extrap_fade'],
                                plot_plane = self.plot_plane)
        else:
            print('Warning, plpython kill programot type',this_plot_data['type'],'not implemented, skipping data')
        return output_dict

    def Animate(self,this_key=0,a_index=0,multi_index=0):
        self.AddNulls()
        self.plot_plane.cla()
        if isinstance(this_key,int):
            this_pd = copy(self.plot_data.iloc[:,0])
            if this_pd['label'] is None:
                this_pd['label'] = self.plot_data.columns[0]
        else:
            this_pd = copy(self.plot_data[this_key])
            if this_pd['label'] is None:
                this_pd['label'] = this_key
        if this_pd['color'] is None:
            this_pd['color'] = colourset8[0]
        if this_pd['symbol'] is None:
            this_pd['symbol'] = markerset[0]
        if this_pd['shift'] is None:
            this_pd['shift'] = shiftset[0]
        if '_vary' not in this_pd['type']:
            raise EnvironmentError('animating something must be _vary type of plot')
        this_pd['color'] = fmt_Qt5_col(this_pd['color'])
        all_keys = self.PlotElement(this_pd)
        if not isinstance(all_keys,pa.MultiIndex):
            raise EnvironmentError('error ploting for animate')
        if this_pd['key_select'] == 'First':
            raise EnvironmentError('Animation not working for First atm')
        if len(this_pd['key_select']) < a_index:
            print('a_index',a_index)
            print('keys',this_pd['key_select'])
            raise IOError('dimensions is smaller than index to animate over')
        if isinstance(this_pd['key_select'][a_index],slice):
            print(a_index)
            raise IOError('cannot animate over dimension that is being plotted')
        if isinstance(a_index,str):
            a_index = list(all_keys.names).index(a_index)
        animate_keys = list(OrderedDict.fromkeys(all_keys.get_level_values(a_index)))
        # for ia_key in list(OrderedDict.fromkeys(animate_keys)):
        this_pd['key_select'] = list(this_pd['key_select'])
        this_pd['key_select'][a_index] = animate_keys[multi_index]
        this_pd['key_select'] = tuple(this_pd['key_select'])
        self.plot_plane.cla()
        self.PlotElement(this_pd,no_key_formatting=True)
        self.Configure_Plane()
        this_dir = '/'.join(self.plot_info['save_file'].split('/')[:-1])
        mkdir_p(this_dir)
        self.plot_info['save_file'] = self.plot_info['save_file'].replace('.pdf','')+'.pdf'
        self.figure.savefig(self.plot_info['save_file'])
        self.WriteData()
        self.figure.show()
        if multi_index < len(animate_keys)-1:
            return multi_index+1
        else:
            return -1

    def Configure_Plane(self):
        self.plot_plane.legend( loc=self.plot_info['leg_loc'],ncol=self.plot_info['leg_ncol'],
                                fontsize=self.plot_info['leg_fontsize'],framealpha=self.plot_info['leg_alpha'])
        if 'title' in self.plot_info:
            if 'title_dict' in self.plot_info:
                if 'title_pad' in self.plot_info:
                    self.plot_plane.set_title(  self.plot_info['title'],
                                                fontdict=self.plot_info['title_dict'],
                                                pad=self.plot_info['title_pad'])
                else:
                    self.plot_plane.set_title(  self.plot_info['title'],
                                                fontdict=self.plot_info['title_dict']
                                                ,pad=tpad)
            else:
                if 'title_pad' in self.plot_info:
                    self.plot_plane.set_title(self.plot_info['title'],pad=self.plot_info['title_pad'])
                else:
                    self.plot_plane.set_title(self.plot_info['title'],{'fontsize':tsize})
        if 'xlabel' in self.plot_info:
            if 'xlabel_dict' in self.plot_info:
                self.plot_plane.set_xlabel(self.plot_info['xlabel'],self.plot_info['xlabel_dict'])
            else:
                self.plot_plane.set_xlabel(self.plot_info['xlabel'],{'fontsize':xsize})
        if 'ylabel' in self.plot_info:
            if 'ylabel_dict' in self.plot_info:
                self.plot_plane.set_ylabel(self.plot_info['ylabel'],self.plot_info['ylabel_dict'])
            else:
                self.plot_plane.set_ylabel(self.plot_info['ylabel'],{'fontsize':ysize})
        if 'x_plot_scale' in self.plot_info and self.plot_info['x_plot_scale'].lower() != 'none':
            self.plot_plane.set_xscale(self.plot_info['x_plot_scale'])
        if 'y_plot_scale' in self.plot_info and self.plot_info['y_plot_scale'].lower() != 'none':
            self.plot_plane.set_yscale(self.plot_info['y_plot_scale'])
        if 'xlims' in self.plot_info:
            self.plot_info['xlims'] = list(self.plot_info['xlims'])
            if self.plot_info['xlims'][0] == None:
                self.plot_info['xlims'][0] = self.plot_plane.get_xlim()[0]
            if self.plot_info['xlims'][1] == None:
                self.plot_info['xlims'][1] = self.plot_plane.get_xlim()[1]
            self.plot_plane.set_xlim(self.plot_info['xlims'])
        if 'ylims' in self.plot_info:
            self.plot_info['ylims'] = list(self.plot_info['ylims'])
            if self.plot_info['ylims'][0] == None:
                self.plot_info['ylims'][0] = self.plot_plane.get_ylim()[0]
            if self.plot_info['ylims'][1] == None:
                self.plot_info['ylims'][1] = self.plot_plane.get_ylim()[1]
            self.plot_plane.set_ylim(self.plot_info['ylims'])

        if 'do_xTicks' in self.plot_info and self.plot_info['do_xTicks']:
            self.plot_plane.set_xticks(np.arange(self.plot_info['xTick_min'],self.plot_info['xTick_max'],step=self.plot_info['xTick_inc']))
        if 'do_yTicks' in self.plot_info and self.plot_info['do_yTicks']:
            self.plot_plane.set_yticks(np.arange(self.plot_info['yTick_min'],self.plot_info['yTick_max'],step=self.plot_info['yTick_inc']))



    def PlotAll(self):
        self.AddNulls()
        self.plot_plane.cla()
        errlist = []
        for icplot,(ikey,iplot) in enumerate(self.plot_data.items()):
            if iplot['color'] is None:
                try:
                    holdcol = next(self.colcyc)
                except StopIteration:
                    self.colcyc = iter(colourset8)
                    holdcol = next(self.colcyc)
                iplot['color'] = holdcol
            elif 'previous' == iplot['color']:
                iplot['color'] = holdcol
            else:
                self.plot_data[ikey]['color'] = fmt_Qt5_col(iplot['color'])
                holdcol = self.plot_data[ikey]['color']

            if iplot['symbol'] is None:
                try:
                    holdsym = next(self.symcyc)
                except StopIteration:
                    self.symcyc = iter(markerset)
                    holdsym = next(self.symcyc)
                iplot['symbol'] = holdsym
            elif 'previous' in iplot['symbol']:
                iplot['symbol'] = holdsym
            else:
                holdsym = iplot['symbol']

            if iplot['shift'] is None:
                try:
                    holdshift = next(self.shiftcyc)
                except StopIteration:
                    self.shiftcyc = iter(shiftset)
                    holdshift = next(self.shiftcyc)
                iplot['shift'] = holdshift
            elif 'previous' in str(iplot['shift']) :
                iplot['shift'] = holdshift
            else:
                holdshift = iplot['shift']

            if iplot['label'] is None:
                iplot['label'] = ikey
            this_out = self.PlotElement(iplot)
            if not isinstance(this_out,int):
                this_out = 1
            errlist.append(this_out)
        if all(np.array(errlist) == -1):
            print('Nothing to plot')
            return
        self.Configure_Plane()
        if 'x_zero_line' in self.plot_info and self.plot_info['x_zero_line']:
            if 'x_zero_line_val' in self.plot_info:
                self.PlotLine('x',self.plot_info['x_zero_line_val'])
            else:
                self.PlotLine('x',0)
        if 'y_zero_line' in self.plot_info and self.plot_info['y_zero_line']:
            if 'y_zero_line_val' in self.plot_info:
                self.PlotLine('y',self.plot_info['y_zero_line_val'])
            else:
                self.PlotLine('y',0)
        # pl.xlim(*defxlim)
        if 'save_file' in self.plot_info:
            this_dir = '/'.join(self.plot_info['save_file'].split('/')[:-1])
            mkdir_p(this_dir)
            self.plot_info['save_file'] = self.plot_info['save_file'].replace('.pdf','')+'.pdf'
            self.figure.savefig(self.plot_info['save_file'])
            self.WriteData()
        if 'show_figure' in self.plot_info and self.plot_info['show_figure']:
            self.figure.show()
            # pl.show()
        # self.ClearFigure()

    def PlotLine(self,axis,value):
        if 'x' in axis:
            self.plot_plane.axhline(value,linestyle='-',color='black')
        elif 'y' in axis:
            self.plot_plane.axvline(value,linestyle='-',color='black')


    def get_xlim(self):
        return self.plot_plane.get_xlim()
    #     data = np.array(self.plot_data.loc['x_data'].values).flatten()
    #     print 'DEBUG',data
    #     print 'DEBUG2',[np.min(data),np.max(data)]
    #     return [np.min(data),np.max(data)]
    #     # return self.plot_plane.get_xlim()
    #
    def get_ylim(self):
        return self.plot_plane.get_ylim()
    #     data = np.array(self.plot_data.loc['y_data'].values).flatten()
    #     return [np.min(data),np.max(data)]
    #     # return self.plot_plane.get_ylim()

    def WriteData(self):
        outDict = ODNested()
        if hasattr(self,'window_size'):
            outDict['window_size_x'] = self.window_size[0]
            outDict['window_size_y'] = self.window_size[1]
        outDict['Info'] = Series_TO_ODict(self.plot_info)
        for iplot_key,iplot_data in self.plot_data.items():
            is_vary = '_vary' in iplot_data['type']
            key_str = str(iplot_key).replace(' ','_')
            for ikey,idata in iplot_data.items():
                dict_key = ikey.replace(' ','_')
                if ikey == 'color':
                    idata = fmt_Qt5_col(idata)
                if '_data' not in ikey and 'key_select' not in ikey:
                    if hasattr(idata,'name'):
                        outDict[key_str][dict_key] = idata.name
                    elif hasattr(idata,'__name__'):
                        outDict[key_str][dict_key] = idata.__name__
                    else:
                        if idata is None:
                            outDict[key_str][dict_key] = 'None'
                        else:
                            outDict[key_str][dict_key] = idata
            if 'plot' in iplot_data['type'] or 'error_bar' in iplot_data['type'] or 'scatter' in  iplot_data['type'] or 'error_band' in iplot_data['type'] :
                if 'plot' in iplot_data['type']:
                    yerr_test,xerr_test = False,False
                else:
                    yerr_test = 'yerr_data' in iplot_data and iplot_data['yerr_data'] is not None and \
                        iplot_data['yerr_data'] is not False
                    xerr_test = 'xerr_data' in iplot_data and iplot_data['xerr_data'] is not None and \
                        iplot_data['xerr_data'] is not False
                if is_vary:
                    for ic,ikey in enumerate(iplot_data['key_select']):
                        if isinstance(ikey,str):
                            outDict[key_str]['key_select_'+str(ic)] = ikey
                    if isinstance(iplot_data['x_data'],pa.Series):
                        this_xdata = iplot_data['x_data'].loc[tuple_wrap(iplot_data['key_select'])].values
                    elif isinstance(iplot_data['x_data'],str) and iplot_data['x_data'] == 'from_keys':
                        # iindex = list(iplot_data['key_select']).index(slice(None))
                        this_key_select = Un_Fix_Key_Select(tuple_wrap(iplot_data['key_select']))
                        # print('DEBUG')
                        # print(this_key_select)
                        # print(iplot_data['y_data'])
                        # print(iplot_data['y_data'].loc[this_key_select])
                        this_xdata = list(iplot_data['y_data'].loc[this_key_select].index)
                        if isinstance(this_xdata[0],(list,tuple,np.ndarray)) and len(this_xdata[0]) > 1:
                            slice_loc = tuple_wrap(iplot_data['key_select']).index(slice(None))
                            if isinstance(this_xdata[0],(list,tuple,np.ndarray)):
                                this_xdata = [ix[slice_loc] for ix in this_xdata]
                            this_xdata = list(map(get_val_float_0,this_xdata))
                        elif isinstance(this_xdata[0],(list,tuple,np.ndarray)):
                            this_xdata = [get_val_float_0(ix[0]) for ix in this_xdata]
                        elif isinstance(this_xdata,(list,tuple,np.ndarray)):
                            this_xdata = [get_val_float_0(ix) for ix in this_xdata]
                    else:
                        this_xdata = iplot_data['x_data']
                    this_key_select = Un_Fix_Key_Select(tuple_wrap(iplot_data['key_select']))
                    this_ydata = iplot_data['y_data'].loc[this_key_select].values
                    if yerr_test:
                        this_yerrdata = iplot_data['yerr_data'].loc[this_key_select].values
                    if xerr_test:
                        this_xerrdata = iplot_data['xerr_data']
                else:
                    this_xdata = iplot_data['x_data']
                    this_ydata = iplot_data['y_data']
                    if yerr_test:
                        this_yerrdata = iplot_data['yerr_data']
                    if xerr_test:
                        this_xerrdata = iplot_data['xerr_data']
                    if isinstance(this_xdata,pa.Series):
                        this_xdata = this_xdata[this_key_select].values
                    elif isinstance(this_xdata,str) and this_xdata == 'from_keys':
                        this_xdata = list(this_ydata.index)
                        if isinstance(this_xdata[0],(tuple,list,np.ndarray)):
                            this_xdata = [ix[0] for ix in this_xdata]
                        this_xdata = list(map(get_val_float_0,this_xdata))
                if iplot_data['x_scale'] != 1:
                    this_xdata = np.array(this_xdata)*iplot_data['x_scale']
                    if xerr_test:
                        this_xerrdata = np.array(this_xerrdata)*iplot_data['x_scale']
                if iplot_data['scale'] != 1:
                    this_ydata = np.array(this_ydata)*iplot_data['scale']
                    if yerr_test:
                        this_yerrdata = np.array(this_yerrdata)*iplot_data['scale']
                if yerr_test:
                    if xerr_test:
                        for ix,ixerr,iy,iyerr in zip(   this_xdata,this_xerrdata,
                                                        this_ydata,this_yerrdata):
                            if minfmt < abs(ix) < maxfmt: frmtflag = 'f'
                            else: frmtflag = 'e'
                            outDict[key_str]['data'][('x_{0:5.10'+frmtflag+'} {1:5.10'+frmtflag+'}').format(ix,ixerr).replace('+','')] = AvgStdToFormat(iy,iyerr)
                    else:
                        for ix,iy,iyerr in zip(this_xdata,this_ydata,this_yerrdata):
                            if minfmt < abs(ix) < maxfmt: frmtflag = 'f'
                            else: frmtflag = 'e'
                            outDict[key_str]['data'][('x_{0:5.10'+frmtflag+'}').format(ix).replace('+','')] = AvgStdToFormat(iy,iyerr)
                else:
                    for ix,iy in zip(this_xdata,this_ydata):
                        if minfmt < abs(ix) < maxfmt: frmtflag = 'f'
                        else: frmtflag = 'e'
                        outDict[key_str]['data'][('x_{0:5.10'+frmtflag+'}').format(ix).replace('+','')] = '{0:10.20e}'.format(iy)
            elif 'vline' in iplot_data['type']:
                if isinstance(iplot_data['x_data'],pa.Series) and is_vary:
                    plotval = iplot_data['x_data'].loc[tuple_wrap(iplot_data['key_select'])]
                elif isinstance(iplot_data['x_data'],(list,tuple,np.ndarray)):
                    plotval = iplot_data['x_data'][0]
                else:
                    plotval = iplot_data['x_data']
                if iplot_data['x_scale'] != 1:
                    plotval = type(plotval)(np.array(plotval)*iplot_data['x_scale'])
                outDict[key_str]['vline_val'] = '{0:10.20e}'.format(plotval)
            elif 'hline' in iplot_data['type']:
                if isinstance(iplot_data['y_data'],pa.Series) and is_vary:
                    plotval = iplot_data['y_data'].loc[tuple_wrap(iplot_data['key_select'])]
                    if isinstance(plotval,pa.Series):
                        plotval = plotval.values[0]
                elif isinstance(iplot_data['y_data'],(list,tuple,np.ndarray)):
                    plotval = iplot_data['y_data'][0]
                else:
                    plotval = iplot_data['y_data']
                if iplot_data['scale'] != 1:
                    plotval = type(plotval)(np.array(plotval)*iplot_data['scale'])
                outDict[key_str]['hline_val'] = '{0:10.20e}'.format(plotval)
            elif 'fit' in iplot_data['type']:
                if is_vary:
                    this_key = tuple_wrap(iplot_data['key_select'])
                    if any([type(jkey) == slice for jkey in this_key]):
                        this_key = iplot_data['fit_class'].index[0]
                    if this_key in iplot_data['fit_class']:
                        this_plot = iplot_data['fit_class'].loc[this_key]
                        if isinstance(this_plot,pa.Series):
                            this_plot = this_plot.iloc[0]
                    else:
                        this_plot = iplot_data['fit_class'].iloc[0]
                else:
                    this_plot = iplot_data['fit_class']
                for iparam,boot_param in this_plot.fit_data['Params'].items():
                    outDict[key_str]['fit_parameters'][iparam] = AvgStdToFormat(boot_param.Avg,boot_param.Std)
                if 'Chi2DoF' in this_plot.fit_data.iloc[0]:
                    outDict[key_str]['chi^{2}_{pdf}'] = AvgStdToFormat(this_plot.fit_data['Chi2DoF'].iloc[0].Avg,this_plot.fit_data.iloc[0]['Chi2DoF'].Std)
        WriteXml(self.HumanFile,{'Results':outDict},Bak=False)
        self.RemoveFuns()
        pickle_out = copy(self.__dict__)
        del pickle_out['figure']
        del pickle_out['plot_plane']
        pickle_out['colcyc'] = list(pickle_out['colcyc'])
        pickle_out['symcyc'] = list(pickle_out['symcyc'])
        pickle_out['shiftcyc'] = list(pickle_out['shiftcyc'])
        WritePickle(self.PickleFile,pickle_out,Bak=False)
        self.GetFuns()

    def FixRead(self,readdict):
        readdict['colcyc'] = iter(colourset8)
        readdict['symcyc'] = iter(markerset)
        readdict['shiftcyc'] = iter(shiftset)
        if 'save_file' in list(self.plot_info.keys()):
            readdict['plot_info']['save_file'] = self.plot_info['save_file']
        if 'show_figure' in list(self.plot_info.keys()):
            readdict['plot_info']['show_figure'] = self.plot_info['show_figure']
        if 'HumanFile' in list(self.plot_info.keys()):
            readdict['HumanFile'] = self.HumanFile
        if 'PickleFile' in list(self.plot_info.keys()):
            readdict['PickleFile'] = self.PickleFile
        return readdict

    def LoadPickle(self,DefWipe=True,ForceWindowSize=False):
        # print 'filename ' , self.PickleFile

        if os.path.isfile(self.PickleFile) and not DefWipe:
            # print 'Loading Pickle for ' , self.name
            loadeddict = ReadPickleWrap(self.PickleFile)
            loadeddict = self.FixRead(loadeddict)
            self.__dict__.update(loadeddict)
            self.GetFuns()
            if ForceWindowSize is not False:
                self.MakePlane(window_size=ForceWindowSize, overwrite=True)
            else:
                self.MakePlane()
            self.AddNulls()
        elif os.path.isfile(self.PickleFile.replace('.py3p','.p')):
            loadeddict = ReadPickleWrap(self.PickleFile.replace('.py3p','.p'))
            loadeddict = self.FixRead(loadeddict)
            self.__dict__.update(loadeddict)
            self.PickleFile = self.PickleFile.replace('.p','.py3p')
            self.GetFuns()
            if ForceWindowSize is not False:
                self.MakePlane(window_size=ForceWindowSize, overwrite=True)
            else:
                self.MakePlane()
            self.AddNulls()
        elif not os.path.isfile(self.PickleFile):
            print('Plot FNF:',self.PickleFile)

    def MakePlane(self,window_size = [12.,12.275],overwrite=False):
        if overwrite or not hasattr(self,'window_size'):
            self.window_size = window_size
        self.figure = pl.figure(figsize=self.window_size)
        self.plot_plane = self.figure.add_subplot(111)


    def ClearFigure(self):
        self.figure.clf()
        del self.plot_plane
        del self.figure
        # self.plot_plane.cla()
        # self.figure.clf()

    def ShowFigure(self):
        self.figure.show()


def Test_Papers_FF():
    import QuantityLists as ql
    this_info = pa.Series()
    mkdir_p('/home/jackdra/LQCD/Scripts/Python_Analysis/TestGraphs/')
    this_info['save_file'] = '/home/jackdra/LQCD/Scripts/Python_Analysis/TestGraphs/Paper_FF.pdf'
    this_info['title'] = r'Neutron EDM Chiral Plot (Not $\alpha$ Rotated)'
    this_info['xlabel'] = r'$m_{\pi}[MeV]$'
    this_info['ylabel'] = r'$d_{n}[efm]$'
    data_plot = Plotting(plot_info=this_info)
    from itertools import cycle
    colcyc = cycle(colourset8)
    symcyc = cycle(markerset)
    shiftcyc = cycle(shiftset)

    for icol,coldata in ql.FF_1512_00566_plot.items():
        coldata['color'] = next(colcyc)
        coldata['symbol'] = next(symcyc)
        coldata['shift'] = next(shiftcyc)
        data_plot.AppendData(coldata)

    for icol,coldata in ql.FF_1510_05823_plot.items():
        coldata['color'] = next(colcyc)
        coldata['symbol'] = next(symcyc)
        coldata['shift'] = next(shiftcyc)
        data_plot.AppendData(coldata)

    for icol,coldata in ql.FF_1502_02295_plot.items():
        coldata['color'] = next(colcyc)
        coldata['symbol'] = next(symcyc)
        coldata['shift'] = next(shiftcyc)
        data_plot.AppendData(coldata)

    for icol,coldata in ql.FF_0808_1428_plot.items():
        coldata['color'] = next(colcyc)
        coldata['symbol'] = next(symcyc)
        coldata['shift'] = next(shiftcyc)
        data_plot.AppendData(coldata)

    for icol,coldata in ql.FF_0512004_plot.items():
        coldata['color'] = next(colcyc)
        coldata['symbol'] = next(symcyc)
        coldata['shift'] = next(shiftcyc)
        data_plot.AppendData(coldata)

    for icol,coldata in ql.FF_0505022_plot.items():
        coldata['color'] = next(colcyc)
        coldata['symbol'] = next(symcyc)
        coldata['shift'] = next(shiftcyc)
        data_plot.AppendData(coldata)
    data_plot.PlotAll()
    data_plot.ShowFigure()
    data_plot.PrintData()
    return data_plot



def Test_Paper_Single():
    import QuantityLists as ql
    this_info = pa.Series()
    mkdir_p('/home/jackdra/LQCD/Scripts/Python_Analysis/TestGraphs')
    this_info['save_file'] = '/home/jackdra/LQCD/Scripts/Python_Analysis/TestGraphs/Paper_FF_Single.pdf'
    this_info['title'] = r'Neutron EDM Chiral Plot (Not $\alpha$ Rotated)'
    this_info['xlabel'] = r'$m_{\pi}[MeV]$'
    this_info['ylabel'] = r'$d_{n}[efm]$'
    data_plot = Plotting(plot_info=this_info)
    from itertools import cycle
    colcyc = cycle(colourset8)
    symcyc = cycle(markerset)
    shiftcyc = cycle(shiftset)

    for icol,coldata in ql.FF_0512004_plot.items():
        coldata['color'] = next(colcyc)
        coldata['symbol'] = next(symcyc)
        coldata['shift'] = next(shiftcyc)
        data_plot.AppendData(coldata)
    data_plot.PlotAll()
    data_plot.ShowFigure()
    data_plot.PrintData()
    return data_plot

def Test_Plot():
    xdata = [1,2,3]
    ydata = [4,5,6]
    yerrdata = [1,0.5,0.6]
    ydata_2 = [5,6,7]
    yerrdata_2 = [2,1.5,0.6]
    columns = list(zip(*[   np.array(['A','B','A','B','A','B']),
                            np.array(list(map(str,[1,1,1,2,2,2])))]))
    indicies = pa.MultiIndex.from_tuples(columns,names=['upper','number'])
    yseries = pa.Series(ydata+ydata_2,index=indicies)
    yerrseries = pa.Series(yerrdata+yerrdata_2,index=indicies)
    # this_data = [[xdata,ydata,yerrdata],[ydata,xdata,yerrdata]]

    this_data = pa.DataFrame()
    hold_series = pa.Series()
    hold_series['type'] = 'error_bar_vary'
    hold_series['key_select'] = (slice(None),'1')
    hold_series['x_data'] = xdata
    hold_series['y_data'] = yseries
    hold_series['yerr_data'] = yerrseries
    hold_series['label'] = 'test data 1'
    this_data['test_1'] = hold_series


    hold_series = pa.Series()
    hold_series['type'] = 'error_bar_vary'
    hold_series['key_select'] = ('A',slice(None))
    hold_series['x_data'] = xdata
    hold_series['y_data'] = yseries
    hold_series['yerr_data'] = yerrseries
    hold_series['label'] = 'test data 2'
    this_data['test_2'] = hold_series


    this_info = pa.Series()
    this_info['save_file'] = '/home/jackdra/LQCD/Scripts/Python_Analysis/TestGraphs/plotting_test.pdf'
    this_info['title'] = 'Test Graph'
    this_info['xlabel'] = 'Test x'
    this_info['ylabel'] = 'Test y'
    # this_info['xlims'] = [0,10]
    # this_info['ylims'] = [0,15]
    data_plot = Plotting(plot_data=this_data,plot_info=this_info)
    data_plot.PlotAll()
    data_plot.ShowFigure()
    data_plot.PrintData()
    return data_plot

def mpi701_plot():
    load_file = '/home/jackdra/LQCD/Scripts/Python_Analysis/TestGraphs/mpi701.csv'
    data = pa.read_csv(load_file)
    xdata = data.iloc[:,0]
    ydata = data['Q']
    yerrdata = data['Qerr']
    # this_data = [[xdata,ydata,yerrdata],[ydata,xdata,yerrdata]]

    this_data = pa.DataFrame()
    hold_series = pa.Series()
    hold_series['type'] = 'error_bar'
    hold_series['x_data'] = xdata
    hold_series['y_data'] = ydata
    hold_series['yerr_data'] = yerrdata
    hold_series['label'] = r'$m_{\pi}=701 MeV$'
    this_data['mpi701'] = hold_series

    load_file = '/home/jackdra/LQCD/Scripts/Python_Analysis/TestGraphs/mpi701_updown.csv'
    data = pa.read_csv(load_file)
    xdata = data['Q2']
    yup,ydown = data['yup'],data['ydown']
    yavg = (yup + ydown)/2.
    yerr = (yup - ydown)/2.
    hold_series = pa.Series()
    hold_series['type'] = 'error_band'
    hold_series['alpha'] = fillalpha
    hold_series['x_data'] = xdata
    hold_series['y_data'] = yavg
    hold_series['yerr_data'] = yerr
    hold_series['label'] = r'$Par2 = 18.6(478)\times 10^{-4}$'
    this_data['mpi701_fitband'] = hold_series

    this_info = pa.Series()
    this_info['save_file'] = load_file.replace('.csv','.pdf')
    this_info['title'] = r'$m_{\pi}=701 MeV$ Form Factor plot'
    this_info['xlabel'] = r'$Q^{2}[GeV]^{2}$'
    this_info['ylabel'] = r'$F(Q^{2})$'
    # this_info['xlims'] = [0,10]
    # this_info['ylims'] = [0,15]
    data_plot = Plotting(plot_data=this_data,plot_info=this_info)
    data_plot.PlotAll()
    data_plot.ShowFigure()
    data_plot.PrintData()
    return data_plot

def Plot_E_Jangho():
    L16_jangho_file = '/home/jackdra/LQCD/Scripts/Python_Analysis/Configs/t2_energy_16X32.dat'
    L20_jangho_file = '/home/jackdra/LQCD/Scripts/Python_Analysis/Configs/t2_energy_20X40.dat'
    L28_jangho_file = '/home/jackdra/LQCD/Scripts/Python_Analysis/Configs/t2_energy_28X56.dat'
    L16_data = pa.read_csv(L16_jangho_file,delimiter=' ')
    L20_data = pa.read_csv(L20_jangho_file,delimiter=' ')
    L28_data = pa.read_csv(L28_jangho_file,delimiter=' ')




    this_data = pa.DataFrame()
    hold_series = pa.Series()
    hold_series['type'] = 'error_bar'
    hold_series['x_data'] = L16_data.iloc[:,0]
    hold_series['y_data'] = L16_data.iloc[:,1]
    hold_series['yerr_data'] = L16_data.iloc[:,2]
    hold_series['label'] = r'Jangho_L16'
    this_data['Jangho_L16'] = hold_series


    hold_series = pa.Series()
    hold_series['type'] = 'error_bar'
    hold_series['x_data'] = L20_data.iloc[:,0]
    hold_series['y_data'] = L20_data.iloc[:,1]
    hold_series['yerr_data'] = L20_data.iloc[:,2]
    hold_series['label'] = r'Jangho_L20'
    this_data['Jangho_L20'] = hold_series


    hold_series = pa.Series()
    hold_series['type'] = 'error_bar'
    hold_series['x_data'] = L28_data.iloc[:,0]
    hold_series['y_data'] = L28_data.iloc[:,1]
    hold_series['yerr_data'] = L28_data.iloc[:,2]
    hold_series['label'] = r'Jangho_L28'
    this_data['Jangho_L28'] = hold_series


    this_info = pa.Series()
    this_info['save_file'] = '/home/jackdra/LQCD/Scripts/Python_Analysis/TestGraphs/Jangho_E.pdf'
    this_info['title'] = r'Jangho E plot'
    this_info['xlabel'] = r'$\sqrt{8t} fm$'
    this_info['ylabel'] = r'$t^2 E(t)$'
    # this_info['xlims'] = [0,10]
    # this_info['ylims'] = [0,15]
    data_plot = Plotting(plot_data=this_data,plot_info=this_info)
    data_plot.PlotAll()
    data_plot.ShowFigure()
    data_plot.PrintData()
    return data_plot


def Plot_one_on_x():
    this_data = pa.DataFrame()
    hold_series = pa.Series()
    hold_series['type'] = 'plot'
    hold_series['x_data'] = np.arange(0,2,0.01)
    hold_series['y_data'] = 1/hold_series['x_data']
    hold_series['label'] = r'$1/x$'
    this_data['one_on_x'] = hold_series

    this_info = pa.Series()
    this_info['save_file'] = '/home/jackdra/LQCD/Scripts/Python_Analysis/TestGraphs/one_on_x.pdf'
    this_info['title'] = r''
    this_info['xlabel'] = r''
    this_info['ylabel'] = r''
    # this_info['xlims'] = [0,10]
    # this_info['ylims'] = [0,15]
    data_plot = Plotting(plot_data=this_data,plot_info=this_info)
    data_plot.PlotAll()
    data_plot.ShowFigure()
    data_plot.PrintData()
    return data_plot


def Plot_one_on_x2():
    this_data = pa.DataFrame()
    hold_series = pa.Series()
    hold_series['type'] = 'plot'
    hold_series['x_data'] = np.arange(0,2,0.01)
    hold_series['y_data'] = 1/(hold_series['x_data']**2)
    hold_series['label'] = r'$1/x^{2}$'
    this_data['one_on_x2'] = hold_series

    this_info = pa.Series()
    this_info['save_file'] = '/home/jackdra/LQCD/Scripts/Python_Analysis/TestGraphs/one_on_x2.pdf'
    this_info['title'] = r''
    this_info['xlabel'] = r''
    this_info['ylabel'] = r''
    # this_info['xlims'] = [0,10]
    # this_info['ylims'] = [0,15]
    data_plot = Plotting(plot_data=this_data,plot_info=this_info)
    data_plot.PlotAll()
    data_plot.ShowFigure()
    data_plot.PrintData()
    return data_plot


if __name__ == '__main__':
    # Test_Plot()

    data_plot = Test_Papers_FF()

    # # %matplotlib inline
    # import sys
    # if len(sys.argv) < 2:
    #     print 'Pease pass file name'
    #     print 'defaulting to test graph'
    #     plot_fn = '/home/jackdra/PHD/CHROMA/TestVar/scratch/resultsChiDistDebug/graphs/AlphaFull/Set1_Alpha_tsum.pdf'
    #     print plot_fn
    # elif not os.path.isfile(sys.arvg[1]):
    #     print 'File not found:'
    #     print sys.arvg[1]
    #     print 'defaulting to test graph'
    #     plot_fn = '/home/jackdra/PHD/CHROMA/TestVar/scratch/resultsChiDistDebug/graphs/AlphaFull/Set1_Alpha_tsum.pdf'
    #     print plot_fn
    # else:
    #     plot_fn = sys.argv[1]
    # this_info = pa.Series()
    # this_info['save_file'] = plot_fn
    # this_info['show_figure'] = True
    # data_plot = Plotting(plot_info=this_info)
    # data_plot.LoadPickle(DefWipe=False)
    # data_plot.PlotAll()
    # data_plot.ShowFigure()
    # data_plot.PrintData()
