#!/usr/bin/env python

from traits.etsconfig.api import ETSConfig
ETSConfig.toolkit = 'qt'

import traits.api as ta
import traitsui.api as tua
import PlotData as jpl
import QuantityLists as ql
import SetsOfFits as sff
import pandas as pad
from copy import deepcopy
from collections import OrderedDict
from BootStrapping import BootStrap
from MiscFuns import convert_to_hex,op_dict,op_dict_std
import pickle as pickle
import os
import matplotlib.pyplot as pl
import numpy as np
import time


# jpl.tpad,tsize,xsize,ysize = 20,40,60,60
# legcol,legsize,legalpha = 1,25,0.7
# linewidth,dotsize,errorbarlen = 2,10,6
# def_tick_size = 30
# def_grid = True
# jpl.linewidth,dotsize,errorbarlen = 4,20,12
# def_tick_size = 55
# jpl.def_grid = False

def_arrow_styles = ['-','->','-[','|-|','-|>','<-','<->','<|-','<|-|>','fancy','simple','wedge']

data_plot = None
backup_plot = []
hold_file = './prev_file.py3p'


x_funs = OrderedDict()
x_funs['None'] = None
x_funs['square'] = np.square
x_funs['sqrt'] = np.sqrt
x_funs['log'] = np.log
x_funs['exp'] = np.exp

def check_index_and_None(this_data,this_index,this_key,def_val = False):
    if this_index not in this_data.index or this_data.loc[this_index,this_key] is None:
        return def_val
    else:
        return this_data.loc[this_index,this_key]


def fmt_Range(this_list):
    this_list = list(OrderedDict.fromkeys(this_list))
    # return map(get_val,this_list)
    return this_list

class DataFrame( ta.HasTraits ):

    Load = ta.Button()
    Apply = ta.Button()
    Duplicate = ta.Button()
    Wipe = ta.Button()
    Undo_Wipe = ta.Button()
    Apply_Move = ta.Button()
    Create_Arrow = ta.Button()
    Combine = ta.Button()


    resample = ta.Str('1000')
    transpose_table = ta.Bool(False)
    chi_min = ta.Float(0)
    chi_max = ta.Float(2)
    fix_min = ta.Int(0)
    cut_max = ta.Int(100)
    Hist_Comb_Min = ta.Button()
    fix_max = ta.Int(100)
    cut_min = ta.Int(0)
    Hist_Comb_Max = ta.Button()

    this_aindex = ta.Int(0)
    animate_dim = ta.Int(0)
    Animate = ta.Button()

    Dump_Fit_Table = ta.Button()
    fit_file_name = ta.File()
    fmt_latex = ta.Bool(True)

    x_range_min = ta.Int(0)
    x_range_max = ta.Int(-1)
    x_increment = ta.Int(1)
    x_scale = ta.Float(1)
    x_function = ta.Enum(list(x_funs.keys()))
    shift = ta.Float(1.0)
    shift_overlap = ta.Float(0.005)
    do_band = ta.Bool(False)
    band_index = ta.Int(-1)
    show_band_eval = ta.Bool(True)

    arrow_x = ta.Float(0.0)
    arrow_y = ta.Float(0.0)
    arrow_x_text = ta.Float(1.0)
    arrow_y_text = ta.Float(1.0)
    arrow_text_size = ta.Float(30)
    arrow_color = ta.Str('black')
    arrow_style = ta.Enum(def_arrow_styles)


    key_list = ta.List([''])
    this_key = ta.Enum(values='key_list')
    this_key_2 = ta.Enum(values='key_list')
    plot_this = ta.Bool(True)
    label = ta.Str('test')
    Physical_Units = ta.Bool(True)
    plot_err = ta.Enum(['No','relative','absolute'])
    supress_legend = ta.Bool(False)
    color = ta.Color(None)
    symbol = ta.Str('o')
    Median = ta.Bool(False)
    scale_y = ta.Float(1.0)
    shift_y = ta.Float(0)
    fit_dim_pick = ta.Int(1)
    ShowEval = ta.List(ta.Float,[0])
    extrap_fade = ta.Float(1.5)
    Hairline_plot = ta.Bool(False)
    line_alpha_scale = ta.Float(20)
    line_alpha_abs = ta.Float(1)
    Force_Param_Name = ta.Str('')
    Force_Param_Units = ta.Str('')
    fit_custom_xrange = ta.Bool(False)
    fit_x_min = ta.Str('FromDataMin')
    fit_x_max = ta.Str('ToDataMax')


    suppress_key_index = ta.Bool(False)
    dim_list = ta.List([''])
    pick_dimension = ta.Enum(values='dim_list')
    test_pick_dim = ta.Bool(False)
    dim_1_list = ta.List([''])
    dim_2_list = ta.List([''])
    dim_3_list = ta.List([''])
    dim_4_list = ta.List([''])
    dim_5_list = ta.List([''])
    dim_6_list = ta.List([''])
    dim_7_list = ta.List([''])
    dim_8_list = ta.List([''])
    dim_1 = ta.Enum([1,2,3,4])
    dim_1 = ta.Enum(values = 'dim_1_list')
    dim_2 = ta.Enum(values = 'dim_2_list')
    dim_3 = ta.Enum(values = 'dim_3_list')
    dim_4 = ta.Enum(values = 'dim_4_list')
    dim_5 = ta.Enum(values = 'dim_5_list')
    dim_6 = ta.Enum(values = 'dim_6_list')
    dim_7 = ta.Enum(values = 'dim_7_list')
    dim_8 = ta.Enum(values = 'dim_8_list')
    test_dim_1 = ta.Bool(False)
    test_dim_2 = ta.Bool(False)
    test_dim_3 = ta.Bool(False)
    test_dim_4 = ta.Bool(False)
    test_dim_5 = ta.Bool(False)
    test_dim_6 = ta.Bool(False)
    test_dim_7 = ta.Bool(False)
    test_dim_8 = ta.Bool(False)
    vary_len = 0
    move_to_index = ta.Int(0)

    combine_operation = ta.Enum(list(op_dict.keys()))
    new_label = ta.Str('')

    view = tua.View(
        tua.Group(
        tua.Item('this_key'),
        tua.Item('plot_this'),
        tua.Item('label'),
        tua.Item('Physical_Units'),
        tua.Item('plot_err'),
        tua.Item('color',style='custom'),
        tua.Item('symbol'),
        tua.Item('Median'),
        tua.Item('scale_y'),
        tua.Item('shift_y'),
        tua.Item('supress_legend'),
        tua.Item('suppress_key_index'),
        tua.Item('pick_dimension',visible_when='test_pick_dim'),
        tua.Item('dim_1',visible_when='test_dim_1',name='pie'),
        tua.Item('dim_2',visible_when='test_dim_2'),
        tua.Item('dim_3',visible_when='test_dim_3'),
        tua.Item('dim_4',visible_when='test_dim_4'),
        tua.Item('dim_5',visible_when='test_dim_5'),
        tua.Item('dim_6',visible_when='test_dim_6'),
        tua.Item('dim_7',visible_when='test_dim_7'),
        tua.Item('dim_8',visible_when='test_dim_8'),
        # tua.Item('Load', show_label=False),
        tua.Item('Apply', show_label=False),
        tua.Item('Duplicate', show_label=False),
        tua.Item('Wipe', show_label=False),
        tua.Item('Undo_Wipe', show_label=False),
        tua.Item('move_to_index'),
        tua.Item('Apply_Move', show_label=False),
        label='Main'
        ),tua.Group(
        tua.Item('this_key'),
        tua.Item('plot_this'),
        tua.Item('x_scale'),
        tua.Item('x_function'),
        tua.Item('shift'),
        tua.Item('shift_overlap'),
        tua.Item('x_range_min'),
        tua.Item('x_range_max'),
        tua.Item('x_increment'),
        tua.Item('do_band'),
        tua.Item('band_index',enabled_when='do_band'),
        tua.Item('show_band_eval',enabled_when='do_band'),
        tua.Item('Apply', show_label=False),
        tua.Item('Duplicate', show_label=False),
        tua.Item('Wipe', show_label=False),
        tua.Item('Undo_Wipe', show_label=False),
        tua.Item('move_to_index'),
        tua.Item('Apply_Move', show_label=False),
        tua.Item('animate_dim'),
        tua.Item('Animate', show_label=False),
        label='x_axis_modify'
        ),tua.Group(
        tua.Item('this_key'),
        tua.Item('plot_this'),
        tua.Item('extrap_fade'),
        tua.Item('Hairline_plot'),
        tua.Item('line_alpha_scale',visible_when='Hairline_plot'),
        tua.Item('line_alpha_abs',visible_when='Hairline_plot'),
        tua.Item('fit_custom_xrange'),
        tua.Item('fit_x_min',enabled_when='fit_custom_xrange'),
        tua.Item('fit_x_max',enabled_when='fit_custom_xrange'),
        tua.Item('ShowEval'),
        tua.Item('fit_dim_pick'),
        tua.Item('Force_Param_Name'),
        tua.Item('Force_Param_Units'),
        tua.Item('Apply', show_label=False),
        tua.Item('Duplicate', show_label=False),
        tua.Item('Wipe', show_label=False),
        tua.Item('Undo_Wipe', show_label=False),
        label='Fitting_stuff'
        ),tua.Group(
        tua.Item('this_key'),
        tua.Item('plot_this'),
        tua.Item('arrow_x'),
        tua.Item('arrow_y'),
        tua.Item('arrow_x_text'),
        tua.Item('arrow_y_text'),
        tua.Item('arrow_text_size'),
        tua.Item('arrow_color'),
        tua.Item('arrow_style'),
        tua.Item('supress_legend'),
        tua.Item('label'),
        tua.Item('Apply', show_label=False),
        tua.Item('Duplicate', show_label=False),
        tua.Item('Wipe', show_label=False),
        tua.Item('Undo_Wipe', show_label=False),
        tua.Item('Create_Arrow',show_label=False),
        tua.Item('fit_file_name'),
        tua.Item('transpose_table'),
        tua.Item('Dump_Fit_Table', show_label=False),
        label='Arrow'
        ),tua.Group(
        tua.Item('this_key'),
        tua.Item('this_key_2'),
        tua.Item('combine_operation'),
        tua.Item('new_label'),
        tua.Item('Combine', show_label=False),
        tua.Item('fit_file_name'),
        tua.Item('fmt_latex'),
        tua.Item('resample'),
        tua.Item('chi_min'),
        tua.Item('chi_max'),
        tua.Item('fix_min'),
        tua.Item('cut_max'),
        tua.Item('Hist_Comb_Min', show_label=False),
        tua.Item('fix_max'),
        tua.Item('cut_min'),
        tua.Item('Hist_Comb_Max', show_label=False),
        tua.Item('fit_file_name'),
        tua.Item('transpose_table'),
        tua.Item('Dump_Fit_Table', show_label=False),
        label='Combine'
        ),
        buttons=['OK'],
        resizable=True
    )

    def _Hist_Comb_Max_fired(self):
        global data_plot
        this_data = data_plot.plot_data.loc[:,self.this_key]
        if 'histogram_vary' == this_data['type']:
            if 'boot_data' in this_data and hasattr(this_data['boot_data'],'Comb_Max_fitr'):
                app_data = deepcopy(this_data)
                app_data['type'] = 'histogram_vary'
                if self.resample.lower() == 'false':
                    app_data['boot_data'] = app_data['boot_data'].Comb_Max_fitr(self.fix_max,cut_min=self.cut_min,
                                                                                chi_min=self.chi_min,chi_max=self.chi_max)
                else:
                    app_data['boot_data'] = app_data['boot_data'].Comb_Max_fitr(self.fix_max,cut_min=self.cut_min,
                                                                                resample=int(self.resample),
                                                                                chi_min=self.chi_min,chi_max=self.chi_max)
                app_data['label'] = app_data['label']+'_fixmax='+str(self.fix_max)+'_cutmin='+str(self.cut_min)
                app_data['key_select'] = app_data['boot_data'].index[0]
                data_plot.AppendData(app_data)
                self.key_list = list(data_plot.plot_data.keys())
                self._Load_fired()
                self._Apply_fired()

    def _Hist_Comb_Min_fired(self):
        global data_plot
        this_data = data_plot.plot_data.loc[:,self.this_key]
        if 'histogram_vary' == this_data['type']:
            if 'boot_data' in this_data and hasattr(this_data['boot_data'],'Comb_Min_fitr'):
                app_data = deepcopy(this_data)
                app_data['type'] = 'histogram_vary'
                if self.resample.lower() == 'false':
                    app_data['boot_data'] = app_data['boot_data'].Comb_Min_fitr(self.fix_min,cut_max=self.cut_max,
                                                                    chi_min=self.chi_min,chi_max=self.chi_max)
                else:
                    app_data['boot_data'] = app_data['boot_data'].Comb_Min_fitr(self.fix_min,cut_max=self.cut_max,
                                                                    resample=int(self.resample),
                                                                    chi_min=self.chi_min,chi_max=self.chi_max)
                app_data['label'] = app_data['label']+'_fixmin='+str(self.fix_min)+'_cutmax='+str(self.fit_max)
                app_data['key_select'] = app_data['boot_data'].index[0]
                data_plot.AppendData(app_data)
                self.key_list = list(data_plot.plot_data.keys())
                self._Load_fired()
                self._Apply_fired()

    def _Dump_Fit_Table_fired(self):
        global data_plot
        this_data = data_plot.plot_data.loc[:,self.this_key]
        if 'type' in this_data and ('fit_vary' == this_data['type'] or 'histogram_vary' == this_data['type']):
            if 'fit_class' in this_data and hasattr(this_data['fit_class'],'Get_Formatted_Table'):
                this_type = 'fit_class'
                do_GFT = True
            elif 'boot_data' in this_data:
                this_type = 'boot_data'
                do_GFT = hasattr(this_data['boot_data'],'Get_Formatted_Table')
            else:
                do_GFT = False
            if do_GFT:
                this_xml = this_data[this_type].Get_Formatted_Table(fmt_latex=self.fmt_latex)
            elif isinstance(this_data[this_type].iloc[0],BootStrap):
                this_xml = this_data[this_type].apply(lambda x : x.MakeValAndErr(Dec=2,latex=True))
            else:
                this_xml = this_data[this_type]
            if self.transpose_table:
                this_xml = this_xml.transpose()
            with open(self.fit_file_name,'w') as f:
                f.write(this_xml.to_latex(escape=False))
            with open(self.fit_file_name+'.py3p','wb') as f:
                if hasattr(this_data[this_type],'RemoveFuns'): this_data[this_type].RemoveFuns()
                pickle.dump(this_data[this_type],f)
                if hasattr(this_data[this_type],'GetFuns'): this_data[this_type].GetFuns()
        else:
            print('Type of data is not fit vary')


    def _Combine_fired(self):
        first_data = deepcopy(data_plot.plot_data.loc[:,self.this_key])
        second_data = deepcopy(data_plot.plot_data.loc[:,self.this_key_2])
        if first_data['type'] != second_data['type']:
            raise TypeError('both data types must be equal!!')
        if 'plot' in first_data['type'] or 'error_bar' in first_data['type']:
            if isinstance(first_data['y_data'],(list,tuple)):
                first_data['y_data'] = np.array(first_data['y_data'])
            if isinstance(second_data['y_data'],(list,tuple)):
                second_data['y_data'] = np.array(second_data['y_data'])
            if first_data['scale'] != 1:
                first_data['y_data'] = first_data['y_data']*first_data['scale']
            if second_data['scale'] != 1:
                second_data['y_data'] = second_data['y_data']*second_data['scale']

            first_data['y_data'] = op_dict[self.combine_operation](first_data['y_data'], second_data['y_data'])
            if 'yerr_data' in first_data and 'yerr_data' in second_data:
                if isinstance(first_data['yerr_data'],(list,tuple)):
                    first_data['yerr_data'] = np.array(first_data['yerr_data'])
                if isinstance(second_data['yerr_data'],(list,tuple)):
                    second_data['yerr_data'] = np.array(second_data['yerr_data'])
                if first_data['scale'] != 1:
                    first_data['yerr_data'] = first_data['y_data']*first_data['scale']
                if second_data['scale'] != 1:
                    second_data['yerr_data'] = second_data['y_data']*second_data['scale']
                first_data['yerr_data'] = op_dict_std[self.combine_operation](first_data['yerr_data'],
                                                                              second_data['yerr_data'],
                                                                              first_data['y_data'],
                                                                              second_data['y_data'])
        else:
            print(first_data['type'])
            print(second_data['type'])
            raise TypeError('data must be plot or error_bar')
        first_data['scale'] = 1.0
        first_data['label'] = self.new_label
        first_data.name = self.new_label
        data_plot.AppendData(first_data)
        self.key_list = list(data_plot.plot_data.keys())
        self._Load_fired()
        self._Apply_fired()



    def _Create_Arrow_fired(self):
        global data_plot
        app_data = deepcopy(data_plot.plot_data.loc[:,self.this_key])
        app_data['type'] = 'arrow'
        app_data['x_data'] = [self.arrow_x,self.arrow_x_text]
        app_data['y_data'] = [self.arrow_y,self.arrow_y_text]
        app_data['label'] = self.label
        app_data['supress_legend'] = self.supress_legend
        app_data['arrow_text_size'] = self.arrow_text_size
        app_data['arrow_color'] = self.arrow_color
        app_data['arrow_style'] = self.arrow_style
        data_plot.AppendData(app_data)
        self.key_list = list(data_plot.plot_data.keys())
        self._Load_fired()
        self._Apply_fired()

    def _Apply_Move_fired(self):
        global data_plot
        this_data = deepcopy(data_plot.plot_data)
        index_list = np.array(this_data.columns)
        if self.move_to_index >= len(data_plot.plot_data.index):
            print('move_to_index is not in range of 0 to ',len(data_plot.plot_data.index)-1)
            return
        from_index = list(index_list).index(self.this_key)
        if from_index > self.move_to_index:
            index_list = np.insert(index_list,self.move_to_index,self.this_key)
            index_list = np.delete(index_list,from_index+1)
        elif from_index < self.move_to_index:
            index_list = np.insert(index_list,self.move_to_index+1,self.this_key)
            index_list = np.delete(index_list,from_index)
        this_data = this_data[index_list]
        data_plot.ImportData(this_data)
        data_plot.PrintData()
        data_plot.PlotAll()
        self._Load_fired()


    def _Undo_Wipe_fired(self):
        global data_plot,backup_plot
        if len(backup_plot) == 0:
            print('Error, no backup present')
        data_plot.plot_data = backup_plot[-1]
        del backup_plot[-1]
        self.key_list = list(data_plot.plot_data.keys())
        self._Load_fired()
        self._Apply_fired()

    def _Wipe_fired(self):
        global data_plot,backup_plot
        # if self.this_key[-2] == '_':
        backup_plot.append(deepcopy(data_plot.plot_data))
        data_plot.plot_data = data_plot.plot_data.drop(self.this_key,axis=1)
        print('removed' , self.this_key)
        self.key_list = list(data_plot.plot_data.keys())
        self._Load_fired()
        self._Apply_fired()
        # else:
        #     print 'Cannot remove original data.'



    def _Duplicate_fired(self):
        global data_plot
        app_data = deepcopy(data_plot.plot_data.loc[:,self.this_key])
        data_plot.AppendData(app_data)
        self.key_list = list(data_plot.plot_data.keys())
        self._Load_fired()
        self._Apply_fired()

    def __init__(self):
        global data_plot
        self.key_list = list(data_plot.plot_data.keys())
        self.this_key = self.key_list[0]
        self.this_key_2 = self.key_list[0]
        self._Load_fired()
        if self.test_pick_dim:
            self.pick_dimension = self.dim_list[data_plot.plot_data.loc['key_select',self.this_key].index(slice(None))]
        self._Apply_fired()

    # def _label_fired(self):
    #     self._Apply_fired()
    # def _color_fired(self):
    #     self._Apply_fired()
    # def _shift_fired(self):
    #     self._Apply_fired()
    # def _symbol_fired(self):
    #     if len(self.symbol)==1:
    #         self._Apply_fired()
    #
    # def _scale_y_fired(self):
    #     self._Apply_fired()

    def _this_key_fired(self):
        self._Load_fired()

    def _this_key_2_fired(self):
        self._Load_fired()

    def _pick_dimension_fired(self):
        self._Load_fired()

    # def _dim_1_fired(self):
    #     self._Load_fired()
    # def _dim_2_fired(self):
    #     self._Load_fired()
    # def _dim_3_fired(self):
    #     self._Load_fired()
    # def _dim_4_fired(self):
    #     self._Load_fired()
    # def _dim_5_fired(self):
    #     self._Load_fired()
    # def _dim_6_fired(self):
    #     self._Load_fired()
    # def _dim_7_fired(self):
    #     self._Load_fired()
    # def _dim_8_fired(self):
    #     self._Load_fired()


    def _Load_fired(self):
        global data_plot
        self.test_pick_dim = False
        for ic in range(1,8):
            setattr(self,'test_dim_'+str(ic),False)

        if data_plot is None:
            raise EnvironmentError('plot data has not been loaded')
        # self.key_list = data_plot.plot_data.keys()
        self.plot_this = bool(check_index_and_None(data_plot.plot_data,'plot_this',self.this_key,def_val=False))


        if 'arrow' in data_plot.plot_data.loc['type',self.this_key]:
            if 'x_data' in data_plot.plot_data.index:
                self.arrow_x,self.arrow_x_text = data_plot.plot_data.loc['x_data',self.this_key][:2]
            else:
                self.arrow_x,self.arrow_x_text = 0,1
            if 'y_data' in data_plot.plot_data.index:
                self.arrow_y,self.arrow_y_text = data_plot.plot_data.loc['y_data',self.this_key][:2]
            else:
                self.arrow_y,self.arrow_y_text = 0,1
            self.arrow_text_size = check_index_and_None(data_plot.plot_data,'arrow_text_size',self.this_key,def_val = 30)
            self.arrow_color = check_index_and_None(data_plot.plot_data,'arrow_color',self.this_key,def_val = 'black')
            self.arrow_style = check_index_and_None(data_plot.plot_data,'arrow_style',self.this_key,def_val = 'fancy')



        self.label = check_index_and_None(data_plot.plot_data,'label',self.this_key,def_val = 'None')
        self.Physical_Units = check_index_and_None(data_plot.plot_data,'phys_axies',self.this_key,def_val = True)
        self.extrap_fade = check_index_and_None(data_plot.plot_data,'extrap_fade',self.this_key,def_val = 1.5)
        self.Hairline_plot = check_index_and_None(data_plot.plot_data,'hairline',self.this_key,def_val = True)
        if 'hair_alpha' in data_plot.plot_data.index:
            if isinstance(data_plot.plot_data.loc['hair_alpha',self.this_key],(tuple,list,np.ndarray)):
                self.line_alpha_scale = data_plot.plot_data.loc['hair_alpha',self.this_key][0]
                self.line_alpha_abs = data_plot.plot_data.loc['hair_alpha',self.this_key][1]
            else:
                self.line_alpha_scale = data_plot.plot_data.loc['hair_alpha',self.this_key]
                self.line_alpha_abs = 1
        else:
            self.line_alpha_scale = 20
            self.line_alpha_abs = 1
        if 'xdatarange' in data_plot.plot_data.index:
            if data_plot.plot_data.loc['xdatarange',self.this_key] == 'Data':
                self.fit_custom_xrange = False
            else:
                self.fit_custom_xrange = True
                try:
                    self.fit_x_min,self.fit_x_max = data_plot.plot_data.loc['xdatarange',self.this_key]
                except Exception as err:
                    pass
        else:
            self.fit_custom_xrange = False
        if 'plot_err' in data_plot.plot_data.index:
            for itype in ['No','relative','absolute']:
                if data_plot.plot_data.loc['plot_err',self.this_key].lower() in itype.lower():
                    self.plot_err = itype
                    break
                else:
                    self.plot_err = 'No'
        else:
            self.plot_err = 'No'
        if 'color' in data_plot.plot_data.index:
            if isinstance(data_plot.plot_data.loc['color',self.this_key],(list,tuple)):
                self.color = convert_to_hex(data_plot.plot_data.loc['color',self.this_key])
            else:
                try:
                    self.color = data_plot.plot_data.loc['color',self.this_key]
                except Exception as err:
                    pass
        else:
            self.color = 'blue'

        self.shift = float(check_index_and_None(data_plot.plot_data,'shift',self.this_key,def_val = 0.0))
        self.supress_legend = bool(check_index_and_None(data_plot.plot_data,'supress_legend',self.this_key,def_val = False))

        xrmin_hold = check_index_and_None(data_plot.plot_data,'x_range_min',self.this_key,def_val = 0)
        if not isinstance(xrmin_hold,bool):
            self.x_range_min = int(xrmin_hold)
        xrmax_hold = check_index_and_None(data_plot.plot_data,'x_range_max',self.this_key,def_val = -1)
        if not isinstance(xrmax_hold,bool):
            self.x_range_max = int(xrmax_hold)
        if self.x_range_max == 0:
            self.x_range_max = -1

        xrinc_hold = check_index_and_None(data_plot.plot_data,'x_increment',self.this_key,def_val = 1)
        if not isinstance(xrmax_hold,bool):
            self.x_increment = int(xrinc_hold)
        if self.x_increment == 0:
            self.x_increment = 1
        this_band = check_index_and_None(data_plot.plot_data,'do_band',self.this_key,def_val = False)
        if isinstance(this_band,int) and not isinstance(this_band,bool):
            self.do_band = True
            self.band_index = this_band
            self.show_band_eval = bool(check_index_and_None(data_plot.plot_data,'show_band_eval',self.this_key,def_val = False))
        else:
            self.do_band = False
            self.show_band_eval = False

        self.symbol = str(check_index_and_None(data_plot.plot_data,'symbol',self.this_key,def_val = 'None'))
        self.Median = bool(check_index_and_None(data_plot.plot_data,'Median',self.this_key,def_val = False))
        self.suppress_key_index = bool(check_index_and_None(data_plot.plot_data,'suppress_key_index',self.this_key,def_val = False))
        hold_xfun = check_index_and_None(data_plot.plot_data,'x_fun',self.this_key,def_val = 'None')
        if hasattr(hold_xfun,'__name__') and not isinstance(hold_xfun.__name__,bool):
            self.x_function = hold_xfun.__name__
        else:
            self.x_function = 'None'

        self.x_scale = float(check_index_and_None(data_plot.plot_data,'x_scale',self.this_key,def_val = 1.0))
        self.scale_y = float(check_index_and_None(data_plot.plot_data,'scale',self.this_key,def_val = 1.0))
        self.shift_y = float(check_index_and_None(data_plot.plot_data,'y_shift',self.this_key,def_val = 0.0))

        # if len(self.ShowEval) == 0:
        #     pass
        hold_ShowEval = check_index_and_None(data_plot.plot_data,'ShowEval',self.this_key,def_val = [])
        if isinstance(hold_ShowEval,(list,tuple,np.ndarray)):
            self.ShowEval = list(hold_ShowEval)
        else:
            self.ShowEval = []

        self.fit_dim_pick = int(check_index_and_None(data_plot.plot_data,'xaxis',self.this_key,def_val = 0))

        self.Force_Param_Name = str(check_index_and_None(data_plot.plot_data,'FPName',self.this_key,def_val = ''))
        self.Force_Param_Units = str(check_index_and_None(data_plot.plot_data,'FPUnits',self.this_key,def_val = ''))


        if '_vary' in data_plot.plot_data.loc['type',self.this_key]:
            this_series = data_plot.plot_data.loc[:,self.this_key]
            update_dim = isinstance(this_series['key_select'],tuple) and not jpl.Py2Read
            if 'fit' in this_series['type']:
                if 'fit_class' in this_series:
                    if isinstance(this_series['fit_class'],sff.SetOfFitFuns):
                        this_fit = this_series['fit_class'].Fit_Stats_fmt['Fit']
                    else:
                        this_fit = this_series['fit_class']
                if isinstance(this_fit.index,pad.MultiIndex):
                    self.test_pick_dim = False
                    self.vary_len = len(this_fit.index.names)
                    for icd,(idim,ikey) in enumerate(zip(this_fit.index.names,this_series['key_select'])):
                        setattr(self,'test_dim_'+str(icd+1),True)
                        this_list = list(this_fit.index.get_level_values(idim))
                        setattr(self,'dim_'+str(icd+1)+'_list',fmt_Range(this_list))
                        if update_dim and not isinstance(ikey,slice):
                            setattr(self,'dim_'+str(icd+1),str(ikey))
                else:
                    self.vary_len = 1
                    self.test_pick_dim = False
                    self.test_dim_1 = True
                    self.dim_1_list = fmt_Range(list(this_fit.index))
            elif 'boot_data' in list(this_series.keys()) and this_series['boot_data'] is not None \
            and not isinstance(this_series['boot_data'],bool):
                if isinstance(this_series['boot_data'],sff.SetOfFitFuns):
                    this_boot = this_series['boot_data'].Get_Extrapolation(fmted=True)
                else:
                    this_boot = this_series['boot_data']
                if isinstance(this_boot.index,pad.MultiIndex):
                    self.test_pick_dim = False
                    self.vary_len = len(this_boot.index.names)
                    for icd,(idim,ikey) in enumerate(zip(this_boot.index.names,this_series['key_select'])):
                        setattr(self,'test_dim_'+str(icd+1),True)
                        this_list = list(this_boot.index.get_level_values(idim))
                        setattr(self,'dim_'+str(icd+1)+'_list',fmt_Range(this_list))
                        if update_dim and not isinstance(ikey,slice):
                            setattr(self,'dim_'+str(icd+1),str(ikey))
                else:
                    self.vary_len = 1
                    self.test_pick_dim = False
                    self.test_dim_1 = True
                    self.dim_1_list = fmt_Range(list(this_boot.index))
            elif 'y_data' in list(this_series.keys()) and isinstance(this_series['y_data'],pad.Series):
                if isinstance(this_series['y_data'].index,pad.MultiIndex):
                    if 'line' in data_plot.plot_data.loc['type',self.this_key]:
                        self.test_pick_dim = False
                        self.vary_len = len(this_series['y_data'].index.names)
                        for icd,idim in enumerate(this_series['y_data'].index.names):
                            setattr(self,'test_dim_'+str(icd+1),True)
                            this_list = list(this_series['y_data'].index.get_level_values(idim))
                            setattr(self,'dim_'+str(icd+1)+'_list',fmt_Range(this_list))
                    else:
                        self.test_pick_dim = True
                        self.dim_list = this_series['y_data'].index.names
                        self.vary_len = len(this_series['y_data'].index.names)
                        for icd,(idim,ikey) in enumerate(zip(this_series['y_data'].index.names,this_series['key_select'])):
                            if idim == str(self.pick_dimension):
                                setattr(self,'test_dim_'+str(icd+1),False)
                            else:
                                setattr(self,'test_dim_'+str(icd+1),True)
                                this_list = list(this_series['y_data'].index.get_level_values(idim))
                                setattr(self,'dim_'+str(icd+1)+'_list',fmt_Range(this_list))
                                if update_dim and not isinstance(ikey,slice):
                                    setattr(self,'dim_'+str(icd+1),str(ikey))
                else:
                    self.test_pick_dim = False
                    self.test_dim_1 = True
                    this_list = list(this_series['y_data'].index)
                    self.dim_1_list = fmt_Range(this_list)
            if 'x_data' in list(this_series.keys()) and isinstance(this_series['x_data'],pad.Series):
                if isinstance(this_series['x_data'].index,pad.MultiIndex):
                    self.test_pick_dim = True
                    self.dim_list = this_series['x_data'].index.names
                    self.vary_len = len(this_series['x_data'].index.names)
                    for icd,(idim,ikey) in enumerate(zip(this_series['x_data'].index.names,this_series['key_select'])):
                        if idim == str(self.pick_dimension):
                            setattr(self,'test_dim_'+str(icd+1),False)
                        else:
                            setattr(self,'test_dim_'+str(icd+1),True)
                            this_list = list(this_series['x_data'].index.get_level_values(idim))
                            setattr(self,'dim_'+str(icd+1)+'_list',fmt_Range(this_list))
                            if update_dim and not isinstance(ikey,slice):
                                setattr(self,'dim_'+str(icd+1),str(ikey))
                else:
                    self.test_pick_dim = False
                    self.test_dim_1 = True
                    self.dim_1_list = fmt_Range(list(this_series['x_data'].index))
        self.fit_file_name = data_plot.plot_info['save_file'].replace('.pdf',self.pick_dimension+'.tex')


    def _Animate_fired(self):
        self._Apply_fired()
        self.this_aindex = data_plot.Animate(this_key=self.this_key,a_index=self.animate_dim,multi_index=self.this_aindex)
        if self.this_aindex is -1:
            self.this_aindex = 0

    def _Apply_fired(self):
        # global customleg,customleglist
        # global wtitle,wxlab,wrlab
        # global legloc
        global data_plot
        this_data = deepcopy(data_plot.plot_data)
        this_data.loc['plot_this',self.this_key] = self.plot_this
        if self.label != 'None':
            this_data.loc['label',self.this_key] = self.label

        if self.Physical_Units != 'None':
            this_data.loc['phys_axies',self.this_key] = self.Physical_Units

        if self.extrap_fade != 'None' and self.extrap_fade != 0:
            this_data.loc['extrap_fade',self.this_key] = self.extrap_fade
        else:
            this_data.loc['extrap_fade',self.this_key] = 1.5
        if self.Hairline_plot != 'None':
            this_data.loc['hairline',self.this_key] = self.Hairline_plot
        alpha_hold = []
        if self.line_alpha_scale != 'None':
            if self.line_alpha_scale < 0:
                print('Warning, no negative line_alpha_scale! breaks alpha scaling... assuming abs')
            alpha_hold.append(abs(self.line_alpha_scale))
        else:
            if isinstance(this_data.loc['hair_alpha',self.this_key],(tuple,list,np.ndarray)):
                alpha_hold.append(this_data.loc['hair_alpha',self.this_key][0])
            else:
                alpha_hold.append(this_data.loc['hair_alpha',self.this_key])
        if self.fit_custom_xrange:
            this_data.loc['xdatarange',self.this_key] = [self.fit_x_min,self.fit_x_max]
        else:
            this_data.loc['xdatarange',self.this_key] = 'Data'
        if self.line_alpha_abs != 'None':
            if 0 < self.line_alpha_abs > 1:
                print('Warning, line_alpha_abs is between 0 and 1! breaks alpha scaling... assuming 1')
                alpha_hold.append(1.0)
            else:
                alpha_hold.append(self.line_alpha_abs)
        else:
            if isinstance(this_data.loc['hair_alpha',self.this_key],(tuple,list,np.ndarray)):
                alpha_hold.append(this_data.loc['hair_alpha',self.this_key][1])
        this_data.loc['hair_alpha',self.this_key] = alpha_hold
        this_data.loc['plot_err',self.this_key] = self.plot_err
        if self.color != 'None':
            this_data.loc['color',self.this_key] = self.color
        if self.supress_legend != 'None':
            this_data.loc['supress_legend',self.this_key] = self.supress_legend
        if self.shift_overlap != 'None':
            this_data.loc['shift_overlap',self.this_key] = self.shift_overlap
        if self.shift != 'None':
            this_data.loc['shift',self.this_key] = self.shift
        this_data.loc['x_range_min',self.this_key] = self.x_range_min
        this_data.loc['x_range_max',self.this_key] = self.x_range_max
        this_data.loc['x_increment',self.this_key] = self.x_increment
        if self.do_band:
            this_data.loc['do_band',self.this_key] = self.band_index
            this_data.loc['show_band_eval',self.this_key] = self.show_band_eval
        else:
            this_data.loc['do_band',self.this_key] = False
            this_data.loc['show_band_eval',self.this_key] = False
        if self.symbol != 'None':
            this_data.loc['symbol',self.this_key] = self.symbol
        if self.Median != 'None':
            this_data.loc['Median',self.this_key] = self.Median
        if self.scale_y != 'None':
            this_data.loc['scale',self.this_key] = self.scale_y
        if self.shift_y != 'None':
            this_data.loc['y_shift',self.this_key] = self.shift_y
        if self.suppress_key_index != 'None':
            this_data.loc['suppress_key_index',self.this_key] = self.suppress_key_index
        if len(list(self.ShowEval)) > 0:
            this_data.loc['ShowEval',self.this_key] = list(self.ShowEval)
        else:
            this_data.loc['ShowEval',self.this_key] = None
        this_data.loc['xaxis',self.this_key] = int(self.fit_dim_pick)
        if self.Force_Param_Name != 'None':
            this_data.loc['FPName',self.this_key] = self.Force_Param_Name
        if self.Force_Param_Units != 'None':
            this_data.loc['FPUnits',self.this_key] = self.Force_Param_Units
        if self.Force_Param_Units != 'None':
            this_data.loc['FPUnits',self.this_key] = self.Force_Param_Units

        if 'arrow' in this_data.loc['type',self.this_key]:
            if self.arrow_x != 'None' and self.arrow_x_text != 'None':
                this_data.loc['x_data',self.this_key] = [self.arrow_x,self.arrow_x_text]
            if self.arrow_y != 'None' and self.arrow_y_text != 'None':
                this_data.loc['y_data',self.this_key] = [self.arrow_y,self.arrow_y_text]
            if self.arrow_text_size != 'None':
                this_data.loc['arrow_text_size',self.this_key] = self.arrow_text_size
            if self.arrow_color != 'None':
                this_data.loc['arrow_color',self.this_key] = self.arrow_color
            if self.arrow_style != 'None':
                this_data.loc['arrow_style',self.this_key] = self.arrow_style

        this_data.loc['x_fun',self.this_key] = x_funs[self.x_function]

        this_data.loc['x_scale',self.this_key] = self.x_scale

        if '_vary' in data_plot.plot_data.loc['type',self.this_key]:
            if self.test_pick_dim:
                this_key = []
                for ic in range(1,self.vary_len+1):
                    if getattr(self,'test_dim_'+str(ic)):
                        this_key.append(str(getattr(self,'dim_'+str(ic))))
                    else:
                        this_key.append(slice(None))
                this_key = tuple(this_key)
            else:
                this_key = []
                for ic in range(1,self.vary_len+1):
                    if getattr(self,'test_dim_'+str(ic)):
                        this_key.append(str(getattr(self,'dim_'+str(ic))))
                if len(this_key) == 0:
                    this_key = slice(None)
            this_data.loc['key_select',self.this_key] = this_key
        data_plot.ImportData(this_data)
        print(data_plot)
        data_plot.PlotAll()
        self._Load_fired()

class InfoFrame( ta.HasTraits ):
    save_file = ta.File('')
    legend_location = ta.Str('')
    legend_location = ta.Enum(ql.LegLocList)
    legend_number_of_columns = ta.Int(jpl.legcol)
    legend_font_size = ta.Int(jpl.legsize)
    legend_alpha = ta.Float(jpl.legalpha)
    title = ta.Str('')
    x_label = ta.Str('')
    y_label = ta.Str('')
    title_pad = ta.Int(jpl.tpad)
    title_size = ta.Int(jpl.tsize)
    x_label_size = ta.Int(jpl.xsize)
    y_label_size = ta.Int(jpl.ysize)
    grid = ta.Bool(jpl.def_grid)
    dot_size = ta.Int(jpl.dotsize)
    line_width = ta.Int(jpl.linewidth)
    error_cap_len = ta.Int(jpl.errorbarlen)
    x_axis_min = ta.Str('None')
    x_axis_max = ta.Str('None')
    y_axis_min = ta.Str('None')
    y_axis_max = ta.Str('None')
    y_axis_max = ta.Str('None')
    custom_y_labels = ta.Bool(False)
    flip_vert = ta.Bool(False)
    mirror_vert = ta.Bool(False)
    y_axis_labels = ta.List([''])
    y_scale = ta.Enum(['None',"linear", "log", "symlog", "logit"])
    x_scale = ta.Enum(['None',"linear", "log", "symlog", "logit"])
    x_zero_line = ta.Bool(False)
    x_zero_line_val = ta.Float(0)
    y_zero_line = ta.Bool(False)
    y_zero_line_val = ta.Float(0)
    Load = ta.Button()
    Apply = ta.Button()
    Reset_Axies = ta.Button()

    tick_size = ta.Int(jpl.def_tick_size)
    xTick_min = ta.Float(0)
    xTick_max = ta.Float(0)
    xTick_inc = ta.Float(0)
    yTick_min = ta.Float(0)
    yTick_max = ta.Float(0)
    yTick_inc = ta.Float(0)


    view = tua.View(
        tua.Group(
        tua.Item('save_file'),
        tua.Item('legend_location'),
        tua.Item('legend_number_of_columns'),
        tua.Item('legend_font_size'),
        tua.Item('legend_alpha'),
        tua.Item('title'),
        tua.Item('x_label'),
        tua.Item('y_label'),
        tua.Item('title_pad'),
        tua.Item('title_size'),
        tua.Item('x_label_size'),
        tua.Item('y_label_size'),
        tua.Item('grid'),
        tua.Item('x_zero_line'),
        tua.Item('x_zero_line_val'),
        tua.Item('y_zero_line'),
        tua.Item('y_zero_line_val'),
        tua.Item('dot_size'),
        tua.Item('line_width'),
        tua.Item('error_cap_len'),
        tua.Item('x_axis_min'),
        tua.Item('x_axis_max'),
        tua.Item('x_scale'),
        tua.Item('y_axis_min'),
        tua.Item('y_axis_max'),
        tua.Item('y_scale'),
        tua.Item('Apply', show_label=False),
        tua.Item('Reset_Axies', show_label=False),
        label='Main'
        ),tua.Group(
        tua.Item('tick_size'),
        tua.Item('xTick_min'),
        tua.Item('xTick_max'),
        tua.Item('xTick_inc'),
        tua.Item('yTick_min'),
        tua.Item('yTick_max'),
        tua.Item('yTick_inc'),
        tua.Item('custom_y_labels'),
        tua.Item('flip_vert',enabled_when='custom_y_labels'),
        tua.Item('mirror_vert',enabled_when='custom_y_labels'),
        tua.Item('y_axis_labels',enabled_when='custom_y_labels'),
        tua.Item('Apply', show_label=False),
        tua.Item('Reset_Axies', show_label=False),
        label='Ticks'
        ),
        buttons=['OK'],
        resizable=True
    )

    def __init__(self):
        global data_plot
        self._Load_fired()
        self.save_file = data_plot.plot_info['save_file']
        self._Apply_fired()

    def _Load_fired(self):
        if data_plot is None:
            raise EnvironmentError('plot data has not been loaded')
        if data_plot.plot_info['leg_ncol'] == 'best':
            self.legend_number_of_columns = 1
        else:
            self.legend_number_of_columns = int(data_plot.plot_info['leg_ncol'])
        self.legend_font_size = int(data_plot.plot_info['leg_fontsize'])
        self.legend_alpha = float(data_plot.plot_info['leg_alpha'])

        if 'xTick_min' in data_plot.plot_info:
            self.xTick_min = data_plot.plot_info['xTick_min']
        else:
            self.xTick_max = 0
        if 'xTick_max' in data_plot.plot_info:
            self.xTick_max = data_plot.plot_info['xTick_max']
        else:
            self.xTick_inc = 0
        if 'xTick_inc' in data_plot.plot_info:
            self.xTick_inc = data_plot.plot_info['xTick_inc']
        else:
            self.xTick_inc = 0
        if 'yTick_min' in data_plot.plot_info:
            self.yTick_min = data_plot.plot_info['yTick_min']
        else:
            self.yTick_min = 0
        if 'yTick_max' in data_plot.plot_info:
            self.yTick_max = data_plot.plot_info['yTick_max']
        else:
            self.yTick_max = 0
        if 'yTick_inc' in data_plot.plot_info:
            self.yTick_inc = data_plot.plot_info['yTick_inc']
        else:
            self.yTick_inc = 0

        if 'mirror_vert' in data_plot.plot_info:
            self.mirror_vert = data_plot.plot_info['mirror_vert']
        else:
            self.mirror_vert = False

        if 'flip_vert' in data_plot.plot_info:
            self.flip_vert = data_plot.plot_info['flip_vert']
        else:
            self.flip_vert = False

        if 'custom_yTicks' in data_plot.plot_info:
            self.custom_y_labels = True
            self.y_axis_labels = data_plot.plot_info['custom_yTicks']
        else:
            self.custom_y_labels = False


        if 'title' in data_plot.plot_info:
            self.title = data_plot.plot_info['title']
        else:
            self.title = ''
        if 'xlabel' in data_plot.plot_info:
            self.x_label = data_plot.plot_info['xlabel']
        else:
            self.x_label = ''
        if 'ylabel' in data_plot.plot_info:
            self.y_label = data_plot.plot_info['ylabel']
        else:
            self.y_label = ''

        if 'title_dict' in data_plot.plot_info and 'fontsize' in list(data_plot.plot_info['title_dict'].keys()):
            self.title_size = data_plot.plot_info['title_dict']['fontsize']
        else:
            self.title_size = jpl.tsize
        if 'title_pad' in data_plot.plot_info:
            if data_plot.plot_info['title_pad'] is None:
                self.title_pad = 0
            else:
                self.title_pad = data_plot.plot_info['title_pad']
        else:
            self.title_pad = jpl.tpad
        if 'xlabel_dict' in data_plot.plot_info and 'fontsize' in list(data_plot.plot_info['title_dict'].keys()):
            self.x_label_size = data_plot.plot_info['xlabel_dict']['fontsize']
        else:
            self.x_label_size = jpl.xsize
        if 'ylabel_dict' in data_plot.plot_info and 'fontsize' in list(data_plot.plot_info['title_dict'].keys()):
            self.y_label_size = data_plot.plot_info['ylabel_dict']['fontsize']
        else:
            self.y_label_size = jpl.ysize
        if 'xlims' in data_plot.plot_info:
            self.x_axis_min = str(data_plot.plot_info['xlims'][0])
            self.x_axis_max = str(data_plot.plot_info['xlims'][1])
        else:
            self.x_axis_min = 'None'
            self.x_axis_max = 'None'
        if 'ylims' in data_plot.plot_info:
            self.y_axis_min = str(data_plot.plot_info['ylims'][0])
            self.y_axis_max = str(data_plot.plot_info['ylims'][1])
        else:
            self.y_axis_min = 'None'
            self.y_axis_max = 'None'
        if 'x_plot_scale' in data_plot.plot_info:
            self.x_scale = data_plot.plot_info['x_plot_scale']
        if 'y_plot_scale' in data_plot.plot_info:
            self.y_scale = data_plot.plot_info['y_plot_scale']

    def _Reset_Axies_fired(self):
        self.x_axis_min = 'None'
        self.x_axis_max = 'None'
        self.y_axis_min = 'None'
        self.y_axis_max = 'None'

    def _Apply_fired(self):
        # global customleg,customleglist
        # global wtitle,wxlab,wrlab
        # global legloc
        global data_plot
        plot_info = pad.Series()
        plot_info['save_file'] = self.save_file.replace('.pdf','')+'.pdf'
        plot_info['leg_loc'] = self.legend_location
        plot_info['leg_ncol'] = self.legend_number_of_columns
        plot_info['leg_fontsize'] = self.legend_font_size
        plot_info['leg_alpha'] = self.legend_alpha
        plot_info['title'] = self.title
        plot_info['title_dict'] = {'fontsize':self.title_size}
        if int(self.title_pad) == 0:
            plot_info['title_pad'] = None
        else:
            plot_info['title_pad'] = self.title_pad
        plot_info['ylabel'] = self.y_label
        plot_info['ylabel_dict'] = {'fontsize':self.y_label_size}
        plot_info['xlabel'] = self.x_label
        plot_info['xlabel_dict'] = {'fontsize':self.x_label_size}
        plot_info['x_zero_line'] = self.x_zero_line
        plot_info['x_zero_line_val'] = self.x_zero_line_val
        plot_info['y_zero_line'] = self.y_zero_line
        plot_info['y_zero_line_val'] = self.y_zero_line_val
        if self.x_axis_min == 'None':
            xmin = None
        else:
            xmin = float(self.x_axis_min)
        if self.x_axis_max == 'None':
            xmax = None
        else:
            xmax = float(self.x_axis_max)
        if self.y_axis_min == 'None':
            ymin = None
        else:
            ymin = float(self.y_axis_min)
        if self.y_axis_max == 'None':
            ymax = None
        else:
            ymax = float(self.y_axis_max)
        plot_info['xlims'] = [xmin,xmax]

        plot_info['ylims'] = [ymin,ymax]
        plot_info['x_plot_scale'] = self.x_scale
        plot_info['y_plot_scale'] = self.y_scale
        if self.xTick_inc > 0:
            plot_info['do_xTicks'] = True
            plot_info['xTick_min'] = self.xTick_min
            plot_info['xTick_max'] = self.xTick_max
            plot_info['xTick_inc'] = self.xTick_inc
        if self.yTick_inc > 0:
            plot_info['do_yTicks'] = True
            plot_info['yTick_min'] = self.yTick_min
            plot_info['yTick_max'] = self.yTick_max
            plot_info['yTick_inc'] = self.yTick_inc
        if self.custom_y_labels:
            plot_info['mirror_vert'] =  self.mirror_vert
            plot_info['flip_vert'] =    self.flip_vert
            plot_info['custom_yTicks'] = self.y_axis_labels
        data_plot.ImportInfo(plot_info)
        this_rc = jpl.params
        this_rc['lines.markersize'] = self.dot_size
        this_rc['axes.grid'] = self.grid
        this_rc['lines.markersize'] = self.dot_size
        this_rc['lines.linewidth'] = self.line_width
        this_rc['errorbar.capsize'] = self.error_cap_len
        this_rc['xtick.labelsize'] = self.tick_size
        this_rc['ytick.labelsize'] = self.tick_size
        # jpl.capwidth = self.line_width
        pl.rcParams.update(this_rc)
        data_plot.PrintData()
        data_plot.PlotAll()
        self._Load_fired()

def Load_Prev_File():
    global hold_file
    if os.path.isfile(hold_file):
        with open(hold_file,'rb') as f:
            file_name = pickle.load(f)
    else:
        print('Warning: ',hold_file,' not found, setting to default')
        file_name  = '/home/jackdra/PHD/CHROMA/TestVar/scratch/resultsSCfgs/HPCCRES/graphs/AlphaFull/Set1_Alpha_tsum.pdf'
    return file_name

def Save_Prev_File(file_name):
    global hold_file
    with open(hold_file,'wb') as f:
        # print 'DEBUG'
        # print 'printing ', file_name,' to ', hold_file
        pickle.dump(file_name,f)


class BeginFrame( ta.HasTraits ):
    """
    starting frame, to select what type of thing you want to do
    options are:
    """

    file_name = ta.File(Load_Prev_File(),filter=['*.pdf'])

    Load_Plot = ta.Button()
    Merge_Into_Plot = ta.Button()
    is_loaded = ta.Bool(False)
    is_not_loaded = ta.Bool(True)
    force_py2 = ta.Bool(False)


    view = tua.View(
        tua.HSplit(
        tua.Item('file_name',style='custom',springy=True),
        tua.VGroup(
        # tua.Item('Load_Plot', show_label=False,enabled_when='is_not_loaded'),
        tua.Item('file_name',springy=True),
        tua.Item('force_py2'),
        tua.Item('Load_Plot', show_label=False),
        tua.Item('Merge_Into_Plot', show_label=False,enabled_when='is_loaded'),
        ),
        ),
        resizable=True,
        height = 1000,
        width = 1500
    )

    def _Load_Plot_fired(self):
        Save_Prev_File(str(self.file_name))
        global data_plot
        if hasattr(data_plot,'ClearFigure'):
            data_plot.ClearFigure()
        plot_info = pad.Series()
        plot_info['save_file'] = self.file_name.replace('.pdf','')+'.pdf'
        data_plot = jpl.Plotting(plot_info=plot_info)
        data_plot.LoadPickle(DefWipe=False,force_py2=self.force_py2)
        self.Info_window = InfoFrame()
        self.Data_window = DataFrame()
        self.Info_window.configure_traits()
        self.Data_window.configure_traits()
        data_plot.ShowFigure()
        data_plot.ShowFigure()
        self.is_loaded = True
        self.is_not_loaded = False
        # Save_Prev_File(str(self.file_name))

    def _Merge_Into_Plot_fired(self):
        global data_plot2,data_plot
        plot_info = pad.Series()
        plot_info['save_file'] = self.file_name.replace('.pdf','')+'.pdf'
        data_plot2 = jpl.Plotting(plot_info=plot_info)
        data_plot2.LoadPickle(DefWipe=False,force_py2=self.force_py2)
        for iname,idata in data_plot2.plot_data.items():
            idata.name = iname
            idata['plot_this'] = False
            data_plot.AppendData(idata)
        self.Info_window._Load_fired()
        self.Data_window.key_list = list(data_plot.plot_data.keys())
        self.Data_window._Load_fired()
        data_plot.ShowFigure()
        data_plot.ShowFigure()
        Save_Prev_File(str(self.file_name))





if  __name__ == "__main__":
    Start = BeginFrame()
    Start.configure_traits()
