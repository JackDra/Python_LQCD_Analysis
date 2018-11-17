#!/usr/bin/env python

from traits.etsconfig.api import ETSConfig
ETSConfig.toolkit = 'qt'

import traits.api as ta
import traitsui.api as tua
import PlotData as jpl
import QuantityLists as ql
import pandas as pad
from copy import deepcopy
from collections import OrderedDict
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

    this_aindex = ta.Int(0)
    animate_dim = ta.Int(0)
    Animate = ta.Button()



    x_range_min = ta.Int(0)
    x_range_max = ta.Int(-1)
    x_increment = ta.Int(1)
    x_scale = ta.Float(1)
    x_function = ta.Enum(list(x_funs.keys()))
    shift = ta.Float(1.0)
    shift_overlap = ta.Float(0.005)

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
        tua.Item('supress_legend'),
        tua.Item('suppress_key_index'),
        tua.Item('pick_dimension',visible_when='test_pick_dim'),
        tua.Item('dim_1',visible_when='test_dim_1'),
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
        label='Arrow'
        ),tua.Group(
        tua.Item('this_key'),
        tua.Item('this_key_2'),
        tua.Item('combine_operation'),
        tua.Item('new_label'),
        tua.Item('Combine', show_label=False),
        label='Combine'
        ),
        buttons=['OK'],
        resizable=True
    )

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

    def _dim_1_fired(self):
        self._Load_fired()
    def _dim_2_fired(self):
        self._Load_fired()
    def _dim_3_fired(self):
        self._Load_fired()
    def _dim_4_fired(self):
        self._Load_fired()
    def _dim_5_fired(self):
        self._Load_fired()
    def _dim_6_fired(self):
        self._Load_fired()
    def _dim_7_fired(self):
        self._Load_fired()
    def _dim_8_fired(self):
        self._Load_fired()


    def _Load_fired(self):
        global data_plot
        self.test_pick_dim = False
        for ic in range(1,8):
            setattr(self,'test_dim_'+str(ic),False)

        if data_plot is None:
            raise EnvironmentError('plot data has not been loaded')
        # self.key_list = data_plot.plot_data.keys()
        if 'plot_this' in data_plot.plot_data.index:
            self.plot_this = data_plot.plot_data.loc['plot_this',self.this_key]
        else:
            self.plot_this = True

        if 'arrow' in data_plot.plot_data.loc['type',self.this_key]:
            if 'x_data' in data_plot.plot_data.index:
                self.arrow_x,self.arrow_x_text = data_plot.plot_data.loc['x_data',self.this_key][:2]
            else:
                self.arrow_x,self.arrow_x_text = 0,1
            if 'y_data' in data_plot.plot_data.index:
                self.arrow_y,self.arrow_y_text = data_plot.plot_data.loc['y_data',self.this_key][:2]
            else:
                self.arrow_y,self.arrow_y_text = 0,1
            if 'arrow_text_size' in data_plot.plot_data.index:
                self.arrow_text_size = data_plot.plot_data.loc['arrow_text_size',self.this_key]
            else:
                self.arrow_text_size = 30
            if 'arrow_color' in data_plot.plot_data.index:
                self.arrow_color = data_plot.plot_data.loc['arrow_color',self.this_key]
            else:
                self.arrow_color = 'black'
            if 'arrow_style' in data_plot.plot_data.index:
                self.arrow_style = data_plot.plot_data.loc['arrow_style',self.this_key]
            else:
                self.arrow_style = 'fancy'

        if 'label' in data_plot.plot_data.index:
            self.label = data_plot.plot_data.loc['label',self.this_key]
        else:
            self.label = 'None'
        if 'phys_axies' in data_plot.plot_data.index:
            self.Physical_Units = data_plot.plot_data.loc['phys_axies',self.this_key]
        else:
            self.Physical_Units = True
        if 'extrap_fade' in data_plot.plot_data.index:
            self.extrap_fade = data_plot.plot_data.loc['extrap_fade',self.this_key]
        else:
            self.extrap_fade = 1.5
        if 'hairline' in data_plot.plot_data.index:
            self.Hairline_plot = data_plot.plot_data.loc['hairline',self.this_key]
        else:
            self.Hairline_plot = True
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
                try:
                    self.fit_custom_xrange = True
                    self.fit_x_min,self.fit_x_max = data_plot.plot_data.loc['xdatarange',self.this_key]
                except:
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
                except:
                    pass
        else:
            self.color = 'blue'

        if self.this_key not in data_plot.plot_data or 'shift' not in data_plot.plot_data.loc[:,self.this_key] or data_plot.plot_data[self.this_key]['shift'] is None:
            self.shift = 0.0
        else:
            self.shift = float(data_plot.plot_data[self.this_key]['shift'])

        if self.this_key not in data_plot.plot_data or 'supress_legend' not in data_plot.plot_data.loc[:,self.this_key] or data_plot.plot_data[self.this_key]['supress_legend'] is None:
            self.supress_legend = False
        else:
            self.supress_legend = bool(data_plot.plot_data[self.this_key]['supress_legend'])



        if self.this_key not in data_plot.plot_data or 'x_range_min' not in data_plot.plot_data.loc[:,self.this_key] or data_plot.plot_data[self.this_key]['x_range_min'] is None:
            self.x_range_min = 0
        else:
            self.x_range_min = int(data_plot.plot_data[self.this_key]['x_range_min'])

        if self.this_key not in data_plot.plot_data or 'x_range_max' not in data_plot.plot_data.loc[:,self.this_key] or data_plot.plot_data[self.this_key]['x_range_max'] is None:
            self.x_range_max = -1
        elif int(data_plot.plot_data[self.this_key]['x_range_max']) != 0:
            self.x_range_max = int(data_plot.plot_data[self.this_key]['x_range_max'])
        else:
            self.x_range_max = -1

        if self.this_key not in data_plot.plot_data or 'x_increment' not in data_plot.plot_data.loc[:,self.this_key] or data_plot.plot_data[self.this_key]['x_increment'] is None:
            self.x_increment = 1
        elif int(data_plot.plot_data[self.this_key]['x_increment']) != 0:
            self.x_increment = int(data_plot.plot_data[self.this_key]['x_increment'])

        if 'symbol' in data_plot.plot_data.index:
            self.symbol = str(data_plot.plot_data.loc['symbol',self.this_key])
        else:
            self.symbol = 'None'
        if 'Median' in data_plot.plot_data.index:
            self.Median = data_plot.plot_data.loc['Median',self.this_key]
        else:
            self.Median = False
        if 'suppress_key_index' in data_plot.plot_data.index:
            self.suppress_key_index = data_plot.plot_data.loc['suppress_key_index',self.this_key]
        else:
            self.suppress_key_index = False
        if 'x_fun' in data_plot.plot_data.index and hasattr(data_plot.plot_data.loc['x_fun',self.this_key],'__name__'):
            if data_plot.plot_data.loc['x_fun',self.this_key].__name__ is not False:
                self.x_function = data_plot.plot_data.loc['x_fun',self.this_key].__name__
        else:
            self.x_function = 'None'
        if 'x_scale' in data_plot.plot_data.index:
            self.x_scale = data_plot.plot_data.loc['x_scale',self.this_key]
        else:
            self.x_scale = 1
        if 'scale' in data_plot.plot_data.index:
            self.scale_y = data_plot.plot_data.loc['scale',self.this_key]
        else:
            self.scale_y = 1.0
        if len(self.ShowEval) == 0:
            pass
        elif 'ShowEval' in data_plot.plot_data.index and isinstance(data_plot.plot_data.loc['ShowEval',self.this_key],(list,tuple,np.ndarray)):
            self.ShowEval = data_plot.plot_data.loc['ShowEval',self.this_key]
        else:
            self.ShowEval = []
        if 'xaxis' in data_plot.plot_data.index:
            self.fit_dim_pick = data_plot.plot_data.loc['xaxis',self.this_key]
        else:
            self.fit_dim_pick = 0

        if 'FPName' in data_plot.plot_data.index:
            self.Force_Param_Name = str(data_plot.plot_data.loc['FPName',self.this_key])
        else:
            self.Force_Param_Name = ''

        if 'FPUnits' in data_plot.plot_data.index:
            self.Force_Param_Units = str(data_plot.plot_data.loc['FPUnits',self.this_key])
        else:
            self.Force_Param_Units = ''



        if '_vary' in data_plot.plot_data.loc['type',self.this_key]:
            this_series = data_plot.plot_data.loc[:,self.this_key]
            if 'fit' in this_series['type']:
                if isinstance(this_series['fit_class'].index,pad.MultiIndex):
                    self.test_pick_dim = False
                    self.vary_len = len(this_series['fit_class'].index.names)
                    for icd,idim in enumerate(this_series['fit_class'].index.names):
                        setattr(self,'test_dim_'+str(icd+1),True)
                        this_list = list(this_series['fit_class'].index.get_level_values(idim))
                        setattr(self,'dim_'+str(icd+1)+'_list',fmt_Range(this_list))
                else:
                    self.test_pick_dim = False
                    self.test_dim_1 = True
                    self.dim_1_list = fmt_Range(list(this_series['fit_class'].index))
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
                        for icd,idim in enumerate(this_series['y_data'].index.names):
                            if idim == str(self.pick_dimension):
                                setattr(self,'test_dim_'+str(icd+1),False)
                            else:
                                setattr(self,'test_dim_'+str(icd+1),True)
                                this_list = list(this_series['y_data'].index.get_level_values(idim))
                                setattr(self,'dim_'+str(icd+1)+'_list',fmt_Range(this_list))
                else:
                    self.test_pick_dim = False
                    self.test_dim_1 = True
                    this_list = list(this_series['y_data'].index)
                    self.dim_1_list = fmt_Range(this_list)
            elif 'boot_data' in list(this_series.keys()) and isinstance(this_series['boot_data'],pad.Series):
                if isinstance(this_series['boot_data'].index,pad.MultiIndex):
                    self.test_pick_dim = False
                    self.vary_len = len(this_series['boot_data'].index.names)
                    for icd,idim in enumerate(this_series['boot_data'].index.names):
                        setattr(self,'test_dim_'+str(icd+1),True)
                        this_list = list(this_series['boot_data'].index.get_level_values(idim))
                        setattr(self,'dim_'+str(icd+1)+'_list',fmt_Range(this_list))
                else:
                    self.test_pick_dim = False
                    self.test_dim_1 = True
                    self.dim_1_list = fmt_Range(list(this_series['boot_data'].index))
            if 'x_data' in list(this_series.keys()) and isinstance(this_series['x_data'],pad.Series):
                if isinstance(this_series['x_data'].index,pad.MultiIndex):
                    self.test_pick_dim = True
                    self.dim_list = this_series['x_data'].index.names
                    self.vary_len = len(this_series['x_data'].index.names)
                    for icd,idim in enumerate(this_series['x_data'].index.names):
                        if idim == str(self.pick_dimension):
                            setattr(self,'test_dim_'+str(icd+1),False)
                        else:
                            setattr(self,'test_dim_'+str(icd+1),True)
                            this_list = list(this_series['x_data'].index.get_level_values(idim))
                            setattr(self,'dim_'+str(icd+1)+'_list',fmt_Range(this_list))
                else:
                    self.test_pick_dim = False
                    self.test_dim_1 = True
                    self.dim_1_list = fmt_Range(list(this_series['x_data'].index))


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
        if self.x_range_min != 'None':
            this_data.loc['x_range_min',self.this_key] = self.x_range_min
        if self.x_range_max != 'None':
            this_data.loc['x_range_max',self.this_key] = self.x_range_max
        if self.x_increment != 'None':
            this_data.loc['x_increment',self.this_key] = self.x_increment
        if self.symbol != 'None':
            this_data.loc['symbol',self.this_key] = self.symbol
        if self.Median != 'None':
            this_data.loc['Median',self.this_key] = self.Median
        if self.scale_y != 'None':
            this_data.loc['scale',self.this_key] = self.scale_y
        if self.suppress_key_index != 'None':
            this_data.loc['suppress_key_index',self.this_key] = self.suppress_key_index
        if len(list(self.ShowEval)) > 0:
            this_data.loc['ShowEval',self.this_key] = list(self.ShowEval)
        else:
            this_data.loc['ShowEval',self.this_key] = None
        if self.fit_dim_pick != 'None':
            this_data.loc['xaxis',self.this_key] = self.fit_dim_pick
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
        data_plot.PrintData()
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
        tua.Item('tick_size'),
        tua.Item('xTick_min'),
        tua.Item('xTick_max'),
        tua.Item('xTick_inc'),
        tua.Item('yTick_min'),
        tua.Item('yTick_max'),
        tua.Item('yTick_inc'),

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


    view = tua.View(
        tua.HSplit(
        tua.Item('file_name',style='custom',springy=True),
        tua.VGroup(
        # tua.Item('Load_Plot', show_label=False,enabled_when='is_not_loaded'),
        tua.Item('file_name',springy=True),
        tua.Item('Load_Plot', show_label=False),
        tua.Item('Merge_Into_Plot', show_label=False,enabled_when='is_loaded'),
        ),
        ),
        resizable=True,
        height = 1000,
        width = 1500
    )

    def _Load_Plot_fired(self):
        global data_plot
        if hasattr(data_plot,'ClearFigure'):
            data_plot.ClearFigure()
        plot_info = pad.Series()
        plot_info['save_file'] = self.file_name.replace('.pdf','')+'.pdf'
        data_plot = jpl.Plotting(plot_info=plot_info)
        data_plot.LoadPickle(DefWipe=False)
        self.Info_window = InfoFrame()
        self.Data_window = DataFrame()
        self.Info_window.configure_traits()
        self.Data_window.configure_traits()
        data_plot.ShowFigure()
        data_plot.ShowFigure()
        self.is_loaded = True
        self.is_not_loaded = False
        Save_Prev_File(str(self.file_name))

    def _Merge_Into_Plot_fired(self):
        global data_plot2,data_plot
        plot_info = pad.Series()
        plot_info['save_file'] = self.file_name.replace('.pdf','')+'.pdf'
        data_plot2 = jpl.Plotting(plot_info=plot_info)
        data_plot2.LoadPickle(DefWipe=False)
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
