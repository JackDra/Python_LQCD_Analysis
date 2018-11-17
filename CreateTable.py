#!/usr/bin/env python

import traits.api as ta
import traitsui.api as tua
from FileIO import ReadXml
from XmlFormatting import MakeValAndErr
from xmltodict import unparse
from collections import OrderedDict
import pandas as pa
import numpy as np
import pickle
# from copy import deepcopy
import os

xml_data = None
latex_data = None
hold_file = './prev_file_table.py3p'

class PlotFrame( ta.HasTraits ):
    # key_list = range(1,21)
    # this_key = ta.Enum(key_list)
    output_file = ta.Str('')
    this_data_type = ta.Enum(['data','fit_parameters'])
    transpose_table = ta.Bool(False)
    include_chi = ta.Bool(True)
    fmt_latex = ta.Bool(True)
    Load = ta.Button()
    Show = ta.Button()
    Show_latex = ta.Button()
    Apply = ta.Button()
    #
    view = tua.View(
        tua.Item('output_file'),
        tua.Item('this_data_type'),
        tua.Item('transpose_table'),
        tua.Item('include_chi'),
        tua.Item('fmt_latex'),
        tua.Item('Load', show_label=False),
        tua.Item('Show', show_label=False),
        tua.Item('Show_latex', show_label=False),
        tua.Item('Apply', show_label=False),
        buttons=['OK'],
        resizable=True
    )
    #
    def _Load_fired(self):
        global xml_data
        if xml_data is None:
            raise EnvironmentError('plot data has not been loaded')
        # self.key_list = xml_data.plot_data.keys()
        self.output_file = xml_data['Results']['save_file']


    def _Show_fired(self):
        global xml_data
        print(unparse(xml_data,pretty=True))

    def _Show_latex_fired(self):
        global latex_data
        if latex_data is None:
            print('latex data not generated yet')
        else:
            print(latex_data.to_string())


    def _Apply_fired(self):
        global xml_data,latex_data
        this_data = xml_data['Results']
        table_data = pa.DataFrame()
        for cfg_key,cfg_data in this_data.items():
            if cfg_key in ['window_size_x','window_size_y','Info','save_file']: continue
            this_data = pa.Series()
            if self.this_data_type in cfg_data.keys():
                if self.this_data_type == 'fit_parameters' and self.include_chi \
                and 'chi_pow_2_pdf' in cfg_data.keys():
                    this_data[r'$\chi^{2}_{pdf}$'] = '$'+MakeValAndErr(cfg_data['chi_pow_2_pdf']['Avg'],cfg_data['chi_pow_2_pdf']['Std'],Dec=2,latex=self.fmt_latex)+'$'
                for data_key,data_data in cfg_data[self.this_data_type].items():
                    if isinstance(data_data,dict):
                        this_data[data_key] = '$'+MakeValAndErr(data_data['Avg'],data_data['Std'],Dec=2,latex=self.fmt_latex)+'$'
                    else:
                        if isinstance(data_data,(list,tuple,np.ndarray)):
                            if len(data_data) == 2:
                                this_data[data_key] = '$'+MakeValAndErr(data_data[0],data_data[1],Dec=2,latex=self.fmt_latex)+'$'
                            elif len(data_data) == 1:
                                this_data[data_key] = f'${data_data[0]:.3f}$'
                        else:
                            this_data[data_key] = f'${data_data[0]:.3f}$'
            if len(this_data) > 0:
                table_data[cfg_key] = this_data
        latex_data = table_data
        if self.transpose_table:
            latex_data = latex_data.transpose()
        with open(xml_data['Results']['save_file'],'w') as f:
            format_output = FormatLatex(latex_data.to_latex(escape=False))
            f.write(format_output)


class AlphaFrame( ta.HasTraits ):
    # key_list = range(1,21)
    # this_key = ta.Enum(key_list)
    output_file = ta.Str('')
    this_momentum = ta.Str('p000')
    this_flow_time = ta.Str('t_f6.01')
    this_parameter = ta.Str(r'\alpha')
    do_chi = ta.Bool(True)
    do_chi_err = ta.Bool(True)
    Load = ta.Button()
    Show = ta.Button()
    Show_latex = ta.Button()
    Apply = ta.Button()
    #
    view = tua.View(
        tua.Item('output_file'),
        tua.Item('this_momentum'),
        tua.Item('this_flow_time'),
        tua.Item('this_parameter'),
        tua.Item('do_chi'),
        tua.Item('do_chi_err',enabled_when='do_chi'),
        tua.Item('Load', show_label=False),
        tua.Item('Show', show_label=False),
        tua.Item('Show_latex', show_label=False),
        tua.Item('Apply', show_label=False),
        buttons=['OK'],
        resizable=True
    )
    #
    def _Load_fired(self):
        global xml_data
        if xml_data is None:
            raise EnvironmentError('plot data has not been loaded')
        # self.key_list = xml_data.plot_data.keys()
        self.output_file = xml_data['Results']['save_file']


    def _Show_fired(self):
        global xml_data
        print(unparse(xml_data,pretty=True))

    def _Show_latex_fired(self):
        global latex_data
        if latex_data is None:
            print('latex data not generated yet')
        else:
            print(latex_data)


    def _Apply_fired(self):
        global xml_data,latex_data
        this_data = xml_data['Results']
        if 'Fits' not in list(this_data.keys()):
            print('Fits not in keys')
            print(list(this_data.keys()))
            return
        this_data = this_data['Fits']
        if self.this_momentum not in list(this_data.keys()):
            print(self.this_momentum , ' not in momentum keys ')
            print(list(this_data.keys()))
            return
        this_data = this_data[self.this_momentum]
        if self.this_flow_time not in list(this_data.keys()):
            print(self.this_flow_time , ' not in flow time keys ')
            print(list(this_data.keys()))
            return
        this_data = this_data[self.this_flow_time]
        vall,ilist,chil,slist1,slist2 = [],[],[],[],[]
        for ikey,idata in this_data.items():
            this_p = self.this_parameter.replace('\\alpha','alpha')
            if this_p in list(idata.keys()):
                if isinstance(idata[this_p],dict) or isinstance(idata[this_p],OrderedDict):
                    vall.append(MakeValAndErr(idata[this_p]['Avg'],idata[this_p]['Std'],latex=True))
                else:
                    vall.append(idata[this_p])
                if 'Chi_pow_2_pdf' in list(idata.keys()):
                    if isinstance(idata['Chi_pow_2_pdf'],dict) or isinstance(idata['Chi_pow_2_pdf'],OrderedDict):
                        chil.append(MakeValAndErr(idata['Chi_pow_2_pdf']['Avg'],idata['Chi_pow_2_pdf']['Std'],latex=True))
                    else:
                        chil.append(idata['Chi_pow_2_pdf'])
                ilist.append(ikey)
                leftval,rightval = list(map(float,ikey[-4:].split('-')))
                slist1.append(1/rightval)
                slist2.append(leftval)
        dump,dump2,ilist,vall,chil =  (list(x) for x in zip(*sorted(zip(slist1,slist2,ilist,vall,chil))))
        latex_data = pa.DataFrame(index=ilist)
        latex_data.loc[:,self.this_parameter+' Avg'] = pa.Series(vall,index=ilist)
        if len(vall) == len(chil) and self.do_chi:
            latex_data.loc[:,self.this_parameter+' chi2pdf'] = pa.Series(chil,index=ilist)
        with open(xml_data['Results']['save_file'],'w') as f:
            format_output = FormatLatex(latex_data.to_latex(escape=False))
            f.write(format_output)


class FlowOpFrame( ta.HasTraits ):
    # key_list = range(1,21)
    # this_key = ta.Enum(key_list)
    output_file = ta.Str('')
    this_flow_time = ta.Str('t_f6.01')
    this_t_sum = ta.Str('ts63')
    this_parameter = ta.Str('A')
    do_chi = ta.Bool(True)
    do_chi_err = ta.Bool(True)
    Load = ta.Button()
    Show = ta.Button()
    Show_latex = ta.Button()
    Apply = ta.Button()
    #
    view = tua.View(
        tua.Item('output_file'),
        tua.Item('this_flow_time'),
        tua.Item('this_t_sum'),
        tua.Item('this_parameter'),
        tua.Item('do_chi'),
        tua.Item('do_chi_err',enabled_when='do_chi'),
        tua.Item('Load', show_label=False),
        tua.Item('Show', show_label=False),
        tua.Item('Show_latex', show_label=False),
        tua.Item('Apply', show_label=False),
        buttons=['OK'],
        resizable=True
    )
    #
    def _Load_fired(self):
        global xml_data
        if xml_data is None:
            raise EnvironmentError('plot data has not been loaded')
        # self.key_list = xml_data.plot_data.keys()
        self.output_file = xml_data['Results']['save_file']


    def _Show_fired(self):
        global xml_data
        print(unparse(xml_data,pretty=True))

    def _Show_latex_fired(self):
        global latex_data
        if latex_data is None:
            print('latex data not generated yet')
        else:
            print(latex_data)


    def _Apply_fired(self):
        global xml_data,latex_data
        this_data = xml_data['Results']
        vall,ilist,chil,slist1,slist2 = [],[],[],[],[]
        if 'Chi_boot' not in list(this_data.keys()):
            print('Integrated FlowOp result not in keys')
        else:
            this_data = this_data['Chi_boot']
            if self.this_flow_time not in list(this_data.keys()):
                print(self.this_flow_time , ' not in flow time keys ')
                print(list(this_data.keys()))
                return
            this_data = this_data[self.this_flow_time]
            if isinstance(this_data,dict) or isinstance(this_data,OrderedDict):
                vall.append(MakeValAndErr(this_data['Avg'],this_data['Std'],latex=True))
            else:
                vall.append(this_data)
            chil.append('N.A.')
            ilist.append('chiint')

        this_data = xml_data['Results']
        if 'Prop_Fit' not in list(this_data.keys()):
            print('Fits not in keys')
            print(list(this_data.keys()))
            return
        this_data = this_data['Prop_Fit']
        if self.this_flow_time not in list(this_data.keys()):
            print(self.this_flow_time , ' not in flow time keys ')
            print(list(this_data.keys()))
            return
        this_data = this_data[self.this_flow_time]
        if self.this_t_sum not in list(this_data.keys()):
            print(self.this_t_sum , ' not in sum source time keys ')
            print(list(this_data.keys()))
            return
        this_data = this_data[self.this_t_sum]
        for ikey,idata in this_data.items():
            this_p = self.this_parameter
            if this_p in list(idata.keys()):
                if isinstance(idata[this_p],dict) or isinstance(idata[this_p],OrderedDict):
                    vall.append(MakeValAndErr(idata[this_p]['Avg'],idata[this_p]['Std'],latex=True))
                else:
                    vall.append(idata[this_p])
                if 'Chi_pow_2_pdf' in list(idata.keys()):
                    if isinstance(idata['Chi_pow_2_pdf'],dict) or isinstance(idata['Chi_pow_2_pdf'],OrderedDict):
                        chil.append(MakeValAndErr(idata['Chi_pow_2_pdf']['Avg'],idata['Chi_pow_2_pdf']['Std'],latex=True))
                    else:
                        chil.append(idata['Chi_pow_2_pdf'])
                ilist.append(ikey)
                leftval,rightval = list(map(float,ikey.replace('fitr','').split('-')))
                slist1.append(1/rightval)
                slist2.append(leftval)
        dump,dump2,ilisttemp,valltemp,chiltemp =  (list(x) for x in zip(*sorted(zip(slist1,slist2,ilist[1:],vall[1:],chil[1:]))))
        ilist = [ilist[0]] + ilisttemp
        vall = [vall[0]] + valltemp
        chil = [chil[0]] + chiltemp
        latex_data = pa.DataFrame(index=ilist)
        latex_data.loc[:,self.this_parameter+' Avg'] = pa.Series(vall,index=ilist)
        if len(vall) == len(chil) and self.do_chi:
            latex_data.loc[:,self.this_parameter+' chi2pdf'] = pa.Series(chil,index=ilist)
        with open(xml_data['Results']['save_file'],'w') as f:
            format_output = FormatLatex(latex_data.to_latex(escape=False))
            f.write(format_output)


def FormatLatex(output):
    this_out = output.replace(r'\textbackslash ','\\')
    this_out = output.replace(r'\$',r'$')
    this_out = output.replace(r'\textasciicircum \{',r'^')
    this_out = output.replace(r'\}',r'}')
    this_out = this_out.replace(r'\_',r' ')
    this_out = this_out.replace('chi2pdf',r'\chi^{2}_{pdf}')
    this_out = this_out.replace('chiint',r'\chi_{int}')
    return this_out

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

    file_name = ta.File(Load_Prev_File(),filter=['*.xml'])
    Alpha_window = AlphaFrame()
    FlowOp_window = FlowOpFrame()
    Plot_window = PlotFrame()
    Load_FlowOp = ta.Button()
    Load_Alpha = ta.Button()
    Load_Plot = ta.Button()

    view = tua.View(
        '_',
        tua.Item('file_name'),
        tua.Item('Load_FlowOp', show_label=False),
        tua.Item('Load_Alpha', show_label=False),
        tua.Item('Load_Plot', show_label=False),
        buttons=['OK'],
        resizable=True
    )

    def _Load_Alpha_fired(self):
        global xml_data
        xml_data = {}
        xml_data['save_file'] = self.file_name
        if os.path.isfile(self.file_name):
            xml_data = ReadXml(self.file_name)
            xml_data['Results']['save_file'] = self.file_name.replace('.xml','.tex')
            self.Alpha_window.configure_traits()
            Save_Prev_File(str(self.file_name))
        else:
            print('File not found:')
            print(xml_data['save_file'])

    def _Load_Plot_fired(self):
        global xml_data
        xml_data = {}
        xml_data['save_file'] = self.file_name
        if os.path.isfile(self.file_name):
            xml_data = ReadXml(self.file_name)
            xml_data['Results']['save_file'] = self.file_name.replace('.xml','.tex')
            self.Plot_window.configure_traits()
            Save_Prev_File(str(self.file_name))
        else:
            print('File not found:')
            print(xml_data['save_file'])


    def _Load_FlowOp_fired(self):
        global xml_data
        xml_data = {}
        xml_data['save_file'] = self.file_name
        if os.path.isfile(self.file_name):
            xml_data = ReadXml(self.file_name)
            xml_data['Results']['save_file'] = self.file_name.replace('.xml','.tex')
            self.FlowOp_window.configure_traits()
            Save_Prev_File(str(self.file_name))
        else:
            print('File not found:')
            print(xml_data['save_file'])





if  __name__ == "__main__":
    Start = BeginFrame()
    Start.configure_traits()
