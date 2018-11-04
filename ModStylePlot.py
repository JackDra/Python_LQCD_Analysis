#!/usr/bin/env python

from traits.etsconfig.api import ETSConfig
ETSConfig.toolkit = 'qt'

import traits.api as ta
import traitsui.api as tua
import PlotData as jpl
import QuantityLists as ql
import pandas as pad
import os

# import pickle as pickle
import glob
import matplotlib.pyplot as pl

global files_selected,select_files

def GetAllPDF(this_dir):
    return glob.glob(this_dir+'/**/*.pdf', recursive=True)


def UpdateFileList(plot_info,rc_params,file_list,window_size=None):
    for ifile in file_list:
        this_plot = jpl.Plotting(plot_info = {'save_file':ifile})
        if not os.path.isfile(this_plot.PickleFile):
            print('FNF:', this_plot.PickleFile)
            pass
        if window_size is not None:
            this_plot.LoadPickle(DefWipe=False,ForceWindowSize=window_size)
        else:
            this_plot.LoadPickle(DefWipe=False)

        this_plot.UpdateInfo(plot_info)
        pl.rcParams.update(rc_params)
        this_plot.PlotAll()



class BeginFrame( ta.HasTraits):

    Start = ta.Button()


    view = tua.View(
        tua.Item('Start',show_label=False)
        )

    def _Start_fired(self):
        global files_selected,select_files
        files_selected = FileSelectedFrame()
        select_files = FileFrame()
        atributes = PlotAtributes()
        files_selected.configure_traits()
        select_files.configure_traits()
        atributes.configure_traits()

class PlotAtributes( ta.HasTraits ):
    """
        Frame for selecting the plot attributes to change
    """
    window_size_x = ta.Float(12.)
    window_size_y = ta.Float(12.275)
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
    x_zero_line = ta.Bool(False)
    y_zero_line = ta.Bool(False)
    Apply = ta.Button()
    Reset_Axies = ta.Button()

    tick_size = ta.Int(jpl.def_tick_size)
    xTick_min = ta.Float(0)
    xTick_max = ta.Float(0)
    xTick_inc = ta.Float(0)
    yTick_min = ta.Float(0)
    yTick_max = ta.Float(0)
    yTick_inc = ta.Float(0)


    Do_window_size = ta.Bool(False)
    Do_legend_location = ta.Bool(False)
    Do_legend_location = ta.Bool(False)
    Do_legend_number_of_columns = ta.Bool(False)
    Do_legend_font_size = ta.Bool(False)
    Do_legend_alpha = ta.Bool(False)
    Do_title = ta.Bool(False)
    Do_x_label = ta.Bool(False)
    Do_y_label = ta.Bool(False)
    Do_title_pad = ta.Bool(False)
    Do_title_size = ta.Bool(False)
    Do_x_label_size = ta.Bool(False)
    Do_y_label_size = ta.Bool(False)
    Do_grid = ta.Bool(False)
    Do_dot_size = ta.Bool(False)
    Do_line_width = ta.Bool(False)
    Do_error_cap_len = ta.Bool(False)
    Do_xlims = ta.Bool(False)
    Do_x_axis_min = ta.Bool(False)
    Do_x_axis_max = ta.Bool(False)
    Do_ylims = ta.Bool(False)
    Do_y_axis_min = ta.Bool(False)
    Do_y_axis_max = ta.Bool(False)
    Do_x_zero_line = ta.Bool(False)
    Do_y_zero_line = ta.Bool(False)

    Do_tick_size = ta.Bool(False)
    Do_xTick = ta.Bool(False)
    Do_yTick = ta.Bool(False)


    view = tua.View(
        tua.Group(
        tua.Item('window_size_x'),
        tua.Item('window_size_y'),
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
        tua.Item('y_zero_line'),
        tua.Item('dot_size'),
        tua.Item('line_width'),
        tua.Item('error_cap_len'),
        tua.Item('x_axis_min'),
        tua.Item('x_axis_max'),
        tua.Item('y_axis_min'),
        tua.Item('y_axis_max'),
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
        ),tua.Group(
        tua.Item('Do_window_size'),
        tua.Item('Do_legend_location'),
        tua.Item('Do_legend_number_of_columns'),
        tua.Item('Do_legend_font_size'),
        tua.Item('Do_legend_alpha'),
        tua.Item('Do_title'),
        tua.Item('Do_x_label'),
        tua.Item('Do_y_label'),
        tua.Item('Do_title_pad'),
        tua.Item('Do_title_size'),
        tua.Item('Do_x_label_size'),
        tua.Item('Do_y_label_size'),
        tua.Item('Do_grid'),
        tua.Item('Do_x_zero_line'),
        tua.Item('Do_y_zero_line'),
        tua.Item('Do_dot_size'),
        tua.Item('Do_line_width'),
        tua.Item('Do_error_cap_len'),
        tua.Item('Do_xlims'),
        tua.Item('Do_x_axis_min'),
        tua.Item('Do_x_axis_max'),
        tua.Item('Do_ylims'),
        tua.Item('Do_y_axis_min'),
        tua.Item('Do_y_axis_max'),
        tua.Item('Do_tick_size'),
        tua.Item('Do_xTick'),
        tua.Item('Do_yTick'),
        tua.Item('Apply', show_label=False),
        tua.Item('Reset_Axies', show_label=False),
        label='Selection'
        ),
        buttons=['OK'],
        resizable=True
    )

    def _Do_xlims_changed(self):
        if not self.Do_xlims_changed:
            self.Do_x_axis_max=False
            self.Do_x_axis_min=False
    def _Do_ylims_changed(self):
        if not self.Do_xlims_changed:
            self.Do_y_axis_max=False
            self.Do_y_axis_min=False


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
        if self.Do_legend_location:
            plot_info['leg_loc'] = self.legend_location
        if self.Do_legend_number_of_columns:
            plot_info['leg_ncol'] = self.legend_number_of_columns
        if self.Do_legend_font_size:
            plot_info['leg_fontsize'] = self.legend_font_size
        if self.Do_legend_alpha:
            plot_info['leg_alpha'] = self.legend_alpha
        if self.Do_title:
            plot_info['title'] = self.title
        if self.Do_title_size:
            plot_info['title_dict'] = {'fontsize':self.title_size}
        if self.Do_title_pad:
            if int(self.title_pad) == 0:
                plot_info['title_pad'] = None
            else:
                plot_info['title_pad'] = self.title_pad
        if self.Do_y_label:
            plot_info['ylabel'] = self.y_label
        if self.Do_y_label_size:
            plot_info['ylabel_dict'] = {'fontsize':self.y_label_size}
        if self.Do_x_label:
            plot_info['xlabel'] = self.x_label
        if self.Do_x_label_size:
            plot_info['xlabel_dict'] = {'fontsize':self.x_label_size}
        if self.Do_x_zero_line:
            plot_info['x_zero_line'] = self.x_zero_line
        if self.Do_y_zero_line:
            plot_info['y_zero_line'] = self.y_zero_line
        if self.Do_x_axis_min:
            if self.x_axis_min == 'None':
                xmin = None
            else:
                xmin = float(self.x_axis_min)
        if self.Do_x_axis_max:
            if self.x_axis_max == 'None':
                xmax = None
            else:
                xmax = float(self.x_axis_max)
        if self.Do_y_axis_min:
            if self.y_axis_min == 'None':
                ymin = None
            else:
                ymin = float(self.y_axis_min)
        if self.Do_y_axis_max:
            if self.y_axis_max == 'None':
                ymax = None
            else:
                ymax = float(self.y_axis_max)
        if self.Do_xlims:
            plot_info['xlims'] = [xmin,xmax]


        if self.Do_ylims:
            plot_info['ylims'] = [ymin,ymax]
        if self.Do_xTick:
            if self.xTick_inc > 0:
                plot_info['do_xTicks'] = True
                plot_info['xTick_min'] = self.xTick_min
                plot_info['xTick_max'] = self.xTick_max
                plot_info['xTick_inc'] = self.xTick_inc
        if self.Do_yTick:
            if self.yTick_inc > 0:
                plot_info['do_yTicks'] = True
                plot_info['yTick_min'] = self.yTick_min
                plot_info['yTick_max'] = self.yTick_max
                plot_info['yTick_inc'] = self.yTick_inc
        this_rc = jpl.params
        if self.Do_dot_size:
            this_rc['lines.markersize'] = self.dot_size
        if self.Do_grid:
            this_rc['axes.grid'] = self.grid
        if self.Do_dot_size:
            this_rc['lines.markersize'] = self.dot_size
        if self.Do_line_width:
            this_rc['lines.linewidth'] = self.line_width
        if self.Do_error_cap_len:
            this_rc['errorbar.capsize'] = self.error_cap_len
        if self.Do_tick_size:
            this_rc['xtick.labelsize'] = self.tick_size
            this_rc['ytick.labelsize'] = self.tick_size
        # if self.Do_line_width:
        #     jpl.capwidth = self.line_width
        # data_plot.UpdateInfo(plot_info)
        # data_plot.PrintData()
        global files_selected
        if self.Do_window_size:
            UpdateFileList(plot_info,this_rc,files_selected.file_list,window_size=[self.window_size_x,self.window_size_y])
        else:
            UpdateFileList(plot_info,this_rc,files_selected.file_list)




class FileFrame( ta.HasTraits ):
    """
    Frame for file selecting
    """

    file_name = ta.File('./',filter=['*.pdf'])
    file_directory = ta.Directory('./')

    Add_File = ta.Button()
    Add_Folder = ta.Button()
    # Undo_Add = ta.Button()


    view = tua.View(
        tua.HSplit(
        tua.Item('file_directory',style='custom',springy=True),
        tua.Item('file_name',style='custom',springy=True),
        tua.VGroup(
        tua.Item('file_directory',springy=True),
        tua.Item('file_name',springy=True),
        tua.Item('Add_File',show_label=False),
        tua.Item('Add_Folder',show_label=False)
        # tua.Item('Undo_Add',show_label=False),
        )),
        resizable=True,
        height = 1000,
        width = 1500
    )


    def _file_name_changed(self):
        self.file_directory = '/'.join(self.file_name.split('/')[:-1])+'/'

    def _file_directory_changed(self):
        file_list = GetAllPDF(self.file_directory)
        if len(file_list) > 0:
            self.file_name = GetAllPDF(self.file_directory)[0]

    def _Add_File_fired(self):
        global files_selected
        files_selected.file_list.append(self.file_name)

    def _Add_Folder_fired(self):
        global files_selected
        files_selected.file_list += GetAllPDF(self.file_directory)



class FileSelectedFrame( ta.HasTraits ):
    """
    Frame for current files selected
    """


    file_list = ta.List(ta.Str,[])

    Add_File = ta.Button()
    Add_Folder = ta.Button()
    Undo_Add = ta.Button()


    view = tua.View(
        tua.Item('file_list'),
        tua.Item('Add_File',show_label=False),
        tua.Item('Add_Folder',show_label=False),
        tua.Item('Undo_Add',show_label=False),
        resizable=True
    )

    def _Add_File_fired(self):
        global select_files
        self.file_list.append(select_files.file_name)

    def _Add_Folder_fired(self):
        global select_files
        self.file_list += GetAllPDF(select_files.file_directory)

    def _Undo_Add_fired(self):
        del self.file_list[-1]

if  __name__ == "__main__":
    beginframe = BeginFrame()
    beginframe.configure_traits()
