#!/usr/bin/env python

# from traits.etsconfig.api import ETSConfig
# ETSConfig.toolkit = 'qt'

# Major library imports
import numpy as np
import PlotData as jpl

# Enthought library imports
from enable.api import Component, ComponentEditor, ColorTrait
from traits.api import HasTraits, Instance, Float, String, List, Enum, Bool, Int
from traitsui.api import Item, HSplit,VGroup, View

# Chaco imports
from chaco.api import ArrayPlotData, HPlotContainer, Plot
from chaco.api import marker_trait
from chaco.tools.api import PanTool, ZoomTool

data_plot = jpl.TestChaco()
simple_col = [ 'blue','red','green','purple', 'gold','black']
def hex_to_rgb(value):
    if not isinstance(value,str):
        return value
    value = value.lstrip('#')
    lv = len(value)
    return tuple(float(int(value[i:i + lv // 3], 16))/255. for i in range(0, lv, lv // 3))

class BeginFrame( HasTraits ):
    """
    starting frame, to select what type of thing you want to do
    options are:
    """
    plot = Instance(Plot)
    color = ColorTrait("blue")
    err_alpha = Float(0.5)
    key_list = List([''])
    index_key = Int(0)
    this_key = Enum(values='key_list')
    plot_this = Bool(True)

    dim_list = List([''])
    plot_dim = Enum(values='dim_list')

    view = View(HSplit(
        Item('plot', editor=ComponentEditor(), show_label=False),
        VGroup(
            Item('this_key'),
            Item('plot_this'),
            Item('color', style="custom"),
            Item('err_alpha'),
            Item('plot_dim')
              )),
        width=500, height=500, resizable=True, title="Chaco Plot")

    def _err_alpha_changed(self):
        self.line_dict['yup_'+str(self.index_key)].alpha = self.err_alpha
        self.line_dict['ydown_'+str(self.index_key)].alpha = self.err_alpha

    def _color_changed(self):
        self.line_dict['y_'+str(self.index_key)].color = self.color
        self.line_dict['yup_'+str(self.index_key)].color = self.color
        self.line_dict['ydown_'+str(self.index_key)].color = self.color
        global data_plot
        data_plot.plot_data[self.this_key]['color'] = self.color

    def _this_key_changed(self):
        global data_plot
        self.plot_this = data_plot.plot_data[self.this_key]['plot_this']
        self.index_key = self.key_list.index(self.this_key) + 1
        self.color = hex_to_rgb(data_plot.plot_data[self.this_key]['color'])
        self.dim_list = data_plot.plot_data[self.this_key]['y_data'].index.names

    def _plot_dim_changed(self):
        global data_plot
        first_key = list(data_plot.plot_data[self.this_key]['y_data'].keys()[0])
        first_key[self.dim_list.index(self.plot_dim)] = slice(None)
        data_plot.plot_data[self.this_key]['key_select'] = first_key
        this_x,this_xerr,this_y,this_yerr,this_fit,this_key_select = data_plot.PlotElement(data_plot.plot_data[self.this_key],IgnorePlot=True)
        if this_x is None or this_x is False: return
        self.plotdata.set_data('x_'+str(self.index_key), this_x)
        self.plotdata.set_data('y_'+str(self.index_key), this_y)
        self.plotdata.set_data('yup_'+str(self.index_key), np.array(this_y) + np.array(this_yerr))
        self.plotdata.set_data('ydown_'+str(self.index_key), np.array(this_y) - np.array(this_yerr))

    def _plot_this_changed(self):
        if self.plot_this:
            self.plot.showplot(self.this_key)
            self.plot.showplot(self.this_key+'_up')
            self.plot.showplot(self.this_key+'_down')
        else:
            self.plot.hideplot(self.this_key)
            self.plot.hideplot(self.this_key+'_up')
            self.plot.hideplot(self.this_key+'_down')
        global data_plot
        data_plot.plot_data[self.this_key]['plot_this'] = self.plot_this

    def _plot_default(self):
        global data_plot
        data_dict = {}
        this_count = 1
        plot_list = []
        for icol,col_data in data_plot.plot_data.items():
            this_x,this_xerr,this_y,this_yerr,this_fit,this_key_select = data_plot.PlotElement(col_data,IgnorePlot=True)
            if this_x is None or this_x is False: continue
            plot_list.append(col_data['label'])
            data_dict['x_'+str(this_count)] = this_x
            data_dict['y_'+str(this_count)] = this_y
            data_dict['yup_'+str(this_count)] = np.array(this_y) + np.array(this_yerr)
            data_dict['ydown_'+str(this_count)] = np.array(this_y) - np.array(this_yerr)
            this_count += 1
        self.plotdata = ArrayPlotData(**data_dict)
        plot = Plot(self.plotdata)
        self.line_dict = {}
        for iplot,plot_name in enumerate(plot_list):
            this_col = simple_col[iplot%len(simple_col)]
            self.line_dict['y_'+str(iplot+1)] = plot.plot(('x_'+str(iplot+1),'y_'+str(iplot+1)),
                                                        type="line",color=this_col,name=plot_name)[0]
            self.line_dict['yup_'+str(iplot+1)] = plot.plot(('x_'+str(iplot+1),'yup_'+str(iplot+1)),
                                                       type="line",color=this_col,alpha=0.5,name=plot_name+'_up')[0]
            self.line_dict['ydown_'+str(iplot+1)] = plot.plot(('x_'+str(iplot+1),'ydown_'+str(iplot+1)),
                                                                type="line",color=this_col,alpha=0.5,
                                                                name=plot_name+'_down')[0]
            if not data_plot.plot_data[plot_name]['plot_this']:
                plot.hideplot(plot_name)
                plot.hideplot(plot_name+'_up')
                plot.hideplot(plot_name+'_down')
        plot.title = 'Test Plot'
        plot.legend.visible = True
        self.key_list = list(plot_list)
        self.index_key = 0
        self.this_key = self.key_list[0]
        self.dim_list = data_plot.plot_data[self.this_key]['y_data'].index.names
        return plot

if  __name__ == "__main__":
    Start = BeginFrame()
    Start.configure_traits()
