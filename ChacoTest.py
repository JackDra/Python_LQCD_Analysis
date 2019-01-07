#!/usr/bin/env python

# from numpy import linspace, pi, sin
# from chaco.shell import *
#
# x = linspace(-2*pi, 2*pi, 100)
# y = sin(x)
#
# plot(x, y, "r-")
# title("First plot")
# ytitle("sin(x)")
# show()

# from numpy import linspace, sin
# from traits.api import HasTraits, Instance
# from traitsui.api import View, Item
# from chaco.api import Plot, ArrayPlotData
# from enable.api import ComponentEditor
#
# class LinePlot(HasTraits):
#     plot = Instance(Plot)
#
#     traits_view = View(
#         Item('plot',editor=ComponentEditor(), show_label=False),
#         width=500, height=500, resizable=True, title="Chaco Plot")
#
#     def _plot_default(self):
#         x = linspace(-14, 14, 100)
#         y = sin(x) * x**3
#         yerr = y/(abs(x)+1)
#         yup = y + yerr
#         ydown = y - yerr
#         plotdata = ArrayPlotData(x=x, yup=yup,ydown=ydown)
#
#         plot = Plot(plotdata)
#         plot.candle_plot(("x", "ydown","yup"))
#         plot.title = "sin(x) * x^3"
#         return plot
#
# if __name__ == "__main__":
#     LinePlot().configure_traits()


# from numpy import linspace, sin
# from traits.api import HasTraits, Instance
# from traitsui.api import Item, View
# from chaco.api import ArrayPlotData, HPlotContainer, Plot
# from enable.api import ComponentEditor
#
# class ContainerExample(HasTraits):
#
#     plot = Instance(HPlotContainer)
#
#     traits_view = View(Item('plot', editor=ComponentEditor(), show_label=False),
#                        width=1000, height=600, resizable=True, title="Chaco Plot")
#
#     def _plot_default(self):
#         x = linspace(-14, 14, 100)
#         y = sin(x) * x**3
#         plotdata = ArrayPlotData(x=x, y=y)
#
#         scatter = Plot(plotdata)
#         scatter.plot(("x", "y"), type="scatter", color="blue")
#
#         line = Plot(plotdata)
#         line.plot(("x", "y"), type="line", color="blue")
#
#         container = HPlotContainer(scatter, line)
#         return container
#
# if __name__ == "__main__":
#     ContainerExample().configure_traits()


# from numpy import linspace, sin
# from traits.api import HasTraits, Instance, Int
# from traitsui.api import Item, Group, View
# from chaco.api import ArrayPlotData, marker_trait, Plot
# from enable.api import ColorTrait, ComponentEditor
#
# class ScatterPlotTraits(HasTraits):
#
#     plot = Instance(Plot)
#     color = ColorTrait("blue")
#     marker = marker_trait
#     marker_size = Int(4)
#
#     traits_view = View(
#         Group(Item('color', label="Color", style="custom"),
#               Item('marker', label="Marker"),
#               Item('marker_size', label="Size"),
#               Item('plot', editor=ComponentEditor(), show_label=False),
#                    orientation = "vertical"),
#               width=800, height=600, resizable=True, title="Chaco Plot")
#
#     def _plot_default(self):
#         x = linspace(-14, 14, 100)
#         y = sin(x) * x**3
#         plotdata = ArrayPlotData(x = x, y = y)
#
#         plot = Plot(plotdata)
#
#         self.renderer = plot.plot(("x", "y"), type="scatter", color="blue")[0]
#         return plot
#
#     def _color_changed(self):
#         self.renderer.color = self.color
#
#     def _marker_changed(self):
#         self.renderer.marker = self.marker
#
#     def _marker_size_changed(self):
#         self.renderer.marker_size = self.marker_size
#
# if __name__ == "__main__":
#     ScatterPlotTraits().configure_traits()

# from scipy.special import jn
# from numpy import linspace
# from traits.api import Enum, HasTraits, Instance
# from traitsui.api import Item, View
# from chaco.api import ArrayPlotData, Plot
# from enable.api import ComponentEditor
#
# class DataChooser(HasTraits):
#
#     plot = Instance(Plot)
#
#     data_name = Enum("jn0", "jn1", "jn2")
#
#     traits_view = View(
#         Item('data_name', label="Y data"),
#         Item('plot', editor=ComponentEditor(), show_label=False),
#         width=800, height=600, resizable=True,
#         title="Data Chooser")
#
#     def _plot_default(self):
#         x = linspace(-5, 10, 100)
#
#         # jn is the Bessel function or order n
#         self.data = {"jn0": jn(0, x),
#                      "jn1": jn(1, x),
#                      "jn2": jn(2, x)}
#
#         self.plotdata = ArrayPlotData(x = x, y = self.data["jn0"])
#
#         plot = Plot(self.plotdata)
#         plot.plot(("x", "y"), type="line", color="blue")
#         return plot
#
#     def _data_name_changed(self):
#         self.plotdata.set_data("y", self.data[self.data_name])
#
# if __name__ == "__main__":
#     DataChooser().configure_traits()


# from numpy import linspace, sin
# from traits.api import Enum, HasTraits, Instance
# from traitsui.api import Item, View
# from chaco.api import ArrayPlotData, Plot
# from chaco.tools.api import PanTool, ZoomTool
# from enable.api import ComponentEditor
#
# class PlotEditor(HasTraits):
#
#     plot = Instance(Plot)
#
#     plot_type = Enum("scatter", "line")
#
#     orientation = Enum("horizontal", "vertical")
#
#     traits_view = View(Item('orientation', label="Orientation"),
#                        Item('plot', editor=ComponentEditor(),
#                             show_label=False),
#                        width=500, height=500, resizable=True,
#                        title="Chaco Plot")
#
#     def _plot_default(self):
#         x = linspace(-14, 14, 100)
#         y = sin(x) * x**3
#         plotdata = ArrayPlotData(x = x, y = y)
#
#         plot = Plot(plotdata)
#         plot.plot(("x", "y"), type=self.plot_type, color="blue")
#
#         plot.tools.append(PanTool(plot))
#         plot.tools.append(ZoomTool(plot))
#         return plot
#
#     def _orientation_changed(self):
#         if self.orientation == "vertical":
#             self.plot.orientation = "v"
#         else:
#             self.plot.orientation = "h"
#
# if __name__ == "__main__":
#
#     # create two plots, one of type "scatter", one of type "line"
#     scatter = PlotEditor(plot_type = "scatter")
#     line = PlotEditor(plot_type = "line")
#
#     # connect the axes of the two plots
#     scatter.plot.range2d = line.plot.range2d
#
#     # open two windows
#     line.edit_traits()
#     scatter.configure_traits()

# from numpy import linspace, sin
# from traits.api import HasTraits, Instance
# from traitsui.api import Item, View
# from chaco.api import ArrayPlotData, Plot
# from chaco.tools.api import DragZoom, PanTool, ZoomTool
# from enable.api import ComponentEditor
#
# class ToolsExample(HasTraits):
#
#     plot = Instance(Plot)
#
#     traits_view = View(
#         Item('plot',editor=ComponentEditor(), show_label=False),
#         width=500, height=500,
#         resizable=True,
#         title="Chaco Plot")
#
#     def _plot_default(self):
#         x = linspace(-14, 14, 100)
#         y = sin(x) * x**3
#         plotdata = ArrayPlotData(x = x, y = y)
#         plot = Plot(plotdata)
#         plot.plot(("x", "y"), type="line", color="blue")
#
#         # append tools to pan, zoom, and drag
#         plot.tools.append(PanTool(plot))
#         plot.tools.append(ZoomTool(plot))
#         plot.tools.append(DragZoom(plot, drag_button="right"))
#         return plot
#
# if __name__ == "__main__":
#     ToolsExample().configure_traits()

# from numpy import linspace, sin
# from traits.api import HasTraits, Instance, List
# from traitsui.api import CheckListEditor, Item, View
# from chaco.api import ArrayPlotData, Plot
# from chaco.tools.api import DragZoom, PanTool, ZoomTool
# from enable.api import ComponentEditor
#
# class ToolsExample2(HasTraits):
#
#     plot = Instance(Plot)
#
#     tools = List(editor=CheckListEditor(values = ["PanTool",
#                                  "SimpleZoom", "DragZoom"]))
#
#     traits_view = View(
#         Item('plot',editor=ComponentEditor(), show_label=False),
#         Item('tools',show_label=False),
#         width=500, height=500,
#         resizable=True,
#         title="Chaco Plot")
#
#     def _plot_default(self):
#         x = linspace(-14, 14, 100)
#         y = sin(x) * x**3
#         plotdata = ArrayPlotData(x = x, y = y)
#         plot = Plot(plotdata)
#         plot.plot(("x", "y"), type="line", color="blue")
#         return plot
#
#     def _tools_changed(self):
#         classes = [eval(class_name) for class_name in self.tools]
#
#         # Remove all tools from the plot
#         plot_tools = self.plot.tools
#         for tool in plot_tools:
#             plot_tools.remove(tool)
#
#         # Create new instances for the selected tool classes
#         for cls in classes:
#             self.plot.tools.append(cls(self.plot))
#
# if __name__ == "__main__":
#     ToolsExample2().configure_traits()

# Major library imports
import numpy as np

# Enthought imports
from traits.api import (Array, Callable, Enum, Float, HasTraits, Instance, Int,
                        Trait)
from traitsui.api import Group, HGroup, Item, View, spring, Handler
from pyface.timer.api import Timer

# Chaco imports
from chaco.chaco_plot_editor import ChacoPlotItem


class Viewer(HasTraits):
    """ This class just contains the two data arrays that will be updated
    by the Controller.  The visualization/editor for this class is a
    Chaco plot.
    """
    index = Array

    data = Array

    plot_type = Enum("line", "scatter")

    view = View(ChacoPlotItem("index", "data",
                              type_trait="plot_type",
                              resizable=True,
                              x_label="Time",
                              y_label="Signal",
                              color="blue",
                              bgcolor="white",
                              border_visible=True,
                              border_width=1,
                              padding_bg_color="lightgray",
                              width=800,
                              height=380,
                              marker_size=2,
                              show_label=False),
                HGroup(spring, Item("plot_type", style='custom'), spring),
                resizable = True,
                buttons = ["OK"],
                width=800, height=500)


class Controller(HasTraits):

    # A reference to the plot viewer object
    viewer = Instance(Viewer)

    # Some parameters controller the random signal that will be generated
    distribution_type = Enum("normal", "lognormal")
    mean = Float(0.0)
    stddev = Float(1.0)

    # The max number of data points to accumulate and show in the plot
    max_num_points = Int(100)

    # The number of data points we have received; we need to keep track of
    # this in order to generate the correct x axis data series.
    num_ticks = Int(0)

    # private reference to the random number generator.  this syntax
    # just means that self._generator should be initialized to
    # random.normal, which is a random number function, and in the future
    # it can be set to any callable object.
    _generator = Trait(np.random.normal, Callable)

    view = View(Group('distribution_type',
                      'mean',
                      'stddev',
                      'max_num_points',
                      orientation="vertical"),
                      buttons=["OK", "Cancel"])

    def timer_tick(self, *args):
        """
        Callback function that should get called based on a timer tick.  This
        will generate a new random data point and set it on the `.data` array
        of our viewer object.
        """
        # Generate a new number and increment the tick count
        new_val = self._generator(self.mean, self.stddev)
        self.num_ticks += 1

        # grab the existing data, truncate it, and append the new point.
        # This isn't the most efficient thing in the world but it works.
        cur_data = self.viewer.data
        new_data = np.hstack((cur_data[-self.max_num_points+1:], [new_val]))
        new_index = np.arange(self.num_ticks - len(new_data) + 1,
                              self.num_ticks + 0.01)

        self.viewer.index = new_index
        self.viewer.data = new_data
        return

    def _distribution_type_changed(self):
        # This listens for a change in the type of distribution to use.
        if self.distribution_type == "normal":
            self._generator = np.random.normal
        else:
            self._generator = np.random.lognormal


class DemoHandler(Handler):

    def closed(self, info, is_ok):
        """ Handles a dialog-based user interface being closed by the user.
        Overridden here to stop the timer once the window is destroyed.
        """

        info.object.timer.Stop()
        return


class Demo(HasTraits):
    controller = Instance(Controller)
    viewer = Instance(Viewer, ())
    timer = Instance(Timer)
    view = View(Item('controller', style='custom', show_label=False),
                Item('viewer', style='custom', show_label=False),
                handler=DemoHandler,
                resizable=True)

    def edit_traits(self, *args, **kws):
        # Start up the timer! We should do this only when the demo actually
        # starts and not when the demo object is created.
        self.timer=Timer(100, self.controller.timer_tick)
        return super(Demo, self).edit_traits(*args, **kws)

    def configure_traits(self, *args, **kws):
        # Start up the timer! We should do this only when the demo actually
        # starts and not when the demo object is created.
        self.timer=Timer(100, self.controller.timer_tick)
        return super(Demo, self).configure_traits(*args, **kws)

    def _controller_default(self):
        return Controller(viewer=self.viewer)


# NOTE: examples/demo/demo.py looks for a 'demo' or 'popup' or 'modal popup'
# keyword when it executes this file, and displays a view for it.
popup=Demo()


if __name__ == "__main__":
    popup.configure_traits()
