#!/usr/bin/env python

import pandas as pa
import tensorflow as tf
import numpy as np
from BootStrapping import BootStrap,FlattenBootstrapDF
from Params import nboot

class Fitting_TF(object):


    def __init__(self,data=None,name='',data_bsed=False,nboot = nboot):
        self.name = name
        self.data = None
        self.data_bs = None
        self.nboot = nboot
        self.ImportData(data)
        if data_bsed:
            self.FormatBootstrap()
        ##TODO add different options for tensorflow


    def ImportData(self,data):
        if data is None:
            self.data = None
        elif isinstance(data,pa.DataFrame):
            if all([idata in list(data.columns) for idata in ['xdata','ydata']]):
                self.data = data
            else:
                print(list(data.columns))
                raise IOError('DataFrame must have a column for xdata and ydata')
        else:
            raise NotImplementedError('formatting data of type '+str(type(data))+' not implemented yet')

    ## can add functionallity to fit over configuration values opposed to bootstrapped values
    def FormatBootstrap(self):
        self.data_bs = self.data
        self.data = FlattenBootstrapDF(self.data,this_nboot=self.nboot)

    def SetupGraph(self,data=None):
        if data is not None:
            self.ImportData(data)
        xvals = self.data['xdata'].values
        yvals = self.data['ydata'].values
        if  any([isinstance(ix,BootStrap) for ix in xvals]) or \
            any([isinstance(iy,BootStrap) for iy in yvals]):
            print('Found bootstrap data in dataframe, formatting now and resetting up')
            self.FormatBootstrap()
            self.SetupGraph()

        xvals = xvals.reshape((len(xvals),1))
        yvals = yvals.reshape((len(yvals),1))

        xdata = tf.constant(xvals, dtype=tf.float32)
        ydata_true = tf.constant(yvals, dtype=tf.float32)


        linear_model = tf.layers.Dense(units=1)

        self.ydata_pred = linear_model(xdata)
        self.loss = tf.losses.mean_squared_error(labels=ydata_true, predictions=self.ydata_pred)

        self.optimizer = tf.train.GradientDescentOptimizer(0.01)
        self.train = self.optimizer.minimize(self.loss)


    def RunGraph(self,nGPU=0):
        init = tf.global_variables_initializer()
        config = tf.ConfigProto(device_count = {'GPU': 0})
        sess = tf.Session(config=config)
        # sess = tf.Session()
        sess.run(init)
        self.loss_values = []
        for i in range(100):
            _, iloss = sess.run((self.train, self.loss))
            self.loss_values.append(iloss)
        y_pred_vals = sess.run(self.ydata_pred)
        if self.data_bs is not None:
            self.data_bs['ydata_pred'] = pa.Series(y_pred_vals.flatten(),index=self.data_bs.index)
            self.data_bs['pred_err'] = (self.data_bs['ydata_pred'] - self.data_bs['ydata'].apply(lambda x : x.Avg))**2
            self.data_bs['pred_err_rel'] = self.data_bs['pred_err']/self.data_bs['ydata'].apply(lambda x : x.Avg)**2
        else:
            self.data['ydata_pred'] = pa.Series(y_pred_vals,index=self.data.index)
            self.data['pred_err'] = (self.data['ydata_pred'] - self.data['ydata'].apply(lambda x : x.Avg))**2
            self.data['pred_err_rel'] = self.data['pred_err']/self.data['ydata'].apply(lambda x : x.Avg)**2


    ## just a wrapper for setup and run
    def Fit(self,nGPU=0,data=None):
        self.SetupGraph(data=data)
        self.RunGraph(nGPU=nGPU)

def TestFit_TF(low=-5,high=5,npoints=10,nboot=200):
    values = np.random.uniform(low=low,high=high,size=(npoints,nboot))
    values = np.array([ival + icv for icv,ival in enumerate(values)])
    testdata = [BootStrap(nboot,name='test_bootstrap',cfgvals=ival,thisDelVal=False) for ival in values]

    testdf = pa.DataFrame()
    testdf.loc[:,'xdata'] = pa.Series(np.arange(npoints))
    testdf.loc[:,'ydata'] = pa.Series(testdata)
    testdf.loc[:,'ydata_Avg'] = pa.Series([idata.Avg for idata in testdata])
    testdf.loc[:,'ydata_Std'] = pa.Series([idata.Std for idata in testdata])

    # testdf_plot = pa.DataFrame()
    # testdf_plot.loc[:,'xdata'] = pa.Series(np.arange(npoints))
    # testdf_plot.loc[:,'ydata'] = pa.Series(testdata)
    # testdf_plot.loc[:,'ydata_Avg'] = pa.Series([idata.Avg for idata in testdata])
    # testdf_plot.loc[:,'ydata_Std'] = pa.Series([idata.Std for idata in testdata])

    testfit = Fitting_TF(data = testdf,name='Test',nboot=nboot)

    testfit.Fit()

    # pl.clf()
    # pl.errorbar(testdf_plot['xdata'],testdf_plot['ydata_Avg'],testdf_plot['ydata_Std'],label='data')
    # testdf_plot.plot(x='xdata',y='ydata_Avg',label='data')
    # this_data = pa.DataFrame()
    # hold_series = pa.Series()
    # hold_series['x_data'] = testdf_plot['xdata']
    # hold_series['y_data'] = testdf_plot['ydata_Avg']
    # hold_series['yerr_data'] = testdf_plot['ydata_Std']
    # hold_series['fit_class'] = None
    # hold_series['type'] = 'error_bar'
    # hold_series['label'] = 'test data'
    # hold_series['color'] = 'Blue'
    # this_data['test_data'] = hold_series
    #
    #
    # hold_series = pa.Series()
    # hold_series['type'] = 'fit'
    # hold_series['fit_class'] = testfit
    # hold_series['label'] = 'test fit'
    # hold_series['color'] = 'previous'
    # this_data['test_fit'] = hold_series
    #
    #
    # this_info = pa.Series()
    # this_info['save_file'] = './TestGraphs/TestFitPlot.pdf'
    # this_info['title'] = 'Test Fitting'
    # this_info['xlabel'] = 'Test x'
    # this_info['ylabel'] = 'Test y'
    # # this_info['xlims'] = [0,10]
    # # this_info['ylims'] = [0,15]
    #
    # import PlotData as jpl
    # data_plot = jpl.Plotting(plot_data=this_data,plot_info=this_info)
    # data_plot.LoadPickle(DefWipe=True)
    # data_plot.PrintData()
    # data_plot.PlotAll()


    # testfit.PlotFunction()
    # pl.legend()
    # pl.savefig('./TestGraphs/TestFuns.pdf')

    #
    return testfit
