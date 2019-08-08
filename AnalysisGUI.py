#!/usr/bin/env python

'''

Controls the GUI interface for running the codebase with selected parameters.

'''

JobFolder = './JobsFolder/'
JobNumber = 1

from traits.etsconfig.api import ETSConfig
ETSConfig.toolkit = 'qt'

if  __name__ == "__main__":
    from MiscFuns import mkdir_p,Remove_Empty_Str
    import sys
    mkdir_p(JobFolder)
    if len(sys.argv) > 1:
        from RunSetCommands import RunJob
        global data
        data = RunJob(sys.argv[1])
        sys.exit()
    else:
        import traits.api as ta
        import traitsui.api as tua
        from RunSetCommands import RunJob
else:
    import traits.api as ta
    import traitsui.api as tua
    from RunSetCommands import RunJob
    from MiscFuns import mkdir_p,Remove_Empty_Str
    import sys

import os

JobFolder = os.getcwd()+JobFolder[1:]

from copy import deepcopy
from XmlFormatting import tstr,untstr
from collections import OrderedDict
from GammaMatricies import GammaTOChroma,ChromaTOGamma
from RunSetCommands import FOpRun,FOpFullRun
from RunSetCommands import TPRun,ARun,AFullRun
from RunSetCommands import RRun,FFRun,FFFinish
from RunSetCommands import Out3ptRun,Out2ptRun,OutBoot2ptRun
from RunSetCommands import WriteCSH
# import SetsOfFlowOps as sfo
# import SetsOfTwoPtCorrs as so
# import SetsOfRatios as sr
# import SetsOfFFs as ffs
# import matplotlib.pyplot as pl
# import cPickle as pickle
import dill as pickle
import os
import numpy as np
from MiscFuns import GenSinkFun,CreateNameAndLL,rreload,op_str_list,DEBUGprint
import QuantityLists as ql
import Params as pa
import FileIO as IO
# from Params import *
# from QuantityLists import *
import pandas as pad

# from traits.api import *
# from traitsui.api import *

thisInfo = deepcopy(pa.defInfo)
thisInfoFlow = deepcopy(pa.defInfoFlow)
job_params = deepcopy(pa.def_job_params)
job_params_solve = deepcopy(pa.def_job_params_solve)
paramfile = pa.this_dir+'/Param_p3.py3p'
prev_file = pa.this_dir+'/PrevRun_p3.py3p'
job_paramfile = pa.this_dir+'/JobParam_p3.py3p'
job_solve_paramfile = pa.this_dir+'/JobParamSolve_p3.py3p'

thisInfoList = [thisInfo]
thisInfoFlowList = [thisInfoFlow]
thisDefWipe = False
plot_info = pad.Series()
info_index = 0


legloc = 'best'

# params = {'legend.fontsize': 7,
#           'legend.numpoints': 1,
#           'axes.labelsize' : 10,
#           'figure.autolayout': True,
#           'axes.grid': True,
#           'errorbar.capsize': 5,
#           'axes.xmargin':0.01,
#           'axes.titlesize' : 20,
#           'axes.ymargin':0.01}
# pl.style.use(pa.plot_style)
# pl.rcParams.update(params)



tsize,xsize,ysize = 20,20,20
wtitle,wxlab,wylab= False,False,False
thetitle,thexlab,theylab = 'Def','Def','Def'
def_checkcfgs = True
def_Corr_Sets = True
# Wipe_Fit = False
customleg = False
customleglist = ['']

Mcore = pa.defMcore

def SavePreviousRun():
    global thisInfoList,thisInfoFlowList
    print('Saving previous run to ',prev_file)
    with open(prev_file,'wb') as pfile:
        pickle.dump( [thisInfoList,thisInfoFlowList], pfile)


def LoadPreviousRun():
    global thisInfoList,thisInfoFlowList
    if os.path.isfile(prev_file):
        print('Reading previous run from ',prev_file)
        with open(prev_file,'rb') as pfile:
            thisInfoList,thisInfoFlowList = pickle.load(pfile)
    else:
        print('No previous pickle file',prev_file)


def PickleParams(paramdict):
    if os.path.isfile(paramfile):
        try:
            with open( paramfile, "rb" ) as pfile:
                currpdict =  pickle.load(pfile)
            currpdict.update(paramdict)
        except Exception as err:
            print('Error reading before pickling ',paramfile,' Ignoring')
            currpdict = paramdict
    else:
        currpdict = paramdict
    with open( paramfile, "wb" ) as pfile:
        pickle.dump( currpdict, pfile)

# def LoadPickle():
#     if os.path.isfile(paramfile):
#         with open( paramfile, "rb" ) as pfile:
#             return pickle.load(pfile)
#     else:
#         print 'params file ' , paramfile
#         raise EnvironmentError('params file not set up')

def JobPickleParams(paramdict):
    if os.path.isfile(job_paramfile):
        with open( job_paramfile, "rb" ) as pfile:
            currpdict =  pickle.load(pfile)
        currpdict.update(paramdict)
    else:
        currpdict = paramdict
    with open( job_paramfile, "wb" ) as pfile:
        pickle.dump( currpdict, pfile)


def JobSolvePickleParams(paramdict):
    if os.path.isfile(job_solve_paramfile):
        with open( job_solve_paramfile, "rb" ) as pfile:
            currpdict =  pickle.load(pfile)
        currpdict.update(paramdict)
    else:
        currpdict = paramdict
    with open( job_solve_paramfile, "wb" ) as pfile:
        pickle.dump( currpdict, pfile)

# def JobLoadPickle():
#     if os.path.isfile(job_paramfile):
#         with open( job_paramfile, "rb" ) as pfile:
#             return pickle.load(pfile)
#     else:
#         print 'job_params file ' , job_paramfile
#         raise EnvironmentError('job_params file not set up')

class CalculationParams( ta.HasTraits ):
    """

    Contains general parameters for analysis

    """

    global thisDefWipe,def_def_checkcfgs,def_Corr_Sets
    scratch_hotfix = ta.Bool(pa.scratch_hotfix)
    Debug_Mode = ta.Bool(pa.Debug_Mode)
    force_wipe = ta.Bool(thisDefWipe)
    Check_Configurations = ta.Bool(def_checkcfgs)
    Correlate_Sets = ta.Bool(def_Corr_Sets)
    Wipe_All_Fits = ta.Bool(pa.Wipe_All_Fits)
    Multicore = ta.Bool(pa.defMcore)
    num_of_procs = ta.Int(pa.defnProc)
    max_job_per_thread = ta.Int(pa.defMaxTasks)
    write_excel = ta.Bool(pa.write_excel)
    this_list = ql.file_formats
    this_list.insert(0, this_list.pop(this_list.index(pa.cfg_file_type)))
    cfg_file_type = ta.Enum(this_list)
    Show_Write = ta.Bool(pa.ShowXmlWrite)
    Show_Read = ta.Bool(pa.ShowRead)
    Max_Iters = ta.Int(pa.defMaxIter)
    Solver_Prec = ta.Float(pa.defPrec)
    data_directory = ta.Str(pa.datadir)
    output_directory = ta.Str(pa.outputdir)
    scratch_directory = ta.Str(pa.scratchdir)
    formatted_cfg_dir = ta.Str(pa.cfgfmtdir)
    cfun_prefix = ta.Str(pa.cfundir.replace(pa.datadir,''))
    flow_new_format = ta.Bool(pa.flow_new_format)
    n_boot = ta.Int(pa.nboot)
    nchi_threshold = ta.Int(pa.nchi_threshold)
    nchi_fit_threshold = ta.Int(pa.nchi_fit_threshold)
    plot_the_ff = ta.Bool(pa.plot_the_ff)
    normed_histograms = ta.Bool(pa.normed_histograms)
    def_min_fitr = ta.Int(pa.def_min_fitr)

    Apply = ta.Button()


    view = tua.View(
        '_',
        tua.Item('scratch_hotfix'),
        tua.Item('Debug_Mode'),
        tua.Item('force_wipe'),
        tua.Item('Check_Configurations'),
        tua.Item('Correlate_Sets'),
        tua.Item('Wipe_All_Fits'),
        tua.Item('write_excel'),
        tua.Item('cfg_file_type'),
        tua.Item('Show_Write'),
        tua.Item('Show_Read'),
        tua.Item('Multicore'),
        tua.Item('num_of_procs',enabled_when='Multicore'),
        tua.Item('max_job_per_thread',enabled_when='Multicore'),
        tua.Item('Max_Iters'),
        tua.Item('Solver_Prec'),
        tua.Item('data_directory'),
        tua.Item('output_directory'),
        tua.Item('scratch_directory'),
        tua.Item('formatted_cfg_dir'),
        tua.Item('flow_new_format'),
        tua.Item('cfun_prefix'),
        tua.Item('n_boot'),
        tua.Item('nchi_threshold'),
        tua.Item('nchi_fit_threshold'),
        tua.Item('plot_the_ff'),
        tua.Item('normed_histograms'),
        tua.Item('def_min_fitr'),
        tua.Item('Apply', show_label=False),
        buttons=['OK'],
        width=0.3,
        resizable=True

    )

    def _Apply_fired(self):
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,Mcore
        thisDefWipe = self.force_wipe
        def_checkcfgs = self.Check_Configurations
        def_Corr_Sets = self.Correlate_Sets
        Mcore = self.Multicore
        PickleParams(self.__dict__)
        rreload(pa)
        for iindex in range(len(thisInfoList)):
            thisInfoList[iindex]['outdir'] = self.output_directory
            thisInfoFlowList[iindex]['outdir'] = self.output_directory
        print('Please restart program for changes to take effect!')

class RunParams( ta.HasTraits ):
    """

    Contains parameters for the particular ensemble to be analysed

    """

    n_space = ta.Int(pa.defInfo['nxyzt'][0]) ## space dimensions are all the same
    n_time = ta.Int(pa.defInfo['nxyzt'][-1]) ## time is last one
    lattice_spacing = ta.Float(pa.defInfo['a']) ## in fermi
    kappa_ud = ta.Float(float('0.'+str(pa.defInfo['kud'])))
    kappa_s = ta.Float(float('0.'+str(pa.defInfo['ks'])))
    t_sink = ta.Int(pa.defInfo['tsink'])
    Csw = ta.Int(pa.defInfo['Csw'])
    Beta = ta.Int(pa.defInfo['Beta'])
    Stats_from_fo = ta.Enum(ql.FOcfglist)
    Pick_streams = ta.Bool(False)
    Stream_list = ta.List(ta.Str,['-a-','-b-'])
    Ensemble = ta.Enum(list(ql.ens_dict.keys()))
    Load = ta.Button()
    Apply = ta.Button()
    Apply_All = ta.Button()
    Setup_Type = ta.Enum(list(ql.setup_dict.keys()))
    Apply_Setup = ta.Button()


    view = tua.View(
        '_',
        tua.Item('n_space'),
        tua.Item('n_time'),
        tua.Item('lattice_spacing'),
        tua.Item('kappa_ud'),
        tua.Item('kappa_s'),
        tua.Item('Csw'),
        tua.Item('Beta'),
        tua.Item('t_sink'),
        tua.Item('Stats_from_fo'),
        tua.Item('Pick_streams'),
        tua.Item('Stream_list',enabled_when='Pick_streams'),
        tua.Item('Ensemble'),
        tua.Item('Load', show_label=False),
        tua.Item('Apply', show_label=False),
        tua.Item('Apply_All', show_label=False),
        tua.Item('Setup_Type', show_label=False),
        tua.Item('Apply_Setup', show_label=False),
        resizable=True,
        buttons=['OK']
    )


    def _Load_fired(self):
        this_dict = ql.ens_dict[self.Ensemble]
        self.n_space,self.n_time = this_dict['nxyzt']
        self.lattice_spacing = this_dict['a']
        self.kappa_ud = this_dict['kud']
        self.kappa_s = this_dict['ks']
        self.Csw = this_dict['Csw']
        self.Beta = this_dict['Beta']
        self.t_sink = this_dict['tsink']
        if 'stream_list' in list(this_dict.keys()):
            self.Pick_streams = True
            self.Stream_list = this_dict['stream_list']
        else:
            self.Pick_streams = False

    def _Apply_fired(self):
        thisdict = self.__dict__
        if '__traits_listener__' in list(thisdict.keys()): del thisdict['__traits_listener__']
        PickleParams(thisdict)
        thisInfoList[info_index]['nxyzt'] = [self.n_space,self.n_space,self.n_space,self.n_time]
        thisInfoList[info_index]['a'] = self.lattice_spacing
        ## this multipilcation is a hack, probably need better way of dealing with this
        thisInfoList[info_index]['kud'] = int(str(self.kappa_ud).replace('0.',''))
        thisInfoList[info_index]['ks'] = int(str(self.kappa_s).replace('0.',''))
        multfact_ud = 7-len(str(thisInfoList[info_index]['kud']))
        multfact_s = 7-len(str(thisInfoList[info_index]['ks']))
        thisInfoList[info_index]['kud'] = thisInfoList[info_index]['kud']*10**multfact_ud
        thisInfoList[info_index]['ks'] = thisInfoList[info_index]['ks']*10**multfact_s
        thisInfoList[info_index]['Csw'] = self.Csw
        thisInfoList[info_index]['Beta'] = self.Beta
        thisInfoList[info_index]['tsink'] = self.t_sink
        thisInfoFlowList[info_index]['Csw'] = self.Csw
        thisInfoFlowList[info_index]['Beta'] = self.Beta
        thisInfoFlowList[info_index]['tsink'] = self.t_sink
        if self.Stats_from_fo == 'None':
            thisInfoList[info_index]['fo_for_cfgs'] = False
            thisInfoFlowList[info_index]['fo_for_cfgs'] = False
        else:
            thisInfoList[info_index]['fo_for_cfgs'] = self.Stats_from_fo
            thisInfoFlowList[info_index]['fo_for_cfgs'] = self.Stats_from_fo
        thisInfoFlowList[info_index]['nxyzt'] = [self.n_space,self.n_space,self.n_space,self.n_time]
        thisInfoFlowList[info_index]['a'] = thisInfoList[info_index]['a']
        thisInfoFlowList[info_index]['kud'] = thisInfoList[info_index]['kud']
        thisInfoFlowList[info_index]['ks'] = thisInfoList[info_index]['ks']
        if bool(self.Pick_streams):
            thisInfoFlowList[info_index]['stream_list'] = Remove_Empty_Str(list(self.Stream_list))
            thisInfoList[info_index]['stream_list'] = Remove_Empty_Str(list(self.Stream_list))
        else:
            thisInfoFlowList[info_index]['stream_list'] = 'All'
            thisInfoList[info_index]['stream_list'] = 'All'
        this_dict = ql.ens_dict[self.Ensemble]
        thisInfoFlowList[info_index]['jsmCoeffs'] = this_dict['jsmCoeffs']
        thisInfoFlowList[info_index]['Fits'] = this_dict['Fits']
        thisInfoFlowList[info_index]['FitsAlpha'] = this_dict['FitsAlpha']
        thisInfoFlowList[info_index]['FitsAlphaFull'] = this_dict['FitsAlphaFull']
        thisInfoFlowList[info_index]['Rat_tsum_range'] = this_dict['Rat_tsum_range']


    def _Apply_All_fired(self):
        thisdict = self.__dict__
        if '__traits_listener__' in list(thisdict.keys()): del thisdict['__traits_listener__']
        PickleParams(thisdict)
        for iindex in range(len(thisInfoList)):
            thisInfoList[iindex]['nxyzt'] = [self.n_space,self.n_space,self.n_space,self.n_time]
            thisInfoList[iindex]['a'] = self.lattice_spacing
            ## this multipilcation is a hack, probably need better way of dealing with this
            thisInfoList[iindex]['kud'] = int(str(self.kappa_ud).replace('0.',''))
            thisInfoList[iindex]['ks'] = int(str(self.kappa_s).replace('0.',''))
            multfact_ud = 7-len(str(thisInfoList[iindex]['kud']))
            multfact_s = 7-len(str(thisInfoList[iindex]['ks']))
            thisInfoList[iindex]['kud'] = thisInfoList[iindex]['kud']*10**multfact_ud
            thisInfoList[iindex]['ks'] = thisInfoList[iindex]['ks']*10**multfact_s
            thisInfoList[iindex]['Csw'] = self.Csw
            thisInfoList[iindex]['Beta'] = self.Beta
            thisInfoList[iindex]['tsink'] = self.t_sink
            thisInfoFlowList[iindex]['Csw'] = self.Csw
            thisInfoFlowList[iindex]['Beta'] = self.Beta
            thisInfoFlowList[iindex]['tsink'] = self.t_sink
            if self.Stats_from_fo == 'None':
                thisInfoList[iindex]['fo_for_cfgs'] = False
                thisInfoFlowList[iindex]['fo_for_cfgs'] = False
            else:
                thisInfoList[iindex]['fo_for_cfgs'] = self.Stats_from_fo
                thisInfoFlowList[iindex]['fo_for_cfgs'] = self.Stats_from_fo
            thisInfoFlowList[iindex]['nxyzt'] = [self.n_space,self.n_space,self.n_space,self.n_time]
            thisInfoFlowList[iindex]['a'] = thisInfoList[iindex]['a']
            thisInfoFlowList[iindex]['kud'] = thisInfoList[iindex]['kud']
            thisInfoFlowList[iindex]['ks'] = thisInfoList[iindex]['ks']
            if bool(self.Pick_streams):
                thisInfoFlowList[iindex]['stream_list'] = Remove_Empty_Str(list(self.Stream_list))
                thisInfoList[iindex]['stream_list'] = Remove_Empty_Str(list(self.Stream_list))
            else:
                thisInfoFlowList[iindex]['stream_list'] = 'All'
                thisInfoList[iindex]['stream_list'] = 'All'
            this_dict = ql.ens_dict[self.Ensemble]
            thisInfoFlowList[iindex]['Fits'] = this_dict['Fits']
            thisInfoFlowList[iindex]['jsmCoeffs'] = this_dict['jsmCoeffs']
            thisInfoFlowList[iindex]['FitsAlpha'] = this_dict['FitsAlpha']
            thisInfoFlowList[iindex]['FitsAlphaFull'] = this_dict['FitsAlphaFull']
            thisInfoFlowList[iindex]['Rat_tsum_range'] = this_dict['Rat_tsum_range']


    def _Apply_Setup_fired(self):
        thisdict = self.__dict__
        if '__traits_listener__' in list(thisdict.keys()): del thisdict['__traits_listener__']
        PickleParams(thisdict)
        # if len(thisInfoList) > 1:
        #     print 'please only have 1 dictionary set up!'
        # else:
        for idict,this_dict in enumerate(ql.setup_dict[str(self.Setup_Type)]):
            self.n_space,self.n_time = this_dict['nxyzt']
            self.lattice_spacing = this_dict['a']
            self.kappa_ud = this_dict['kud']
            self.kappa_s = this_dict['ks']
            self.Csw = this_dict['Csw']
            self.Beta = this_dict['Beta']
            self.t_sink = this_dict['tsink']
            if 'stream_list' in list(this_dict.keys()):
                self.Pick_streams = True
                self.Stream_list = this_dict['stream_list']
            else:
                self.Pick_streams = False
            thisInfoList[idict]['nxyzt'] = [self.n_space,self.n_space,self.n_space,self.n_time]
            thisInfoList[idict]['a'] = self.lattice_spacing
            ## this multipilcation is a hack, probably need better way of dealing with this
            thisInfoList[idict]['kud'] = int(str(self.kappa_ud).replace('0.',''))
            thisInfoList[idict]['ks'] = int(str(self.kappa_s).replace('0.',''))
            multfact_ud = 7-len(str(thisInfoList[idict]['kud']))
            multfact_s = 7-len(str(thisInfoList[idict]['ks']))
            thisInfoList[idict]['kud'] = thisInfoList[idict]['kud']*10**multfact_ud
            thisInfoList[idict]['ks'] = thisInfoList[idict]['ks']*10**multfact_s
            thisInfoList[idict]['Csw'] = self.Csw
            thisInfoList[idict]['Beta'] = self.Beta
            thisInfoList[idict]['tsink'] = self.t_sink
            thisInfoFlowList[idict]['Csw'] = self.Csw
            thisInfoFlowList[idict]['Beta'] = self.Beta
            thisInfoFlowList[idict]['tsink'] = self.t_sink
            if self.Stats_from_fo == 'None':
                thisInfoList[idict]['fo_for_cfgs'] = False
                thisInfoFlowList[idict]['fo_for_cfgs'] = False
            else:
                thisInfoList[idict]['fo_for_cfgs'] = self.Stats_from_fo
                thisInfoFlowList[idict]['fo_for_cfgs'] = self.Stats_from_fo
            thisInfoFlowList[idict]['nxyzt'] = [self.n_space,self.n_space,self.n_space,self.n_time]
            thisInfoFlowList[idict]['a'] = thisInfoList[idict]['a']
            thisInfoFlowList[idict]['kud'] = thisInfoList[idict]['kud']
            thisInfoFlowList[idict]['ks'] = thisInfoList[idict]['ks']
            if bool(self.Pick_streams):
                thisInfoFlowList[idict]['stream_list'] = Remove_Empty_Str(list(self.Stream_list))
                thisInfoList[idict]['stream_list'] = Remove_Empty_Str(list(self.Stream_list))
            else:
                thisInfoFlowList[idict]['stream_list'] = 'All'
                thisInfoList[idict]['stream_list'] = 'All'
            thisInfoFlowList[idict]['jsmCoeffs'] = this_dict['jsmCoeffs']
            thisInfoFlowList[idict]['Fits'] = this_dict['Fits']
            thisInfoFlowList[idict]['FitsAlpha'] = this_dict['FitsAlpha']
            thisInfoFlowList[idict]['FitsAlphaFull'] = this_dict['FitsAlphaFull']
            thisInfoFlowList[idict]['Rat_tsum_range'] = this_dict['Rat_tsum_range']
            if idict < len(thisInfoList) and len(thisInfoList) < len(ql.setup_dict[str(self.Setup_Type)]):
                thisInfoList.append(deepcopy(thisInfoList[idict]))
                thisInfoFlowList.append(deepcopy(thisInfoFlowList[idict]))

class PlotParams( ta.HasTraits ):
    """

    DEPRECIATED

    """
    legend_location = ta.Enum(ql.LegLocList)
    force_title = ta.Bool(False)
    force_x_label = ta.Bool(False)
    force_y_label = ta.Bool(False)
    title = ta.Str()
    x_label = ta.Str()
    y_label = ta.Str()
    title_size = ta.Int(tsize)
    x_label_size = ta.Int(xsize)
    y_label_size = ta.Int(ysize)
    x_axis_min = ta.Str('None')
    x_axis_max = ta.Str('None')
    y_axis_min = ta.Str('None')
    y_axis_max = ta.Str('None')
    # Custom_legend = ta.Bool(customleg)
    # Custom_legend_list = ta.List(ta.Str,customleglist)
    Apply = ta.Button()


    view = tua.View(
        '_',
        tua.Item('legend_location'),
        tua.Item('force_title'),
        tua.Item('title',enabled_when='force_title'),
        tua.Item('force_x_label'),
        tua.Item('x_label',enabled_when='force_x_label'),
        tua.Item('force_y_label'),
        tua.Item('y_label',enabled_when='force_y_label'),
        tua.Item('title_size'),
        tua.Item('x_label_size'),
        tua.Item('y_label_size'),
        tua.Item('x_axis_min'),
        tua.Item('x_axis_max'),
        tua.Item('y_axis_min'),
        tua.Item('y_axis_max'),
        # tua.Item('Custom_legend'),
        # tua.Item('Custom_legend_list',enabled_when='Custom_legend'),
        tua.Item('Apply', show_label=False),
        resizable=True,
        buttons=['OK']
    )

    def _Apply_fired(self):
        global plot_info
        # global customleg,customleglist
        # global wtitle,wxlab,wrlab
        # global legloc
        plot_info['leg_loc'] = self.legend_location
        if self.force_title:
            plot_info['title'] = self.title
        plot_info['title_dict'] = {'fontsize':self.title_size}
        if self.force_y_label:
            plot_info['ylabel'] = self.y_label
        plot_info['ylabel_dict'] = {'fontsize':self.y_label_size}
        if self.force_x_label:
            plot_info['xlabel'] = self.x_label
        plot_info['xlabel_dict'] = {'fontsize':self.x_label_size}
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


class JobParams( ta.HasTraits ):
    """

    When calling 'Create job' or 'Create separate job', the job will have these parameters

    """
    global JobNumber
    def_hour,def_min,def_seconds = list(map(int,job_params['time'].split(':')))
    machine_name = ta.Str(job_params['machine'])
    python_executable = ta.Str(job_params['python_exe'])
    submit_command = ta.Str(job_params['Scom'])
    run_time_note = ta.Str('Run_Time')
    hours = ta.Int(def_hour)
    minutes = ta.Int(def_min)
    seconds = ta.Int(def_seconds)
    memory = ta.Str(job_params['memory'])
    if 'email' in job_params:
        email = ta.Str(job_params['email'])
    else:
        email = ta.Str('')
    nproc = ta.Int(job_params['nproc'])
    run_ptg = ta.Bool(job_params['RunPTG'])
    queue_type = ta.Str(job_params['quetype'])
    Job_Number = ta.Int(JobNumber)
    Apply = ta.Button()


    view = tua.View(
        '_',
        tua.Item('machine_name'),
        tua.Item('submit_command'),
        tua.Item('python_executable'),
        tua.Item('run_time_note',show_label=False,style='readonly'),
        tua.Item('hours'),
        tua.Item('minutes'),
        tua.Item('seconds'),
        tua.Item('memory'),
        tua.Item('email'),
        tua.Item('nproc'),
        tua.Item('run_ptg'),
        tua.Item('queue_type'),
        tua.Item('Job_Number'),
        tua.Item('Apply', show_label=False),
        resizable=True,
        buttons=['OK']
    )

    def _Apply_fired(self):
        this_job_params = OrderedDict()
        this_job_params['machine'] = str(self.machine_name)
        this_job_params['Scom'] = str(self.submit_command)
        this_job_params['python_exe'] = str(self.python_executable)
        str_hours = str(self.hours).zfill(2)
        str_minutes = str(self.minutes).zfill(2)
        str_seconds = str(self.seconds).zfill(2)
        this_job_params['time'] = ':'.join([str_hours,str_minutes,str_seconds])
        if len(str(self.email)) > 0:
            this_job_params['email'] = str(self.email)
        this_job_params['nproc'] = int(self.nproc)
        this_job_params['RunPTG'] = bool(self.run_ptg)
        this_job_params['quetype'] = str(self.queue_type)
        this_job_params['memory'] = str(self.memory)
        JobPickleParams(this_job_params)
        rreload(pa)
        global job_params,JobNumber
        JobNumber = self.Job_Number
        job_params = deepcopy(this_job_params)


class JobParamsSolve( ta.HasTraits ):
    """

    When calling 'Create job' or 'Create separate job', the job will have these parameters
    when running the 'solve' part of the system of equations.


    """
    global JobNumber
    def_hour,def_min,def_seconds = list(map(int,job_params_solve['time'].split(':')))
    machine_name = ta.Str(job_params_solve['machine'])
    python_executable = ta.Str(job_params_solve['python_exe'])
    submit_command = ta.Str(job_params_solve['Scom'])
    run_time_note = ta.Str('Run_Time')
    hours = ta.Int(def_hour)
    minutes = ta.Int(def_min)
    seconds = ta.Int(def_seconds)
    memory = ta.Str(job_params_solve['memory'])
    if 'email' in job_params_solve:
        email = ta.Str(job_params_solve['email'])
    else:
        email = ta.Str('')
    nproc = ta.Int(job_params_solve['nproc'])
    run_ptg = ta.Bool(job_params_solve['RunPTG'])
    queue_type = ta.Str(job_params_solve['quetype'])
    Job_Number = ta.Int(JobNumber)
    Apply = ta.Button()


    view = tua.View(
        '_',
        tua.Item('machine_name'),
        tua.Item('submit_command'),
        tua.Item('python_executable'),
        tua.Item('run_time_note',show_label=False,style='readonly'),
        tua.Item('hours'),
        tua.Item('minutes'),
        tua.Item('seconds'),
        tua.Item('memory'),
        tua.Item('email'),
        tua.Item('nproc'),
        tua.Item('run_ptg'),
        tua.Item('queue_type'),
        tua.Item('Job_Number'),
        tua.Item('Apply', show_label=False),
        resizable=True,
        buttons=['OK']
    )

    def _Apply_fired(self):
        this_job_params_solve = OrderedDict()
        this_job_params_solve['machine'] = str(self.machine_name)
        this_job_params_solve['Scom'] = str(self.submit_command)
        this_job_params_solve['python_exe'] = str(self.python_executable)
        str_hours = str(self.hours).zfill(2)
        str_minutes = str(self.minutes).zfill(2)
        str_seconds = str(self.seconds).zfill(2)
        this_job_params_solve['time'] = ':'.join([str_hours,str_minutes,str_seconds])
        if len(str(self.email)) > 0:
            this_job_params_solve['email'] = str(self.email)
        this_job_params_solve['nproc'] = int(self.nproc)
        this_job_params_solve['RunPTG'] = bool(self.run_ptg)
        this_job_params_solve['quetype'] = str(self.queue_type)
        this_job_params_solve['memory'] = str(self.memory)
        JobSolvePickleParams(this_job_params_solve)
        rreload(pa)
        global job_params_solve,JobNumber
        JobNumber = self.Job_Number
        job_params_solve = deepcopy(this_job_params_solve)


class TPCorr( ta.HasTraits ):
    """

    Controler for running the two-point correlation function analysis

    """
    reformat_directory = ta.Str(pa.outputdir+'/Reformat2pt/')
    boot_directory = ta.Str(pa.outputdir+'/Boot2pt/')
    var_extra_mpi = []
    test_subtraction = ta.Bool(False)
    mass_ratio_guess = ta.Float(3.1)
    sub_fit_state = ta.Str('state2')
    sub_fit_range = ta.Str('fitr3-20')
    sub_fit_guess = ta.List(ta.Float,[1/77000.,50,0.5,0.03])
    norm_plot = ta.Bool(True)
    Only_Var_plot = ta.Bool(False)
    Physical_units = ta.Bool(True)
    Mes_or_Bar = ta.Enum(ql.BMList)
    Bar_Mes_List = ta.List(ql.BarNumbList)
    Bar_Mes_Number = ta.Enum(values = 'Bar_Mes_List')
    combine_all_corrs = ta.Enum(ql.CombineList)
    Power_of_corr = ta.Str('None')
    Coeff_of_corr = ta.Float(1.0)
    source_smearing = ta.Int(pa.defInfo['ism'])
    t_source = ta.Int(pa.defInfo['tsrc'])
    combine_sinks = ta.Bool('sm2ptRat' in list(pa.defInfo.keys()))
    not_combine_sinks = ta.Bool(not ('sm2ptRat' in list(pa.defInfo.keys())))
    sink_smearing = ta.Int(pa.defInfo['jsm'])
    smlist = list(map(str,pa.defInfo['sm2ptRat']))
    smear_list = ta.Str(', '.join(['sm'+ival for ival in smlist]))
    smear_comment = ta.Str(' Set coeff to 0 to ignore smearing ')
    sink_coeffs = ta.Array(np.float64,(len(smlist),),pa.defInfo['jsmCoeffs'])
    sink_label = ta.Str(pa.defInfo['sinklab'])
    EffM_state_fit = ta.Str(pa.defInfo['Fits'][0][0])
    EffM_fit_range = ta.Str(pa.defInfo['Fits'][0][1])
    momentum_list = ta.List(ta.Str,pa.defpmomlist)
    momentum_note = ta.Str(' Only the first momentum will be plotted, order correctly ')
    do_var_method = ta.Bool(True)
    prony_ism = ta.Int(0)
    var_symetrize = ta.Bool(True)
    var_t0 = ta.Int(2)
    var_dt = ta.Int(2)
    # pformlist = ['p'+imom.replace(' ','') for imom in momentum_list]



    title_size = ta.Int(20)
    ylab_size = ta.Int(20)
    xlab_size = ta.Int(20)
    RUN = ta.Button()
    CreateJob = ta.Button()
    CreateSeparateJob = ta.Button()
    Reformat = ta.Button()
    WriteBoot = ta.Button()
    Load = ta.Button()
    Apply = ta.Button()


    view = tua.View(
        tua.Group(
        tua.Item('reformat_directory'),
        tua.Item('boot_directory'),
        tua.Item('norm_plot'),
        tua.Item('Only_Var_plot'),
        tua.Item('Physical_units'),
        tua.Item('Mes_or_Bar'),
        tua.Item('Bar_Mes_Number'),
        tua.Item('source_smearing'),
        tua.Item('t_source'),
        tua.Item('combine_sinks'),
        tua.Item('sink_smearing',enabled_when='not_combine_sinks'),
        tua.Item('smear_list',enabled_when='combine_sinks',style='readonly'),
        tua.Item('smear_comment',show_label=False,enabled_when='combine_sinks',style='readonly'),
        tua.Item('sink_coeffs',enabled_when='combine_sinks',style='custom'),
        tua.Item('sink_label',enabled_when='combine_sinks'),
        tua.Item('EffM_state_fit'),
        tua.Item('EffM_fit_range'),
        tua.Item('momentum_note',show_label=False,style='readonly'),
        tua.Item('momentum_list'),
        # tua.Item('Load', show_label=False),
        tua.Item('Apply', show_label=False),
        tua.Item('CreateJob', show_label=False ),
        tua.Item('CreateSeparateJob', show_label=False ),
        tua.Item('Reformat',show_label=False),
        tua.Item('WriteBoot',show_label=False),
        tua.Item('RUN', show_label=False ),
        label='Page 1'
        ),tua.Group(
        tua.Item('combine_all_corrs'),
        tua.Item('Power_of_corr'),
        tua.Item('Coeff_of_corr'),
        tua.Item('test_subtraction'),
        tua.Item('mass_ratio_guess',enabled_when='test_subtraction'),
        tua.Item('sub_fit_state',enabled_when='test_subtraction'),
        tua.Item('sub_fit_range',enabled_when='test_subtraction'),
        tua.Item('sub_fit_guess',enabled_when='test_subtraction'),
        # tua.Item('Load', show_label=False),
        tua.Item('Apply', show_label=False),
        tua.Item('CreateJob', show_label=False ),
        tua.Item('CreateSeparateJob', show_label=False ),
        tua.Item('Reformat',show_label=False),
        tua.Item('WriteBoot',show_label=False),
        tua.Item('RUN', show_label=False ),
        label='Page 2'
        ),tua.Group(
        tua.Item('do_var_method'),
        tua.Item('prony_ism',enabled_when='do_var_method'),
        tua.Item('var_symetrize',enabled_when='do_var_method'),
        tua.Item('var_t0',enabled_when='do_var_method'),
        tua.Item('var_dt',enabled_when='do_var_method'),
        # tua.Item('Load', show_label=False),
        tua.Item('Apply', show_label=False),
        tua.Item('CreateJob', show_label=False ),
        tua.Item('CreateSeparateJob', show_label=False ),
        tua.Item('Reformat',show_label=False),
        tua.Item('WriteBoot',show_label=False),
        tua.Item('RUN', show_label=False ),
        label='Variational Method'
        ),
        resizable=True
    )

    def _Mes_or_Bar_changed(self):
        if self.Mes_or_Bar == 'Baryon':
            self.Bar_Mes_List = ql.BarNumbList
        elif self.Mes_or_Bar == 'Meson':
            self.Bar_Mes_List = ql.MesNumbList


    def _combine_sinks_changed(self):
        self.not_combine_sinks = not self.combine_sinks

    # def __init__(self):
    #     self._Load_fired()

    def _Load_fired(self):
        if 'MpiPlotParams' in list(thisInfoList[info_index].keys()):
            self.EffM_state_fit = thisInfoList[info_index]['MpiPlotParams'][0]
            self.EffM_fit_range = thisInfoList[info_index]['MpiPlotParams'][-1]
        if 'MesOrBar' in list(thisInfoList[info_index].keys()):
            self.Mes_or_Bar = thisInfoList[info_index]['MesOrBar']
            if self.Mes_or_Bar == 'Meson':
                self.Bar_Mes_Number = ChromaTOGamma(thisInfoList[info_index]['Interp'])
            elif self.Mes_or_Bar == 'Baryon':
                self.Bar_Mes_Number = thisInfoList[info_index]['Interp']


        if 'Popp' in list(thisInfoList[info_index].keys()):
            self.Power_of_corr      = thisInfoList[info_index]['Popp']
        if 'coeffopp' in list(thisInfoList[info_index].keys()):
            self.Coeff_of_corr      = thisInfoList[info_index]['coeffopp']
        if 'ism' in list(thisInfoList[info_index].keys()):
            self.source_smearing    = thisInfoList[info_index]['ism']
        if 'tsrc' in list(thisInfoList[info_index].keys()):
            self.t_source           = thisInfoList[info_index]['tsrc']
        if 'do_var_method' in list(thisInfoList[info_index].keys()):
            self.do_var_method           = thisInfoList[info_index]['do_var_method']
        if 'prony_ism' in list(thisInfoList[info_index].keys()):
            self.prony_ism           = thisInfoList[info_index]['prony_ism']
        if 'var_symetrize' in list(thisInfoList[info_index].keys()):
            self.var_symetrize           = thisInfoList[info_index]['var_symetrize']
        if 'var_t0' in list(thisInfoList[info_index].keys()):
            self.var_t0           = thisInfoList[info_index]['var_t0']
        if 'var_dt' in list(thisInfoList[info_index].keys()):
            self.var_dt           = thisInfoList[info_index]['var_dt']
        if 'jsm' in list(thisInfoList[info_index].keys()):
            self.sink_smearing      = thisInfoList[info_index]['jsm']
        if 'combine_sinks' in list(thisInfoList[info_index].keys()):
            self.combine_sinks =      thisInfoList[info_index]['combine_sinks']
        if 'sinklab' in list(thisInfoList[info_index].keys()):
            self.sink_label    =      thisInfoList[info_index]['sinklab']
        if 'pmom' in list(thisInfoList[info_index].keys()):
            if thisInfoList[info_index]['pmom'] == 'All':
                self.momentum_list = ['0 0 0','1 0 0','1 1 0','1 1 1','2 0 0']
            else:
                self.momentum_list = thisInfoList[info_index]['pmom']
        self.pformlist = [str('p'+imom.replace(' ','')) for imom in thisInfoList[info_index]['pmom']]
        if 'ComCoffs' in list(thisInfoList[info_index].keys()):
            self.Coeff_of_corr = thisInfoList[info_index]['CombCoffs']
        if 'CombOpp'  in list(thisInfoList[info_index].keys()):
            self.combine_all_corrs = thisInfoList[info_index]['CombOpp']


    def _Apply_fired(self):
        thisInfoList[info_index]['MpiPlotParams'] = (str(self.EffM_state_fit),str(self.momentum_list[0]),str(self.EffM_fit_range))
        thisInfoList[info_index]['MesOrBar'] = self.Mes_or_Bar
        thisInfoList[info_index]['Popp'] = self.Power_of_corr
        thisInfoList[info_index]['coeffopp'] = self.Coeff_of_corr
        if self.Mes_or_Bar == 'Meson':
            thisInfoList[info_index]['Interp'] = GammaTOChroma(self.Bar_Mes_Number)
        else:
            thisInfoList[info_index]['Interp'] = self.Bar_Mes_Number

        thisInfoList[info_index]['ism'] = self.source_smearing

        thisInfoList[info_index]['tsrc'] = self.t_source
        thisInfoList[info_index]['jsm'] = self.sink_smearing
        thisInfoList[info_index]['do_var_method'] = self.do_var_method
        if self.do_var_method:
            thisInfoList[info_index]['prony_ism'] = self.prony_ism
            thisInfoList[info_index]['var_symetrize'] = self.var_symetrize
            thisInfoList[info_index]['var_t0'] = self.var_t0
            thisInfoList[info_index]['var_dt'] = self.var_dt

        thisInfoList[info_index]['Fits'] = deepcopy([(str(self.EffM_state_fit),str(self.EffM_fit_range))])
        thisInfoList[info_index]['pmom'] = Remove_Empty_Str(list(map(str,self.momentum_list)))
        self.pformlist = [str('p'+imom.replace(' ','')) for imom in thisInfoList[info_index]['pmom']]
        thisInfoList[info_index]['combine_sinks'] = self.combine_sinks
        thisInfoList[info_index]['sinklab'] = self.sink_label
        if self.combine_sinks:
            thisInfoList[info_index]['sinklen'] = len(self.smlist)
            smint,scoeffs = list(map(int,self.smlist)),list(map(np.float64,self.sink_coeffs))

            # if not all([ilist['jsm'] == jsm for ilist,jsm in zip(thisInfoList[info_index:info_index+3],smint)]):
            for icount,(jsm,icoeff) in enumerate(zip(smint,scoeffs)):
                if len(thisInfoList) < info_index+icount+1:
                    thisInfoList.append(deepcopy(thisInfoList[info_index]))
                    thisInfoFlowList.append(deepcopy(thisInfoFlowList[info_index]))
                thisInfoList[info_index+icount]['combine_sinks'] = True
                thisInfoList[info_index+icount]['jsm'] = jsm
                thisInfoList[info_index+icount]['jsmCoeffs'] = [icoeff]
                self.var_extra_mpi.append((icount,thisInfoList[info_index]['MpiPlotParams']))
                # thisInfoList[info_index]['sm2ptRat'] = map(int,self.smlist)
                # thisInfoList[info_index]['jsmCoeffs'] = map(np.float64,self.sink_coeffs)
        else:
            thisInfoList[info_index]['sinklen'] = 1
            if 'sm2ptRat' in list(thisInfoList[info_index].keys()):
                del thisInfoList[info_index]['sm2ptRat']

        if self.combine_all_corrs != 'No':
            ## TODO, add coeffs as input parameter
            thisInfoList[info_index]['CombCoffs'] = self.Coeff_of_corr
            thisInfoList[info_index]['CombOpp'] = self.combine_all_corrs

    def _Reformat_fired(self):
        print()
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,Mcore
        print('Begining reformatting of 2 point correlation functions')
        Out2ptRun(thisInfoList,def_Corr_Sets,self.reformat_directory)
        # data = sr.SetOfRat(InfoDict=thisInfoFlowList,corr_cfglists=def_Corr_Sets)
        # data.Reformat3pt(self.reformat_directory)
        # data = None
        print('Reformatting complete')


    def _WriteBoot_fired(self):
        print()
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,Mcore
        print('Begining writing boots of 2 point correlation functions')
        OutBoot2ptRun(thisInfoList,def_Corr_Sets,self.boot_directory)
        # data = sr.SetOfRat(InfoDict=thisInfoFlowList,corr_cfglists=def_Corr_Sets)
        # data.Reformat3pt(self.reformat_directory)
        # data = None
        print('Writing boots complete')

    def _CreateSeparateJob_fired(self):
        print()
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,JobNumber
        SavePreviousRun()
        for iinfo in thisInfoList:
            this_job_file = JobFolder+'TwoPt'+str(JobNumber)+'.py3p'
            IO.WritePickle(this_job_file,[[iinfo],def_Corr_Sets,def_checkcfgs,thisDefWipe])
            WriteCSH('Ana2pt',this_job_file,job_params,JobNumber)
            print('Created job script for Two-point correlator run', this_job_file.replace('.py3p','.csh'))
            JobNumber += 1


    def _CreateJob_fired(self):
        print()
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,JobNumber
        SavePreviousRun()
        this_job_file = JobFolder+'TwoPt'+str(JobNumber)+'.py3p'
        IO.WritePickle(this_job_file,[thisInfoList,def_Corr_Sets,def_checkcfgs,thisDefWipe])
        WriteCSH('Ana2pt',this_job_file,job_params,JobNumber)
        print('Created job script for Two-point correlator run', this_job_file.replace('.py3p','.csh'))
        JobNumber += 1

    def _RUN_fired(self):

        print()
        print('Begining Two-point correlator run')
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,plot_info
        SavePreviousRun()
        data = TPRun(thisInfoList,def_Corr_Sets,def_checkcfgs,thisDefWipe)
        # data = so.SetOfTwoPt(InfoDict=thisInfoList,corr_cfglists=def_Corr_Sets)
        # data.LoadPickle(DefWipe=thisDefWipe,CheckCfgs=def_checkcfgs)
        # data.DoPowers()
        # data.DoCoeff()
        fitlist = [idict['Fits'][0] for idict in thisInfoList]
        momplotlist = [['p'+idict['pmom'][0].replace(' ','')] for idict in thisInfoList]
        if self.do_var_method:
            n_dim = len(thisInfoList)**(1/2)
            data.VariationalMethod(t0=self.var_t0,dt=self.var_dt,symetrize=self.var_symetrize,n_states=n_dim)
            data.PronyMethod(t0=self.var_t0,dt=self.var_dt,symetrize=self.var_symetrize,
                             ism=self.prony_ism,n_states=n_dim)
        for initlist,idict in enumerate(thisInfoList):
            if idict['combine_sinks']:
                sinkcoeffs = []
                these_indicies = []
                for icount in range(idict['sinklen']):
                    these_indicies.append(icount+initlist)
                    sinkcoeffs.append(thisInfoList[icount+initlist]['jsmCoeffs'][0])

                thisname,thisleglab,jsmlab = CreateNameAndLL(data.SetC2,thisInfoList[initlist]['sinklab'],from_index=initlist)
                thisfun = GenSinkFun(list(map(np.float64,sinkcoeffs)))
                print('combining indicies ', ','.join(map(str,these_indicies)))
                data.AppendCCorrs(thisfun,filename=thisname,LegLab=thisleglab,jsm = thisInfoList[initlist]['sinklab'],indicies = these_indicies)
                fitlist += [fitlist[0]]
                momplotlist += [momplotlist[0]]
                break
        data.EffMassPlot(plot_info,fitDict=fitlist,momlist=momplotlist,OnlyVar=self.Only_Var_plot,Phys=self.Physical_units)
        if self.test_subtraction:
            data.ExcitedSubtraction(self.mass_ratio_guess,[str(self.sub_fit_state),
                                    str(self.sub_fit_range),Remove_Empty_Str(list(self.sub_fit_guess))])
            fitlist += [[str(self.sub_fit_state),str(self.sub_fit_range)]]
            momplotlist += [momplotlist[0]]
        elif self.combine_all_corrs != 'No':
            thisname = ('_'+self.combine_all_corrs+'_').join([idict.filename for idict in data.SetC2.values()])
            thisleglab = ('_'+self.combine_all_corrs+'_').join([idict.LegLab for idict in data.SetC2.values()])
            data.AppendCCorrs(ql.Getnpfun(self.combine_all_corrs),filename = thisname, LegLab = thisleglab)
            momplotlist += [momplotlist[0]]
            fitlist += [fitlist[0]]

        # global xmoremin,xmoremax,ymoremin,ymoremax
        data.EffMassPlot(plot_info,fitDict=fitlist,momlist=momplotlist,OnlyVar=self.Only_Var_plot,Phys=self.Physical_units)
        # data.LogPlot(plot_info,momlist=momplotlist,norm=self.norm_plot,OnlyVar=self.Only_Var_plot)
        data.Plot(plot_info,momlist=momplotlist,norm=self.norm_plot,OnlyVar=self.Only_Var_plot)
        input_plot_params = [iplotparam['MpiPlotParams'] for iplotparam in thisInfoList]
        [input_plot_params.insert(*ivar) for ivar in self.var_extra_mpi]
        data.PlotMassVsPionMass(input_plot_params,plot_info,OnlyVar=self.Only_Var_plot)
        data.PlotEnergyVSRecRel(input_plot_params,plot_info,OnlyVar=self.Only_Var_plot,Phys=self.Physical_units)
        data = None
        print('Two-point correlator run complete')


class FlowedOp( ta.HasTraits ):
    """

    Controler for running the Flowed operator (Q or W) analysis

    """

    Calculate_Chi = ta.Bool(True)
    Improved_Overlap = ta.Bool(False)
    Ahmed_comparison = ta.Bool(False)
    Auto_correlation_analysis = ta.Bool(False)
    S_parameter = ta.Float(pa.defSparam)
    tflow_plot_min = ta.Float(pa.defxlimOp[0])
    tflow_plot_max = ta.Float(pa.defxlimOp[-1])
    tflow_plot_inc = ta.Float(pa.defxlimOp[1]-pa.defxlimOp[0])
    observable = ta.Enum(ql.ObsListExt)
    flow_fit_range = ta.Str(pa.defInfoFlow['FlowFits'][0])
    # flow_fit_fun = ta.Str(pa.defInfoFlow['FlowFitFun'])
    tflow_fit = ta.Str(pa.defInfoFlow['tflowfit'][0])
    tsink_fit = ta.Str(pa.defInfoFlow['tsinkfit'][0])
    plot_FO_times_delta = ta.Bool(False)
    epsilon_cut = ta.Float(pa.Qeps)
    # Random_Times = ta.Bool(False)
    # Rand_Times_Per_GF = ta.Str('1-64')
    # randt_fit_max = ta.Str('t63')
    RUN = ta.Button()
    CreateJob = ta.Button()
    CreateSeparateJob = ta.Button()
    Apply = ta.Button()
    Load = ta.Button()

    view = tua.View(
        '_',
        tua.Item('Calculate_Chi'),
        tua.Item('Improved_Overlap'),
        tua.Item('Ahmed_comparison'),
        tua.Item('Auto_correlation_analysis'),
        tua.Item('S_parameter'),
        tua.Item('tflow_plot_min'),
        tua.Item('tflow_plot_max'),
        tua.Item('tflow_plot_inc'),
        tua.Item('observable'),
        tua.Item('flow_fit_range'),
        # tua.Item('flow_fit_fun'),
        tua.Item('tflow_fit'),
        tua.Item('tsink_fit'),
        tua.Item('plot_FO_times_delta'),
        tua.Item('epsilon_cut'),
        # tua.Item('Random_Times'),
        # tua.Item('Rand_Times_Per_GF',enabled_when='Random_Times'),
        # tua.Item('randt_fit_max',enabled_when='Random_Times'),
        # tua.Item('Load',show_label=False),
        tua.Item('Apply',show_label=False),
        tua.Item('CreateJob', show_label=False ),
        tua.Item('CreateSeparateJob', show_label=False ),
        tua.Item('RUN', show_label=False ),
        resizable=True
    )

    def _Load_fired(self):
        self.tflow_plot_min = thisInfoFlowList[info_index]['tflowlist'][0]
        self.tflow_plot_max = thisInfoFlowList[info_index]['tflowlist'][-1]
        if len(thisInfoFlowList[info_index]['tflowlist']) > 1:
            self.tflow_plot_inc = thisInfoFlowList[info_index]['tflowlist'][1] - self.tflow_plot_min
        else:
            self.tflow_plot_inc = 1
        if 'Observable' in list(thisInfoFlowList[info_index].keys()):
            if thisInfoFlowList[info_index]['Observable'] == 'TopCharge':
                self.observable = 'Q'
            if thisInfoFlowList[info_index]['Observable'] == 'Weinberg':
                self.observable = 'W'
            if thisInfoFlowList[info_index]['Observable'] == 'TopCharge_Weinberg':
                self.observable = 'QW'
            if thisInfoFlowList[info_index]['Observable'] == 'Weinberg_Weinberg':
                self.observable = 'WW'
            if thisInfoFlowList[info_index]['Observable'] == 'TopCharge_TopCharge':
                self.observable = 'QQ'
        self.flow_fit_range           = thisInfoFlowList[info_index]['FlowFits'][0]
        self.tsink_fit                = thisInfoFlowList[info_index]['tsinkfit'][0]
        self.Auto_correlation_analysis  = thisInfoFlowList[info_index]['autocorr']
        self.tflow_fit                = thisInfoFlowList[info_index]['tflowfit'][0]
        self.S_parameter                = thisInfoFlowList[info_index]['Sparam']
        if 'IOver' in list(thisInfoFlowList[info_index].keys()):
            self.Improved_Overlap           = thisInfoFlowList[info_index]['IOver']
        if 'deltaEps' in list(thisInfoFlowList[info_index].keys()):
            self.epsilon_cut                = thisInfoFlowList[info_index]['deltaEps']
        if 'Do_Chi' in list(thisInfoFlowList[info_index].keys()):
            self.Calculate_Chi              = thisInfoFlowList[info_index]['Do_Chi']

    def _Apply_fired(self):
        thisInfoFlowList[info_index]['tflowlist'] = np.arange(max(self.tflow_plot_min,0.01),min(self.tflow_plot_max+self.tflow_plot_inc,10),self.tflow_plot_inc)
        if not (1 < self.S_parameter < 2):
            print('Warning, paper suggests choosing 1<S<2')
        if 'TopCharge' in self.observable or 'Q' == self.observable:
            thisInfoFlowList[info_index]['Observable'] = 'TopCharge'
            thisInfoFlowList[info_index]['ObsComb'] = []
        elif 'Weinberg' in self.observable or 'W' == self.observable:
            thisInfoFlowList[info_index]['Observable'] = 'Weinberg'
            thisInfoFlowList[info_index]['ObsComb'] = []
        elif 'Q' in self.observable and 'W' in self.observable:
            thisInfoFlowList[info_index]['Observable'] = 'TopCharge_Weinberg'
            thisInfoFlowList[info_index]['ObsComb'] = ['TopCharge','Weinberg']
        elif 'W' in self.observable:
            thisInfoFlowList[info_index]['Observable'] = 'Weinberg_Weinberg'
            thisInfoFlowList[info_index]['ObsComb'] = ['Weinberg','Weinberg']
        elif 'Q' in self.observable:
            thisInfoFlowList[info_index]['Observable'] = 'TopCharge_TopCharge'
            thisInfoFlowList[info_index]['ObsComb'] = ['TopCharge','TopCharge']
        thisInfoFlowList[info_index]['FlowFits'] = [self.flow_fit_range]
        thisInfoFlowList[info_index]['tsinkfit']= [self.tsink_fit]
        thisInfoFlowList[info_index]['autocorr']= self.Auto_correlation_analysis
        thisInfoFlowList[info_index]['tflowfit']= [self.tflow_fit]
        thisInfoFlowList[info_index]['Sparam'] = self.S_parameter
        thisInfoFlowList[info_index]['IOver'] = self.Improved_Overlap
        thisInfoFlowList[info_index]['deltaEps']= self.epsilon_cut
        thisInfoFlowList[info_index]['Do_Chi'] = self.Calculate_Chi
        # thisInfoFlowList[info_index]['FlowOpRandTimes'] = self.Random_Times
        # if self.Random_Times:
        #     thisInfoFlowList[info_index]['RandTFitM']= self.randt_fit_max
        #     try:
        #         thisInfoFlowList[info_index]['nrandt'] = [int(self.Rand_Times_Per_GF)]
        #     except Exception as err:
        #         thisInfoFlowList[info_index]['nrandt'] = str(self.Rand_Times_Per_GF)
        # else:
            # thisInfoFlowList[info_index]['deltaEps']= pa.Qeps
            # thisInfoFlowList[info_index]['RandTFitM']= 't63'

    def _CreateSeparateJob_fired(self):
        print()
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,JobNumber
        SavePreviousRun()
        for iinfo in thisInfoFlowList:
            this_job_file = JobFolder+'FOp'+str(JobNumber)+'.py3p'
            IO.WritePickle(this_job_file,[[iinfo],bool(self.Ahmed_comparison),def_Corr_Sets,def_checkcfgs,thisDefWipe])
            WriteCSH('AnaFOp',this_job_file,job_params,JobNumber)
            print('Created job script for Flowed Operator run', this_job_file.replace('.py3p','.csh'))
            JobNumber += 1


    def _CreateJob_fired(self):
        print()
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,JobNumber
        SavePreviousRun()
        this_job_file = JobFolder+'FOp'+str(JobNumber)+'.py3p'
        IO.WritePickle(this_job_file,[thisInfoFlowList,bool(self.Ahmed_comparison),def_Corr_Sets,def_checkcfgs,thisDefWipe])
        WriteCSH('AnaFOp',this_job_file,job_params,JobNumber)
        print('Created job script for Flowed Operator run', this_job_file.replace('.py3p','.csh'))
        JobNumber += 1

    def _RUN_fired(self):
        print()
        print('Begining Flow run')
        fitlist = [str(idict['FlowFits'][0]) for idict in thisInfoFlowList]
        global thisDefWipe,def_checkcfgs,def_Corr_Sets
        SavePreviousRun()
        data = FOpRun(thisInfoFlowList,bool(self.Ahmed_comparison),def_Corr_Sets,def_checkcfgs,thisDefWipe)
        # data = sfo.SetOfFO(InfoDict=thisInfoFlowList,CompAhmed=bool(self.Ahmed_comparison),corr_cfglists=def_Corr_Sets)
        # data.LoadPickle(DefWipe=thisDefWipe,CheckCfgs=def_checkcfgs)
        # data.FlowMakeChit(DefWipe=thisDefWipe)
        # randtfitmax = [itflow['RandTFitM'] for itflow in thisInfoFlowList]
        tflowpicked = [itflow['tflowfit'][-1] for itflow in thisInfoFlowList]
        epslist = [itflow['deltaEps'] for itflow in thisInfoFlowList]
        global plot_info
        data.PlotMonteTime(tflowpicked,plot_info)
        data.PlotMonteTimeDeltaEps(tflowpicked,plot_info,epslist,MultQ=bool(self.plot_FO_times_delta))
        data.Plot(plot_info,fitlist)
        data.PlotTauIntOverFlow(plot_info)
        data.PlotTauInt(tflowpicked,plot_info)
        data.PlotWopt(tflowpicked,plot_info)
        # data.PlotRandTimes(tflowpicked,randt_fitr=randtfitmax,tsize=tsize,xsize=xsize,ysize=ysize,ftitle=thetitle,fxlab=thexlab,fylab=theylab,legloc=legloc)
        data.PlotVsPionMass(fitlist,plot_info)
        data.PlotVsQuarkMass(tflowpicked,plot_info)
        data = None
        print('Flow run complete')



class FlowedOpFull( ta.HasTraits ):
    """

    Controler for running the code for analysing the Euclidean time dependnace of
    the Flowed operators.

    """

    Calculate_Chi = ta.Bool(True)
    Only_Q = ta.Bool(False)
    Auto_correlation_analysis = ta.Bool(False)
    S_parameter = ta.Float(pa.defSparam)
    Sub_Lattice_Method = ta.Bool(False)
    Shift_Corr = ta.Bool(False)
    Number_of_Random_Times = ta.Int(1)
    Sum_Tsink = ta.Bool(True)
    Op_fold = ta.Bool(True)
    Number_of_Sink_Times = ta.Int(33)
    Not_Only_Picked_Tflow = ta.Bool(False)
    tflow_read_min = ta.Float(pa.defxlimOp[0])
    tflow_read_max = ta.Float(pa.defxlimOp[-1])
    tflow_read_inc = ta.Float(pa.defxlimOp[1]-pa.defxlimOp[0])
    Not_Sum_Whole_Lat = ta.Bool(True)
    tsum_read_All = ta.Bool(True)
    not_tsum_read_All = ta.Bool(not bool(tsum_read_All))
    tsum_read_min = ta.Int(63)
    tsum_read_max = ta.Int(63)
    observable = ta.Enum(ql.ObsList)
    tflow_picked = ta.Str(pa.defInfoFlow['tflowfit'][0])
    fit_par_picked = ta.Str('First')
    tsink_picked = ta.Str('tr32')
    tsum_picked = ta.Str('ts63')
    t_fit_range = ta.Str('fitr2-8')
    do_multiple_fits = ta.Bool(False)
    t_fit_min_min = ta.Int(2)
    t_fit_min_max = ta.Int(8)
    t_fit_max_min = ta.Int(8)
    t_fit_max_max = ta.Int(20)
    fit_min_picked = ta.Int(2)
    fit_max_picked = ta.Int(20)
    RUN = ta.Button()
    CreateJob = ta.Button()
    CreateSeparateJob = ta.Button()
    Apply = ta.Button()
    Load = ta.Button()

    view = tua.View(
        tua.Group(
        tua.Item('Only_Q'),
        tua.Item('Calculate_Chi'),
        tua.Item('Auto_correlation_analysis'),
        tua.Item('S_parameter'),
        tua.Item('Sub_Lattice_Method'),
        tua.Item('Shift_Corr'),
        tua.Item('Not_Sum_Whole_Lat'),
        tua.Item('Number_of_Random_Times',enabled_when='Not_Sum_Whole_Lat'),
        tua.Item('Number_of_Sink_Times',enabled_when='Not_Sum_Whole_Lat'),
        tua.Item('Sum_Tsink'),
        tua.Item('Op_fold'),
        tua.Item('Not_Only_Picked_Tflow'),
        tua.Item('tflow_read_min',enabled_when='Not_Only_Picked_Tflow'),
        tua.Item('tflow_read_max',enabled_when='Not_Only_Picked_Tflow'),
        tua.Item('tflow_read_inc',enabled_when='Not_Only_Picked_Tflow'),
        tua.Item('tsum_read_All',enabled_when='Not_Sum_Whole_Lat'),
        tua.Item('tsum_read_min',enabled_when='not_tsum_read_All'),
        tua.Item('tsum_read_max',enabled_when='not_tsum_read_All'),
        tua.Item('observable'),
        tua.Item('tflow_picked'),
        tua.Item('fit_par_picked'),
        tua.Item('tsink_picked',enabled_when='Not_Sum_Whole_Lat'),
        tua.Item('tsum_picked',enabled_when='Not_Sum_Whole_Lat'),
        tua.Item('t_fit_range'),
        tua.Item('Apply',show_label=False),
        tua.Item('CreateJob', show_label=False ),
        tua.Item('CreateSeparateJob', show_label=False ),
        tua.Item('RUN', show_label=False ),
        label='Page 1'
        ),tua.Group(
        tua.Item('do_multiple_fits'),
        tua.Item('t_fit_min_min',enabled_when='do_multiple_fits'),
        tua.Item('t_fit_min_max',enabled_when='do_multiple_fits'),
        tua.Item('t_fit_max_min',enabled_when='do_multiple_fits'),
        tua.Item('t_fit_max_max',enabled_when='do_multiple_fits'),
        tua.Item('fit_min_picked',enabled_when='do_multiple_fits'),
        tua.Item('fit_max_picked',enabled_when='do_multiple_fits'),
        # tua.Item('Load',show_label=False),
        tua.Item('Apply',show_label=False),
        tua.Item('CreateJob', show_label=False ),
        tua.Item('CreateSeparateJob', show_label=False ),
        tua.Item('RUN', show_label=False ),
        label='Mult_Fits'
        ),
        resizable=True
    )

    def _tsum_read_All_fired(self):
        self.not_tsum_read_All = not self.tsum_read_All

    def _Load_fired(self):
        if len(thisInfoFlowList[info_index]['tflowlist']) == 1:
            self.Not_Only_Picked_Tflow = True
        else:
            self.Not_Only_Picked_Tflow = False
            self.tflow_read_min = thisInfoFlowList[info_index]['tflowlist'][0]
            self.tflow_read_inc = thisInfoFlowList[info_index]['tflowlist'][1] - self.tflow_read_min
            self.tflow_read_max = thisInfoFlowList[info_index]['tflowlist'][-1]

        self.tsum_read_All = str(thisInfoFlowList[info_index]['tsum_list']) == 'All'
        if not self.tsum_read_All :
            self.tsum_read_min = int(thisInfoFlowList[info_index]['tsum_list'][0].replace('ts',''))
            self.tsum_read_max = int(thisInfoFlowList[info_index]['tsum_list'][-1].replace('ts',''))
        if not (1 < self.S_parameter < 2):
            print('Warning, paper suggests choosing 1<S<2')


        if 'tsum_picked' in list(thisInfoFlowList[info_index].keys()):
            if isinstance(thisInfoFlowList[info_index]['tsum_picked'],(tuple,list,np.ndarray)):
                self.tsum_picked  = thisInfoFlowList[info_index]['tsum_picked'][0]
            else:
                self.tsum_picked  = thisInfoFlowList[info_index]['tsum_picked']
        if 'nrandt' in list(thisInfoFlowList[info_index].keys()):
            self.Number_of_Random_Times=  thisInfoFlowList[info_index]['nrandt']
        if 'tsink_picked' in list(thisInfoFlowList[info_index].keys()):
            self.tsink_picked=                 thisInfoFlowList[info_index]['tsink_picked']
        if 'nranget' in list(thisInfoFlowList[info_index].keys()):
            self.Number_of_Sink_Times=    thisInfoFlowList[info_index]['nranget']
        if 'autocorr' in list(thisInfoFlowList[info_index].keys()):
            self.Auto_correlation_analysis=    thisInfoFlowList[info_index]['autocorr']
        if 'Do_Chi' in list(thisInfoFlowList[info_index].keys()):
            self.Calculate_Chi=                thisInfoFlowList[info_index]['Do_Chi']
        if 'Sparam' in list(thisInfoFlowList[info_index].keys()):
            self.S_parameter=           thisInfoFlowList[info_index]['Sparam']
        if 'Observable' in list(thisInfoFlowList[info_index].keys()):
            if 'TopCharge' in thisInfoFlowList[info_index]['Observable']:
                self.observable=                   'TopCharge'
            elif 'Weinberg' in thisInfoFlowList[info_index]['Observable']:
                self.observable=                   'Weinberg'
        if 'Sum_Range' in list(thisInfoFlowList[info_index].keys()):
            self.Sum_Tsink=                    thisInfoFlowList[info_index]['Sum_Range']
        if 'Op_fold' in list(thisInfoFlowList[info_index].keys()):
            self.Op_fold =                    thisInfoFlowList[info_index]['Op_fold']
        if 'tflow_picked' in list(thisInfoFlowList[info_index].keys()):
            self.tflow_picked=                 thisInfoFlowList[info_index]['tflow_picked']
        # if 'fit_min_picked' in thisInfoFlowList[info_index].keys():
        #     self.do_multiple_fits = thisInfoFlowList[info_index]['fit_min_picked'] == float('NaN')
        if self.do_multiple_fits:
            split_list = [list(map(int,ix.replace('fitr','').split('-'))) for ix in thisInfoFlowList[info_index]['t_fit_range']]
            start_list,fin_list = list(zip(*split_list))
            self.t_fit_min_min = np.min(start_list)
            self.t_fit_min_max = np.max(start_list)
            self.t_fit_max_min = np.min(fin_list)
            self.t_fit_max_max = np.max(fin_list)
            self.fit_min_picked =    thisInfoFlowList[info_index]['fit_min_picked']
            self.fit_max_picked =    thisInfoFlowList[info_index]['fit_max_picked']
        else:
            if 't_fit_range' in list(thisInfoFlowList[info_index].keys()):
                self.t_fit_range = thisInfoFlowList[info_index]['t_fit_range'][0]
        if 'fit_par_picked' in list(thisInfoFlowList[info_index].keys()):
            self.fit_par_picked = thisInfoFlowList[info_index]['fit_par_picked']

    def _Apply_fired(self):
        if self.Not_Only_Picked_Tflow:
            thisInfoFlowList[info_index]['tflowlist'] = np.arange(max(self.tflow_read_min,0.01),min(self.tflow_read_max+self.tflow_read_inc,10),self.tflow_read_inc)
        else:
            thisInfoFlowList[info_index]['tflowlist'] = np.array([float(self.tflow_picked.replace('t_f',''))])

        if self.Sub_Lattice_Method:
            thisInfoFlowList[info_index]['Sublat_Method'] =  True
            thisInfoFlowList[info_index]['tsum_list'] =  ['ts0']
            thisInfoFlowList[info_index]['tsum_picked'] = 'ts0'
            thisInfoFlowList[info_index]['nrandt']= 64
            thisInfoFlowList[info_index]['tsink_picked'] = 'tr63'
            thisInfoFlowList[info_index]['nranget']= 64
        elif not self.Not_Sum_Whole_Lat:
            thisInfoFlowList[info_index]['tsum_list'] =  ['ts63']
            thisInfoFlowList[info_index]['tsum_picked'] = 'ts63'
            thisInfoFlowList[info_index]['nrandt']= 1
            thisInfoFlowList[info_index]['tsink_picked'] = 'tr32'
            thisInfoFlowList[info_index]['nranget']= 33
        elif self.tsum_read_All:
            thisInfoFlowList[info_index]['tsum_list'] =  'All'
            thisInfoFlowList[info_index]['tsum_picked'] = self.tsum_picked
            thisInfoFlowList[info_index]['nrandt']= int(self.Number_of_Random_Times)
            thisInfoFlowList[info_index]['tsink_picked'] = self.tsink_picked
            thisInfoFlowList[info_index]['nranget']= int(self.Number_of_Sink_Times)
        else:
            thisInfoFlowList[info_index]['tsum_list'] =  np.arange(max(self.tsum_read_min,0),min(self.tsum_read_max+1,64))
            thisInfoFlowList[info_index]['tsum_list'] = ['ts'+str(it) for it in thisInfoFlowList[info_index]['tsum_list']]
            thisInfoFlowList[info_index]['tsum_picked'] = self.tsum_picked
            thisInfoFlowList[info_index]['nrandt']= int(self.Number_of_Random_Times)
            thisInfoFlowList[info_index]['tsink_picked'] = self.tsink_picked
            thisInfoFlowList[info_index]['nranget']= int(self.Number_of_Sink_Times)
        thisInfoFlowList[info_index]['autocorr']= self.Auto_correlation_analysis
        if not (1 < self.S_parameter < 2):
            print('Warning, paper suggests choosing 1<S<2')
        thisInfoFlowList[info_index]['Do_Chi'] = self.Calculate_Chi
        thisInfoFlowList[info_index]['Sparam'] = float(self.S_parameter)
        thisInfoFlowList[info_index]['Observable'] = self.observable
        thisInfoFlowList[info_index]['Sum_Range'] = self.Sum_Tsink
        thisInfoFlowList[info_index]['Op_fold'] = self.Op_fold
        thisInfoFlowList[info_index]['tflow_picked'] = self.tflow_picked
        if self.do_multiple_fits:
            thisInfoFlowList[info_index]['t_fit_range'] = []
            for fitmin in range(self.t_fit_min_min,self.t_fit_min_max+1):
                for fitmax in range(self.t_fit_max_min,self.t_fit_max_max+1):
                    thisInfoFlowList[info_index]['t_fit_range'].append('fitr'+str(fitmin)+'-'+str(fitmax))
            thisInfoFlowList[info_index]['fit_min_picked'] = str(self.fit_min_picked)
            thisInfoFlowList[info_index]['fit_max_picked'] = str(self.fit_max_picked)
        else:
            thisInfoFlowList[info_index]['t_fit_range'] = [str(self.t_fit_range)]
            thisInfoFlowList[info_index]['fit_min_picked'] = str('NaN')
            thisInfoFlowList[info_index]['fit_max_picked'] = str('NaN')
        thisInfoFlowList[info_index]['fit_par_picked'] = str(self.fit_par_picked)

    def _CreateSeparateJob_fired(self):
        print()
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,JobNumber,plot_info
        SavePreviousRun()
        for iinfo in thisInfoFlowList:
            this_job_file = JobFolder+'FOpFull'+str(JobNumber)+'.py3p'
            IO.WritePickle(this_job_file,[[iinfo],def_Corr_Sets,def_checkcfgs,thisDefWipe,self.Only_Q])
            WriteCSH('AnaFOpFull',this_job_file,job_params,JobNumber)
            print('Created job script for Flowed Operator Full Analysis run', this_job_file.replace('.py3p','.csh'))
            JobNumber += 1


    def _CreateJob_fired(self):
        print()
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,JobNumber,plot_info
        SavePreviousRun()
        this_job_file = JobFolder+'FOpFull'+str(JobNumber)+'.py3p'
        IO.WritePickle(this_job_file,[thisInfoFlowList,def_Corr_Sets,def_checkcfgs,thisDefWipe,self.Only_Q])
        WriteCSH('AnaFOpFull',this_job_file,job_params,JobNumber)
        print('Created job script for Flowed Operator Full Analysis run', this_job_file.replace('.py3p','.csh'))
        JobNumber += 1

    def _RUN_fired(self):
        print()
        print('Begining Flow Full run')
        global thisDefWipe,def_checkcfgs,def_Corr_Sets
        SavePreviousRun()
        data = FOpFullRun(thisInfoFlowList,def_Corr_Sets,def_checkcfgs,thisDefWipe,self.Only_Q)
        # data = sfo.SetOfFOFull(InfoDict=thisInfoFlowList,corr_cfglists=def_Corr_Sets)
        # data.LoadPickle(DefWipe=thisDefWipe,CheckCfgs=def_checkcfgs)
        tflow_picked = [iinfo['tflow_picked'] for iinfo in thisInfoFlowList]
        tsink_picked = [iinfo['tsink_picked'] for iinfo in thisInfoFlowList]
        tsum_picked = [iinfo['tsum_picked'] for iinfo in thisInfoFlowList]
        t_fit_range_picked = [iinfo['t_fit_range'] for iinfo in thisInfoFlowList]
        fit_par_picked = [iinfo['fit_par_picked'] for iinfo in thisInfoFlowList]
        fit_min_picked = [iinfo['fit_min_picked'] for iinfo in thisInfoFlowList]
        fit_max_picked = [iinfo['fit_max_picked'] for iinfo in thisInfoFlowList]
        if self.Only_Q:
            data.PlotOverT(tflow_picked,plot_info)
        else:
            data.PlotOverTsink(tflow_picked,tsum_picked,plot_info,fitrlist=t_fit_range_picked,shifted=self.Shift_Corr)
            data.PlotTauIntOverTsink(tflow_picked,tsum_picked,plot_info)
            data.PlotOverTsum(tflow_picked,tsink_picked,plot_info)
            data.PlotTauIntOverTsum(tflow_picked,tsink_picked,plot_info)
            # data.PlotOverTflow(tsum_picked,tsink_picked,plot_info)
            data.PlotFitOverTflow(t_fit_range_picked,tsum_picked,fit_par_picked,plot_info)
            data.PlotTauIntOverTflow(tsum_picked,tsink_picked,plot_info)
            data.PlotTauInt(tflow_picked,tsink_picked,tsum_picked,plot_info)
            data.PlotWopt(tflow_picked,tsink_picked,tsum_picked,plot_info)
            data.PlotVsPionMass(tflow_picked,t_fit_range_picked,tsum_picked,fit_par_picked,plot_info)
            if bool(self.do_multiple_fits):
                data.PlotVstfFitr(tflow_picked,tsum_picked,fit_min_picked,plot_info,parlist=fit_par_picked)
                data.PlotVstiFitr(tflow_picked,tsum_picked,fit_max_picked,plot_info,parlist=fit_par_picked)
        data = None
        print('Flow run complete')







class RatioFun( ta.HasTraits ):
    """

    Controler for the Ratio function analysis of the three- and two-point correlation functions

    """
    reformat_directory = ta.Str(pa.outputdir+'/Reformat3pt/')
    Rat_alpha = ta.Bool(True)
    source_smearing = ta.Int(pa.defInfo['ism'])
    t_source = ta.Int(pa.defInfo['tsrc'])
    combine_sinks = ta.Bool('sm2ptRat' in list(pa.defInfo.keys()))
    not_combine_sinks = ta.Bool(not ('sm2ptRat' in list(pa.defInfo.keys())))
    sink_smearing = ta.Int(pa.defInfo['jsm'])
    smlist = list(map(str,pa.defInfo['sm2ptRat']))
    smear_list = ta.Str(', '.join(['sm'+ival for ival in smlist]))
    smear_comment = ta.Str(' Set coeff to 0 to ignore smearing ')
    sink_coeffs = ta.Array(np.float64,(len(smlist),),pa.defInfo['jsmCoeffs'])
    sink_label = ta.Str(pa.defInfoFlow['sinklab'])
    rf_fit_range = ta.Str(pa.defInfoFlow['FitsRF'][0])
    current_momentum_list = ta.List(ta.Str,pa.defInfo['qmom'])
    t_sink = ta.Int(pa.defInfo['tsink'])
    sink_momentum = ta.Str(pa.defInfo['ppmom'])
    doublet_singlet_comb = ta.Enum(ql.DSList)
    projector = ta.Enum(ql.ProjList)
    current_gamma_list = ta.List(ta.Enum(ql.all_gammas),pa.defInfo['gammas'])
    # pformlist = ['p'+imom.replace(' ','') for imom in momentum_list]

    njnq_observable = ta.Bool(False)
    tflow_plot_min = ta.Float(pa.defxlimOp[0])
    tflow_plot_max = ta.Float(pa.defxlimOp[-1])
    tflow_plot_inc = ta.Float(pa.defxlimOp[1]-pa.defxlimOp[0])
    observable = ta.Enum(ql.ObsList)
    tflow_fit = ta.Str(pa.defInfoFlow['tflowfit'][0])
    Rat_tsum_range = ta.Str(pa.defInfoFlow['FitsAlphaFull'][0][1])
    Rat_sum_type = ta.Enum(ql.rat_sum_list)
    Rat_comb_alpha = ta.Bool(False)
    Rat_comb_alpha_fun = ta.Enum(op_str_list)
    frat_bool = ta.Bool(False)
    full_ratios = ta.Bool(False)
    tsum_picked = ta.Str(pa.defInfoFlow['tsum_fit_list'][0])
    tcurr_picked = ta.Str(tstr(untstr(pa.defInfo['tsink'])//2))
    RUN = ta.Button()
    CreateJob = ta.Button()
    CreateSeparateJob = ta.Button()
    Reformat = ta.Button()
    Apply = ta.Button()
    Load = ta.Button()


    view = tua.View(tua.Group(
        tua.Item('reformat_directory'),
        tua.Item('Rat_alpha'),
        tua.Item('source_smearing'),
        tua.Item('t_source'),
        tua.Item('combine_sinks'),
        tua.Item('sink_smearing',enabled_when='not_combine_sinks'),
        tua.Item('smear_list',enabled_when='combine_sinks',style='readonly'),
        tua.Item('smear_comment',show_label=False,enabled_when='combine_sinks',style='readonly'),
        tua.Item('sink_coeffs',enabled_when='combine_sinks',style='custom'),
        tua.Item('sink_label',enabled_when='combine_sinks'),
        tua.Item('rf_fit_range'),
        tua.Item('current_momentum_list'),
        tua.Item('t_sink'),
        tua.Item('sink_momentum'),
        tua.Item('doublet_singlet_comb'),
        tua.Item('projector'),
        tua.Item('current_gamma_list'),
        tua.Item('Apply', show_label=False),
        tua.Item('CreateJob', show_label=False ),
        tua.Item('CreateSeparateJob', show_label=False ),
        tua.Item('Reformat',show_label=False),
        tua.Item('RUN', show_label=False ),
        label='Page 1'
        ),tua.Group(
        tua.Item('njnq_observable'),
        tua.Item('tflow_plot_min',enabled_when='njnq_observable'),
        tua.Item('tflow_plot_max',enabled_when='njnq_observable'),
        tua.Item('tflow_plot_inc',enabled_when='njnq_observable'),
        tua.Item('observable',enabled_when='njnq_observable'),
        tua.Item('tflow_fit',enabled_when='njnq_observable'),
        tua.Item('Rat_comb_alpha',enabled_when='njnq_observable'),
        tua.Item('Rat_comb_alpha_fun',enabled_when='Rat_comb_alpha'),
        tua.Item('full_ratios',enabled_when='njnq_observable'),
        tua.Item('Rat_sum_type',enabled_when='frat_bool'),
        tua.Item('Rat_tsum_range',enabled_when='frat_bool'),
        tua.Item('tsum_picked',enabled_when='frat_bool'),
        tua.Item('tcurr_picked',enabled_when='frat_bool'),
        tua.Item('Apply', show_label=False),
        tua.Item('CreateJob', show_label=False ),
        tua.Item('CreateSeparateJob', show_label=False ),
        tua.Item('Reformat',show_label=False),
        tua.Item('RUN', show_label=False ),
    label='Page 2'),        resizable=True)


    def _Load_fired(self):
        self.combine_sinks = 'sm2ptRat' in list(thisInfoFlowList[info_index].keys())
        try:
            if self.combine_sinks:
                self.smlist = list(map(np.float64,thisInfoFlowList[info_index]['sm2ptRat']))
                self.sink_coeffs = list(map(np.float64,thisInfoFlowList[info_index]['jsmCoeffs']))
                self.sink_label = thisInfoFlowList[info_index]['sinklab']
        except Exception as err:
            print(str(err))
        if 'Rat_alpha' in list(thisInfoFlowList[info_index].keys()):
            self.Rat_alpha     =   thisInfoFlowList[info_index]['Rat_alpha']
        self.source_smearing     =   thisInfoFlowList[info_index]['ism']
        self.t_source          =      thisInfoFlowList[info_index]['tsrc']
        self.sink_smearing     =   thisInfoFlowList[info_index]['jsm']
        self.current_momentum_list   =   thisInfoFlowList[info_index]['qmom']
        self.t_sink      =      int(thisInfoFlowList[info_index]['tsink'])
        self.sink_momentum     =         thisInfoFlowList[info_index]['ppmom']
        self.doublet_singlet_comb     =       thisInfoFlowList[info_index]['DoubSing']
        self.projector          =          thisInfoFlowList[info_index]['Projector']
        self.current_gamma_list  = thisInfoFlowList[info_index]['gammas']
        if 'Full_Rat' in list(thisInfoFlowList[info_index].keys()):
            self.frat_bool         =     thisInfoFlowList[info_index]['Full_Rat']
        else:
            self.frat_bool = False
        if self.frat_bool:
            if isinstance(thisInfoFlowList[info_index]['op_trange_sum_type_3pt'],tuple):
                thisop = '_'.join(thisInfoFlowList[info_index]['op_trange_sum_type_3pt'])
            else:
                thisop = thisInfoFlowList[info_index]['op_trange_sum_type_3pt']
            if '_none' in thisop.lower():
                self.Rat_sum_type = '_'.join(thisop.split('_')[:-1])
            else:
                self.Rat_sum_type = thisop
            if isinstance(thisInfoFlowList[info_index]['tsum_picked'],(tuple,list,np.ndarray)):
                self.tsum_picked  = thisInfoFlowList[info_index]['tsum_picked'][0]
            else:
                self.tsum_picked  = thisInfoFlowList[info_index]['tsum_picked']
            if 'tcurr_picked' in thisInfoFlowList[info_index].keys():
                self.tcurr_picked = thisInfoFlowList[info_index]['tcurr_picked'][0]
            self.rf_fit_range,self.Rat_tsum_range = thisInfoFlowList[info_index]['FitsRFFull'][0][0:2]
        else:
            self.rf_fit_range = thisInfoFlowList[info_index]['FitsRF'][0]

        self.qformlist = [str('q'+imom.replace(' ','')) for imom in Remove_Empty_Str(self.current_momentum_list)]
        self.ppform = str('pp'+self.sink_momentum.replace(' ',''))

        self.njnq_observable = 'Observable' in list(thisInfoFlowList[info_index].keys())
        if self.njnq_observable:
            if 'Rat_comb_alpha_fun' in thisInfoFlowList[info_index]:
                if thisInfoFlowList[info_index]['Rat_comb_alpha_fun'] is False:
                    self.Rat_comb_alpha = False
                else:
                    self.Rat_comb_alpha = True
                    if thisInfoFlowList[info_index]['Rat_comb_alpha_fun'] in op_str_list:
                        self.Rat_comb_alpha_fun = thisInfoFlowList[info_index]['Rat_comb_alpha_fun']
            self.tflow_plot_min = thisInfoFlowList[info_index]['tflowlist'][0]
            self.tflow_plot_max = thisInfoFlowList[info_index]['tflowlist'][-1]
            if len(thisInfoFlowList[info_index]['tflowlist']) == 1:
                self.tflow_plot_inc = 1
            else:
                self.tflow_plot_inc = thisInfoFlowList[info_index]['tflowlist'][1] - self.tflow_plot_min
            if 'TopCharge' in thisInfoFlowList[info_index]['Observable']:
                self.observable=                   'TopCharge'
            elif 'Weinberg' in thisInfoFlowList[info_index]['Observable']:
                self.observable=                   'Weinberg'
            self.tflow_fit = thisInfoFlowList[info_index]['tflowfit'][0]

    def _njnq_observable_changed(self):
        self.frat_bool = self.full_ratios and self.njnq_observable
        if self.frat_bool:
            self.Rat_comb_alpha = False

    def _full_ratios_changed(self):
        self.frat_bool = self.full_ratios and self.njnq_observable
        if self.frat_bool:
            self.Rat_comb_alpha = False

    def _combine_sinks_changed(self):
        self.not_combine_sinks = not self.combine_sinks

    def _Reformat_fired(self):
        print()
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,Mcore
        if any([iflowlist['DoubSing'] not in ['doub','sing'] for iflowlist in thisInfoFlowList]):
            print('Reformatting only supported for doublet and singlet quantities, please select doub or sing')
        else:
            print('Begining reformatting of 3 point correlation functions')
            Out3ptRun(thisInfoFlowList,def_Corr_Sets,self.reformat_directory)
            # data = sr.SetOfRat(InfoDict=thisInfoFlowList,corr_cfglists=def_Corr_Sets)
            # data.Reformat3pt(self.reformat_directory)
            # data = None
            print('Reformatting complete complete')


    def _Apply_fired(self):
        thisInfoFlowList[info_index]['Rat_alpha'] = bool(self.Rat_alpha)
        thisInfoFlowList[info_index]['ism'] = int(self.source_smearing)
        thisInfoFlowList[info_index]['tsrc'] = self.t_source
        thisInfoFlowList[info_index]['jsm'] = self.sink_smearing
        if self.combine_sinks:
            thisInfoFlowList[info_index]['sm2ptRat'] = list(map(int,self.smlist))
            thisInfoFlowList[info_index]['jsmCoeffs'] = list(map(np.float64,self.sink_coeffs))
            thisInfoFlowList[info_index]['sinklab'] = self.sink_label
        else:
            if 'sm2ptRat' in list(thisInfoFlowList[info_index].keys()): del thisInfoFlowList[info_index]['sm2ptRat']

        thisInfoFlowList[info_index]['qmom'] = Remove_Empty_Str(list(map(str,self.current_momentum_list)))
        thisInfoFlowList[info_index]['tsink'] = str(self.t_sink)
        thisInfoFlowList[info_index]['ppmom'] = str(self.sink_momentum)
        thisInfoFlowList[info_index]['DoubSing'] = str(self.doublet_singlet_comb)
        thisInfoFlowList[info_index]['Projector'] = str(self.projector)
        thisInfoFlowList[info_index]['gammas'] = Remove_Empty_Str(list(map(str,self.current_gamma_list)))
        thisInfoFlowList[info_index]['Full_Rat'] = bool(self.frat_bool)
        if self.frat_bool:
            thisInfoFlowList[info_index]['op_trange_sum_type_3pt'] = str(self.Rat_sum_type)+'_None'
            if 'sym' in str(self.Rat_sum_type):
                thisInfoFlowList[info_index]['boot_tsum_list'] = list(range(thisInfoFlowList[info_index]['nxyzt'][-1]//2))
                thisInfoFlowList[info_index]['tsum_fit_list'] = list(range(thisInfoFlowList[info_index]['nxyzt'][-1]//2))
            else:
                thisInfoFlowList[info_index]['boot_tsum_list'] = list(range(thisInfoFlowList[info_index]['nxyzt'][-1]))
                thisInfoFlowList[info_index]['tsum_fit_list'] = list(range(thisInfoFlowList[info_index]['nxyzt'][-1]))
            thisInfoFlowList[info_index]['tsum_picked'] = [str(self.tsum_picked)]
            thisInfoFlowList[info_index]['tcurr_picked'] = [str(self.tcurr_picked)]
            thisInfoFlowList[info_index]['FitsRFFull'] = [(str(self.rf_fit_range),str(self.Rat_tsum_range))]
        else:
            thisInfoFlowList[info_index]['op_trange_sum_type_3pt'] = 'None'
            thisInfoFlowList[info_index]['tsum_picked'] = ['ts32']
            thisInfoFlowList[info_index]['tcurr_picked'] = ['t6']
            thisInfoFlowList[info_index]['FitsRF'] = [str(self.rf_fit_range)]

        thisInfoFlowList[info_index]['Interp'] = 'CPEven'
        thisInfoFlowList[info_index]['MesOrBar'] = 'Baryon'


        self.qformlist = [str('q'+imom.replace(' ','')) for imom in Remove_Empty_Str(self.current_momentum_list)]
        self.ppform = str('pp'+self.sink_momentum.replace(' ',''))

        if self.njnq_observable:
            if self.Rat_comb_alpha:
                thisInfoFlowList[info_index]['Rat_comb_alpha_fun'] = self.Rat_comb_alpha_fun
            else:
                thisInfoFlowList[info_index]['Rat_comb_alpha_fun'] = 'None'
            thisInfoFlowList[info_index]['tflowlist'] = np.arange(max(self.tflow_plot_min,0.01),min(self.tflow_plot_max+self.tflow_plot_inc,10),self.tflow_plot_inc)
            thisInfoFlowList[info_index]['Observable'] = self.observable
            thisInfoFlowList[info_index]['tflowfit']= [str(self.tflow_fit)]
        else:
            if 'Observable' in list(thisInfoFlowList[info_index].keys()):
                del thisInfoFlowList[info_index]['Observable']

    def _CreateSeparateJob_fired(self):
        print()
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,Mcore,JobNumber
        SavePreviousRun()
        for iinfo in thisInfoFlowList:
            this_job_file = JobFolder+'Rat'+str(JobNumber)+'.py3p'
            IO.WritePickle(this_job_file,[[iinfo],def_Corr_Sets,def_checkcfgs,thisDefWipe,Mcore])
            WriteCSH('AnaRat',this_job_file,job_params,JobNumber)
            print('Created job script for Ratio Function run', this_job_file.replace('.py3p','.csh'))
            JobNumber += 1


    def _CreateJob_fired(self):
        print()
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,Mcore,JobNumber
        SavePreviousRun()
        this_job_file = JobFolder+'Rat'+str(JobNumber)+'.py3p'
        IO.WritePickle(this_job_file,[thisInfoFlowList,def_Corr_Sets,def_checkcfgs,thisDefWipe,Mcore])
        WriteCSH('AnaRat',this_job_file,job_params,JobNumber)
        print('Created job script for Ratio Function run', this_job_file.replace('.py3p','.csh'))
        JobNumber += 1

    def _RUN_fired(self):
        # for ikey,ival in thisInfoFlowList[0].iteritems():
        #     print ikey,ival
        print()
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,Mcore,plot_info
        print('Begining Ratio correlator run')
        SavePreviousRun()
        data = RRun(thisInfoFlowList,def_Corr_Sets,def_checkcfgs,thisDefWipe,Mcore)
        # data = sr.SetOfRat(InfoDict=thisInfoFlowList,corr_cfglists=def_Corr_Sets)
        # data.LoadPickle(DefWipe=thisDefWipe,CheckCfgs=def_checkcfgs,Mcore=Mcore)
        fitlist = [str(idict['FitsRF'][0]) for idict in thisInfoFlowList]
        fitlistFull = [str('_'.join(idict['FitsRFFull'][0])) for idict in thisInfoFlowList]
        tflowlist = [str(idict['tflowfit'][0]) for idict in thisInfoFlowList]
        tsumlist = [str(idict['tsum_picked'][0]) for idict in thisInfoFlowList]
        gammamomlist = [{idict['gammas'][0]:[idict['qmom'][0]]} for idict in thisInfoFlowList]
        data.Plot(  plot_info,fitDict=fitlist,fitDictFull=fitlistFull,
                    tflowlist=tflowlist,tsum_list=tsumlist,gammamomlist=gammamomlist)
        if self.full_ratios:
            tcurrlist = [str(idict['tcurr_picked'][0]) for idict in thisInfoFlowList]
            data.PlotTsum(  plot_info,fitDict=fitlistFull,tflowlist=tflowlist,
                            tcurr_list=tcurrlist,gammamomlist=gammamomlist)
        data = None
        print('Ratio correlator run complete')



class Alpha( ta.HasTraits ):
    """

    Controler for the 'Alpha' analysis which uses two-point correlators and flowed operators
    (combined).

    """

    Alpha_Noise_Subtraction = ta.Bool(False)
    # Ahmed_comparison = ta.Bool(False)
    Auto_correlation_analysis = ta.Bool(True)
    S_parameter = ta.Float(pa.defSparam)
    Only_Var_plot = ta.Bool(False)
    Mes_or_Bar = ta.Enum(ql.BMList)
    Bar_Mes_List = ta.List(ql.BarNumbListCPOdd)
    Bar_Mes_Number = ta.Enum(values = 'Bar_Mes_List')
    source_smearing = ta.Int(pa.defInfo['ism'])
    t_source = ta.Int(pa.defInfo['tsrc'])
    # combine_sinks = ta.Bool('sm2ptRat' in pa.defInfo.keys())
    # not_combine_sinks = ta.Bool(not ('sm2ptRat' in pa.defInfo.keys()))
    combine_sinks = ta.Bool(False)
    not_combine_sinks = ta.Bool(True)
    sink_smearing = ta.Int(pa.defInfo['jsm'])
    smlist = list(map(str,pa.defInfo['sm2ptRat']))
    smear_list = ta.Str(', '.join(['sm'+ival for ival in smlist]))
    smear_comment = ta.Str(' Set coeff to 0 to ignore smearing ')
    sink_coeffs = ta.Array(np.float64,(len(smlist),),pa.defInfo['jsmCoeffs'])
    sink_label = ta.Str(pa.defInfo['sinklab'])
    NNQ_fit = ta.Bool(False)
    alpha_fit_range = ta.Str(pa.defInfoFlow['FitsAlpha'][0])
    alpha_fit_comment = ta.Str('fit range is also what range is plotted')
    momenta_comment = ta.Str('First momenta in list is plotted')
    momentum_list = ta.List(ta.Str,pa.defpmomlist)
    momentum_note = ta.Str(' Only the first momentum will be plotted, order correctly ')
    # momentum_list = ta.List(ta.Str,['0 0 0'])
    # pformlist = ['p'+imom.replace(' ','') for imom in momentum_list]

    tflow_plot_min = ta.Float(pa.defxlimOpMore[0])
    tflow_plot_max = ta.Float(pa.defxlimOpMore[-1])
    tflow_plot_inc = ta.Float(pa.defxlimOpMore[1]-pa.defxlimOpMore[0])
    observable = ta.Enum(ql.ObsList)
    tflow_fit = ta.Str(pa.defInfoFlow['tflowfit'][0])
    tsink_fit = ta.Str(pa.defInfoFlow['tsinkfit'][0])
    RUN = ta.Button()
    CreateJob = ta.Button()
    CreateSeparateJob = ta.Button()
    Apply = ta.Button()
    Load = ta.Button()


    view = tua.View(
        '_',
        tua.Item('Alpha_Noise_Subtraction'),
        # tua.Item('Improved_Overlap',enabled_when='Do_Improved'),
        # tua.Item('Ahmed_comparison'),
        tua.Item('Auto_correlation_analysis'),
        tua.Item('S_parameter'),
        tua.Item('Only_Var_plot'),
        tua.Item('Mes_or_Bar'),
        tua.Item('Bar_Mes_Number'),
        tua.Item('source_smearing'),
        tua.Item('t_source'),
        tua.Item('combine_sinks'),
        tua.Item('sink_smearing',enabled_when='not_combine_sinks'),
        tua.Item('smear_list',enabled_when='combine_sinks',style='readonly'),
        tua.Item('smear_comment',show_label=False,enabled_when='combine_sinks',style='readonly'),
        tua.Item('sink_coeffs',enabled_when='combine_sinks',style='custom'),
        tua.Item('sink_label',enabled_when='combine_sinks'),
        tua.Item('NNQ_fit'),
        tua.Item('alpha_fit_range'),
        tua.Item('alpha_fit_comment',show_label=False,style='readonly'),
        tua.Item('momenta_comment',show_label=False,style='readonly'),
        tua.Item('momentum_list'),
        tua.Item('tflow_plot_min'),
        tua.Item('tflow_plot_max'),
        tua.Item('tflow_plot_inc'),
        tua.Item('observable'),
        tua.Item('tflow_fit'),
        tua.Item('tsink_fit'),
        tua.Item('Apply', show_label=False),
        tua.Item('CreateJob', show_label=False ),
        tua.Item('CreateSeparateJob', show_label=False ),
        tua.Item('RUN', show_label=False ),
        resizable=True
    )

    def _Mes_or_Bar_changed(self):
        if self.Mes_or_Bar == 'Baryon':
            self.Bar_Mes_List = ql.BarNumbListCPOdd
        elif self.Mes_or_Bar == 'Meson':
            self.Bar_Mes_List = ql.MesNumbList


    def _Load_fired(self):
        self.tflow_plot_min = thisInfoFlowList[info_index]['tflowlist'][0]
        self.tflow_plot_max = thisInfoFlowList[info_index]['tflowlist'][-1]
        if len(thisInfoFlowList[info_index]['tflowlist']) == 1:
            self.tflow_plot_inc = 1
        else:
            self.tflow_plot_inc = thisInfoFlowList[info_index]['tflowlist'][1] - self.tflow_plot_min

        self.pformlist = [str('p'+imom.replace(' ','')) for imom in thisInfoFlowList[info_index]['pmom']]


        self.Mes_or_Bar          = thisInfoFlowList[info_index]['MesOrBar']
        if self.Mes_or_Bar == 'Meson':
            self.Bar_Mes_Number = ChromaTOGamma(thisInfoFlowList[info_index]['Interp'])
        elif self.Mes_or_Bar == 'Baryon':
            self.Bar_Mes_Number = thisInfoFlowList[info_index]['Interp']
        self.source_smearing     = thisInfoFlowList[info_index]['ism']
        self.t_source            = thisInfoFlowList[info_index]['tsrc']
        self.sink_smearing       = thisInfoFlowList[info_index]['jsm']
        self.alpha_fit_range   = thisInfoFlowList[info_index]['FitsAlpha'][0]
        if isinstance(thisInfoFlowList[info_index]['pmom'],str):
            self.momentum_list       = [thisInfoFlowList[info_index]['pmom']]
        else:
            self.momentum_list       = list(thisInfoFlowList[info_index]['pmom'])
        if 'Observable' in list(thisInfoFlowList[info_index].keys()):
            if 'TopCharge' in thisInfoFlowList[info_index]['Observable']:
                self.observable=                   'TopCharge'
            elif 'Weinberg' in thisInfoFlowList[info_index]['Observable']:
                self.observable=                   'Weinberg'
        self.tflow_fit     = thisInfoFlowList[info_index]['tflowfit'][0]
        self.tsink_fit      = thisInfoFlowList[info_index]['tsinkfit'][0]
        self.combine_sinks         = thisInfoFlowList[info_index]['combine_sinks']
        self.Auto_correlation_analysis  = thisInfoFlowList[info_index]['autocorr']
        self.S_parameter    = thisInfoFlowList[info_index]['Sparam']

        self.combine_sinks = 'sm2ptRat' in list(thisInfoFlowList[info_index].keys())
        if self.combine_sinks:
            self.sink_label = thisInfoFlowList[info_index]['sinklab']

    def _combine_sinks_changed(self):
        self.not_combine_sinks = not self.combine_sinks

    def _Apply_fired(self):
        thisInfoFlowList[info_index]['MesOrBar'] = self.Mes_or_Bar
        if self.Mes_or_Bar == 'Meson':
            thisInfoFlowList[info_index]['Interp'] = GammaTOChroma(self.Bar_Mes_Number)
        else:
            thisInfoFlowList[info_index]['Interp'] = self.Bar_Mes_Number
        thisInfoFlowList[info_index]['ism'] = self.source_smearing
        thisInfoFlowList[info_index]['tsrc'] = self.t_source
        thisInfoFlowList[info_index]['jsm'] = self.sink_smearing
        thisInfoFlowList[info_index]['FitsAlpha'] = [self.alpha_fit_range]
        thisInfoFlowList[info_index]['pmom'] = Remove_Empty_Str(list(map(str,self.momentum_list)))
        self.pformlist = [str('p'+imom.replace(' ','')) for imom in thisInfoFlowList[info_index]['pmom']]

        # if self.Do_Improved:
        #     thisInfoFlowList[info_index]['IExpect'] = float(self.Improved_Overlap)
        thisInfoFlowList[info_index]['tflowlist'] = np.arange(max(self.tflow_plot_min,0.01),min(self.tflow_plot_max+self.tflow_plot_inc,10),self.tflow_plot_inc)
        thisInfoFlowList[info_index]['Observable'] = self.observable
        thisInfoFlowList[info_index]['tflowfit']= [str(self.tflow_fit)]
        thisInfoFlowList[info_index]['tsinkfit']= [str(self.tsink_fit)]
        ## TODO: this currently does not work, need to me moved into a function and called in _RUN_fired
        thisInfoFlowList[info_index]['combine_sinks'] = self.combine_sinks
        thisInfoFlowList[info_index]['autocorr']= self.Auto_correlation_analysis
        thisInfoFlowList[info_index]['Sparam'] = self.S_parameter
        if self.combine_sinks:
            thisInfoFlowList[info_index]['sinklab'] = self.sink_label
            thisInfoFlowList[info_index]['sinklen'] = len(self.smlist)
            smint,scoeffs = list(map(int,self.smlist)),list(map(np.float64,self.sink_coeffs))
            # if not all([ilist['jsm'] == jsm for ilist,jsm in zip(thisInfoFlowList[-3:],smint)]):
            for icount,(jsm,icoeff) in enumerate(zip(smint,scoeffs)):
                if len(thisInfoFlowList) < info_index+icount+1:
                    thisInfoFlowList.append(deepcopy(thisInfoFlowList[info_index]))
                    thisInfoList.append(deepcopy(thisInfoList[info_index]))
                thisInfoFlowList[info_index+icount]['combine_sinks'] = True
                thisInfoFlowList[info_index+icount]['jsm'] = jsm
                thisInfoFlowList[info_index+icount]['jsmCoeffs'] = [icoeff]
                # thisInfoFlowList[info_index]['sm2ptRat'] = map(int,self.smlist)
                # thisInfoFlowList[info_index]['jsmCoeffs'] = map(np.float64,self.sink_coeffs)
            # else:
            #     thisInfoFlowList[info_index]['combine_sinks'] = False
        else:
            thisInfoFlowList[info_index]['sinklen'] = 1
            if 'sm2ptRat' in list(thisInfoFlowList[info_index].keys()):
                del thisInfoFlowList[info_index]['sm2ptRat']

    def _CreateJob_fired(self):
        print()
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,JobNumber
        SavePreviousRun()
        this_job_file = JobFolder+'Alpha'+str(JobNumber)+'.py3p'
        IO.WritePickle( this_job_file,[thisInfoFlowList,bool(self.Alpha_Noise_Subtraction),
                        def_Corr_Sets,def_checkcfgs,thisDefWipe])
        WriteCSH('AAlpha',this_job_file,job_params,JobNumber)
        print('Created job script for Alpha run', this_job_file.replace('.py3p','.csh'))
        JobNumber += 1

    def _CreateSeparateJob_fired(self):
        print()
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,JobNumber
        SavePreviousRun()
        for iinfo in thisInfoFlowList:
            this_job_file = JobFolder+'Alpha'+str(JobNumber)+'.py3p'
            IO.WritePickle( this_job_file,[[iinfo],bool(self.Alpha_Noise_Subtraction),
                            def_Corr_Sets,def_checkcfgs,thisDefWipe])
            WriteCSH('AAlpha',this_job_file,job_params,JobNumber)
            print('Created job script for Alpha run', this_job_file.replace('.py3p','.csh'))
            JobNumber += 1


    def _RUN_fired(self):
        print()
        print('Begining CP Odd Two-point correlator run')
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,plot_info
        SavePreviousRun()
        data = ARun(thisInfoFlowList,bool(self.Alpha_Noise_Subtraction),def_Corr_Sets,def_checkcfgs,thisDefWipe)
        # data = so.SetOfNNQ(InfoDict=thisInfoFlowList,corr_cfglists=def_Corr_Sets,CompAhmed=bool(self.Ahmed_comparison))
        # data.LoadPickle(DefWipe=thisDefWipe,CheckCfgs=def_checkcfgs,CheckMom=True)
        # thisname,thisleglab,jsmlab = CreateNameAndLL(data.GetSet('First'),thisInfoFlowList[info_index]['sinklab'])
        fitlist = [idict['FitsAlpha'][0] for idict in thisInfoFlowList]
        tflowlist = [idict['tflowfit'][0] for idict in thisInfoFlowList]
        tsinklist = [idict['tsinkfit'][0] for idict in thisInfoFlowList]
        alltflow = [ilist['tflowlist'] for ilist in thisInfoFlowList]
        ## will plot first momenta
        momplotlist = [['p'+idict['pmom'][0].replace(' ','')] for idict in thisInfoFlowList]
        allmomlist = [['p'+jdict.replace(' ','') for jdict in idict['pmom']] for idict in thisInfoFlowList]
        for initlist,idict in enumerate(thisInfoFlowList):
            if idict['combine_sinks']:
                sinkcoeffs = []
                these_indicies = []
                for icount in range(idict['sinklen']):
                    these_indicies.append(icount+initlist)
                    sinkcoeffs.append(thisInfoFlowList[icount+initlist]['jsmCoeffs'][0])

                thisname,thisleglab,jsmlab = CreateNameAndLL(data.SetC2,thisInfoFlowList[initlist]['sinklab'],from_index=initlist)
                thisfun = GenSinkFun(list(map(np.float64,sinkcoeffs)))
                print('combining indicies ', ','.join(map(str,these_indicies)))
                data.AppendCCorrs(thisfun,filename=thisname,jsmcoeffs=list(map(np.float64,sinkcoeffs)),LegLab=thisleglab,jsm = thisInfoFlowList[initlist]['sinklab'],indicies = these_indicies)
                fitlist += [fitlist[0]]
                tflowlist += [tflowlist[0]]
                tsinklist += [tsinklist[0]]
                alltflow += [alltflow[0]]
                momplotlist += [momplotlist[0]]
                allmomlist += [allmomlist[0]]
                break
        if bool(self.Alpha_Noise_Subtraction):
            fitlist += [fitlist[0]]
            tflowlist += [tflowlist[0]]
            tsinklist += [tsinklist[0]]
            alltflow += [alltflow[0]]
            momplotlist += [momplotlist[0]]
            allmomlist += [allmomlist[0]]


        # if self.combine_sinks:
        #     thisname,thisleglab,jsmlab = CreateNameAndLL(data.SetC2,thisInfoFlowList[info_index]['sinklab'])
        #     thisfun = GenSinkFun(map(np.float64,self.sink_coeffs))
        #     data.AppendCCorrs(thisfun,filename=thisname,LegLab=thisleglab,jsm = thisInfoFlowList[info_index]['sinklab'])
        #     fitlist += [fitlist[0]]
        #     tflowlist += [tflowlist[0]]
        #     tsinklist += [tsinklist[0]]
        #     alltflow += [alltflow[0]]
        #     momplotlist += [momplotlist[0]]
        data.AlphaPlot(plot_info,tflowlist=tflowlist,fitDict=fitlist,momlist=momplotlist,NNQfit=self.NNQ_fit,OnlyVar=self.Only_Var_plot)
        data.AlphaPlotTflow(plot_info,tlist=tsinklist,tflowlist=alltflow,momlist=momplotlist,OnlyVar=self.Only_Var_plot)
        data.AlphaPlotTauIntWopt(plot_info,momlist=momplotlist,OnlyVar=self.Only_Var_plot)
        data.AlphaPlotTflowTauIntWopt(plot_info,tlist=tsinklist,tflowlist=alltflow,momlist=momplotlist,OnlyVar=self.Only_Var_plot)
        data.AlphaPlotVsEnergy(plot_info,fitDict=fitlist,momlist=allmomlist,OnlyVar=self.Only_Var_plot)
        data.AlphaPlotVsLatSpace(plot_info,fitlist,momlist=momplotlist,plottflowlist=tflowlist,OnlyVar=self.Only_Var_plot)
        data.AlphaPlotVsMpi(plot_info,fitlist,momlist=momplotlist,plottflowlist=tflowlist,OnlyVar=self.Only_Var_plot)
        data = None
        print('CP Odd Two-point correlator run complete')




class AlphaFull( ta.HasTraits ):
    """

    Controler for the 'Alpha' analysis where the Euclidean time dependance of the
    flowed operators is symmetically summed around some location relative to the
    two-point correlator.

    """
    alpha_save_mem = ta.Bool(False)
    alpha_improve_type = ta.Enum(ql.alpha_improv_types)
    alpha_imp_NNQ2Sub = ta.Bool(False)
    Plot_Vs_Energy = ta.Bool(False)
    alpha_sum_type = ta.Enum(ql.alpha_sum_list)
    Auto_correlation_analysis = ta.Bool(True)
    S_parameter = ta.Float(pa.defSparam)
    Only_Var_plot = ta.Bool(False)
    Mes_or_Bar = ta.Enum(ql.BMList)
    Bar_Mes_List = ta.List(ql.BarNumbListCPOdd)
    Bar_Mes_Number = ta.Enum(values = 'Bar_Mes_List')
    source_smearing = ta.Int(pa.defInfo['ism'])
    t_source = ta.Int(pa.defInfo['tsrc'])
    # combine_sinks = ta.Bool('sm2ptRat' in pa.defInfo.keys())
    # not_combine_sinks = ta.Bool(not ('sm2ptRat' in pa.defInfo.keys()))
    combine_sinks = ta.Bool(False)
    not_combine_sinks = ta.Bool(True)
    sink_smearing = ta.Int(pa.defInfo['jsm'])
    smlist = list(map(str,pa.defInfo['sm2ptRat']))
    smear_list = ta.Str(', '.join(['sm'+ival for ival in smlist]))
    smear_comment = ta.Str(' Set coeff to 0 to ignore smearing ')
    sink_coeffs = ta.Array(np.float64,(len(smlist),),pa.defInfo['jsmCoeffs'])
    sink_label = ta.Str(pa.defInfo['sinklab'])
    # alpha_fit_range = ta.Str(pa.defInfoFlow['FitsAlphaFull'][0][0])
    alpha_fit_range = ta.Str('fitr5-15')
    alpha_tsum_fit_range = ta.Str(pa.defInfoFlow['FitsAlphaFull'][0][1])
    momenta_comment = ta.Str('First momenta in list is plotted')
    # momentum_list = ta.List(ta.Str,pa.defpmomlist)
    momentum_list = ta.List(ta.Str,['0 0 0'])
    momentum_note = ta.Str(' Only the first momentum will be plotted, order correctly ')
    fit_par_picked = ta.Str('First')
    # momentum_list = ta.List(ta.Str,['0 0 0'])
    # pformlist = ['p'+imom.replace(' ','') for imom in momentum_list]

    # tflow_plot_min = ta.Float(pa.defxlimOp[0])
    # tflow_plot_max = ta.Float(pa.defxlimOp[-1])
    # tflow_plot_inc = ta.Float(pa.defxlimOp[1]-pa.defxlimOp[0])
    tflow_plot_min = ta.Float(0.01)
    tflow_plot_max = ta.Float(6.01)
    tflow_plot_inc = ta.Float(1)
    tsum_boot_min = ta.Int(0)
    tsum_boot_max = ta.Int(pa.defInfo['nxyzt'][-1]-1)
    tsum_spec_boot_min = ta.Int(0)
    tsum_spec_boot_max = ta.Int(pa.defInfo['nxyzt'][-1]-1)
    # tsum_spec_boot_min = ta.Int(20)
    # tsum_spec_boot_max = ta.Int(pa.defInfo['nxyzt'][-1]-21)
    # t_boot_min = ta.Int(0)
    # t_boot_max = ta.Int(pa.defInfo['nxyzt'][-1]-1)
    # t_boot_inc = ta.Int(1)
    t_boot_min = ta.Int(4)
    t_boot_max = ta.Int(14)
    t_boot_inc = ta.Int(5)
    observable = ta.Enum(ql.ObsList)
    tflow_picked = ta.Str(pa.defInfoFlow['tflowfit'][0])
    tsink_picked = ta.Str(pa.defInfoFlow['tsinkfit'][0])
    # tsum_picked = ta.Str(pa.defInfoFlow['tsum_fit_list'][0])
    tsum_picked = ta.Str('ts7')
    tsum_spec_picked = ta.Str('tss32')
    do_tsink_fit = ta.Bool(False)
    do_multiple_fits = ta.Bool(False)
    t_fit_min_min = ta.Int(2)
    t_fit_min_max = ta.Int(8)
    t_fit_max_min = ta.Int(15)
    t_fit_max_max = ta.Int(20)
    fit_min_picked = ta.Int(2)
    fit_max_picked = ta.Int(20)
    RUN = ta.Button()
    CreateJob = ta.Button()
    CreateSeparateJob = ta.Button()
    Apply = ta.Button()
    Load = ta.Button()


    view = tua.View(tua.Group(
        tua.Item('alpha_save_mem'),
        tua.Item('Plot_Vs_Energy'),
        tua.Item('alpha_improve_type'),
        tua.Item('alpha_sum_type'),
        tua.Item('Auto_correlation_analysis'),
        tua.Item('S_parameter'),
        tua.Item('Only_Var_plot'),
        tua.Item('Mes_or_Bar'),
        tua.Item('Bar_Mes_Number'),
        tua.Item('source_smearing'),
        tua.Item('t_source'),
        tua.Item('combine_sinks'),
        tua.Item('sink_smearing',enabled_when='not_combine_sinks'),
        tua.Item('smear_list',enabled_when='combine_sinks',style='readonly'),
        tua.Item('smear_comment',show_label=False,enabled_when='combine_sinks',style='readonly'),
        tua.Item('sink_coeffs',enabled_when='combine_sinks',style='custom'),
        tua.Item('sink_label',enabled_when='combine_sinks'),
        tua.Item('alpha_fit_range'),
        tua.Item('alpha_tsum_fit_range'),
        tua.Item('Apply', show_label=False),
        tua.Item('CreateJob', show_label=False ),
        tua.Item('CreateSeparateJob', show_label=False ),
        tua.Item('RUN', show_label=False ),
        label='Page 1'
        ),tua.Group(
        tua.Item('momenta_comment',show_label=False,style='readonly'),
        tua.Item('momentum_list'),
        tua.Item('tflow_plot_min'),
        tua.Item('tflow_plot_max'),
        tua.Item('tflow_plot_inc'),
        tua.Item('tsum_boot_min'),
        tua.Item('tsum_boot_max'),
        tua.Item('tsum_spec_boot_min',enabled_when='alpha_imp_NNQ2Sub'),
        tua.Item('tsum_spec_boot_max',enabled_when='alpha_imp_NNQ2Sub'),
        tua.Item('t_boot_min'),
        tua.Item('t_boot_max'),
        tua.Item('t_boot_inc'),
        tua.Item('observable'),
        tua.Item('tflow_picked'),
        tua.Item('tsink_picked'),
        tua.Item('tsum_picked'),
        tua.Item('tsum_spec_picked'),
        tua.Item('fit_par_picked'),
        tua.Item('do_tsink_fit'),
        tua.Item('Apply', show_label=False),
        tua.Item('CreateJob', show_label=False ),
        tua.Item('CreateSeparateJob', show_label=False ),
        tua.Item('RUN', show_label=False ),
        label='Page 2'
        ),
        tua.Group(
        tua.Item('do_multiple_fits'),
        tua.Item('t_fit_min_min',enabled_when='do_multiple_fits'),
        tua.Item('t_fit_min_max',enabled_when='do_multiple_fits'),
        tua.Item('t_fit_max_min',enabled_when='do_multiple_fits'),
        tua.Item('t_fit_max_max',enabled_when='do_multiple_fits'),
        tua.Item('fit_min_picked',enabled_when='do_multiple_fits'),
        tua.Item('fit_max_picked',enabled_when='do_multiple_fits'),
        tua.Item('Apply', show_label=False),
        tua.Item('CreateJob', show_label=False ),
        tua.Item('CreateSeparateJob', show_label=False ),
        tua.Item('RUN', show_label=False ),
        label='Mult_Fits'
        ),
    resizable=True)

    def _Mes_or_Bar_changed(self):
        if self.Mes_or_Bar == 'Baryon':
            self.Bar_Mes_List = ql.BarNumbListCPOdd
        elif self.Mes_or_Bar == 'Meson':
            self.Bar_Mes_List = ql.MesNumbList


    def _Load_fired(self):
        self.tflow_plot_min = thisInfoFlowList[info_index]['tflowlist'][0]
        self.tflow_plot_max = thisInfoFlowList[info_index]['tflowlist'][-1]
        if len(thisInfoFlowList[info_index]['tflowlist']) == 1:
            self.tflow_plot_inc = 1
        else:
            self.tflow_plot_inc = thisInfoFlowList[info_index]['tflowlist'][1] - self.tflow_plot_min

        self.t_boot_min = thisInfoFlowList[info_index]['boot_nt_list'][0]
        self.t_boot_max = thisInfoFlowList[info_index]['boot_nt_list'][-1]
        if len(thisInfoFlowList[info_index]['boot_nt_list']) == 1:
            self.t_boot_inc = 1
        else:
            self.t_boot_inc = thisInfoFlowList[info_index]['boot_nt_list'][1] - self.t_boot_min

        if isinstance(thisInfoFlowList[info_index]['boot_tsum_list'][0],str):
            self.tsum_boot_min = int(thisInfoFlowList[info_index]['boot_tsum_list'][0].replace('ts',''))
            self.tsum_boot_max = int(thisInfoFlowList[info_index]['boot_tsum_list'][-1].replace('ts',''))
        else:
            self.tsum_boot_min = thisInfoFlowList[info_index]['boot_tsum_list'][0]
            self.tsum_boot_max = thisInfoFlowList[info_index]['boot_tsum_list'][-1]

        self.alpha_fit_range = thisInfoFlowList[info_index]['FitsAlphaFull'][0][0]
        self.alpha_tsum_fit_range = thisInfoFlowList[info_index]['FitsAlphaFull'][0][1]
        if thisInfoFlowList[info_index]['op_trange_sum_type'] != 'None':
            self.alpha_sum_type = '_'.join(thisInfoFlowList[info_index]['op_trange_sum_type'].split('_')[:-1])
            self.alpha_improve_type = thisInfoFlowList[info_index]['op_trange_sum_type'].split('_')[-1]

        self.Auto_correlation_analysis               = thisInfoFlowList[info_index]['autocorr']
        self.S_parameter                             = thisInfoFlowList[info_index]['Sparam']
        self.fit_par_picked                     = thisInfoFlowList[info_index]['fit_par_picked']
        self.Mes_or_Bar                              = thisInfoFlowList[info_index]['MesOrBar']
        if self.Mes_or_Bar == 'Meson':
            self.Bar_Mes_Number = ChromaTOGamma(thisInfoFlowList[info_index]['Interp'])
        elif self.Mes_or_Bar == 'Baryon':
            self.Bar_Mes_Number = thisInfoFlowList[info_index]['Interp']
        self.source_smearing                         = thisInfoFlowList[info_index]['ism']
        self.t_source                                = thisInfoFlowList[info_index]['tsrc']
        self.sink_smearing                           = thisInfoFlowList[info_index]['jsm']
        # self.alpha_fit_range                    = thisInfoFlowList[info_index]['alpha_tfit_range']
        if isinstance(thisInfoFlowList[info_index]['pmom'],str):
            self.momentum_list = [thisInfoFlowList[info_index]['pmom']]
        else:
            self.momentum_list = list(thisInfoFlowList[info_index]['pmom'])
        if 'Observable' in list(thisInfoFlowList[info_index].keys()):
            if 'TopCharge' in thisInfoFlowList[info_index]['Observable']:
                self.observable=                   'TopCharge'
            elif 'Weinberg' in thisInfoFlowList[info_index]['Observable']:
                self.observable=                   'Weinberg'
        self.tflow_picked                     = thisInfoFlowList[info_index]['tflowfit'][0]
        self.tsink_picked                    = thisInfoFlowList[info_index]['tsinkfit'][0]
        if isinstance(thisInfoFlowList[info_index]['tsum_fit_list'][0],int):
            self.tsum_picked                      = 'ts'+str(thisInfoFlowList[info_index]['tsum_fit_list'][0])
        else:
            self.tsum_picked                      = thisInfoFlowList[info_index]['tsum_fit_list'][0]
        self.tsum_spec_picked                 = thisInfoFlowList[info_index]['tsum_spec_fit_list'][0]
        self.combine_sinks                           = thisInfoFlowList[info_index]['combine_sinks']

        self.pformlist = [str('p'+imom.replace(' ','')) for imom in thisInfoFlowList[info_index]['pmom']]

        if 'fit_min_picked' in list(thisInfoFlowList[info_index].keys()):
            self.do_multiple_fits = thisInfoFlowList[info_index]['fit_min_picked'] == float('NaN')
        if self.do_multiple_fits:
            split_list = [list(map(int,ix[1].replace('tsumfitr','').split('-'))) for ix in thisInfoFlowList[info_index]['fit_range']]
            start_list,fin_list = list(zip(*split_list))
            self.t_fit_min_min = np.min(start_list)
            self.t_fit_min_max = np.max(start_list)
            self.t_fit_max_min = np.min(fin_list)
            self.t_fit_max_max = np.max(fin_list)
            self.fit_min_picked =    thisInfoFlowList[info_index]['fit_min_picked']
            self.fit_max_picked =    thisInfoFlowList[info_index]['fit_max_picked']

        ## TODO: this currently does not work, need to me moved into a function and called in _RUN_fired
        if self.combine_sinks:
            self.sink_label =  thisInfoFlowList[info_index]['sinklab']


    def _combine_sinks_changed(self):
        self.not_combine_sinks = not self.combine_sinks

    def _alpha_improve_type_fired(self):
        self.alpha_imp_NNQ2Sub = 'NNQ2sub' == str(self.alpha_improve_type)

    def _Apply_fired(self):
        thisInfoFlowList[info_index]['autocorr']= self.Auto_correlation_analysis
        thisInfoFlowList[info_index]['Sparam'] = self.S_parameter
        thisInfoFlowList[info_index]['fit_par_picked'] = str(self.fit_par_picked)
        thisInfoFlowList[info_index]['op_trange_sum_type'] = str(self.alpha_sum_type)
        thisInfoFlowList[info_index]['op_trange_sum_type'] += '_'+str(self.alpha_improve_type)
        thisInfoFlowList[info_index]['MesOrBar'] = self.Mes_or_Bar
        if self.Mes_or_Bar == 'Meson':
            thisInfoFlowList[info_index]['Interp'] = GammaTOChroma(self.Bar_Mes_Number)
        else:
            thisInfoFlowList[info_index]['Interp'] = self.Bar_Mes_Number
        thisInfoFlowList[info_index]['ism'] = self.source_smearing
        thisInfoFlowList[info_index]['tsrc'] = self.t_source
        thisInfoFlowList[info_index]['jsm'] = self.sink_smearing
        thisInfoFlowList[info_index]['alpha_tfit_range'] = str(self.alpha_fit_range)
        thisInfoFlowList[info_index]['FitsAlphaFull'] = [[str(self.alpha_fit_range),str(self.alpha_tsum_fit_range)]]
        thisInfoFlowList[info_index]['pmom'] = Remove_Empty_Str(list(map(str,self.momentum_list)))
        if thisInfoFlowList[info_index]['pmom'][0] == 'All':
            thisInfoFlowList[info_index]['pmom'] = ['0 0 0','1 0 0','1 1 0','1 1 1','2 0 0']
        self.pformlist = [str('p'+imom.replace(' ','')) for imom in thisInfoFlowList[info_index]['pmom']]

        thisInfoFlowList[info_index]['tflowlist'] = np.arange(max(self.tflow_plot_min,0.01),min(self.tflow_plot_max+self.tflow_plot_inc,10),self.tflow_plot_inc)
        thisInfoFlowList[info_index]['boot_tsum_list'] = list(range(self.tsum_boot_min,self.tsum_boot_max+1))
        thisInfoFlowList[info_index]['boot_spec_tsum_list'] = list(range(self.tsum_spec_boot_min,self.tsum_spec_boot_max+1))

        thisInfoFlowList[info_index]['boot_nt_list'] = np.arange(max(self.t_boot_min,0),min(self.t_boot_max+1,pa.defInfo['nxyzt'][-1]),self.t_boot_inc)
        thisInfoFlowList[info_index]['Observable'] = self.observable
        thisInfoFlowList[info_index]['tflowfit']= [str(self.tflow_picked)]
        thisInfoFlowList[info_index]['tsinkfit']= [str(self.tsink_picked)]
        thisInfoFlowList[info_index]['tsum_fit_list'] = [str(self.tsum_picked)]
        thisInfoFlowList[info_index]['tsum_spec_fit_list'] = [str(self.tsum_spec_picked)]
        thisInfoFlowList[info_index]['fit_list']= thisInfoFlowList[info_index]['FitsAlphaFull']
        if self.do_multiple_fits:
            thisInfoFlowList[info_index]['fit_range'] = []
            for fitmin in range(self.t_fit_min_min,self.t_fit_min_max+1):
                for fitmax in range(self.t_fit_max_min,self.t_fit_max_max+1):
                    thisInfoFlowList[info_index]['fit_range'].append([str(self.alpha_fit_range),'tsumfitr'+str(fitmin)+'-'+str(fitmax)])

            thisInfoFlowList[info_index]['fit_min_picked'] = str(self.fit_min_picked)
            thisInfoFlowList[info_index]['fit_max_picked'] = str(self.fit_max_picked)
        else:
            thisInfoFlowList[info_index]['fit_range']= [str(self.alpha_fit_range),str(self.alpha_tsum_fit_range)]
            thisInfoFlowList[info_index]['fit_min_picked'] = str('NaN')
            thisInfoFlowList[info_index]['fit_max_picked'] = str('NaN')
        ## TODO: this currently does not work, need to me moved into a function and called in _RUN_fired
        thisInfoFlowList[info_index]['combine_sinks'] = self.combine_sinks
        if self.combine_sinks:
            thisInfoFlowList[info_index]['sinklab'] = self.sink_label
            thisInfoFlowList[info_index]['sinklen'] = len(self.smlist)
            smint,scoeffs = list(map(int,self.smlist)),list(map(np.float64,self.sink_coeffs))
            # if not all([ilist['jsm'] == jsm for ilist,jsm in zip(thisInfoFlowList[-3:],smint)]):
            for icount,(jsm,icoeff) in enumerate(zip(smint,scoeffs)):
                if len(thisInfoFlowList) < info_index+icount+1:
                    thisInfoFlowList.append(deepcopy(thisInfoFlowList[info_index]))
                    thisInfoList.append(deepcopy(thisInfoList[info_index]))
                thisInfoFlowList[info_index+icount]['combine_sinks'] = True
                thisInfoFlowList[info_index+icount]['jsm'] = jsm
                thisInfoFlowList[info_index+icount]['jsmCoeffs'] = [icoeff]
        else:
            thisInfoFlowList[info_index]['sinklen'] = 1
            if 'sm2ptRat' in list(thisInfoFlowList[info_index].keys()):
                del thisInfoFlowList[info_index]['sm2ptRat']

    def _CreateJob_fired(self):
        print()
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,JobNumber
        SavePreviousRun()
        this_job_file = JobFolder+'AlphaFull'+str(JobNumber)+'.py3p'
        IO.WritePickle(this_job_file,[thisInfoFlowList,def_Corr_Sets,def_checkcfgs,thisDefWipe,bool(self.alpha_save_mem)])
        WriteCSH('AAlphaFull',this_job_file,job_params,JobNumber)
        print('Created job script for Alpha Full Analysis run', this_job_file.replace('.py3p','.csh'))
        JobNumber += 1


    def _CreateSeparateJob_fired(self):
        print()
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,JobNumber
        SavePreviousRun()
        for iInfo in thisInfoFlowList:
            this_job_file = JobFolder+'AlphaFull'+str(JobNumber)+'.py3p'
            IO.WritePickle(this_job_file,[[iInfo],def_Corr_Sets,def_checkcfgs,thisDefWipe,bool(self.alpha_save_mem)])
            WriteCSH('AAlphaFull',this_job_file,job_params,JobNumber)
            print('Created job script for Alpha Full Analysis run', this_job_file.replace('.py3p','.csh'))
            JobNumber += 1


    def _RUN_fired(self):
        print()
        print('Begining CP Odd Two-point Full correlator run')
        global thisDefWipe,def_checkcfgs,def_Corr_Sets,plot_info
        SavePreviousRun()
        data = AFullRun(thisInfoFlowList,def_Corr_Sets,def_checkcfgs,thisDefWipe,bool(self.alpha_save_mem))
        # data = so.SetOfNNQFull(InfoDict=thisInfoFlowList,corr_cfglists=def_Corr_Sets)
        # data.LoadPickle(DefWipe=thisDefWipe,CheckCfgs=def_checkcfgs,CheckMom=True)
        # thisname,thisleglab,jsmlab = CreateNameAndLL(data.GetSet('First'),thisInfoFlowList[info_index]['sinklab'])
        fitlist = [idict['FitsAlphaFull'][0] for idict in thisInfoFlowList]
        tflowlist = [idict['tflowfit'][0] for idict in thisInfoFlowList]
        tsinklist = [idict['tsinkfit'][0] for idict in thisInfoFlowList]
        tsumlist = [idict['tsum_fit_list'][0] for idict in thisInfoFlowList]
        tsum_spec_list = [idict['tsum_spec_fit_list'][0] for idict in thisInfoFlowList]
        fit_list = [idict['fit_range'] for idict in thisInfoFlowList]
        fit_list_comb = ['_'.join(idict['fit_range']) for idict in thisInfoFlowList]
        alltflow = [ilist['tflowlist'] for ilist in thisInfoFlowList]
        fit_par_picked = [iinfo['fit_par_picked'] for iinfo in thisInfoFlowList]
        ## will plot first momenta
        momplotlist = [['p'+idict['pmom'][0].replace(' ','')] for idict in thisInfoFlowList]
        allmomlist = [['p'+jdict.replace(' ','') for jdict in idict['pmom']] for idict in thisInfoFlowList]
        fit_min_picked = [iinfo['fit_min_picked'] for iinfo in thisInfoFlowList]
        fit_max_picked = [iinfo['fit_max_picked'] for iinfo in thisInfoFlowList]
        alpha_tfit_range = [iinfo['alpha_tfit_range'] for iinfo in thisInfoFlowList]

        for initlist,idict in enumerate(thisInfoFlowList):
            if idict['combine_sinks']:
                sinkcoeffs = []
                these_indicies = []
                for icount in range(idict['sinklen']):
                    these_indicies.append(icount+initlist)
                    sinkcoeffs.append(thisInfoFlowList[icount+initlist]['jsmCoeffs'][0])

                thisname,thisleglab,jsmlab = CreateNameAndLL(data.SetC2,thisInfoFlowList[initlist]['sinklab'],from_index=initlist)
                thisfun = GenSinkFun(list(map(np.float64,sinkcoeffs)))
                print('combining indicies ', ','.join(map(str,these_indicies)))
                data.AppendCCorrs(thisfun,filename=thisname,jsmcoeffs=list(map(np.float64,sinkcoeffs)),LegLab=thisleglab,jsm = thisInfoFlowList[initlist]['sinklab'],indicies = these_indicies)
                fitlist += [fitlist[0]]
                tflowlist += [tflowlist[0]]
                tsinklist += [tsinklist[0]]
                tsumlist += [tsumlist[0]]
                alltflow += [alltflow[0]]
                momplotlist += [momplotlist[0]]
                allmomlist += [allmomlist[0]]
                fit_list += [fit_list[0]]
                fit_list_comb += [fit_list_comb[0]]
                fit_par_picked += [fit_par_picked[0]]
                fit_min_picked += [fit_min_picked[0]]
                fit_max_picked += [fit_max_picked[0]]
                alpha_tfit_range += [alpha_tfit_range[0]]
                tsum_spec_list += [tsum_spec_list[0]]
                break
        # if self.combine_sinks:
        #     thisname,thisleglab,jsmlab = CreateNameAndLL(data.SetC2,thisInfoFlowList[info_index]['sinklab'])
        #     thisfun = GenSinkFun(map(np.float64,self.sink_coeffs))
        #     data.AppendCCorrs(thisfun,filename=thisname,LegLab=thisleglab,jsm = thisInfoFlowList[info_index]['sinklab'])
        #     fitlist += [fitlist[0]]
        #     tflowlist += [tflowlist[0]]
        #     tsinklist += [tsinklist[0]]
        #     alltflow += [alltflow[0]]
        #     momplotlist += [momplotlist[0]]
        if not self.do_tsink_fit:
            print('Not doing fits for each tsink')
        data.AlphaPlot( plot_info,fitrlist=fit_list_comb,momlist=momplotlist,tsum_list=tsumlist,
                            spec_tsum_list=tsum_spec_list,tflowlist=tflowlist)
        data.AlphaPlotTsum( plot_info,momlist=momplotlist,tlist=tsinklist,
                            spec_tsum_list=tsum_spec_list,tflowlist=tflowlist,
                            fit_list=fit_list)
        data.AlphaPlotCovar( plot_info,momlist=momplotlist,tlist=tsinklist,
                            spec_tsum_list=tsum_spec_list,tflowlist=tflowlist)
        data.AlphaPlotTflow(plot_info,tlist=tsinklist,tsum_list=tsumlist,
                            tflowlist=alltflow,momlist=momplotlist)
        data.AlphaPlotDist( plot_info,tsinklist,tflowlist,tsumlist,tsum_spec_list=tsum_spec_list,
                            momlist=momplotlist)
        data.AlphaPlotVsMpi( plot_info,fit_list_comb,plottflowlist = tflowlist,momlist=momplotlist)
        data.AlphaPlotVsLatSpace( plot_info,fit_list_comb,plottflowlist = tflowlist,momlist=momplotlist)
        if self.Plot_Vs_Energy:
            data.AlphaPlotVsEnergy( plot_info,fitDict=fitlist,
                                    momlist=allmomlist,tflowlist=tflowlist)
        if bool(self.do_multiple_fits):
            data.AlphaPlotVstiFitr( tflowlist,alpha_tfit_range,fit_max_picked,
                                    plot_info,momlist=momplotlist,
                                    OnlyVar=False)
            data.AlphaPlotVstfFitr( tflowlist,alpha_tfit_range,fit_min_picked,
                                    plot_info,momlist=momplotlist,
                                    OnlyVar=False)
        data = None
        print('CP Odd Two-point correlator run complete')


class FormFacts( ta.HasTraits ):
    """

    Controler for the full Form Factor analysis.

    """


    CommonRatios = ta.Bool(False)
    Fit_Schiff = ta.Bool(False)
    not_Fit_Schiff = ta.Bool(True)
    Vanish_Chiral_Limit = ta.Bool(True)
    mixed_term  = ta.Bool(False)
    twod_fit = ta.Bool(True)
    Order_2_Fit = ta.Bool(True)
    Rat_alpha = ta.Bool(True)
    recalculate_all = ta.Bool(False)
    config_list_from_file = ta.Bool(True)
    all_cfgs_for_alpha = ta.Bool(True)
    source_smearing = ta.Int(pa.defInfo['ism'])
    t_source = ta.Int(pa.defInfo['tsrc'])
    combine_sinks = ta.Bool('sm2ptRat' in list(pa.defInfo.keys()))
    not_combine_sinks = ta.Bool(not ('sm2ptRat' in list(pa.defInfo.keys())))
    sink_smearing = ta.Int(pa.defInfo['jsm'])
    smlist = list(map(str,pa.defInfo['sm2ptRat']))
    smear_list = ta.Str(', '.join(['sm'+ival for ival in smlist]))
    smear_comment = ta.Str(' Set coeff to 0 to ignore smearing ')
    sink_coeffs = ta.Array(np.float64,(len(smlist),),pa.defInfo['jsmCoeffs'])
    sink_label = ta.Str(pa.defInfoFlow['sinklab'])
    t_sink = ta.Int(pa.defInfo['tsink'])
    sink_momentum = ta.Str(pa.defInfo['ppmom'])
    doublet_singlet_comb = ta.Enum(ql.DSList)
    current_type = ta.Enum(ql.CurrTypes)
    Rat_comb_alpha = ta.Bool(False)
    Rat_comb_alpha_fun = ta.Enum(op_str_list)
    form_factor = ta.Enum(ql.FFTypes)
    q2_momentum_list = ta.List(ta.Int,[0,1,2,3,4])
    pp2_momentum = ta.Int(0)
    EffM_state_fit = ta.Str(pa.defInfo['Fits'][0][0])
    EffM_fit_range = ta.Str(pa.defInfo['Fits'][0][1])
    alpha_fit_range = ta.Str(pa.defInfoFlow['FitsAlpha'][0])
    # rf_fit_range = ta.Str(pa.defInfoFlow['FitsRF'][0])
    Chi_Cutoff = ta.Bool(True)
    Only_Important_Ratios = ta.Bool(True)
    Reduce_System = ta.Bool(True)
    Chi2pdf_Min = ta.Float(0.3)
    Chi2pdf_Max = ta.Float(0.8)
    Chi2pdf_number_picked = ta.Str('Chi0-Chi100')
    RF_fit_range_minimum = ta.Int(3)
    RF_fit_min = ta.Int(2)
    RF_fit_max = ta.Int(pa.defInfo['tsink']-2)
    Colapse_Chi_Dist = ta.Bool(True)

    njnq_observable = ta.Bool(False)
    tflow_plot_min = ta.Float(pa.defxlimOp[0])
    tflow_plot_max = ta.Float(pa.defxlimOp[-1])
    tflow_plot_inc = ta.Float(pa.defxlimOp[1]-pa.defxlimOp[0])
    observable = ta.Enum(ql.ObsList)
    tflow_fit = ta.Str(pa.defInfoFlow['tflowfit'][0])
    full_ratios = ta.Bool(False)
    frat_bool = ta.Bool(False)
    Rat_tsum_range = ta.Str(pa.defInfoFlow['FitsAlphaFull'][0][1])
    Rat_sum_type = ta.Enum(ql.rat_sum_list)
    fal_bool = ta.Bool(False)
    full_alpha = ta.Bool(False)
    alpha_tsum_range = ta.Str(pa.defInfoFlow['FitsAlphaFull'][0][1])
    alpha_sum_type = ta.Enum(ql.alpha_sum_list)
    FF_combine = ta.Bool(False)
    FF_combine_type = ta.Enum(ql.CombineList_FF)
    FF_index_list = ta.List(int,[0,1])
    FF_coeff_list = ta.List(float,[1.0,1.0])
    RUN = ta.Button()
    Plot = ta.Button()
    PlotOnlyDist = ta.Button()
    CreateJob = ta.Button()
    CreateSeparateJob = ta.Button()
    Apply = ta.Button()
    Load = ta.Button()


    view = tua.View(
        tua.Group(
        tua.Item('CommonRatios'),
        tua.Item('recalculate_all'),
        tua.Item('Rat_alpha'),
        tua.Item('config_list_from_file'),
        tua.Item('all_cfgs_for_alpha'),
        tua.Item('source_smearing'),
        tua.Item('t_source'),
        tua.Item('combine_sinks'),
        tua.Item('sink_smearing',enabled_when='not_combine_sinks'),
        tua.Item('smear_list',enabled_when='combine_sinks',style='readonly'),
        tua.Item('smear_comment',show_label=False,enabled_when='combine_sinks',style='readonly'),
        tua.Item('sink_coeffs',enabled_when='combine_sinks',style='custom'),
        tua.Item('sink_label',enabled_when='combine_sinks'),
        tua.Item('doublet_singlet_comb'),
        tua.Item('current_type'),
        tua.Item('Rat_comb_alpha_fun',enabled_when='Rat_comb_alpha'),
        tua.Item('form_factor'),
        tua.Item('q2_momentum_list'),
        tua.Item('pp2_momentum'),
        tua.Item('Apply', show_label=False),
        tua.Item('CreateJob', show_label=False ),
        tua.Item('CreateSeparateJob', show_label=False ),
        tua.Item('RUN', show_label=False ),
        tua.Item('Plot', show_label=False ),
        tua.Item('PlotOnlyDist', show_label=False ),
        label='Page 1'
        ),tua.Group(
        tua.Item('EffM_state_fit'),
        tua.Item('EffM_fit_range'),
        tua.Item('alpha_fit_range'),
        tua.Item('Reduce_System'),
        tua.Item('Chi_Cutoff'),
        tua.Item('Only_Important_Ratios',enabled_when='Chi_Cutoff'),
        tua.Item('Chi2pdf_Min',enabled_when='Chi_Cutoff'),
        tua.Item('Chi2pdf_Max',enabled_when='Chi_Cutoff'),
        tua.Item('Chi2pdf_number_picked',enabled_when='Chi_Cutoff'),
        tua.Item('Colapse_Chi_Dist',enabled_when='Chi_Cutoff'),
        tua.Item('RF_fit_range_minimum'),
        tua.Item('RF_fit_min'),
        tua.Item('RF_fit_max'),
        tua.Item('t_sink'),
        tua.Item('Apply', show_label=False),
        tua.Item('CreateJob', show_label=False ),
        tua.Item('CreateSeparateJob', show_label=False ),
        tua.Item('RUN', show_label=False ),
        tua.Item('Plot', show_label=False ),
        tua.Item('PlotOnlyDist', show_label=False ),
        # tua.Item('njnq_observable',style='readonly'),
        label='Page 2'
        ),tua.Group(
        tua.Item('tflow_plot_min',enabled_when='njnq_observable'),
        tua.Item('tflow_plot_max',enabled_when='njnq_observable'),
        tua.Item('tflow_plot_inc',enabled_when='njnq_observable'),
        tua.Item('tflow_fit',enabled_when='njnq_observable'),
        tua.Item('observable',enabled_when='njnq_observable',style='readonly'),
        tua.Item('full_ratios',enabled_when='njnq_observable'),
        tua.Item('Rat_sum_type',enabled_when='frat_bool'),
        tua.Item('Rat_tsum_range',enabled_when='frat_bool'),
        tua.Item('full_alpha',enabled_when='njnq_observable'),
        tua.Item('alpha_sum_type',enabled_when='fal_bool'),
        tua.Item('alpha_tsum_range',enabled_when='fal_bool'),
        tua.Item('FF_combine'),
        tua.Item('FF_combine_type',enabled_when='FF_combine'),
        tua.Item('FF_index_list',enabled_when='FF_combine',style='custom'),
        tua.Item('FF_coeff_list',enabled_when='FF_combine',style='custom'),
        tua.Item('Apply', show_label=False),
        tua.Item('CreateJob', show_label=False ),
        tua.Item('CreateSeparateJob', show_label=False ),
        tua.Item('RUN', show_label=False ),
        tua.Item('Plot', show_label=False ),
        tua.Item('PlotOnlyDist', show_label=False ),
        label='Page 3'
        ),
        tua.Group(
        tua.Item('Fit_Schiff'),
        tua.Item('twod_fit',enabled_when='not_Fit_Schiff'),
        tua.Item('mixed_term',enabled_when='twod_fit'),
        tua.Item('Vanish_Chiral_Limit',enabled_when='not_Fit_Schiff'),
        tua.Item('Order_2_Fit',enabled_when='not_Fit_Schiff'),
        tua.Item('Apply', show_label=False),
        tua.Item('CreateJob', show_label=False ),
        tua.Item('CreateSeparateJob', show_label=False ),
        tua.Item('RUN', show_label=False ),
        tua.Item('Plot', show_label=False ),
        tua.Item('PlotOnlyDist', show_label=False ),
        label='Extrapolation Type'
        ),
        resizable=True
    )
    def _Fit_Schiff_changed(self):
        self.not_Fit_Schiff = not self.Fit_Schiff
        if self.Fit_Schiff:
            self.twod_fit = False
        self._twod_fit_changed()

    def _twod_fit_changed(self):
        if not self.twod_fit:
            self.mixed_term = False

    def _combine_sinks_changed(self):
        self.not_combine_sinks = not self.combine_sinks


    def _full_ratios_changed(self):
        self.frat_bool = self.full_ratios

    def _full_alpha_changed(self):
        self.fal_bool = self.full_alpha


    def _current_type_changed(self):
        if 'Top' in self.current_type or 'Wein' in self.current_type:
            if 'Top' in self.current_type:
                self.observable = 'TopCharge'
            else:
                self.observable = 'Weinberg'
            self.njnq_observable = True
            if self.full_ratios:
                self.frat_bool = True
            if self.full_alpha:
                self.fal_bool = True
        else:
            self.njnq_observable = False
            self.frat_bool = False
            self.fal_bool = False
        self.Rat_comb_alpha = '_mal' in self.current_type



    def _Load_fired(self):
        self.combine_sinks = 'sm2ptRat' in list(thisInfoFlowList[info_index].keys())
        try:
            if self.combine_sinks:
                self.smlist       = thisInfoFlowList[info_index]['sm2ptRat']
                self.sink_coeffs  =  thisInfoFlowList[info_index]['jsmCoeffs']
                self.sink_label   = thisInfoFlowList[info_index]['sinklab']
        except Exception as err:
            print(str(err))

        self.EffM_state_fit = thisInfoFlowList[info_index]['Fits'][0][0]
        self.EffM_fit_range = thisInfoFlowList[info_index]['Fits'][0][1]

        self.source_smearing           = thisInfoFlowList[info_index]['ism']
        if 'Rat_alpha' in thisInfoFlowList[info_index]:
            self.Rat_alpha                 = thisInfoFlowList[info_index]['Rat_alpha']
        self.t_source                  = thisInfoFlowList[info_index]['tsrc']
        self.sink_smearing             = thisInfoFlowList[info_index]['jsm']
        if 'AlphaFree' in list(thisInfoFlowList[info_index].keys()):
            self.all_cfgs_for_alpha        = thisInfoFlowList[info_index]['AlphaFree']
        self.t_sink               = int(thisInfoFlowList[info_index]['tsink'])
        self.doublet_singlet_comb = thisInfoFlowList[info_index]['DoubSing']
        if 'CurrType' in list(thisInfoFlowList[info_index].keys()):
            self.current_type         = thisInfoFlowList[info_index]['CurrType']
        if 'Rat_comb_alpha_fun' in list(thisInfoFlowList[info_index].keys()):
            if str(thisInfoFlowList[info_index]['Rat_comb_alpha_fun']) in op_str_list:
                self.Rat_comb_alpha_fun         = thisInfoFlowList[info_index]['Rat_comb_alpha_fun']

        if 'FFnumb' in list(thisInfoFlowList[info_index].keys()):
            this_ff = thisInfoFlowList[info_index]['FFnumb'].replace('frac{F_{3}}{2m_{N}}','F_{3}/2m_{N}')
            this_ff = this_ff.replace('frac{tilde{F}_{3}}{2m_{N}}','tilde{F}_{3}/2m_{N}')
            self.form_factor          = this_ff
        if 'q2list' in list(thisInfoFlowList[info_index].keys()):
            self.q2_momentum_list = thisInfoFlowList[info_index]['q2list']
        if 'pp2' in list(thisInfoFlowList[info_index].keys()):
            self.pp2_momentum         = thisInfoFlowList[info_index]['pp2']
        if 'Alpha_Full' in list(thisInfoFlowList[info_index].keys()):
            self.fal_bool            = thisInfoFlowList[info_index]['Alpha_Full']
        if 'Full_Rat' in list(thisInfoFlowList[info_index].keys()):
            self.frat_bool           = thisInfoFlowList[info_index]['Full_Rat']
        if 'Chi_Cutoff' in list(thisInfoFlowList[info_index].keys()):
            self.Chi_Cutoff          = thisInfoFlowList[info_index]['Chi_Cutoff']
        if 'RF_fit_range_minimum' in list(thisInfoFlowList[info_index].keys()):
            self.RF_fit_range_minimum = thisInfoFlowList[info_index]['RF_fit_range_minimum']
        if 'RF_fit_min' in list(thisInfoFlowList[info_index].keys()):
            self.RF_fit_min           = thisInfoFlowList[info_index]['RF_fit_min']
        if 'RF_fit_max' in list(thisInfoFlowList[info_index].keys()):
            self.RF_fit_max           = thisInfoFlowList[info_index]['RF_fit_max']

        self.alpha_fit_range = thisInfoFlowList[info_index]['FitsAlpha'][0]
        if self.fal_bool:
            if 'op_trange_sum_type' in list(thisInfoFlowList[info_index].keys()):
                self.alpha_sum_type = '_'.join(thisInfoFlowList[info_index]['op_trange_sum_type'].split('_')[:-1])
            if 'FitsAlphaFull' in list(thisInfoFlowList[info_index].keys()):
                self.alpha_tsum_range = thisInfoFlowList[info_index]['FitsAlphaFull'][0][1]

        if self.frat_bool:
            if 'op_trange_sum_type_3pt' in list(thisInfoFlowList[info_index].keys()):
                self.Rat_sum_type = thisInfoFlowList[info_index]['op_trange_sum_type_3pt'].replace('_None','')
            if 'Rat_tsum_range' in list(thisInfoFlowList[info_index].keys()):
                self.Rat_tsum_range = thisInfoFlowList[info_index]['Rat_tsum_range']
            # thisInfoFlowList[info_index]['FitsRF'] = [str(self.rf_fit_range)]

        if bool(self.Chi_Cutoff):
            if 'Colapse_Chi_Dist' in list(thisInfoFlowList[info_index].keys()):
                self.Colapse_Chi_Dist = thisInfoFlowList[info_index]['Colapse_Chi_Dist']
            if 'Important_FF' in list(thisInfoFlowList[info_index].keys()):
                self.Only_Important_Ratios = thisInfoFlowList[info_index]['Important_FF']
            if 'Reduce_System' in list(thisInfoFlowList[info_index].keys()):
                self.Reduce_System = thisInfoFlowList[info_index]['Reduce_System']
            if 'Chi2pdf_Min' in list(thisInfoFlowList[info_index].keys()):
                self.Chi2pdf_Min = thisInfoFlowList[info_index]['Chi2pdf_Min']
            if 'Chi2pdf_Max' in list(thisInfoFlowList[info_index].keys()):
                self.Chi2pdf_Max = thisInfoFlowList[info_index]['Chi2pdf_Max']
            if 'Chi2pdf_number_picked' in list(thisInfoFlowList[info_index].keys()) :
                if '-' in thisInfoFlowList[info_index]['Chi2pdf_number_picked']:
                    self.Chi2pdf_number_picked = thisInfoFlowList[info_index]['Chi2pdf_number_picked']

        self.njnq_observable = 'Observable' in list(thisInfoFlowList[info_index].keys())
        if self.njnq_observable:
            self.tflow_plot_min = thisInfoFlowList[info_index]['tflowlist'][0]
            self.tflow_plot_max = thisInfoFlowList[info_index]['tflowlist'][-1]
            if len(thisInfoFlowList[info_index]['tflowlist']) == 1:
                self.tflow_plot_inc = 1
            else:
                self.tflow_plot_inc = thisInfoFlowList[info_index]['tflowlist'][1] - self.tflow_plot_min
            if 'TopCharge' in thisInfoFlowList[info_index]['Observable']:
                self.observable=                   'TopCharge'
            elif 'Weinberg' in thisInfoFlowList[info_index]['Observable']:
                self.observable=                   'Weinberg'
            self.tflow_fit = thisInfoFlowList[info_index]['tflowfit'][0]


    def _Apply_fired(self):
        thisInfoFlowList[info_index]['ism'] = self.source_smearing
        thisInfoFlowList[info_index]['Rat_alpha'] = self.Rat_alpha
        thisInfoFlowList[info_index]['tsrc'] = self.t_source
        thisInfoFlowList[info_index]['jsm'] = self.sink_smearing
        thisInfoFlowList[info_index]['AlphaFree'] = self.all_cfgs_for_alpha
        if self.combine_sinks:
            thisInfoFlowList[info_index]['sm2ptRat'] = list(map(int,self.smlist))
            thisInfoFlowList[info_index]['jsmCoeffs'] = list(map(float,self.sink_coeffs))
            thisInfoFlowList[info_index]['sinklab'] = str(self.sink_label)
        else:
            if 'sm2ptRat' in list(thisInfoFlowList[info_index].keys()):
                del thisInfoFlowList[info_index]['sm2ptRat']

        thisInfoFlowList[info_index]['Fits'] = [(str(self.EffM_state_fit),str(self.EffM_fit_range))]
        thisInfoFlowList[info_index]['tsink'] = str(self.t_sink)
        thisInfoFlowList[info_index]['DoubSing'] = str(self.doublet_singlet_comb)
        thisInfoFlowList[info_index]['Interp'] = 'CPEven'
        thisInfoFlowList[info_index]['CurrType'] = str(self.current_type)
        thisInfoFlowList[info_index]['Rat_comb_alpha_fun'] = str(self.Rat_comb_alpha_fun)
        thisInfoFlowList[info_index]['FFnumb'] = str(self.form_factor)
        if '/2m_{N}' in self.form_factor:
            thisInfoFlowList[info_index]['FFnumb'] = r'frac{'+thisInfoFlowList[info_index]['FFnumb'].replace('/',r'}{')+'}'
        thisInfoFlowList[info_index]['q2list'] = list(map(int,self.q2_momentum_list))
        thisInfoFlowList[info_index]['pp2'] = int(self.pp2_momentum)

        thisInfoFlowList[info_index]['Alpha_Full'] = bool(self.fal_bool)
        if self.fal_bool:
            thisInfoFlowList[info_index]['op_trange_sum_type'] = str(self.alpha_sum_type)+'_None'
            thisInfoFlowList[info_index]['FitsAlphaFull'] = [(str(self.alpha_fit_range),str(self.alpha_tsum_range))]
            thisInfoFlowList[info_index]['FitsAlpha'] = [str(self.alpha_fit_range)]
        else:
            thisInfoFlowList[info_index]['op_trange_sum_type'] = 'None'
            thisInfoFlowList[info_index]['FitsAlpha'] = [str(self.alpha_fit_range)]

        thisInfoFlowList[info_index]['Full_Rat'] = bool(self.frat_bool)
        if self.frat_bool:
            thisInfoFlowList[info_index]['op_trange_sum_type_3pt'] = str(self.Rat_sum_type)+'_None'
            thisInfoFlowList[info_index]['Rat_tsum_range'] = str(self.Rat_tsum_range)
        else:
            thisInfoFlowList[info_index]['op_trange_sum_type_3pt'] = 'None'
            thisInfoFlowList[info_index]['Rat_tsum_range'] = 'None'
            # thisInfoFlowList[info_index]['FitsRF'] = [str(self.rf_fit_range)]

        thisInfoFlowList[info_index]['Chi_Cutoff'] = bool(self.Chi_Cutoff)
        if bool(self.Chi_Cutoff):

            thisInfoFlowList[info_index]['Colapse_Chi_Dist'] = bool(self.Colapse_Chi_Dist)
            thisInfoFlowList[info_index]['Important_FF'] = bool(self.Only_Important_Ratios)
            thisInfoFlowList[info_index]['Chi2pdf_Min'] = float(self.Chi2pdf_Min)
            thisInfoFlowList[info_index]['Chi2pdf_Max'] = float(self.Chi2pdf_Max)
            thisInfoFlowList[info_index]['Chi2pdf_number_picked'] = self.Chi2pdf_number_picked
        else:
            thisInfoFlowList[info_index]['Important_FF'] = False
            thisInfoFlowList[info_index]['Colapse_Chi_Dist'] = False
            thisInfoFlowList[info_index]['Chi2pdf_number_picked'] = 'Chi0-Chi100'
        thisInfoFlowList[info_index]['Reduce_System'] = bool(self.Reduce_System)
        thisInfoFlowList[info_index]['RF_fit_range_minimum'] = int(self.RF_fit_range_minimum)
        thisInfoFlowList[info_index]['RF_fit_min'] = int(self.RF_fit_min)
        thisInfoFlowList[info_index]['RF_fit_max'] = int(self.RF_fit_max)


        if self.njnq_observable:
            thisInfoFlowList[info_index]['tflowlist'] = np.arange(max(self.tflow_plot_min,0.01),min(self.tflow_plot_max+self.tflow_plot_inc,10),self.tflow_plot_inc)
            thisInfoFlowList[info_index]['Observable'] = self.observable
            thisInfoFlowList[info_index]['tflowfit']= [str(self.tflow_fit)]
        else:
            if 'Observable' in list(thisInfoFlowList[info_index].keys()):
                del thisInfoFlowList[info_index]['Observable']

    def _CreateJob_fired(self):
        print()
        global thisDefWipe,Mcore,JobNumber
        SavePreviousRun()

        this_job_file = JobFolder+'FFSetup'+str(JobNumber)+'.py3p'
        IO.WritePickle(this_job_file,[  thisInfoFlowList,thisDefWipe,
                                        self.recalculate_all,self.config_list_from_file,
                                        Mcore,self.CommonRatios])
        next_job_file = JobFolder+'FFSolve'+str(JobNumber)+'.csh'
        WriteCSH('AnaFFSetup',this_job_file,job_params,JobNumber,Next_Job=next_job_file)
        print('Created job script for Form Factor setup', this_job_file.replace('.py3p','.csh'))

        this_job_file = JobFolder+'FFSolve'+str(JobNumber)+'.py3p'
        IO.WritePickle( this_job_file,[thisInfoFlowList,thisDefWipe,self.recalculate_all,
                        self.config_list_from_file,Mcore,self.CommonRatios])
        next_job_file = JobFolder+'FFFinish'+str(JobNumber)+'.csh'
        WriteCSH('AnaFFSolve',this_job_file,job_params_solve,JobNumber,Next_Job=next_job_file)
        print('Created job script for Form Factor solve', this_job_file.replace('.py3p','.csh'))

        this_job_file = JobFolder+'FFFinish'+str(JobNumber)+'.py3p'
        IO.WritePickle(this_job_file,[  thisInfoFlowList,thisDefWipe,
                                        self.recalculate_all,self.config_list_from_file,
                                        Mcore,False,False,True,self.CommonRatios])
        WriteCSH('AnaFFFinish',this_job_file,job_params,JobNumber)
        print('Created job script for Form Factor run', this_job_file.replace('.py3p','.csh'))


        JobNumber += 1

    def _CreateSeparateJob_fired(self):
        print()
        global thisDefWipe,Mcore,JobNumber
        SavePreviousRun()
        for iinfo in thisInfoFlowList:


            this_job_file = JobFolder+'FFSetup'+str(JobNumber)+'.py3p'
            IO.WritePickle(this_job_file,[  [iinfo],thisDefWipe,self.recalculate_all,
                                            self.config_list_from_file,Mcore,self.CommonRatios])
            next_job_file = JobFolder+'FFSolve'+str(JobNumber)+'.csh'
            WriteCSH('AnaFFSetup',this_job_file,job_params,JobNumber,Next_Job=next_job_file)
            print('Created job script for Form Factor setup', this_job_file.replace('.py3p','.csh'))

            this_job_file = JobFolder+'FFSolve'+str(JobNumber)+'.py3p'
            IO.WritePickle(this_job_file,[  [iinfo],thisDefWipe,self.recalculate_all,
                                            self.config_list_from_file,Mcore,self.CommonRatios])
            next_job_file = JobFolder+'FFFinish'+str(JobNumber)+'.csh'
            WriteCSH('AnaFFSolve',this_job_file,job_params_solve,JobNumber,Next_Job=next_job_file)
            print('Created job script for Form Factor solve', this_job_file.replace('.py3p','.csh'))


            this_job_file = JobFolder+'FFFinish'+str(JobNumber)+'.py3p'
            IO.WritePickle( this_job_file,[[iinfo],thisDefWipe,self.recalculate_all,
                            self.config_list_from_file,Mcore,False,False,True,self.CommonRatios])
            WriteCSH('AnaFFFinish',this_job_file,job_params,JobNumber)
            print('Created job script for Form Factor run', this_job_file.replace('.py3p','.csh'))


            JobNumber += 1
            # this_job_file = JobFolder+'FF'+str(JobNumber)+'.py3p'
            # IO.WritePickle(this_job_file,[[iinfo],thisDefWipe,self.recalculate_all,self.config_list_from_file,Mcore])
            # WriteCSH('AnaFF',this_job_file,job_params,JobNumber)
            # print 'Created job script for Form Factor run', this_job_file.replace('.py3p','.csh')
            # JobNumber += 1




    def _RUN_fired(self):
        # for ikey,ival in thisInfoFlowList[0].iteritems():
        #     print ikey,ival
        print()
        print('Begining Form Factor run')
        global thisDefWipe,Mcore,data,plot_info
        SavePreviousRun()
        if any([iFF > len(thisInfoFlowList) for iFF in self.FF_index_list]):
            raise EnvironmentError('FF_index_list has index which is out side the number of sets')
        this_chilist = [ival['Chi2pdf_number_picked'] for ival in thisInfoFlowList]
        data = FFRun(   thisInfoFlowList,thisDefWipe,self.recalculate_all,self.config_list_from_file,
                        Mcore,self.CommonRatios)
        # if self.config_list_from_file:
        #     data = ffs.SetOfFF(cfglist='fromfile',cfglist2pt='fromfile',InfoDict=thisInfoFlowList)
        # else:
        #     data = ffs.SetOfFF(InfoDict=thisInfoFlowList)
        #
        # data.LoadPickle(DefWipe=thisDefWipe,HardWipe=self.recalculate_all,Mcore=Mcore)
        tflowlist = [str(idict['tflowfit'][0]) for idict in thisInfoFlowList]
        # data.Plot(self.form_factor,plottflowlist=[[itflow] for itflow in tflowlist],tsize=tsize,xsize=xsize,ysize=ysize)
        FFlist = [str(idict['FFnumb']) for idict in thisInfoFlowList]
        Colape_Chi_list = [idict['Colapse_Chi_Dist'] for idict in thisInfoFlowList]

        if self.FF_combine:
            data.CombineFFs(self.FF_index_list,list(map(np.float64,self.FF_coeff_list)),this_fun=self.FF_combine_type)

        # if Wipe_Fit:
        #     data.WipeAllFits()
        #     data.Write()
        data.FitFFs(thisFFList=FFlist)
        print('Plotting Form Factors')
        data.Plot(FFlist,plot_info,plottflowlist=[[itflow] for itflow in tflowlist],chi_list=this_chilist)
        print('Plotting Form Factor fits vs Mpi')
        data.PlotFitVsMpi(FFlist,plot_info,plottflowlist=tflowlist,Vanish_Clim=self.Vanish_Chiral_Limit,
                            Ord_Two=self.Order_2_Fit,latspace_fit=(self.twod_fit,self.mixed_term),fitSchiff=self.Fit_Schiff)
        print('Plotting Form Factor fits vs UDMass')
        data.PlotFitVsUDMass(FFlist,plot_info,plottflowlist=tflowlist,chi_list=this_chilist)
        print('Plotting Form Factor fits vs lattice spacing')
        data.PlotFitVsLatSpace(FFlist,plot_info,plottflowlist=tflowlist,mpi_fit=(self.twod_fit,self.mixed_term,
                                                                                 self.Vanish_Chiral_Limit,self.Order_2_Fit)
                                                                                 ,fitSchiff=self.Fit_Schiff)
        print('Plotting Form Factor fits vs box size')
        data.PlotFitVsBoxSize(FFlist,plot_info,plottflowlist=tflowlist)
        print('Plotting Form Factor distribution of fits')
        data.PlotDist(  FFlist,plot_info,plottflowlist=tflowlist,
                        chi_colapse=Colape_Chi_list,Vanish_Clim=self.Vanish_Chiral_Limit,Ord_Two=self.Order_2_Fit,
                        latspace_fit=(self.twod_fit,self.mixed_term))
        # data = None
        print('Form Factor run complete')

    def _Plot_fired(self):
        # for ikey,ival in thisInfoFlowList[0].iteritems():
        #     print ikey,ival
        print()
        print('Plotting Form Factor run')
        global thisDefWipe,Mcore,data,plot_info
        SavePreviousRun()
        this_chilist = [ival['Chi2pdf_number_picked'] for ival in thisInfoFlowList]
        if any([iFF > len(thisInfoFlowList) for iFF in self.FF_index_list]):
            raise EnvironmentError('FF_index_list has index which is out side the number of sets')
        data = FFFinish(thisInfoFlowList,thisDefWipe,self.recalculate_all,self.config_list_from_file,Mcore,
                        False,False,False,self.CommonRatios)
        # if self.config_list_from_file:
        #     data = ffs.SetOfFF(cfglist='fromfile',cfglist2pt='fromfile',InfoDict=thisInfoFlowList)
        # else:
        #     data = ffs.SetOfFF(InfoDict=thisInfoFlowList)
        #
        # data.LoadPickle(DefWipe=thisDefWipe,HardWipe=self.recalculate_all,Mcore=Mcore)
        tflowlist = [str(idict['tflowfit'][0]) for idict in thisInfoFlowList]
        # data.Plot(self.form_factor,plottflowlist=[[itflow] for itflow in tflowlist],tsize=tsize,xsize=xsize,ysize=ysize)
        FFlist = [str(idict['FFnumb']) for idict in thisInfoFlowList]
        Colape_Chi_list = [idict['Colapse_Chi_Dist'] for idict in thisInfoFlowList]
        if self.FF_combine:
            data.CombineFFs(self.FF_index_list,list(map(np.float64,self.FF_coeff_list)),this_fun=self.FF_combine_type)
        # if Wipe_Fit:
        #     data.WipeAllFits()
        #     data.Write()
        # data.FitFFs(thisFFList=FFlist,WipeFit=Wipe_Fit)
        print('Plotting Form Factors')
        data.Plot(FFlist,plot_info,plottflowlist=[[itflow] for itflow in tflowlist],chi_list=this_chilist)
        print('Plotting Form Factor fits vs Mpi')

        data.PlotFitVsMpi(FFlist,plot_info,plottflowlist=tflowlist,Vanish_Clim=self.Vanish_Chiral_Limit,
                            Ord_Two=self.Order_2_Fit,latspace_fit=(self.twod_fit,self.mixed_term),fitSchiff=self.Fit_Schiff)
        print('Plotting Form Factor fits vs UDMass')
        data.PlotFitVsUDMass(FFlist,plot_info,plottflowlist=tflowlist,chi_list=this_chilist)
        print('Plotting Form Factor fits vs lattice spacing')
        data.PlotFitVsLatSpace(FFlist,plot_info,plottflowlist=tflowlist,mpi_fit=(self.twod_fit,self.mixed_term,
                                                                                 self.Vanish_Chiral_Limit,self.Order_2_Fit)
                                                                                 ,fitSchiff=self.Fit_Schiff)
        print('Plotting Form Factor fits vs box size')
        data.PlotFitVsBoxSize(FFlist,plot_info,plottflowlist=tflowlist)
        print('Plotting Form Factor distribution of fits')
        data.PlotDist(  FFlist,plot_info,plottflowlist=tflowlist,
                        chi_colapse=Colape_Chi_list,Vanish_Clim=self.Vanish_Chiral_Limit,Ord_Two=self.Order_2_Fit,
                        latspace_fit=(self.twod_fit,self.mixed_term))
        # data = None
        print('Form Factor plot complete')


    def _PlotOnlyDist_fired(self):
        # for ikey,ival in thisInfoFlowList[0].iteritems():
        #     print ikey,ival
        print()
        print('Plotting Form Factor run')
        global thisDefWipe,Mcore,data,plot_info
        SavePreviousRun()
        if self.FF_combine:
            raise EnvironmentError('Cannot combine form factors when loading only fits, please choose Plot or Run')
        # if any([iFF > len(thisInfoFlowList) for iFF in self.FF_index_list]):
        #     raise EnvironmentError('FF_index_list has index which is out side the number of sets')
        this_chilist = [ival['Chi2pdf_number_picked'] for ival in thisInfoFlowList]
        data = FFFinish(thisInfoFlowList,thisDefWipe,self.recalculate_all,self.config_list_from_file,
                        Mcore,False,True,False,self.CommonRatios)
        # if self.config_list_from_file:
        #     data = ffs.SetOfFF(cfglist='fromfile',cfglist2pt='fromfile',InfoDict=thisInfoFlowList)
        # else:
        #     data = ffs.SetOfFF(InfoDict=thisInfoFlowList)
        #
        # data.LoadPickle(DefWipe=thisDefWipe,HardWipe=self.recalculate_all,Mcore=Mcore)
        tflowlist = [str(idict['tflowfit'][0]) for idict in thisInfoFlowList]
        # data.Plot(self.form_factor,plottflowlist=[[itflow] for itflow in tflowlist],tsize=tsize,xsize=xsize,ysize=ysize)
        FFlist = [str(idict['FFnumb']) for idict in thisInfoFlowList]
        Colape_Chi_list = [idict['Colapse_Chi_Dist'] for idict in thisInfoFlowList]
        # if self.FF_combine:
        #     data.CombineFFs(self.FF_index_list,map(np.float64,self.FF_coeff_list),this_fun=self.FF_combine_type)
        # if Wipe_Fit:
        #     data.WipeAllFits()
        #     data.Write()
        # data.FitFFs(thisFFList=FFlist,WipeFit=Wipe_Fit)
        # print 'Plotting Form Factors'
        # data.Plot(FFlist,plot_info,plottflowlist=[[itflow] for itflow in tflowlist],chi_list=this_chilist)
        data.PlotFitVsMpi(FFlist,plot_info,plottflowlist=tflowlist,Vanish_Clim=self.Vanish_Chiral_Limit,
                            Ord_Two=self.Order_2_Fit,latspace_fit=(self.twod_fit,self.mixed_term),fitSchiff=self.Fit_Schiff)
        print('Plotting Form Factor fits vs UDMass')
        data.PlotFitVsUDMass(FFlist,plot_info,plottflowlist=tflowlist,chi_list=this_chilist)
        print('Plotting Form Factor fits vs lattice spacing')
        data.PlotFitVsLatSpace(FFlist,plot_info,plottflowlist=tflowlist,mpi_fit=(self.twod_fit,self.mixed_term,
                                                                                 self.Vanish_Chiral_Limit,self.Order_2_Fit)
                                                                                 ,fitSchiff=self.Fit_Schiff)
        print('Plotting Form Factor fits vs box size')
        data.PlotFitVsBoxSize(FFlist,plot_info,plottflowlist=tflowlist)
        print('Plotting Form Factor distribution of fits')
        data.PlotDist(  FFlist,plot_info,plottflowlist=tflowlist,
                        chi_colapse=Colape_Chi_list,Vanish_Clim=self.Vanish_Chiral_Limit,Ord_Two=self.Order_2_Fit,
                        latspace_fit=(self.twod_fit,self.mixed_term))
        # print 'Plotting Form Factor fits vs Mpi'
        # data.PlotFitVsMpi(FFlist,plot_info,plottflowlist=tflowlist,Vanish_Clim=self.Vanish_Chiral_Limit,Ord_Two=self.Order_2_Fit)
        # print 'Plotting Form Factor fits vs UDMass'
        # data.PlotFitVsUDMass(FFlist,plot_info,plottflowlist=tflowlist,chi_list=this_chilist)
        # print 'Plotting Form Factor fits vs lattice spacing'
        # data.PlotFitVsLatSpace(FFlist,plot_info,plottflowlist=tflowlist)
        # print 'Plotting Form Factor fits vs box size'
        # data.PlotFitVsBoxSize(FFlist,plot_info,plottflowlist=tflowlist)
        # print 'Plotting Form Factor distribution of fits'
        # data.PlotDist(  FFlist,plot_info,plottflowlist=tflowlist,
        #                 chi_colapse=Colape_Chi_list,Vanish_Clim=self.Vanish_Chiral_Limit,Ord_Two=self.Order_2_Fit)
        # data = None
        print('Form Factor plot complete')


class DictFlowEditor(ta.HasTraits):
    '''

    Manual display of all the parameters passed into the analysis (for flowed analysis).

    '''

    Object = ta.Instance( object )
    this_index = ta.Int(0)
    def __init__(self, **traits):
        global thisInfoFlowList
        super(DictFlowEditor, self).__init__(**traits)
        self.Object = thisInfoFlowList[self.this_index]

    def _this_index_changed(self):
        global thisInfoFlowList
        self.Object = thisInfoFlowList[self.this_index]

    def trait_view(self, name=None, view_elements=None):
        return tua.View(
          tua.VGroup(
            tua.Item('Object',
                  label      = 'Debug',
                  id         = 'debug',
                  editor     = tua.ValueEditor(), #ValueEditor()
                  # style      = 'custom',
                  dock       = 'horizontal',
                  show_label = False),),
          title     = 'Flow Information Dictionary Number '+str(self.this_index+1),
          width     = 400,
          height    = 800,
          resizable = True)


class DictEditor(ta.HasTraits):
    '''

    Manual display of all the parameters passed into the analysis.

    '''
    Object = ta.Instance( object )
    this_index = ta.Int(0)
    def __init__(self, **traits):
        global thisInfoList
        super(DictEditor, self).__init__(**traits)
        self.Object = thisInfoList[self.this_index]

    def _this_index_changed(self):
        global thisInfoList
        self.Object = thisInfoList[self.this_index]
        self.trait_view()

    def trait_view(self, name=None, view_elements=None):
        return tua.View(
          tua.VGroup(
            tua.Item('Object',
                  label      = 'Debug',
                  id         = 'debug',
                  editor     = tua.ValueEditor(), #ValueEditor()
                  style      = 'custom',
                  dock       = 'horizontal',
                  show_label = False),),
          title     = 'Information Dictionary Number '+str(self.this_index+1),
          width     = 400,
          height    = 800,
          resizable = True)

class BeginFrame( ta.HasTraits ):
    """

    Master frame.

    Internal state of 'run parameters' is stored in thisInfoList and thisInfoFlowList

    Create and Delete parameter sets elements via 'new info' and 'delete info' buttons

    'show info diff' will display all the parameters that vary between all the
    different parameter sets set up for the analysis.

    'load previous run' will load the previously calculated run again.

    """

    Info_List = ta.List([1])
    Info_Index = ta.Enum(values='Info_List')
    General_Params = ta.Button()
    Run_Params = ta.Button()
    Plot_Params = ta.Button()
    Job_Params = ta.Button()
    Job_Params_Solve = ta.Button()
    Flowed_Opps = ta.Button()
    Flowed_Opps_Full = ta.Button()
    Two_Point_Corrs = ta.Button()
    Alpha_Ratio = ta.Button()
    Alpha_Ratio_Full = ta.Button()
    Ratio_Function = ta.Button()
    Form_Factor = ta.Button()
    Load_Previous_Run = ta.Button()
    New_Info = ta.Button()
    Delete_Info = ta.Button()
    Show_Info_Diff = ta.Button()
    Show_Info = ta.Button()
    Show_Info_Flow = ta.Button()
    Info_Length = ta.Str('Number of Infos: '+str(len(thisInfoList)))
    display = ta.String(' Please check parameters, then select analysis. ')

    Info_window = DictEditor()
    Info_Flow_window = DictFlowEditor()
    GP_window = CalculationParams()
    Run_window = RunParams()
    Plot_window = PlotParams()
    Job_window = JobParams()
    Job_Solve_window = JobParamsSolve()
    TPC_window = TPCorr()
    FO_window = FlowedOp()
    FOF_window = FlowedOpFull()
    Alpha_window = Alpha()
    AlphaFull_window = AlphaFull()
    RF_window = RatioFun()
    FF_window = FormFacts()


    traits_view = tua.View(
        '_',
        tua.Item('display', show_label=False,style='readonly' ),
        '_',
        tua.Item('Info_Length', show_label=False, style='readonly' ),
        tua.Item('Info_Index'),
        tua.Item('New_Info', show_label=False),
        tua.Item('Delete_Info', show_label=False),
        tua.Item('Show_Info_Diff', show_label=False),
        tua.Item('Show_Info', show_label=False),
        tua.Item('Show_Info_Flow', show_label=False),
        '_',
        '_',
        tua.Item('General_Params', show_label=False ),
        tua.Item('Run_Params', show_label=False ),
        tua.Item('Plot_Params', show_label=False ),
        tua.Item('Job_Params', show_label=False ),
        tua.Item('Job_Params_Solve', show_label=False ),
        '_',
        tua.Item('Flowed_Opps', show_label=False ),
        tua.Item('Flowed_Opps_Full', show_label=False ),
        tua.Item('Two_Point_Corrs', show_label=False ),
        tua.Item('Alpha_Ratio', show_label=False ),
        tua.Item('Alpha_Ratio_Full', show_label=False ),
        tua.Item('Ratio_Function', show_label=False ),
        tua.Item('Form_Factor', show_label=False ),
        '_',
        tua.Item('Load_Previous_Run', show_label=False ),
        '_',
        buttons = ['Cancel'],
        title = 'Analysis',
        resizable=True
            )

    def _Info_Index_changed(self):
        self.Info_List = list(range(1,len(thisInfoList)+1))
        if self.Info_Index not in self.Info_List:
            self.Info_Index = self.Info_List[-1]
        if self.Info_Index < 1:
            self.Info_Index = 1
        global info_index
        info_index = int(self.Info_Index)-1
        self.Info_window.this_index = info_index
        self.Info_Flow_window.this_index = info_index
        self.TPC_window._Load_fired()
        self.FO_window._Load_fired()
        self.FOF_window._Load_fired()
        self.RF_window._Load_fired()
        self.Alpha_window._Load_fired()
        self.AlphaFull_window._Load_fired()
        self.FF_window._Load_fired()


    def _Info_Length_changed(self):
        self._Info_Index_changed()

    def __init__(self):
        self._Load_Previous_Run_fired()
        self._Info_Index_changed()
        # self._Show_Info_fired()
        # self._Show_Info_Flow_fired()

    def _New_Info_fired(self):
        thisInfoList.append(deepcopy(thisInfoList[-1]))
        thisInfoFlowList.append(deepcopy(thisInfoFlowList[-1]))
        # self.Info_Index = self.Info_Index + 1
        self.Info_Length = 'Number of Infos: '+str(len(thisInfoList))
        self.Info_Index = self.Info_List[-1]

    def _Delete_Info_fired(self):
        global info_index
        if len(thisInfoList) > 1:
            del thisInfoList[info_index]
            self.Info_Length = 'Number of Infos: '+str(len(thisInfoList))
            self.Info_Index = self.Info_List[-1]
        if len(thisInfoFlowList) > 1:
            del thisInfoFlowList[info_index]
            self.Info_Length = 'Number of Infos Flows: '+str(len(thisInfoFlowList))
            self.Info_Index = self.Info_List[-1]

    def _Show_Info_Diff_fired(self):
        if len(thisInfoList) > 1:
            print()
            print('Info Dict Differences:')
            print()
            diffkeylist = []
            keylist = list(thisInfoList[0].keys())
            for idict in thisInfoList[1:]:
                for ikey,ival in idict.items():
                    if ikey in keylist and ikey not in diffkeylist:
                        thischeck = ival != thisInfoList[0][ikey]
                        if isinstance(thischeck,(list,tuple,np.ndarray)):
                            thischeck = any(thischeck)
                        if thischeck:
                            diffkeylist.append(ikey)
            for ikey in diffkeylist:
                print()
                print(ikey)
                for ic,idict in enumerate(thisInfoList):
                    print('    ', ic , idict[ikey])
        if len(thisInfoFlowList) > 1:
            print()
            print('Flow Info Dict Differences:')
            print()
            diffkeylist = []
            keylist = list(thisInfoFlowList[0].keys())
            for idict in thisInfoFlowList[1:]:
                for ikey,ival in idict.items():
                    if ikey in keylist and ikey not in diffkeylist:
                        thischeck = ival != thisInfoFlowList[0][ikey]
                        if isinstance(thischeck,(list,tuple,np.ndarray)):
                            thischeck = any(thischeck)
                        if thischeck:
                            diffkeylist.append(ikey)
            for ikey in diffkeylist:
                print()
                print(ikey)
                for ic,idict in enumerate(thisInfoFlowList):
                    print('    ', ic , idict[ikey])
        self.Info_Length = 'Number of Infos: '+str(len(thisInfoList))


    def _Show_Info_fired(self):
        # for icd,idict in enumerate(thisInfoList):
        # global info_index
        # print
        # print 'Info Dict #',info_index+1
        # print
        # for ikey,ival in thisInfoList[info_index].iteritems():
        #     print '  ',ikey, ReductList(ival)
        self.Info_window.configure_traits()
        self.Info_Length = 'Number of Infos: '+str(len(thisInfoList))


    def _Show_Info_Flow_fired(self):
        # global info_index
        # print
        # print 'Info Dict #',info_index+1
        # print
        # for ikey,ival in thisInfoFlowList[info_index].iteritems():
        #     print '  ',ikey, ReductList(ival)
        self.Info_Flow_window.configure_traits()
        self.Info_Length = 'Number of Infos: '+str(len(thisInfoList))


    def _General_Params_fired(self):
        self.GP_window.configure_traits()

    def _Run_Params_fired(self):
        self.Run_window.configure_traits()

    def _Plot_Params_fired(self):
        self.Plot_window.configure_traits()

    def _Job_Params_fired(self):
        self.Job_window.configure_traits()

    def _Job_Params_Solve_fired(self):
        self.Job_Solve_window.configure_traits()

    def _Flowed_Opps_fired(self):
        self.FO_window.configure_traits()

    def _Flowed_Opps_Full_fired(self):
        self.FOF_window.configure_traits()


    def _Two_Point_Corrs_fired(self):
        self.TPC_window.configure_traits()

    def _Alpha_Ratio_fired(self):
        self.Alpha_window.configure_traits()

    def _Alpha_Ratio_Full_fired(self):
        self.AlphaFull_window.configure_traits()

    def _Ratio_Function_fired(self):
        self.RF_window.configure_traits()

    def _Form_Factor_fired(self):
        self.FF_window.configure_traits()

    def _Load_Previous_Run_fired(self):
        LoadPreviousRun()
        self.Info_Length = 'Number of Infos: '+str(len(thisInfoList))
        # rreload(sfo)
        # rreload(so)
        # rreload(sr)
        # rreload(ffs)
        # rreload(ql)


if  __name__ == "__main__":
    Start = BeginFrame()
    Start.configure_traits()






# t
# comp_view = tua.View(
#     Group(
#         Group(
#             tua.Item('h1.display', show_label=False,style='readonly' ),
#             '_',
#             tua.Item('h1.Gen_Params', show_label=False ),
#             tua.Item('h1.Flowed_Opps', show_label=False ),
#             tua.Item('h1.Two_Point_Corrs', show_label=False ),
#             tua.Item('h1.Alpha', show_label=False ),
#             tua.Item('h1.Ratio_Function', show_label=False ),
#             tua.Item('h1.Form_Factor', show_label=False ),
#             '_',
#             show_border=True
#         ),
#         Group(
#             tua.Item('h2.MultiCore'),
#             show_border=True
#         ),
#         orientation = 'horizontal'
#     ),
#     title = 'Analysis',
#     buttons = ['Cancel','Run'],
# )
