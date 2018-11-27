#!/usr/bin/env python


#IMPORT THIS FIRST, put in base import file
# import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
# import matplotlib.pyplot as pl
##

from Params import nboot,defSparam,datadir,outputdir,flow_new_format,Wipe_All_Fits
from FileIO import ReadFuns,WriteFuns,ReadWithMeta,WriteWithMeta
from Params import chitcoeff,defxlimOp,TflowToPhys,Qeps,cfg_file_type,cfgfmtdir
import MomParams as mp
from BootStrapping import BootStrap
from MiscFuns import ODNested,mkdir_p,CheckClass,Series_TO_ODict,fmt_file_type,StreamToNumb
# from XmlFormatting import tstr,untstr,tflowstr,untflowstr,unxmlfitr,AvgStdToFormat,tflowTOLeg
from XmlFormatting import tstr,untstr,tflowstr,untflowstr
from XmlFormatting import unxmlfitr_int,AvgStdToFormat,tflowTOLeg,KeyForamtting
# from XmlFormatting import tstr,tflowstr,untflowstr,unxmlfitr,AvgStdToFormat,tflowTOLeg
from XmlFormatting import MakeValAndErr
from FileIO import WriteXml,WritePickle,ReadPickleWrap,WriteExcel
from TimeStuff import Timer
from copy import deepcopy
from Autocorr import AutoCorrelate
from PredefFitFuns import Chi_Prop_Fit_2exp,LinearFitFun
from PlotData import null_series
# from PredefFitFuns import RandTFitFun
# import cPickle as pickle
# import dill as pickle
import itertools
from copy import copy

import numpy as np
import pandas as pa
import re,os
# import FitFunctions as ff
import SetsOfFits as sff
import glob

from collections import OrderedDict
import operator as ope



"""
How to use:

"""
def GetOppFromStr(thisOp):
    if thisOp == '+':
        return ope.add
    elif thisOp == '-':
        return ope.sub
    elif thisOp == '*' or thisOp == 'x':
        return ope.mul
    elif thisOp == '/':
        return ope.div
    elif thisOp == 'leftI':
        def returnopp(a,b):
            return a
        return returnopp
    else:
        raise IOError('Operation combination not supported ' + thisOp)

def CombineInfo(Info,Info2,Opp):
    InfoOut = copy(Info)
    ## Info can contain fit ranges to be calculated
    ## must bein form Info['FlowFits'] = [ ('state#' , 'fit#-#' ) , ...]


    if 'FlowFits' in list(Info.keys()):
        if 'FlowFits' in list(Info2.keys()):
            InfoOut['FlowFits'] = Info['FlowFits'] + Info2['FlowFits']
        else:
            InfoOut['FlowFits'] = Info['FlowFits']
    elif 'FlowFits' in list(Info2.keys()):
        InfoOut['FlowFits'] = Info2['FlowFits']


    ## Info can contain fit ranges to be calculated
    ## must bein form Info['FlowFitFun'] = (Function Object , number of parameters )
    if 'FlowFitFun' in list(Info.keys()):
        if 'FlowFitFun' in list(Info2.keys()):
            if Info['FlowFitFun'] != Info2['FlowFitFun']: print('Warning, different fit functions present, using host fit funtion')
        InfoOut['FlowFitFun'] = Info['FlowFitFun']
    elif 'FlowFitFun' in list(Info2.keys()):
        InfoOut['FlowFitFun'] = Info2['FlowFitFun']

    ## kappa up down (in integer form)
    if 'kud' in list(Info.keys()):
        if 'kud' in list(Info2.keys()):
            if Info['kud'] != Info2['kud']:
                InfoOut['kud'] = Info['kud'] +Opp+ Info2['kud']
        InfoOut['kud'] = Info['kud']
    elif 'kud' in list(Info2.keys()):
        InfoOut['kud'] = Info2['kud']

    ## kappa strange (in integer form)
    if 'ks' in list(Info.keys()):
        if 'ks' in list(Info2.keys()):
            if Info['ks'] != Info2['ks']:
                InfoOut['ks'] = Info['ks'] +Opp+ Info2['ks']
        InfoOut['ks'] = Info['ks']
    elif 'ks' in list(Info2.keys()):
        InfoOut['ks'] = Info2['ks']

    ## Current obersvables are TopCharge and Weinberg
    if 'Observable' in list(Info.keys()):
        if 'Observable' in list(Info2.keys()):
            InfoOut['Observable'] = Info['Observable'] +Opp+ Info2['Observable']
        else:
            InfoOut['Observable'] = Info['Observable']
    elif 'Observable' in list(Info2.keys()):
        InfoOut['Observable'] = Info2['Observable']


    ## this is used display only these flowtimes, Still reads and writes all
    if 'tflowlist' in list(Info.keys()):
        if 'tflowlist' in list(Info2.keys()):
            if any(np.array(Info['tflowlist']) != np.array(Info2['tflowlist'])): print('Warning, different flow time lists present, using host fit funtion')
        InfoOut['tflowlist'] = Info['tflowlist']
    elif 'tflowlist' in list(Info2.keys()):
        InfoOut['tflowlist'] = Info2['tflowlist']

    ## Info can contain initial guess for fits, make sure it array of correct length
    if 'iGuess' in list(Info.keys()):
        if 'iGuess' in list(Info2.keys()):
            if Info['iGuess'] != Info2['iGuess']: print('Warning, different initial condition found, using host initail conidition')
        InfoOut['iGuess'] = Info['iGuess']
    elif 'iGuess' in list(Info2.keys()):
        InfoOut['iGuess'] = Info2['iGuess']

    if 'Sparam' in list(Info.keys()):
        if 'Sparam' in list(Info2.keys()):
            if Info['Sparam'] != Info2['Sparam']: print('Warning, different Sparams, using host initail conidition')
        InfoOut['Sparam'] = Info['Sparam']
    elif 'Sparam' in list(Info2.keys()):
        InfoOut['Sparam'] = Info2['Sparam']


    if 'stream_list' in list(Info.keys()):
        if 'stream_list' in list(Info2.keys()):
            if isinstance(Info['stream_list'],(list,tuple,np.ndarray)):
                if any(np.array(Info['stream_list']) != Info2['stream_list']):
                    raise EnvironmentError('stream list missmatch')
                else:
                    InfoOut['stream_list'] = Info['stream_list']
            elif isinstance(Info2['stream_list'],(list,tuple,np.ndarray)):
                if any(np.array(Info2['stream_list']) != Info['stream_list']):
                    raise EnvironmentError('stream list missmatch')
                else:
                    InfoOut['stream_list'] = Info['stream_list']
            else:
                if Info2['stream_list'] != Info['stream_list']:
                    raise EnvironmentError('stream list missmatch')
                else:
                    InfoOut['stream_list'] = Info['stream_list']
        else:
            InfoOut['stream_list'] = Info['stream_list']
    elif 'stream_list' in list(Info2.keys()):
        InfoOut['stream_list'] = Info2['stream_list']
    return InfoOut



class FlowOp(object):
    """ Op(tf) uses bootstrap class
        tf = flow time
    """

    ## Info is a dictionary containing information for the correlator:
    ## see comments below to see what is required.
    ## missing elements just ignores that field

    ## TODO: create variational flag and interface with this function.

    ## array of all classes that inherit this class. Append this when adding new classes please!
    parentlist = [  'FlowOpps.FlowOpFullSquared',
                    'TwoPtCorrelators.TwoPointCorr','TwoPtCorrelators.NNQCorr',
                    'TwoPtCorrelators.NNQFullCorr','VarMethod.VariationalTwoPt',
                    'RatioCorrelators.RatioCorr',
                    'RatioCorrelators.RatioFOCorr','RatioCorrelators.RatioFOFullCorr']

    def __init__(self, thisnboot=nboot, cfglist={},Info={}, thissym='Not Set',thiscol='Not Set',thisshift=0.0,name='',LLname='',
                man_load_cfgs=False):
        def GetFolderFileOp(thistype):
            if self.FlowNewForm:
                if  'TopCharge' in thistype  :
                    return 'FlowNewForm','q'
                elif 'Weinberg' in thistype :
                    return 'FlowNewForm','W'
                elif 'Energy' in thistype :
                    return 'FlowNewForm','E'
                else:
                    print('thistype ' , thistype)
                    raise IOError('Type of operator not known, please select TopCharge or Weinberg')
            else:
                if  'TopCharge' in thistype  :
                    return 'topcharge','q'
                elif 'Weinberg' in thistype :
                    return 'weinopp','W'
                else:
                    print('thistype ' , thistype)
                    raise IOError('Type of operator not known, please select TopCharge or Weinberg')

        ## these class elemnts are checked when loading in pickled file

        # self.thischecklist = ['nboot','kud','ks','tflowlist','cfglist']
        self.thischecklist = ['nboot','kud','ks','tflowlist','Sparam','FlowNewForm']

        self.flowcurr = 0
        ## Initialising plotting parameters
        self.thisshift = thisshift
        self.thiscol = thiscol
        self.thissym = thissym



        self.Op_cfgs  = pa.DataFrame()
        '''
        self.Op_cfgs DataFrame:

        columns =   config_numbers , Op[itflow]

        rows =      stream , configuration
        '''

        self.Op_Stats = pa.DataFrame()

        '''
        self.Op_Stats DataFrame:

        columns =   boot, Avg,      Std
                    Auto, AutoAvg,  AutoStd

        rows =      flow time
        '''
        self.Op_Stats.loc[:,'boot'] = pa.Series()

        # # Opboot { itf } bs
        # self.Opboot = ODNested()
        # self.OpAvg  = ODNested()
        # self.OpStd  = ODNested()
        #
        # # Opboot { itf } bs
        # self.OpAuto = ODNested()
        # self.OpAutoAvg  = ODNested()
        # self.OpAutoStd  = ODNested()

        ## will be of type SetOfFitFuns from SetOfFits.py
        # self.flowFit_Stats = None
        # '''
        # self.Op_Stats DataFrame:
        #
        # columns =   boot, Avg,      Std
        #
        # rows =      fit range
        # '''
        self.flowFit_Stats = None

        # # FlowFitDict [ Fit_range ,  Fit_Parameter ] bs
        # self.flowFitDict    = ODNested()
        # self.flowFitDictAvg = ODNested()
        # self.flowFitDictStd = ODNested()
        self.latparams = mp.LatticeParameters(Info=Info)
        self.latparams.LoadPickle()

        self.nt = self.latparams.nt
        self.nxyz = self.latparams.nxyz
        self.dim_label = 'RC'+str(self.nxyz)+'x'+str(self.nt)
        self.dim_ll = str(self.nxyz)+r'^{3}\times ' + str(self.nt)
        self.nboot = thisnboot

        if 'Do_Chi' in list(Info.keys()): self.do_chi = Info['Do_Chi']
        else: self.do_chi = False

        if 'Sparam' in list(Info.keys()): self.Sparam = Info['Sparam']
        else: self.Sparam = defSparam

        if 'show_cfgs' in list(Info.keys()): self.show_cfgs = Info['show_cfgs']
        else: self.show_cfgs = False


        ## Info can contain fit ranges to be calculated
        ## must bein form Info['FlowFits'] = [ ('state#' , 'fit#-#' ) , ...]
        if 'FlowFits' in list(Info.keys()): self.flowfit_range = Info['FlowFits']
        else: self.flowfit_range = 'None'


        ## Info can contain initial guess for fits, make sure it array of correct length
        if 'FlowiGuess' in list(Info.keys()): self.flowiGuess = np.array(Info['FlowiGuess'])
        else: self.flowiGuess = 'None' ## Defaults to value set for function in FitFunctions.py

        if 'cfg_file_type' in list(Info.keys()): self.cfg_file_type = Info['cfg_file_type']
        else: self.cfg_file_type = cfg_file_type


        ## Info can contain fit ranges to be calculated
        ## must bein form Info['FlowFitFun'] = (Function Object , number of parameters )
        if 'FlowFitFun' in list(Info.keys()):
            thisFun,thisnpar = Info['FlowFitFun']
            self.FlowSetFunction(thisFun,thisnpar)
        else:
            self.FlowSetFunction('None',0)


        ## kappa up down (in integer form)
        if 'kud' in list(Info.keys()): self.kud = 'kud'+str(Info['kud'])
        else: self.kud = ''

        ## kappa strange (in integer form)
        if 'ks' in list(Info.keys()): self.ks = 'ks'+str(Info['ks'])
        else: self.ks = ''

        ## Current obersvables are TopCharge and Weinberg
        if 'Observable' in list(Info.keys()): self.Observ = Info['Observable']
        else: self.Observ = 'TopCharge' ## default to TopCharge
        self.is_comb = False

        ## this is used display only these flowtimes, Still reads and writes all
        if 'tflowlist' in list(Info.keys()): self.tflowlist = Info['tflowlist']
        else: self.tflowlist = defxlimOp ## default to read all flow times


        ## this is used display only these flowtimes, Still reads and writes all
        self.FlowNewForm = flow_new_format
        # if 'FlowNewForm' in Info.keys(): self.FlowNewForm = Info['FlowNewForm']
        # else: self.FlowNewForm = False ## default to read all flow times

        ## will attempt to find all replica streams in directory.
        ## set self.nstream to number of streams you want to find.
        if 'stream_list' in list(Info.keys()):
            if isinstance(Info['stream_list'],str):
                if Info['stream_list'].lower() != 'all':
                    self.stream_list = [Info['stream_list']]
                else:
                    self.stream_list = 'all'
            else:
                self.stream_list = [iinfo.lower() for iinfo in Info['stream_list']]
        else:
            self.stream_list = 'all'


        listout = []
        for ict,itflow in enumerate(self.tflowlist):
            if 't_f' not in str(itflow):
                listout.append(tflowstr(itflow))
            else:
                listout.append(itflow)
        self.tflowfloat = [untflowstr(itflow) for itflow in listout]
        self.tflowlist = np.array(listout)

        thisfolder,self.flowfilepref = GetFolderFileOp(self.Observ)
        self.flowfilename = [self.flowfilepref+'_flow_b1.90_ng','.out']

        kappafolder = 'Kud0'+self.kud.replace('kud','')+'Ks0'+self.ks.replace('ks','')
        # self.flowfilename = ['RC32x64_B1900'+kappafolder+'C1715','_','_k'+self.kud.replace('kud','')+'_'+self.tsrc+filejsm+fileism+'_nucleon.2cf.xml']
        self.flowInfo = Info
        self.flowfiledir = datadir+thisfolder+'/'+self.dim_label+'_'+kappafolder+'-*-/PerGF/'
        self.FlowSetCustomName(name,LLname)
        if not man_load_cfgs:
            self.LoadCfgs(cfglist,name,LLname)

    def LoadCfgs(self,cfglist,name='',LLname=''):
        self.Create_Streams_Dirs()

        ## cfglist { configuration : [ x-source numbers ] }
        if isinstance(cfglist,(list,np.ndarray,pa.DataFrame)) and len(cfglist) > 0:
            self.FlowImportCfgList(cfglist)
        else:
            self.FlowGetCfgList()

        self.FlowSetCustomName(name,LLname)

    def Create_Streams_Dirs(self):
        if isinstance(self.stream_list,str) and self.stream_list == 'all':
            self.All_Streams()
        else:
            self.Check_Streams()
        self.stream_name = '-'+''.join([istr.replace('-','') for istr in self.stream_list])+'-'
        self.flowdir_list = {}
        for istream in self.stream_list:
            thisdir = self.flowfiledir.replace('-*-',istream)
            if os.path.isdir(thisdir):
                self.flowdir_list[istream] = thisdir
            else:
                self.flowdir_list[istream] = thisdir.replace(self.dim_label+'_','')




    def All_Streams(self):
        self.stream_list = []
        for ifolder in glob.glob(self.flowfiledir):
            self.stream_list.append(re.findall('-.*-',ifolder)[-1])
        for ifolder in glob.glob(self.flowfiledir.replace(self.dim_label+'_','')):
            this_stream = re.findall('-.*-',ifolder)[-1]
            if this_stream not in self.stream_list:
                self.stream_list.append(this_stream)

        if len(self.stream_list) == 0:
            print(self.flowfiledir+'\n No Streams found')

    def Check_Streams(self):
        for istream in self.stream_list:
            this_check = os.path.isdir(self.flowfiledir.replace('-*-',istream))
            this_check = this_check or os.path.isdir(self.flowfiledir.replace('-*-',istream).replace(self.dim_label+'_',''))
            if not this_check:
                print((self.flowfiledir.replace('-*-',istream) + '\n stream directory not found, as well as \n'+
                        self.flowfiledir.replace('-*-',istream).replace(self.dim_label+'_','')))


    def GetInterp(self,this_to):
        fit_info = {}
        fit_info['Funs'] = (LinearFitFun,2)
        fit_info['name'] = self.flowname + ' t0fit'
        left_t0,right_t0 = 100,100
        for it0 in map(untflowstr,list(self.Op_Stats.index)):
            if it0-this_to < 0 and np.abs(it0-this_to) < np.abs(left_t0-this_to):
                left_t0 = it0
            elif it0-this_to > 0 and np.abs(it0-this_to) < np.abs(right_t0-this_to):
                right_t0 = it0
        # fit_min = list(self.Op_Stats.index).index(tflowstr(left_t0))
        # fit_max = list(self.Op_Stats.index).index(tflowstr(right_t0))
        # fit_loc = fit_min + ((this_to-left_t0)*((right_t0-left_t0)/(fit_max-fit_min)))
        # print('DEBUG',fit_min,fit_max)
        # print(self.Op_Stats['boot'])
        self.t0_fit = sff.SetOfFitFuns(data=self.Op_Stats['boot'])
        this_test = self.t0_fit.ScanRange(left_t0,right_t0,fit_info=fit_info,min_fit_len=0,from_xdata=True)
        if this_test:
            self.t0_fit.DoFits()
            # self.t0_fit.SortChi()
            self.FlowWrite()
        else:
            self.t0_fit = None
            return float('NaN')
        return self.t0_fit.Fit_Stats['Fit'].iloc[0].Eval_Function(this_to)


    def Compute_t0(self,nfit=6,fix_energy=True,DefWipe=False):
        ## THIS FUNCTION ASSUMES DATA IS <E(t)>
        if hasattr(self,'t0_values') and not DefWipe:
            print('t0_values alread computed, skipping')
            return
        fit_info = {}
        fit_info['Funs'] = (LinearFitFun,2)
        fit_info['name'] = self.flowname + ' t0fit'
        if fix_energy:
            self.Op_Stats.loc[:,'t2boot'] = self.Op_Stats.apply(lambda x : (untflowstr(x.name)**2)*6*x['boot']/(self.latparams.nt*self.latparams.nxyz**3),axis=1)
        else:
            self.Op_Stats.loc[:,'t2boot'] = self.Op_Stats.apply(lambda x : (untflowstr(x.name)**2)*x['boot'],axis=1)
        self.Op_Stats.loc[:,'t2boot'].apply(lambda x : x.Stats())
        this_df = deepcopy(self.Op_Stats)
        this_df.loc[:,'t2boot_min_p3'] = np.abs(this_df.loc[:,'t2boot'].apply(lambda x : x.Avg)-0.3)
        this_df = this_df.sort_values(by=['t2boot_min_p3'])
        t0_list = list(map(untflowstr,list(this_df.index[:nfit])))
        t0_min,t0_max = tflowstr(min(t0_list)),tflowstr(max(t0_list))
        fit_min = list(self.Op_Stats.index).index(t0_min)
        fit_max = list(self.Op_Stats.index).index(t0_max)
        # print('DEBUG',fit_min,fit_max)
        # print(self.Op_Stats['boot'])
        self.t0_fit = sff.SetOfFitFuns(data=self.Op_Stats['t2boot'])
        this_test = self.t0_fit.ScanRange(fit_min,fit_max,fit_info=fit_info,min_fit_len=0)
        if this_test:
            self.t0_fit.DoFits()
            # self.t0_fit.SortChi()
            self.FlowWrite()
        else:
            self.t0_fit = None
            return
        t0_list,tf_min_list,tf_max_list,fitr_list,fit_len = [],[],[],[],[]
        for ifit,fit_data in self.t0_fit.Fit_Stats_fmt['Fit'].items():
            linear_term = fit_data.fit_data['Params'].iloc[0]
            const_term = fit_data.fit_data['Params'].iloc[1]
            t0 = (0.3-const_term)/linear_term
            t0.Stats()
            t0_list.append(t0)
            fitr_list.append(ifit[1])
            fit_min,fit_max = unxmlfitr_int(ifit[1])
            tf_fitr_min = list(self.Op_Stats.index)[fit_min]
            tf_fitr_max = list(self.Op_Stats.index)[fit_max]
            tf_min_list.append(tf_fitr_min)
            tf_max_list.append(tf_fitr_max)
            fit_len.append(fit_max - fit_min)
        if len(fit_len)>0:
            self.t0_values = pa.DataFrame()
            self.t0_values['boot'] = pa.Series(t0_list)
            self.t0_values['tf_min'] = pa.Series(tf_min_list)
            self.t0_values['tf_max'] = pa.Series(tf_max_list)
            self.t0_values['fitr'] = pa.Series(fitr_list)
            self.t0_values['fit_len'] = pa.Series(fit_len)
            self.t0_values = self.t0_values.sort_values(by=['fit_len'],ascending=False)
            self.t0_values['sqrt8tf'] = self.t0_values['boot'].apply(lambda x : (8*x).Sqrt())
            self.t0_values['sqrt8tf'].apply(lambda x : x.Stats())
            self.latparams.Set_t0(self.t0_values,self.Op_cfgs.loc[:,self.Op_cfgs.columns != 'Op'])
            self.FlowWrite()

    def FlowGetCfgList(self):
        def SplitCfg(ifile):
            icfg = re.findall('_ng.*',ifile)
            if len(icfg) == 0:
                print()
                print('WARNING: file does not have cfg number')
                print(ifile)
                print()
                return ''
            icfg = icfg[0].replace('_ng','').replace('.out','')
            return icfg

        cfglist,ilist = [],[]
        low_test,upper_test = self.flowfilepref.upper(),self.flowfilepref.lower()
        for istream in self.stream_list:
            idir = self.flowdir_list[istream]
            for icfg,ifile in enumerate(glob.glob(idir+'/'+self.flowfilepref+'_*')):
                if low_test in ifile or upper_test in ifile:
                    ilist.append((istream,icfg))
                    cfglist.append(SplitCfg(ifile))
        if len(cfglist) == 0:
            print('No configs found for directory:')
            print(' , '.join(self.flowdir_list))
        else:
            indicies = pa.MultiIndex.from_tuples(ilist,names=['stream','config_number'])
            cfg_df = pa.DataFrame(cfglist,columns=['configs'],index=indicies)
            cfg_df.loc[:,'int_configs'] = pa.Series(list(map(int,cfglist)),index=indicies)
            cfg_df = cfg_df.sort_values('int_configs').sort_index(level='stream',sort_remaining=False)
            del cfg_df['int_configs']
            self.FlowImportCfgList(cfg_df)


    # def MakeOutputFile(self):
    #     self.name = '_'.join([self.tsrc,self.ism,self.jsm])
    #     self.HumanFile = '/'+kappafolder+'/G2/'+self.name+'.xml'
    #     self.PickleFile = '/'+kappafolder+'/G2/Pickle/'+self.name+'.py3p'


    ## Routines for Plotting ###############################################

    def FlowCheckCol(self,thiscol):
        # if 'PreDefine' == thiscol :
            # if 'Not Set' == self.thiscol:
                # raise IOError('Pass in color to initialize it')
        # else:
        if 'PreDefine' not in thiscol: self.thiscol = thiscol

    def FlowCheckSym(self,thissym):
        # if 'PreDefine' == thissym :
            # if 'Not Set' == self.thissym:
                # raise IOError('Pass in symbol to initialize it')
        # else:
        if 'PreDefine' not in thissym: self.thissym = thissym


    def FlowGetShift(self,xlims,thisshift):
        if 'PreDefine' not in thisshift: self.thisshift = thisshift
        xlen = np.abs(xlims[1]-xlims[0])
        return self.thisshift*xlen


    def PlotMonteTime(self,plot_class,flowtime,stream='first',thiscol='PreDefine',thisshift='PreDefine'):
        if 'Op' not in self.Op_cfgs:
            # print 'Warning: no configurations in memory for '+self.flowname + ' attempting to read again'
            # self.FlowRead()
            print('Warning: no configurations in memory for '+self.flowname + ' not doing monteplot')
            return plot_class
        if stream == 'first':
            stream = self.stream_list[0]
        if thiscol == 'PreDefine': thiscol = self.thiscol
        if isinstance(self.tflowlist,np.ndarray): self.tflowlist = self.tflowlist.tolist()
        hold_series = null_series
        hold_series['x_data'] = 'from_keys'
        vall,indexl = [],[]
        for ((istream,dump),cfg_vals),icfg in zip(iter(self.Op_cfgs['Op'].items()),self.Op_cfgs['configs'].values):
            for itflow,flow_val in zip(self.tflowlist,cfg_vals):
                vall.append(flow_val)
                indexl.append((StreamToNumb(istream),icfg,itflow))
        if len(indexl) == 0:
            print('No data to plot for montetime for ',self.name)
            return plot_class
        indicies = pa.MultiIndex.from_tuples(indexl,names=['Stream','Config Number','Flow Time'])
        ploty = pa.Series(vall,index=indicies)
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = ploty
        if thisshift == 'PreDefine': thisshift = self.thisshift
        if isinstance(thisshift,str): thisshift = 0.0
        hold_series['key_select'] = (StreamToNumb(stream),slice(None),flowtime)
        # thisshift = np.abs(hold_series['x_data'][0]-hold_series['x_data'][-1])*thisshift
        hold_series['shift'] = thisshift
        hold_series['label'] = self.flowLegLab
        if 'Not Set' not in thiscol: hold_series['color'] = thiscol
        hold_series['type'] = 'plot_vary'
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class
        # pl.plot(plotx+thisshift,ploty,label=,color=thiscol)
        # pl.scatter(plotx+thisshift,ploty,label=self.flowLegLab+' $'+tflowTOLeg(flowtime)+'$',color=self.thiscol)



    ## MultQ means Q * deltaQ<eps
    def PlotMonteTimeDeltaEps(self,plot_class,flowtime,stream='first',thiseps=Qeps,MultQ=False,thiscol='PreDefine',thisshift='PreDefine'):
        if 'Op' not in self.Op_cfgs:
            # print 'Warning: no configurations in memory for '+self.flowname + ' attempting to read again'
            # self.FlowRead()
            print('Warning: no configurations in memory for '+self.flowname + ' not doing monteplot')
            return plot_class
        if stream == 'first':
            stream = self.stream_list[0]
        if isinstance(self.tflowlist,np.ndarray): self.tflowlist = self.tflowlist.tolist()
        if thiscol == 'PreDefine': thiscol = self.thiscol
        thisindex = self.tflowlist.index(flowtime)

        hold_series = null_series
        hold_series['x_data'] = self.Op_cfgs['configs'][stream].astype(int).values
        if MultQ:
            hold_series['y_data'] = [iOp[thisindex]*int(np.abs(iOp[thisindex])<thiseps) for iOp in self.Op_cfgs['Op'][stream].values]
        else:
            hold_series['y_data'] = [int(np.abs(iOp[thisindex])<thiseps) for iOp in self.Op_cfgs['Op'][stream].values]
        if thisshift == 'PreDefine': thisshift = self.thisshift
        if isinstance(thisshift,str): thisshift = 0.0
        thisshift = np.abs(hold_series['x_data'][0]-hold_series['x_data'][-1])*thisshift
        hold_series['key_select'] = None
        hold_series['shift'] = thisshift
        hold_series['label'] = self.flowLegLab.replace(self.stream_name,stream)+' $'+tflowTOLeg(flowtime)+'$ '
        if 'Not Set' not in thiscol: hold_series['color'] = thiscol
        hold_series['type'] = 'plot'
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class

    def PlotTauIntWOpt(self,plot_class,xlims='data',thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',FlowRad=True):
        if xlims == 'data':
            xlims = list(map(untflowstr,self.tflowlist))
        self.FlowCheckCol(thiscol)
        self.FlowCheckSym(thissym)
        # if not hasattr(self,'OpAuto'):
        #     print 'Warning, Depreciated result has no autocorrelation results, skipping TauIntWOpt plot'
        #     return
        dataavg,dataerr,tlistplot = [],[],[]
        for it,idata in self.Op_Stats['Auto'].items():
            if str(float(untflowstr(it))) in list(map(str,xlims)):
                tlistplot.append(it)
                dataavg.append(idata.tau)
                dataerr.append(idata.tauerr)
        if len(tlistplot) == 0:
            print('Warning, no flow times found for plot')
            return
        tflowphys = np.array(list(map(untflowstr,tlistplot)))
        tflowphys = TflowToPhys(np.array(tflowphys),self.latparams.latspace)
        thisshift = self.FlowGetShift([tflowphys[0],tflowphys[-1]],thisshift)
        hold_series = null_series
        hold_series['x_data'] = tflowphys
        hold_series['y_data'] = dataavg
        hold_series['yerr_data'] = dataerr
        hold_series['key_select'] = None
        hold_series['shift'] = thisshift
        hold_series['label'] = self.flowLegLab
        if 'Not Set' not in self.thiscol: hold_series['color'] = self.thiscol
        if 'Not Set' not in self.thissym: hold_series['symbol'] = self.thissym
        hold_series['type'] = 'error_bar'
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class


        # pl.errorbar(tflowphys+thisshift,dataavg,dataerr,label=,fmt=,color=self.thiscol)
        # pl.xlim(xlims[0],xlims[-1])




    ## make sure all xlims are the same, because this sets the window of xlim (needed for calculating shifts)
    def FlowPlotAuto(self,plot_class,xlims='data',thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',
                     FlowRad=True,t0_xvals=True):
        if xlims == 'data':
            xlims = list(map(untflowstr,self.tflowlist))
        self.FlowCheckCol(thiscol)
        self.FlowCheckSym(thissym)
        # if not hasattr(self,'OpAuto'):
        #     print 'Warning, Depreciated result has no autocorrelation results, skipping Auto plot'
        #     return
        dataavg,dataerr,tlistplot = [],[],[]
        for it,idata in self.Op_Stats['Auto'].items():
            if str(float(untflowstr(it))) in list(map(str,xlims)):
                tlistplot.append(it)
                dataavg.append(idata.Avg)
                dataerr.append(idata.Std)
        if len(tlistplot) == 0:
            print('Warning, no flow times found for plot')
            return plot_class
        if t0_xvals:
            try:
                this_s8t0 = (8*self.latparams.Get_t0()['boot'].iloc[0]).Sqrt()
            except Exception as err:
                print('t0_xvals error:',err)
                t0_xvals = False
        if t0_xvals:
            tflowphys,tflowphys_err = [],[]
            for itflow in map(untflowstr,tlistplot):
                it = np.sqrt(8*itflow)/this_s8t0
                it.Stats()
                tflowphys.append(it.Avg)
                tflowphys_err.append(it.Std)
            tflowphys = np.array(tflowphys)
            tflowphys_err = np.array(tflowphys_err)
        if not t0_xvals:
            tflowphys = np.array(list(map(untflowstr,tlistplot)))
            tflowphys = TflowToPhys(np.array(tflowphys),self.latparams.latspace)
        # tflowphys = np.array(list(map(untflowstr,tlistplot)))
        # tflowphys = TflowToPhys(np.array(tflowphys),self.latparams.latspace)
        thisshift = self.FlowGetShift([tflowphys[0],tflowphys[-1]],thisshift)
        hold_series = null_series
        hold_series['x_data'] = tflowphys
        if t0_xvals:
            hold_series['xerr_data'] = tflowphys_err
        hold_series['y_data'] = dataavg
        hold_series['yerr_data'] = dataerr
        hold_series['key_select'] = None
        hold_series['fit_class'] = None
        hold_series['shift'] = thisshift
        hold_series['label'] = self.flowLegLab+' Auto'
        if 'Not Set' not in self.thiscol: hold_series['color'] = self.thiscol
        if 'Not Set' not in self.thissym: hold_series['symbol'] = self.thissym
        hold_series['type'] = 'error_bar'
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class
        # pl.errorbar(tflowphys+thisshift,dataavg,dataerr,label=self.flowLegLab+' Auto',fmt=self.thissym,color=self.thiscol)
        # pl.xlim(xlims[0],xlims[-1])

    def FlowPlot_mul_tf2(self,plot_class,xlims='data',thiscol='PreDefine',
                         thissym='PreDefine',thisshift='PreDefine',FlowRad=True,BorA='boot',
                         mul=False,tpow=1,t0_xvals=True,just_E=False,lab_append=''):
        if xlims == 'data':
            xlims = list(map(untflowstr,self.tflowlist))
        self.FlowCheckCol(thiscol)
        self.FlowCheckSym(thissym)
        dataavg,dataerr,tlistplot = [],[],[]
        if BorA not in self.Op_Stats.columns:
            raise IOError(BorA+' \nNot Present in Op_Stats: \n'+','.join(self.Op_Stats.columns))
        for it,idata in self.Op_Stats[BorA].items():
            if str(float(untflowstr(it))) in list(map(str,xlims)):
                tlistplot.append(it)
                if just_E:
                    plot_val =  (idata * 6*float(untflowstr(it))**2)/(self.latparams.nt*self.latparams.nxyz**3)
                elif mul:
                    plot_val =  idata * float(untflowstr(it))**tpow
                    # plot_val = plot_val * self.latparams.nxyz**3
                else:
                    plot_val =  idata*(float(untflowstr(it))**tpow / (6*float(untflowstr(it))**2))
                    plot_val = plot_val *(self.latparams.nt*self.latparams.nxyz**3)
                plot_val.Stats()
                dataavg.append(plot_val.Avg)
                dataerr.append(plot_val.Std)
        if len(tlistplot) == 0:
            print('Warning, no flow times found for plot')
            return plot_class
        if t0_xvals:
            try:
                this_s8t0 = (8*self.latparams.Get_t0()['boot'].iloc[0]).Sqrt()
            except Exception as err:
                print('t0_xvals error:',err)
                t0_xvals = False
        if t0_xvals:
            tflowphys,tflowphys_err = [],[]
            for itflow in map(untflowstr,tlistplot):
                it = np.sqrt(8*itflow)/this_s8t0
                it.Stats()
                tflowphys.append(it.Avg)
                tflowphys_err.append(it.Std)
            tflowphys = np.array(tflowphys)
            tflowphys_err = np.array(tflowphys_err)
        else:
            tflowphys = np.array(list(map(untflowstr,tlistplot)))
            tflowphys = TflowToPhys(np.array(tflowphys),self.latparams.latspace)
        thisshift = self.FlowGetShift([tflowphys[0],tflowphys[-1]],thisshift)
        hold_series = null_series
        hold_series['x_data'] = tflowphys
        if t0_xvals:
            hold_series['xerr_data'] = tflowphys_err
        hold_series['y_data'] = dataavg
        hold_series['yerr_data'] = dataerr
        hold_series['shift'] = thisshift
        hold_series['flowphys'] = self.latparams.latspace
        hold_series['label'] = self.flowLegLab + ' '+lab_append
        if 'Auto' in BorA:
            hold_series['label'] += ' Auto'
        if 'Not Set' not in self.thiscol: hold_series['color'] = self.thiscol
        if 'Not Set' not in self.thissym: hold_series['symbol'] = self.thissym
        hold_series['type'] = 'error_bar'
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class

    ## make sure all xlims are the same, because this sets the window of xlim (needed for calculating shifts)
    def FlowPlot(self,plot_class,xlims='data',thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',
                 FlowRad=True,t0_xvals=True,lab_append=''):
        if xlims == 'data':
            xlims = list(map(untflowstr,self.tflowlist))
        self.FlowCheckCol(thiscol)
        self.FlowCheckSym(thissym)
        dataavg,dataerr,tlistplot = [],[],[]
        for it,idata in self.Op_Stats['boot'].items():
            if str(float(untflowstr(it))) in list(map(str,xlims)):
                tlistplot.append(it)
                dataavg.append(idata.Avg)
                dataerr.append(idata.Std)
        if len(tlistplot) == 0:
            print('Warning, no flow times found for plot')
            return plot_class
        tflowphys = np.array(list(map(untflowstr,tlistplot)))
        if t0_xvals:
            try:
                this_s8t0 = (8*self.latparams.Get_t0()['boot'].iloc[0]).Sqrt()
            except Exception as err:
                print('t0_xvals error:',err)
                t0_xvals = False
        if t0_xvals:
            tflowphys,tflowphys_err = [],[]
            for itflow in map(untflowstr,tlistplot):
                it = np.sqrt(8*itflow)/this_s8t0
                it.Stats()
                tflowphys.append(it.Avg)
                tflowphys_err.append(it.Std)
            tflowphys = np.array(tflowphys)
            tflowphys_err = np.array(tflowphys_err)
        else:
            tflowphys = np.array(list(map(untflowstr,tlistplot)))
            tflowphys = TflowToPhys(np.array(tflowphys),self.latparams.latspace)
        # tflowphys = np.array(list(map(untflowstr,tlistplot)))
        # tflowphys = TflowToPhys(np.array(tflowphys),self.latparams.latspace)
        thisshift = self.FlowGetShift([tflowphys[0],tflowphys[-1]],thisshift)
        hold_series = null_series
        hold_series['x_data'] = tflowphys
        if t0_xvals:
            hold_series['xerr_data'] = tflowphys_err
        hold_series['y_data'] = dataavg
        hold_series['yerr_data'] = dataerr
        hold_series['key_select'] = None
        hold_series['fit_class'] = None
        hold_series['shift'] = thisshift
        hold_series['flowphys'] = self.latparams.latspace
        hold_series['label'] = self.flowLegLab + ' '+ lab_append
        if 'Not Set' not in self.thiscol: hold_series['color'] = self.thiscol
        if 'Not Set' not in self.thissym: hold_series['symbol'] = self.thissym
        hold_series['type'] = 'error_bar'
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class
        # pl.errorbar(tflowphys+thisshift,dataavg,dataerr,label=self.flowLegLab,fmt=self.thissym,color=self.thiscol)
        # pl.xlim(xlims[0],xlims[-1])


    ## state is 1 or 2 state fit
    ## fitr is fit range e.g. fitr5-10
    ## fit range of 'Data' is just the original fit range
    def FlowFitPlot(self,plot_class,fitr,xlims = 'Data',thiscol='PreDefine',thisshift='PreDefine',legparindex=0):
        if fitr == 'None': return plot_class
        self.FlowCheckCol(thiscol)
        thisshift = self.FlowGetShift(TflowToPhys(plot_class.get_xlim(),self.latparams.latspace),thisshift)
        if self.flowFit_Stats is None:
            print(fitr , ' has not been done, performing fit now')
            self.FlowFit(fitr)
        if self.flowFit_Stats is not None:
            hold_series = null_series
            hold_series['type'] = 'fit_vary'
            hold_series['fit_class'] = self.flowFit_Stats
            hold_series['label'] = self.flowFit_Stats[fitr].name
            hold_series['key_select'] = 'First'
            hold_series['color'] = self.thiscol
            hold_series['shift'] = thisshift
            hold_series['xdatarange'] = xlims
            hold_series['flowphys'] = self.latparams.latspace
            hold_series['ShowPar'] = self.thisparam[0]
            hold_series['fmt_class'] = KeyForamtting(self.latparams)
            plot_class.AppendData(hold_series)
        return plot_class


    def FlowObsToName(self,obs):
        if obs == 'TopChargexTopCharge':
            return '<Q^{2}>'
        elif obs == 'WeinbergxWeinberg':
            return '<W^{2}>'
        elif 'TopCharge' in obs and 'Weinberg' in obs:
            return '<WQ>'
        elif 'TopCharge' == obs:
            return '<Q>'
        elif 'Weinberg' == obs:
            return '<W>'
        else:
            return obs

    def FlowSetCustomName(self,string='',stringLL='',RTDir=''):
        if not hasattr(self,'stream_name') or self.stream_name == '--':
            self.stream_name = '-def-'
        if string == '':
            self.flowname = '_'.join([self.dim_label,self.kud,self.ks,self.stream_name,self.Observ])
        else:
            self.flowname = string
        if stringLL == '':
            self.flowLegLab = '$'+'\ '.join([self.dim_ll,self.latparams.GetPionMassLab(),self.stream_name,self.FlowObsToName(self.Observ)]) + '$'
        else:
            # print 'setting Legend to ',stringLL
            self.flowLegLab = stringLL
        if self.FlowNewForm:
            self.OpDir = outputdir +'/FlowOps_NewForm/' + RTDir + '/'
            self.OpCfgsDir =  cfgfmtdir +'/FlowOps_NewForm/' + RTDir + '/'
        else:
            self.OpDir = outputdir +'/FlowOps/' + RTDir + '/'
            self.OpCfgsDir =  cfgfmtdir +'/FlowOps/' + RTDir + '/'
        mkdir_p(self.OpCfgsDir)
        mkdir_p(self.OpDir+'/Pickle/')
        mkdir_p(self.OpDir+'/Excel/')
        self.flowHumanFile = self.OpDir+self.flowname+'.xml'
        self.flowExcelFile = self.OpDir+'/Excel/'+self.flowname+'.xlsx'
        self.flowPickleFile = self.OpDir+'/Pickle/'+self.flowname+'.py3p'
        self.flowCfgFile = self.OpCfgsDir+'/'+self.flowname+'.cfgs'



    def Get_Flattened_Cfgs(self):
        def index_fun(val):
            return list(self.stream_list).index(val)

        col_list = list(self.Op_cfgs.index.names) + self.Op_cfgs.columns.tolist() + ['flow_time']
        red_col_list = self.Op_cfgs.columns[self.Op_cfgs.columns != 'Op'].tolist()
        new_df = pa.DataFrame(columns = col_list)
        new_df.reset_index(inplace=True)
        del new_df['index']
        for ikey in self.Op_cfgs.index:
            for itf,tfdata in enumerate(self.Op_cfgs.loc[ikey,'Op']):
                this_tf = untflowstr(self.tflowlist[itf])
                new_df = new_df.append(pa.Series(list(ikey) + self.Op_cfgs.loc[ikey,red_col_list].tolist()+[tfdata,this_tf],index=col_list),ignore_index=True)
        new_df['stream'] = new_df['stream'].apply(index_fun)
        del new_df['stream_cfgs']
        del new_df['config_number']
        new_df['stream'].astype(np.int)
        new_df.loc[:,'configs'] = new_df['configs'].apply(lambda x : np.int(x)).astype(np.int)
        new_df['Op'].astype(np.float64)
        new_df['flow_time'].astype(np.float16)
        return new_df,list(self.stream_list)


    def Format_Flattened_Cfgs(self,this_df):
        indexl = []
        data_dict = OrderedDict()
        for icol in ['configs','stream_cfgs','Op','flow_time']:
            data_dict[icol] = []
        for ic,(ikey,ival) in enumerate(this_df.iterrows()):
            this_stream_cfg = self.stream_list[ival['stream']]+str(int(ival['configs'])).zfill(6)
            if this_stream_cfg not in data_dict['stream_cfgs']:
                indexl.append((self.stream_list[ival['stream']],ic))
                data_dict['configs'].append(str(int(ival['configs'])).zfill(4))
                data_dict['stream_cfgs'].append(this_stream_cfg)
                data_dict['Op'].append([])
            if tflowstr(ival['flow_time']) not in data_dict['flow_time']:
                data_dict['flow_time'].append(tflowstr(ival['flow_time']))
            data_dict['Op'][-1].append(ival['Op'])
        if len(indexl) > 0 and set(self.tflowlist) <= set(data_dict['flow_time']):
            self.tflowlist = data_dict['flow_time']
            del data_dict['flow_time']
            indicies = pa.MultiIndex.from_tuples(indexl,names=['stream','config_number'])
            new_df = pa.DataFrame(data_dict,index=indicies)
        else:
            new_df = None
        return new_df


    def Write_Cfgs(self,file_type='PreDef',show_timer=True,CheckCfgs=True):
        if hasattr(self,'read_cfgs_from_file') and self.read_cfgs_from_file:
            return
        if file_type == 'PreDef':
            file_type = self.cfg_file_type
        if 'to_' not in file_type:
            file_type = 'to_'+file_type
        this_file = self.flowCfgFile+fmt_file_type(file_type.replace('to_',''))
        if show_timer: thistimer = Timer(name='Writing Cfgs '+this_file)
        if len(self.Op_cfgs.index) == 0:
            self.FlowRead(show_timer=show_timer,CheckCfgs=CheckCfgs)
        # this_df,dump = self.Get_Flattened_Cfgs()
        WriteWithMeta(self.Op_cfgs,np.array(self.tflowlist),this_file,file_type=file_type)
        if show_timer: thistimer.Stop()

    def ChopDepVarList(self,this_df,read_meta):
        isin = set(self.tflowlist).issubset(read_meta)
        if isin:
            same_list = len(self.tflowlist) == len(read_meta)
            if not same_list:
                index_list = [np.where(read_meta == it)[0][0] for it in self.tflowlist]
                def Chop_fun(val):
                    return np.array(val)[index_list].tolist()
                this_df.loc[:,'Op'] = this_df.loc[:,'Op'].apply(Chop_fun)
        return this_df,isin

    def Read_Cfgs(self,file_type='PreDef',show_timer=True,CheckCfgs=True):
        def CheckCfgList(cfg_data,read_tflowlist):
            if cfg_data is None or len(cfg_data) == 0:
                return None,False
            if CheckCfgs and 'stream_cfgs' in self.Op_cfgs:
                second_class = pa.Series([cfg_data],index=['Op_cfgs'])
                cfg_bool = CheckClass(self.__dict__,second_class,['Op_cfgs'])
                # cfg_bool = np.array(self.Op_cfgs['stream_cfgs'].values) == np.array(cfg_data['stream_cfgs'].values)
                # if isinstance(cfg_bool,np.ndarray):
                #     cfg_bool = all(cfg_bool)
                # if show_timer and not cfg_bool:
                #     print '\n'+GetCfgMissmatch(  np.array(self.Op_cfgs['stream_cfgs'].values),
                #                             np.array(cfg_data['stream_cfgs'].values))
            else:
                cfg_bool = True
            if cfg_bool:
                cfg_data,outbool = self.ChopDepVarList(cfg_data,read_tflowlist)
                return cfg_data,outbool
            else:
                return None,False

        if file_type == 'PreDef':
            file_type = self.cfg_file_type
        if 'read_' not in file_type:
            file_type = 'read_'+file_type
        this_file = self.flowCfgFile+fmt_file_type(file_type.replace('read_',''))
        if not os.path.isfile(this_file):
            if '-def-' in this_file:
                this_file = this_file.replace('-def-','-*-')
                file_list = glob.glob(this_file)
                if len(file_list) == 0:
                    print('Cfg FNF: \n', this_file)
                    self.read_cfgs_from_file = False
                    return False
                else:
                    this_file = file_list[0]
        if not os.path.isfile(this_file):
            self.read_cfgs_from_file = False
            if show_timer: print('cfg FNF: ' + this_file)
            return False
        if show_timer: thistimer = Timer(name='Reading Cfgs '+this_file)
        meta_data,cfg_data = ReadWithMeta(this_file,file_type=file_type)
        # cfg_data = self.Format_Flattened_Cfgs(cfg_data)
        cfg_data = self.MakeStreamCfgs(cfg_data)
        cfg_data,check_bool = CheckCfgList(cfg_data,meta_data)
        self.read_cfgs_from_file = check_bool
        if check_bool:
            self.FlowImportCfgList(cfg_data)
            if show_timer: thistimer = thistimer.Stop()
            return True
        else:
            if show_timer: print('\n file incorrect,', end=' ')
            self.flowCfgFile = self.flowCfgFile.replace('.cfgs','.cfgs.bak')
            if os.path.isfile(this_file.replace('.cfgs','.cfgs.bak')):
                print(' found backup file')
                return self.Read_Cfgs(file_type=file_type,show_timer=show_timer,CheckCfgs=CheckCfgs)
            else:
                print(' reading unformatted')
                return False

        # return cfg_data is not None and all(self.Op_cfgs['stream_cfgs'].values == cfg_data['stream_cfgs'].values)
    def MakeStreamCfgs(self,cfglist):
        if 'stream_cfgs' not in cfglist.columns and 'configs' in cfglist:
            cfglist.loc[:,'stream_cfgs'] = pa.Series(['00'.join((istream,icfg)) for (istream,iccfg),icfg in cfglist['configs'].items()],index=cfglist.index)
        return cfglist

    def FlowImportCfgList(self,cfglist,stream_list='PreDefine'):
        # if isinstance(cfglist,(dict,OrderedDict,pa.Series)): cfglist = cfglist.keys()
        if not isinstance(cfglist,pa.DataFrame):
            if not isinstance(stream_list,str):
                self.stream_list = stream_list
            ## cfglist must be 2-d array, first dimension is stream, second dimension is configuraiton number
            ilist = []
            for istream,is_cfg in enumerate(cfglist):
                for icf in range(is_cfg):
                    ilist.append((istream,is_cfg))
            indicies = pa.MultiIndex.from_tuples(ilist,names=['stream','config_number'])
            cfglist = pa.DataFrame(np.array(cfglist).flatten(),columns=['configs'],index=indicies)
        self.Op_cfgs = self.MakeStreamCfgs(cfglist)
        self.Update_Stream_List()
        if 'configs' in self.Op_cfgs:
            self.flowncfg_list = [icfg.size for istream,icfg in self.Op_cfgs['configs'].groupby(level='stream')]

    def Update_Stream_List(self):
        if 'Op' in self.Op_cfgs:
            self.stream_list = self.Op_cfgs.index.levels[0]
            self.nstream = len(self.stream_list)
            old_sn = self.stream_name
            self.stream_name = '-'+''.join([istr.replace('-','') for istr in self.stream_list])+'-'
            if hasattr(self,'flowname'):
                self.flowname = self.flowname.replace(old_sn,self.stream_name)
            if hasattr(self,'flowLegLab'):
                self.flowLegLab = self.flowLegLab.replace(old_sn,self.stream_name)

    def FlowStats(self):
        if 'boot' not in self.Op_Stats: return
        Avglist,Stdlist = [],[]
        for iOp in self.Op_Stats['boot']:
            iOp.Stats()
            Avglist.append(iOp.Avg)
            Stdlist.append(iOp.Std)
        self.Op_Stats.loc[:,'Avg'] = pa.Series(Avglist,index=self.Op_Stats.index)
        self.Op_Stats.loc[:,'Std'] = pa.Series(Stdlist,index=self.Op_Stats.index)


    def FlowRemoveVals(self):
        if 'Op' in self.Op_cfgs: del self.Op_cfgs['Op']


    def FlowRead(self,show_timer=True,cfg_file_type= 'PreDef',CheckCfgs=True):
        if self.is_comb:
            return

        def FlowReadTopCharge(thisfile):
            if os.path.isfile(thisfile):
                readfile = thisfile
            elif os.path.isfile(thisfile.replace('ng00','ng')):
                readfile = thisfile.replace('ng00','ng')
            else:
                print('Warning, file not found',thisfile)
                return [],[]

            with open(readfile,'r') as ifile:
                tflowlist,topchargelist = [],[]
                for iline,line in enumerate(ifile):
                    val = list(map(np.float64,line.split()))
                    if len(val) == 3:
                        tflow,topcharge,cmplxdump = val
                    elif len(val) == 2:
                        tflow,topcharge = val
                    else:
                        raise IOError('error with file: \n'+readfile+'\n on line: '+str(iline+1))
                    if tflowstr(tflow) in self.tflowlist:
                        tflowlist.append(np.float(tflow))
                        topchargelist.append(topcharge)
            return tflowlist,topchargelist


        def FlowReadTopCharge_New(thisfile):
            if os.path.isfile(thisfile):
                readfile = thisfile
            elif os.path.isfile(thisfile.replace('ng00','ng')):
                readfile = thisfile.replace('ng00','ng')
            else:
                print('Warning, file not found',thisfile)
                return [],[]
            with open(readfile,'r') as ifile:
                tflowlist,topchargelist = [],[]
                flow_line = True
                for iline,line in enumerate(ifile):
                    if flow_line and len(line.strip()) != 0:

                        val = list(map(np.float64,line.split()))
                        if len(val) == 3:
                            tflow,topcharge,cmplxdump = val
                        elif len(val) == 2:
                            tflow,topcharge = val
                        else:
                            raise IOError('error with file: \n'+readfile+'\n on line: '+str(iline+1))
                        ## this is fix for 10.0 problem
                        if tflowstr(tflow) in self.tflowlist:
                            tflowlist.append(np.float(tflow))
                            topchargelist.append(topcharge)
                        flow_line = False
                    elif len(line.strip()) == 0:
                        flow_line = True
            tflowlist,topchargelist = zip(*[[y,x] for y,x in sorted(zip(tflowlist,topchargelist))])
            return list(tflowlist),list(topchargelist)

        ## making interpolator number more human readable
        # if len(self.flowname.split('_')) > 3:
        #     raise IOError('Do not read after combining correlators, please recreate the correlator')
        # ntflow = 0


        if 'file_names' in self.Op_cfgs.columns:
            if show_timer: thistimer = Timer(linklist=list(self.Op_cfgs['file_names'].values),name='Read '+self.flowname)
            if len(self.Op_cfgs.index) == 0:
                raise IOError('No values found for ' +self.Observ)
            thisdata = []
            for ifile in self.Op_cfgs['file_names'].values:
                # ifile = self.flowdir_list[istream]+icfg.join(self.flowfilename)
                if self.FlowNewForm:
                    tflow,topcharge = FlowReadTopCharge_New(ifile)
                else:
                    tflow,topcharge = FlowReadTopCharge(ifile)
                for itflow in self.tflowfloat:
                    if itflow not in tflow:
                        print()
                        print(self.tflowfloat)
                        print(tflow)
                        print(itflow)
                        print(ifile)
                        raise IOError('flowed operator does not contain all flow times, check flow time list or remove file if corrupt')
            # print
            # print 'Debug'
            # print ifile
            # for iflow,itop in zip(tflow,topcharge):
            #     print iflow,itop
            # if ntflow != len(tflow) and ntflow != 0:
            #     print 'Warning, tflow list has changed for file'
            #     print self.flowfiledir+ifile
            #     print 'base comparision ntflow is ' , ntflow , ' , this file has ' , len(tflow)
            #     continue
            # ntflow = max(len(tflow),ntflow)
                thisdata.append(topcharge)
                if show_timer: thistimer.Lap(ifile)
            self.Op_cfgs.loc[:,'Op'] = pa.Series(thisdata,index=self.Op_cfgs.index)
        else:
            if show_timer: thistimer = Timer(linklist=self.Op_cfgs['configs'].values,name='Read '+self.flowname)
            if len(self.Op_cfgs['configs'].values) == 0:
                raise IOError('No values found for ' +self.Observ)
            thisdata = []
            if self.Read_Cfgs(file_type=cfg_file_type,show_timer=show_timer,CheckCfgs=CheckCfgs):
                return
            else:
                if 'configs' not in self.Op_cfgs:
                    raise EnvironmentError('No config files found anywhere.')
            for (istream,iccfg),icfg in self.Op_cfgs['configs'].items():
                ifile = self.flowdir_list[istream]+icfg.join(self.flowfilename)
                if self.FlowNewForm:
                    tflow,topcharge = FlowReadTopCharge_New(ifile)
                else:
                    tflow,topcharge = FlowReadTopCharge(ifile)
                if any([itflow not in tflow  for itflow in self.tflowfloat]):
                    print()
                    print(self.tflowfloat)
                    print(tflow)
                    print(ifile)
                    raise IOError('flowed operator does not contain all flow times, check flow time list or remove file if corrupt')
            # print
            # print 'Debug'
            # print ifile
            # for iflow,itop in zip(tflow,topcharge):
            #     print iflow,itop
            # if ntflow != len(tflow) and ntflow != 0:
            #     print 'Warning, tflow list has changed for file'
            #     print self.flowfiledir+ifile
            #     print 'base comparision ntflow is ' , ntflow , ' , this file has ' , len(tflow)
            #     continue
            # ntflow = max(len(tflow),ntflow)
                thisdata.append(topcharge)
                if show_timer: thistimer.Lap()
            self.Op_cfgs.loc[:,'Op'] = pa.Series(thisdata,index=self.Op_cfgs.index)
        # if self.flowncfg*len(tflow) != len(thisdata):
        #     print 'ncfg',self.flowncfg,'tflowlen', len(tflow)
        #     print 'thidata',len(thisdata)
        #     raise IOError('read in different ncfg * tflow than in file')
        # thisdata = thisdata.reshape(self.flowncfg,len(tflow))
        # thisdata = np.array(thisdata)
        # self.Op = np.rollaxis(thisdata,0,thisdata.ndim)
            self.Write_Cfgs(show_timer=show_timer,file_type=cfg_file_type,CheckCfgs=CheckCfgs)

    def FlowBootAndWrite(self):
        self.FlowBootstrap()
        self.FlowFit()
        self.FlowWrite()

    def FlowFitAndWrite(self):
        self.FlowFit()
        self.FlowWrite()


    ## ease of use function to create <Q^2> or <W^2>
    def FlowMakeOppSquared(self):
        output = self.FlowCombinePreBS(self,'x')
        # output = output.CombinePreBS(output,'x') ## For Fun
        output.BootAndWrite()
        return output

    def GetCombOpp(self,opp1,opp2):
        if 'TopCharge' in opp1:
            if 'Weinberg' in opp2:
                return 'QW'
            elif 'TopCharge' in opp2:
                return 'Q2'
        elif 'Weinberg' in opp1:
            if 'TopCharge' in opp2:
                return 'QW'
            elif 'Weinberg' in opp2:
                return 'W2'

    def GetChitFuns(self,thiscoeff,thispow):
        def chit(*vals):
            return thiscoeff*(vals[0]**thispow)
        def chitder(*vals):
            return [thiscoeff*thispow*(vals[0]**(thispow-1))]
        return (chit,chitder)






    ## ease of use function to create chit = (hbarc/(latspace*nxyz**(0.75)*nt**(0.25)))
    ## or chitW = (h barc/((latspace*nxyz**(0.75)*nt**(0.25)))**0.5)
    def FlowMakeChit(self,OtherOpp,DefWipe=True,Improve=False,ShowPrint=True,RandT='None'):
        if self.do_chi:
            combfile = self.flowPickleFile.replace(self.Observ,'Chit'+self.GetCombOpp(self.Observ,OtherOpp.Observ))
            combfile = combfile.replace(OtherOpp.Observ,'Chit'+self.GetCombOpp(self.Observ,OtherOpp.Observ))
        else:
            combfile = self.flowPickleFile.replace(self.Observ,self.GetCombOpp(self.Observ,OtherOpp.Observ))
            combfile = combfile.replace(OtherOpp.Observ,self.GetCombOpp(self.Observ,OtherOpp.Observ))
        if Improve: combfile.replace('.py3p','Improv.py3p')
        if os.path.isfile(combfile) and not DefWipe:
            if ShowPrint:
                print('Loading Pickle for ', combfile)
                print()
            output = deepcopy(self)
            # output.flowPickleFile = combfile
            loadeddict = ReadPickleWrap(combfile)
            loadeddict = self.FixRead(loadeddict)
            if self.do_chi:
                loadeddict['flowHumanFile'] = loadeddict['flowHumanFile'].replace(self.Observ,'Chit'+self.GetCombOpp(self.Observ,OtherOpp.Observ))
                loadeddict['flowExcelFile'] = loadeddict['flowExcelFile'].replace(self.Observ,'Chit'+self.GetCombOpp(self.Observ,OtherOpp.Observ))
                loadeddict['flowPickleFile'] = loadeddict['flowPickleFile'].replace(self.Observ,'Chit'+self.GetCombOpp(self.Observ,OtherOpp.Observ))
            else:
                loadeddict['flowHumanFile'] = loadeddict['flowHumanFile'].replace(self.Observ,self.GetCombOpp(self.Observ,OtherOpp.Observ))
                loadeddict['flowExcelFile'] = loadeddict['flowExcelFile'].replace(self.Observ,self.GetCombOpp(self.Observ,OtherOpp.Observ))
                loadeddict['flowPickleFile'] = loadeddict['flowPickleFile'].replace(self.Observ,self.GetCombOpp(self.Observ,OtherOpp.Observ))

            thischeck = self.thischecklist
            if CheckClass(output.__dict__,loadeddict,thischeck):
                output.__dict__.update(loadeddict)

                # output.FlowLoadPickle()
                output.thiscol = self.thiscol
                output.thissym = self.thissym
                output.thisshift = self.thisshift
                # if RandT == 'None':
                #     output.FlowSetCustomName(string=output.flowname+'_RandT'+str(RandT),stringLL=output.flowLegLab+'_RandT'+str(RandT),RTDir='RandT')
                output.FlowStats()
                output.FlowFitAndWrite()
                return output
        if ShowPrint:
            print('Combinging operators ')
            print(self)
            print(OtherOpp)
            print()

        if Improve:
            if 'TopCharge' in self.Observ and 'TopCharge' in OtherOpp.Observ:
                output = self.FlowCombinePreBS(OtherOpp,'leftI')
                output.FlowBootstrap(WipeData=False,Improve=True)
            else:
                raise EnvironmentError('Improved is only currently for Q^2')
        else:
            output = self.FlowCombinePreBS(OtherOpp,'x')
            output.FlowBootstrap(WipeData=False)
        Ophold = deepcopy(output.Op_cfgs['Op'])
        # print 'debug'
        # for ikey,iout in output.Opboot.iteritems():
        #     print ikey, iout.Avg,iout.Std
        if 'TopCharge' in self.Observ:
            if 'Weinberg' in OtherOpp.Observ:
                # output = chitcoeff(self.nranget,self.latparams.nxyz,self.latparams.latspace,'WQ') * output**(1/6.)
                # output.FlowSetCustomName('_'.join([self.kud,self.ks,'ChitQW']),'$'+'\ '.join([self.latparams.GetPionMassLab(),
                #                                                                               r'\frac{1}{V^{1/6}}<QW>^{1/6}'])+'$')
                if self.do_chi:
                    self.chitcoeff = chitcoeff(self.latparams.nt,self.latparams.nxyz,self.latparams.latspace,'WQ')
                    output = self.chitcoeff * output
                    output.FlowSetCustomName(   '_'.join([self.dim_label,self.kud,self.ks,self.stream_name,'ChitQW']),
                                                '$'+'\ '.join([self.latparams.GetPionMassLab(),self.stream_name,
                                                 r'\frac{1}{V}<QW>'])+'$')
                else:
                    self.chitcoeff = 1.0
                    output.FlowSetCustomName(   '_'.join([self.dim_label,self.kud,self.ks,self.stream_name,'QW']),
                                                '$'+'\ '.join([self.latparams.GetPionMassLab(),self.stream_name,
                                                 r'<QW>'])+'$')
                thispow = 1.0

                if ShowPrint: print('The Coefficient of <QW> is', self.chitcoeff, ' giving a result in GeV^(6)')
            elif 'TopCharge' in OtherOpp.Observ:
                if self.do_chi:
                    self.chitcoeff = chitcoeff(self.latparams.nt,self.latparams.nxyz,self.latparams.latspace,'QQ')
                    thispow = 0.25
                    output = self.chitcoeff * output**thispow
                    output.FlowSetCustomName(   '_'.join([self.dim_label,self.kud,self.ks,self.stream_name,'ChitQ2']),
                                                '$'+'\ '.join([self.latparams.GetPionMassLab(),self.stream_name,
                                                 r'\frac{1}{V^{1/4}}<Q^{2}>^{1/4}'])+'$')
                else:
                    self.chitcoeff = 1.0
                    thispow = 1.0
                    output.FlowSetCustomName(   '_'.join([self.dim_label,self.kud,self.ks,self.stream_name,'Q2']),
                                                '$'+'\ '.join([self.latparams.GetPionMassLab(),self.stream_name,
                                                 r'<Q^{2}>'])+'$')
                if Improve: output.flowLegLab += ' Improve'
                if ShowPrint: print('The Coefficient of <QQ> is', self.chitcoeff ,' giving a result in GeV')
                # self.chitcoeff = thispow = 1
                # output.FlowSetCustomName('_'.join([self.kud,self.ks,'ChitQ2']),'$'+'\ '.join([self.latparams.GetPionMassLab(),
                #                                                                               r'<Q^{2}>'])+'$')
        elif 'Weinberg' in self.Observ:
            if 'Weinberg' in OtherOpp.Observ:
                if self.do_chi:
                    self.chitcoeff = chitcoeff(self.latparams.nt,self.latparams.nxyz,self.latparams.latspace,'WW')
                    thispow = 0.125
                    output = self.chitcoeff * output**thispow
                    output.FlowSetCustomName(   '_'.join([self.dim_label,self.kud,self.ks,self.stream_name,'ChitW2']),
                                                '$'+'\ '.join([self.latparams.GetPionMassLab(),self.stream_name,
                                                    r'\frac{1}{V^{1/8}}<W^{2}>^{1/8}'])+'$')
                else:
                    self.chitcoeff = 1.0
                    thispow = 1.0
                    output.FlowSetCustomName(   '_'.join([self.dim_label,self.kud,self.ks,self.stream_name,'W2']),
                                                '$'+'\ '.join([self.latparams.GetPionMassLab(),self.stream_name,
                                                    r'<W^{2}>'])+'$')
                if ShowPrint: print('The Coefficient of <WW> is', self.chitcoeff ,' giving a result in GeV')
            elif 'TopCharge' in OtherOpp.Observ:
                if self.do_chi:
                    self.chitcoeff = chitcoeff(self.latparams.nt,self.latparams.nxyz,self.latparams.latspace,'QW')
                    thispow = 1
                    output = self.chitcoeff * output
                    output.FlowSetCustomName(   '_'.join([self.dim_label,self.kud,self.ks,self.stream_name,'ChitWQ']),
                                                '$'+'\ '.join([self.latparams.GetPionMassLab(),self.stream_name,
                                                              r'\frac{1}{V}<WQ>'])+'$')
                else:
                    self.chitcoeff = 1.0
                    thispow = 1
                    output.FlowSetCustomName(   '_'.join([self.dim_label,self.kud,self.ks,self.stream_name,'WQ']),
                                                '$'+'\ '.join([self.latparams.GetPionMassLab(),self.stream_name,
                                                              r'<WQ>'])+'$')

                if ShowPrint: print('The Coefficient of <WQ> is', self.chitcoeff ,' giving a result in GeV^(6)')
        ## Improve does not work for this yet, so skip.
        if not Improve:
            autol,avgl,stdl = [],[],[]
            keyl = []
            if ShowPrint: thistimer = Timer(linklist=self.tflowfloat,name='Autocorrelating '+self.flowname)
            for ictflow,itflow in enumerate(self.tflowlist):
                # print itflow, hold, np.mean(output.Opboot[itflow].bootvals)**(1/thispow)/self.chitcoeff,output.OpAuto[itflow].Avg, np.mean(output.Opboot[itflow].bootvals)
                tdata = Ophold.apply(lambda x: x[ictflow]).to_frame('Op')
                autol.append(AutoCorrelate(Fun=self.GetChitFuns(self.chitcoeff,thispow),Sparam=self.Sparam,
                                                      name=output.flowname + ' $' +tflowTOLeg(itflow)+'$',data=tdata))
                keyl.append(itflow)
                avgl.append(autol[-1].Avg)
                stdl.append(autol[-1].Std)
                if ShowPrint: thistimer.Lap()
            output.tflowlist = keyl
            index_list = pa.Index(keyl,name='Flow Times')
            output.Op_Stats.loc[:,'Auto'] = pa.Series(autol,index=index_list)
            output.Op_Stats.loc[:,'AutoAvg'] = pa.Series(avgl,index=index_list)
            output.Op_Stats.loc[:,'AutoStd'] = pa.Series(stdl,index=index_list)
        if RandT != 'None':
            output.FlowSetCustomName(string=output.flowname+'_RandT'+str(RandT),stringLL=output.flowLegLab+'_RandT'+str(RandT),RTDir='RandT')
        output.FlowStats()
        output.chitcoeff = self.chitcoeff
        if RandT == 'None': output.FlowFitAndWrite()
        return output


    ## this function is used for creating operators like <Q^2> or <QW>
    ## Process is , create 2class, class1/2.Read() , class_comb = class1.CombinePreBS(class2) , class_comb.Bootstrap() , class_comb.Write() etc..

    def FlowCombinePreBS_FixFT_Interp(self,FlowOps2,Operation='x',flow_time=6.01,Phys=False):
        if Phys:
            flow_time=flow_time/(self.latparams.latspace**2)
        flow_list = list(map(untflowstr,list(self.Op_Stats.index)))
        def SortFun_pos(val):
            val = flow_time-val
            if val > 0:
                return abs(val)
            else:
                return 1000
        def SortFun_neg(val):
            val = flow_time-val
            if val < 0:
                return abs(val)
            else:
                return 1000
        left_flow = min(flow_list,key=lambda x :SortFun_pos(x))
        right_flow = min(flow_list,key=lambda x :SortFun_neg(x))
        left_flow = 't_f'+str(left_flow)
        right_flow = 't_f'+str(right_flow)
        left_data = self.FlowCombinePreBS_FixFT(FlowOps2,Operation=Operation,flow_time=left_flow)
        right_data = self.FlowCombinePreBS_FixFT(FlowOps2,Operation=Operation,flow_time=right_flow)
        return left_data,right_data,left_flow,right_flow

    def FlowCombinePreBS_FixFT(self,FlowOps2,Operation='x',flow_time=0):
        thisOpFun = GetOppFromStr(Operation)
        if not isinstance(FlowOps2,FlowOp): raise IOError('FlowCombinePreBS_FixFT only takes FlowOp types, not '+FlowOps2.__class__.__name__)
        output = FlowOp(thisnboot=self.nboot, cfglist=self.Op_cfgs,Info=CombineInfo(self.flowInfo,FlowOps2.flowInfo,Operation.replace('*','x')),
                        thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        this_Op = []
        this_index = []
        # print 'DEBUG 3, '
        # print self, len(self.Op),
        # print FlowOps2,len(FlowOps2.Op)
        # print
        if isinstance(flow_time,str):
            flow_time = tflowstr(untflowstr(flow_time))
            if flow_time not in self.tflowlist:
                out_str = 'flow_time passed to function is: '+flow_time + '\n'
                out_str += 'not found in list \n' + str(self.tflowlist)
                raise IOError(out_str)
            flow_time = list(self.tflowlist).index(flow_time)
        for ikey in list(self.Op_cfgs['Op'].keys()):
            if ikey not in list(FlowOps2.Op_cfgs['Op'].keys()):
                raise IOError('Configuration lists differ, make sure .cfglist are identical')
            this_index.append(ikey)
            this_Op.append(thisOpFun(np.array(self.Op_cfgs['Op'][ikey]),FlowOps2.Op_cfgs['Op'][ikey][flow_time]))
        this_index = pa.MultiIndex.from_tuples(this_index,names=self.Op_cfgs.index.names)
        output.Op_cfgs = deepcopy(self.Op_cfgs)
        output.Op_cfgs.loc[:,'Op'] = pa.Series(this_Op,index=this_index)
        return output



    def FlowCombinePreBS(self,FlowOps2,Operation='x'):
        thisOpFun = GetOppFromStr(Operation)
        if not isinstance(FlowOps2,FlowOp): raise IOError('FlowCombinePreBS only takes FlowOp types, not '+FlowOps2.__class__.__name__)
        output = FlowOp(thisnboot=self.nboot, cfglist=self.Op_cfgs,Info=CombineInfo(self.flowInfo,FlowOps2.flowInfo,Operation.replace('*','x')),
                        thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        this_Op = []
        this_index = []
        # print 'DEBUG 3, '
        # print self, len(self.Op),
        # print FlowOps2,len(FlowOps2.Op)
        # print
        for ikey in list(self.Op_cfgs['Op'].keys()):
            if ikey not in list(FlowOps2.Op_cfgs['Op'].keys()):
                raise IOError('Configuration lists differ, make sure .cfglist are identical')
            this_index.append(ikey)
            this_Op.append(thisOpFun(np.array(self.Op_cfgs['Op'][ikey]),np.array(FlowOps2.Op_cfgs['Op'][ikey])))
        this_index = pa.MultiIndex.from_tuples(this_index,names=self.Op_cfgs.index.names)
        output.Op_cfgs = deepcopy(self.Op_cfgs)
        output.Op_cfgs.loc[:,'Op'] = pa.Series(this_Op,index=this_index)
        return output




    def FlowImportPlotParams(self,thiscol,thissym,thisshift):
        self.FlowCheckCol(thiscol)
        self.FlowCheckSym(thissym)
        if thisshift != 'PreDefine': self.thisshift = thisshift

    def FlowBootstrap(self,tflowlist='PreDef',WipeData=True,show_timer=True,Improve=False,CheckCfgs=True):
        if tflowlist != 'PreDef': self.tflowlist = tflowlist
        ReadAll = (len(self.tflowlist) == 0)
        if 'Op' not in self.Op_cfgs:  self.FlowRead(show_timer=show_timer,CheckCfgs=CheckCfgs)
        bootl,avgl,stdl = [],[],[]
        keyl = []
        if show_timer: thistimer = Timer(linklist=self.tflowfloat,name='Bootstrapping '+self.flowname)
        rlist = None
        for ictflow,tflows in enumerate(self.tflowfloat):
            strit = tflowstr(tflows) ## XmlFormatting.py
            if strit in self.tflowlist or tflows in self.tflowlist or ReadAll:
                tdata = self.Op_cfgs['Op'].apply(lambda x: x[ictflow])
                ##ONLY FOR CHIT
                if Improve:
                    bootl.append(BootStrap(self.nboot, name=self.flowfilepref+'('+strit+')',cfgvals=tdata,QforImprov=tdata,rand_list=rlist))
                else:
                    bootl.append(BootStrap(self.nboot, name=self.flowfilepref+'('+strit+')',cfgvals=tdata,rand_list=rlist))
                rlist = bootl[-1].Get_Rand_List()
                keyl.append(strit)
                avgl.append(bootl[-1].Avg)
                stdl.append(bootl[-1].Std)
                if show_timer: thistimer.Lap()
        if len(keyl) > 0:
            self.tflowlist = keyl
            index_list = pa.Index(keyl,name='Flow Times')
            self.Op_Stats['boot'] = pa.Series(bootl,index=index_list)
            self.Op_Stats['Avg'] = pa.Series(avgl,index=index_list)
            self.Op_Stats['Std'] = pa.Series(stdl,index=index_list)
        self.FlowAutocorr(tflowlist=tflowlist,show_timer=show_timer)
        if WipeData:
            del self.Op_cfgs['Op']

    def GetObsFun(self):
        def Ident(*val):
            return val[0]
        def IdentDer(*val):
            return [1.0]+[0.0]*(len(val)-1)
        return (Ident,IdentDer)

    def FlowAutocorr(self,tflowlist='PreDef',show_timer=True,CheckCfgs=True):
        # if not hasattr(self,'OpAuto'):
        #     return
        if tflowlist != 'PreDef': self.tflowlist = tflowlist
        ReadAll = (len(self.tflowlist) == 0)
        autol,aavgl,astdl = [],[],[]
        keyl = []
        if 'Op' not in self.Op_cfgs:  self.FlowRead(show_timer=show_timer,CheckCfgs=CheckCfgs)
        if show_timer: thistimer = Timer(linklist=self.tflowfloat,name='Autocorrelating '+self.flowname)
        for ictflow,tflows in enumerate(self.tflowfloat):
            strit = tflowstr(tflows) ## XmlFormatting.py
            if strit in self.tflowlist or tflows in self.tflowlist or ReadAll:
                tdata = self.Op_cfgs['Op'].apply(lambda x: x[ictflow]).to_frame('Op')
                keyl.append(strit)
                autol.append(AutoCorrelate( Fun=self.GetObsFun(),Sparam=self.Sparam,
                                            name=self.flowname + ' $' +tflowTOLeg(strit)+'$',data=tdata))
                aavgl.append(autol[-1].Avg)
                astdl.append(autol[-1].Std)
                if show_timer: thistimer.Lap()
        if len(keyl) > 0:
            self.tflowlist = keyl
            index_list = pa.Index(keyl,name='Flow Times')
            self.Op_Stats.loc[:,'Auto'] = pa.Series(autol,index=index_list)
            self.Op_Stats.loc[:,'AutoAvg'] = pa.Series(aavgl,index=index_list)
            self.Op_Stats.loc[:,'AutoStd'] = pa.Series(astdl,index=index_list)


    def Get_Extrapolation(self,fmted=True):
        if self.flowFit_Stats is None:
            # print 'No fits found for extrapolation, attempting to fit now.'
            self.FlowFit()
        if self.flowFit_Stats is None:
            # raise EnvironmentError('Flowed Operator Fit Failed')
            return None
        else:
            return self.flowFit_Stats.Get_Extrapolation(fmted=fmted)



    ## knoledge of number of parameters is needed
    def FlowSetFunction(self,Fun,npar):
        self.FlowFun = Fun
        if hasattr(Fun,'__name__'):
            self.FFun_name = Fun.__name__
        self.flownpar = npar


    ## Fitting the correlator
    ## flowfit_range is list of fit ranges to fit the correlator to
    ## formatted as fitr#-# where # are the min and max fit values.
    ## states correspond to how many states to fit to (currenlty only working for states = 1 or 2)
    def FlowFit(self,flowfit_range='PreDef',FlowiGuess='PreDef',EstDir=False,WipeFit=True):

        if self.FlowFun == 'None': raise IOError('Please define funciton using .SetFunction(Fun,npar) for fitting before calling .Fit()')

        if 'PreDef' != flowfit_range: self.flowfit_range = flowfit_range
        if isinstance(self.flowfit_range,(tuple,list,np.ndarray)):
            self.flowfit_range = self.flowfit_range[0]
        if str(self.flowfit_range) == 'None': return
        if FlowiGuess != 'PreDef': self.flowiGuess = FlowiGuess
        if self.flownpar == 1:
            self.thisparam = ['Value']
        else:
            self.thisparam = 'PreDef'

        fit_info = {}
        if EstDir:
            fit_info['Funs'] = (self.FlowFun,self.flownpar,'Estimate')
        else:
            fit_info['Funs'] = (self.FlowFun,self.flownpar)
        fit_info['iGuess'] = self.flowiGuess
        fit_info['name'] = self.flowname + ' Fits'
        fit_info['paramunits'] = 'GeV'
        fit_min,fit_max = unxmlfitr_int(self.flowfit_range)
        # print('DEBUG',fit_min,fit_max)
        # print(self.Op_Stats['boot'])
        self.flowFit_Stats = sff.SetOfFitFuns(data=self.Op_Stats['boot'])
        this_test = self.flowFit_Stats.ScanRange(fit_min,fit_max,fit_info=fit_info,min_fit_len=0)
        if this_test:
            self.flowFit_Stats.DoFits()
            self.flowFit_Stats.SortChi()
            self.FlowWrite()
        else:
            self.flowFit_Stats = None

    ## Comparisons

    def __str__(self):
        return self.flowname


    def __iter__(self):
        return self

    def Flowitems(self):
        return list(self.Op_Stats['boot'].items())



    def FlowitemsAvgStd(self):
        outlist = []
        for it,tdata in self.Op_Stats['boot'].items():
            outlist.append((it,tdata.Avg,tdata.Std))
        return outlist


    def Flowvalues(self):
        return self.Op_Stats['boot'].values

    def FlowvaluesAvgStd(self):
        return list(zip(self.Op_Stats['Avg'].values,self.Op_Stats['Std'].values))


    def Flowkeys(self):
        return self.Op_Stats['boot'].index


    def __setitem__(self, key, value ):
        self.Op_Stats.loc[key,'boot'] = value

    def __getitem__(self, key):
        return self.Op_Stats.loc[key,'boot']

    def __reversed__(self):
        self.Op_Stats.iiloc[::-1]
        return self



    def Flownext(self):
        if self.flowcurr >= len(self.Flowkeys()):
            self.flowcurr = 0
            raise StopIteration
        else:
            self.flowcurr += 1
            if self.flowcurr >= len(self.Flowkeys())+1:
                return 0.0
            else:
                thiskey = self.Flowkeys()[self.flowcurr-1]
                return self[thiskey]

    ## incase next gets overloaded
    def __next__(self):
        return self.Flownext(self)

    def FlowWrite(self,ExtraInfo='None'):
        def FixDictArray(thislist,ikey):
            outDict = OrderedDict()
            for ic,ilist in enumerate(thislist):
                outDict[ikey+str(ic+1)] = ilist
            return outDict

        outDict = ODNested()
        if not hasattr(self,'show_cfgs'):
            self.show_cfgs = False



        outDict['name'] = self.flowname
        outDict['kud'] = self.kud
        outDict['ks'] = self.ks
        outDict['nt'] = self.nt
        outDict['nxyz'] = self.nxyz

        if hasattr(self,'flowncfg_list'):
            for istream,incfg in zip(self.stream_list,self.flowncfg_list):
                outDict['ncfg_'+istream] = incfg
        outDict['nboot'] = self.nboot
        if hasattr(self,'chitcoeff'):
            outDict['chit_coeff'] = self.chitcoeff

        if hasattr(self,'t0_values'):
            for irow,rowdata in self.t0_values.iterrows():
                this_key = rowdata['tf_min'] + '_to_'+rowdata['tf_max']
                outDict[this_key]['t0'] = AvgStdToFormat(rowdata['boot'].Avg,rowdata['boot'].Std,frmtflag='f')
                outDict[this_key]['sqrt8tf'] = AvgStdToFormat(rowdata['sqrt8tf'].Avg,rowdata['sqrt8tf'].Std,frmtflag='f')


        excel_params = pa.Series(deepcopy(outDict))




        if 'boot' in self.Op_Stats.columns:
            for it,tdata in self.Op_Stats['boot'].items():
                outDict['Opboot'][it] = AvgStdToFormat(tdata.Avg,tdata.Std,frmtflag='e')
        # if hasattr(self,'FlowFun'):
        #     outDict['FlowFits']['Function'] = self.FlowFun.__name__
        # elif hasattr(self,'FlowFun_name'):
        #     outDict['FlowFits']['Function'] = self.FlowFun_name
        # if self.flowiGuess == 'None':
        #     outDict['FlowFits']['Initial_Guess'] = 'Default (see FitFunctions.py)'
        # else:
        #     for iparam in xrange(1,self.flownpar+1):
        #         thispar = 'Par'+str(iparam)
        #         outDict['FlowFits']['Initial_Guess'][thispar] = self.flowiGuess[iparam]
            # for thispar,parbs in fitdict.Params.iteritems():
            #     outDict['FlowFits'][ifit][thispar] = AvgStdToFormat(parbs.Avg,parbs.Std)
            # outDict['FlowFits'][ifit]['Chi^2_pdf'] = AvgStdToFormat(fitdict.Chi2DoF.Avg,fitdict.Chi2DoF.Std)
            # outDict['FlowFits'][ifit]['MaxIter'] = fitdict.MaxIter
            # outDict['FlowFits'][ifit]['Prec'] = fitdict.Prec
        if 'Auto' in self.Op_Stats.columns:
            outDict['S_param'] = self.Sparam
            for it,tdata in self.Op_Stats['Auto'].items():
                outDict['OpAuto'][it] = AvgStdToFormat(tdata.Avg,tdata.Std,frmtflag='e')
                outDict['OpAuto_Wopt'][it] = tdata.Wopt
                outDict['OpAuto_tau'][it] = AvgStdToFormat(tdata.tau,tdata.tauerr,frmtflag='e')

        if self.flowFit_Stats is not None:
            outDict['FlowFits'] = self.flowFit_Stats.GetOutputDict()


        if self.show_cfgs and hasattr(self,'flowncfg_list') and any(np.array(self.flowncfg_list) > 0) and 'configs' in self.Op_cfgs.columns:
            outDict['cfglist'] = Series_TO_ODict(self.Op_cfgs['configs'])

        if ExtraInfo != 'None':
            for ikey,ival in ExtraInfo.items():
                if ikey not in list(outDict.keys()):
                    outDict[ikey] = ival



        ## Write human readable file with data
        WriteXml(self.flowHumanFile,{'Results':outDict})

        ## write to excel file
        if 'configs' in self.Op_cfgs.columns:
            WriteExcel(self.flowExcelFile,{'Op_results':deepcopy(self.Op_Stats)},flowcfglist=self.Op_cfgs['configs'],params=excel_params)

            if hasattr(self,'flowncfg_list') and any(np.array(self.flowncfg_list) > 0):
                outDict['cfglist'] = Series_TO_ODict(self.Op_cfgs['configs'])

        ## pickles rest of data for reading
        self.RemoveFuns()
        WritePickle(self.flowPickleFile,self.__dict__)
        self.GetFuns()

    def GetFuns(self):
        if hasattr(self.flowFit_Stats,'GetFuns'):
            self.flowFit_Stats.GetFuns()
        if not hasattr(self,'FlowFun'):
            self.FlowFun = ReadFuns(self.FlowFun_name)[0]
            del self.FlowFun_name

    def RemoveFuns(self):
        if hasattr(self.flowFit_Stats,'RemoveFuns'):
            self.flowFit_Stats.RemoveFuns()
        if hasattr(self,'FlowFun') and hasattr(self.FlowFun,'__name__'):
            self.FlowFun_name = self.FlowFun.__name__
            WriteFuns(self.FlowFun)
            del self.FlowFun

    ##################################

    def FlowReadAndWrite(self,tflowlist='PreDef',WipeData=False,show_timer=True,CheckCfgs=True):
        self.FlowRead(show_timer=show_timer,CheckCfgs=CheckCfgs)
        self.FlowBootstrap(tflowlist=tflowlist,WipeData=WipeData,CheckCfgs=CheckCfgs)
        self.FlowFit()
        self.FlowWrite()

    def FixRead(self,readdict):
        readdict['flowInfo']['outdir'] = self.flowInfo['outdir']
        readdict['flowfit_range'] = self.flowfit_range
        readdict['flowPickleFile'] = self.flowPickleFile
        readdict['flowHumanFile'] = self.flowHumanFile
        readdict['flowExcelFile'] = self.flowExcelFile
        readdict['flowCfgFile'] = self.flowCfgFile
        readdict['flowLegLab'] = self.flowLegLab
        readdict['cfg_file_type'] = self.cfg_file_type
        # readdict['read_cfgs_from_file'] = self.read_cfgs_from_file
        if 'is_comb' not in list(readdict.keys()):
            readdict['is_comb'] = self.is_comb
        if 'cfg_file_type' not in list(readdict.keys()):
            readdict['cfg_file_type'] = self.cfg_file_type
        if 'FFun_name' in list(readdict.keys()) and readdict['FFun_name'] != self.FFun_name:
            print('Warning, fit function in file')
            print(readdict['FFun_name'])
            print('is not what was asked for')
            print(self.FFun_name)
            readdict['FFun_name'] = self.FFun_name
        if 'FlowNewForm' not in list(readdict.keys()):
            readdict['FlowNewForm'] = False
        return readdict


    def FlowLoadPickle(self,DefWipe = False,WipeData=False,CheckCfgs=True,show_timer=True):
        # print 'Loading Pickle for ' , self.flowPickleFile
        if os.path.isfile(self.flowPickleFile) and not DefWipe:
            loadeddict = ReadPickleWrap(self.flowPickleFile)
            loadeddict = self.FixRead(loadeddict)
            thischeck = self.thischecklist
            ## TODO: implement cfglist check with new style
            if CheckCfgs: thischeck = thischeck + ['Op_cfgs']
            # checkdict = deepcopy(self.__dict__)
            # checkdict['flowPickleFile'] = self.flowPickleFile.replace('.bak','')
            # checkdict['flowExcelFile'] = self.flowExcelFile.replace('.bak','')
            # checkdict['flowHumanFile'] = self.flowHumanFile.replace('.bak','')
            if CheckClass(self.__dict__,loadeddict,thischeck):
                self.read_cfgs_from_file = True
                self.__dict__.update(loadeddict)
                if Wipe_All_Fits or not isinstance(self.flowFit_Stats,sff.SetOfFitFuns):
                    self.flowFit_Stats = None
            else:
                if os.path.isfile(self.flowPickleFile+'.bak'):
                    print('using backupfile for '+self.flowname)
                    self.flowPickleFile = self.flowPickleFile+'.bak'
                    self.flowHumanFile = self.flowHumanFile+'.bak'
                    self.flowExcelFile = self.flowExcelFile.replace('.xlsx','.bak.xlsx')
                    self.FlowLoadPickle(DefWipe=DefWipe,WipeData=WipeData,CheckCfgs=CheckCfgs,show_timer=show_timer)
                    return
                print('Warning, file ' , self.flowHumanFile , ' has different parameters than this instance, reading and writing over file now')
                print()
                self.FlowReadAndWrite(WipeData=WipeData,show_timer=show_timer,CheckCfgs=CheckCfgs)
            # self.FlowWrite()
        else:
            self.FlowReadAndWrite(WipeData=WipeData,show_timer=show_timer,CheckCfgs=CheckCfgs)
    ## Operator overloading

    #Too Much work, only works for one operator. Please use SetCustomName after doing all algebraic manipulation

    ## this could to with some improving.... or some more generallity somewhere

    def FlowUpdateName(self,opp,Op,Op2):
        def GetOppSelf(opp):
            if opp == '+':
                return 'x2'
            elif opp == '-':
                return 'x0'
            elif opp == 'x':
                return 'pow2'
            elif opp == 'div':
                return 'divself'
            elif opp == 'pow':
                return 'powself'
        def GetLegSelf(opp):
            if opp == '+':
                return '\times 2'
            elif opp == '-':
                return '\times 0'
            elif opp == 'x':
                return '^{2}'
            elif opp == 'div':
                return 'divself'
            elif opp == 'pow':
                return 'powself'

        def LLFormatOp(opp):
            oppout = opp.replace('x',r'\times ')
            oppout = oppout.replace('pow','^{}')
            return oppout

        outname = []
        outLL = []
        oppLL = LLFormatOp(opp) ## \times looks better :)
        if isinstance(Op2,FlowOp) and isinstance(Op,FlowOp):
            op1LL = Op.flowLegLab[1:-1] ## removes $ symbols, re adds at end
            op2LL = Op2.flowLegLab[1:-1] ## removes $ symbols, re adds at end
            if Op.flowname == Op2.flowname:
                outname = Op.flowname,GetOppSelf(opp)
                outLL = op1LL,' Comb ',GetLegSelf(oppLL)
            else:
                for iname,iname2 in zip(Op.flowname.split('_'),Op2.flowname.split('_')):
                    if iname in iname2:
                        outname.append(iname2)
                    elif iname2 in iname:
                        outname.append(iname)
                    else:
                        outname.append(iname+opp+iname2)
                outname = '_'.join(outname)
                for iflowLegLab,iflowLegLab2 in zip(op1LL.split('_'),op2LL.split('\ ')):
                    if iflowLegLab in iflowLegLab2:
                        outLL.append(iflowLegLab2)
                    elif iflowLegLab2 in iflowLegLab:
                        outLL.append(iflowLegLab)
                    else:
                        outLL.append(iflowLegLab+oppLL+iflowLegLab2)
                outLL = '_'.join(outLL)
        elif isinstance(Op,FlowOp):
            op1LL = Op.flowLegLab[1:-1] ## removes $ symbols, re adds at end
            if isinstance(Op2,BootStrap):
                outname = Op.flowname,opp+Op2.name
                outLL  =  op1LL,oppLL+Op2.name
            elif isinstance(Op2,int):
                outname  =  Op.flowname,opp+str(Op2)
                outLL  =  op1LL,oppLL+str(Op2)
            else:
                try:
                    if Op2 > 0.01:
                        outname = Op.flowname,opp+'{:.2f}'.format(Op2)
                        outLL = op1LL,oppLL+'{:.2f}'.format(Op2)
                except Exception as err:
                    print(type(Op2))
                    raise EnvironmentError('Invalid value to combine with FlowOp class \n'+str(err))
        elif isinstance(Op2,FlowOp):
            op2LL = Op2.flowLegLab[1:-1] ## removes $ symbols, re adds at end
            if isinstance(Op,BootStrap):
                outname  =  Op.name + opp ,Op2.flowname
                outLL  =  Op.name + oppLL,op2LL
            elif isinstance(Op,int):
                outname  =  str(Op) + opp ,Op2.flowname
                outLL  =  str(Op) + oppLL,op2LL
            else:
                if Op > 0.01:
                    outname = '{:.2f}'.format(Op) + opp,Op2.flowname
                    outLL = '{:.2f}'.format(Op) + oppLL,op2LL
        else:
            if Op > 0.01:
                outname = '{:.2f}'.format(float(Op)) + opp,'{:.2f}'.format(float(Op2))
                outLL = '{:.2f}'.format(float(Op)) + oppLL,'{:.2f}'.format(float(Op2))

        # print r'$'+r'\ '.join(outLL)+'$'
        outLL = list(outLL)
        for ic,iLL in enumerate(outLL):
            if iLL is not None and '^{}' in iLL:
                outLL[ic] = iLL.replace('^{}','^{')+'}'
        if outname is not None:
            if outLL is not None:
                outLL = np.array(outLL)
                outLL = outLL[outLL is None]
                self.FlowSetCustomName('_'.join(outname),r'$'+r'\ '.join(outLL)+'$')
            else:
                self.FlowSetCustomName('_'.join(outname))

    ## Numerical operatoions

    def __add__(self,Op2):
        if any([ipar in str(type(Op2)) for ipar in FlowOp.parentlist]):
            return Op2.__radd__(self)
        else:
            result = FlowOp(thisnboot=self.nboot, cfglist=self.Op_cfgs,
                                  Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
            result.FlowUpdateName('+',self,Op2)
            result.Op_Stats = pa.DataFrame(columns=['boot'],index=self.Op_Stats.index)
            if isinstance(Op2,FlowOp):
                for (it,iOp) in self.Flowitems():
                    result[it] = iOp + Op2[it]
            else:
                try:
                    for it,iOp in self.Flowitems():
                        result[it] = iOp + Op2
                except:
                    print(type(Op2))
                    raise EnvironmentError('Invalid value to combine with FlowOp class')
            result.is_comb = True
            return result

    def __sub__(self,Op2):
        if any([ipar in str(type(Op2)) for ipar in FlowOp.parentlist]):
            return Op2.__rsub__(self)
        else:
            result = FlowOp(thisnboot=self.nboot, cfglist=self.Op_cfgs,
                                  Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
            result.FlowUpdateName('-',self,Op2)
            result.Op_Stats = pa.DataFrame(columns=['boot'],index=self.Op_Stats.index)
            if isinstance(Op2,FlowOp):
                for (it,iOp) in self.Flowitems():
                    result[it] = iOp - Op2[it]
            else:
                try:
                    for it,iOp in self.Flowitems():
                        result[it] = iOp - Op2
                except:
                    print(type(Op2))
                    raise EnvironmentError('Invalid value to combine with FlowOp class')
            result.is_comb = True
            return result

    def __mul__(self,Op2):
        if any([ipar in str(type(Op2)) for ipar in FlowOp.parentlist]):
            return Op2.__rmul__(self)
        else:
            result = FlowOp(thisnboot=self.nboot, cfglist=self.Op_cfgs,
                                  Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
            result.FlowUpdateName('x',self,Op2)
            result.Op_Stats = pa.DataFrame(columns=['boot'],index=self.Op_Stats.index)
            if isinstance(Op2,FlowOp):
                for (it,iOp) in self.Flowitems():
                    result[it] = iOp * Op2[it]
            else:
                try:
                    for it,iOp in self.Flowitems():
                        result[it] = iOp * Op2
                except:
                    print(type(Op2))
                    raise EnvironmentError('Invalid value to combine with FlowOp class')
            result.is_comb = True
            return result

    def __div__(self,Op2):
        if any([ipar in str(type(Op2)) for ipar in FlowOp.parentlist]):
            return Op2.__rdiv__(self)
        else:
            result = FlowOp(thisnboot=self.nboot, cfglist=self.Op_cfgs,
                                  Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
            result.FlowUpdateName('div',self,Op2)
            result.Op_Stats = pa.DataFrame(columns=['boot'],index=self.Op_Stats.index)
            if isinstance(Op2,FlowOp):
                for (it,iOp) in self.Flowitems():
                    try:
                        result[it] = iOp / Op2[it]
                    except:
                        if any([ibs == 0 for ibs in Op2[it].bootvals]):
                            raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+Op2.flowname + ' '+Op2[it].flowname)

            else:
                # if Op2 == 0:
                #     raise ZeroDivisionError('Dividing by zero found in bootstrap division')
                try:
                    for it,iOp in self.Flowitems():
                        result[it] = iOp / Op2
                except:
                    print(type(Op2))
                    raise EnvironmentError('Invalid value to combine with FlowOp class')
            result.is_comb = True
            return result


    def __pow__(self,Op2):
        if any([ipar in str(type(Op2)) for ipar in FlowOp.parentlist]):
            return Op2.__rpow__(self)
        else:
            result = FlowOp(thisnboot=self.nboot, cfglist=self.Op_cfgs,
                                  Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
            result.FlowUpdateName('pow',self,Op2)
            result.Op_Stats = pa.DataFrame(columns=['boot'],index=self.Op_Stats.index)
            if isinstance(Op2,FlowOp):
                for (it,iOp) in self.Flowitems():
                    result[it] = iOp ** Op2[it]
            else:
                try:
                    result.Op_Stats.loc[:,'boot'] = pa.Series(self.Op_Stats['boot'] ** Op2,index=self.tflowlist)
                except:
                    print(type(Op2))
                    raise EnvironmentError('Invalid value to combine with FlowOp class')

            result.is_comb = True
            return result


    ## Right multiplication functions


    def __radd__(self,Op2):
        result = FlowOp(thisnboot=self.nboot, cfglist=self.Op_cfgs,
                              Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.FlowUpdateName('+',Op2,self)
        result.Op_Stats = pa.DataFrame(columns=['boot'],index=self.Op_Stats.index)
        if isinstance(Op2,FlowOp):
            for (it,iOp) in self.Flowitems():
                result[it] = Op2[it] + iOp
        else:
            try:
                for it,iOp in self.Flowitems():
                    result[it] = Op2 + iOp
            except:
                print(type(Op2))
                raise EnvironmentError('Invalid value to combine with FlowOp class')
        result.is_comb = True
        return result

    def __rsub__(self,Op2):
        result = FlowOp(thisnboot=self.nboot, cfglist=self.Op_cfgs,
                              Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.FlowUpdateName('-',Op2,self)
        result.Op_Stats = pa.DataFrame(columns=['boot'],index=self.Op_Stats.index)
        if isinstance(Op2,FlowOp):
            for (it,iOp) in self.Flowitems():
                result[it] =  Op2[it] - iOp
        else:
            try:
                for it,iOp in self.Flowitems():
                    result[it] =  Op2 - iOp
            except:
                print(type(Op2))
                raise EnvironmentError('Invalid value to combine with FlowOp class')
        result.is_comb = True
        return result

    def __rmul__(self,Op2):
        result = FlowOp(thisnboot=self.nboot, cfglist=self.Op_cfgs,
                              Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.FlowUpdateName('x',Op2,self)
        result.Op_Stats = pa.DataFrame(columns=['boot'],index=self.Op_Stats.index)
        if isinstance(Op2,FlowOp):
            for (it,iOp) in self.Flowitems():
                result[it] =  Op2[it] * iOp
        else:
            try:
                result.Op_Stats.loc[:,'boot'] = pa.Series(self.Op_Stats['boot'] * Op2,index=self.tflowlist)
            except:
                print(type(Op2))
                raise EnvironmentError('Invalid value to combine with FlowOp class')
        result.is_comb = True
        return result

    def __rdiv__(self,Op2):
        result = FlowOp(thisnboot=self.nboot, cfglist=self.Op_cfgs,
                              Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.FlowUpdateName('div',Op2,self)
        result.Op_Stats = pa.DataFrame(columns=['boot'],index=self.Op_Stats.index)
        if isinstance(Op2,FlowOp):
            for (it,iOp) in self.Flowitems():
                try:
                    result[it] = Op2[it] / iOp
                except:
                    if any([ibs == 0 for ibs in iOp.bootvals]):
                        raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+self.flowname + ' '+iOp.flowname)

        else:
            for it,iOp in self.Flowitems():
                try:
                    result[it] = Op2 / iOp
                except:
                    if iOp == 0:
                        raise ZeroDivisionError('Dividing by zero found in bootstrap division')
                    else:
                        print(type(Op2))
                        raise EnvironmentError('Invalid value to combine with FlowOp class')
        result.is_comb = True
        return result


    def __rpow__(self,Op2):
        result = FlowOp(thisnboot=self.nboot, cfglist=self.Op_cfgs,
                              Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.FlowUpdateName('pow',Op2,self)
        result.Op_Stats = pa.DataFrame(columns=['boot'],index=self.Op_Stats.index)
        if isinstance(Op2,FlowOp):
            for (it,iOp) in self.Flowitems():
                result[it] =  Op2[it] ** iOp
        else:
            try:
                for it,iOp in self.Flowitems():
                    result[it] = Op2 ** iOp
            except:
                print(type(Op2))
                raise EnvironmentError('Invalid value to combine with FlowOp class')
        result.is_comb = True
        return result


    ## unitary arithmetic operations

    def __neg__(self):
        for it in self.Flowkeys():
            self[it] = -self[it]
        return self


    def __pos__(self):
        return self

    def __abs__(self):
        for it in self.Flowkeys():
            self[it] = abs(self[it])
        return self

    ## Casting, just casts the bootstrap values, ignore return value (only needed since return has to be same type).


    def __complex__(self):
        for it in self.Flowkeys():
            complex(self[it])
        self.FlowStats()
        return 1+0j

    def __int__(self):
        for it in self.Flowkeys():
            int(self[it])
        self.FlowStats()
        return 1

    def __float__(self):
        for it in self.Flowkeys():
            float(self[it])
        self.FlowStats()
        return 1.0






######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################
######################################################################################################################################################







class FlowOpFullSquared(object):
    """ Op(tf,t) uses bootstrap class
        tf = flow time
        t = time slice
    """

    ## Info is a dictionary containing information for the correlator:
    ## see comments below to see what is required.
    ## missing elements just ignores that field

    ## TODO: create variational flag and interface with this function.

    ## array of all classes that inherit this class. Append this when adding new classes please!
    parentlist = [  'TwoPtCorrelators.TwoPointCorr','TwoPtCorrelators.NNQCorr',
                    'TwoPtCorrelators.NNQFullCorr','VarMethod.VariationalTwoPt',
                    'RatioCorrelators.RatioCorr',
                    'RatioCorrelators.RatioFOCorr','RatioCorrelators.RatioFOFullCorr']

    def __init__(self, thisnboot=nboot, cfglist={},Info={}, thissym='Not Set',thiscol='Not Set',thisshift=0.0,name='',LLname=''
                 ,man_load_cfgs=False):
        def GetFolderFileOp(thistype):
            if self.FlowNewForm:
                if  'TopCharge' in thistype  :
                    return 'FlowNewForm','q'
                elif  'Weinberg' in thistype :
                    return 'FlowNewForm','W'
                else:
                    print('thistype ' , thistype)
                    raise IOError('Type of operator not known, please select TopCharge or Weinberg')
            else:
                if  'TopCharge' in thistype  :
                    return 'Flow','Q'
                elif  'Weinberg' in thistype :
                    return 'Flow','W'
                else:
                    print('thistype ' , thistype)
                    raise IOError('Type of operator not known, please select TopCharge or Weinberg')

        ## these class elemnts are checked when loading in pickled file

        # self.thischecklist = ['nboot','kud','ks','tflowlist','cfglist']
        self.thischecklist = ['nboot','kud','ks','FlowNewForm']

        self.latparams = mp.LatticeParameters(Info=Info)
        self.latparams.LoadPickle()
        self.nt = self.latparams.nt
        self.nxyz = self.latparams.nxyz
        self.dim_label = 'RC'+str(self.nxyz)+'x'+str(self.nt)
        self.dim_ll = str(self.nxyz)+r'^{3}\times ' + str(self.nt)

        self.flowcurr = 0
        ## Initialising plotting parameters
        self.thisshift = thisshift
        self.thiscol = thiscol
        self.thissym = thissym


        ## Op.... are the individual results for each time and flow time
        # #self.Op = [ itf , it , iconf ]
        # self.Op     = np.array([])

        self.Op_cfgs  = pa.DataFrame()
        '''
        self.Op_cfgs DataFrame:

        columns =   config_numbers , Op[itflow,it]

        rows =      stream , configuration
        '''


        ## when bootstrapping, the random source selection is done as well
        self.Op_Stats = pa.DataFrame()

        '''
        self.Op_Stats DataFrame:

        columns =   boot, Avg,      Std

        rows =      flow time, it_src
        '''

        self.Op_avg_Stats = pa.DataFrame()

        '''
        self.Op_avg_Stats DataFrame:

        columns =   boot, Avg,      Std
                    Auto, AutoAvg,  AutoStd

        rows =      flow time
        '''

        ## <Q^2> data
        self.Chi_Stats = pa.DataFrame()
        '''
        self.Chi_Stats DataFrame:

        columns =   boot, Avg,      Std
                    Auto, AutoAvg,  AutoStd

        rows =      flow time
        '''



        ## <Q(0)Q(t)> data
        # #self.Op2 = [ itf , it_distance , itsum , it_src , iconf ]
        # self.Op2     = np.array([])

        self.Op2_Stats = pa.DataFrame()
        '''
        self.Op2_Stats DataFrame:

        columns =   boot, Avg,      Std

        rows =      flow time , time separation , sum time, source time
        '''

        ## averaged over source locations
        self.Op2_avg_Stats = pa.DataFrame()
        '''
        self.Op2_avg_Stats DataFrame:

        columns =   boot, Avg,      Std
                    Auto, AutoAvg,  AutoStd

        rows =      flow time , time separation , sum time

        '''

        ## averaged over source locations
        self.Prop_Fit_Col_Names = ['sum_time','flow_time']
        self.Prop_Fit_Stats = pa.DataFrame()
        '''
        self.Prop_Fit_Stats DataFrame:

        columns =   boot, Avg,      Std

        rows =      sum time , flow time

        type = SetOfFitFuns

        '''

        self.nboot = thisnboot

        ## for autocorrelation analysis
        if 'Sparam' in list(Info.keys()): self.Sparam = Info['Sparam']
        else: self.Sparam = defSparam

        if 'show_cfgs' in list(Info.keys()): self.show_cfgs = Info['show_cfgs']
        else: self.show_cfgs = False


        ## kappa up down (in integer form)
        if 'kud' in list(Info.keys()): self.kud = 'kud'+str(Info['kud'])
        else: self.kud = ''

        ## kappa strange (in integer form)
        if 'ks' in list(Info.keys()): self.ks = 'ks'+str(Info['ks'])
        else: self.ks = ''

        ## Current obersvables are TopCharge and Weinberg
        if 'Observable' in list(Info.keys()): self.Observ = Info['Observable']
        else: self.Observ = 'TopCharge' ## default to TopCharge

        if 'cfg_file_type' in list(Info.keys()): self.cfg_file_type = Info['cfg_file_type']
        else: self.cfg_file_type = cfg_file_type

        ## this is used display only these flowtimes, Still reads and writes all
        if 'tflowlist' in list(Info.keys()): self.tflowlist = Info['tflowlist']
        else: self.tflowlist = defxlimOp ## default to read all flow times



        ## number of random sources per gauge field to average over
        if 'nrandt' in list(Info.keys()): self.nrandt = Info['nrandt']
        else: self.nrandt = self.latparams.nt

        if 'Op_fold' in list(Info.keys()): self.Op_fold = Info['Op_fold']
        else: self.Op_fold = True

        ## number of random sources per gauge field to average over
        if 'nranget' in list(Info.keys()): self.nranget = Info['nranget']
        else: self.nranget = self.latparams.nt/2


        ## number of values summed after to.
        ## default is 1 (for single source) and all (to give back original Q^2 result)
        if 'tsum_list' in list(Info.keys()) and Info['tsum_list'] != 'All': self.tsum_list = Info['tsum_list']
        else: self.tsum_list = ['ts'+str(it) for it in range(self.latparams.nt)]

        ## include the calculation of <Q^2>
        if 'Do_Chi' in list(Info.keys()): self.do_chi = Info['Do_Chi']
        else: self.do_chi = False

        ## selects if the sum is taken or not.
        if 'Sum_Range' in list(Info.keys()): self.do_sum_range = Info['Sum_Range']
        else: self.do_sum_range = True

        ## selects if the sum is taken or not.
        if 'Sublat_Method' in list(Info.keys()): self.do_sublat = Info['Sublat_Method']
        else: self.do_sublat = False

        if self.do_sublat:
            self.prop_fit_fun = [LinearFitFun,2]
        else:
            self.prop_fit_fun = [LinearFitFun,2]
            # self.prop_fit_fun = [Chi_Prop_Fit_2exp,5]


        ## will attempt to find all replica streams in directory.
        ## set self.nstream to number of streams you want to find.
        if 'stream_list' in list(Info.keys()):
            if isinstance(Info['stream_list'],str):
                if Info['stream_list'].lower() != 'all':
                    self.stream_list = [Info['stream_list']]
                else:
                    self.stream_list = 'all'
            else:
                self.stream_list = [iinfo.lower() for iinfo in Info['stream_list']]
        else:
            self.stream_list = 'all'

        self.FlowNewForm = flow_new_format
        # if 'FlowNewForm' in Info.keys(): self.FlowNewForm = Info['FlowNewForm']
        # else: self.FlowNewForm = False ## default to read all flow times


        # print 'DEBUG'
        # print 'tsumlist',self.tsum_list
        # print 'nrandt',self.nrandt

        listout = []
        for ict,itflow in enumerate(self.tflowlist):
            if 't_f' not in str(itflow):
                listout.append(tflowstr(itflow))
            else:
                listout.append(itflow)
        self.tflowfloat = [untflowstr(itflow) for itflow in listout]
        self.tflowlist = np.array(listout)

        thisfolder,self.flowfilepref = GetFolderFileOp(self.Observ)
        self.flowfilename = [self.flowfilepref+'_flow_b1.90_ng','.out']

        kappafolder = 'Kud0'+self.kud.replace('kud','')+'Ks0'+self.ks.replace('ks','')
        self.flowfilelist = ['_'+self.flowfilepref+'t_t{:2.6f}'.format(untflowstr(itflow)) for itflow in self.tflowlist]
        # self.flowfilename = ['RC32x64_B1900'+kappafolder+'C1715','_','_k'+self.kud.replace('kud','')+'_'+self.tsrc+filejsm+fileism+'_nucleon.2cf.xml']
        self.flowfiledir = datadir+thisfolder+'/'+self.dim_label+'_'+kappafolder+'-*-/'

        if self.do_chi:
            if self.Observ == 'TopCharge':
                self.chi_coeff = chitcoeff(self.latparams.nt,self.latparams.nxyz,self.latparams.latspace,'QQ')
                self.thispow = 0.25
            elif self.Observ == 'Weinberg':
                self.chi_coeff = chitcoeff(self.latparams.nt,self.latparams.nxyz,self.latparams.latspace,'WW')
                self.thispow = 0.125
            else:
                raise ImportError('Observable not recognised')
        else:
            self.chi_coeff = 1.
            self.thispow = 1.

        self.Create_Streams_Dirs()
        self.flowInfo = Info

        self.FlowSetCustomName(name,LLname)
        if not man_load_cfgs:
            self.LoadCfgs(cfglist,name,LLname)

    def LoadCfgs(self,cfglist,name='',LLname=''):
        ## cfglist { configuration : [ x-source numbers ] }
        if isinstance(cfglist,(list,np.ndarray,pa.DataFrame)):
            self.FlowImportCfgList(cfglist)
        else:
            self.FlowGetCfgList()

        ## predefine a time picking
        ## self.timelist [ it0_rand , istream, icfg]
        if 'PickedTimes' in list(self.flowInfo.keys()): self.ImportTimeList(self.flowInfo['PickedTimes'])
        else: self.MakeRandomTimeList()
        self.FlowSetCustomName(name,LLname)


    def GetFuns(self):
        if self.Prop_Fit_Stats is not None and 'boot' in self.Prop_Fit_Stats.columns:
            for ival in self.Prop_Fit_Stats['boot'].values:
                ival.GetFuns()
        if not hasattr(self,'prop_fit_fun'):
            self.prop_fit_fun = [ReadFuns(self.prop_fit_fun_name)[0],self.prop_fit_npar]
            if hasattr(self,'prop_fit_der_name'):
                self.prop_fit_fun += [ReadFuns(self.prop_fit_der_name)[0]]

    def RemoveFuns(self):
        if self.Prop_Fit_Stats is not None and 'boot' in self.Prop_Fit_Stats.columns:
            for ival in self.Prop_Fit_Stats['boot'].values:
                ival.RemoveFuns()
        if hasattr(self,'prop_fit_fun'):
            self.prop_fit_fun_name = self.prop_fit_fun[0].__name__
            self.prop_fit_npar = self.prop_fit_fun[1]
            if len(self.prop_fit_fun) == 3:
                self.prop_fit_der_name = self.prop_fit_fun[2].__name__
                WriteFuns(self.prop_fit_fun[0],self.prop_fit_fun[3])
            else:
                WriteFuns(self.prop_fit_fun[0])
            del self.prop_fit_fun


    def Create_Streams_Dirs(self):
        if isinstance(self.stream_list,str) and self.stream_list == 'all':
            self.All_Streams()
        else:
            self.Check_Streams()
        self.stream_name = '-'+''.join([istr.replace('-','') for istr in self.stream_list])+'-'
        self.flowdir_list = {}
        for istream in self.stream_list:
            thisdir = self.flowfiledir.replace('-*-',istream)
            if os.path.isdir(thisdir):
                self.flowdir_list[istream] = thisdir
            else:
                self.flowdir_list[istream] = thisdir.replace(self.dim_label+'_','')

    def All_Streams(self):
        self.stream_list = []
        for ifolder in glob.glob(self.flowfiledir):
            self.stream_list.append(re.findall('-.*-',ifolder)[-1])
        for ifolder in glob.glob(self.flowfiledir.replace(self.dim_label+'_','')):
            this_stream = re.findall('-.*-',ifolder)[-1]
            if this_stream not in self.stream_list:
                self.stream_list.append(this_stream)
        # if len(self.stream_list) == 0:
        #     raise IOError(self.flowfiledir+'\n No Streams found')

    def Check_Streams(self):
        for istream in self.stream_list:
            this_check = os.path.isdir(self.flowfiledir.replace('-*-',istream))
            this_check = this_check or os.path.isdir(self.flowfiledir.replace('-*-',istream).replace(self.dim_label+'_',''))
            if not this_check:
                raise IOError(  self.flowfiledir.replace('-*-',istream) + '\n stream directory not found, as well as \n'+
                                self.flowfiledir.replace('-*-',istream).replace(self.dim_label+'_',''))




    def FlowReadAndWrite(self,DefWipe=False,WipeData=False,show_timer=True,OnlyQ=False,CheckCfgs=True):
        self.FlowRead(show_timer=show_timer,CheckCfgs=CheckCfgs)
        self.FlowBootstrap(DefWipe=DefWipe,OnlyQ=OnlyQ,CheckCfgs=CheckCfgs)
        if not OnlyQ:
            self.FlowWrite()

    def FixRead(self,readdict):
        readdict['flowInfo']['outdir'] = self.flowInfo['outdir']
        readdict['latparams'].outputdir = self.latparams.outputdir
        readdict['latparams'].PickleFile = self.latparams.PickleFile
        readdict['latparams'].HumanFile = self.latparams.HumanFile
        readdict['tflowlist'] = self.tflowlist
        readdict['nrandt'] = self.nrandt
        readdict['nranget'] = self.nranget
        readdict['tsum_list'] = self.tsum_list
        if hasattr(self,'flowncfg_list'):
            readdict['flowncfg_list'] = self.flowncfg_list
        readdict['Op_cfgs'] = self.Op_cfgs
        readdict['Op_Stats'] = self.Op_Stats
        readdict['Op_avg_Stats'] = self.Op_avg_Stats
        readdict['Chi_Stats'] = self.Chi_Stats
        readdict['Op2_Stats'] = self.Op2_Stats
        readdict['Op2_avg_Stats'] = self.Op2_avg_Stats
        # print 'DEBUGGING PERPOUSES, please comment line below'
        # if self.prop_fit_fun[0].__name__ != readdict['prop_fit_fun'][0].__name__:
        #     readdict['Prop_Fit_Stats'] = self.Prop_Fit_Stats
        #     readdict['prop_fit_fun'] = self.prop_fit_fun
        readdict['flowfilelist'] = self.flowfilelist
        readdict['flowfilepref'] = self.flowfilepref
        readdict['flowdir_list'] = self.flowdir_list
        readdict['flowPickleFile'] = self.flowPickleFile
        readdict['flowHumanFile'] = self.flowHumanFile
        readdict['flowExcelFile'] = self.flowExcelFile
        readdict['cfg_file_type'] = self.cfg_file_type
        readdict['FlowNewForm'] = self.FlowNewForm
        return readdict

    def FlowLoadPickle(self,DefWipe = False,WipeData=False,CheckCfgs=False,OnlyQ=False,show_timer=True):
        if os.path.isfile(self.flowPickleFile) and not DefWipe:
            print('Loading Pickle for ' , self.flowPickleFile)
            loadeddict = ReadPickleWrap(self.flowPickleFile)
            loadeddict = self.FixRead(loadeddict)
            # thischeck = self.thischecklist
            # if CheckCfgs: thischeck = thischeck + ['flowcfglist']
            # if CheckClass(self.__dict__,loadeddict,thischeck):
            self.read_cfgs_from_file = True
            self.__dict__.update(loadeddict)
            if Wipe_All_Fits:
                self.Prop_Fit_Stats = pa.DataFrame()
            elif len(self.Prop_Fit_Stats.values) > 0 and not isinstance(self.Prop_Fit_Stats.values[0][0],sff.SetOfFitFuns):
                print('Legacy file, fits need to be recomputed')
                self.Prop_Fit_Stats = pa.DataFrame()
            # print 'DEBUG'
            # print self.Prop_Fit_Stats
            # else:
            #     print 'Warning, file ' , self.flowHumanFile , ' has different parameters than this instance, reading and writing over file now'
            #     print
            # self.FlowWrite()
        self.FlowReadAndWrite(  DefWipe=DefWipe,WipeData=WipeData,show_timer=show_timer,
                                OnlyQ=OnlyQ,CheckCfgs=CheckCfgs)

    def FlowGetCfgList_New(self):
        def SplitCfg(ifile):
            icfg = re.findall('_ng.*',ifile)
            if len(icfg) == 0:
                print()
                print('WARNING: file does not have cfg number')
                print(ifile)
                print()
                return ''
            icfg = icfg[0].replace('_ng','').replace('.out','')
            return icfg

        cfglist,ilist = [],[]
        low_test,upper_test = self.flowfilepref.upper(),self.flowfilepref.lower()
        for istream in self.stream_list:
            idir = self.flowdir_list[istream]+'PerGF/'
            for icfg,ifile in enumerate(glob.glob(idir+'/'+self.flowfilepref+'*')):
                if low_test in ifile or upper_test in ifile:
                    ilist.append((istream,icfg))
                    cfglist.append(SplitCfg(ifile))
        if len(cfglist) == 0:
            print('No configs found for directory:')
            print(' , '.join(self.flowdir_list))
        else:
            indicies = pa.MultiIndex.from_tuples(ilist,names=['stream','config_number'])
            cfg_df = pa.DataFrame(cfglist,columns=['configs'],index=indicies)
            cfg_df.loc[:,'int_configs'] = pa.Series(list(map(int,cfglist)),index=indicies)
            cfg_df = cfg_df.sort_values('int_configs').sort_index(level='stream',sort_remaining=False)
            del cfg_df['int_configs']
            self.FlowImportCfgList(cfg_df)



    def FlowGetCfgList(self):
        if self.FlowNewForm:
            self.FlowGetCfgList_New()
            return
        def CheckCfgFolder(thisdir):
            for icheck in self.flowfilelist:
                if not os.path.isfile(thisdir+'/'+icheck):
                    return False
            return True

        cfglist,ilist = [],[]
        for istream in self.stream_list:
            idir = self.flowdir_list[istream]
            if len(glob.glob(idir+'/*')) == 1:
                idir = idir + 'Comb/'
            for icfg,ifile in enumerate(glob.glob(idir+'/*')):
                if 'cfg' in ifile and CheckCfgFolder(ifile+'/'):
                    ilist.append((istream,icfg))
                    cfglist.append(ifile.replace(idir+'cfg',''))
        if len(cfglist) == 0:
            print('No configs found for directory:')
            print(' , '.join(self.flowdir_list))
        else:
            indicies = pa.MultiIndex.from_tuples(ilist,names=['stream','config_number'])
            cfg_df = pa.DataFrame(cfglist,columns=['configs'],index=indicies)
            cfg_df.loc[:,'int_configs'] = pa.Series(list(map(int,cfglist)),index=indicies)
            cfg_df = cfg_df.sort_values('int_configs').sort_index(level='stream',sort_remaining=False)
            del cfg_df['int_configs']
            self.FlowImportCfgList(cfg_df)

    def Write_Cfgs(self,file_type='PreDef',show_timer=True,CheckCfgs=True):
        if hasattr(self,'read_cfgs_from_file') and self.read_cfgs_from_file:
            return
        if file_type == 'PreDef':
            file_type = self.cfg_file_type
        if 'to_' not in file_type:
            file_type = 'to_'+file_type
        this_file = self.flowCfgFile+fmt_file_type(file_type.replace('to_',''))
        if show_timer: thistimer = Timer(name='Writing Cfgs '+this_file)
        if len(self.Op_cfgs.index) == 0:
            self.FlowRead(show_timer=show_timer,CheckCfgs=CheckCfgs)
        # this_df,dump = self.Get_Flattened_Cfgs()
        this_t_list = list(range(len(self.Op_cfgs['Op'].iloc[0][0])))
        this_meta = np.array([self.tflowlist,this_t_list])
        WriteWithMeta(self.Op_cfgs,this_meta,this_file,file_type=file_type)
        if show_timer: thistimer.Stop()

    def ChopDepVarList(self,this_df,read_meta,show_err=True):
        read_tflow,read_t = read_meta
        this_t = list(range(self.nt))
        isin = True
        this_index = [slice(None),slice(None)]
        name_list = ['flow_times','euclidian_times']
        for ic,(ithis_list,icheck_list) in enumerate([(this_t,read_t),(self.tflowlist,read_tflow)]):
            # isin = isin and set(ithis_list).issubset(icheck_list)
            isin = isin and set(ithis_list).issubset(icheck_list)
            if isin:
                same_list = len(ithis_list) == len(icheck_list)
                if not same_list:
                    index_list = [np.where(np.array(icheck_list) == it)[0][0] for it in ithis_list]
                    this_index[ic] = index_list
            elif show_err:
                print('\n missmatch in ',name_list[ic])
                print(ithis_list)
                print(icheck_list)
                break
            else:
                break
        if isin:
            def Chop_fun(val):
                return np.array(val)[this_index].tolist()
            this_df.loc[:,'Op'] = this_df.loc[:,'Op'].apply(Chop_fun)
        return this_df,isin

    def Read_Cfgs(self,file_type='PreDef',show_timer=True,CheckCfgs=True):
        def CheckCfgList(cfg_data,read_tflowlist):
            if cfg_data is None  or len(cfg_data) == 0:
                return None,False
            if CheckCfgs and 'stream_cfgs' in self.Op_cfgs:
                second_class = pa.Series([cfg_data],index=['Op_cfgs'])
                cfg_bool = CheckClass(self.__dict__,second_class,['Op_cfgs'])
                # cfg_bool = np.array(self.Op_cfgs['stream_cfgs'].values) == np.array(cfg_data['stream_cfgs'].values)
                # if isinstance(cfg_bool,np.ndarray):
                #     cfg_bool = all(cfg_bool)
                # if show_timer and not cfg_bool:
                #     print '\n'+GetCfgMissmatch(  np.array(self.Op_cfgs['stream_cfgs'].values),
                #                             np.array(cfg_data['stream_cfgs'].values))
            else:
                cfg_bool = True
            if cfg_bool:
                cfg_data,outbool = self.ChopDepVarList(cfg_data,read_tflowlist,show_err=show_timer)
                return cfg_data,outbool
            else:
                return None,False


        if file_type == 'PreDef':
            file_type = self.cfg_file_type
        if 'read_' not in file_type:
            file_type = 'read_'+file_type
        this_file = self.flowCfgFile+fmt_file_type(file_type.replace('read_',''))
        if not os.path.isfile(this_file):
            if '-def-' in this_file or '--' in this_file:
                this_file = this_file.replace('-def-','-*-').replace('--','-*-')
                file_list = glob.glob(this_file)
                if len(file_list) == 0:
                    print('Cfg FNF: \n', this_file)
                    self.read_cfgs_from_file = False
                    return False
                else:
                    this_file = file_list[0]
        if not os.path.isfile(this_file):
            self.read_cfgs_from_file = False
            if show_timer: print('cfg FNF: ' + this_file)
            return False
        if show_timer: thistimer = Timer(name='Reading Cfgs '+this_file)
        meta_data,cfg_data = ReadWithMeta(this_file,file_type=file_type)
        # cfg_data = self.Format_Flattened_Cfgs(cfg_data)
        cfg_data = self.MakeStreamCfgs(cfg_data)
        cfg_data,check_bool = CheckCfgList(cfg_data,meta_data)
        if check_bool:
            self.FlowImportCfgList(cfg_data)
            if show_timer: thistimer = thistimer.Stop()
            self.read_cfgs_from_file = True
            return True
        else:
            if show_timer: print('\n file incorrect,', end=' ')
            self.flowCfgFile = self.flowCfgFile.replace('.cfgs','.cfgs.bak')
            if os.path.isfile(this_file.replace('.cfgs','.cfgs.bak')):
                print(' found backup file')
                return self.Read_Cfgs(file_type=file_type,show_timer=show_timer,CheckCfgs=CheckCfgs)
            else:
                print(' reading unformatted')
                self.read_cfgs_from_file = False
                return False

    def MakeRandomTimeList(self):
        if 'configs' not in self.Op_cfgs:
            # raise EnvironmentError('configuration list need to be read in before getting random time list')
            print('Warning: configuration list need to be read in before getting random time list')
        else:
            readnt = self.latparams.nt
            ## timelist [  irandt , istream, icfg ]
            timelist = []
            for istream in self.flowncfg_list:
                for icfg in range(istream):
                    timelist.append(np.random.randint(0,readnt,size=self.nrandt).tolist())
            self.Op_cfgs.loc[:,'random_time_list'] = pa.Series(timelist[:len(self.Op_cfgs.index)],index=self.Op_cfgs.index)

        # self.timelist = np.random.randint(0,readnt,(self.flowncfg,self.nrandt))

    def ImportTimeList(self,timelist_in):
        if 'configs' not in self.Op_cfgs:
            raise EnvironmentError('configuration list need to be read in before import random time list')
        if isinstance(timelist_in,pa.Series):
            self.Op_cfgs['random_time_list'] = timelist_in
        else:
            ## timelist_in must be of shape [ itsrc_number, istream, icfg]
            timelist = []
            for isrc in timelist_in:
                this_timelist = []
                for istream in isrc:
                    this_timelist += istream
                timelist.append(this_timelist)
            self.Op_cfgs.loc[:,'random_time_list'] = pa.Series(timelist,index=self.Op_cfgs.index)



    def FlowRead_New(self,show_timer=True,cfg_file_type= 'PreDef',CheckCfgs=True):
        def FlowReadTopCharge_New(thisfile):
            if os.path.isfile(thisfile):
                readfile = thisfile
            elif os.path.isfile(thisfile.replace('ng00','ng')):
                readfile = thisfile.replace('ng00','ng')
            else:
                print('Warning, file not found',thisfile)
                return [],[],[]
            # print 'DEBUG reading: ',readfile
            with open(readfile,'r') as ifile:
                tflowlist,oplist = [],[]
                euclid_tlist = []
                flow_line,skip_t = True,False
                for iline,line in enumerate(ifile):
                    if flow_line and len(line.strip()) != 0:
                        val = list(map(np.float64,line.split()))
                        if len(val) == 3:
                            tflow,topcharge,cmplxdump = val
                        elif len(val) == 2:
                            tflow,topcharge = val
                        else:
                            raise IOError('error with file: \n'+readfile+'\n on line: '+str(iline+1))
                        if tflowstr(tflow) in self.tflowlist:
                            tflowlist.append(np.float(tflow))
                            oplist.append([])
                            euclid_tlist.append([])
                            skip_t = False
                        else:
                            skip_t = True
                            # oplist.append(topcharge)
                        flow_line = False
                    elif len(line.strip()) == 0:
                        flow_line = True
                        skip_t = False
                    elif not skip_t:
                        val = list(map(np.float64,line.split()))
                        if len(val) == 3:
                            et,topcharge,cmplxdump = val
                        elif len(val) == 2:
                            et,topcharge = val
                        else:
                            raise IOError('error with file: \n'+readfile+'\n on line: '+str(iline+1))
                        oplist[-1].append(topcharge)
                        euclid_tlist[-1].append(et)
            tflowlist,oplist,euclid_tlist = zip(*[[y,x,z] for y,x,z in sorted(zip(tflowlist,oplist,euclid_tlist))])
            return list(tflowlist),list(oplist),list(euclid_tlist)

        if self.Read_Cfgs(file_type=cfg_file_type,show_timer=show_timer,CheckCfgs=CheckCfgs):
            # if 'random_time_list' not in self.Op_cfgs :
            self.MakeRandomTimeList()
            return

        thisdata = []
        if show_timer: thistimer = Timer(linklist=self.Op_cfgs['configs'],name='Read '+self.flowname)
        if len(self.Op_cfgs['configs'].values) == 0:
            raise IOError('No values found for ' +self.Observ)
        for (istream,iccfg),icfg in self.Op_cfgs['configs'].items():
            ifile = self.flowdir_list[istream]+'PerGF/'+icfg.join(self.flowfilename)
            tflow,topcharge,euclid_t = FlowReadTopCharge_New(ifile)

            missing_tf = []
            for itflow,et_list in zip(tflow,euclid_t):
                if itflow not in tflow or len(et_list) != self.nt:
                    missing_tf.append(itflow)
            if len(missing_tf) > 0:
                print()
                # print self.tflowfloat
                print(itflow)
                print(ifile)
                raise IOError('flowed operator does not contain all flow times, check flow time list or remove file if corrupt')
            thisdata.append(topcharge)
            if show_timer: thistimer.Lap()
        # thisdata = np.array(thisdata)
        if len(thisdata) != len(self.Op_cfgs.index):
            raise IOError('Failed reading (see above) caused error. fix files please.')
        ## self.Op_cfgs['Op'][istream,icfg] = [ itflow , it ]
        self.Op_cfgs.loc[:,'Op'] = pa.Series(thisdata,index=self.Op_cfgs.index)
        # self.Op = np.rollaxis(thisdata,0,thisdata.ndim)
        # ## self.Op = [ it , icfg , itflow ]
        # self.CreateOpCorr(show_timer=show_timer)
        ## self.Op = [ itflow , it , icfg  ]
        print('read raw cofgs complete:')
        print(str(self.Op_cfgs))
        if 'random_time_list' not in self.Op_cfgs :
            self.MakeRandomTimeList()
        self.Write_Cfgs(file_type=cfg_file_type,show_timer=show_timer,CheckCfgs=CheckCfgs)


    def FlowRead(self,show_timer=True,cfg_file_type= 'PreDef',CheckCfgs=True):
        def ReadTsInFlowOp(thisfile):
            if os.path.isfile(thisfile):
                with open(thisfile,'r') as ifile:
                    tlist,topchargelist = [],[]
                    for line in ifile:
                        hold = list(map(np.float64,line.split()))
                        if len(hold) < 3: continue
                        tval,topcharge,cmplxdump = hold
                        tlist.append(tval)
                        topchargelist.append(topcharge)
                return tlist,topchargelist,False
            else:
                print('Warning, file not found',thisfile)
                return [],[],True


        ## making interpolator number more human readable
        # if len(self.flowname.split('_')) > 3:
        #     raise IOError('Do not read after combining correlators, please recreate the correlator')
        # ntflow = 0
        if self.FlowNewForm:
            self.FlowRead_New(show_timer=show_timer,cfg_file_type=cfg_file_type,CheckCfgs=CheckCfgs)
            return

        if self.Read_Cfgs(file_type=cfg_file_type,show_timer=show_timer,CheckCfgs=CheckCfgs):
            # if 'random_time_list' not in self.Op_cfgs :
            self.MakeRandomTimeList()
            return

        thisdata = []
        if show_timer: thistimer = Timer(linklist=self.Op_cfgs['configs'],name='Read '+self.flowname)
        if len(self.Op_cfgs['configs'].values) == 0:
            raise IOError('No values found for ' +self.Observ)
        for (istream,iccfg),icfg in self.Op_cfgs['configs'].items():
            istream_dir = self.flowdir_list[istream]
            idir = 'cfg'+icfg
            tflow,thisop = [],[]
            for itflowfile in self.flowfilelist:
                holdtflow,holdtop,holdError = ReadTsInFlowOp(istream_dir+'/'+idir+'/'+itflowfile)
                if holdError:
                    print(istream_dir+'/'+idir+'/'+itflowfile)
                    print('Failed to read in this file')
                    break
                tflow.append(holdtflow)
                thisop.append(holdtop)
            if not holdError: thisdata.append(thisop)
            if show_timer: thistimer.Lap()
        # thisdata = np.array(thisdata)
        if len(thisdata) != len(self.Op_cfgs.index):
            raise IOError('Failed reading (see above) caused error. fix files please.')
        ## self.Op_cfgs['Op'][istream,icfg] = [ itflow , it ]
        self.Op_cfgs.loc[:,'Op'] = pa.Series(thisdata,index=self.Op_cfgs.index)
        # self.Op = np.rollaxis(thisdata,0,thisdata.ndim)
        # ## self.Op = [ it , icfg , itflow ]
        # self.CreateOpCorr(show_timer=show_timer)
        ## self.Op = [ itflow , it , icfg  ]
        if 'random_time_list' not in self.Op_cfgs :
            self.MakeRandomTimeList()
        self.Write_Cfgs(file_type=cfg_file_type,show_timer=show_timer,CheckCfgs=CheckCfgs)

    def FlowBootAndWrite(self):
        self.FlowBootstrap()
        self.FlowWrite()


    def CreateOpCorrDummy(self):
        for tflow in self.tflowlist:
            for it_range in range(self.nranget):
                for itsum in self.tsum_list:
                    yield tflow,it_range,itsum

    def CreateChiTinv(self):
        for ictflow,tflow in enumerate(self.tflowlist):
            tsrc_data = self.Op_cfgs['Op'].apply(lambda x: np.sum(x[ictflow])**2)
            '''
            tsrc_data   pa.Series()
            rows =      stream , configuration
            '''
            yield tflow, tsrc_data

    def SublatMethod(self,DefWipe=False):
        # if PostRead: self.Op = np.rollaxis(self.Op,0,self.Op.ndim)
        ## self.Op = [ itflow , it , icfg  ]
        Oplen = self.latparams.nt
        for ictflow,tflow in enumerate(self.tflowlist):
            flowdata = self.Op_cfgs['Op'].apply(lambda x: x[ictflow])
            for it_range in range(self.nranget):

                # degenerate_props = (it_range == 0 or it_range == self.latparams.nt/2.)
                # for itsum in self.tsum_list:
                for itsum in self.tsum_list:
                    strit_range = tstr(it_range).replace('t','tr')
                    if self.FlowNewForm:
                        thisdir = self.OpDir+'/Pickle/SplitUp_NewForm/'+'/'.join([tflow,itsum,'trand'+str(self.nrandt)])
                    else:
                        thisdir = self.OpDir+'/Pickle/SplitUp/'+'/'.join([tflow,itsum,'trand'+str(self.nrandt)])
                    mkdir_p(thisdir)
                    thisfile = thisdir + '/'+self.flowname+'_'.join([tflow,itsum,strit_range])+'.py3p'
                    if (not os.path.isfile(thisfile)) or DefWipe:
                        def This_Apply_Fun(flow_val):
                            tsrc_data = []
                            for itsrc in range(Oplen):
                                sublat_data = 0.0
                                for ix,iy in itertools.product(range(it_range+1),range(it_range+1)):
                                    sublat_data += flow_val[(itsrc+ix)%Oplen]*flow_val[(itsrc+iy)%Oplen]
                                    # if not degenerate_props:
                                    #     sublat_data += flow_val[(Oplen-(itsrc+ix))%Oplen]*flow_val[(Oplen-(itsrc+iy))%Oplen]
                                tsrc_data.append(sublat_data)
                            return tsrc_data
                        tsrc_data = flowdata.apply(This_Apply_Fun)
                        ## yielded value [ itsrc]
                        yield tflow, it_range, itsum, tsrc_data,thisfile
                    else:
                        yield None,None,None,None,False


    def CreateOpCorrSumR(self,DefWipe=False):
        index_list = np.arange(self.latparams.nt)
        def GetIndicies(to,ts,tr,thislen):
            to_tr_TO_to_ts_tr   = []
            to_tr_TO_to_ts_mintr   = []
            to_TO_ts = []
            np.random.shuffle(index_list)
            to_ts = (to+ts+1)%thislen
            if to_ts <= to:
                this_ilist = np.append(index_list[to:],index_list[:to_ts])
            else:
                this_ilist = index_list[to:to_ts]
            for itr in range(tr+1):
                to_TO_ts.append(this_ilist)
                to_tr_TO_to_ts_tr.append((this_ilist+itr)%thislen)
                to_tr_TO_to_ts_mintr.append((this_ilist-itr)%thislen)
            return to_TO_ts,to_tr_TO_to_ts_tr,to_tr_TO_to_ts_mintr


        # if PostRead: self.Op = np.rollaxis(self.Op,0,self.Op.ndim)
        ## self.Op = [ itflow , it , icfg  ]
        Oplen = self.latparams.nt
        for ictflow,tflow in enumerate(self.tflowlist):
            flowdata = self.Op_cfgs['Op'].apply(lambda x: x[ictflow])
            for it_range in range(self.nranget):
                # degenerate_props = (it_range == 0 or it_range == self.latparams.nt/2.)
                for itsum in self.tsum_list:
                    strit_range = tstr(it_range).replace('t','tr')
                    its = untstr(itsum)
                    if self.FlowNewForm:
                        thisdir = self.OpDir+'/Pickle/SplitUp_NewForm/'+'/'.join([tflow,itsum,'trand'+str(self.nrandt)])
                    else:
                        thisdir = self.OpDir+'/Pickle/SplitUp/'+'/'.join([tflow,itsum,'trand'+str(self.nrandt)])
                    mkdir_p(thisdir)
                    thisfile = thisdir + '/'+self.flowname+'_'.join([tflow,itsum,strit_range])+'.py3p'
                    if (not os.path.isfile(thisfile)) or DefWipe:
                        def This_Apply_Fun(flow_val):
                            ileft,iright,irightmin = GetIndicies(0,its,it_range,Oplen)
                            ## first value is tr = 0, so dont double count for +- tr
                            flow_val = np.array(flow_val)

                            tsrc_data = np.array([np.sum(flow_val[ileft[0]]*flow_val[iright[0]])])
                            for iileft,iiright,iirightmin in zip(ileft[1:],iright[1:],irightmin[1:]):
                                tsrc_data[0] = tsrc_data[0] + np.sum(flow_val[iileft]*flow_val[iiright])
                                if iiright[0] != iirightmin[0]: tsrc_data[0] = tsrc_data[0] + np.sum(flow_val[iileft]*flow_val[iirightmin])

                            for itsrc in range(1,Oplen):
                                ileft,iright,irightmin = GetIndicies(itsrc,its,it_range,Oplen)
                                ## first value is tr = 0, so dont double count for +- tr
                                tsrc_data = np.append(tsrc_data,[np.sum(flow_val[ileft[0]]*flow_val[iright[0]])])
                                for iileft,iiright,iirightmin in zip(ileft[1:],iright[1:],irightmin[1:]):
                                    tsrc_data[itsrc] = tsrc_data[itsrc] + np.sum(flow_val[iileft]*flow_val[iiright])
                                    if iiright[0] != iirightmin[0]: tsrc_data[itsrc] = tsrc_data[itsrc] + np.sum(flow_val[iileft]*flow_val[iirightmin])
                                #iterates over [flattened itflow, it_range,itsum,itsrc]
                            return tsrc_data
                        tsrc_data = flowdata.apply(This_Apply_Fun)
                        ## yielded value [ itsrc]
                        yield tflow, it_range, itsum, tsrc_data,thisfile
                    else:
                        yield None,None,None,None,False


    def CreateOpCorr(self,DefWipe=False):
        index_list = np.arange(self.latparams.nt)
        def GetIndicies(to,ts,tr,thislen):
            np.random.shuffle(index_list)
            to_ts = (to+ts+1)%thislen
            if to_ts <= to:
                this_ilist = np.append(index_list[to:],index_list[:to_ts])
            else:
                this_ilist = index_list[to:to_ts]
            to_TO_ts = this_ilist
            to_tr_TO_to_ts_tr = (this_ilist+tr)%thislen
            to_tr_TO_to_ts_mintr = (this_ilist-tr)%thislen
            is_deg = to_tr_TO_to_ts_tr[0] == to_tr_TO_to_ts_mintr[0]
            return to_TO_ts,to_tr_TO_to_ts_tr,to_tr_TO_to_ts_mintr,is_deg


        # if PostRead: self.Op = np.rollaxis(self.Op,0,self.Op.ndim)
        ## self.Op = [ itflow , it , icfg  ]
        Oplen = self.latparams.nt
        for ictflow,tflow in enumerate(self.tflowlist):
            flowdata = self.Op_cfgs['Op'].apply(lambda x: x[ictflow])
            for it_range in range(self.nranget):
                for itsum in self.tsum_list:
                    its = untstr(itsum)
                    strit_range = tstr(it_range).replace('t','tr')
                    if self.FlowNewForm:
                        thisdir = self.OpDir+'/Pickle/SplitUp_NewForm/'+'/'.join([tflow,itsum,'trand'+str(self.nrandt)])
                    else:
                        thisdir = self.OpDir+'/Pickle/SplitUp/'+'/'.join([tflow,itsum,'trand'+str(self.nrandt)])
                    mkdir_p(thisdir)
                    thisfile = thisdir + '/'+self.flowname+'_'.join([tflow,itsum,strit_range])+'.py3p'
                    if (not os.path.isfile(thisfile)) or DefWipe:
                        def This_Apply_Fun(flow_val):
                            flow_val = np.array(flow_val)
                            ileft,iright,irightneg,itr_deg = GetIndicies(0,its,it_range,Oplen)
                            rightval = np.sum(flow_val[ileft]*flow_val[iright])
                            leftval = np.sum(flow_val[ileft]*flow_val[irightneg])
                            if itr_deg: tsrc_data = np.array([rightval])
                            else: tsrc_data = np.array([leftval+rightval])
                            for itsrc in range(1,Oplen):
                                ileft,iright,irightneg,itr_deg = GetIndicies(itsrc,its,it_range,Oplen)
                                rightval = np.sum(flow_val[ileft]*flow_val[iright])
                                leftval = np.sum(flow_val[ileft]*flow_val[irightneg])
                                if itr_deg or not self.Op_fold: tsrc_data = np.append(tsrc_data,[rightval])
                                else: tsrc_data = np.append(tsrc_data,[(leftval+rightval)/2.])
                                #iterates over [flattened itflow, it_range,itsum,itsrc]
                            return tsrc_data
                        ## yielded value [ itsrc ]
                        tsrc_data = flowdata.apply(This_Apply_Fun)
                        yield tflow, it_range, itsum, tsrc_data,thisfile
                    else:
                        yield None,None,None,None,False
        ## self.Op2 = [ itflow , it_range , itsum , ito ,  icfg ]
        # if PostRead: self.Op = np.rollaxis(self.Op,-1)

    def FlowBootO2(self,DefWipe=False,show_timer=True):
        # mkdir_p(self.OpDir+'/Pickle//temp/')
        if self.do_sublat:
            this_gen = self.SublatMethod
        else:
            if self.do_sum_range:
                this_gen = self.CreateOpCorrSumR
            else:
                this_gen = self.CreateOpCorr

        thislambda_Auto = lambda x : np.mean([x['values'][ival] for ival in x['rlist']])
        thislambda_boot = [(lambda x : x['values'][x['rlist'][icv]]) for icv in range(self.nrandt)]
        if show_timer: thistimer = Timer(   linklist=list(range(len(list(self.CreateOpCorrDummy())))),
                                            name='Booting And Autoing Op2 '+self.flowname)

        for tflows, it_range, strit_sum, tdata_sum,thisfile in this_gen(DefWipe=DefWipe):
            if thisfile is not False:
                print('FNF',thisfile, end=' ')
                strit_range = tstr(it_range).replace('t','tr')
                it_sum = untstr(strit_sum)
                if self.FlowNewForm:
                    thisdir = self.OpDir+'/Pickle/SplitUp_NewForm/'+'/'.join([tflows,strit_sum,'trand'+str(self.nrandt)])
                else:
                    thisdir = self.OpDir+'/Pickle/SplitUp/'+'/'.join([tflows,strit_sum,'trand'+str(self.nrandt)])
                mkdir_p(thisdir)
                # thisfile = thisdir + '/'+self.flowname+'_'.join([tflows,strit_sum,strit_range])+'.py3p'
                pickledata = []
                tdata = tdata_sum.to_frame(name = 'values')
                tdata['rlist'] = self.Op_cfgs['random_time_list']
                tdata_Auto = tdata.apply(thislambda_Auto,axis=1).to_frame(name='OpVals')
                if self.do_chi:
                    iTr,iTs = 1+2*it_range,it_sum+1
                    if iTr > 64:
                        iTr = iTr - 1
                    if self.do_sublat:
                        iTs = self.latparams.nt
                    # print 'DEBUG',iTr , iTs, (self.latparams.nt/float(iTr))**self.thispow
                    # thiscoeff = self.chi_coeff*(self.latparams.nt/np.sqrt(iTr*iTs))**self.thispow
                    # thiscoeff = self.chi_coeff*(self.latparams.nt/float(iTr))**self.thispow
                    thiscoeff = self.chi_coeff*(self.latparams.nt/float(iTs))**self.thispow
                    # thiscoeff = self.chi_coeff
                    thisfun_input = self.GetChitFuns(thiscoeff,self.thispow)
                else:
                    thisfun_input = self.GetChitFuns(1.,1.)
                thisAuto = AutoCorrelate(   Fun=thisfun_input,Sparam=self.Sparam,
                                            name=self.flowname + ' $' +','.join([strit_range,strit_sum,tflowTOLeg(tflows)])+'$',
                                            data=tdata_Auto)
                pickledata.append(thisAuto)
                    # pickledata[0].Power(self.thispow)
                    # pickledata[0].MultCoeff(thischi)
                    # pickledata[0].MultCoeff(self.chi_coeff)
                pickledata.append(thisAuto.Avg)
                pickledata.append(thisAuto.Std)
                rlist = None
                for it_src in range(self.nrandt):
                    # thisfile = self.OpDir+'/Pickle/temp/'+self.flowname+tflows+strit_range+'.py3p'
                    # if not os.path.isfile(thisfile):
                    # thisdict,thisdictAvg,thisdictStd = ODNested(),ODNested(),ODNested()
                    boothold = []
                    strit_src = tstr(it_src+1).replace('t','to')
                    tdata_boot = tdata.apply(thislambda_boot[it_src],axis=1)
                    thisboot = BootStrap(self.nboot,
                                        name=self.flowfilepref+'('+','.join([strit_range,strit_src,tflows])+')',
                                        cfgvals=tdata_boot,rand_list=rlist)
                    rlist = thisboot.Get_Rand_List()
                    if self.do_chi:
                        iTr,iTs = 1+2*it_range,it_sum+1
                        if iTr > 64:
                            iTr = iTr - 1
                        if self.do_sublat:
                            iTs = self.latparams.nt
                        # thiscoeff = self.chi_coeff*(self.latparams.nt/np.sqrt(iTr*iTs))**self.thispow
                        # thiscoeff = self.chi_coeff*(self.latparams.nt/float(iTr))**self.thispow
                        thiscoeff = self.chi_coeff*(self.latparams.nt/float(iTs))**self.thispow
                        # thiscoeff = self.chi_coeff
                        thisboot = thiscoeff * thisboot ** self.thispow
                        # boothold[0] = self.chi_coeff * boothold[0] ** self.thispow
                        thisboot.Stats()
                    boothold.append(thisboot)
                    boothold.append(thisboot.Avg)
                    boothold.append(thisboot.Std)
                    pickledata.append(boothold)
                WritePickle(thisfile,pickledata)
            if show_timer: thistimer.Lap()

        if show_timer: thistimer = Timer(linklist=list(self.CreateOpCorrDummy()),name='Reading Pickled Op2 '+self.flowname)
        lOp2,lOp2Avg,lOp2Std = [],[],[]
        lOp2A,lOp2AAvg,lOp2AStd = [],[],[]
        ilist,ilist_src = [],[]
        for tflows, it_range, strit_sum  in self.CreateOpCorrDummy():
            strit_range = tstr(it_range).replace('t','tr')
            if self.FlowNewForm:
                thisdir = self.OpDir+'/Pickle/SplitUp_NewForm/'+'/'.join([tflows,strit_sum,'trand'+str(self.nrandt)])
            else:
                thisdir = self.OpDir+'/Pickle/SplitUp/'+'/'.join([tflows,strit_sum,'trand'+str(self.nrandt)])
            mkdir_p(thisdir)
            thisfile = thisdir + '/'+self.flowname+'_'.join([tflows,strit_sum,strit_range])+'.py3p'
            data = iter(ReadPickleWrap(thisfile,thisShowRead=False))
            ilist.append((tflows,strit_range,strit_sum))
            lOp2A.append(next(data))
            lOp2AAvg.append(next(data))
            lOp2AStd.append(next(data))
            for it_src in range(self.nrandt):
                bootdata = next(data)
                strit_src = tstr(it_src+1).replace('t','to')
                ilist_src.append((tflows,strit_range,strit_sum,strit_src))
                lOp2.append(bootdata[0])
                lOp2Avg.append(bootdata[1])
                lOp2Std.append(bootdata[2])
        indicies = pa.MultiIndex.from_tuples(ilist,names=['flow_times','time_ranges','sum_time'])
        indicies_src = pa.MultiIndex.from_tuples(ilist_src,names=['flow_times','time_ranges','sum_time','it_src'])
        self.Op2_avg_Stats.loc[:,'Auto'] = pa.Series(lOp2A,index=indicies)
        self.Op2_avg_Stats.loc[:,'AutoAvg'] = pa.Series(lOp2AAvg,index=indicies)
        self.Op2_avg_Stats.loc[:,'AutoStd'] = pa.Series(lOp2AStd,index=indicies)
        self.Op2_Stats.loc[:,'boot'] = pa.Series(lOp2,index=indicies_src)
        self.Op2_Stats.loc[:,'Avg'] = pa.Series(lOp2Avg,index=indicies_src)
        self.Op2_Stats.loc[:,'Std'] = pa.Series(lOp2Std,index=indicies_src)


    def GetChitFuns(self,thiscoeff,thispow):
        def chit(*vals):
            return thiscoeff*(np.power(vals[0],thispow))
        def chitder(*vals):
            return [thiscoeff*thispow*(np.power(vals[0],(thispow-1)))]
        return (chit,chitder)



    def FlowBootChi(self,show_timer=True):
        if show_timer: thistimer = Timer(linklist=[self.tflowlist],name='Boot And Auto Chi '+self.flowname)
        lChi,lChiAvg,lChiStd = [],[],[]
        lChiA,lChiAAvg,lChiAStd = [],[],[]
        indicies = []
        rlist = None
        for tflows,tflowdata in self.CreateChiTinv():
            thisChi = BootStrap(self.nboot, name=self.flowfilepref+'('+','.join([tflows])+')',cfgvals=tflowdata,rand_list=rlist)
            # self.Chiboot[tflows] = BootStrap(self.nboot, name=self.flowfilepref+'('+','.join([tflows])+')',cfgvals=tflowdata)
            rlist = thisChi.Get_Rand_List()
            if self.do_chi:
                thiscoeff = self.chi_coeff
                thispow = self.thispow
                thisChi = thiscoeff*thisChi**thispow
                thisChi.Stats()
            else:
                thiscoeff = 1.
                thispow = 1.
            thisfun_input = self.GetChitFuns(thiscoeff,thispow)
            thisAuto = AutoCorrelate(Fun=thisfun_input,Sparam=self.Sparam,name=self.flowname + ' $' +','.join([tflowTOLeg(tflows)])+'$',data=tflowdata.to_frame('OpVals'))
            indicies.append(tflows)
            lChi.append(thisChi)
            lChiAvg.append(thisChi.Avg)
            lChiStd.append(thisChi.Std)
            lChiA.append(thisAuto)
            lChiAAvg.append(thisAuto.Avg)
            lChiAStd.append(thisAuto.Std)
            if show_timer: thistimer.Lap()
        indicies = pa.Index(indicies,name='Flow Times')
        self.Chi_Stats.loc[:,'boot'] = pa.Series(lChi,index=indicies)
        self.Chi_Stats.loc[:,'Avg'] = pa.Series(lChiAvg,index=indicies)
        self.Chi_Stats.loc[:,'Std'] = pa.Series(lChiStd,index=indicies)
        self.Chi_Stats.loc[:,'Auto'] = pa.Series(lChiA,index=indicies)
        self.Chi_Stats.loc[:,'AutoAvg'] = pa.Series(lChiAAvg,index=indicies)
        self.Chi_Stats.loc[:,'AutoStd'] = pa.Series(lChiAStd,index=indicies)


    def FlowBootstrap(self,DefWipe=False,show_timer=True,WipeData=True,OnlyQ=False,CheckCfgs=True):
        if 'Op' not in self.Op_cfgs or  len(self.tflowlist) != len(self.Op_cfgs['Op'].iloc[0]):
            self.FlowRead(show_timer=show_timer,CheckCfgs=CheckCfgs)
        if show_timer: thistimer = Timer(linklist=[self.tflowlist,list(range(self.nrandt))],name='Bootstrapping Op '+self.flowname)
        lOp,lOpAvg,lOpStd = [],[],[]
        lOpA,lOpAAvg,lOpAStd = [],[],[]
        indicies = []
        rlist = None
        thisfun_input = self.GetAutoFuns()
        for itcf,tflows in enumerate(self.tflowlist):
            for itsrc in range(self.nrandt):
                def temp_fun(x):
                    return x['Op'][itcf][itsrc]
                    # return x['Op'][itcf][x['random_time_list'][itsrc]]
                strit = tstr(itsrc+1).replace('t','to')
                indicies.append((tflows,strit))
                tdata = self.Op_cfgs.apply(temp_fun,axis=1).to_frame('OpVals')
                thisAuto = AutoCorrelate(   Fun=thisfun_input,Sparam=self.Sparam,
                                            name=self.flowname + ' $' +','.join([strit,tflows])+'$',data=tdata)
                tdata = self.Op_cfgs.apply(temp_fun ,axis=1)
                thisboot = BootStrap(   self.nboot, name=self.flowfilepref+'('+','.join([strit,tflows])+')',
                                        cfgvals=tdata,rand_list=rlist)
                rlist = thisboot.Get_Rand_List()
                lOp.append(thisboot)
                lOpAvg.append(thisboot.Avg)
                lOpStd.append(thisboot.Std)
                lOpA.append(thisAuto)
                lOpAAvg.append(thisAuto.Avg)
                lOpAStd.append(thisAuto.Std)
                if show_timer: thistimer.Lap()
        indicies = pa.MultiIndex.from_tuples(indicies,names=['flow_times','it_src'])
        self.Op_Stats.loc[:,'boot'] = pa.Series(lOp,index=indicies)
        self.Op_Stats.loc[:,'Avg'] = pa.Series(lOpAvg,index=indicies)
        self.Op_Stats.loc[:,'Std'] = pa.Series(lOpStd,index=indicies)
        self.Op_Stats.loc[:,'Auto'] = pa.Series(lOpA,index=indicies)
        self.Op_Stats.loc[:,'AutoAvg'] = pa.Series(lOpAAvg,index=indicies)
        self.Op_Stats.loc[:,'AutoStd'] = pa.Series(lOpAStd,index=indicies)

        if show_timer: thistimer = Timer(linklist=[self.tflowlist],name='Autocorrelating Op '+self.flowname)
        # thisfun_input = self.GetAvgAutoFuns()
        lOp,lOpAvg,lOpStd = [],[],[]
        indicies = []
        for itcf,tflows in enumerate(self.tflowlist):
            # tdata = self.Op_cfgs.apply(lambda x : np.mean([x['Op'][itcf][isrc] for isrc in x['random_time_list']]),axis=1).to_frame('OpVals')
            # tdata = self.Op_cfgs.apply(lambda x : np.mean(x['Op'][itcf]),axis=1).to_frame('OpVals')
            tdata = self.Op_cfgs.apply(lambda x : np.sum(x['Op'][itcf]),axis=1).to_frame('OpVals')
            thisAuto = AutoCorrelate(   Fun=thisfun_input,Sparam=self.Sparam,
                                        name=self.flowname + ' $' +','.join([tflowTOLeg(tflows)])+'$',data=tdata)
            indicies.append(tflows)
            lOp.append(thisAuto)
            lOpAvg.append(thisAuto.Avg)
            lOpStd.append(thisAuto.Std)
            if show_timer: thistimer.Lap()
        indicies = pa.Index(indicies,name='Flow Times')
        self.Op_avg_Stats.loc[:,'Auto'] = pa.Series(lOp,index=indicies)
        self.Op_avg_Stats.loc[:,'AutoAvg'] = pa.Series(lOpAvg,index=indicies)
        self.Op_avg_Stats.loc[:,'AutoStd'] = pa.Series(lOpStd,index=indicies)
        self.FlowBootChi(show_timer=show_timer)
        if not OnlyQ:
            self.FlowBootO2(DefWipe=DefWipe,show_timer=show_timer)
        self.Do_toAvg(show_timer=show_timer)
        if WipeData: del self.Op_cfgs['Op']

    def Do_toAvg(self,show_timer=True):
        if show_timer: thistimer = Timer(linklist=self.Op_Stats['boot'].index.get_level_values(0),name='Creating Opboot_toavg '+self.flowname)

        lOp,lOpAvg,lOpStd = [],[],[]
        indicies = []
        for itflow,flowdata in self.Op_Stats['boot'].groupby(level='flow_times'):
            indicies.append(itflow)
            # thisboot = flowdata.mean()
            thisboot = np.sum(flowdata.values)
            thisboot.Stats()
            lOp.append(thisboot)
            lOpAvg.append(thisboot.Avg)
            lOpStd.append(thisboot.Std)
            if show_timer: thistimer.Lap()
        indicies = pa.Index(indicies,name='Flow Times')
        self.Op_avg_Stats.loc[:,'boot'] = pa.Series(lOp,index=indicies)
        self.Op_avg_Stats.loc[:,'Avg'] = pa.Series(lOpAvg,index=indicies)
        self.Op_avg_Stats.loc[:,'Std'] = pa.Series(lOpStd,index=indicies)
        if show_timer: thistimer = Timer(linklist=[list(range(len(self.tflowlist))),list(range(self.nranget)),list(range(len(self.tsum_list)))],name='Creating Op2boot_toavg '+self.flowname)

        if 'boot' in self.Op2_Stats:
            lOp2,lOp2Avg,lOp2Std = [],[],[]
            indicies = []
            for (itflow,it_range,it_sum),tsumdata in self.Op2_Stats['boot'].groupby(level=['flow_times','time_ranges','sum_time']):
                indicies.append((itflow,it_range,it_sum))
                # thisboot = tsumdata.mean()
                thisboot = np.mean(tsumdata.values)
                thisboot.Stats()
                lOp2.append(thisboot)
                lOp2Avg.append(thisboot.Avg)
                lOp2Std.append(thisboot.Std)
                if show_timer: thistimer.Lap()
            indicies = pa.MultiIndex.from_tuples(indicies,names=['flow_times','time_ranges','sum_time'])
            self.Op2_avg_Stats.loc[:,'boot'] = pa.Series(lOp2,index=indicies)
            self.Op2_avg_Stats.loc[:,'Avg'] = pa.Series(lOp2Avg,index=indicies)
            self.Op2_avg_Stats.loc[:,'Std'] = pa.Series(lOp2Std,index=indicies)




    def GetAvgAutoFuns(self):
        def AvgFun(*vals):
            return vals[0]
        def AvgFunDer(*vals):
            return [1.]
        return AvgFun,AvgFunDer


    def GetAutoFuns(self):
        def AvgFun(*vals):
            return vals[0]
        def AvgFunDer(*vals):
            return [1.]
        return AvgFun,AvgFunDer

    def GetAvgChiAutoFuns(self,thiscoeff,thispow):
        def AvgFun(*vals):
            return thiscoeff*np.mean(vals)**thispow
        def AvgFunDer(*vals):
            thismean = np.mean(vals)
            return thiscoeff*thispow*thismean**(thispow-1)*np.array([((np.arange(self.nrandt) == jc).astype(float) /self.nrandt) for jc in range(self.nrandt)])
        return AvgFun,AvgFunDer






    def FlowObsToName(self,obs):
        if obs == 'TopChargexTopCharge':
            return '<Q^{2}>'
        elif obs == 'WeinbergxWeinberg':
            return '<W^{2}>'
        elif 'TopCharge' in obs and 'Weinberg' in obs:
            return '<WQ>'
        elif 'TopCharge' == obs:
            return '<Q>'
        elif 'Weinberg' == obs:
            return '<W>'
        else:
            return obs

    def FlowObsToChi(self,obs):
        if obs == 'TopChargexTopCharge':
            return r'\chi_{Q}'
        elif obs == 'WeinbergxWeinberg':
            return r'\chi_{W}'
        else:
            return r'\chi'

    def FlowSetCustomName(self,string='',stringLL=''):
        if not hasattr(self,'stream_name'):
            self.stream_name = '-def-'
        if string == '':
            self.flowname = '_'.join([self.dim_label,self.kud,self.ks,self.stream_name,self.Observ,'Full'])
            cfgfilename = self.flowname
            if self.do_chi: self.flowname +='chit'
            if self.do_sublat: self.flowname +='_sublat'
            elif self.do_sum_range: self.flowname +='_sumrange'
        else:
            self.flowname = string
            cfgfilename = self.flowname
        if stringLL == '':
            self.flowLegLab = '$'+'\ '.join([self.dim_ll,self.latparams.GetPionMassLab(),self.stream_name,self.FlowObsToName(self.Observ),'nto='+str(self.nrandt),'Full']) + '$'
            self.flowLegLab2 = '$'+'\ '.join([self.dim_ll,self.latparams.GetPionMassLab(),self.stream_name,self.FlowObsToName(self.Observ+'x'+self.Observ),'nto='+str(self.nrandt),'Full']) + '$'
            if self.do_chi:
                self.flowLegLab2 = self.flowLegLab2.replace(self.FlowObsToName(self.Observ+'x'+self.Observ),self.FlowObsToChi(self.Observ+'x'+self.Observ))
            if not self.Op_fold:
                self.flowname = self.flowname + '_NoFold'
            else:
                self.flowLegLab2 = self.flowLegLab2 + r'$\ Opfold$'
                self.flowLegLab = self.flowLegLab + r'$\ Opfold$'
            if self.do_sublat:
                self.flowLegLab2 = self.flowLegLab2 + r'$\ sublat$'
            elif self.do_sum_range:
                self.flowLegLab2 = self.flowLegLab2 + r'$\ sumtr$'
        else:
            # print 'setting Legend to ',stringLL
            self.flowLegLab = stringLL
        if self.FlowNewForm:
            self.OpDir = outputdir +'/FlowOps_NewForm/'
            self.OpCfgsDir = cfgfmtdir +'/FlowOps_NewForm/'
        else:
            self.OpDir = outputdir +'/FlowOps/'
            self.OpCfgsDir = cfgfmtdir +'/FlowOps/'
        mkdir_p(self.OpCfgsDir)
        mkdir_p(self.OpDir+'/Pickle/')
        mkdir_p(self.OpDir+'/Excel/')
        self.flowHumanFile = self.OpDir+self.flowname+'.xml'
        self.flowExcelFile = self.OpDir+'/Excel/'+self.flowname+'.xlsx'
        self.flowPickleFile = self.OpDir+'/Pickle/'+self.flowname+'.py3p'
        self.flowCfgFile = self.OpCfgsDir+'/'+cfgfilename+'.cfgs'


    def MakeStreamCfgs(self,cfglist):
        if 'stream_cfgs' not in cfglist.columns and 'configs' in cfglist:
            cfglist.loc[:,'stream_cfgs'] = pa.Series(['00'.join((istream,icfg)) for (istream,iccfg),icfg in cfglist['configs'].items()],index=cfglist.index)
        return cfglist

    def FlowImportCfgList(self,cfglist,stream_list='PreDefine'):
        # if isinstance(cfglist,(dict,OrderedDict,pa.Series)): cfglist = cfglist.keys()
        if not isinstance(stream_list,str):
            self.stream_list = stream_list
        if not isinstance(cfglist,pa.DataFrame):
            ## cfglist must be 2-d array, first dimension is stream, second dimension is configuraiton number
            ilist = []
            for istream,is_cfg in enumerate(cfglist):
                for icf in range(is_cfg):
                    ilist.append((istream,is_cfg))
            indicies = pa.MultiIndex.from_tuples(ilist,names=['stream','config_number'])
            cfglist = pa.DataFrame(np.array(cfglist).flatten(),columns=['configs'],index=indicies)
        self.Op_cfgs = self.MakeStreamCfgs(cfglist)
        self.Update_Stream_List()
        if 'configs' in self.Op_cfgs:
            self.flowncfg_list = [icfg.size for istream,icfg in self.Op_cfgs['configs'].groupby(level='stream')]

    def Update_Stream_List(self):
        if 'configs' in self.Op_cfgs:
            self.stream_list = self.Op_cfgs.index.levels[0]
            self.nstream = len(self.stream_list)
            old_sn = self.stream_name
            self.stream_name = '-'+''.join([istr.replace('-','') for istr in self.stream_list])+'-'
            if hasattr(self,'flowname'):
                self.flowname = self.flowname.replace(old_sn,self.stream_name)
            if hasattr(self,'flowLegLab'):
                self.flowLegLab = self.flowLegLab.replace(old_sn,self.stream_name)
            if hasattr(self,'flowLegLab2'):
                self.flowLegLab2 = self.flowLegLab2.replace(old_sn,self.stream_name)


    def FlowStats(self):
        if 'boot' not in self.Op_Stats: return
        for itf,it,iOp in self.Flowitems():
            iOp.Stats()
            self.OpAvg[itf][it] = iOp.Avg
            self.OpStd[itf][it] = iOp.Std


    def FlowRemoveVals(self):
        if 'Op' in self.Op_cfgs: del self.Op_cfgs['Op']


    def FlowImportPlotParams(self,thiscol,thissym,thisshift):
        self.FlowCheckCol(thiscol)
        self.FlowCheckSym(thissym)
        if thisshift != 'PreDefine':self.thisshift = thisshift


    def GetObsFun(self):
        def Ident(*val):
            return val[0]
        def IdentDer(*val):
            return [1.0]+[0.0]*(len(val)-1)
        return (Ident,IdentDer)


    def FlowCheckCol(self,thiscol):
        if 'PreDefine' == thiscol :
            if 'Not Set' == self.thiscol:
                raise IOError('Pass in color to initialize it')
        else:
            self.thiscol = thiscol

    def FlowCheckSym(self,thissym):
        if 'PreDefine' == thissym :
            if 'Not Set' == self.thissym:
                raise IOError('Pass in symbol to initialize it')
        else:
            self.thissym = thissym


    def FlowGetShift(self,xlims,thisshift):
        if thisshift != 'PreDefine': self.thisshift = thisshift
        xlen = np.abs(xlims[1]-xlims[0])
        return self.thisshift*xlen


    # def Fit_Propagator(self,t_fit_min,t_fit_max,t_sum,t_flow,fit_fun_info='PreDef'):
    #     fitr = 'fitr'+str(untstr(t_fit_min))+'-'+str(untstr(t_fit_max))
    #     if fit_fun_info != 'PreDef': self.prop_fit_fun = fit_fun_info
    #     if untstr(t_fit_max)-untstr(t_fit_min) < self.prop_fit_fun[1]:
    #         print 'Insufficient degrees of fredom for fit',fitr, 'for fit function', self.prop_fit_fun
    #         return None
    #     prop_fit_name = '_'.join([fitr,t_sum,t_flow])
    #     fit_data = self.Op2_avg_Stats['boot'].loc[t_flow,:,t_sum][t_fit_min:t_fit_max]
    #     parse_data = pa.DataFrame(columns=['ydata','xdata'],index=[0])
    #     parse_data['ydata'][0] = fit_data.values
    #     fitxdata = list(fit_data.index)
    #     if isinstance(fitxdata[0],tuple):
    #         fitxdata = [iOp2[1] for iOp2 in fitxdata]
    #     parse_data['xdata'][0] = map(untstr,fitxdata)
    #     fit_class = ff.Fitting(Funs=self.prop_fit_fun,data=parse_data,name=prop_fit_name)
    #     ikey = fit_class.name.split('_')
    #     ikey = ikey[:-2]+ [ikey[-2]+'_'+ikey[-1]]
    #     if tuple(ikey) in list(self.Prop_Fit_Stats.index):
    #         return None
    #     fit_class.FitBoots()
    #     return fit_class

    # def Fit_All_Props(  self,fit_range,sum_list='All',tflow_list='All',
    #                     fun_info='PreDef',WipeFit=True,show_timer=False):
    def Fit_All_Props(  self,fit_range,sum_list='All',tflow_list='All',
                        fun_info='PreDef',WipeFit=False,show_timer=True):
        if fun_info != 'PreDef': self.prop_fit_fun = fun_info
        if sum_list == 'All': sum_list = self.tsum_list
        if tflow_list == 'All': tflow_list = self.tflowlist
        if isinstance(fit_range[0],(list,tuple,np.ndarray)) or (isinstance(fit_range[0],str) and 'fitr' in fit_range[0]):
            print('Warning, new version with implementation of SetsOfFits scans all fit ranges')
            print('no need for list of fit_rangegs in Fit_All_Props')
            fit_range = fit_range[0]

        # if show_timer: thistimer = Timer(   linklist=[ival[0] for ival in itertools.product(sum_list,tflow_list)],name='Setting up Fits ')
        flist,index_list = [],[]
        for isum,itflow in itertools.product(sum_list,tflow_list):
            ikey = (isum,itflow)
            if tuple(ikey) in list(self.Prop_Fit_Stats.index):
                if WipeFit:
                    # print 'DEBUG'
                    # print ikey
                    # print self.Prop_Fit_Stats
                    self.Prop_Fit_Stats.drop([tuple(ikey)],inplace=True)
                else:
                    continue
            fit_out = sff.SetOfFitFuns(data=self.Op2_avg_Stats['boot'].loc[itflow,:,isum],name = '_'.join(ikey))
            susc = fit_out.ScanRange_fitr(fit_range,fit_info = {'Funs':self.prop_fit_fun})
            if susc:
                flist.append(fit_out)
                index_list.append(ikey)
            else:
                print('problem constructing setoffits for ',ikey)
        if len(index_list) > 0:
            indicies = pa.MultiIndex.from_tuples(index_list,names=self.Prop_Fit_Col_Names)
            if 'boot' in self.Prop_Fit_Stats.columns:
                this_df = pa.DataFrame()
                this_df.loc[:,'boot'] = pa.Series(flist,index=indicies)
                # this_df.loc[:,'Avg'] = pa.Series(aflist,index=indicies)
                # this_df.loc[:,'Std'] = pa.Series(sflist,index=indicies)
                self.Prop_Fit_Stats = self.Prop_Fit_Stats.append(this_df)
            else:
                self.Prop_Fit_Stats.loc[:,'boot'] = pa.Series(flist,index=indicies)
                # self.Prop_Fit_Stats.loc[:,'Avg'] = pa.Series(aflist,index=indicies)
                # self.Prop_Fit_Stats.loc[:,'Std'] = pa.Series(sflist,index=indicies)
            if show_timer: thistimer = Timer(linklist=list(range(len(self.Prop_Fit_Stats['boot']))),name='Fitting')
            def DoFit_Fun(val):
                print()
                val.DoFits(show_timer=show_timer)
                print()
                if show_timer: thistimer.Lap(val.name)
                return val
            self.Prop_Fit_Stats['boot'].apply(DoFit_Fun)
            self.FlowWrite()

### BLOCK OF PLOTTING FUNCTIONS ################################################


    def Plot_Prop_Tflow(self,plot_class,thiscol='PreDefine',thissym='PreDefine',
                        thisshift='PreDefine',tfitr='First',tsum='First',thisparam='First'):
        self.FlowImportPlotParams(thiscol,thissym,thisshift)
        if any(np.array([tfitr,tsum]) == 'First') and ('boot' not in self.Prop_Fit_Stats \
                                                        or len(self.Prop_Fit_Stats.index) == 0):
            print('No Fits found, not plotting any')
            return plot_class
        elif ('boot' not in self.Prop_Fit_Stats or len(self.Prop_Fit_Stats.index) == 0):
            print(tfitr,tsum,'Not found, attempting to fit all flow times now')
            self.Fit_All_Props(tfitr,sum_list=[tsum],show_timer=True)
        ## assuming result is the first fit parameter.
        ## keys are [itsum,itflow,fit_fun,fit_range,parameter]
        fit_data = sff.PullSOFSeries_Wpar(self.Prop_Fit_Stats.loc[:,'boot'])
        avg_data = fit_data.apply(lambda x : x.Avg)
        std_data = fit_data.apply(lambda x : x.Std)
        if tsum == 'First':
            tsum = list(fit_data.index)[0][0]
        fit_fun = fit_data.index[0][2]
        if tfitr == 'First':
            tfitr = list(fit_data.index)[0][3]
        elif isinstance(tfitr,(list,tuple,np.ndarray)) and 'fitr' in tfitr[0]:
            tfitr = tfitr[0]
        if thisparam == 'First':
            thisparam = fit_data.index[0][4]
        elif isinstance(thisparam,int):
            thisparam = fit_data.index.levels[thisparam]

        updated = False
        for itflow in self.tflowlist:
            if 'boot' not in self.Prop_Fit_Stats.columns or (tsum,itflow) not in self.Prop_Fit_Stats['boot'].index:
                updated = True
                print(tfitr,tsum,itflow , 'Not found, attempting to fit all flow times now')
                self.Fit_All_Props(tfitr,sum_list=[tsum],show_timer=True)
                if (tsum,itflow) not in self.Prop_Fit_Stats['boot'].index:
                    print(list(self.Prop_Fit_Stats['boot'].index))
                    print(tsum,itflow)
                    raise IOError('error with fitting')
        if updated:
            fit_data = sff.PullSOFSeries_Wpar(self.Prop_Fit_Stats.loc[:,'boot'])
            avg_data = fit_data.apply(lambda x : x.Avg)
            std_data = fit_data.apply(lambda x : x.Std)
        # fit_data = self.Prop_Fit_Stats.loc[(tfitr,tsum),'boot']
        this_key = (tsum,slice(None),fit_fun,tfitr,thisparam)
        # tflowplot = list(fit_data.loc[tfitr,tsum,slice(None),fit_parameter_label].index)
        # if isinstance(tflowplot[0],tuple):
        #     tflowplot = [iOp[2] for iOp in tflowplot]
        # tflowphys = np.array(map(untflowstr,tflowplot))
        # tflowphys = TflowToPhys(np.array(tflowphys),self.latparams.latspace)
        # plotshift = self.thisshift*np.abs(tflowphys[-1]-tflowphys[0])
        # print 'DEBUG'
        # print self.tflowlist
        # print self.Op2_avg_Stats
        # print self.Prop_Fit_Stats
        # for ix,iy,iyerr in zip(np.array(tflowphys),plot_data.apply(lambda x : x.Avg),plot_data.apply(lambda x : x.Std)):
        #     print ix,iy,iyerr
        hold_series = null_series
        hold_series['type'] = 'error_bar_vary'
        hold_series['label'] = 'Prop '# self.flowLegLab # fit_data.name + ' ' +fit_parameter_label
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = self.thisshift
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = avg_data
        hold_series['key_select'] = this_key
        hold_series['yerr_data'] = std_data
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        # pl.errorbar(np.array(plotx)+plotshift,ploty,plotyerr,label=self.flowLegLab2+' $'+str(flowtime)+'$ '+thistsum+' Auto',fmt=self.thissym,color=self.thiscol)
        return plot_class

    def Plot_Prop(self,plot_class,thiscol,thisshift,tfitr='First',tsum='First',tflow='First'):
        if 'boot' not in self.Prop_Fit_Stats:
            if any(np.array([tfitr,tsum,tflow]) == 'First'):
                print('No Fits found, not plotting any')
                return plot_class
            else:
                self.Fit_All_Props(tfitr,sum_list=[tsum],tflow_list=[tflow])
        this_fit = sff.PullSOFSeries(self.Prop_Fit_Stats['boot'],fmted=True)
        if any(np.array([tfitr,tsum,tflow]) == 'First') and len(this_fit.index) == 0:
            print('No Fits found, not plotting any')
            return plot_class
        if tsum == 'First':
            tsum = list(this_fit.index)[0][0]
        if tflow == 'First':
            tflow = list(this_fit.index)[0][1]
        fit_fun = list(this_fit.index)[0][2]
        if tfitr == 'First':
            tfitr = list(this_fit.index)[0][3]
        this_key = (tsum,tflow,fit_fun,tfitr)
        if this_key not in this_fit:
            print(this_key , 'Not found, attempting to fit now')
            self.Fit_All_Props(tfitr,sum_list=[tsum],tflow_list=[tflow])
            this_fit = sff.PullSOFSeries(self.Prop_Fit_Stats['boot'],fmted=True)
        if this_key not in this_fit:
            print(list(this_fit.index))
            print(this_key)
            raise IOError('error with fitting')
        # fit_plot = this_fit['boot'][tfitr,tsum,tflow]
        hold_series = null_series
        hold_series['type'] = 'fit_vary'
        hold_series['fit_class'] = this_fit
        hold_series['key_select'] = this_key
        hold_series['label'] = this_fit.name
        if self.do_sublat:
            final_result = this_fit.Eval_Function(self.latparams.nt/2)
            hold_series['label'] += r'$\chi ='+ MakeValAndErr(final_result.Avg,final_result.Std,Dec=3)+'$'
        hold_series['color'] = thiscol
        hold_series['symbol'] = 'o'
        hold_series['shift'] = thisshift
        if not self.do_sublat:
            hold_series['ShowPar'] = r'A'
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class
        # fit_plot.PlotFunction(color=thiscol,shift=thisshift,ShowPar=r'A')


    def PlotChiLine(self,plot_class,flowtime,thiscol,Vary=True):
        if flowtime in list(self.Chi_Stats['boot'].keys()):
            # mid,err = self.Chi_Stats['boot'][flowtime].Avg,self.Chi_Stats['boot'][flowtime].Std
            mid,err = self.Chi_Stats['boot'].apply(lambda x :x.Avg),self.Chi_Stats['boot'].apply(lambda x :x.Std)
            # up,down = mid+err,mid-err
            # linelab = r' $'+r' '.join([r'\chi =',MakeValAndErr(mid[flowtime],err[flowtime],Dec=3)])+r'$'
            hold_series = null_series
            hold_series['label'] = r'ACTUAL_VALUE_\chi'
            hold_series['color'] = thiscol
            hold_series['symbol'] = 'o'
            if not Vary:
                hold_series['type'] = 'hline'
                mid,err = mid[flowtime],err[flowtime]
            else:
                hold_series['type'] = 'hline_vary'
                hold_series['key_select'] = flowtime
            hold_series['y_data'] = mid
            hold_series['yerr_data'] = err
            hold_series['fmt_class'] = KeyForamtting(self.latparams)
            plot_class.AppendData(hold_series)
            # pl.axhline(mid,color=thiscol,label = linelab)
            # pl.axhspan(up,down, alpha=fillalpha, color=thiscol)
        return plot_class

    def PlotChiLineAuto(self,plot_class,flowtime,thiscol,Vary=True):
        if flowtime in list(self.Chi_Stats['Auto'].keys()):
            mid,err = self.Chi_Stats['Auto'].apply(lambda x :x.Avg),self.Chi_Stats['Auto'].apply(lambda x :x.Std)
            # up,down = mid+err,mid-err
            # linelab = r' $'+r' '.join([r'\chi =',MakeValAndErr(mid[flowtime],err[flowtime],Dec=3)])+r'$'
            hold_series = null_series
            hold_series['label'] = r'ACTUAL_VALUE_\chi'
            hold_series['color'] = thiscol
            hold_series['key_select'] = flowtime
            if not Vary:
                hold_series['type'] = 'hline'
                mid,err = mid[flowtime],err[flowtime]
            else:
                hold_series['type'] = 'hline_vary'
                hold_series['key_select'] = flowtime
            hold_series['y_data'] = mid
            hold_series['yerr_data'] = err
            hold_series['fmt_class'] = KeyForamtting(self.latparams)
            plot_class.AppendData(hold_series)
            # pl.axhline(mid,color=thiscol,label = linelab)
            # pl.axhspan(up,down, alpha=fillalpha, color=thiscol)
        return plot_class


    def PlotFOLineAuto(self,plot_class,flowtime,thiscol,Vary=True):
        if flowtime in list(self.Op_avg_Stats['Auto'].keys()):
            mid,err = self.Op_avg_Stats['Auto'].apply(lambda x : x.Avg),self.Op_avg_Stats['Auto'].apply(lambda x : x.Std)
            # up,down = mid+err,mid-err
            # linelab = r' $'+r' '.join([r'\chi =',MakeValAndErr(mid[flowtime],err[flowtime],Dec=3)])+r'$'
            hold_series = null_series
            hold_series['label'] = r'ACTUAL_VALUE_\chi'
            if not Vary:
                hold_series['type'] = 'hline'
                mid,err = mid[flowtime],err[flowtime]
            else:
                hold_series['type'] = 'hline_vary'
                hold_series['key_select'] = flowtime
            hold_series['color'] = thiscol
            hold_series['y_data'] = mid
            hold_series['yerr_data'] = err
            hold_series['fmt_class'] = KeyForamtting(self.latparams)
            plot_class.AppendData(hold_series)
            # pl.axhline(mid,color=thiscol,label = linelab)
            # pl.axhspan(up,down, alpha=fillalpha, color=thiscol)
        return plot_class

    def PlotFOLine(self,plot_class,flowtime,thiscol,Vary=True):
        if flowtime in list(self.Op_avg_Stats['boot'].keys()):
            mid,err = self.Op_avg_Stats['boot'].apply(lambda x : x.Avg),self.Op_avg_Stats['boot'].apply(lambda x : x.Std)
            # mid,err = self.Op_avg_Stats['boot'][flowtime].Avg,self.Op_avg_Stats['boot'][flowtime].Std
            # up,down = mid+err,mid-err
            # linelab = r' $'+r' '.join([r'\chi =',MakeValAndErr(mid[flowtime],err[flowtime],Dec=3)])+r'$'
            hold_series = null_series
            if not Vary:
                hold_series['type'] = 'hline'
                mid,err = mid[flowtime],err[flowtime]
            else:
                hold_series['type'] = 'hline_vary'
                hold_series['key_select'] = flowtime
            hold_series['label'] = r'ACTUAL_VALUE_\chi'
            hold_series['color'] = thiscol
            hold_series['y_data'] = mid
            hold_series['yerr_data'] = err
            hold_series['fmt_class'] = KeyForamtting(self.latparams)
            plot_class.AppendData(hold_series)
            # pl.axhline(mid,color=thiscol,label = linelab)
            # pl.axhspan(up,down, alpha=fillalpha, color=thiscol)
        return plot_class


    def PlotChiTflow(   self,plot_class,thiscol='PreDefine',
                        thissym='PreDefine',thisshift='PreDefine'):
        self.FlowImportPlotParams(thiscol,thissym,thisshift)
        plotx = list(self.Chi_Stats['boot'].keys())
        plotx = np.array(list(map(untflowstr,plotx)))
        plotx = TflowToPhys(np.array(plotx),self.latparams.latspace)
        plotshift = self.thisshift*np.abs(plotx[-1]-plotx[0])
        ploty = self.Chi_Stats['boot'].apply(lambda x : x.Avg).values
        plotyerr = self.Chi_Stats['boot'].apply(lambda x : x.Std).values
        linelab = self.flowLegLab2+r'$ \chi_{int}$'
        hold_series = null_series
        hold_series['type'] = 'error_bar'
        hold_series['fit_class'] = False
        hold_series['label'] = linelab
        hold_series['color'] = self.thiscol
        hold_series['symbol'] = self.thissym
        hold_series['shift'] = plotshift
        hold_series['x_data'] = plotx
        hold_series['y_data'] = ploty
        hold_series['yerr_data'] = plotyerr
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        # pl.axhline(mid,color=thiscol,label = linelab)
        # pl.axhspan(up,down, alpha=fillalpha, color=thiscol)
        return plot_class


    # TODO: leave these untill we implement SetsOfFits
    def PlotVstiFitr(   self,plot_class,flowtime,thistsum,fit_max,thiscol='PreDefine',
                            thissym='PreDefine',thisshift='PreDefine',fit_min_list='All',
                            thisparam='First',x_is_min=True):
        fit_max = int(fit_max)
        if isinstance(fit_min_list,(tuple,list,np.ndarray)):
            print('SetsOfFits now does not require lst for fit_min_list, defaulting to minimum value of list')
            fit_min_list = min(list(map(int,fit_min_list)))
        elif fit_min_list == 'All':
            fit_min_list = 0
        else:
            fit_min_list = int(fit_min_list)
        fitr_range = 'fitr'+str(fit_min_list)+'-'+str(fit_max)
        self.FlowImportPlotParams(thiscol,thissym,thisshift)
        if 'boot' not in self.Prop_Fit_Stats or (thistsum,flowtime) not in self.Prop_Fit_Stats['boot']:
            print(thistsum,flowtime,'has no fits, fitting now')
            self.Fit_All_Props(fitr_range,sum_list=[thistsum],tflow_list=[flowtime])
            # raise IOError('invalid flow/tsum time picked to plot')
        else:
            fit_data = self.Prop_Fit_Stats.loc[(thistsum,flowtime),'boot'].Fit_Stats
            for imin in range(fit_min_list,fit_max+1):
                for imax in range(imin+self.prop_fit_fun[1],fit_max+1):
                    if (self.prop_fit_fun[0].__name__,str(imax),str(imin)) not in fit_data.swaplevel(1,3):
                        self.Fit_All_Props(fitr_range,sum_list=[thistsum],tflow_list=[flowtime])
                        break
        fit_data = sff.PullSOFSeries_Wpar(self.Prop_Fit_Stats.loc[:,'boot'],fmted=False)
        avg_data = fit_data.apply(lambda x : x.Avg)
        std_data = fit_data.apply(lambda x : x.Std)


        if thisparam == 'First':
            thisparam = fit_data.index[0][-1]
        elif isinstance(thisparam,int):
            thisparam = fit_data.index.levels[thisparam]
        if x_is_min:
            this_key = (thistsum,flowtime,self.prop_fit_fun[0].__name__,slice(None),str(imax),thisparam)
        else:
            this_key = (thistsum,flowtime,self.prop_fit_fun[0].__name__,str(imin),slice(None),thisparam)
        hold_series = null_series
        hold_series['type'] = 'error_bar_vary'
        hold_series['label'] = 'Fits'
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = self.thisshift
        hold_series['key_select'] = this_key
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = avg_data
        hold_series['yerr_data'] = std_data
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        # pl.errorbar(np.array(plotx)+plotshift,ploty,plotyerr,label=self.flowLegLab2+' $'+str(flowtime)+'$ '+thistsum,
        #             fmt=self.thissym,color=self.thiscol)
        plot_class = self.PlotChiLine(plot_class,flowtime,self.thiscol)
        return plot_class


    def PlotVstfFitr(   self,plot_class,flowtime,thistsum,fit_min,thiscol='PreDefine',
                            thissym='PreDefine',thisshift='PreDefine',fit_max_list='All',thisparam='First'):
        return self.PlotVstiFitr(plot_class,flowtime,thistsum,fit_max_list,thiscol=thiscol,
                            thissym=thissym,thisshift=thisshift,fit_min_list=fit_min,
                            thisparam=thisparam,x_is_min=False)

    def PlotOpT(   self,plot_class,flowtime,thiscol='PreDefine',
                    thissym='PreDefine',thisshift='PreDefine'):
        if 'boot' in self.Op_Stats and (flowtime,) not in self.Op_Stats['boot']:
            print(list(self.Op_Stats['boot'].index))
            print(flowtime)
            raise IOError('invalid flow time picked to plot')
        self.FlowImportPlotParams(thiscol,thissym,thisshift)
        thisOp = self.Op_Stats.loc[:,'boot']
        Opilist = list(thisOp.loc[(flowtime,slice(None))].index)
        if isinstance(Opilist[0],tuple):
            Opilist = [iOp[1] for iOp in Opilist]
        plotx = list(map(untstr,Opilist))
        ploty,plotyerr = thisOp.apply(lambda x : x.Avg),thisOp.apply(lambda x : x.Std),
        # ploty,plotyerr = zip(*[(iOp.Avg,iOp.Std) for iOp in thisOp.values])
        plotshift = self.thisshift*np.abs(plotx[-1]-plotx[0])
        hold_series = null_series
        hold_series['type'] = 'error_bar_vary'
        hold_series['key_select'] = (flowtime,slice(None))
        hold_series['fit_class'] = False
        hold_series['label'] = self.flowLegLab2 #+' $'+str(flowtime)+'$ '
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = plotshift
        # hold_series['x_data'] = np.array(plotx)
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = ploty
        hold_series['yerr_data'] = plotyerr
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        plot_class = self.PlotFOLine(plot_class,flowtime,self.thiscol)
        # pl.errorbar(np.array(plotx)+plotshift,ploty,plotyerr,label=self.flowLegLab2+' $'+str(flowtime)+'$ '+thistsum,
        #             fmt=self.thissym,color=self.thiscol)
        return plot_class


    def PlotOpTAuto(   self,plot_class,flowtime,thiscol='PreDefine',
                    thissym='PreDefine',thisshift='PreDefine'):
        if 'boot' in self.Op_Stats and (flowtime,) not in self.Op_Stats['Auto']:
            print(list(self.Op_Stats['Auto'].index))
            print(flowtime)
            raise IOError('invalid flow time picked to plot')
        self.FlowImportPlotParams(thiscol,thissym,thisshift)
        thisOp = self.Op_Stats.loc[:,'Auto']
        Opilist = list(thisOp.loc[(flowtime,slice(None))].index)
        if isinstance(Opilist[0],tuple):
            Opilist = [iOp[1] for iOp in Opilist]
        plotx = list(map(untstr,Opilist))
        # ploty,plotyerr = zip(*[(iOp.Avg,iOp.Std) for iOp in thisOp.values])
        ploty,plotyerr = thisOp.apply(lambda x : x.Avg),thisOp.apply(lambda x : x.Std),
        plotshift = self.thisshift*np.abs(plotx[-1]-plotx[0])
        hold_series = null_series
        hold_series['type'] = 'error_bar_vary'
        hold_series['fit_class'] = False
        hold_series['key_select'] = (flowtime,slice(None))
        hold_series['label'] = self.flowLegLab2+' Auto'
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = plotshift
        hold_series['x_data'] = 'from_keys'
        # hold_series['x_data'] = np.array(plotx)
        hold_series['y_data'] = ploty
        hold_series['yerr_data'] = plotyerr
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        plot_class = self.PlotFOLineAuto(plot_class,flowtime,self.thiscol)
        # pl.errorbar(np.array(plotx)+plotshift,ploty,plotyerr,label=self.flowLegLab2+' $'+str(flowtime)+'$ '+thistsum,
        #             fmt=self.thissym,color=self.thiscol)
        return plot_class

    def CreateShiftCorr(self,shift_num=1,Force=True):
        if 'boot_shift' in self.Op2_avg_Stats.columns and not Force:
            print('Op2 shift already computed')
            return
        self.Op2_avg_Stats['boot_shift'] = self.Op2_avg_Stats['boot']
        for (itflow,itsum),iseries in self.Op2_avg_Stats['boot'].groupby(level=('flow_times','sum_time')):
            this_data = iseries.values.flatten()
            self.Op2_avg_Stats['boot_shift'].loc[itflow,:,itsum] = this_data - np.roll(this_data,-1)
        self.Op2_avg_Stats['boot_shift'].apply(lambda x : x.Stats())
        # print('DEBUG')
        # print(self.Op2_avg_Stats['boot_shift'].iloc[0].Avg)
        # print(self.Op2_avg_Stats['boot_shift'].iloc[0].Std)

    ## Plots the observable against
    def PlotOp2CorrTsink(   self,plot_class,flowtime,thistsum,thiscol='PreDefine',
                            thissym='PreDefine',thisshift='PreDefine',fitr='First',shifted=False):
        this_col = 'boot'
        if shifted:
            this_col = 'boot_shift'
            if this_col not in self.Op2_avg_Stats.columns:
                self.CreateShiftCorr()
        if this_col in self.Op2_avg_Stats and (flowtime,thistsum) not in self.Op2_avg_Stats[this_col].swaplevel(1,2):
            print(list(self.Op2_avg_Stats[this_col].swaplevel(1,2).index))
            print(flowtime, thistsum)
            raise IOError('invalid flow/tsum time picked to plot')
        self.FlowImportPlotParams(thiscol,thissym,thisshift)
        # thisOp2 = self.Op2_avg_Stats.loc[(flowtime,slice(None),thistsum),this_col]
        thisOp2 = self.Op2_avg_Stats.loc[:,this_col]
        Op2ilist = list(thisOp2.loc[(flowtime,slice(None),thistsum)].index)
        if isinstance(Op2ilist[0],tuple):
            Op2ilist = [iOp2[1] for iOp2 in Op2ilist]
        plotx = list(map(untstr,Op2ilist))
        ploty,plotyerr = thisOp2.apply(lambda x : x.Avg),thisOp2.apply(lambda x : x.Std)
        # ploty,plotyerr = zip(*[(iOp.Avg,iOp.Std) for iOp in thisOp2.values])
        plotshift = self.thisshift*np.abs(plotx[-1]-plotx[0])
        hold_series = null_series
        hold_series['type'] = 'error_bar_vary'
        hold_series['fit_class'] = False
        # hold_series['label'] = self.flowLegLab2+' $'+str(flowtime)+'$ '+thistsum
        if shifted:
            hold_series['label'] = self.flowLegLab2 + ' Csh'
        else:
            hold_series['label'] = self.flowLegLab2
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = plotshift
        hold_series['key_select'] = (flowtime,slice(None),thistsum)
        # hold_series['x_data'] = np.array(plotx)
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = ploty
        hold_series['yerr_data'] = plotyerr
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        # pl.errorbar(np.array(plotx)+plotshift,ploty,plotyerr,label=self.flowLegLab2+' $'+str(flowtime)+'$ '+thistsum,
        #             fmt=self.thissym,color=self.thiscol)
        if not shifted:
            plot_class = self.PlotChiLine(plot_class,flowtime,self.thiscol,Vary=True)
            if fitr != 'None':
                plot_class = self.Plot_Prop(plot_class,self.thiscol,plotshift,tfitr=fitr,tsum=thistsum,tflow=flowtime)
        return plot_class

    def PlotOp2CorrTsinkAuto(self,plot_class,flowtime,thistsum,thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine'):
        if 'Auto' in self.Op2_avg_Stats and (flowtime,thistsum) not in self.Op2_avg_Stats['Auto'].swaplevel(1,2):
            print(flowtime, thistsum)
            print(list(self.Op2_avg_Stats.index))
            raise IOError('invalid flow/tsum time picked to plot')
        self.FlowImportPlotParams(thiscol,thissym,thisshift)
        thisOp2 = self.Op2_avg_Stats.loc[:,'Auto']
        Op2ilist = list(thisOp2.loc[(flowtime,slice(None),thistsum)].index)
        if isinstance(Op2ilist[0],tuple):
            Op2ilist = [iOp2[1] for iOp2 in Op2ilist]
        plotx = list(map(untstr,Op2ilist))
        # ploty,plotyerr = zip(*[(iOp.Avg,iOp.Std) for iOp in thisOp2.values])
        ploty,plotyerr = thisOp2.apply(lambda x : x.Avg),thisOp2.apply(lambda x : x.Std)
        plotshift = self.thisshift*np.abs(plotx[-1]-plotx[0])
        hold_series = null_series
        hold_series['type'] = 'error_bar_vary'
        hold_series['label'] = self.flowLegLab2+' Auto'
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['key_select'] = (flowtime,slice(None),thistsum)
        hold_series['shift'] = plotshift
        hold_series['x_data'] = 'from_keys'
        # hold_series['x_data'] = np.array(plotx)
        hold_series['y_data'] = ploty
        hold_series['yerr_data'] = plotyerr
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        # pl.errorbar(np.array(plotx)+plotshift,ploty,plotyerr,label=self.flowLegLab2+' $'+str(flowtime)+'$ '+thistsum+' Auto',fmt=self.thissym,color=self.thiscol)
        plot_class = self.PlotChiLineAuto(plot_class,flowtime,self.thiscol,Vary=True)
        return plot_class

    def PlotTsinkAuto_TauIntWOpt(self,plot_class,flowtime,thistsum,thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine'):
        if 'Auto' in self.Op2_avg_Stats and (flowtime,thistsum) not in self.Op2_avg_Stats['Auto'].swaplevel(1,2):
            print(flowtime, thistsum)
            raise IOError('invalid flow/tsum time picked to plot')
        self.FlowImportPlotParams(thiscol,thissym,thisshift)
        thisOp2 = self.Op2_avg_Stats['Auto']
        Op2ilist = list(thisOp2.loc[(flowtime,slice(None),thistsum)].index)
        # Op2ilist = list(thisOp2.index)
        if isinstance(Op2ilist[0],tuple):
            Op2ilist = [iOp2[1] for iOp2 in Op2ilist]
        plotx = list(map(untstr,Op2ilist))
        ploty,plotyerr = thisOp2.apply(lambda x : x.tau),thisOp2.apply(lambda x : x.tauerr)
        # ploty,plotyerr = zip(*[(iOp.tau,iOp.tauerr) for iOp in thisOp2.values])
        plotshift = self.thisshift*np.abs(plotx[-1]-plotx[0])
        hold_series = null_series
        hold_series['type'] = 'error_bar_vary'
        hold_series['key_select'] = (flowtime,slice(None),thistsum)
        hold_series['label'] = self.flowLegLab2#+' $'+str(flowtime)+'$ '+thistsum
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = plotshift
        # hold_series['x_data'] = np.array(plotx)
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = ploty
        hold_series['yerr_data'] = plotyerr
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        # pl.errorbar(np.array(plotx)+plotshift,ploty,plotyerr,label=self.flowLegLab2+' $'+str(flowtime)+'$ '+thistsum,fmt=self.thissym,color=self.thiscol)
        return plot_class

    ## Plots the observable against
    def PlotOp2CorrTsum(self,plot_class,flowtime,thistrange,thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine'):
        if 'boot' in self.Op2_avg_Stats and (flowtime,thistrange) not in self.Op2_avg_Stats['boot']:
            print(flowtime, thistrange)
            raise IOError('invalid flow/range time picked to plot')
        self.FlowImportPlotParams(thiscol,thissym,thisshift)
        thisOp2 = self.Op2_avg_Stats['boot']
        Op2ilist = list(thisOp2.loc[(flowtime,thistrange,slice(None))].index)
        if isinstance(Op2ilist[0],tuple):
            Op2ilist = [iOp2[2] for iOp2 in Op2ilist]
        plotx = list(map(untstr,Op2ilist))
        ploty,plotyerr = thisOp2.apply(lambda x : x.Avg),thisOp2.apply(lambda x : x.Std)
        # ploty,plotyerr = zip(*[(iOp.Avg,iOp.Std) for iOp in thisOp2.values])
        plotshift = self.thisshift*np.abs(plotx[-1]-plotx[0])
        hold_series = null_series
        hold_series['type'] = 'error_bar_vary'
        hold_series['key_select'] = (flowtime,thistrange,slice(None))
        hold_series['label'] = self.flowLegLab2#+' $'+str(flowtime)+'$ '+thistrange
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = plotshift
        # hold_series['x_data'] = np.array(plotx)
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = ploty
        hold_series['yerr_data'] = plotyerr
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        # pl.errorbar(np.array(plotx)+plotshift,ploty,plotyerr,label=self.flowLegLab2+' $'+str(flowtime)+'$ '+thistrange,fmt=self.thissym,color=self.thiscol)
        plot_class = self.PlotChiLine(plot_class,flowtime,self.thiscol)
        return plot_class

    def PlotOp2CorrTsumAuto(self,plot_class,flowtime,thistrange,thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine'):
        if 'Auto' in self.Op2_avg_Stats and (flowtime,thistrange) not in self.Op2_avg_Stats['Auto']:
            print(flowtime, thistrange)
            raise IOError('invalid flow/range time picked to plot')
        self.FlowImportPlotParams(thiscol,thissym,thisshift)
        thisOp2 = self.Op2_avg_Stats['Auto']
        Op2ilist = list(thisOp2.loc[(flowtime,thistrange,slice(None))].index)
        if isinstance(Op2ilist[0],tuple):
            Op2ilist = [iOp2[2] for iOp2 in Op2ilist]
        plotx = list(map(untstr,Op2ilist))
        ploty,plotyerr = thisOp2.apply(lambda x : x.Avg),thisOp2.apply(lambda x : x.Std)
        # ploty,plotyerr = zip(*[(iOp.Avg,iOp.Std) for iOp in thisOp2.values])
        plotshift = self.thisshift*np.abs(plotx[-1]-plotx[0])
        hold_series = null_series
        hold_series['type'] = 'error_bar_vary'
        hold_series['key_select'] = (flowtime,thistrange,slice(None))
        hold_series['label'] = self.flowLegLab2+' Auto'
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = plotshift
        # hold_series['x_data'] = np.array(plotx)
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = ploty
        hold_series['yerr_data'] = plotyerr
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        # pl.errorbar(np.array(plotx)+plotshift,ploty,plotyerr,label=self.flowLegLab2+' $'+str(flowtime)+'$ '+thistrange+' Auto',fmt=self.thissym,color=self.thiscol)
        plot_class = self.PlotChiLineAuto(plot_class,flowtime,self.thiscol)
        return plot_class

    def PlotTsumAuto_TauIntWOpt(self,plot_class,flowtime,thistrange,thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine'):
        if 'Auto' in self.Op2_avg_Stats and (flowtime,thistrange) not in self.Op2_avg_Stats['Auto']:
            print(flowtime, thistrange)
            raise IOError('invalid flow/range time picked to plot')
        self.FlowImportPlotParams(thiscol,thissym,thisshift)
        thisOp2 = self.Op2_avg_Stats['Auto']
        Op2ilist = list(thisOp2.loc[(flowtime,thistrange,slice(None))].index)
        if isinstance(Op2ilist[0],tuple):
            Op2ilist = [iOp2[2] for iOp2 in Op2ilist]
        plotx = list(map(untstr,Op2ilist))
        ploty,plotyerr = thisOp2.apply(lambda x : x.tau),thisOp2.apply(lambda x : x.tauerr)
        # ploty,plotyerr = zip(*[(iOp.tau,iOp.tauerr) for iOp in thisOp2.values])
        plotshift = self.thisshift*np.abs(plotx[-1]-plotx[0])
        hold_series = null_series
        hold_series['type'] = 'error_bar_vary'
        hold_series['label'] = self.flowLegLab2#+' $'+str(flowtime)+'$ '+thistrange
        hold_series['key_select'] = (flowtime,thistrange,slice(None))
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = plotshift
        # hold_series['x_data'] = np.array(plotx)
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = ploty
        hold_series['yerr_data'] = plotyerr
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class
        # pl.errorbar(np.array(plotx)+plotshift,ploty,plotyerr,label=self.flowLegLab2+' $'+str(flowtime)+'$ '+thistrange,fmt=self.thissym,color=self.thiscol)

    ## Plots the observable against
    def PlotOp2CorrTflow(self,plot_class,thistrange,thistsum,thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine'):
        if 'boot' in self.Op2_avg_Stats and (thistsum,thistrange) not in self.Op2_avg_Stats['boot'].swaplevel(0,2):
            print(thistrange, thistsum)
            raise IOError('invalid range/sum time picked to plot')
        self.FlowImportPlotParams(thiscol,thissym,thisshift)
        thisOp2 = self.Op2_avg_Stats['boot']
        Op2ilist = list(thisOp2.loc[(slice(None),thistrange,thistsum)].index)
        if isinstance(Op2ilist[0],tuple):
            Op2ilist = [iOp2[0] for iOp2 in Op2ilist]
        plotx = list(map(untflowstr,Op2ilist))
        # ploty,plotyerr = zip(*[(iOp.Avg,iOp.Std) for iOp in thisOp2.values])
        ploty,plotyerr = thisOp2.apply(lambda x : x.Avg),thisOp2.apply(lambda x : x.Std)
        plotshift = self.thisshift*np.abs(plotx[-1]-plotx[0])
        hold_series = null_series
        hold_series['type'] = 'error_bar_vary'
        hold_series['key_select'] = (slice(None),thistrange,thistsum)
        hold_series['label'] = self.flowLegLab2#+' '+thistrange+ ' '+thistsum
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = plotshift
        hold_series['x_data'] = np.array(plotx)
        hold_series['y_data'] = ploty
        hold_series['yerr_data'] = plotyerr
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class
        # pl.errorbar(np.array(plotx)+plotshift,ploty,plotyerr,label=self.flowLegLab2+' '+thistrange+ ' '+thistsum,fmt=self.thissym,color=self.thiscol)

    def PlotOp2CorrTflowAuto(self,plot_class,thistrange,thistsum,thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine'):
        if 'Auto' in self.Op2_avg_Stats and (thistsum,thistrange) not in self.Op2_avg_Stats['Auto'].swaplevel(0,2):
            print(thistrange, thistsum)
            raise IOError('invalid range/sum time picked to plot')
        self.FlowImportPlotParams(thiscol,thissym,thisshift)
        thisOp2 = self.Op2_avg_Stats['Auto']
        Op2ilist = list(thisOp2.loc[(slice(None),thistrange,thistsum)].index)
        if isinstance(Op2ilist[0],tuple):
            Op2ilist = [iOp2[0] for iOp2 in Op2ilist]
        plotx = list(map(untflowstr,Op2ilist))
        ploty,plotyerr = thisOp2.apply(lambda x : x.Avg),thisOp2.apply(lambda x : x.Std)
        # ploty,plotyerr = zip(*[(iOp.Avg,iOp.Std) for iOp in thisOp2.values])
        plotshift = self.thisshift*np.abs(plotx[-1]-plotx[0])
        hold_series = null_series
        hold_series['type'] = 'error_bar_vary'
        hold_series['label'] = self.flowLegLab2+' Auto'
        hold_series['key_select'] = (slice(None),thistrange,thistsum)
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = plotshift
        hold_series['x_data'] = np.array(plotx)
        hold_series['y_data'] = ploty
        hold_series['yerr_data'] = plotyerr
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class
        # pl.errorbar(np.array(plotx)+plotshift,ploty,plotyerr,label=self.flowLegLab2+' '+thistrange+ ' '+thistsum+' Auto',fmt=self.thissym,color=self.thiscol)

    def PlotTflowAuto_TauIntWOpt(self,plot_class,thistrange,thistsum,thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine'):
        if 'Auto' in self.Op2_avg_Stats and (thistsum,thistrange) not in self.Op2_avg_Stats['Auto'].swaplevel(0,2):
            print(thistrange, thistsum)
            raise IOError('invalid range/sum time picked to plot')
        self.FlowImportPlotParams(thiscol,thissym,thisshift)
        thisOp2 = self.Op2_avg_Stats['Auto']
        Op2ilist = list(thisOp2.loc[(slice(None),thistrange,thistsum)].index)
        if isinstance(Op2ilist[0],tuple):
            Op2ilist = [iOp2[0] for iOp2 in Op2ilist]
        plotx = list(map(untflowstr,Op2ilist))
        ploty,plotyerr = thisOp2.apply(lambda x : x.tau),thisOp2.apply(lambda x : x.tauerr)
        # ploty,plotyerr = zip(*[(iOp.tau,iOp.tauerr) for iOp in thisOp2.values])
        plotshift = self.thisshift*np.abs(plotx[-1]-plotx[0])
        hold_series = null_series
        hold_series['type'] = 'error_bar_vary'
        hold_series['key_select'] = (slice(None),thistrange,thistsum)
        hold_series['label'] = self.flowLegLab2+' '+thistrange+ ' '+thistsum
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = plotshift
        hold_series['x_data'] = np.array(plotx)
        hold_series['y_data'] = ploty
        hold_series['yerr_data'] = plotyerr
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class
        # pl.errorbar(np.array(plotx)+plotshift,ploty,plotyerr,label=self.flowLegLab2+' '+thistrange+ ' '+thistsum,fmt=self.thissym,color=self.thiscol)


### BLOCK OF PLOTTING FUNCTIONS ################################################








    ## Comparisons

    def __str__(self):
        return self.flowname


    def __iter__(self):
        return self

    def Flowitems(self):
        outlist = []
        for ikey,data in self.Op_Stats['boot'].items():
            if isinstance(ikey,tuple):
                outlist.append(ikey+(data,))
            else:
                outlist.append((ikey,data))
        return outlist

    def FlowitemsAuto(self):
        outlist = []
        for ikey,data in self.Op_avg_Stats['Auto'].items():
            if isinstance(ikey,tuple):
                outlist.append(ikey+(data,))
            else:
                outlist.append((ikey,data))
        return outlist

    def Flowitems2(self):
        outlist = []
        for ikey,data in self.Op2_Stats['boot'].items():
            if isinstance(ikey,tuple):
                outlist.append(ikey+(data,))
            else:
                outlist.append((ikey,data))
        return outlist

    def Flowitems2Auto(self):
        outlist = []
        for ikey,data in self.Op2_avg_Stats['Auto'].items():
            if isinstance(ikey,tuple):
                outlist.append(ikey+(data,))
            else:
                outlist.append((ikey,data))
        return outlist

    def Flowitems_toavg(self):
        outlist = []
        for ikey,data in self.Op_avg_Stats['boot'].items():
            if isinstance(ikey,tuple):
                outlist.append(ikey+(data,))
            else:
                outlist.append((ikey,data))
        return outlist


    def Flowitems2_toavg(self):
        outlist = []
        for ikey,data in self.Op2_avg_Stats['boot'].items():
            if isinstance(ikey,tuple):
                outlist.append(ikey+(data,))
            else:
                outlist.append((ikey,data))
        return outlist


    def FlowitemsAvgStd(self):
        outlist = []
        for ikey,data in self.Op_Stats['boot'].items():
            if isinstance(ikey,tuple):
                outlist.append(ikey+(data.Avg,data.Std))
            else:
                outlist.append((ikey,data.Avg,data.Std))
        return outlist


    def Flowvalues(self):
        return self.Op_Stats['boot'].values

    def FlowvaluesAuto(self):
        return self.Op_avg_Stats['Auto'].values

    def FlowvaluesAvgStd(self):
        return list(zip(self.Op_Stats['Avg'].values,self.Op_Stats['Std']))

    def Flowkeys(self):
        return self.Op_Stats['boot'].index


    def __setitem__(self, key, value ):
        self.Op_Stats.loc[key,'boot'] = value

    def __getitem__(self, key):
        return self.Op_Stats.loc[key,'boot']

    def __reversed__(self):
        self.Op_Stats.iiloc[::-1]
        self.Op_avg_Stats.iiloc[::-1]
        self.Op2_Stats.iiloc[::-1]
        self.Op2_avg_Stats.iiloc[::-1]
        self.Chi_Stats.iiloc[::-1]
        return self



    def Flownext(self):
        if self.flowcurr >= len(self.Flowkeys()):
            self.flowcurr = 0
            raise StopIteration
        else:
            self.flowcurr += 1
            if self.flowcurr >= len(self.Flowkeys())+1:
                return 0.0
            else:
                thiskey = self.Flowkeys()[self.flowcurr-1]
                return self[thiskey]

    ## incase next gets overloaded
    def __next__(self):
        return self.Flownext(self)

    def FlowWrite(self):
        def FixDictArray(thislist,ikey):
            outDict = OrderedDict()
            for ic,ilist in enumerate(thislist):
                outDict[ikey+str(ic+1)] = ilist
            return outDict

        outDict = ODNested()
        if not hasattr(self,'show_cfgs'):
            self.show_cfgs = False

        outDict['name'] = self.flowname
        outDict['kud'] = self.kud
        outDict['ks'] = self.ks
        outDict['nt'] = self.nt
        outDict['nxyz'] = self.nxyz

        for istream,incfg in zip(self.stream_list,self.flowncfg_list):
            outDict['ncfg_'+istream] = incfg
        outDict['nboot'] = self.nboot

        excel_params = pa.Series(deepcopy(outDict))

        ## uncomment this if you want to actually see every time, rather than random sums
        if 'boot' in self.Prop_Fit_Stats.columns:
            for (itsum,itflow),idata in self.Prop_Fit_Stats['boot'].items():
                outDict['Prop_Fit'][itflow][itsum] = idata.GetOutputDict()

        if 'boot' in self.Op_Stats.columns:
            for itflow,it,tdata in self.Flowitems():
                if hasattr(tdata,'Avg'):
                    outDict['Opboot'][itflow][it] = AvgStdToFormat(tdata.Avg,tdata.Std,frmtflag='e')
                else:
                    outDict['Opboot'][itflow][it] = tdata
        if 'boot' in self.Op_avg_Stats.columns:
            for itflow,tdata in self.Flowitems_toavg():
                outDict['Opboot_toavg'][itflow] = AvgStdToFormat(tdata.Avg,tdata.Std,frmtflag='e')
        if 'Auto' in self.Op_avg_Stats.columns:
            outDict['S_param'] = self.Sparam
            for itflow,tdata in self.FlowitemsAuto():
                outDict['OpAuto'][itflow] = AvgStdToFormat(tdata.Avg,tdata.Std,frmtflag='e')
                outDict['OpAuto_Wopt'][itflow] = tdata.Wopt
                outDict['OpAuto_tau'][itflow] = AvgStdToFormat(tdata.tau,tdata.tauerr,frmtflag='e')
        # if len(self.Op2boot.keys()) > 0:
        #     for itflow,itrange,itsum,it,tdata in self.Flowitems2():
        #         outDict['Op2boot'][itflow][itrange][itsum][it] = AvgStdToFormat(tdata.Avg,tdata.Std,frmtflag='e')
        if 'boot' in self.Chi_Stats.columns:
            for itflow,tdata in self.Chi_Stats['boot'].items():
                outDict['Chi_boot'][itflow] = AvgStdToFormat(tdata.Avg,tdata.Std,frmtflag='e')
        if 'Auto' in self.Chi_Stats.columns:
            for itflow,tdata in self.Chi_Stats['Auto'].items():
                outDict['Chi_Auto'][itflow] = AvgStdToFormat(tdata.Avg,tdata.Std,frmtflag='e')
        if 'boot' in self.Op2_avg_Stats.columns:
            for itflow,itrange,itsum,tdata in self.Flowitems2_toavg():
                outDict['Op2boot'][itflow][itrange][itsum] = AvgStdToFormat(tdata.Avg,tdata.Std,frmtflag='e')
        if 'Auto' in self.Op2_avg_Stats.columns:
            outDict['S_param'] = self.Sparam
            for itflow,itrange,itsum,tdata in self.Flowitems2Auto():
                outDict['Op2Auto'][itflow][itrange][itsum] = AvgStdToFormat(tdata.Avg,tdata.Std,frmtflag='e')
                outDict['Op2Auto_Wopt'][itflow][itrange][itsum] = tdata.Wopt
                outDict['Op2Auto_tau'][itflow][itrange][itsum] = AvgStdToFormat(tdata.tau,tdata.tauerr,frmtflag='e')


        if self.show_cfgs and 'configs' in self.Op_cfgs.columns:
            outDict['cfglist'] = Series_TO_ODict(self.Op_cfgs['configs'])
            for (istream,icfg),ivals in self.Op_cfgs.iterrows():
                for itsrc in range(self.nrandt):
                    thistsrc = 't_src'+str(itsrc)
                    thiscfg = ivals['configs']
                    outDict['random_sources'][istream][thiscfg][thistsrc] = ivals['random_time_list']



        ## Write human readable file with data
        WriteXml(self.flowHumanFile,{'Results':outDict})

        stats_dict = OrderedDict()
        keylist = ['Op_Stats','Op2_Stats','Op_avg_Stats','Op2_avg_Stats','Chi_Stats']
        for ikey in keylist:
            stats_dict[ikey] = deepcopy(getattr(self,ikey))

        if 'configs' in self.Op_cfgs.columns:
            WriteExcel(self.flowExcelFile,stats_dict,flowcfglist=self.Op_cfgs['configs'],params=excel_params)

            outDict['cfglist'] = Series_TO_ODict(self.Op_cfgs['configs'])

        ## pickles rest of data for reading
        # WritePickle(self.flowPickleFile,self.RemoveOp2(deepcopy(self.__dict__)))
        self.RemoveFuns()
        WritePickle(self.flowPickleFile,self.__dict__)
        self.GetFuns()

    ##################################

    def RemoveOp2(self,thisdict):
        del thisdict['Op2_Stats']
        del thisdict['Op2_avg_Stats']
        return thisdict


    ## Operator overloading

    #Too Much work, only works for one operator. Please use SetCustomName after doing all algebraic manipulation

    ## this could to with some improving.... or some more generallity somewhere

    # def FlowUpdateName(self,opp,Op,Op2):
    #     def GetOppSelf(opp):
    #         if opp == '+':
    #             return 'x2'
    #         elif opp == '-':
    #             return 'x0'
    #         elif opp == 'x':
    #             return 'pow2'
    #         elif opp == 'div':
    #             return 'divself'
    #         elif opp == 'pow':
    #             return 'powself'
    #     def GetLegSelf(opp):
    #         if opp == '+':
    #             return '\times 2'
    #         elif opp == '-':
    #             return '\times 0'
    #         elif opp == 'x':
    #             return '^{2}'
    #         elif opp == 'div':
    #             return 'divself'
    #         elif opp == 'pow':
    #             return 'powself'

    #     def LLFormatOp(opp):
    #         oppout = opp.replace('x',r'\times ')
    #         oppout = oppout.replace('pow','^{}')
    #         return oppout

    #     outname = []
    #     outLL = []
    #     oppLL = LLFormatOp(opp) ## \times looks better :)
    #     if isinstance(Op2,FlowOp) and isinstance(Op,FlowOp):
    #         op1LL = Op.flowLegLab[1:-1] ## removes $ symbols, re adds at end
    #         op2LL = Op2.flowLegLab[1:-1] ## removes $ symbols, re adds at end
    #         if Op.flowname == Op2.flowname:
    #             outname = Op.flowname,GetOppSelf(opp)
    #             outLL = op1LL,' Comb ',GetLegSelf(oppLL)
    #         else:
    #             for iname,iname2 in zip(Op.flowname.split('_'),Op2.flowname.split('_')):
    #                 if iname in iname2:
    #                     outname.append(iname2)
    #                 elif iname2 in iname:
    #                     outname.append(iname)
    #                 else:
    #                     outname.append(iname+opp+iname2)
    #             outname = '_'.join(outname)
    #             for iflowLegLab,iflowLegLab2 in zip(op1LL.split('_'),op2LL.split('\ ')):
    #                 if iflowLegLab in iflowLegLab2:
    #                     outLL.append(iflowLegLab2)
    #                 elif iflowLegLab2 in iflowLegLab:
    #                     outLL.append(iflowLegLab)
    #                 else:
    #                     outLL.append(iflowLegLab+oppLL+iflowLegLab2)
    #             outLL = '_'.join(outLL)
    #     elif isinstance(Op,FlowOp):
    #         op1LL = Op.flowLegLab[1:-1] ## removes $ symbols, re adds at end
    #         if isinstance(Op2,int):
    #             outname  =  Op.flowname,opp+str(Op2)
    #             outLL  =  op1LL,oppLL+str(Op2)
    #         else:
    #             if Op2 > 0.01:
    #                 outname = Op.flowname,opp+'{:.2f}'.format(Op2)
    #                 outLL = op1LL,oppLL+'{:.2f}'.format(Op2)
    #     else:
    #         op2LL = Op2.flowLegLab[1:-1] ## removes $ symbols, re adds at end
    #         if isinstance(Op,int):
    #             outname  =  str(Op) + opp ,Op2.flowname
    #             outLL  =  str(Op) + oppLL,op2LL
    #         else:
    #             if Op > 0.01:
    #                 outname = '{:.2f}'.format(Op) + opp,Op2.flowname
    #                 outLL = '{:.2f}'.format(Op) + oppLL,op2LL
    #     # print r'$'+r'\ '.join(outLL)+'$'
    #     outLL = list(outLL)
    #     for ic,iLL in enumerate(outLL):
    #         if '^{}' in iLL:
    #             outLL[ic] = iLL.replace('^{}','^{')+'}'
    #     self.FlowSetCustomName('_'.join(outname),r'$'+r'\ '.join(outLL)+'$')


    # ## Numerical operatoions

    # def __add__(self,Op2):
    #     if any([ipar in str(type(Op2)) for ipar in FlowOp.parentlist]):
    #         return Op2.__radd__(self)
    #     else:
    #         result = FlowOp(thisnboot=self.nboot, cfglist=self.flowcfglist,
    #                               Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
    #         result.FlowUpdateName('+',self,Op2)
    #         if isinstance(Op2,FlowOp):
    #             for (it,iOp) in self.Flowitems():
    #                 result[it] = iOp + Op2[it]
    #         else:
    #             try:
    #                 for it,iOp in self.Flowitems():
    #                     result[it] = iOp + Op2
    #             except:
    #                 print type(Op2)
    #                 raise EnvironmentError('Invalid value to combine with FlowOp class')
    #         return result

    # def __sub__(self,Op2):
    #     if any([ipar in str(type(Op2)) for ipar in FlowOp.parentlist]):
    #         return Op2.__rsub__(self)
    #     else:
    #         result = FlowOp(thisnboot=self.nboot, cfglist=self.flowcfglist,
    #                               Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
    #         result.FlowUpdateName('-',self,Op2)
    #         if isinstance(Op2,FlowOp):
    #             for (it,iOp) in self.Flowitems():
    #                 result[it] = iOp - Op2[it]
    #         else:
    #             try:
    #                 for it,iOp in self.Flowitems():
    #                     result[it] = iOp - Op2
    #             except:
    #                 print type(Op2)
    #                 raise EnvironmentError('Invalid value to combine with FlowOp class')
    #         return result

    # def __mul__(self,Op2):
    #     if any([ipar in str(type(Op2)) for ipar in FlowOp.parentlist]):
    #         return Op2.__rmul__(self)
    #     else:
    #         result = FlowOp(thisnboot=self.nboot, cfglist=self.flowcfglist,
    #                               Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
    #         result.FlowUpdateName('x',self,Op2)
    #         if isinstance(Op2,FlowOp):
    #             for (it,iOp) in self.Flowitems():
    #                 result[it] = iOp * Op2[it]
    #         else:
    #             try:
    #                 for it,iOp in self.Flowitems():
    #                     result[it] = iOp * Op2
    #             except:
    #                 print type(Op2)
    #                 raise EnvironmentError('Invalid value to combine with FlowOp class')
    #         return result

    # def __div__(self,Op2):
    #     if any([ipar in str(type(Op2)) for ipar in FlowOp.parentlist]):
    #         return Op2.__rdiv__(self)
    #     else:
    #         result = FlowOp(thisnboot=self.nboot, cfglist=self.flowcfglist,
    #                               Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
    #         result.FlowUpdateName('div',self,Op2)
    #         if isinstance(Op2,FlowOp):
    #             for (it,iOp) in self.Flowitems():
    #                 try:
    #                     result[it] = iOp / Op2[it]
    #                 except:
    #                     if any([ibs == 0 for ibs in Op2[it].bootvals]):
    #                         raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+Op2.flowname + ' '+Op2[it].flowname)

    #         else:
    #             if Op2 == 0:
    #                 raise ZeroDivisionError('Dividing by zero found in bootstrap division')
    #             try:
    #                 for it,iOp in self.Flowitems():
    #                     result[it] = iOp / Op2
    #             except:
    #                 print type(Op2)
    #                 raise EnvironmentError('Invalid value to combine with FlowOp class')
    #         return result


    # def __pow__(self,Op2):
    #     if any([ipar in str(type(Op2)) for ipar in FlowOp.parentlist]):
    #         return Op2.__rpow__(self)
    #     else:
    #         result = FlowOp(thisnboot=self.nboot, cfglist=self.flowcfglist,
    #                               Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
    #         result.FlowUpdateName('pow',self,Op2)
    #         if isinstance(Op2,FlowOp):
    #             for (it,iOp) in self.Flowitems():
    #                 result[it] = iOp ** Op2[it]
    #         else:
    #             try:
    #                 for it,iOp in self.Flowitems():
    #                     result[it] = iOp ** Op2
    #             except:
    #                 print type(Op2)
    #                 raise EnvironmentError('Invalid value to combine with FlowOp class')
    #         return result


    # ## Right multiplication functions


    # def __radd__(self,Op2):
    #     result = FlowOp(thisnboot=self.nboot, cfglist=self.flowcfglist,
    #                           Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
    #     result.FlowUpdateName('+',Op2,self)
    #     if isinstance(Op2,FlowOp):
    #         for (it,iOp) in self.Flowitems():
    #             result[it] = Op2[it] + iOp
    #     else:
    #         try:
    #             for it,iOp in self.Flowitems():
    #                 result[it] = Op2 + iOp
    #         except:
    #             print type(Op2)
    #             raise EnvironmentError('Invalid value to combine with FlowOp class')
    #     return result

    # def __rsub__(self,Op2):
    #     result = FlowOp(thisnboot=self.nboot, cfglist=self.flowcfglist,
    #                           Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
    #     result.FlowUpdateName('-',Op2,self)
    #     if isinstance(Op2,FlowOp):
    #         for (it,iOp) in self.Flowitems():
    #             result[it] =  Op2[it] - iOp
    #     else:
    #         try:
    #             for it,iOp in self.Flowitems():
    #                 result[it] =  Op2 - iOp
    #         except:
    #             print type(Op2)
    #             raise EnvironmentError('Invalid value to combine with FlowOp class')
    #     return result

    # def __rmul__(self,Op2):
    #     result = FlowOp(thisnboot=self.nboot, cfglist=self.flowcfglist,
    #                           Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
    #     result.FlowUpdateName('x',Op2,self)
    #     if isinstance(Op2,FlowOp):
    #         for (it,iOp) in self.Flowitems():
    #             result[it] =  Op2[it] * iOp
    #     else:
    #         try:
    #             for it,iOp in self.Flowitems():
    #                 result[it] = Op2 * iOp
    #         except:
    #             print type(Op2)
    #             raise EnvironmentError('Invalid value to combine with FlowOp class')
    #     return result

    # def __rdiv__(self,Op2):
    #     result = FlowOp(thisnboot=self.nboot, cfglist=self.flowcfglist,
    #                           Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
    #     result.FlowUpdateName('div',Op2,self)
    #     if isinstance(Op2,FlowOp):
    #         for (it,iOp) in self.Flowitems():
    #             try:
    #                 result[it] = Op2[it] / iOp
    #             except:
    #                 if any([ibs == 0 for ibs in iOp.bootvals]):
    #                     raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+self.flowname + ' '+iOp.flowname)

    #     else:
    #         for it,iOp in self.Flowitems():
    #             try:
    #                 result[it] = Op2 / iOp
    #             except:
    #                 if iOp == 0:
    #                     raise ZeroDivisionError('Dividing by zero found in bootstrap division')
    #                 else:
    #                     print type(Op2)
    #                     raise EnvironmentError('Invalid value to combine with FlowOp class')
    #     return result


    # def __rpow__(self,Op2):
    #     result = FlowOp(thisnboot=self.nboot, cfglist=self.flowcfglist,
    #                           Info=self.flowInfo,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
    #     result.FlowUpdateName('pow',Op2,self)
    #     if isinstance(Op2,FlowOp):
    #         for (it,iOp) in self.Flowitems():
    #             result[it] =  Op2[it] ** iOp
    #     else:
    #         try:
    #             for it,iOp in self.Flowitems():
    #                 result[it] = Op2 ** iOp
    #         except:
    #             print type(Op2)
    #             raise EnvironmentError('Invalid value to combine with FlowOp class')
    #     return result


    ## unitary arithmetic operations

    def __neg__(self):
        for itflow,it in self.Flowkeys():
            self[(itflow,it)] = -self[(itflow,it)]
        return self


    def __pos__(self):
        return self

    def __abs__(self):
        for itflow,it in self.Flowkeys():
            self[(itflow,it)] = abs(self[(itflow,it)])
        return self

    ## Casting, just casts the bootstrap values, ignore return value (only needed since return has to be same type).


    def __complex__(self):
        for itflow,it in self.Flowkeys():
            complex(self[(itflow,it)])
        self.FlowStats()
        return 1+0j

    def __int__(self):
        for itflow,it in self.Flowkeys():
            int(self[(itflow,it)])
        self.FlowStats()
        return 1

    def __float__(self):
        for itflow,it in self.Flowkeys():
            float(self[(itflow,it)])
        self.FlowStats()
        return 1.0



def TestFO(DefWipe=False):
    from Params import defInfoFlow
    defInfoFlow['tflowlist'] = np.arange(0.01,10,0.1)
    defInfoFlow['FlowNewForm'] = True
    # defInfoFlow['cfg_file_type'] = 'html'
    data = FlowOp(Info=defInfoFlow)
    data.FlowLoadPickle(DefWipe=DefWipe)
    data.FlowWrite()
    this_info = pa.Series()
    this_info['save_file'] = './TestGraphs/flow_plot_test.pdf'
    this_info['title'] = 'Test FlowOps'
    import PlotData as jpl
    data_plot = jpl.Plotting(plot_info=this_info)
    # data_plot.LoadPickle(DefWipe=False)
    # data_plot = data.PlotMonteTime(data_plot,'t_f6.01')
    # data_plot = data.PlotMonteTimeDeltaEps(data_plot,'t_f6.01',thiseps=1,thiscol='Red')
    # data_plot = data.PlotTauIntWOpt(data_plot)
    data_plot = data.FlowPlot(data_plot)
    datacomb = data.FlowMakeChit(data)
    data_plot = datacomb.FlowPlot(data_plot,thiscol='Blue',thissym='o')
    data_plot = datacomb.FlowFitPlot(data_plot,'fitr4-9')
    datacomb.FlowWrite()
    data_plot.PrintData()
    data_plot.PlotAll()
    return data,datacomb

def TestFOFull(DefWipe=False):
    from Params import defInfoFlow
    # defInfoFlow['tflowlist'] = np.array([4.01,5.01,6.01])
    # defInfoFlow['tsum_fit_list'] = ['ts28','ts30','ts32']
    # defInfoFlow['tsum_list'] = ['ts28','ts30','ts32']
    defInfoFlow['tflowlist'] = np.array([6.01])
    defInfoFlow['tsum_fit_list'] = ['ts32']
    defInfoFlow['tsum_list'] = ['ts32']
    defInfoFlow['nrandt'] = 1
    defInfoFlow['FlowNewForm'] = True
    dataFull = FlowOpFullSquared(Info=defInfoFlow)
    dataFull.FlowLoadPickle(DefWipe=DefWipe)
    # dataFull.ExpFitTime(thistflow='t_f6.01',funDebug=True)
    this_info = pa.Series()
    this_info['save_file'] = './TestGraphs/flowfull_plot_test.pdf'
    this_info['title'] = 'Test FlowOpsFull'
    import PlotData as jpl
    data_plot = jpl.Plotting(plot_info=this_info)
    dataFull.Fit_All_Props('fitr4-15')
    data_plot = dataFull.Plot_Prop_Tflow(data_plot,thiscol='Blue',thissym='o',thisshift=0.0)
    # data_plot = dataFull.PlotOp2CorrTsink(data_plot,'t_f6.01','ts32',thiscol='Blue',thissym='o',thisshift=0.0)
    data_plot.PrintData()
    data_plot.PlotAll()

    print('Look at data and datacomb and dataFull')
    return dataFull


## runing this file will do a test run!
if __name__ == '__main__':
    # %matplotlib inline
    data,datacomb = TestFO()
    dataFull = TestFOFull()
