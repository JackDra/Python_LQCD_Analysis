#!/usr/bin/env python


#IMPORT THIS FIRST, put in base import file
# import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
# import matplotlib.pyplot as pl
##


from Params import nboot,defxlimAlpha,defxlimOp,datadir,cfgfmtdir,Debug_Mode,this_dir
from FileIO import ReadFuns,WriteFuns,Construct_File_Object
from Params import cfundir,defxlim,outputdir,TflowToPhys,Qeps,cfg_file_type,defSparam,Wipe_All_Fits
from MiscFuns import ODNested,mkdir_p,CheckClass,GetInfoList,RemoveAllBoots,fmt_file_type
from MiscFuns import Series_fix_key,CombineCfgs
from MiscFuns import CreateNameAndLL,GenSinkFun,flip_fun_2arg,Series_TO_ODict
from XmlFormatting import tstr,unxmlfitr,untstr,AvgStdToFormat,epsstr,unxmlfitr_int
from XmlFormatting import tflowTOLeg,tflowstr,tsinkTOLeg,untflowstr,tsumstr,KeyForamtting
# from XmlFormatting import MakeValAndErr
from PredefFitFuns import C2OneStateFitFunNoExp,C2OneStateFitFun,C2TwoStateFitFun,C2TwoStateFitFunNoExp,C2OSFAntiper
from PredefFitFuns import ConstantFitFun,Alpha_Prop_Fit_plin,C2OSFAntiperDer
from Autocorr import AutoCorrelate
from NullPlotData import null_series
from warnings import warn
import SetsOfFits as sff
# import traceback


import MomParams as mp
from ReadBinaryCfuns import R2CChromaXMLFileList,RC2Full
from BootStrapping import BootStrap
from FileIO import WriteXml,WritePickle,ReadPickleWrap,WriteExcel,ReadWithMeta,WriteWithMeta
from copy import deepcopy
# import cPickle as pickle
from TimeStuff import Timer

import numpy as np
import pandas as pa
import operator as op
import re,glob,os
import FitFunctions as ff
from functools import partial
import itertools

from collections import OrderedDict
from functools import partial

from FlowOpps import FlowOp,FlowOpFullSquared


"""
How to use:

TwoPointCorr(cfglist=[], Info={})
creates class using:
 configuration list cfglist
 Information Dictionary Info (See below for what can be in Info)


.LoadPickle()
Does everything basically.
Will ReadAndWrite() if pickle file is not found
i.e. will run analysis if not already done before.


NNQCorr class is below, (search for string! this is big file)

WARNING: only 1 stream is implemented now. only pass in config lists or have directories
         with 1 steam (-a- for example)

"""

def GetInterpFlag(thisInterp):
    if thisInterp == 'P4' or thisInterp == 'CPEven':
        return '9'
    elif thisInterp == 'g5P4' or thisInterp == 'CPOdd':
        return '17'
    elif thisInterp == 'g5' or thisInterp == 'g5div2':
        return '1'
    elif thisInterp == 'P4g5':
        return '0'
    elif thisInterp == 'Pion':
        return '15'
    else:
        return thisInterp



class TwoPointCorr(object):
    """ C2(p,t) uses bootstrap class
        p =  momentum
        t = time
    """
    def Construct_2pt_File(this_file):
        return Construct_File_Object(this_file,TwoPointCorr)


    ## Info is a dictionary containing information for the correlator:
    ## see comments below to see what is required.
    ## missing elements just ignores that field

    ## TODO: create variational flag and interface with this function.

    ## array of all classes that inherit this class. Append this when adding new classes please!
    parentlist = ['SetOfCorrs.SetOfTwoPt','TwoPtCorrelators.NNQCorr','TwoPtCorrelators.NNQFullCorr','VarMethod.VariationalTwoPt']

    def __init__(self, thisnboot=nboot, cfglist={},Info={},thissym='Not Set',thiscol='Not Set',thisshift=0.0,name='',
                 man_load_cfgs=False):

        ## these class elemnts are checked when loading in pickled file
        self.thischecklist = ['nboot','nt','ProjStr','ism','jsm','tsrc','kud','ks']

        self.latparams = mp.LatticeParameters(Info=Info)
        self.latparams.LoadPickle()
        # self.checkmomlist = self.latparams.momstrlistAvg

        self.nt = self.latparams.nt
        self.nxyz = self.latparams.nxyz
        self.current = 0
        ## Initialising plotting parameters
        self.thisshift = thisshift
        self.thiscol = thiscol
        self.thissym = thissym
        ## the actual shift ammount needs to be scaled by the x axis size
        ## see GetShift
        self.thisshiftscale = 'Not Set'

        #self.C2 = [ ip , it , ,istream , iconf ]
        # self.C2     = np.array([])

        self.C2_cfgs = pa.DataFrame()
        '''
        self.C2_cfgs DataFrame:

        columns =   config_numbers , Op[ip, it]

        rows =      stream , configuration
        '''

        self.C2_Col_Names = ['momentum','source_sink_separation']
        self.C2_Stats = pa.DataFrame()

        '''
        self.C2_Stats DataFrame:

        columns =   boot,   Avg,        Std
                    EffM,   EffMAvg,    EffMStd

        rows =      momenta, time separation
        '''

        # # C2boot [ ip , it ] bs
        # self.C2boot = ODNested()
        # self.C2Avg  = ODNested()
        # self.C2Std  = ODNested()
        #
        # # EffM [ ip , it ] bs
        # self.EffM    = ODNested()
        # self.EffMAvg = ODNested()
        # self.EffMStd = ODNested()


        self.C2_Fit_Col_Names = ['states','momentum']
        self.C2_Fit_Stats = pa.DataFrame()

        '''
        self.C2_Fit_Stats DataFrame:

        columns =   boot

        rows =      states fitted to, momenta


        boot is SetOfFits object
        Avg/Std elements are pa.Series instances with:
            index = fit parameters
            values = Avg/Std values
        '''


        # # FitDict = [ States_Fitted_To, ip, Fit_range ] fitting class
        # self.FitDict    = ODNested()
        # # FitDictAvg/Std [ States_Fitted_To , ip ,  Fit_range , Fit_Parameter ] bs
        # self.FitDictAvg = ODNested()
        # self.FitDictStd = ODNested()



        self.nboot = thisnboot


        ## toggle to show exact configurations and sources used in xml outputfiles
        if 'show_cfgs' in list(Info.keys()): self.show_cfgs = Info['show_cfgs']
        else: self.show_cfgs = False


        ## used to pick out different configuraiton/source lists.
        ## for example, if there are less flowed operators than 2 point correlators, passing 'TopCharge' will read that result
        if 'fo_for_cfgs' in list(Info.keys()): self.fo_for_cfgs = Info['fo_for_cfgs']
        else: self.fo_for_cfgs = False



        ## Info can contain fit ranges to be calculated
        ## must bein form Info['Fits'] = [ ('state#' , 'fit#-#' ) , ...]
        if 'Fits' in list(Info.keys()): self.PredefFits = Info['Fits']
        else: self.PredefFits = []


        ## Info can contain projector used (or type of meson in psibar Gamma psi)
        if 'Interp' in list(Info.keys()):
            try:
                self.ProjStr = 'P'+str(int(Info['Interp']))
            except Exception as err:
                self.ProjStr = Info['Interp']
            self.Interp = GetInterpFlag(Info['Interp'])
        else:
            self.ProjStr = 'CPEven'
            self.Interp = GetInterpFlag('CPEven') ## default to regular CP Even projector

        if self.ProjStr == 'g5':
            self.ProjStr = 'g5div2'
        ## Info can contain projector used (or type of meson in psibar Gamma psi)

        if 'MesOrBar' in list(Info.keys()): self.MesOrBar = Info['MesOrBar']
        else: self.MesOrBar = 'Baryon' ## Default to baryon

        if 'cfg_file_type' in list(Info.keys()): self.cfg_file_type = Info['cfg_file_type']
        else: self.cfg_file_type = cfg_file_type

        ## Info can contain initial guess for fits, make sure it array of correct length
        if 'iGuess' in list(Info.keys()): self.iGuess = np.array(Info['iGuess'])
        else: self.iGuess = 'None' ## Defaults to value set for function in FitFunctions.py

        ## Source Smearing ism value ( in integer form)
        if 'ism' in list(Info.keys()): self.ism = 'ism'+str(Info['ism'])
        else: self.ism = ''

        ## Sink Smearing jsm value ( in integer form)
        if 'jsm' in list(Info.keys()): self.jsm = 'jsm'+str(Info['jsm'])
        else: self.jsm = ''

        ## t_source location ( in integer form)
        if 'tsrc' in list(Info.keys()): self.tsrc = 'tsrc'+str(Info['tsrc'])
        else: self.tsrc = ''

        if 'Sparam' in list(Info.keys()): self.Sparam = Info['Sparam']
        else: self.Sparam = defSparam


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


        # ## x_source list has been pushed to cfglist, this is depriciated
        # if 'xsrc' in Info.keys(): self.xsrc = ['xsrc'+str(isrc) for isrc in Info['xsrc']]
        # else: self.xsrc = ['xsrc1']

        ## momentum ( in form of [px py pz] )
        if 'pmom' in list(Info.keys()):
            if 'All' in Info['pmom']:
                self.pmom = self.latparams.GetMomList(list_type='str_Avg').tolist() ## aparenlty ReadBinary likes lists, not numpy arrays
            else:
                self.pmom = list(map(self.latparams.GetAvgmomstr,Info['pmom']))
        else: self.pmom = self.latparams.GetMomList(list_type='str_Avg').tolist()
        self.pmom = self.latparams.RemoveDup(self.pmom)
        self.checkmomlist = self.pmom

        self.pform = ['p'+ip.replace(' ','') for ip in self.pmom]


        ## kappa up down (in integer form)
        if 'kud' in list(Info.keys()): self.kud = 'kud'+str(Info['kud'])
        else: self.kud = ''

        ## kappa strange (in integer form)
        if 'ks' in list(Info.keys()): self.ks = 'ks'+str(Info['ks'])
        else: self.ks = ''

        ## kappa strange (in integer form)
        if 'Csw' in list(Info.keys()): self.Csw = 'C'+str(Info['Csw'])
        else: self.Csw = 'C1715'
        ## kappa strange (in integer form)
        if 'Beta' in list(Info.keys()): self.Beta = 'B'+str(Info['Beta'])
        else: self.Beta = 'B1900'

        fileism = self.ism.replace('ism','sm')
        filejsm = self.jsm.replace('jsm','si')
        self.kappafolder = 'Kud0'+self.kud.replace('kud','')+'Ks0'+self.ks.replace('ks','')
        self.kappafolder2 = 'Kud0'+self.kud.replace('kud','')[:-2]+'Ks0'+self.ks.replace('ks','')[:-2]
        self.dim_label = 'RC'+str(self.nxyz)+'x'+str(self.nt)
        self.dim_ll = str(self.nxyz)+r'^{3}\times ' + str(self.nt)
        self.readfilename = [self.dim_label+'_'+self.Beta+self.kappafolder+self.Csw,'_',
                             '_k'+self.kud.replace('kud','')+'_'+self.tsrc+fileism+filejsm+'_nucleon.2cf.xml']
        self.readfilename2 = [self.dim_label+'_'+self.Beta+self.kappafolder2+self.Csw,'_',
                             '_k1382500_k1370000_'+self.tsrc+fileism+filejsm+'_nucleon.2cf.xml']
        self.readfilename3 = [self.dim_label+'_'+self.Beta+self.kappafolder2+self.Csw,'_',
                             '_k1370000_k1364000_'+self.tsrc+fileism+filejsm+'_nucleon.2cf.xml']
        if os.path.isdir(cfundir+self.dim_label+'_'+self.kappafolder+'_AY'):
            self.filedir = cfundir+self.dim_label+'_'+self.kappafolder+'_AY/twoptRandT/twopt'+fileism+filejsm+'/'
        elif os.path.isdir(cfundir+self.dim_label+'_'+self.kappafolder):
            self.filedir = cfundir+self.dim_label+'_'+self.kappafolder+'/twoptRandT/twopt'+fileism+filejsm+'/'
        else:
            self.filedir = cfundir+self.kappafolder+'/twoptRandT/twopt'+fileism+filejsm+'/'

        ## cfglist { configuration : [ x-source numbers ] }
        self.Info = Info


        self.SetCustomName(name)
        self.LogScale = 1.0
        if not man_load_cfgs:
            self.LoadCfgs(cfglist,name)

    def LoadCfgs(self,cfglist,name=''):

        if len(cfglist) == 0:
            self.GetCfgList()
        else:
            self.ImportCfgList(cfglist)
        self.SetCustomName(name)


    def GetFuns(self):
        if self.C2_Fit_Stats is not None and 'boot' in self.C2_Fit_Stats.columns:
            for ival in self.C2_Fit_Stats['boot'].values:
                ival.GetFuns()

    def RemoveFuns(self):
        if self.C2_Fit_Stats is not None and 'boot' in self.C2_Fit_Stats.columns:
            for ival in self.C2_Fit_Stats['boot'].values:
                ival.RemoveFuns()

    # ## thisInfo must contain correct information for combining the correlators (the smearings, the eigenvectors, etc...)
    # def GetCombCorr(self,thisInfo):
    #     from SetsOfTwoPtCorrs import SetOfTwoPt
    #     sinklab,jsmcoeffs,sm2ptList,InfoList = GetInfoList(thisInfo)

    #     C2Set = SetOfTwoPt(cfglist=self.NNQ_cfgs[['configs','xsrc_list']],InfoDict = InfoList)
    #     C2Set.LoadPickle(CheckCfgs=CheckCfgs)
    #     thisname,thisleglab,jsmlab = CreateNameAndLL(C2Set.SetC2,sinklab)
    #     thisfun = GenSinkFun(jsmcoeffs)
    #     thisfile = self.filename.replace(self.jsm,sinklab)
    #     thisLL = self.LegLab.replace(self.jsm,sinklab)
    #     # self.SetCustomName(thisfile,thisLL)
    #     # self.jsm = jsmlab
    #     return C2Set.CombineCorrs(thisfun,filename=thisfile,LegLab=thisLL,jsmlab=jsmlab)



    def GetCfgList(self,name=''):
        def SplitCfgSrc(ifile):
            isrc = re.findall('_xsrc.*_k',ifile)
            icfg = re.findall('-.-00.*_xsrc',ifile)
            istream = re.findall('-.-',ifile)
            if len(isrc) == 0 or len(icfg) == 0:
                print()
                print('WARNING: file does not have xsrc or -.-####### (source or cfg number)')
                print(ifile)
                print()
                return '','',0,0
            isrc = isrc[0].split('_')[1]
            icfg = icfg[0].replace('_xsrc','')[5:]
            istream = istream[0]
            # return istream,isrc,icfg,int(icfg),int(isrc.replace('xsrc',''))
            return istream,isrc,icfg

        cfg_dict = ODNested()
        this_check = isinstance(self.stream_list,str) and self.stream_list == 'all'
        for ifile in glob.glob(self.filedir+'/*'):
            istream,ixsrc,icfg = SplitCfgSrc(ifile)
            if this_check or istream in self.stream_list:
                if istream not in list(cfg_dict.keys()) or icfg not in list(cfg_dict[istream].keys()):
                    cfg_dict[istream][icfg] = [ixsrc]
                else:
                    cfg_dict[istream][icfg].append(ixsrc)
        if len(list(cfg_dict.keys())) == 0:
            print('no configs found in ')
            print(self.filedir)


        cfglist,xsrclist,ilist = [],[],[]
        this_lambda = lambda ix : int(ix.replace('xsrc',''))
        for istream,stream_dict in cfg_dict.items():
            for icfg,(the_cfg,ixsrc_list) in enumerate(stream_dict.items()):
                ilist.append((istream,icfg))
                cfglist.append(the_cfg)
                ## sorting xsrc list with respect to its INTEGER value (not string sort!)
                ixsrc_int_list = list(map(this_lambda,ixsrc_list))
                xsrclist.append([ix for _,ix in sorted(zip(ixsrc_int_list,ixsrc_list))])
        ## finds cfgs and xsrcs conbinations that are in all streams.
        ##MAKE SURE TO ORDER CFG AND SRC

        if len(ilist) == 0:
            # raise EnvironmentError('No configs found for \n'+name)
            print('Warning, No configs found for \n',name)
        else:
            indicies = pa.MultiIndex.from_tuples(ilist,names=['stream','config_number'])
            cfg_df = pa.DataFrame(cfglist,columns=['configs'],index=indicies)
            cfg_df.loc[:,'xsrc_list'] = pa.Series(xsrclist,index=indicies)
            cfg_df.loc[:,'int_configs'] = pa.Series(list(map(int,cfglist)),index=indicies)
            cfg_df = cfg_df.sort_values('int_configs').sort_index(level='stream',sort_remaining=False)
            del cfg_df['int_configs']
            self.ImportCfgList(cfg_df)


    # def MakeOutputFile(self):
    #     self.name = '_'.join([self.tsrc,self.ism,self.jsm])
    #     self.HumanFile = '/'+self.kappafolder+'/G2/'+self.name+'.xml'
    #     self.PickleFile = '/'+self.kappafolder+'/G2/Pickle/'+self.name+'.py3p'


    ## Routines for Plotting ###############################################

    def ImportPlotParams(self,thiscol,thissym,thisshift):
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        self.thisshift = thisshift

    def CheckCol(self,thiscol):
        if 'PreDefine' == thiscol :
            if 'Not Set' == self.thiscol:
                raise IOError('Pass in color to initialize it')
        else:
            self.thiscol = thiscol

    def CheckSym(self,thissym):
        if 'PreDefine' == thissym :
            if 'Not Set' == self.thissym:
                raise IOError('Pass in symbol to initialize it')
        else:
            self.thissym = thissym


    def GetShift(self,xlims,thisshift):
        if thisshift != 'PreDefine': self.thisshift = thisshift
        # if self.thisshiftscale == 'Not Set':
        xlen = np.abs(xlims[1]-xlims[0])
        self.thisshiftscale = self.thisshift*xlen
        return self.thisshiftscale
        # else:
        #     return self.thisshiftscale




    ## make sure all xlims are the same, because this sets the window of xlim (needed for calculating shifts)
    def LogPlot(self,plot_class,xlims='All',momlist=['p000'],
                thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',norm=False):
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        self.thisshiftscale = 'Not Set'
        if xlims == 'All':
            xlims = range(self.nt+1)
        thisshift = self.GetShift(xlims,thisshift)
        if len(momlist) == 0:
            momlist = self.C2_Stats.index.get_level_values('momentum')
        for ip in momlist:
            if 'All' in ip: ip = 'p000'
            if ip not in self.C2_Stats.index.get_level_values('momentum'):
                print(ip, ' not found in effective mass list')
                continue
            pdata = self.C2_Stats['boot'].xs(ip, level='momentum')
            tlist,databoot = pdata.index,pdata.values
            dataavg,dataerr = [],[]
            if norm:
                scale = databoot[xlims[0]].Log()
                scale.Stats()
                scale = scale.Avg
            else:
                scale = 1.0
            for idata in databoot[xlims[0]:xlims[1]]:
                logdata = idata.Log().__div__(scale)
                logdata.Stats()
                dataavg.append(logdata.Avg)
                dataerr.append(logdata.Std)
            tlist = tlist[xlims[0]:xlims[1]]
            hold_series = pa.Series()
            hold_series['x_data'] = np.array(list(map(untstr,tlist)))
            hold_series['y_data'] = dataavg
            hold_series['yerr_data'] = dataerr
            hold_series['xerr_data'] = None
            hold_series['type'] = 'error_bar'
            hold_series['fit_class'] = None
            hold_series['label'] = self.LegLab
            hold_series['symbol'] = self.thissym
            hold_series['color'] = self.thiscol
            hold_series['shift'] = thisshift
            hold_series['xdatarange'] = xlims
            hold_series['scale'] = self.LogScale
            hold_series['Phys'] = None
            hold_series['fmt_class'] = KeyForamtting(self.latparams)
            plot_class.AppendData(hold_series)
            self.LogScale=scale
        return plot_class
        # plot_class.get_xlim(*xlims)

    ## make sure all xlims are the same, because this sets the window of xlim (needed for calculating shifts)
    def Plot(   self,plot_class,momlist=['p000'],
                thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',norm=False,log=False):

        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        ## resets thisshiftscale so that you can plot EffMass and Log plot next to eachother,

        self.thisshiftscale = 'Not Set'
        # thisshift = self.GetShift(xlims,thisshift)
        if len(momlist) == 0:
            momlist = self.pform
        # m_getattr = lambda x, y: getattr(y,x)
        ip = momlist[0]
        if 'All' in ip: ip = 'p000'
        if ip not in self.pform:
            print(ip, ' not found in effective mass list')
            return plot_class
        # pdata = self.C2_Stats['EffM'].xs(ip, level='momentum')
        pdata = self.C2_Stats['boot']
        if log:
            pdata = pdata.apply(lambda x : x.Log())
            pdata.apply(lambda x : x.Stats())
        dataavg = pdata.apply(lambda x : x.Avg)
        dataerr = pdata.apply(lambda x : x.Std)
        # tlist = pdata.index[xlims[0]:xlims[1]]
        # dataavg = dataavg.iloc[xlims[0]:xlims[1]]
        # dataerr = dataerr.iloc[xlims[0]:xlims[1]]
        ## np.abs makes negative effective masses coming from G2^-1 be positive.
        hold_series = pa.Series()
        hold_series['x_data'] = 'from_keys'
        # hold_series['x_data'] = np.array(map(untstr,tlist))
        hold_series['y_data'] = dataavg
        hold_series['key_select'] = (ip,slice(None))
        hold_series['yerr_data'] = dataerr
        hold_series['xerr_data'] = None
        hold_series['type'] = 'error_bar_vary'
        hold_series['fit_class'] = None
        hold_series['label'] = self.LegLab
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        # hold_series['xdatarange'] = xlims
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class
        # plot_class.get_xlim(*xlims)


    ## make sure all xlims are the same, because this sets the window of xlim (needed for calculating shifts)
    def EffMassPlot(self,plot_class,xlims = defxlim,momlist=['p000'],
                    thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',Phys=True):
        def DoPhys(value):
            thisval = value*self.latparams.hbarcdivlat
            thisval.Stats()
            return thisval

        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        ## resets thisshiftscale so that you can plot EffMass and Log plot next to eachother,

        self.thisshiftscale = 'Not Set'
        thisshift = self.GetShift(xlims,thisshift)
        if len(momlist) == 0:
            momlist = self.pform
        # m_getattr = lambda x, y: getattr(y,x)
        ip = momlist[0]
        if 'All' in ip: ip = 'p000'
        if ip not in self.pform:
            print(ip, ' not found in effective mass list')
            return plot_class
        # pdata = self.C2_Stats['EffM'].xs(ip, level='momentum')
        pdata = self.C2_Stats['EffM']
        if Phys: pdata = pdata.apply(DoPhys)
        dataavg = pdata.apply(lambda x : x.Avg)
        dataerr = pdata.apply(lambda x : x.Std)
        # tlist = pdata.index[xlims[0]:xlims[1]]
        # dataavg = dataavg.iloc[xlims[0]:xlims[1]]
        # dataerr = dataerr.iloc[xlims[0]:xlims[1]]
        ## np.abs makes negative effective masses coming from G2^-1 be positive.
        hold_series = pa.Series()
        hold_series['x_data'] = 'from_keys'
        # hold_series['x_data'] = np.array(map(untstr,tlist))
        hold_series['y_data'] = dataavg
        hold_series['key_select'] = (ip,slice(None))
        hold_series['yerr_data'] = dataerr
        hold_series['xerr_data'] = None
        hold_series['type'] = 'error_bar_vary'
        hold_series['fit_class'] = None
        hold_series['label'] = self.LegLab
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        hold_series['xdatarange'] = xlims
        hold_series['Phys'] = Phys
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class
        # plot_class.get_xlim(*xlims)

    ## state is 1 or 2 state fit
    ## fitr is fit range e.g. fitr5-10
    ## fit range of 'Data' is just the original fit range
    ## make sure to plot this after LogPlot or the xlims wont work properly
    def LogFitPlot( self,plot_class,state,fitr,
                    xlims = 'Data',momlist=['p000'],thiscol='PreDefine',thisshift='PreDefine'):
        return self.EffMassFitPlot( plot_class,state,fitr,
                        xlims = xlims,momlist=momlist,thiscol=thiscol,thisshift=thisshift,
                        plot_type='logplot')
    ## state is 1 or 2 state fit
    ## fitr is fit range e.g. fitr5-10
    ## fit range of 'Data' is just the original fit range
    ## make sure to plot this after EffMassPlot or the xlims wont work properly
    def EffMassFitPlot( self,plot_class,state,fitr,
                        xlims = 'Data',momlist=['p000'],thiscol='PreDefine',thisshift='PreDefine',
                        Phys=True,plot_type='effmass'):
        self.CheckCol(thiscol)
        thisshift = self.GetShift(plot_class.get_xlim(),thisshift)
        if isinstance(momlist,(list,tuple,np.ndarray)):
            print('momlist is not needed anymore, picking p000')
            momlist = 'p000'
        if isinstance(state,int): state = 'state'+str(state)
        if 'boot' not in self.C2_Fit_Stats.columns:
            print('no fits done, starting')
            self.Fit(int(state.replace('state','')),fitr)
        fit_data = sff.PullSOFSeries(self.C2_Fit_Stats['boot'],fmted=True)
        if len(fit_data) == 0: return plot_class
        # fit_data = Series_fix_key(fit_data,4,'fittwor','tsumfitr')
        this_fun = fit_data.index[0][2]
        this_key = (state,momlist,this_fun,fitr)
        if this_key not in fit_data.index:
            print(this_key, ' not found in FitDict, performing fit now')
            print(fit_data)
            self.Fit(int(state.replace('state','')),fitr)
        hold_series = null_series
        if 'effmass' in plot_type.lower():
            hold_series['type'] = 'effm_fit_vary'
        elif 'log' in plot_type.lower():
            hold_series['type'] = 'log_fit_vary'

        hold_series['fit_class'] = fit_data
        hold_series['key_select'] = this_key
        hold_series['label'] = 'Exp Fit'
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['Phys'] = Phys
        hold_series['shift'] = thisshift
        hold_series['xdatarange'] = xlims
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class
    ## ##################### ###############################################

    ## state is 1 or 2 state fit
    ## fitr is fit range e.g. fitr5-10
    ## fit range of 'Data' is just the original fit range
    ## make sure to plot this after EffMassPlot or the xlims wont work properly
    ## TODO, fix this
    def RecRelPlot(self,plot_class,state,fitr,thiscol='PreDefine',thisshift='PreDefine',thissym='PreDef',Phys=True):
        print('RecRelPlot Not implemeted for ImpSOF, ON THE TODO LIST')
        return plot_class
        # def DoPhys(value):
        #     thisval = value*self.latparams.hbarcdivlat
        #     thisval.Stats()
        #     return thisval
        #
        # self.CheckCol(thiscol)
        # self.CheckSym(thissym)
        # if isinstance(state,int): state = 'state'+str(state)
        # # fitr_params = self.FitDict[state][fitr]
        # for imom in self.pform:
        #     if (state,imom,fitr) not in self.C2_Fit_Stats.index:
        #         print state,imom,fitr , ' has not been done, performing fit now'
        #         self.Fit(int(state.replace('state','')),[fitr])
        # this_fit = self.C2_Fit_Stats['boot']
        # varl,errl,indexl,bootl = [],[],[],[]
        # for (istate,imom,ifitr),ifit in this_fit.iteritems():
        #     for ipar,parval in ifit.fit_data['Params'].iteritems():
        #         indexl.append((istate,imom,ifitr,ipar))
        #         if ipar == 'Energy' and Phys:
        #             parval = parval*self.latparams.hbarcdivlat
        #             parval.Stats()
        #         bootl.append(parval)
        #         varl.append(parval.Avg)
        #         errl.append(parval.Std)
        # if len(indexl) > 0:
        #     indicies = pa.MultiIndex.from_tuples(indexl,names=list(self.C2_Fit_Stats.index.names) + ['Parameter'])
        #     fit_boot = pa.Series(bootl,index=indicies)
        #     fit_data = pa.Series(varl,index=indicies)
        #     fiterr_data = pa.Series(errl,index=indicies)
        #
        #
        # xboot = [self.latparams.EnergyFromRecRel(imom,zboot) for imom,zboot in fit_boot.loc[(istate,slice(None),ifitr,ipar)].iteritems()]
        # ploty = fit_data.loc[(istate,slice(None),ifitr,ipar)].values
        # plotx = map(lambda x : x.Avg,xboot)
        # plotxerr = map(lambda x : x.Std,xboot)
        # thisshift = self.GetShift((min(plotx),max(plotx)),thisshift)
        # hold_series = pa.Series()
        # hold_series['x_data'] = np.array(plotx)
        # hold_series['y_data'] = fit_data
        # hold_series['xerr_data'] = plotxerr
        # hold_series['yerr_data'] = fiterr_data
        # hold_series['type'] = 'error_bar_vary'
        # hold_series['key_select'] = (istate,slice(None),ifitr,ipar)
        # hold_series['fit_class'] = None
        # hold_series['label'] = self.LegLab
        # hold_series['symbol'] = self.thissym
        # hold_series['color'] = self.thiscol
        # hold_series['shift'] = thisshift
        # plot_class.AppendData(hold_series)
        # return plot_class,max([max(ploty),max(plotx)]),min([min(ploty),min(plotx)])



    def SetCustomName(self,string='',stringLL='',fo_for_cfgs='PreDef'):
        if not hasattr(self,'stream_name'):
            self.stream_name = '-def-'
        if string == '':
            self.name = '_'.join([self.dim_label,self.kappafolder,self.stream_name,self.MesOrBar,self.ProjStr,self.tsrc,self.ism,self.jsm])
            self.filename = '_'.join([self.MesOrBar,self.stream_name,self.ProjStr,self.tsrc,self.ism,self.jsm])
            # self.filename = self.name
        else:
            self.filename = self.name = string
        if stringLL == '':
            self.LegLab = '$'+'\ '.join([self.dim_ll,self.latparams.GetPionMassLab(),self.stream_name,self.MesOrBar,self.ProjStr,self.ism,self.jsm])+'$' ## customise this how you like
        else:
            self.LegLab = stringLL


        if isinstance(fo_for_cfgs,str):
            if fo_for_cfgs == 'PreDef':
                fo_for_cfgs = self.fo_for_cfgs
        if isinstance(fo_for_cfgs,str):
            self.G2dir =        outputdir + '/'+self.dim_label+self.kappafolder+'/G2/'+fo_for_cfgs+'/'
            self.G2Cfgsdir =    cfgfmtdir + '/'+self.dim_label+self.kappafolder+'/G2/'+fo_for_cfgs+'/'
        else:
            self.G2dir = outputdir + '/'+self.dim_label+self.kappafolder+'/G2/'
            self.G2Cfgsdir = cfgfmtdir + '/'+self.dim_label+self.kappafolder+'/G2/'
        mkdir_p(self.G2dir+'/Pickle/')
        mkdir_p(self.G2dir+'/Excel/')
        mkdir_p(self.G2Cfgsdir+'/')
        self.HumanFile = self.G2dir+self.filename+'.xml'
        self.ExcelFile = self.G2dir+'/Excel/'+self.filename+'.xlsx'
        self.PickleFile = self.G2dir+'/Pickle/'+self.filename+'.py3p'
        self.CfgFile = self.G2Cfgsdir+'/'+self.filename+'.cfgs'

    def ImportCfgList(self,cfgsrclist,stream_list = 'PreDefine'):
        if not isinstance(cfgsrclist,pa.DataFrame):
            if not isinstance(stream_list,str):
                self.stream_list = stream_list
            ## cfgsrclist is [cfglist , xsrclist]
            ## cfglist must be 2-d array, first dimension is stream, second dimension is configuraiton number
            ## xsrclist is 3-d array, first dimension is stream, second dimension is cfg number, third is xsrc number
            ## xsrclist must have values 'xsrc##' with no zero pre buffering
            ilist = []
            xsrc_list = []
            for istream,(is_cfg,is_xsrc) in enumerate(zip(cfgsrclist[0],cfgsrclist[1])):
                for icf in range(is_cfg):
                    xsrc_list.append(is_xsrc)
                    ilist.append((istream,is_cfg))
            indicies = pa.MultiIndex.from_tuples(ilist,names=['stream','config_number'])
            cfg_df = pa.DataFrame(np.array(cfgsrclist[0]).flatten(),columns=['configs'],index=indicies)
            cfg_df.loc[:,'xsrc_list'] = pa.Series(xsrc_list,index=indicies)
            self.C2_cfgs = cfg_df
        else:
            self.C2_cfgs = cfgsrclist
        if 'configs' in self.C2_cfgs.columns:
            self.C2_cfgs.loc[:,'stream_cfgs'] = pa.Series(['00'.join((istream,icfg)) for (istream,iccfg),icfg in self.C2_cfgs['configs'].items()],index=self.C2_cfgs.index)
            # print
            # print 'DEBUG importing cfgs'
            # print self.C2_cfgs
            self.Update_Stream_List()
            self.ncfg_list = [icfg.size for istream,icfg in self.C2_cfgs['configs'].groupby(level='stream')]
            self.xsrcAvg = np.mean(list(map(len,self.C2_cfgs['xsrc_list'].values)))
            self.nMeas = self.xsrcAvg*sum(self.ncfg_list)
        ## filecfglist [ istream , icfg , isrc ]

    def Update_Stream_List(self):
        if 'configs' in self.C2_cfgs:
            self.stream_list = self.C2_cfgs.index.levels[0]
            self.nstream = len(self.stream_list)
            if hasattr(self,'stream_name'):
                old_sn = self.stream_name
                self.stream_name = '-'+''.join([istr.replace('-','') for istr in self.stream_list])+'-'
                if hasattr(self,'filename'):
                    self.filename = self.filename.replace(old_sn,self.stream_name)
                if hasattr(self,'name'):
                    self.name = self.name.replace(old_sn,self.stream_name)
                if hasattr(self,'LegLab'):
                    self.LegLab = self.LegLab.replace(old_sn,self.stream_name)
            else:
                self.stream_name = '-'+''.join([istr.replace('-','') for istr in self.stream_list])+'-'


    def Stats(self):
        if 'boot' not in self.C2_Stats: return
        lAvg,lStd = [],[]
        for (ip,it),iC2 in self.items():
            iC2.Stats()
            lAvg.append(iC2.Avg)
            lStd.append(iC2.Std)
        indicies = pa.MultiIndex.from_tuples(self.C2_Stats.index,names=self.C2_Col_Names)
        self.C2_Stats.loc[:,'Avg'] = pa.Series(lAvg,index=indicies)
        self.C2_Stats.loc[:,'Std'] = pa.Series(lStd,index=indicies)
        self.EffMass()


    def RemoveVals(self):
        if 'C2' in self.C2_cfgs: del self.C2_cfgs['C2']


    #self.C3 = [ igamma , ip , it , iconf ]
    def WriteCfgsToFile(self,thisdir,show_timer=True):
        thisoutdir = thisdir + '/'+self.dim_label+'_'+self.kappafolder+'/'+self.ProjStr+'_'.join([self.ism,self.jsm])+'/'
        # for igamma,gammadata in zip(self.gammalist,self.C3_cfgs['C3']):
        mkdir_p(thisoutdir)
        if show_timer: thistimer = Timer(linklist=self.C2_cfgs['configs'].values,name='Writing 2-ptCorrs ')
        for ((istream,iccfg),cfgdata),icfg in zip(iter(self.C2_cfgs['C2'].items()),self.C2_cfgs['configs'].values):
            cfgstr = 'C2'+istream+'cfg'+str(icfg).zfill(5)
            thisfile = thisoutdir + '/'+cfgstr+'.xml'
            dictout = ODNested()
            for ip,qdata in zip(self.pform,cfgdata):
                for it,tdata in enumerate(qdata):
                    itstr = tstr(it+1)
                    dictout[ip][itstr] = tdata
            WriteXml(thisfile,{'Results':dictout})
            if show_timer: thistimer.Lap()

    #self.C3 = [ igamma , ip , it , iconf ]
    def WriteBootToFile(self,thisdir,show_timer=True):
        thisoutdir = thisdir + '/'+self.dim_label+'_'+self.kappafolder+'/'+self.ProjStr+'_'.join([self.ism,self.jsm])+'/'
        # for igamma,gammadata in zip(self.gammalist,self.C3_cfgs['C3']):
        mkdir_p(thisoutdir)
        if show_timer: thistimer = Timer(linklist=list(range(1,self.nboot+1)),name='Writing booted 2-ptCorrs ')
        this_nboot = self.C2_Stats.loc[:,'boot'].iloc[0].bootvals.index
        for iboot in this_nboot:
            boot_val = self.C2_Stats.loc[:,'boot'].apply(lambda ival : ival.bootvals.loc[iboot])
            cfgstr = 'C2_iboot'+str(iboot).zfill(5)
            thisfile = thisoutdir + '/'+cfgstr+'.xml'
            dictout = ODNested()
            for (ip,it),ival in boot_val.items():
                dictout[ip][it] = ival
            WriteXml(thisfile,{'Results':dictout})
            if show_timer: thistimer.Lap()


    ## self.C2 = [ ip, it , iconf ]
    def Read(self,show_timer=True,Full=False,file_type='PreDef'):
        ## making interpolator number more human readable
        if Full:
            self.CfgFile = self.CfgFile.replace('.cfgs','_Full.cfgs')
        if self.Read_Cfgs(file_type=file_type,show_timer=show_timer):
            return
        thisdata = []
        prop_src = []
        # if len(self.filename.split('_')) > 5:
        #     raise IOError('Do not read after combining correlators, please recreate the correlator')
        thistfix = False
        # if self.nt == 40 and self.nxyz == 20: thistfix = True
        if show_timer: thistimer = Timer(linklist=self.C2_cfgs['stream_cfgs'].values,name='Read '+self.name)
        for i_stream_cfg,ixsrc_list in zip(self.C2_cfgs['stream_cfgs'].values,self.C2_cfgs['xsrc_list'].values):
            thisxsrclist = [self.filedir+ self.readfilename[0]+i_stream_cfg+self.readfilename[1]+ixsrc+self.readfilename[2] for ixsrc in ixsrc_list]
            if not os.path.isfile(thisxsrclist[0]):
                thisxsrclist = [self.filedir+ self.readfilename2[0]+i_stream_cfg+self.readfilename2[1]+ixsrc+self.readfilename2[2] for ixsrc in ixsrc_list]
            if not os.path.isfile(thisxsrclist[0]):
                thisxsrclist = [self.filedir+ self.readfilename3[0]+i_stream_cfg+self.readfilename3[1]+ixsrc+self.readfilename3[2] for ixsrc in ixsrc_list]

            if Full:
                idata = RC2Full(thisxsrclist,self.pmom,InterpNumb=self.Interp,MesOrBar=self.MesOrBar)
            else:
                idata = R2CChromaXMLFileList(thisxsrclist,self.pmom,InterpNumb=self.Interp,MesOrBar=self.MesOrBar,tfix=thistfix)

            thisdata.append(idata.data)
            prop_src.append(idata.tshiftlist)

            if show_timer: thistimer.Lap()
        if Full: self.C2_cfgs.loc[:,'tsrc_list'] = pa.Series(prop_src,index=self.C2_cfgs.index)

        self.C2_cfgs.loc[:,'C2'] = pa.Series(thisdata,index=self.C2_cfgs.index)
        if self.Interp == '1':
            self.C2_cfgs.loc[:,'C2'] = self.C2_cfgs.loc[:,'C2'].apply(lambda x : (np.array(x)/2.).tolist())
        self.Write_Cfgs(show_timer=show_timer,file_type=file_type)

        # self.C2_cfgs.to_csv('./Debug2pt_cfgs.csv')
        # pa.Series(self.Info).to_csv('./Debug2pt_Info.csv')

        # print 'DEBUG mean'
        # for ival in range(self.nt):
        #     print 't'+str(ival),self.C2_cfgs.loc[:,'C2'].apply(lambda x : x[0][ival]).mean()
        #     # print self.C2_cfgs.loc[:,'C2'].apply(lambda x : x[0][ival]).to_csv('./Debug2pt_C2.csv')
        # print
        #
        # print 'DEBUG t1 in stream -2-'
        # for ikey,ival in self.C2_cfgs.loc[:,'C2'].iteritems():
        #     print ikey,ival[0][0]
        #     # self.C2_cfgs.loc[:,'stream_cfgs'].to_csv('./Debug2pt_cfgs.csv')
        #     # print self.C2_cfgs.loc[:,'C2'].apply(lambda x : x[0][ival]).to_csv('./Debug2pt_C2.csv')
        # print

    def Bootstrap(self,WipeData=True):
        if 'C2' not in self.C2_cfgs:
            self.Read()
        lC2,lC2Avg,lC2Std = [],[],[]
        lCA2,lCA2Avg,lCA2Std = [],[],[]
        ilist = []
        rlist = None
        # print 'DEBUG',self.name
        for icp,ip in enumerate(self.pform):
            for it in range(self.latparams.nt):
                strit = tstr(it+1)
                ilist.append((ip,strit))
                this_lambda = lambda x : x[icp][it]
                tdata = self.C2_cfgs['C2'].apply(this_lambda)

                thisboot = BootStrap(self.nboot, name='G_{2}('+ip+','+strit+')',cfgvals=tdata,rand_list=rlist)
                rlist = thisboot.Get_Rand_List()
                # print ip,it,'bootval=',thisboot.Avg,'bootstd',thisboot.Std,'diff=', np.abs(thisboot.Avg-np.mean(tdata))/np.mean(tdata)
                lC2.append(thisboot)
                lC2Avg.append(thisboot.Avg)
                lC2Std.append(thisboot.Std)
                def eye_fun(*x):
                    return x[0]
                def one_fun(*x):
                    return [x[0]/x[0]]
                tdata = tdata.to_frame('C2')
                thisAuto =  AutoCorrelate(  Fun=(eye_fun,one_fun),Sparam=self.Sparam,
                                            name=self.name + ip+' $t='+str(it)+'$',data=tdata)
                lCA2.append(thisAuto)
                lCA2Avg.append(thisAuto.Avg)
                lCA2Std.append(thisAuto.Std)
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.C2_Col_Names)
            self.C2_Stats.loc[:,'boot'] = pa.Series(lC2,index=indicies)
            self.C2_Stats.loc[:,'Avg'] = pa.Series(lC2Avg,index=indicies)
            self.C2_Stats.loc[:,'Std'] = pa.Series(lC2Std,index=indicies)
            self.C2_Stats.loc[:,'Auto'] = pa.Series(lCA2,index=indicies)
            self.C2_Stats.loc[:,'AutoAvg'] = pa.Series(lCA2Avg,index=indicies)
            self.C2_Stats.loc[:,'AutoStd'] = pa.Series(lCA2Std,index=indicies)

        if WipeData: self.RemoveVals()


    # def ReadAndBoot(self):
    #     self.Read()
    #     self.Bootstrap()

    def ImportFitRanges(self,fit_range):
        if fit_range == 'PreDef':
            self.fit_range = self.PredefFits[0][1]
        else:
            self.fit_range = fit_range
            if isinstance(self.fit_range,(list,tuple,np.ndarray)):
                print('SetsOfFits implementation only takes 1 fit range to vary over, choosing first')
                self.fit_range = self.fit_range[0]


    ## Fitting the correlator
    ## fit_range is either list of fit ranges to fit the correlator to, or a single fit range
    ## formatted as fitr#-# where # are the min and max fit values.
    ## states correspond to how many states to fit to (currenlty only working for states = 1 or 2)
    ## NB, iGuess here is not implemented yet, needs to be passed into the function
    def Fit(self,state='PreDef',fit_range='PreDef',iGuess='PreDef',EstDir=False,Ratio=False,WipeFit=False):

        if iGuess != 'PreDef': self.iGuess = iGuess

        self.ImportFitRanges(fit_range)
        if state == 'PreDef':
            state = self.PredefFits[0][0]

        if isinstance(state,str):
            if 'state' in state:
                state = state.replace('state','')
            state = int(state)
        if self.MesOrBar == 'Baryon':
            if state == 1:
                this_fun = [C2OneStateFitFunNoExp,2]
            elif state == 2:
                this_fun =[C2TwoStateFitFunNoExp,4]
            elif state > 2:
                raise IOError('fitting state for Bayrons is only implemented up to 2 state so far. Cannot do '+str(state))
            this_fit_info = {'Funs':this_fun}
        elif self.MesOrBar == 'Meson':
            if state == 1:
                def C2atT(t,p):
                    return C2OSFAntiper(t,p,self.latparams.nt+1)
                C2atT.file_name = 'C2atT'+str(self.latparams.nt+1)
                this_fun = [C2atT,2]
                def C2DeratT(t,p):
                    return C2OSFAntiperDer(t,p,self.latparams.nt+1)
                C2DeratT.file_name = 'C2DeratT'+str(self.latparams.nt+1)
            elif state > 1:
                raise IOError('fitting state for Mesons is only implemented up to 1 state so far. Cannot do '+str(state))
            this_fit_info = {'Funs':this_fun,'FunDer':C2DeratT}

        ststr = 'state'+str(state)
        lFit,ilist = [],[]
        for ip,pdata in self.C2_Stats['boot'].groupby(level='momentum'):
            # print 'Fitting state',str(state) , ' ' , ifit
            this_key = (ststr,ip)
            if 'boot' in self.C2_Fit_Stats and this_key in self.C2_Fit_Stats['boot']:
                this_sof = self.C2_Fit_Stats.loc[this_key,'boot']
                if isinstance(this_sof,pa.Series):
                    this_sof = this_sof.iloc[0]
                if this_sof.CheckRange_fitr(self.fit_range,fit_info=this_fit_info):
                    if WipeFit:
                        self.C2_Fit_Stats.drop([(ststr,ip)],inplace=True)
                    else:
                        continue
            # print '     Fitting ', ip
            thisff = sff.SetOfFitFuns(data=pdata[ip],name='_'.join(this_key))
            if self.MesOrBar == 'Baryon':
                thisff.ScanRange_fitr(self.fit_range,fit_info=this_fit_info)
            elif self.MesOrBar == 'Meson':
                thisff.ScanRange_fitr_Sym(self.fit_range,fit_info=this_fit_info)
            lFit.append(thisff)
            ilist.append(this_key)
        if len(ilist)>0:
            # if Debug_Mode:
            #     print 'DEBUG'
            #     print self.C2_Fit_Col_Names
            #     print ilist
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.C2_Fit_Col_Names)
            if 'boot' in self.C2_Fit_Stats.columns:
                this_df = pa.DataFrame()
                this_df.loc[:,'boot'] = pa.Series(lFit,index=indicies)
                self.C2_Fit_Stats = self.C2_Fit_Stats.append(this_df)
            else:
                self.C2_Fit_Stats.loc[:,'boot'] = pa.Series(lFit,index=indicies)
            # def DoFits(val):
            #     val.DoFits()
            #     # if state == 1:
            #     #     val.UnExpParams(paramkeys=('Energy',))
            #     #     val.ImportFunction( C2OneStateFitFunNoExp,2)
            #     # elif state == 2:
            #     #     val.UnExpParams(paramkeys=('Energy','DE'))
            #     #     val.ImportFunction( C2TwoStateFitFunNoExp,4)
            #     return val
            self.C2_Fit_Stats['boot'] = self.C2_Fit_Stats['boot'].apply(lambda x : x.DoFits())
            self.Write()



    ## Comparisons

    def __str__(self):
        return self.name


    def __iter__(self):
        return self

    def items(self):
        return list(self.C2_Stats['boot'].items())

    def iteritems(self):
        return iter(self.C2_Stats['boot'].items())

    def itemsAvgStd(self):
        outlist = []
        for ivals in self.C2_Stats.itertuples():
            outlist.append((ivals[0],ivals[2],ivals[3]))
        return outlist


    def values(self):
        return self.C2_Stats['boot'].values

    def valuesAvgStd(self):
        outlist = []
        for ivals in self.C2_Stats.itertuples(index=False):
            outlist.append((ivals[1],ivals[2]))
        return outlist


    def keys(self):
        return self.C2_stats.index


    def __setitem__(self, key, value ):
        self.C2_Stats.loc[key,'boot'] = value

    def __getitem__(self, key):
        return self.C2_Stats.loc[key,'boot']

    def __reversed__(self):
        ## TODO this is incorrect, fix if required
        self.C2_Stats = self.C2_Stats.iiloc[::-1]
        return self


    def __next__(self):
        if self.current >= len(list(self.keys())):
            self.current = 0
            raise StopIteration
        else:
            self.current += 1
            if self.current >= len(list(self.keys()))+1:
                return 0.0
            else:
                thiskey = list(self.keys())[self.current-1]
                return self[thiskey]


    def EffMass(self):
        if not hasattr(self,'MessEffMass'):
            self.MessEffMass = False
        overwrite_effm = False
        if self.MesOrBar == 'Meson':
            try:
                from pynverse import inversefunc
                overwrite_effm = not self.MessEffMass
                self.MesEffMass = True
            except Exception as err:
                warn('pynverse needs to be install to analytically solve the effective mass \
                     function for Mesons, using incorrect effective mass funciton. Run pip install \
                     pynverse to get it')
                self.MesEffMass = False
        else:
            self.MesEffMass = False
        if 'EffM' not in self.C2_Stats or overwrite_effm:
            lEM,lEMA,lEMS = [],[],[]
            if self.MesEffMass:
                for (ip,it),idata in self.items():
                    nextt = (untstr(it) % (self.nt-1)) + 1
                    idata_nextt = self.C2_Stats['boot'][ip,tstr(nextt)]
                    Tshift = self.latparams.nt/2-untstr(it) + 1
                    def ThisEffMFun(mass):
                            return np.cosh(-mass*Tshift)/np.cosh(-mass*(Tshift-1))
                    thisEmass = BootStrap(bootvals=np.abs(inversefunc(ThisEffMFun,
                                                               y_values=(idata/idata_nextt).bootvals)))
                    thisEmass.Stats()
                    lEM.append(thisEmass)
                    lEMA.append(thisEmass.Avg)
                    lEMS.append(thisEmass.Std)
            else:
                for (ip,it),idata in self.items():
                        nextt = (untstr(it) % (self.nt-1)) + 1
                        thisEmass = (idata/self.C2_Stats['boot'][ip,tstr(nextt)]).Log()
                        thisEmass.Stats()
                        lEM.append(thisEmass)
                        lEMA.append(thisEmass.Avg)
                        lEMS.append(thisEmass.Std)
            self.C2_Stats.loc[:,'EffM'] = pa.Series(lEM,index=self.C2_Stats.index)
            self.C2_Stats.loc[:,'EffMAvg'] = pa.Series(lEMA,index=self.C2_Stats.index)
            self.C2_Stats.loc[:,'EffMStd'] = pa.Series(lEMS,index=self.C2_Stats.index)

    # def ReadAndBootAndEffM(self,DefWipe=False):
    #     if os.path.isfile(self.PickleFile) and not DefWipe:
    #         self.LoadPickle()
    #     else:
    #         self.ReadAndBoot()
    #         self.EffMass()

    def Write(self):
        def FixDictArray(thislist,ikey):
            outDict = OrderedDict()
            for ic,ilist in enumerate(thislist):
                outDict[ikey+str(ic+1)] = ilist
            return outDict

        if not hasattr(self,'show_cfgs'):
            self.show_cfgs = False

        outDict = ODNested()
        outDict['name'] = self.name
        outDict['kud'] = self.kud
        outDict['ks'] = self.ks
        outDict['nt'] = self.nt
        outDict['nxyz'] = self.nxyz

        for istream,incfg in zip(self.stream_list,self.ncfg_list):
            outDict['ncfg_'+istream] = incfg
        outDict['nstream'] = self.nstream
        outDict['nxsrc_avg'] = self.xsrcAvg
        outDict['nMeas'] = self.nMeas
        outDict['nboot'] = self.nboot
        outDict['tsrc'] = self.tsrc
        # outDict['xsrc'] = FixDictArray(self.xsrc,'xsrc')
        outDict['ism'] = self.ism
        outDict['jsm'] = self.jsm
        if self.MesOrBar == 'Meson':
            if self.MesEffMass:
                outDict['EffMassFun'] = 'Cosh_ratio_inverse'
            else:
                outDict['EffMassFun'] = 'log_of_ratio'

        excel_params = pa.Series(deepcopy(outDict))

        if 'boot' in self.C2_Fit_Stats.columns:
            # print 'debug'
            # print self.C2_Fit_Stats['boot']
            for (istate,imom),fitdict in self.C2_Fit_Stats['boot'].items():
                if 'Fit' in fitdict.Fit_Stats.columns:
                    outDict['Fits'][istate][imom] = fitdict.GetOutputDict()


        # for ststr,statedict in self.FitDict.iteritems():
        #     if len(statedict.keys()) == 0: continue
        #     for ip,pdict in statedict.iteritems():
        #         for ifit,fitdict in pdict.iteritems():
        #             for ipar,pardict in fitdict.Params.iteritems():
        #                 thisfmtflag = 'f'
        #                 if ipar == 'A_Ep': thisfmtflag = 'e'
        #                 outDict['Fits'][ststr][ip][ifit][ipar] = AvgStdToFormat(pardict.Avg,pardict.Std,frmtflag=thisfmtflag)
        #             outDict['Fits'][ststr][ip][ifit]['Chi^2_pdf'] = AvgStdToFormat(fitdict.Chi2DoF.Avg,fitdict.Chi2DoF.Std)
        #             outDict['Fits'][ststr][ip][ifit]['MaxIter'] = fitdict.MaxIter
        #             outDict['Fits'][ststr][ip][ifit]['Prec'] = fitdict.Prec

        # print 'DEBUG',self.name
        # print self.C2_Stats
        if 'boot' in self.C2_Stats.columns:
            for (ip,it),pdata in self.C2_Stats['boot'].items():
                outDict['C2boot'][ip][it] = AvgStdToFormat(pdata.Avg,pdata.Std,frmtflag='e')
        if 'Auto' in self.C2_Stats.columns:
            for (ip,it),pdata in self.C2_Stats['Auto'].items():
                outDict['C2Auto'][ip][it] = AvgStdToFormat(pdata.Avg,pdata.Std,frmtflag='e')
        if 'EffM' in self.C2_Stats.columns:
            for (ip,it),pdata in self.C2_Stats['EffM'].items():
                outDict['EffM'][ip][it] = AvgStdToFormat(pdata.Avg,pdata.Std,frmtflag='e')

        if self.show_cfgs and any(np.array(self.ncfg_list) > 0):
            for ((istream,iccfg),icfg),ixsrc_list in zip(iter(self.C2_cfgs['configs'].items()),self.C2_cfgs['xsrc_list'].values):
                outDict['cfglist'][istream][icfg] = 'n_xsrc_'+str(len(ixsrc_list))
            if 'tsrc_list' in self.C2_cfgs:
                for ((istream,iccfg),icfg),itsrc_list in zip(iter(self.C2_cfgs['configs'].items()),self.C2_cfgs['tsrc_list'].values):
                    for ictsrc,itsrc in enumerate(itsrc_list):
                        outDict['tsrc_list'][istream][icfg]['itsrc_'+str(ictsrc)] = str(itsrc)



        ## Write human readable file with data
        WriteXml(self.HumanFile,{'Results':outDict})

        ## TODO, include fitting in the excell file
        WriteExcel(self.ExcelFile,{'C2_results':deepcopy(self.C2_Stats)},cfglist=self.C2_cfgs,params=excel_params)

        if any(np.array(self.ncfg_list) > 0):
            for ((istream,iccfg),icfg),ixsrc_list in zip(iter(self.C2_cfgs['configs'].items()),self.C2_cfgs['xsrc_list'].values):
                outDict['cfglist'][istream][icfg] = 'n_xsrc_'+str(len(ixsrc_list))
            if 'tsrc_list' in self.C2_cfgs:
                for ((istream,iccfg),icfg),itsrc_list in zip(iter(self.C2_cfgs['configs'].items()),self.C2_cfgs['tsrc_list'].values):
                    for ictsrc,itsrc in enumerate(itsrc_list):
                        outDict['tsrc_list'][istream][icfg]['itsrc_'+str(ictsrc)] = str(itsrc)

        ## pickles rest of data for reading
        self.RemoveFuns()
        WritePickle(self.PickleFile,self.__dict__)
        self.GetFuns()
    ##################################

    def ReadAndWrite(self,WipeData=True,NoFit=False):
        self.Read()
        self.Bootstrap(WipeData=WipeData)
        self.EffMass()
        if NoFit:
            self.Fit()
        self.Write()


    def FixRead(self,readdict):
        readdict['Info']['outdir'] = self.Info['outdir']
        readdict['latparams'].outputdir = self.latparams.outputdir
        readdict['latparams'].PickleFile = self.latparams.PickleFile
        readdict['latparams'].HumanFile = self.latparams.HumanFile
        readdict['PickleFile'] = self.PickleFile
        readdict['ExcelFile'] = self.ExcelFile
        readdict['HumanFile'] = self.HumanFile
        readdict['CfgFile'] = self.CfgFile
        readdict['C2_Fit_Col_Names'] = self.C2_Fit_Col_Names
        if 'cfg_file_type' not in list(readdict.keys()):
            readdict['cfg_file_type'] = self.cfg_file_type
        return readdict

    def Write_Cfgs(self,file_type='PreDef',show_timer=True):
        if hasattr(self,'read_cfgs_from_file') and self.read_cfgs_from_file:
            return
        if file_type == 'PreDef':
            file_type = self.cfg_file_type
        if 'to_' not in file_type:
            file_type = 'to_'+file_type
        this_file = self.CfgFile+fmt_file_type(file_type.replace('to_',''))
        if show_timer: thistimer = Timer(name='Writing Cfgs '+this_file)
        if len(self.C2_cfgs.index) == 0:
            self.Read(show_timer=show_timer)
        # this_df,dump = self.Get_Flattened_Cfgs()
        this_meta = [self.pmom,self.nt]
        WriteWithMeta(self.C2_cfgs,this_meta,this_file,file_type=file_type)
        if show_timer: thistimer.Stop()

    def ChopDepVarList(self,this_df,read_meta,show_err=True):
        isin = True
        name_list = ['momentum_list']
        index = [slice(None),slice(None)]
        if self.nt != read_meta[1]:
            if show_err:
                print('\n missmatch in euclidian time, read_nt='+str(read_meta[1])+', setup_nt='+str(self.nt))
            return this_df,False
        mom_list = read_meta[0]
        for ic,(ithis_list,icheck_list) in enumerate([(self.pmom,mom_list)]):
            isin = isin and set(ithis_list).issubset(icheck_list)
            if isin:
                same_list = len(ithis_list) == len(icheck_list)
                if not same_list:
                    index_list = [np.where(np.array(icheck_list) == it)[0][0] for it in ithis_list]
                    index[ic] = index_list
            elif show_err:
                print('\n missmatch in ',name_list[ic])
                print(ithis_list)
                print(icheck_list)
                break
        if isin:
            def Chop_fun(val):
                return np.array(val)[tuple(index)].tolist()
            this_df.loc[:,'C2'] = this_df.loc[:,'C2'].apply(Chop_fun)
        return this_df,isin




    def Read_Cfgs(self,file_type='PreDef',show_timer=True,CheckCfgs=True):
        def CheckCfgList(cfg_data,read_meta):
            if cfg_data is None or len(cfg_data) == 0:
                return None,False
            if CheckCfgs and 'stream_cfgs' in self.C2_cfgs:
                second_class = pa.Series([cfg_data],index=['C2_cfgs'])
                cfg_bool = CheckClass(self.__dict__,second_class,['C2_cfgs'])
                # if show_timer and not this_bool:
                #     print '\n'+GetCfgMissmatch( np.array(self.C2_cfgs['stream_cfgs'].values),
                #                                 np.array(cfg_data['stream_cfgs'].values))
                #
            else:
                cfg_bool = True
            if cfg_bool:
                cfg_data,outbool = self.ChopDepVarList(cfg_data,read_meta)
                return cfg_data,outbool
            else:
                return None,False

        if file_type == 'PreDef':
            file_type = self.cfg_file_type
        if 'read_' not in file_type:
            file_type = 'read_'+file_type
        this_file = self.CfgFile+fmt_file_type(file_type.replace('read_',''))
        if not os.path.isfile(this_file):
            if '-def-' in this_file:
                this_file = this_file.replace('-def-','-*-')
                file_list = glob.glob(this_file)
                if len(file_list) == 0:
                    print('Cfg FNF: \n', this_file)
                    return False
                else:
                    this_file = file_list[0]
        if not os.path.isfile(this_file):
            if show_timer: print('cfg FNF: ' + this_file)
            return False
        if show_timer: thistimer = Timer(name='Reading Cfgs '+this_file)
        meta_data,cfg_data = ReadWithMeta(this_file,file_type=file_type)
        cfg_data,check_bool = CheckCfgList(cfg_data,meta_data)
        # check_bool = check_bool and CheckMeta(meta_data)
        # cfg_data = self.Format_Flattened_Cfgs(cfg_data)
        self.read_cfgs_from_file = check_bool
        if check_bool:
            self.ImportCfgList(cfg_data)
            if show_timer: thistimer = thistimer.Stop()
            self.ncfg_list = [icfg.size for istream,icfg in self.C2_cfgs['configs'].groupby(level='stream')]
            return True
        else:
            if show_timer: print('\n file incorrect,', end=' ')
            self.CfgFile = self.CfgFile.replace('.cfgs','.cfgs.bak')
            if os.path.isfile(this_file.replace('.cfgs','.cfgs.bak')):
                print(' found backup file')
                return self.Read_Cfgs()
            else:
                print(' reading unformatted')
                return False


    def LoadPickle(self,DefWipe=False,WipeData = True,CheckCfgs=False,CheckMom=True,NoFit=False):
        # print 'filename ' , self.PickleFile
        if os.path.isfile(self.PickleFile) and not DefWipe:
            # print 'Loading Pickle for ' , self.name
            loadeddict = ReadPickleWrap(self.PickleFile)
            loadeddict = self.FixRead(loadeddict)
            checklist = self.thischecklist
            if CheckCfgs: checklist = checklist + ['C2_cfgs']
            if CheckMom: checklist = checklist + ['checkmomlist']

            loadeddict['checkmomlist'] = [self.latparams.pformTOpstr(ip,nopref=True) for ip in loadeddict['C2_Stats'].index.get_level_values('momentum')]
            if CheckClass(self.__dict__,loadeddict,checklist):
                self.__dict__.update(loadeddict)
                if Wipe_All_Fits:
                    self.C2_Fit_Stats = pa.DataFrame()
                elif len(self.C2_Fit_Stats.values)> 0 and not isinstance(self.C2_Fit_Stats.values[0][0],sff.SetOfFitFuns):
                    print('Legacy file, fits need to be recomputed')
                    self.C2_Fit_Stats = pa.DataFrame()
                self.read_cfgs_from_file = True
            else:
                if os.path.isfile(self.PickleFile+'.bak'):
                    print('using backupfile for '+self.name)
                    self.PickleFile = self.PickleFile+'.bak'
                    self.HumanFile = self.HumanFile+'.bak'
                    self.ExcelFile = self.ExcelFile.replace('.xlsx','.bak.xlsx')
                    self.LoadPickle(DefWipe=DefWipe,WipeData=WipeData,CheckCfgs=CheckCfgs,CheckMom=CheckMom)
                    return
                print('Warning, file ' , self.HumanFile , ' has different parameters than this instance, reading and writing over file now')
                self.ReadAndWrite(WipeData=WipeData,NoFit=NoFit)
                # self.Write()
        else:
            self.ReadAndWrite(WipeData=WipeData,NoFit=NoFit)
        self.EffMass()
    ## Operator overloading

    #Too Much work, only works for one operator. Please use SetCustomName after doing all algebraic manipulation


    def UpdateName(self,opp,C2,C22):
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
        if isinstance(C22,TwoPointCorr) and isinstance(C2,TwoPointCorr):
            op1LL = C2.LegLab[1:-1] ## removes $ symbols, re adds at end
            op2LL = C22.LegLab[1:-1] ## removes $ symbols, re adds at end
            if C2.name == C22.name:
                outname = C2.name,GetOppSelf(opp)
                outLL = op1LL,' Comb ',GetLegSelf(oppLL)
            else:
                for iname,iname2 in zip(C2.name.split('_'),C22.name.split('_')):
                    if iname in iname2:
                        outname.append(iname2)
                    elif iname2 in iname:
                        outname.append(iname)
                    else:
                        outname.append(iname+opp+iname2)
                outname = '_'.join(outname)
                for iLegLab,iLegLab2 in zip(op1LL.split('\ '),op2LL.split('\ ')):
                    if iLegLab in iLegLab2:
                        outLL.append(iLegLab2)
                    elif iLegLab2 in iLegLab:
                        outLL.append(iLegLab)
                    else:
                        outLL.append(iLegLab+oppLL+iLegLab2)
                outLL = '_'.join(outLL)
        elif isinstance(C2,TwoPointCorr):
            op1LL = C2.LegLab[1:-1] ## removes $ symbols, re adds at end
            if isinstance(C22,int):
                outname  =  C2.name,opp+str(C22)
                outLL  =  op1LL,oppLL+str(C22)
            else:
                if C22 > 0.01:
                    outname = C2.name,opp+'{:.2f}'.format(C22)
                    outLL = op1LL,oppLL+'{:.2f}'.format(C22)
        else:
            op2LL = C22.LegLab[1:-1] ## removes $ symbols, re adds at end
            if isinstance(C2,int):
                outname  =  str(C2) + opp ,C22.name
                outLL  =  str(C2) + oppLL,op2LL
            else:
                if C2 > 0.01:
                    outname = '{:.2f}'.format(C2) + opp,C22.name
                    outLL = '{:.2f}'.format(C2) + oppLL,op2LL
        # print r'$'+r'\ '.join(outLL)+'$'
        outLL = list(outLL)
        for ic,iLL in enumerate(outLL):
            if '^{}' in iLL:
                outLL[ic] = iLL.replace('^{}','^{')+'}'
        self.SetCustomName('_'.join(outname),r'$'+r'\ '.join(outLL)+'$')


    def Overload_wrap(self,C22,this_fun):
        result = TwoPointCorr(thisnboot=self.nboot, cfglist=self.C2_cfgs,
                              Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.SetCustomName(string=self.filename,stringLL=self.LegLab,fo_for_cfgs=self.fo_for_cfgs)
        # result.UpdateName('+',self,C22)
        if isinstance(C22,TwoPointCorr):
            result.C2_Stats.loc[:,'boot'] = this_fun(self.C2_Stats.ix[:, 'boot'] ,C22.C2_Stats.ix[:, 'boot'])
        else:
            try:
                result.C2_Stats.loc[:,'boot'] = this_fun(self.C2_Stats.loc[:,'boot'], C22)
            except Exception as e:
                print(type(C22), C22)
                print(self.C2_Stats)
                print(e)
                raise EnvironmentError('Invalid value to combine with TwoPointCorr class')
        return result


    def __add__(self,C22):
        if any([ipar in str(type(C22)) for ipar in TwoPointCorr.parentlist]):
            return C22.__radd__(self)
        else:
            return self.Overload_wrap(C22,op.add)

    def __sub__(self,C22):
        if any([ipar in str(type(C22)) for ipar in TwoPointCorr.parentlist]):
            return C22.__rsub__(self)
        else:
            return self.Overload_wrap(C22,op.sub)

    def __mul__(self,C22):
        if any([ipar in str(type(C22)) for ipar in TwoPointCorr.parentlist]):
            return C22.__rmul__(self)
        else:
            return self.Overload_wrap(C22,op.mul)

    def __div__(self,C22):
        if any([ipar in str(type(C22)) for ipar in TwoPointCorr.parentlist]):
            return C22.__rdiv__(self)
        else:
            return self.Overload_wrap(C22,op.truediv)

    def __pow__(self,C22):
        if any([ipar in str(type(C22)) for ipar in TwoPointCorr.parentlist]):
            return C22.__rpow__(self)
        else:
            return self.Overload_wrap(C22,op.pow)


    ## Right multiplication functions


    def __radd__(self,C22):
        if any([ipar in str(type(C22)) for ipar in TwoPointCorr.parentlist]):
            return C22.__add__(self)
        else:
            return self.Overload_wrap(C22,flip_fun_2arg(op.add))

    def __rsub__(self,C22):
        if any([ipar in str(type(C22)) for ipar in TwoPointCorr.parentlist]):
            return C22.__sub__(self)
        else:
            return self.Overload_wrap(C22,flip_fun_2arg(op.sub))

    def __rmul__(self,C22):
        if any([ipar in str(type(C22)) for ipar in TwoPointCorr.parentlist]):
            return C22.__mul__(self)
        else:
            return self.Overload_wrap(C22,flip_fun_2arg(op.mul))

    def __rdiv__(self,C22):
        if any([ipar in str(type(C22)) for ipar in TwoPointCorr.parentlist]):
            return C22.__div__(self)
        else:
            return self.Overload_wrap(C22,flip_fun_2arg(op.truediv))


    def __rpow__(self,C22):
        if any([ipar in str(type(C22)) for ipar in TwoPointCorr.parentlist]):
            return C22.__pow__(self)
        else:
            return self.Overload_wrap(C22,flip_fun_2arg(op.pow))

    ## Numerical operatoions

    ## unitary arithmetic operations

    def __neg__(self):
        for ival in list(self.keys()):
            self[ival] = -self[ival]
        return self


    def __pos__(self):
        return self

    def __abs__(self):
        for ival in list(self.keys()):
            self[ival] = abs(self[ival])
        return self

    ## Casting, just casts the bootstrap values, ignore return value (only needed since return has to be same type).


    def __complex__(self):
        for ival in list(self.keys()):
            complex(self[ival])
        self.Stats()
        return 1+0j

    def __int__(self):
        for ival in list(self.keys()):
            int(self[ival])
        self.Stats()
        return 1

    def __float__(self):
        for ival in list(self.keys()):
            np.float64(self[ival])
        self.Stats()
        return 1.0

####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################


class NNQCorr(object):
    """ NNQ(p,t,tf) uses bootstrap class
        p  = momentum
        t  = time
        tf = tflow
    """
    def Construct_NNQ_File(this_file):
        return Construct_File_Object(this_file,NNQCorr)

    ## Info is a dictionary containing information for the correlator:
    ## see comments below to see what is required.
    ## missing elements just ignores that field
    ## All conflicts between TwoPointCorr and FlowOp were resolved by adding flow to variable names, and Flow to function names.
    ## All new functions and variables where I do not want to overwrite TwoPointCorr have the keyword NNQ and nnq respectivly.

    ## TODO: create variational flag and interface with this function.

    ## array of all classes that inherit this class. Append this when adding new classes please!
    parentlist = ['SetOfCorrs.SetOfNNQ','TwoPtCorrelators.NNQFullCorr']

    def __init__(self, thisnboot=nboot, cfglist={},
                 Info={},thissym='Not Set',thiscol='Not Set',thisshift=0.0,name='',man_load_cfgs=False):
        ## initialises all the data dictionaries


        ## Setting to g5P+NNQ (or g5P+NNW) by forceing Info to read correct projector and BarOrMes

        # if 'Interp' not in Info.keys() or 'g5' not in Info['Interp']:
        #     Info['Interp'] = 'CPOdd'
        Info['MesOrBar'] = 'Baryon'
        ## Instance of FlowOp
        self.FO = FlowOp(thisnboot,cfglist,Info,thissym,thiscol,thisshift,name)
        # ## Instance of TwoPointCorr
        # Info['fo_for_cfgs'] = self.FO.Observ

        self.NN = TwoPointCorr(thisnboot,cfglist,Info,thissym,thiscol,thisshift,name)

        ## these class elemnts are checked when loading in pickled file
        self.thischecklist = ['nboot','nt','ProjStr','ism','jsm','tsrc','kud','ks','NNfile','FOfile','tflowlist']

        self.checkmomlist = self.NN.checkmomlist

        ## pull stuff from the NN and FO instances
        self.IsComb = self.IsSubComb = False ## if it is a combination of correlators
        self.NNfile=self.NN.PickleFile
        self.FOfile=self.FO.flowPickleFile
        self.tflowlist = self.FO.tflowlist
        self.ProjStr = self.NN.ProjStr
        self.ism = self.NN.ism
        self.tsrc = self.NN.tsrc
        self.kud = self.NN.kud
        self.ks = self.NN.ks

        ## Grab some stuff (from .NN and .FO)
        self.stream_list = self.NN.stream_list
        self.nboot = thisnboot
        self.nt = self.NN.nt
        self.nxyz = self.NN.nxyz
        self.CPEvenFile = 'None'
        self.Info = Info
        self.thissym = thissym
        self.thiscol = thiscol
        self.thisshift = thisshift
        self.thisshiftscale = 'Not Set'
        self.name = name

        self.jsm = self.NN.jsm

        ## used to pick out different configuraiton/source lists.
        ## for example, if there are less flowed operators than 2 point correlators, passing 'TopCharge' will read that result
        if 'fo_for_cfgs' in list(Info.keys()): self.fo_for_cfgs = Info['fo_for_cfgs']
        else: self.fo_for_cfgs = False


        if 'Sparam' in list(Info.keys()): self.Sparam = Info['Sparam']
        else: self.Sparam = defSparam


        if 'show_cfgs' in list(Info.keys()): self.show_cfgs = Info['show_cfgs']
        else: self.show_cfgs = False


        ## Info can contain initial guess for fits, make sure it array of correct length
        if 'AlphaiGuess' in list(Info.keys()): self.AlphaiGuess = np.array(Info['AlphaiGuess'])
        else: self.AlphaiGuess = 'None' ## Defaults to value set for function in FitFunctions.py

        ## Info can contain initial guess for fits, make sure it array of correct length
        if 'NNQiGuess' in list(Info.keys()): self.NNQiGuess = np.array(Info['NNQiGuess'])
        else: self.NNQiGuess = 'None' ## Defaults to value set for function in FitFunctions.py


        ## Info can contain fit ranges for alpha to be calculated
        if 'FitsAlpha' in list(Info.keys()): self.PredefFits = Info['FitsAlpha']
        else: self.PredefFits = []

        ## Info can contain default tflow picked, used for plotting
        if 'tflowfit' in list(Info.keys()): self.SetTflowFit(Info['tflowfit'])
        else: self.SetTflowFit('Not Set')


        ## Info can contain default t sink picked, used for plotting
        if 'tsinkfit' in list(Info.keys()): self.SetTFit(Info['tsinkfit'])
        else: self.SetTFit('Not Set')

        ## Info can contain fit ranges to be calculated
        ## must bein form Info['FitFunAlpha'] = (Function Object , number of parameters )
        if 'FitFunAlpha' in list(Info.keys()):
            thisFun,thisnpar = Info['FitFunAlpha']
            self.SetFunction(thisFun,thisnpar)
        else:
            self.SetFunction('None',0)

        ## selects the improved overlap, where the noise is subtracted from Q in NNQ
        if 'IExpect' in list(Info.keys()):
            self.IExpect = Info['IExpect']
        else:
            self.IExpect = False


        # self.NNQOp = ODNested()
        ## Series [ momentum , source_sink_sep , flow_time , configurations]
        self.NNQ_cfgs = pa.DataFrame()
        '''
        self.NNQ_cfgs DataFrame:

        columns =   xsrc_list , configs , C2, Op, NNQ

        rows =      istream, iccfg
        '''

        self.NNQ_col_names = ['momentum','source_sink_separation','flow_time']
        self.NNQ_Stats = pa.DataFrame()

        '''
        self.NNQ_Stats DataFrame:

        columns =   NNQboot,    NNQAvg,         NNQStd
                    EffM,       EffMAvg,        EffMStd
                    Alphaboot,  AlphaAvg,       AlphaStd
                    AlphaAuto,  AlphaAutoAvg,   AlphaAutoStd
        rows =      momentum,   source_sink_sep,    flow_time
        '''




        self.NNQ_Fit_Stats = pa.DataFrame()
        self.NNQ_Fit_col_names = ['momentum','flow_time']
        '''
        self.NNQ_Stats DataFrame:

        columns =   NNQboot,    NNQAvg,         NNQStd
                    boot

        rows =      momentum,   flow_time,  fit_range
        boot elements are SetsOfFitFuns
        '''

        # self.LoadPickle(DefWipe=DefWipe)
        # self.Read()
        self.SetCustomName(name)
        if not man_load_cfgs:
            self.LoadCfgs(cfglist,name)

    def LoadCfgs(self,cfglist,name=''):
        if len(cfglist) == 0:
            self.GetCfgList()
        else:
            self.ImportCfgList(cfglist)
        self.SetCustomName(name)

    ## to be consistant with TwoPtCorr;
    def GetCfgList(self):
        self.CombineCfgLists()

    def ImportPlotParams(self,thiscol,thissym,thisshift):
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        self.thisshift = thisshift


    def CombineCfgLists(self):
        if 'configs' not in self.FO.Op_cfgs.columns:
            if 'configs' in self.NN.C2_cfgs.columns:
                self.ImportCfgList(self.NN.C2_cfgs)
            return
        if 'configs' not in self.NN.C2_cfgs.columns:
            self.ImportCfgList(self.FO.Op_cfgs['configs'])
            return
        self.ImportCfgList(CombineCfgs(self.FO.Op_cfgs,self.NN.C2_cfgs))
        # FO_cfgs = self.FO.Op_cfgs['configs']
        # CombineCfgs
        # comb_cfgs,comb_xsrc = [],[]
        # ilist = []
        # for (istream,iccfg),NN_cfg_data in self.NN.C2_cfgs.iterrows():
        #     if NN_cfg_data['configs'] in FO_cfgs[istream].values:
        #         comb_cfgs.append(NN_cfg_data['configs'])
        #         comb_xsrc.append(NN_cfg_data['xsrc_list'])
        #         ilist.append((istream,iccfg))
        # indicies = pa.MultiIndex.from_tuples(ilist,names=['stream','config_number'])
        # cfglist = pa.Series(comb_cfgs,index=indicies).to_frame(name='configs')
        # cfglist.loc[:,'xsrc_list'] = pa.Series(comb_xsrc,index=indicies)

    def Stats(self,DefWipeG2=False,DoAuto=True):
        def DoStats(val):
            if isinstance(val,BootStrap): val.Stats()
            return val
        def DoAvg(val):
            if isinstance(val,BootStrap):
                return val.Avg
            else:
                return val
        def DoStd(val):
            if isinstance(val,BootStrap):
                return val.Std
            else:
                return val
        if not self.IsComb:
            self.NN.Stats()
            self.FO.FlowStats()
        ## maybe effmass and alpha not here???
        self.NNQ_Stats.loc[:,'boot'] = self.NNQ_Stats['boot'].apply(DoStats)
        self.NNQ_Stats.loc[:,'Avg'] = self.NNQ_Stats['boot'].apply(DoAvg)
        self.NNQ_Stats.loc[:,'Std'] = self.NNQ_Stats['boot'].apply(DoStd)
        if 'Alphaboot' in self.NNQ_Stats:
            self.NNQ_Stats.loc[:,'Alphaboot'] = self.NNQ_Stats['Alphaboot'].apply(DoStats)
            self.NNQ_Stats.loc[:,'AlphaAvg'] = self.NNQ_Stats['Alphaboot'].apply(DoAvg)
            self.NNQ_Stats.loc[:,'AlphaStd'] = self.NNQ_Stats['Alphaboot'].apply(DoStd)
        self.EffMass()
        if not hasattr(self,'IsSubComb'):
            self.IsSubComb = False
        if not self.IsSubComb:
            self.AlphaRatio(DefWipeG2=DefWipeG2,DoAuto=DoAuto)

    ## not needed, taken care of when bootstrapping
    # def NNQCombPreBS(self):
    #     self.NNQ = np.array([np.multiply.outer(iC2,Op) for iC2,iOp in zip(self.C2,self.Op)])

    ## DefWipe can be 'C2' or 'Flow' to only delete that specific data
    ## You can even call self.FlowCombinePreBS(OtherOpp) if you wanted to look at other flowed operators too!
    def Bootstrap(self,WipeData = False,tflowlist = 'PreDef',show_timer=True):
        # self.C2 = [ ip , it , iconf ]
        # self.Op = [ itflow , iconf ]
        if 'boot' in self.NNQ_Stats:
            self.RemoveVals(WipeData = WipeData)
            return
        if tflowlist != 'PreDef': self.tflowlist = tflowlist
        ## Do not wipe the data, we need it for combining on the configuration level!
        self.NN.Bootstrap(WipeData=False)
        self.FO.FlowBootstrap(tflowlist=tflowlist,WipeData=False,show_timer=False)
        # ReadAll = (len(self.tflowlist) == 0)
        NNQlcfg = []
        NNQl,NNQlAvg,NNQlStd = [],[],[]

        indexl = []
        if show_timer: thistimer = Timer(linklist=[self.NN.pform,list(range(self.nt)),self.FO.tflowfloat],name='Booting '+self.name)
        for icp,ip in enumerate(self.NN.pform):
            NNQlcfg.append([])
            for it in range(self.nt):
                NNQlcfg[-1].append([])
                tdata = self.NNQ_cfgs['C2'].apply(lambda x : x[icp][it])
                strit = tstr(it+1)
                for ictflow,tflows in enumerate(self.FO.tflowfloat):
                    tflowdata = self.NNQ_cfgs['Op'].apply(lambda x : x[ictflow])
                    strtflows = tflowstr(tflows)
                    # print self.tflowlist , strtflows
                    # print 'test',ip,strit,strtflows
                    thisdata = tdata * tflowdata
                    indexl.append((ip,strit,strtflows))
                    # for icfg,icfg_data in zip(self.FO.flowcfglist,thisdata):
                    NNQlcfg[-1][-1].append(thisdata.values)

                    if hasattr(self,'IExpect') and self.IExpect:
                        raise EnvironmentError('IExpect not implemented after pandas rework. on the TODO list!')
                        # if isinstance(self.IExpect,float):
                        #     self.NNQboot[ip][strit][strtflows] = BootStrap(self.nboot, name='G_{2}('+','.join([ip,strit,strtflows])+')',cfgvals=tdata,
                        #                                                    thiseps=self.IExpect,QforImprov=tflowdata)
                        #     # thisfile = './temp/'+'_'([ip,strit,strtflows,epsstr(self.IExpect)])+'.py3p'
                        #     # if os.path.isfile(thisfile):
                        #     #     with open(thisfile,'rb') as f:
                        #     #         self.NNQboot[ip][strit][strtflows] = pickle.load(f)
                        #     # else:
                        #     #     for ieps in np.arange(0.1,self.IExpect,0.1):
                        #     #         mkdir_p('./temp')
                        #     #         with open(thisfile.replace(epsstr(self.IExpect),epsstr(ieps)),'wb') as f:
                        #     #             pickle.dump(BootStrap(self.nboot, name='G_{2}('+','.join([ip,strit,strtflows])+')',data=tdata,
                        #     #                                                thiseps=ieps,QforImprov=tflowdata),f)
                        # elif self.IExpect == True:
                        #     self.NNQboot[ip][strit][strtflows] = BootStrap(self.nboot, name='G_{2}('+','.join([ip,strit,strtflows])+')',cfgvals=tdata,QforImprov=tflowdata)
                        # elif self.IExpect == False:
                        #     self.NNQboot[ip][strit][strtflows] = BootStrap(self.nboot, name='G_{2}('+','.join([ip,strit,strtflows])+')',cfgvals=thisdata)
                    else:
                        thisboot = BootStrap(self.nboot, name='G_{2}('+','.join([ip,strit,strtflows])+')',cfgvals=thisdata)
                        NNQl.append(thisboot)
                    NNQlAvg.append(thisboot.Avg)
                    NNQlStd.append(thisboot.Std)
                    # thistimer.Lap()
                    if show_timer: thistimer.Lap((ip,strit,strtflows))
        indicies = pa.MultiIndex.from_tuples(indexl,names=self.NNQ_col_names)
        try:
            self.NNQ_cfgs.loc[:,'NNQ'] = pa.Series(np.rollaxis(np.array(NNQlcfg),-1).tolist(),index=self.NNQ_cfgs.index)
        except Exception as err:
            out_str = '\n missmatch in configuration sizes:\n'
            out_str += '    number of configs in NNQ:'+str(len(self.NNQ_cfgs))+'\n'
            out_str += '    number of configs passed into NNQ:'+str(len(np.rollaxis(np.array(NNQlcfg),-1).tolist()))
            raise EnvironmentError(str(err) + out_str)
        self.NNQ_Stats.loc[:,'boot'] = pa.Series(NNQl,index=indicies)
        self.NNQ_Stats.loc[:,'Avg'] = pa.Series(NNQlAvg,index=indicies)
        self.NNQ_Stats.loc[:,'Std'] = pa.Series(NNQlStd,index=indicies)
        self.RemoveVals(WipeData = WipeData)

    def ObsToName(self,thisop):
        thisobs = self.FO.FlowObsToName(thisop)
        return thisobs.replace('<',r'\alpha_{').replace('>',r'}')

    # ## thisInfo must contain correct information for combining the correlators (the smearings, the eigenvectors, etc...)
    # def GetCombCorr(self,thisInfo):
    #     from SetsOfTwoPtCorrs import SetOfNNQ
    #     sinklab,jsmcoeffs,sm2ptList,InfoList = GetInfoList(thisInfo)

    #     C2Set = SetOfNNQ(cfglist=self.NNQ_cfgs[['configs','xsrc_list']],InfoDict = InfoList)
    #     C2Set.LoadPickle(CheckCfgs=CheckCfgs)
    #     thisname,thisleglab,jsmlab = CreateNameAndLL(C2Set.SetC2,sinklab)
    #     thisfun = GenSinkFun(jsmcoeffs)
    #     thisfile = self.filename.replace(self.jsm,sinklab)
    #     thisLL = self.LegLab.replace(self.jsm,sinklab)
    #     # self.SetCustomName(thisfile,thisLL)
    #     # self.jsm = jsmlab
    #     return C2Set.CombineCorrs(thisfun,filename=thisfile,LegLab=thisLL,jsmlab=jsmlab)


    def SetCustomName(self,string='',stringLL='',fo_for_cfgs='PreDef'):
        if string == '':
            self.name = '_'.join([self.NN.dim_label,self.NN.kappafolder,self.NN.stream_name,self.NN.tsrc,self.NN.ism,self.jsm,self.FO.Observ])
            self.filename = '_'.join([self.NN.MesOrBar,self.NN.stream_name,self.ProjStr,self.NN.tsrc,self.NN.ism,self.jsm,self.FO.Observ])
            if hasattr(self,'IExpect'):
                if isinstance(self.IExpect,float):
                    self.filename += '_Ieps'+epsstr(self.IExpect)
                elif self.IExpect:
                    self.filename += '_Ieps'+epsstr(Qeps)
        else:
            self.filename = self.name = string
        if stringLL == '':
            self.LegLab = '$'+'\ '.join([self.NN.dim_ll,self.NN.latparams.GetPionMassLab(),self.NN.stream_name,self.ProjStr,self.NN.ism,self.jsm,self.ObsToName(self.FO.Observ)])+'$' ## customise this how you like
            if hasattr(self,'IExpect'):
                if isinstance(self.IExpect,float):
                    self.LegLab += ' Ieps'+epsstr(self.IExpect)
                elif self.IExpect:
                    self.LegLab += ' Ieps'+epsstr(Qeps)
        else:
            self.LegLab = stringLL
        if isinstance(fo_for_cfgs,str):
            if fo_for_cfgs == 'PreDef':
                fo_for_cfgs = self.fo_for_cfgs
        if isinstance(fo_for_cfgs,str):
            self.NNQdir = outputdir + '/'+self.NN.dim_label+self.NN.kappafolder+'/alpha/'+fo_for_cfgs+'/'
        else:
            self.NNQdir = outputdir + '/'+self.NN.dim_label+self.NN.kappafolder+'/alpha/'
        mkdir_p(self.NNQdir+'/Pickle/')
        mkdir_p(self.NNQdir+'/Excel/')
        self.HumanFile = self.NNQdir+self.filename+'.xml'
        self.ExcelFile = self.NNQdir+'/Excel/'+self.filename+'.xlsx'
        self.PickleFile = self.NNQdir+'/Pickle/'+self.filename+'.py3p'


    def ImportCfgList(self,cfglist):
        self.FO.FlowImportCfgList(cfglist)
        self.NN.ImportCfgList(cfglist)
        self.NNQ_cfgs = self.NN.C2_cfgs
        # self.NNQ_cfgs['Op'] = pa.Series(self.FO.Op_cfgs.loc[:,'Op'].values,index=self.NNQ_cfgs.index)
        self.Update_Stream_List()

    def Update_Stream_List(self):
        self.stream_list = self.NN.stream_list
        self.ncfg_list = self.NN.ncfg_list
        self.xsrcAvg = self.NN.xsrcAvg
        self.nMeas = self.NN.nMeas
        self.nstream = len(self.stream_list)
        if hasattr(self,'stream_name'):
            old_sn = self.stream_name
            self.stream_name = '-'+''.join([istr.replace('-','') for istr in self.stream_list])+'-'
            if hasattr(self,'filename'):
                self.filename = self.filename.replace(old_sn,self.stream_name)
            if hasattr(self,'name'):
                self.name = self.name.replace(old_sn,self.stream_name)
            if hasattr(self,'LegLab'):
                self.LegLab = self.LegLab.replace(old_sn,self.stream_name)
        else:
            self.stream_name = '-'+''.join([istr.replace('-','') for istr in self.stream_list])+'-'



    def RemoveVals(self,WipeData = True):
        if 'C2' == WipeData:
            self.NN.RemoveVals()
        elif 'Flow' ==  WipeData:
            self.FO.FlowRemoveVals()
        elif WipeData:
            self.RemoveVals('C2')
            self.RemoveVals('Flow')

    ## self.C2 = [ ip, it , iconf ]
    def Read(self):
        self.NN.Read()
        self.FO.FlowRead()
        self.NNQ_cfgs = self.NN.C2_cfgs
        self.NNQ_cfgs.loc[:,'Op'] = pa.Series(self.FO.Op_cfgs.loc[:,'Op'].values,index=self.NNQ_cfgs.index)

    # def ReadAndBoot(self):
    #     self.Read()
    #     self.Bootstrap()

    def ImportFitRanges(self,fit_range):
        self.fit_range = fit_range
        if isinstance(self.fit_range,(list,tuple,np.ndarray)):
            print('SetsOfFits implementation only takes 1 fit range to vary over, choosing first')
            self.fit_range = self.fit_range[0]


    def CheckCol(self,thiscol):
        if 'PreDefine' == thiscol :
            if 'Not Set' == self.thiscol:
                raise IOError('Pass in color to initialize it')
        else:
            self.thiscol = thiscol

    def CheckSym(self,thissym):
        if 'PreDefine' == thissym :
            if 'Not Set' == self.thissym:
                raise IOError('Pass in symbol to initialize it')
        else:
            self.thissym = thissym

    def GetShift(self,xlims,thisshift):
        if thisshift != 'PreDefine': self.thisshift = thisshift
        # if self.thisshiftscale == 'Not Set':
        xlen = np.abs(xlims[1]-xlims[0])
        self.thisshiftscale = self.thisshift*xlen
        return self.thisshiftscale
        # else:
        #     return self.thisshiftscale

    def SetTflowFit(self,tflow):
        # if tflow not in self.tflowlist and tflow != 'Not Set':
        #     print 'Warning: ', tflow, ' not in flow list, tflowfit not added'
        # else:
        self.tflowfit = tflow

    def SetTFit(self,tsink):
        # if tsink not in ['t'+str(it) for it in xrange(1,self.nt+1)] and tsink != 'Not Set':
        #     print 'Warning: ', tsink, ' not in  list, tsinkfit not added'
        # else:
        self.tsinkfit = tsink

    def SetFunction(self,Fun,npar):
        self.Fun = Fun
        self.npar = npar

    def GetFuns(self):
        # self.NN.GetFuns()
        if self.NNQ_Fit_Stats is not None and 'boot' in self.NNQ_Fit_Stats.columns:
            for ival in self.NNQ_Fit_Stats['boot'].values:
                ival.GetFuns()
        self.FO.GetFuns()
        if not hasattr(self,'Fun'):
            self.Fun = ReadFuns(self.Fun_name)[0]
        if hasattr(self,'NNCPEven'):
            self.NNCPEven.GetFuns()

    def RemoveFuns(self):
        # self.NN.RemoveFuns()
        if self.NNQ_Fit_Stats is not None and 'boot' in self.NNQ_Fit_Stats.columns:
            for ival in self.NNQ_Fit_Stats['boot'].values:
                ival.RemoveFuns()
        self.FO.RemoveFuns()
        if hasattr(self,'Fun'):
            self.Fun_name = self.Fun.__name__
            WriteFuns(self.Fun)
            del self.Fun
        if hasattr(self,'NNCPEven'):
            self.NNCPEven.RemoveFuns()

    def FitNNQ(self,NNfit='First',fit_range='PreDef',iGuess='PreDef',EstDir=False,thistflowlist = 'PreDef',WipeFit=True):
        print('NNQ fit not implemented for setsoffits, skipping')
        # if self.Fun == 'None': raise IOError('Please define funciton using .SetFunction(Fun,npar) for fitting before calling .FitNNQ()')
        # if 'PreDef' == fit_range:
        #     self.fit_range = self.PredefFits
        # else:
        #     self.fit_range = fit_range
        # if iGuess != 'PreDef': self.NNQiGuess = iGuess
        # if thistflowlist == 'PreDef':
        #     if self.tflowfit == 'Not Set': raise IOError('tflowfit not set, please set by using .SetTflowFit()')
        # else:
        #     if thistflowlist == 'All':
        #         self.tflowfit = self.tflowlist
        #     else:
        #         self.tflowfit = thistflowlist
        #
        # lFit,lFitA,lFitS = [],[],[]
        # ilist = []
        # for ip,pdata in self.NNQ_Stats['boot'].groupby(level='momentum'):
        #     if isinstance(NNfit,(list,tuple,np.ndarray)):
        #         NNfit = NNfit[0]
        #     indexlist = self.NNCPEven.C2_Fit_Stats['boot'].index.levels
        #     if NNfit == 'First':
        #         if 'state1' in indexlist[0] and ip in indexlist[1] and len(self.NNCPEven.C2_Fit_Stats['boot'].loc['state1',ip,:]) == 0:
        #             self.NNCPEven.Fit(state='state1')
        #             if 'state1' in indexlist[0] and ip in indexlist[1] and len(self.NNCPEven.C2_Fit_Stats['boot'].loc['state1',ip,:]) == 0:
        #                 raise EnvironmentError('NN fit dictonary has not fits, please run self.NN.Fit()')
        #         thisindex = self.NNCPEven.C2_Fit_Stats['boot'].loc['state1',ip,:].index[0]
        #     else:
        #         thisindex = ('state1',ip,NNfit)
        #         if thisindex not in self.NNCPEven.C2_Fit_Stats['boot']:
        #             self.NNCPEven.Fit(state='state1',fit_range=NNfit)
        #     NNfit = self.NNCPEven.C2_Fit_Stats['boot'][thisindex]
        #     if isinstance(NNfit,pa.Series):
        #         NNfit = NNfit.iloc[0]
        #     NNparams = NNfit.fit_data['Params']
        #     NNmass = NNparams['Energy']
        #     NNCoeff = NNparams['A_Ep']
        #
        #
        #     # print 'debug', NNmass.Avg,NNCoeff.Avg
        #     for itflow in self.tflowfit:
        #         if (ip,itflow) not in pdata.swaplevel(1,2): continue
        #         for ifit in self.fit_range:
        #             if (ip,itflow,ifit) in self.NNQ_Fit_Stats.index:
        #                 if WipeFit:
        #                     self.NNQ_Fit_Stats.drop([(ip,itflow,ifit)],inplace=True)
        #                 else:
        #                     continue
        #             ilist.append((ip,itflow,ifit))
        #             ifitmin,ifitmax = np.array(unxmlfitr(ifit),dtype=np.float64)
        #             tdatarange,ydata = [[]],[]
        #             for it,idata in pdata.swaplevel(1,2)[ip,itflow].iteritems():
        #                 tint = untstr(it)
        #                 if ifitmin <= tint <= ifitmax:
        #                     tdatarange[0].append(tint)
        #                     # print
        #                     # print 't, flowdata',tint,tflowdata.Avg
        #                     # print 'NNmass',NNmass.Avg
        #                     # print 'NNCoeff',NNCoeff.Avg
        #                     # print '(inverse) exp value Avg',np.exp(NNmass.Avg*tint)
        #                     # print 'denominator for ratio Avg',np.exp(NNmass.Avg*tint)/NNCoeff.Avg
        #                     # print '(inverse) exp value',np.mean((NNmass*tint).Exp().bootvals)
        #                     # print 'denominator for ratio',np.mean(((NNmass*tint).Exp()/NNCoeff).bootvals)
        #                     thisdata = idata*(NNmass*tint).Exp()/NNCoeff
        #                     thisdata.Stats()
        #                     ydata.append(thisdata)
        #                     # print 'result value',np.mean(ydata[-1].bootvals)
        #                     # print
        #
        #             # ##debugging
        #             # print
        #             # print 'FitData:'
        #             # for it,iy in zip(tdatarange[0],ydata):
        #             #     iy.Stats()
        #             #     print it, iy.Avg
        #             # print
        #
        #             if EstDir:
        #                 thisff = ff.Fitting(Funs=[self.Fun,self.npar,'Estimate'],data=[np.array(tdatarange),np.array(ydata)],name='NNQ'+ifit,iGuess=self.NNQiGuess)
        #             else:
        #                 thisff = ff.Fitting(Funs=[self.Fun,self.npar],data=[np.array(tdatarange),np.array(ydata)],name='NNQ'+ifit,iGuess=self.NNQiGuess)
        #
        #             thistest = thisff.FitBoots()
        #             if not np.isnan(thistest):
        #                 ilist.append((ip,itflow,ifit))
        #                 lFit.append(thisff)
        #                 lFitA.append(thisff.fit_data['ParamsAvg'])
        #                 lFitS.append(thisff.fit_data['ParamsStd'])
        #                 # print 'results after fit'
        #                 # for iparam,paramval in thisff.Params.iteritems():
        #                 #     print iparam,paramval.Avg
        #                 # print
        # if len(lFit) != 0:
        #     indicies = pa.MultiIndex.from_tuples(ilist,names=self.NNQ_Fit_col_names)
        #     if 'NNQboot' in self.NNQ_Fit_Stats.columns:
        #         this_df = pa.DataFrame()
        #         this_df.loc[:,'NNQboot'] = pa.Series(lFit,index=indicies)
        #         this_df.loc[:,'NNQAvg'] = pa.Series(lFitA,index=indicies)
        #         this_df.loc[:,'NNQStd'] = pa.Series(lFitS,index=indicies)
        #         self.NNQ_Fit_Stats = self.NNQ_Fit_Stats.append(this_df)
        #     else:
        #         self.NNQ_Fit_Stats.loc[:,'NNQboot'] = pa.Series(lFit,index=indicies)
        #         self.NNQ_Fit_Stats.loc[:,'NNQAvg'] = pa.Series(lFitA,index=indicies)
        #         self.NNQ_Fit_Stats.loc[:,'NNQStd'] = pa.Series(lFitS,index=indicies)
        #     self.Write()


    def FitAlpha(   self,fit_range='PreDef',iGuess='PreDef',
                    EstDir=False,thistflowlist = 'PreDef',WipeFit=True,show_timer=False):
        if self.Fun == 'None':
            raise IOError('Please define funciton using .SetFunction(Fun,npar) for fitting before calling .FitAlpha()')
        if 'PreDef' == fit_range:
            self.fit_range = self.PredefFits
        else:
            self.fit_range = fit_range
        if not isinstance(self.fit_range,str):
            print('Fit ranges for setsoffits implementation now only takes 1 fit range to scan with.')
            print('selecting first: ',self.fit_range[0])
            self.fit_range = self.fit_range[0]
        else:
            print('Fitting Alpha scanned in ',str(self.fit_range))
        if iGuess != 'PreDef': self.AlphaiGuess = iGuess
        if thistflowlist == 'PreDef':
            if self.tflowfit == 'Not Set':
                raise IOError('tflowfit not set, please set by using .SetTflowFit()')
        else:
            if thistflowlist == 'All':
                self.tflowfit = self.tflowlist
            else:
                self.tflowfit = thistflowlist
        lFit,ilist = [],[]
        for ip_key,pdata in self.NNQ_Stats['Alphaboot'].groupby(level=('momentum','flow_time')):
            if ip_key[1] not in self.tflowfit:
                print(ip_key[1], 'not in flow list, skipping')
                # print('DEBUG',self.tflowlist)
                continue
            ## checks and sets for if data is already present.
            if ip_key in self.NNQ_Fit_Stats.index:
                this_sof = self.NNQ_Fit_Stats.loc[ip_key,'boot']
                if isinstance(this_sof,pa.Series):
                    this_sof = this_sof.iloc[0]
                if EstDir:
                    this_fun = (self.Fun,self.npar,'Estimate')
                else:
                    this_fun = (self.Fun,self.npar)
                if this_sof.CheckRange_fitr(self.fit_range,fit_info={'Funs':this_fun}):
                    if WipeFit:
                        self.NNQ_Fit_Stats.drop([ip_key],inplace=True)
                    else:
                        continue

            # ##debugging
            # print
            # print 'FitData:'
            # for it,iy in zip(tdatarange[0],ydata):
            #     print it, iy
            # print

            fit_info = {}
            if EstDir:
                fit_info['Funs'] = (self.Fun,self.npar,'Estimate')
            else:
                fit_info['Funs'] = (self.Fun,self.npar)
            fit_info['iGuess'] = self.AlphaiGuess
            fit_info['name'] = self.name + ' Fits'
            fit_info['paramlist'] = [r'\alpha']
            this_fit = sff.SetOfFitFuns(data=pdata.loc[ip_key[0],:,ip_key[1]])
            this_test = this_fit.ScanRange_fitr(self.fit_range,fit_info=fit_info)
            if this_test:
                ilist.append(ip_key)
                lFit.append(this_fit)
            else:
                print('Alpha fit for ', ip_key,' failed to setup')
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.NNQ_Fit_col_names)
            if 'boot' in self.NNQ_Fit_Stats.columns:
                this_df = pa.DataFrame()
                this_df.loc[:,'boot'] = pa.Series(lFit,index=indicies)
                self.NNQ_Fit_Stats = self.NNQ_Fit_Stats.append(this_df)
            else:
                self.NNQ_Fit_Stats.loc[:,'boot'] = pa.Series(lFit,index=indicies)

            def DoFit(val):
                val.DoFits(show_timer=show_timer)
                val.SortChi()
                return val
            self.NNQ_Fit_Stats.loc[:,'boot'] = self.NNQ_Fit_Stats['boot'].apply(DoFit)
            self.Write()
        else:
            print(str(self))
            print('Alpha found nothing to fit!')


    def Get_Extrapolation(self,keep_keys=slice(None),fmted=True):
        if 'boot' in self.NNQ_Fit_Stats:
            return sff.PullSOFSeries_Wpar(self.NNQ_Fit_Stats['boot'],keep_keys=keep_keys,fmted=fmted)
        else:
            return pa.Series()

    def AlphaPlotTauIntWopt(    self,plot_class,xlims=defxlimAlpha,
                                momlist=['p000'],thiscol='PreDefine',
                                thissym='PreDefine',thisshift='PreDefine',
                                thistflowlist = 'PreDef'):
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        self.thisshiftscale = 'Not Set'
        if thistflowlist == 'PreDef': thistflowlist = self.tflowfit
        if isinstance(xlims[0],str): xlims = list(map(untstr,xlims))
        thisshift = self.GetShift(xlims,thisshift)
        if len(momlist) == 0: momlist = self.NNQ_Stats.index.get_level_values('momentum')
        ip = momlist[0]
        iflow = thistflowlist[0]
        # for ip in momlist:
        #     for iflow in thistflowlist:
        ## cannot figure out a good solution for this test. It will just throw exception
        if 'All' in ip: ip = 'p000'
        if 'AlphaAuto' not in self.NNQ_Stats.columns:
            print('Autocorrlation result not present (might be combined result?), skipping tauintWopt plot')
            return plot_class
        if (ip,iflow) not in self.NNQ_Stats['AlphaAuto'].swaplevel(1,2):
            print(ip, iflow,' not found in AlphaAuto keys')
            return plot_class
        idata = self.NNQ_Stats['AlphaAuto']
        dataavg = idata.apply(lambda x : x.tau)
        dataerr = idata.apply(lambda x : x.tauerr)
        # tlist,values = idata.keys()[xlims[0]:xlims[1]],idata.values[xlims[0]:xlims[1]]
        # dataavg,dataerr = zip(*[(ival.,ival.) for ival in values])
        hold_series = pa.Series()
        hold_series['x_data'] = 'from_keys'
        # hold_series['x_data'] = np.array(map(untstr,tlist))
        hold_series['y_data'] = dataavg
        hold_series['yerr_data'] = dataerr
        hold_series['key_select'] = (ip,slice(None),iflow)
        hold_series['xerr_data'] = None
        hold_series['type'] = 'error_bar_vary'
        hold_series['fit_class'] = None
        hold_series['label'] = self.LegLab+' Auto'
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        hold_series['xdatarange'] = None
        hold_series['ShowPar'] = r'\alpha'
        hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
        plot_class.AppendData(hold_series)
        return plot_class
        # plot_class.get_xlim(*xlims)





    def AlphaPlotTflowTauIntWopt(   self,plot_class,tflowplot = defxlimOp,momlist=['p000'],
                                    thiscol='PreDefine',thissym='PreDefine',
                                    thisshift='PreDefine',thistsinklist = 'PreDef'):
        # def WithinBounds(tlist):
        #     output,indicies = [],[]
        #     for ic,ilist in enumerate(tlist):
        #         if untflowstr(ilist) in tflowplot:
        #             output.append(ilist)
        #             indicies.append(ic)
        #     return output,indicies
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        tflowplot = list(map(untflowstr,tflowplot))
        if thistsinklist == 'PreDef': thistsinklist = self.tsinkfit
        if len(momlist) == 0: momlist = self.NNQ_Stats.index.get_level_values('momentum')
        ip = momlist[0]
        itsink = thistsinklist[0]
        # for ip in momlist:
        #     for itsink in thistsinklist:
        if 'All' in ip: ip = 'p000'
        if 'AlphaAuto' not in self.NNQ_Stats.columns:
            print('Autocorrlation result not present (might be combined result?), skipping tauintwop plot')
            return plot_class
        if (ip,itsink) not in self.NNQ_Stats['AlphaAuto']:
            print(ip, itsink ,' not found in AlphaAuto keys')
            return plot_class
        tflowdata = self.NNQ_Stats['AlphaAuto']
        tflowlist = list(tflowdata.loc[(ip,itsink,slice(None))].keys())
        # tflowlist,indicies = WithinBounds(np.array(tflowlist))
        if len(tflowlist) > 0 and isinstance(tflowlist[0],(list,tuple,np.ndarray)):
            tflowlist = [itflow[2] for itflow in tflowlist]
        tflowlist = TflowToPhys(np.array(tflowlist),self.NN.latparams.latspace)
        plotshift = self.GetShift([tflowlist[0],tflowlist[-1]],thisshift)
        # print
        # print thisshift
        # print plotshift
        # print [tflowlist[0],tflowlist[-1]]
        # print np.array(tflowlist)+plotshift
        dataavg = tflowdata.apply(lambda x : x.tau)
        dataerr = tflowdata.apply(lambda x : x.tauerr)
        hold_series = pa.Series()
        hold_series['x_data'] = np.array(tflowlist)
        hold_series['y_data'] = dataavg
        hold_series['key_select'] = (ip,itsink,slice(None))
        hold_series['yerr_data'] = dataerr
        hold_series['xerr_data'] = None
        hold_series['type'] = 'error_bar_vary'
        hold_series['fit_class'] = None
        hold_series['label'] = self.LegLab+' Auto'
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = plotshift
        hold_series['xdatarange'] = None
        hold_series['ShowPar'] = r'\alpha'
        hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
        plot_class.AppendData(hold_series)
        return plot_class
        # plot_class.get_xlim(*xlims)


    ## make sure all xlims are the same, because this sets the window of xlim (needed for calculating shifts)
    ## thistflowlist is array of tflows to plot over (currently supports length 1 arrays).
    ## ToDo: make color and symbol and shift iterate over this. Currently handeled externally (SetsOfCorrs.py) with repeated calls
    def AlphaPlot(  self,plot_class,xlims=defxlimAlpha,momlist=['p000'],
                    thiscol='PreDefine',thissym='PreDefine',
                    thisshift='PreDefine',thistflowlist = 'PreDef',Auto=False):
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        self.thisshiftscale = 'Not Set'
        if thistflowlist == 'PreDef': thistflowlist = self.tflowfit
        if isinstance(xlims[0],str): xlims = list(map(untstr,xlims))
        thisshift = self.GetShift(xlims,thisshift)
        if Auto:
            if 'AlphaAuto' not in self.NNQ_Stats.columns:
                print('Autocorrlation result not present (might be combined result?), skipping alphaplot')
                return plot_class
            thisdata = self.NNQ_Stats['AlphaAuto']
        else:
            thisdata = self.NNQ_Stats['Alphaboot']
        if len(momlist) == 0: momlist = self.NNQ_Stats.index.get_level_values('momentum')
        ip = momlist[0]
        iflow = thistflowlist[0]
        if 'All' in ip: ip = 'p000'
        if (ip,iflow) not in thisdata.swaplevel(1,2):
            print(ip,iflow, ' not found in Alpha keys')
            return plot_class
        thisleglab = self.LegLab
        if Auto: thisleglab  = thisleglab + ' Auto'

        # xlims[0] = max(xlims[0],self.boot_tsum_list[0])
        # xlims[1] = min(xlims[1],self.boot_tsum_list[-1])

        dataavg = thisdata.apply(lambda x : x.Avg)
        dataerr = thisdata.apply(lambda x : x.Std)
        hold_series = pa.Series()
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = dataavg
        hold_series['yerr_data'] = dataerr
        hold_series['key_select'] = (ip,slice(None),iflow)
        hold_series['xerr_data'] = None
        hold_series['type'] = 'error_bar_vary'
        hold_series['fit_class'] = None
        hold_series['label'] = thisleglab
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        hold_series['xdatarange'] = None
        hold_series['ShowPar'] = r'\alpha'
        hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
        plot_class.AppendData(hold_series)
        return plot_class
        # plot_class.get_xlim(*xlims)



    def AlphaPlotTflow( self,plot_class,tflowplot = defxlimOp,momlist=['p000'],
                        thiscol='PreDefine',thissym='PreDefine',
                        thisshift='PreDefine',thistsinklist = 'PreDef',Auto=False):
        # def WithinBounds(tlist):
        #     output,indicies = [],[]
        #     for ic,ilist in enumerate(tlist):
        #         if untflowstr(ilist) in tflowplot:
        #             output.append(ilist)
        #             indicies.append(ic)
        #     return output,indicies
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        tflowplot = list(map(untflowstr,tflowplot))
        if thistsinklist == 'PreDef': thistsinklist = self.tsinkfit
        if Auto:
            if 'AlphaAuto' not in self.NNQ_Stats.columns:
                print('Autocorrlation result not present (might be combined result?), skipping tflowplot')
                return plot_class
            thisdata = self.NNQ_Stats['AlphaAuto']
        else:
            thisdata = self.NNQ_Stats['Alphaboot']
        if len(momlist) == 0: momlist = self.NNQ_Stats.index.get_level_values('momentum')
        ip = momlist[0]
        itsink = thistsinklist[0]
        if 'All' in ip: ip = 'p000'
        if (ip,itsink) not in thisdata:
            print(ip,itsink, ' not found in Alpha keys')
            return plot_class
        thisleglab = self.LegLab
        if Auto: thisleglab = thisleglab + ' Auto'
        tflowdata = thisdata[ip,itsink]
        tflowlist = list(tflowdata.keys())
        # tflowlist,indicies = WithinBounds(np.array(tflowlist))
        if len(tflowlist) > 0 and isinstance(tflowlist[0],(list,tuple,np.ndarray)):
            tflowlist = [itflow[2] for itflow in tflowlist]
        tflowlist = TflowToPhys(np.array(tflowlist),self.NN.latparams.latspace)
        plotshift = self.GetShift([tflowlist[0],tflowlist[-1]],thisshift)
        # print
        # print thisshift
        # print plotshift
        # print [tflowlist[0],tflowlist[-1]]
        # print np.array(tflowlist)+plotshift
        dataavg = thisdata.apply(lambda x : x.Avg)
        dataerr = thisdata.apply(lambda x : x.Std)
        hold_series = pa.Series()
        hold_series['x_data'] = np.array(tflowlist)
        hold_series['y_data'] = dataavg
        hold_series['yerr_data'] = dataerr
        hold_series['xerr_data'] = None
        hold_series['type'] = 'error_bar_vary'
        hold_series['key_select'] = (ip,itsink,slice(None))
        hold_series['fit_class'] = None
        hold_series['label'] = thisleglab
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = plotshift
        hold_series['xdatarange'] = None
        hold_series['ShowPar'] = r'\alpha'
        hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
        plot_class.AppendData(hold_series)
        return plot_class
        # plot_class.get_xlim(*xlims)



    #state is 1 or 2 state fit
    ## fitr is fit range e.g. fitr5-10
    ## fit range of 'Data' is just the original fit range
    def AlphaFitPlot(   self,plot_class,fitr,xlims = 'Data',momlist=['p000'],
                        thiscol='PreDefine',thisshift='PreDefine',
                        thistflowlist = 'PreDef',Auto=False,WipeFit=False):
        # autocorrelation has not fits (currenlty)
        if Auto: return plot_class
        if fitr == 'None': return plot_class
        self.CheckCol(thiscol)
        if isinstance(fitr,(tuple,list,np.ndarray)):
            fitr = '_'.join(fitr)
        thisshift = self.GetShift(plot_class.get_xlim(),thisshift)
        if len(momlist) == 0: momlist = self.NNQ_Stats.index.get_level_values('momentum')
        if thistflowlist == 'PreDef': thistflowlist = self.tflowfit
        imom = momlist[0]
        iflow = thistflowlist[0]
        if 'All' in imom: imom = 'p000'
        if WipeFit or 'boot' not in self.NNQ_Fit_Stats or (imom,iflow) not in self.NNQ_Fit_Stats['boot']:
            print(imom,iflow ,' not found in FitDict, fitting now')
            self.FitAlpha(fitr,thistflowlist = [iflow],WipeFit=WipeFit)

        fit_data = sff.PullSOFSeries(self.NNQ_Fit_Stats['boot'],fmted=True)
        this_fun = fit_data.index[0][2]
        this_key = (imom,iflow,this_fun,fitr)
        if this_key not in fit_data.index:
            print(this_key, ' not found in FitDict, performing fit now')
            print(fit_data)
            self.FitAlpha(fitr,thistflowlist = [iflow],WipeFit=WipeFit)
        hold_series = null_series
        hold_series['type'] = 'fit_vary'
        hold_series['fit_class'] = fit_data
        hold_series['key_select'] = this_key
        hold_series['label'] = 'Fit '
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        hold_series['xdatarange'] = xlims
        hold_series['ShowPar'] = r'\alpha'
        hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
        plot_class.AppendData(hold_series)
        return plot_class

    #state is 1 or 2 state fit
    ## fitr is fit range e.g. fitr5-10
    ## fit range of 'Data' is just the original fit range
    def AlphaNNQFitPlot(self,plot_class,fitr,xlims = 'Data',momlist=['p000'],
                        thiscol='PreDefine',thisshift='PreDefine',
                        thistflowlist = 'PreDef',Auto=False):
        print('NNQFitPlot not implemented with setoffitfuns, skipping')
        # if Auto: return plot_class
        # if fitr == 'None': return plot_class
        # self.CheckCol(thiscol)
        # thisshift = self.GetShift(plot_class.get_xlim(),thisshift)
        # if len(momlist) == 0: momlist = self.NNQ_Stats.index.get_level_values('momentum')
        # if thistflowlist == 'PreDef': thistflowlist = self.tflowfit
        # imom = momlist[0]
        # iflow = thistflowlist[0]
        # if 'All' in imom: imom = 'p000'
        # if 'NNQboot' not in self.NNQ_Fit_Stats or (imom,iflow,fitr) not in self.NNQ_Fit_Stats['NNQboot']:
        #     print imom,iflow,fitr, ' not found in NNQFitDict, fitting now'
        #     self.FitNNQ([fitr],thistflowlist = [iflow])
        # # print 'debug',self.NNQFitDict[imom][iflow][fitr]
        # # print self.HumanFile
        # # print 'results'
        # # for iparam,paramval in self.NNQFitDict[imom][iflow][fitr].Params.iteritems():
        # #     print imom,iflow,fitr,iparam,paramval.Avg
        # # print
        # hold_series = pa.Series()
        # hold_series['x_data'] = None
        # hold_series['y_data'] = None
        # hold_series['yerr_data'] = None
        # hold_series['xerr_data'] = None
        # hold_series['type'] = 'fit_vary'
        # hold_series['key_select'] = (imom,iflow,fitr)
        # hold_series['fit_class'] = self.NNQ_Fit_Stats['NNQboot']
        # hold_series['label'] = 'NNQ Fit'
        # hold_series['symbol'] = self.thissym
        # hold_series['color'] = self.thiscol
        # hold_series['shift'] = thisshift
        # hold_series['xdatarange'] = xlims
        # # hold_series['ShowPar'] = r'\alpha'
        # hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
        # plot_class.AppendData(hold_series)
        # return plot_class

    #state is 1 or 2 state fit
    ## fitr is fit range e.g. fitr5-10
    ## fit range of 'Data' is just the original fit range
    def AlphaPlotVsEnergy(  self,plot_class,fitr,momlist=['p000'],
                            thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',
                            thistflowlist = 'PreDef',Auto=False,PhysEnergy=True,WipeFit=False):
        # autocorrelation has not fits (currenlty)
        if Auto: return plot_class
        if fitr == 'None': return plot_class
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        thisshift = self.GetShift(plot_class.get_xlim(),thisshift)
        if len(momlist) == 0: momlist = self.NNQ_Stats.index.get_level_values('momentum')
        if thistflowlist == 'PreDef': thistflowlist = self.tflowfit

        iflow = thistflowlist[0]
        fit_data = self.Get_Extrapolation()
        fit_data = self.NN.latparams.Fmt_pmom_index(fit_data,0)
        this_fun = fit_data.index[0][2]
        fit_param = fit_data.index[0][4]
        this_key = (slice(None),iflow,this_fun,fitr,fit_param)
        test_key = (fit_param,iflow,this_fun,fitr)
        if WipeFit or 'boot' not in self.NNQ_Fit_Stats or test_key not in fit_data.swaplevel(0,4):
            print(test_key, ' is not in FitDictionary, fitting now')
            self.FitAlpha(fitr,thistflowlist = [iflow],WipeFit=WipeFit)

        # this_x = self.NN.latparams.EnergyFromRecRel(momlist,self.NN.latparams.GetNucleonMass(Phys=False))
        # if PhysEnergy: this_x = this_x*self.NN.latparams.hbarcdivlat
        # plotx.append(this_x)

        ## into GeV
        # plotshift = self.GetShift([plotx[0],plotx[-1]],thisshift)
        hold_series = null_series
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = fit_data.apply(lambda x : x.Avg)
        hold_series['yerr_data'] = fit_data.apply(lambda x : x.Std)
        hold_series['key_select'] = this_key
        hold_series['type'] = 'error_bar_vary'
        hold_series['label'] = self.LegLab
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
        plot_class.AppendData(hold_series)
        return plot_class



    # ## make sure all xlims are the same, because this sets the window of xlim (needed for calculating shifts)
    # ## thistflowlist is array of tflows to plot over (currently supports length 1 arrays).
    # ## ToDo: make color and symbol and shift iterate over this. Currently handeled externally (SetsOfCorrs.py) with repeated calls
    # def AlphaPlotAuto(self,xlims=defxlimAlpha,momlist=['p000'],thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',thistflowlist = 'PreDef'):
    #     self.CheckCol(thiscol)
    #     self.CheckSym(thissym)
    #     self.thisshiftscale = 'Not Set'
    #     if thistflowlist == 'PreDef': thistflowlist = self.tflowfit
    #     if len(momlist) == 0: momlist = self.AlphaAuto.keys()
    #     if isinstance(xlims[0],basestring): xlims = map(untstr,xlims)
    #     thisshift = self.GetShift(xlims,thisshift)
    #     for ip in momlist:
    #         if ip not in self.AlphaAuto.keys():
    #             print ip, ' not found in AlphaAuto keys'
    #             continue
    #         for iflow in thistflowlist:
    #             thisleglab = self.LegLab +'$\ '+ tflowTOLeg(iflow) + '$'
    #             if len(momlist) > 1: thisleglab = thisleglab+'\ '+ip.replace('p','p=')
    #             pdata = self.AlphaAuto[ip]
    #             tlist,dataAuto = zip(*pdata.items())
    #             if any([iflow not in iAuto.keys() for iAuto in dataAuto]):
    #                 print iflow, ' not found in AlphaAuto'
    #                 continue
    #             dataAuto = [idata[iflow] for idata in dataAuto[xlims[0]:xlims[1]]]
    #             dataavg,dataerr = zip(*[(idata.Avg,idata.Std) for idata in dataAuto])
    #             tlist = tlist[xlims[0]:xlims[1]]
    #             pl.errorbar(np.array(map(untstr,tlist))+thisshift,dataavg,dataerr,label=thisleglab,fmt=self.thissym,color=self.thiscol)
    #     # plot_class.get_xlim(*xlims)



    # def AlphaPlotTflowAuto(self,tflowplot = defxlimOp,momlist=['p000'],thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',thistsinklist = 'PreDef'):
    #     def WithinBounds(tlist):
    #         output,indicies = [],[]
    #         for ic,ilist in enumerate(tlist):
    #             if untflowstr(ilist) in tflowplot:
    #                 output.append(ilist)
    #                 indicies.append(ic)
    #         return output,indicies
    #     self.CheckCol(thiscol)
    #     self.CheckSym(thissym)
    #     tflowplot = map(untflowstr,tflowplot)
    #     if thistsinklist == 'PreDef': thistsinklist = self.tsinkfit
    #     if len(momlist) == 0: momlist = self.AlphaAuto.keys()
    #     for ip in momlist:
    #         if ip not in self.AlphaAuto.keys():
    #             print ip, ' not found in AlphaAuto keys'
    #             continue
    #         for itsink in thistsinklist:
    #             if itsink not in self.AlphaAuto[ip].keys():
    #                 print itsink, ' not found in AlphaAuto'
    #                 continue
    #             thisleglab = self.LegLab +'$\ '+ tsinkTOLeg(itsink) + '$'
    #             if len(momlist) > 1: thisleglab = thisleglab+'\ '+ip.replace('p','p=')
    #             tflowdata = self.AlphaAuto[ip][itsink]
    #             tflowlist,dataAuto = zip(*tflowdata.items())
    #             tflowlist,indicies = WithinBounds(np.array(tflowlist))
    #             tflowlist = TflowToPhys(np.array(tflowlist),self.NN.latparams.latspace)
    #             plotshift = self.GetShift([tflowlist[0],tflowlist[-1]],thisshift)
    #             # print
    #             # print thisshift
    #             # print plotshift
    #             # print [tflowlist[0],tflowlist[-1]]
    #             # print np.array(tflowlist)+plotshift
    #             dataAuto = np.array(dataAuto)[indicies]
    #             dataavg,dataerr = zip(*[(idata.Avg,idata.Std) for idata in dataAuto])
    #             pl.errorbar(np.array(tflowlist)+plotshift,dataavg,dataerr,label=thisleglab,fmt=self.thissym,color=self.thiscol)
    #     # plot_class.get_xlim(*xlims)



    # TODO: implement FitEffMass function with same layout as TwoPtCorr

    # ## make sure all xlims are the same, because this sets the window of xlim (needed for calculating shifts)
    # def EffMassPlot(self,xlims = defxlim,momlist=['p000'],thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine'):
    #     self.CheckCol(thiscol)
    #     self.CheckSym(thissym)
    #     thisshift = self.GetShift(xlims,thisshift)
    #     if len(momlist) == 0:
    #         momlist = self.EffMAvg.keys()
    #     for ip in momlist:
    #         if ip not in self.EffMAvg.keys():
    #             print ip, ' not found in effective mass list'
    #             continue
    #         pdata = self.EffMAvg[ip]
    #         tlist,dataavg = zip(*pdata.items())
    #         tlist,dataerr = zip(*self.EffMStd[ip].items())
    #         tlist = tlist[xlims[0]:xlims[1]]
    #         dataavg = dataavg[xlims[0]:xlims[1]]
    #         dataerr = dataerr[xlims[0]:xlims[1]]
    #         pl.errorbar(np.array(map(untstr,tlist))+thisshift,dataavg,dataerr,label=self.LegLab,fmt=self.thissym,color=self.thiscol)
    #     plot_class.get_xlim(*xlims)


    ## state is 1 or 2 state fit
    # ## fitr is fit range e.g. fitr5-10
    # ## fit range of 'Data' is just the original fit range
    # def EffMassFitPlot(self,state,fitr,xlims = 'Data',momlist=['p000'],thiscol='PreDefine',thisshift='PreDefine'):
    #     self.CheckCol(thiscol)
    #     thisshift = self.GetShift(plot_class.get_xlim(),thisshift)
    #     if len(momlist) == 0:
    #         momlist = self.EffMAvg.keys()
    #     if isinstance(state,int): state = 'state'+str(state)
    #     if fitr not in self.FitDict[state].keys():
    #         print fitr , ' has not been done, performing fit now'
    #         self.Fit(int(state.replace('state','')),[fitr])
    #     fitr_params = self.FitDict[state][fitr]
    #     for imom in momlist:
    #         if imom not in fitr_params.keys():
    #             print imom, ' not found in FitDict'
    #             continue
    #         fitr_params[imom].PlotEffMass(xdatarange=xlims,color=self.thiscol,shift=thisshift)



    ## Comparisons

    def __str__(self):
        return self.name


    def __iter__(self):
        return self

    def items(self):
        return list(self.NNQ_Stats['boot'].items())


    def itemsAvgStd(self):
        outlist = []
        for ivals in self.NNQ_Stats.itertuples():
            outlist.append((ivals[0],ivals[2],ivals[3]))
        return outlist


    def values(self):
        return self.NNQ_Stats['boot'].values

    def valuesAvgStd(self):
        outlist = []
        for ivals in self.NNQ_Stats.itertuples(index=False):
            outlist.append((ivals[1],ivals[2]))
        return outlist


    def keys(self):
        return self.NNQ_Stats.index


    def __setitem__(self, key, value ):
        self.NNQ_Stats.loc[key,'boot'] = value

    def __getitem__(self, key):
        return self.NNQ_Stats.loc[key,'boot']

    def __reversed__(self):
        self.NNQ_Stats = self.NNQ_Stats.iiloc[::-1]
        return self

    def __next__(self):
        if self.current >= len(list(self.keys())):
            self.current = 0
            raise StopIteration
        else:
            self.current += 1
            if self.current >= len(list(self.keys()))+1:
                return 0.0
            else:
                thiskey = list(self.keys())[self.current-1]
                return self[thiskey]



    def EffMass(self):
        if 'boot' not in self.NNQ_Stats: return
        lEff,lEffA,lEffS = [],[],[]
        ilist = []
        for (ip,itflow),tdata in self.NNQ_Stats['boot'].groupby(level=('momentum','flow_time')):
            next_tdata = np.roll(tdata.swaplevel(1,2)[ip,itflow],1)
            for (it,idata),idata_shift in zip(iter(tdata.swaplevel(1,2)[ip,itflow].items()),next_tdata):
                thisEff = (idata/idata_shift).Log()
                thisEff.Stats()
                ilist.append((ip,it,itflow))
                lEff.append(thisEff)
                lEffA.append(thisEff.Avg)
                lEffS.append(thisEff.Std)
        indicies = pa.MultiIndex.from_tuples(ilist,names=self.NNQ_col_names)
        self.NNQ_Stats.loc[:,'EffM'] = pa.Series(lEff,index=indicies)
        self.NNQ_Stats.loc[:,'EffMAvg'] = pa.Series(lEffA,index=indicies)
        self.NNQ_Stats.loc[:,'EffMStd'] = pa.Series(lEffS,index=indicies)

    def RatFunAndDer(self):
        def chit(*vals):
            return vals[0]/vals[1]
        def chitder(*vals):
            return [1.0/vals[1],-vals[0]/(vals[1]**2.)]
        return (chit,chitder)

    def AlphaRatio(self,DefWipeG2=False,DoAuto=True,show_timer= True):
        if self.IsSubComb: return
        cpevenInfo = deepcopy(self.Info)
        cpevenInfo['Interp'] = 'CPEven'


        # print 'Debug iscomb , ' , self.IsComb
        if self.IsComb:
            ## TODO, autocorrelation needs to be done for combining correlators...
            ## NOT FEASABLE, too much memory me thinks...
            from SetsOfTwoPtCorrs import SetOfTwoPt
            sinklab,jsmcoeffs,sm2ptList,InfoList = GetInfoList(cpevenInfo)
            C2Set = SetOfTwoPt(cfglist=self.NNQ_cfgs[['configs','xsrc_list']],InfoDict = InfoList)
            C2Set.LoadPickle(DefWipe=DefWipeG2,CheckCfgs=True,CheckMom=True)
            # for iset in C2Set.SetC2.values():
            #     if len(iset.C2) == 0 and DoAuto:
            #         iset.LoadPickle(DefWipe=True,CheckCfgs=True,CheckMom=True,WipeData=not DoAuto)
            thisname,thisleglab,jsmlab = CreateNameAndLL(C2Set.SetC2,sinklab)
            thisfun = GenSinkFun(jsmcoeffs)
            thisfile = self.filename.replace(self.jsm,sinklab)
            thisLL = self.LegLab.replace(self.jsm,sinklab)
            self.SetCustomName(thisfile,thisLL)
            self.jsm = jsmlab
            self.NNCPEven = C2Set.CombineCorrs(thisfun,filename=thisfile,LegLab=thisLL,jsmlab=jsmlab)
            thisname = self.NNCPEven.name + '_ForAlpha'
            self.NNCPEven.SetCustomName(string=thisname)
            self.NNCPEven.Stats()
        else:
            self.NNCPEven = TwoPointCorr(   self.nboot,self.NNQ_cfgs[['configs','xsrc_list']],cpevenInfo,
                                            self.thissym,self.thiscol,self.thisshift)
            # self.NNCPEven = self.NNCPEven.GetCombCorr(cpevenInfo)
            thisname = self.NNCPEven.name + '_ForAlpha'
            self.NNCPEven.SetCustomName(string=thisname)
            self.NNCPEven.LoadPickle(DefWipe=DefWipeG2,CheckCfgs=True,CheckMom=True,WipeData=not DoAuto)
            if 'C2' not in self.NNCPEven.C2_cfgs or not DoAuto:
                self.NNCPEven.LoadPickle(DefWipe=True,CheckCfgs=True,CheckMom=True,WipeData=not DoAuto)



        self.CPEvenFile = self.NNCPEven.HumanFile
        if DoAuto and 'NNQ' in self.NNQ_cfgs and not self.IsComb:
            thisfuns = self.RatFunAndDer()
            lAl,lAlAvg,lAlStd = [],[],[]
            ilist = []
            this_gen = itertools.product(enumerate(self.NN.pform),list(range(self.nt)),enumerate(self.FO.tflowlist))
            if show_timer: thistimer = Timer(linklist=list(this_gen),name='Autcor '+self.name)
            this_gen = itertools.product(enumerate(self.NN.pform),list(range(self.nt)),enumerate(self.FO.tflowlist))
            for ((icp,ip),ict,(ictflow,itflow)) in this_gen:
                it = tstr(ict+1)
                # cant do these checks XO.... dictionaries are better!
                # if ip not in self.NNCPEven.C2.keys(): continue
                # if it not in self.NNCPEven.C2[ip].keys(): continue
                # if ('C2' not in self.NNCPEven.C2_cfgs
                #     or icp*self.nt+ict >= self.NNCPEven.C2_cfgs['C2'].size): continue
                auto_df = self.NNQ_cfgs['NNQ'].apply(lambda x : x[icp][ict][ictflow]).to_frame(name='NNQ')
                auto_df.loc[:,'CPEven'] = pa.Series(self.NNCPEven.C2_cfgs['C2'].apply(lambda x : x[icp][ict]).values,index=auto_df.index)
                # print 'debug',idata[5],twoptdata[5]
                # if len(idata) != len(twoptdata):
                #     errorstr = '''
                #     NNQ and NN have different configuration lists, \n
                #     maybe force delete the pickled files for NN, \n
                #     or select a different -stats from- directory corresponding to operator Q or W
                #     '''
                #     raise EnvironmentError(errorstr)
                ilist.append((ip,it,itflow))

                thisAuto =  AutoCorrelate(  Fun=thisfuns,Sparam=self.Sparam,
                                            name=self.name + ip+' $t='+str(it)+'\ '+tflowTOLeg(itflow)+'$',data=auto_df)
                lAl.append(thisAuto)
                lAlAvg.append(thisAuto.Avg)
                lAlStd.append(thisAuto.Std)
                if show_timer: thistimer.Lap((ip,it,itflow))
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.NNQ_col_names)
            self.NNQ_Stats.loc[:,'AlphaAuto'] = pa.Series(lAl,index=indicies)
            self.NNQ_Stats.loc[:,'AlphaAutoAvg'] = pa.Series(lAlAvg,index=indicies)
            self.NNQ_Stats.loc[:,'AlphaAutoStd'] = pa.Series(lAlStd,index=indicies)
            del self.NNQ_cfgs['NNQ']

        lAl,lAlAvg,lAlStd = [],[],[]
        ilist = []
        for (ip,it,itflow),idata in list(self.items()):
            if (ip,it) not in self.NNCPEven.C2_Stats['boot']: continue
            ilist.append((ip,it,itflow))
            idata.Stats()
            thisboot = idata/self.NNCPEven.C2_Stats['boot'][ip,it]
            thisboot.Stats()
            # print 'DEBUG'
            # print 'iscomb',self.IsComb
            # print 'nncpeven',self.NNCPEven
            # print 'nncpodd',self.NN
            # print 'flow',self.FO
            # print (ip,it,itflow, idata.Avg,idata.Std,
            #                     self.NNCPEven.C2_Stats['boot'][ip,it].Avg,self.NNCPEven.C2_Stats['boot'][ip,it].Std,
            #                     thisboot.Avg, thisboot.Std)
            lAl.append(thisboot)
            lAlAvg.append(thisboot.Avg)
            lAlStd.append(thisboot.Std)
        indicies = pa.MultiIndex.from_tuples(ilist,names=self.NNQ_col_names)
        self.NNQ_Stats.loc[:,'Alphaboot'] = pa.Series(lAl,index=indicies)
        self.NNQ_Stats.loc[:,'AlphaAvg'] = pa.Series(lAlAvg,index=indicies)
        self.NNQ_Stats.loc[:,'AlphaStd'] = pa.Series(lAlStd,index=indicies)


    # def ReadAndBootAndEffM(self,DefWipe=False):
    #     if os.path.isfile(self.PickleFile) and not DefWipe:
    #         self.LoadPickle()
    #     else:
    #         self.ReadAndBoot()
    #         self.EffMass()


    def Write(self,WipeData=True):
        def FixDictArray(thislist,ikey):
            outDict = OrderedDict()
            for ic,ilist in enumerate(thislist):
                outDict[ikey+str(ic+1)] = ilist
            return outDict

        # self.FO.FlowWrite()
        # if not self.IsComb:
        #     self.NN.Write()

        if not hasattr(self,'show_cfgs'):
            self.show_cfgs = False

        outDict = ODNested()
        if self.IsComb:
            outDict['TwoPtCorr_File'] = 'Not outputted'
        #     self.NN.HumanFile = self.NN.HumanFile.replace(self.NN.jsm,self.jsm)
        else:
            outDict['TwoPtCorr_File'] = self.NN.HumanFile
        outDict['FlowOp_File'] = self.FO.flowHumanFile
        if self.CPEvenFile != 'None':
            outDict['TwoPtCorr_CPEven_file'] = self.CPEvenFile
        excel_params = pa.Series(deepcopy(outDict))
        outDict['name'] = self.name
        outDict['kud'] = self.NN.kud
        outDict['ks'] = self.NN.ks
        outDict['nt'] = self.NN.nt
        outDict['nxyz'] = self.NN.nxyz

        for istream,incfg in zip(self.stream_list,self.ncfg_list):
            outDict['ncfg_'+istream] = incfg
        outDict['nstream'] = self.NN.nstream
        outDict['nxsrc_avg'] = self.NN.xsrcAvg
        outDict['nMeas'] = self.NN.nMeas
        outDict['nboot'] = self.NN.nboot

        outDict['Fits']['Function'] = self.Fun.__name__
        if self.AlphaiGuess == 'None':
            outDict['Fits']['Initial_Guess'] = 'Default (see FitFunctions.py)'
        else:
            for iparam in range(1,self.npar+1):
                thispar = 'Par'+str(iparam)
                outDict['Fits']['Initial_Guess'][thispar] = self.AlphaiGuess[iparam]
        if 'boot' in self.NNQ_Fit_Stats:
            for (ip,itflow),fitdict in self.NNQ_Fit_Stats['boot'].items():
                if not isinstance(fitdict,sff.SetOfFitFuns): continue
                outDict['Fits'][ip][itflow] = fitdict.GetOutputDict()


        outDict['NNQFits']['Function'] = self.Fun.__name__
        if self.NNQiGuess == 'None':
            outDict['NNQFits']['Initial_Guess'] = 'Default (see FitFunctions.py)'
        else:
            for iparam in range(1,self.npar+1):
                thispar = 'Par'+str(iparam)
                outDict['NNQFits']['Initial_Guess'][thispar] = self.NNQiGuess[iparam]
        # if 'NNQboot' in self.NNQ_Fit_Stats:
        #     for (ip,itflow,ifit),fitdict in self.NNQ_Fit_Stats['NNQboot'].iteritems():
        #         if not isinstance(fitdict,ff.Fitting): continue
        #         outDict['Fits'][ip][itflow][ifit] = fitdict.GetOutputDict()


        for col_key,col_data in self.NNQ_Stats.items():
            if 'Avg' in col_key or 'Std' in col_key: continue
            for (ip,it,itflow),idata in col_data.items():
                outDict[col_key][ip][it][itflow] = AvgStdToFormat(idata.Avg,idata.Std)

        if any(np.array(self.ncfg_list) > 0) and self.show_cfgs:
            for ((istream,iccfg),icfg),ixsrc_list in zip(iter(self.NNQ_cfgs['configs'].items()),self.NNQ_cfgs['xsrc_list'].values):
                outDict['cfglist'][istream][icfg] = 'n_xsrc_'+str(len(ixsrc_list))

        if WipeData: self.RemoveVals()

        ## Write human readable file with data
        WriteXml(self.HumanFile,{'Results':outDict})
        ## pickles rest of data for reading

        if any(np.array(self.ncfg_list) > 0):
            for ((istream,iccfg),icfg),ixsrc_list in zip(iter(self.NNQ_cfgs['configs'].items()),self.NNQ_cfgs['xsrc_list'].values):
                outDict['cfglist'][istream][icfg] = 'n_xsrc_'+str(len(ixsrc_list))
        ## TODO, include fitting in the excell file
        WriteExcel(self.ExcelFile,{'NNQ_results':deepcopy(self.NNQ_Stats)},cfglist=self.NNQ_cfgs,params=excel_params)

        self.RemoveFuns()
        WritePickle(self.PickleFile,self.__dict__)
        self.GetFuns()
        ##################################

    def ReadAndWrite(self,WipeData=True,NoFit=False):
        self.Read()
        self.Bootstrap(WipeData=WipeData)
        self.Stats()
        if NoFit:
            self.FitAlpha(WipeFit=WipeData)
        # self.FitNNQ()
        self.Write()


    def FixRead(self,readdict):
        readdict['Info']['outdir'] = self.Info['outdir']
        readdict['FO'].latparams.Info = self.FO.latparams.Info
        readdict['FO'].latparams.outputdir = self.FO.latparams.outputdir
        readdict['FO'].latparams.PickleFile = self.FO.latparams.PickleFile
        readdict['FO'].latparams.HumanFile = self.FO.latparams.HumanFile
        readdict['NN'].latparams.outputdir = self.NN.latparams.outputdir
        readdict['NN'].latparams.PickleFile = self.NN.latparams.PickleFile
        readdict['NN'].latparams.HumanFile = self.NN.latparams.HumanFile
        readdict['FO'].flowHumanFile = self.FO.flowHumanFile
        readdict['FO'].flowPickleFile = self.FO.flowPickleFile
        readdict['NN'].HumanFile = self.NN.HumanFile
        readdict['NN'].PickleFile = self.NN.PickleFile

        readdict['PickleFile'] = self.PickleFile
        readdict['HumanFile'] = self.HumanFile
        readdict['NNQ_Fit_col_names'] = self.NNQ_Fit_col_names

        readdict['IsSubComb'] = self.IsSubComb
        if not os.path.isfile(readdict['NNfile']):
            readdict['NNfile'] = self.NNfile
        if not os.path.isfile(readdict['FOfile']):
            readdict['FOfile'] = self.FOfile
        # if not os.path.isfile(readdict['ppparams'].PickleFile):
        #     readdict['ppparams'].PickleFile = self.ppparams.PickleFile
        #     readdict['ppparams'].HumanFile = self.ppparams.HumanFile
        #     readdict['qparams'].PickleFile = self.qparams.PickleFile
        #     readdict['qparams'].HumanFile = self.qparams.HumanFile
        #     readdict['Info']['outdir'] = self.Info['outdir']
        return readdict




    ## TODO: Pickle load and dump self.NN and self.FO from their respective directories.
    ##       At this point, I am just dumping this whole class instead. (doubles up data)
    def LoadPickle(self,WipeData=True,DefWipe = False,CheckCfgs=False,CheckMom=True,legacy=False,NoFit=False):
        # print 'Loading Pickle',self.PickleFile
        read_pickle = self.PickleFile.replace('_'.join([self.NN.MesOrBar,self.NN.stream_name])+'_','')
        if os.path.isfile(self.PickleFile) and not DefWipe and not legacy:
            loadeddict = ReadPickleWrap(self.PickleFile)
            loadeddict = self.FixRead(loadeddict)
            checklist = self.thischecklist
            if CheckCfgs: checklist = checklist + ['NNQ_cfgs']
            if CheckMom: checklist = checklist + ['checkmomlist']
            loadeddict['checkmomlist'] = [self.NN.latparams.pformTOpstr(ip,nopref=True) for ip in loadeddict['NNQ_Stats'].index.get_level_values('momentum')]
            if CheckClass(self.__dict__,loadeddict,checklist):
                print('loading pickle:')
                print(self.PickleFile)
                print()
                self.__dict__.update(loadeddict)
                if Wipe_All_Fits:
                    self.NNQ_Fit_Stats = pa.DataFrame()
                elif len(self.NNQ_Fit_Stats.values)> 0 and not isinstance(self.NNQ_Fit_Stats.values[0][0],sff.SetOfFitFuns):
                    print('Legacy file, fits need to be recomputed')
                    self.NNQ_Fit_Stats = pa.DataFrame()
            else:
                if os.path.isfile(self.PickleFile+'.bak'):
                    print('using backupfile for '+self.name)
                    self.PickleFile = self.PickleFile+'.bak'
                    self.HumanFile = self.HumanFile+'.bak'
                    self.ExcelFile = self.ExcelFile.replace('.xlsx','.bak.xlsx')
                    self.LoadPickle(DefWipe=DefWipe,WipeData=WipeData,CheckCfgs=CheckCfgs,CheckMom=CheckMom)
                    return
                print('Warning, file ' , self.PickleFile , ' has different parameters than this instance, reading and writing over file now')
                print('Attempting to read legacy file')
                print()
                self.PickleFile = self.PickleFile.replace('.bak','')
                self.HumanFile = self.HumanFile.replace('.bak','')
                self.ExcelFile = self.ExcelFile.replace('.bak','')
                self.LoadPickle(DefWipe=DefWipe,WipeData=WipeData,CheckCfgs=CheckCfgs,CheckMom=CheckMom,legacy=True)
                return
                # self.ReadAndWrite(WipeData=WipeData)
            # self.Write()
        elif os.path.isfile(read_pickle) and not DefWipe:
            loadeddict = ReadPickleWrap(read_pickle)
            loadeddict = self.FixRead(loadeddict)
            checklist = self.thischecklist
            if CheckCfgs: checklist = checklist + ['NNQ_cfgs']
            if CheckMom: checklist = checklist + ['checkmomlist']
            loadeddict['checkmomlist'] = [self.NN.latparams.pformTOpstr(ip,nopref=True) for ip in loadeddict['NNQ_Stats'].index.get_level_values('momentum')]
            if CheckClass(self.__dict__,loadeddict,checklist):
                print('Legacy file found, reading and writing to proper file')
                print(read_pickle)
                print()
                self.__dict__.update(loadeddict)
                self.Write()
            else:
                if os.path.isfile(self.PickleFile+'.bak'):
                    print('using backupfile for '+self.name)
                    self.PickleFile = self.PickleFile+'.bak'
                    self.HumanFile = self.HumanFile+'.bak'
                    self.ExcelFile = self.ExcelFile.replace('.xlsx','.bak.xlsx')
                    self.LoadPickle(DefWipe=DefWipe,WipeData=WipeData,CheckCfgs=CheckCfgs,CheckMom=CheckMom,legacy=True)
                    return
                print('Warning, legacy file ' , read_pickle , ' has different parameters than this instance, reading and writing over file now')
                print()
                self.ReadAndWrite(WipeData=WipeData,NoFit=NoFit)
        else:
            if not DefWipe:
                print('Pickle file not found:')
                print(self.PickleFile)
                print('or legacy file')
                print(read_pickle)
                print()
            self.ReadAndWrite(WipeData=WipeData,NoFit=NoFit)
    ## Operator overloading


    ## unitary arithmetic operations

    def __neg__(self):
        for ival in list(self.keys()):
            self[ival] = -self[ival]
        return self


    def __pos__(self):
        return self

    def __abs__(self):
        for ival in list(self.keys()):
            self[ival] = abs(self[ival])
        return self

    ## Casting, just casts the bootstrap values, ignore return value (only needed since return has to be same type).


    def __complex__(self):
        for ival in list(self.keys()):
            complex(self[ival])
        self.Stats()
        return 1+0j

    def __int__(self):
        for ival in list(self.keys()):
            int(self[ival])
        self.Stats()
        return 1

    def __float__(self):
        for ival in list(self.keys()):
            np.float64(self[ival])
        self.Stats()
        return 1.0

    def sin(self):
        result = NNQCorr(thisnboot=self.nboot, cfglist=self.NNQ_cfgs[['configs','xsrc_list']],
                         Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.NNQ_Stats.loc[:,'boot'] = np.sin(self.NNQ_Stats.loc[:,'boot'])
        return result

    def cos(self):
        result = NNQCorr(thisnboot=self.nboot, cfglist=self.NNQ_cfgs[['configs','xsrc_list']],
                         Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.NNQ_Stats = pa.DataFrame(columns = ['boot'],index=self.NNQ_Stats.index)
        result.NNQ_Stats.loc[:,'boot'] = np.cos(self.NNQ_Stats.loc[:,'boot'])
        return result


    def Overload_wrap(self,NNQ2,this_fun):
        result = NNQCorr(thisnboot=self.nboot, cfglist=self.NNQ_cfgs[['configs','xsrc_list']],
                              Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.NNQ_Stats = pa.DataFrame(columns = ['boot'],index=self.NNQ_Stats.index)
        # result.UpdateName('+',self,NNQ2)
        if isinstance(NNQ2,NNQCorr):
            result.NNQ_Stats.loc[:,'boot'] = this_fun(self.NNQ_Stats.loc[:,'boot'],NNQ2.NNQ_Stats.loc[:,'boot'])
            if 'Alphaboot' in self.NNQ_Stats and 'Alphaboot' in NNQ2.NNQ_Stats:
                result.NNQ_Stats.loc[:,'Alphaboot'] = this_fun(self.NNQ_Stats.loc[:,'Alphaboot'],NNQ2.NNQ_Stats.loc[:,'Alphaboot'])
        elif isinstance(NNQ2,TwoPointCorr):
            for (ip,it,itflow),iNNQ in self.NNQ_Stats['boot'].items():
                if (ip,it) not in NNQ2.C2_Stats.index: continue
                result[(ip,it,itflow)] = this_fun(iNNQ, NNQ2[(ip,it)])
            if 'Alphaboot' in self.NNQ_Stats:
                for (ip,it,itflow),iNNQ in self.NNQ_Stats['Alphaboot'].items():
                    if (ip,it) not in NNQ2.C2_Stats.index: continue
                    result.NNQ_Stats.loc[(ip,it,itflow),'Alphaboot'] = this_fun(iNNQ, NNQ2[(ip,it)])
        elif isinstance(NNQ2,FlowOp):
            for (ip,it,itflow),iNNQ in self.NNQ_Stats['boot'].items():
                if itflow not in NNQ2.Op_Stats.index: continue
                result[(ip,it,itflow)] = this_fun(iNNQ, NNQ2[itflow])
            if 'Alphaboot' in self.NNQ_Stats:
                for (ip,it,itflow),iNNQ in self.NNQ_Stats['Alphaboot'].items():
                    if itflow not in NNQ2.Op_Stats.index: continue
                    result.NNQ_Stats.loc[(ip,it,itflow),'Alphaboot'] = this_fun(iNNQ, NNQ2[itflow])
        else:
            try:
                result.NNQ_Stats.loc[:,'boot'] = this_fun(self.NNQ_Stats.loc[:,'boot'], NNQ2)
            except Exception as err:
                print(type(NNQ2), NNQ2)
                print(self.NNQ_Stats)
                raise EnvironmentError('Invalid value to combine with NNQCorr class')
            if 'Alphaboot' in self.NNQ_Stats:
                try:
                    result.NNQ_Stats.loc[:,'Alphaboot'] = this_fun(self.NNQ_Stats.loc[:,'Alphaboot'], NNQ2)
                except Exception as err:
                    print(type(NNQ2), NNQ2)
                    print(self.NNQ_Stats)
                    raise EnvironmentError('Invalid value to combine with NNQCorr class')
        return result


    def __add__(self,NNQ2):
        if any([ipar in str(type(NNQ2)) for ipar in NNQCorr.parentlist]):
            return NNQ2.__radd__(self)
        else:
            return self.Overload_wrap(NNQ2,op.add)

    def __sub__(self,NNQ2):
        if any([ipar in str(type(NNQ2)) for ipar in NNQCorr.parentlist]):
            return NNQ2.__rsub__(self)
        else:
            return self.Overload_wrap(NNQ2,op.sub)

    def __mul__(self,NNQ2):
        if any([ipar in str(type(NNQ2)) for ipar in NNQCorr.parentlist]):
            return NNQ2.__rmul__(self)
        else:
            return self.Overload_wrap(NNQ2,op.mul)

    def __div__(self,NNQ2):
        if any([ipar in str(type(NNQ2)) for ipar in NNQCorr.parentlist]):
            return NNQ2.__rdiv__(self)
        else:
            return self.Overload_wrap(NNQ2,op.truediv)

    def __pow__(self,NNQ2):
        if any([ipar in str(type(NNQ2)) for ipar in NNQCorr.parentlist]):
            return NNQ2.__rpow__(self)
        else:
            return self.Overload_wrap(NNQ2,op.pow)


    ## Right multiplication functions


    def __radd__(self,NNQ2):
        if any([ipar in str(type(NNQ2)) for ipar in NNQCorr.parentlist]):
            return NNQ2.__add__(self)
        else:
            return self.Overload_wrap(NNQ2,flip_fun_2arg(op.add))

    def __rsub__(self,NNQ2):
        if any([ipar in str(type(NNQ2)) for ipar in NNQCorr.parentlist]):
            return NNQ2.__sub__(self)
        else:
            return self.Overload_wrap(NNQ2,flip_fun_2arg(op.sub))

    def __rmul__(self,NNQ2):
        if any([ipar in str(type(NNQ2)) for ipar in NNQCorr.parentlist]):
            return NNQ2.__mul__(self)
        else:
            return self.Overload_wrap(NNQ2,flip_fun_2arg(op.mul))

    def __rdiv__(self,NNQ2):
        if any([ipar in str(type(NNQ2)) for ipar in NNQCorr.parentlist]):
            return NNQ2.__div__(self)
        else:
            return self.Overload_wrap(NNQ2,flip_fun_2arg(op.truediv))


    def __rpow__(self,NNQ2):
        if any([ipar in str(type(NNQ2)) for ipar in NNQCorr.parentlist]):
            return NNQ2.__pow__(self)
        else:
            return self.Overload_wrap(NNQ2,flip_fun_2arg(op.pow))



####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################


class NNQFullCorr(object):
    """ NNQFull(p,t,tf) uses bootstrap class
        p  = momentum
        t  = time
        tf = tflow
    """
    def Construct_NNQFull_File(this_file):
        return Construct_File_Object(this_file,NNQCorr)

    ## Info is a dictionary containing information for the correlator:
    ## see comments below to see what is required.
    ## missing elements just ignores that field
    ## All conflicts between TwoPointCorr and FlowOp were resolved by adding flow to variable names, and Flow to function names.
    ## All new functions and variables where I do not want to overwrite TwoPointCorr have the keyword NNQFull and NNQFull respectivly.

    ## TODO: create variational flag and interface with this function.

    ## array of all classes that inherit this class. Append this when adding new classes please!
    parentlist = ['SetOfCorrs.SetOfNNQFull']

    def __init__(self, thisnboot=nboot, cfglist={},
                 Info={},thissym='Not Set',thiscol='Not Set',thisshift=0.0,name='',man_load_cfgs=False):
        ## initialises all the data dictionaries


        ## Setting to g5P+NNQ (or g5P+NNW) by forceing Info to read correct projector and BarOrMes
        Info['Do_Chi'] = True
        # if 'g5' not in Info['Interp']: Info['Interp'] = 'CPOdd'
        Info['MesOrBar'] = 'Baryon'
        # c2_info = deepcopy(Info)
        # c2_info['pmom'] = 'All'
        ## Instance of FlowOp
        self.FOFull = FlowOpFullSquared(thisnboot,cfglist,Info,thissym,thiscol,thisshift,name)
        # ## Instance of TwoPointCorr
        # Info['fo_for_cfgs'] = self.FOFull.Observ
        self.NN = TwoPointCorr(thisnboot,cfglist,Info,thissym,thiscol,thisshift,name)

        ## these class elemnts are checked when loading in pickled file
        # self.thischecklist = [  'nboot','nt','ProjStr','ism','jsm','tsrc','kud','ks',
        #                         'NNfile','FOFullfile','tflowlist','boot_tsum_list']
        self.thischecklist = [  'nboot','boot_nt_list','ProjStr','ism','jsm','tsrc','kud','ks',
                                'tflowlist','boot_tsum_list','Fun']



        ## pull stuff from the NN and FO instances
        self.IsComb = False ## if it is a combination of correlators
        self.NNfile = self.NN.PickleFile
        self.FOFullfile = self.FOFull.flowPickleFile
        self.tflowlist = self.FOFull.tflowlist
        self.ProjStr = self.NN.ProjStr
        self.ism = self.NN.ism
        self.tsrc = self.NN.tsrc
        self.kud = self.NN.kud
        self.ks = self.NN.ks

        ## Grab some stuff (from .NN and .FOfull)
        self.stream_list = self.NN.stream_list
        self.nboot = thisnboot
        self.nt = self.NN.nt
        self.CPEvenFile = 'None'
        self.Info = Info
        self.thissym = thissym
        self.thiscol = thiscol
        self.thisshift = thisshift
        self.thisshiftscale = 'Not Set'
        self.name = name

        self.jsm = self.NN.jsm

        ## used to pick out different configuraiton/source lists.
        ## for example, if there are less flowed operators than 2 point correlators, passing 'TopCharge' will read that result
        if 'fo_for_cfgs' in list(Info.keys()): self.fo_for_cfgs = Info['fo_for_cfgs']
        else: self.fo_for_cfgs = False


        if 'Sparam' in list(Info.keys()): self.Sparam = Info['Sparam']
        else: self.Sparam = defSparam


        if 'show_cfgs' in list(Info.keys()): self.show_cfgs = Info['show_cfgs']
        else: self.show_cfgs = False


        ## Info can contain initial guess for fits, make sure it array of correct length
        if 'AlphaiGuess' in list(Info.keys()): self.AlphaiGuess = np.array(Info['AlphaiGuess'])
        else: self.AlphaiGuess = 'None' ## Defaults to value set for function in FitFunctions.py

        ## Info can contain initial guess for fits, make sure it array of correct length
        if 'NNQFulliGuess' in list(Info.keys()): self.NNQFulliGuess = np.array(Info['NNQFulliGuess'])
        else: self.NNQFulliGuess = 'None' ## Defaults to value set for function in FitFunctions.py


        ## Info can contain fit ranges for alpha to be calculated
        if 'FitsAlphaFull' in list(Info.keys()): self.PredefFits = Info['FitsAlphaFull'][0]
        else: self.PredefFits = []

        ## Info can contain default tflow picked, used for plotting
        if 'tflowfit' in list(Info.keys()): self.SetTflowFit(Info['tflowfit'])
        else: self.SetTflowFit('Not Set')


        ## Info can contain default t sink picked, used for plotting
        if 'tsinkfit' in list(Info.keys()): self.SetTFit(Info['tsinkfit'])
        else: self.SetTFit('Not Set')


        ## Info can contain default t sink picked, used for plotting
        if 'min_fit_len_tsum' in list(Info.keys()): self.min_fit_len_tsum = Info['min_fit_len_tsum']
        else: self.min_fit_len_tsum = 'Default'

        ## Info can contain fit ranges to be calculated
        ## must bein form Info['FitFunAlpha'] = (Function Object , number of parameters )
        # if 'FitFunAlpha' in Info.keys():
        #     thisFun,thisnpar = Info['FitFunAlpha']
        #     self.SetFunction(thisFun,thisnpar)
        # else:
        #     self.SetFunction('None',0)

        ## momentum ( in form of [px py pz] )
        self.checkmomlist = self.NN.checkmomlist
        self.pform = self.NN.pform

        ## selects the improved overlap, where the noise is subtracted from Q in NNQFull
        if 'IExpect' in list(Info.keys()):
            self.IExpect = Info['IExpect']
        else:
            self.IExpect = False

        if 'tsum_fit_list' in list(Info.keys()):
            self.tsum_fit_list = list(map(tsumstr,Info['tsum_fit_list']))
        else:
            self.tsum_fit_list = list(map(tsumstr,list(range(self.nt))))



        if 'boot_nt_list' in list(Info.keys()):
            self.boot_nt_list = list(map(untstr,Info['boot_nt_list']))
            hold_list = []
            for its in self.boot_nt_list:
                if -1 < its < self.nt:
                    hold_list.append(its)
            self.boot_nt_list = hold_list
        else:
            self.boot_nt_list = list(range(self.nt))





        '''
        I'm thinking types could be:
        from_src_sink_sum       is summing from source and sink, outwards
        from_src_sum            is summing from source
        from_sink_sum           is summing from sink
        from_center_sum         is summing from central point (tsink+tsrc)/2

        from_src_sym_sum                is summing symmetrically from source
        from_sink_sym_sum               is summing symmetrically from sink
        from_center_sym_sum             is summing symmetrically from central point (tsink+tsrc)/2

        from_src_sink           is values from source and sink, outwards
        from_src                is values from source
        from_sink               is values from sink
        from_center             is values from central point (tsink+tsrc)/2

        from_src_sym                is values symmetrically from source
        from_sink_sym               is values symmetrically from sink
        from_center_sym             is values symmetrically from central point (tsink+tsrc)/2

        mult   multiplys result by t (when doing non summed)
        etc...
        '''
        if 'op_trange_sum_type' in list(Info.keys()):
            self.tr_sum_type = Info['op_trange_sum_type']
        else:
            self.tr_sum_type = 'from_src_sink'
        self.tr_sum_name = self.tr_sum_type.replace('src','sr')
        self.tr_sum_name = self.tr_sum_name.replace('sink','si')
        self.tr_sum_name = self.tr_sum_name.replace('sum','su')
        self.tr_sum_name = self.tr_sum_name.replace('center','ce')
        self.tr_sum_name = self.tr_sum_name.replace('sym','sy')
        self.tr_sum_name = self.tr_sum_name.replace('mult','mt')
        self.tr_sum_name = self.tr_sum_name.replace('from','')

        if 'boot_tsum_list' in list(Info.keys()):
            self.boot_tsum_list = list(map(untstr,Info['boot_tsum_list']))
            hold_list = []
            for its in self.boot_tsum_list:
                if '_sym' in self.tr_sum_type and '_sum' in self.tr_sum_type:
                    if -1 < its < self.nt//2 + 2:
                        hold_list.append(its)
                else:
                    if -1 < its < self.nt:
                        hold_list.append(its)
            self.boot_tsum_list = hold_list
        else:
            self.boot_tsum_list = list(range(self.nt))

        self.is_Qsub = '_Qsub' in self.tr_sum_name
        self.is_Q2sub = '_Q2sub' in self.tr_sum_name

        self.is_NNQ2sub = '_NNQ2sub' in self.tr_sum_name
        self.is_NNQsub = '_NNQsub' in self.tr_sum_name
        self.is_P5NNQsub = '_P5NNQsub' in self.tr_sum_name
        self.is_NNSingQ2sub = '_NNSingQ2sub' in self.tr_sum_name

        self.is_improved =  self.is_P5NNQsub or self.is_NNQsub or self.is_Qsub or \
                            self.is_NNQ2sub or self.is_Q2sub or self.is_NNSingQ2sub



        if 'boot_spec_tsum_list' in list(Info.keys()) and self.is_NNQ2sub:
            self.boot_spec_tsum_list = list(map(untstr,Info['boot_spec_tsum_list']))
            hold_list = []
            for its in self.boot_spec_tsum_list:
                if '_sym' in self.tr_sum_type and '_sum' in self.tr_sum_type:
                    print('Warning, not sure how spectator Q is dealt with in sym and sum')
                    if -1 < its < self.nt//2 + 2:
                        hold_list.append(its)
                else:
                    if -1 < its < self.nt:
                        hold_list.append(its)
            self.boot_spec_tsum_list = hold_list
        else:
            self.boot_spec_tsum_list = list(range(self.nt))


        if self.is_NNQsub:
            this_info = deepcopy(Info)
            this_info['Interp'] = 'CPEven'
            self.NN2 = TwoPointCorr(thisnboot,cfglist,this_info,thissym,thiscol,thisshift,name)

        self.tr_sum_name = self.tr_sum_name.replace('_None','')
        self.tr_sum_type = self.tr_sum_type.replace('_None','')



        # self.NNQFullOp = ODNested()
        ## Series [ momentum , source_sink_sep , flow_time , configurations]
        self.NNQFull_cfgs = pa.DataFrame()
        '''
        self.NNQFull_cfgs DataFrame:

        columns =   xsrc_list , configs

        rows =      istream, iccfg
        '''

        if self.is_NNQ2sub:
            self.NNQ_subNNQ2_col_names = [ 'momentum','source_sink_separation','flow_time',
                                                'flow_op_trange','flow_op_2_trange']
            self.NNQ_subNNQ2_Stats = pa.DataFrame()

            '''
            self.NNQ_subNNQ2_Stats DataFrame:

            columns =   Alphaboot,  AlphaAvg,       AlphaStd
                        AlphaAuto,  AlphaAutoAvg,   AlphaAutoStd

            rows =      momentum,   source_sink_sep,    flow_time,  flow_op_trange, flow_op_2_trange
            '''

        self.NNQFull_col_names = [  'momentum','source_sink_separation','flow_time',
                                    'flow_op_trange']
        self.NNQFull_Stats = pa.DataFrame()

        '''
        self.NNQFull_Stats DataFrame:

        columns =   NNQFullboot,    NNQFullAvg,         NNQFullStd
                    Alphaboot,      AlphaAvg,           AlphaStd
                    AlphaAuto,      AlphaAutoAvg,       AlphaAutoStd

        rows =      momentum,   source_sink_sep,    flow_time,  flow_op_trange
        '''

        self.NNQInt_col_names = ['momentum','source_sink_separation','flow_time']
        self.NNQInt_Stats = pa.DataFrame()

        '''
        self.NNQInt_Stats DataFrame:

        columns =   NNQIntboot,     NNQIntAvg,          NNQIntStd
                    AlphaIntboot,   AlphaIntAvg,        AlphaIntStd
                    AlphaAuto,      AlphaAutoAvg,       AlphaAutoStd

        rows =      momentum,   source_sink_sep,    flow_time
        '''




        self.NNQFull_Fit_Stats = pa.DataFrame()
        self.NNQFull_Fit_col_names = ['momentum','flow_time']
        '''
        self.NNQFull_Fit_Stats DataFrame:

        columns =   NNQFullboot,    NNQFullAvg,         NNQFullStd
                    boot

        rows =      momentum,   flow_time, fit_range
        boot elements are SetsOfFitFuns instances
        '''

        # # FitDict = [ ip, itflow , Fit_range ] fitting class
        # self.FitDict    = ODNested()
        # # FitDictAvg/Std [ ip ,  itflow , Fit_range , Fit_Parameter ] bs
        # self.FitDictAvg = ODNested()
        # self.FitDictStd = ODNested()
        #
        # # FitDict = [ ip, itflow , Fit_range ] fitting class
        # self.NNQFullFitDict    = ODNested()
        # # FitDictAvg/Std [ ip ,  itflow , Fit_range , Fit_Parameter ] bs
        # self.NNQFullFitDictAvg = ODNested()
        # self.NNQFullFitDictStd = ODNested()

        ##
        # self.Tau_Fit_Stats = pa.DataFrame()
        # '''
        # self.Tau_Fit_Stats DataFrame:
        #
        # columns =   boot, Avg,      Std
        #
        # rows =      fit ranges , momentum , source sink sep , flow time
        #
        # '''
        self.SetFunction(ConstantFitFun,1)
        # self.SetFunction(Alpha_Prop_Fit_plin,4)
        # self.SetFunction(Alpha_Prop_Fit,3)
        # self.SetFunction(Alpha_Prop_Fit_2exp,5)


        # self.LoadPickle(DefWipe=DefWipe)
        # self.Read()
        self.SetCustomName(name)
        if not man_load_cfgs:
            self.LoadCfgs(cfglist,name)

    def LoadCfgs(self,cfglist,name=''):
        if len(cfglist) == 0:
            self.GetCfgList()
        else:
            self.ImportCfgList(cfglist)
        self.SetCustomName(name)

    ## to be consistant with TwoPtCorr;
    def GetCfgList(self):
        self.CombineCfgLists()

    def ImportPlotParams(self,thiscol,thissym,thisshift):
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        self.thisshift = thisshift

    def GetFuns(self):
        # self.NN.GetFuns()
        if self.NNQFull_Fit_Stats is not None and 'boot' in self.NNQFull_Fit_Stats.columns:
            for ival in self.NNQFull_Fit_Stats['boot'].values:
                ival.GetFuns()
        self.FOFull.GetFuns()
        if not hasattr(self,'Fun'):
            self.Fun = ReadFuns(self.Fun_name)[0]
        if hasattr(self,'NNCPEven'):
            self.NNCPEven.GetFuns()

    def RemoveFuns(self):
        # self.NN.RemoveFuns()
        if self.NNQFull_Fit_Stats is not None and 'boot' in self.NNQFull_Fit_Stats.columns:
            for ival in self.NNQFull_Fit_Stats['boot'].values:
                ival.RemoveFuns()
        self.FOFull.RemoveFuns()
        if hasattr(self,'Fun'):
            self.Fun_name = self.Fun.__name__
            WriteFuns(self.Fun)
            del self.Fun
        if hasattr(self,'NNCPEven'):
            self.NNCPEven.RemoveFuns()

    def Get_Extrapolation(self,keep_keys=slice(None),fmted=True):
        if 'boot' in self.NNQFull_Fit_Stats:
            output =  sff.PullSOFSeries_Wpar(  self.NNQFull_Fit_Stats['boot'],
                                                keep_keys=keep_keys,fmted=fmted)
            return Series_fix_key(output,4,'fittwor','tsumfitr')
        else:
            return pa.Series()

    def CombineCfgLists(self):
        if 'configs' not in self.FOFull.Op_cfgs:
            if 'configs' in self.NN.C2_cfgs:
                self.ImportCfgList(self.NN.C2_cfgs)
            return
        if 'configs' not in self.NN.C2_cfgs:
            self.ImportCfgList(self.FOFull.Op_cfgs['configs'])
            return
        self.ImportCfgList(CombineCfgs(self.FOFull.Op_cfgs,self.NN.C2_cfgs))

        # FO_cfgs = self.FOFull.Op_cfgs['configs']
        # comb_cfgs,comb_xsrc = [],[]
        # ilist = []
        # for (istream,iccfg),NN_cfg_data in self.NN.C2_cfgs.iterrows():
        #     if NN_cfg_data['configs'] in FO_cfgs[istream].values:
        #         comb_cfgs.append(NN_cfg_data['configs'])
        #         comb_xsrc.append(NN_cfg_data['xsrc_list'])
        #         ilist.append((istream,iccfg))
        # indicies = pa.MultiIndex.from_tuples(ilist,names=['stream','config_number'])
        # cfglist = pa.Series(comb_cfgs,index=indicies).to_frame(name='configs')
        # cfglist.loc[:,'xsrc_list'] = pa.Series(comb_xsrc,index=indicies)
        # self.ImportCfgList(cfglist)

    def Stats(self,DefWipeG2=False,DoAuto=True):
        def DoStats(val):
            if isinstance(val,BootStrap): val.Stats()
            return val
        def DoAvg(val):
            if isinstance(val,BootStrap):
                return val.Avg
            else:
                return val
        def DoStd(val):
            if isinstance(val,BootStrap):
                return val.Std
            else:
                return val
        if not self.IsComb:
            self.NN.Stats()
            self.FOFull.FlowStats()
        ## maybe effmass and alpha not here???
        self.NNQFull_Stats.loc[:,'boot'] = self.NNQFull_Stats['boot'].apply(DoStats)
        self.NNQFull_Stats.loc[:,'Avg'] = self.NNQFull_Stats['boot'].apply(DoAvg)
        self.NNQFull_Stats.loc[:,'Std'] = self.NNQFull_Stats['boot'].apply(DoStd)
        if 'Alphaboot' in self.NNQFull_Stats:
            self.NNQFull_Stats.loc[:,'Alphaboot'] = self.NNQFull_Stats['Alphaboot'].apply(DoStats)
            self.NNQFull_Stats.loc[:,'AlphaAvg'] = self.NNQFull_Stats['Alphaboot'].apply(DoAvg)
            self.NNQFull_Stats.loc[:,'AlphaStd'] = self.NNQFull_Stats['Alphaboot'].apply(DoStd)
        self.NNQInt_Stats.loc[:,'boot'] = self.NNQInt_Stats['boot'].apply(DoStats)
        self.NNQInt_Stats.loc[:,'Avg'] = self.NNQInt_Stats['boot'].apply(DoAvg)
        self.NNQInt_Stats.loc[:,'Std'] = self.NNQInt_Stats['boot'].apply(DoStd)
        # self.EffMass()
        self.AlphaRatio(DefWipeG2=DefWipeG2,DoAuto=DoAuto)

    ## not needed, taken care of when bootstrapping
    # def NNQFullCombPreBS(self):
    #     self.NNQFull = np.array([np.multiply.outer(iC2,Op) for iC2,iOp in zip(self.C2,self.Op)])

    '''
    I'm thinking types could be:
        from_src_sink_sum       is summing from source and sink, outwards
        from_src_sum            is summing from source
        from_sink_sum           is summing from sink
        from_center_sum         is summing from central point (tsink+tsrc)/2

        from_src_sym_sum                is summing symmetrically from source
        from_sink_sym_sum               is summing symmetrically from sink
        from_center_sym_sum             is summing symmetrically from central point (tsink+tsrc)/2

        from_src_sink           is values from source and sink, outwards
        from_src                is values from source
        from_sink               is values from sink
        from_center             is values from central point (tsink+tsrc)/2

        from_src_sym                is values symmetrically from source
        from_sink_sym               is values symmetrically from sink
        from_center_sym             is values symmetrically from central point (tsink+tsrc)/2

        etc...
    '''
    def Create_Tsrc_Fun(self,ictflow,it_sink):
        if isinstance(it_sink,str):
            it_sink = untstr(it_sink.replace('tsink','t'))
        # if it_sink > self.nt/2:
        #     t_mid = (self.nt-it_sink)/2
        # else:
        #     t_mid = it_sink/2
        # if 'from_src_sink' in self.tr_sum_type:
        #     def Tsrc_Fun(this_series):
        #         ## tsrc_list = [ ixsrc ]
        #         out_list = []
        #         for itsrc in this_series['tsrc_list']:
        #             ## op_series = [ tflow, t_src ]
        #             flow_op = this_series['Op'][ictflow]
        #             tflow_src = np.roll(flow_op,-itsrc)
        #             tflow_sink = np.roll(tflow_src,-it_sink)
        #             tflow_src_left = tflow_src[::-1]
        #             tflow_sink_left = tflow_sink[::-1]
        #             tflow_src_left[0] = np.float64(0.0)
        #             tflow_sink_left[0] = np.float64(0.0)
        #             if it_sink > self.nt/2:
        #                 # tflow_src_left[:t_mid] = tflow_sink[:t_mid]= np.float64(0.0)
        #                 this_out = (np.lib.pad(tflow_src_left[:t_mid],(0,self.nt-t_mid),'constant',constant_values=(0,0.))
        #                             + tflow_src + tflow_sink_left +
        #                             np.lib.pad(tflow_sink[:t_mid],(0,self.nt-t_mid),'constant',constant_values=(0,0.)))
        #             else:
        #                 # tflow_src[:t_mid] = tflow_sink_left[:t_mid] = np.float64(0.0)
        #                 this_out = (np.lib.pad(tflow_src[:t_mid],(0,self.nt-t_mid),'constant',constant_values=(0,0.))
        #                             + tflow_src_left + tflow_sink +
        #                             np.lib.pad(tflow_sink_left[:t_mid],(0,self.nt-t_mid),'constant',constant_values=(0,0.)))
        #             if '_sum' in self.tr_sum_type:
        #                 denom_list = [2.] + [2.+ival*4 for ival in range(1,t_mid+1)]
        #                 last_denom = denom_list[-1]
        #                 denom_list = denom_list + [ival*2 + last_denom for ival in range(1,self.nt-t_mid)]
        #                 this_out = np.cumsum(this_out)
        #                 out_list.append(this_out*(self.nt/np.array(denom_list,dtype=np.float64)))
        #             else:
        #                 out_list.append(this_out*self.nt)
        #         return pa.Series([out_list],index=['iFO'])
        # else:
        def Tsrc_Fun(this_series):
            ## tsrc_list = [ ixsrc ]
            out_list = []
            for itsrc in this_series['tsrc_list']:
                if 'from_center' in self.tr_sum_type:
                    n_roll = (itsrc+it_sink//2)%self.nt
                elif 'src' in self.tr_sum_type:
                    n_roll = itsrc
                elif 'from_sink' in self.tr_sum_type:
                    n_roll = (itsrc+it_sink)%self.nt
                else:
                    raise EnvironmentError(self.tr_sum_type + ' sum time not supported.')
                ## op_series = [ tflow, t_src ]
                flow_op = this_series['Op'][ictflow]
                tflow_rolled = np.roll(flow_op,-n_roll)
                if '_src_sink' in self.tr_sum_type:
                    tflow_rolled_tsink = np.roll(np.roll(tflow_rolled,-it_sink)[::-1],1)
                    if '_sym' in self.tr_sum_type:
                        tflow_rolled = tflow_rolled+np.roll(tflow_rolled[::-1],1)
                        tflow_rolled_tsink = tflow_rolled_tsink+tflow_rolled_tsink[::-1]
                    elif '_rev' in self.tr_sum_type:
                        tflow_rolled = tflow_rolled[::-1]
                        tflow_rolled_tsink = tflow_rolled_tsink[::-1]
                    tflow_rolled = (tflow_rolled + tflow_rolled_tsink)/2.
                elif '_rev' in self.tr_sum_type:
                    tflow_rolled = np.roll(tflow_rolled[::-1],1)
                elif '_sym' in self.tr_sum_type:
                    tflow_rolled = tflow_rolled+np.roll(tflow_rolled[::-1],1)
                if '_sum' in self.tr_sum_type:
                    if '_sym' in self.tr_sum_type and '_src_sink' not in self.tr_sum_type:
                        tflow_rolled[0] = tflow_rolled[0]/2.
                        if self.nt//2 == self.nt/2:
                            tflow_rolled[int(self.nt//2-1)] = tflow_rolled[int(self.nt//2-1)]/2.
                    out_list.append(np.cumsum(tflow_rolled))
                elif '_mult' in self.tr_sum_type:
                    out_list.append(tflow_rolled*self.nt)
                else:
                    out_list.append(tflow_rolled)
            return pa.Series([out_list],index=['iFO'])
        return Tsrc_Fun

    def Create_NNQ_fun(self,it_sum):
        if '_FO2' in self.tr_sum_type:
            if '_mult' in self.tr_sum_type:
                scale_coeff = np.float(self.nt)
            elif '_sum' in self.tr_sum_type:
                if '_sym' in self.tr_sum_type:
                    scale_coeff = np.float(it_sum)*2+1
                    if it_sum == self.nt//2-1:
                        scale_coeff -= 1.0
                else:
                    scale_coeff = np.float(it_sum)+1
            else:
                scale_coeff = 1.0
            def NNQ_fun(this_df):
                return np.mean(this_df['iC2']*(np.array(this_df['iFO'])[:,it_sum]**2)/scale_coeff)
        else:
            def NNQ_fun(this_df):
                return np.mean(this_df['iC2']*np.array(this_df['iFO'])[:,it_sum])
        return NNQ_fun

    def Create_NNQ_CPEven_fun(self,it_sum):
        def NNQ_fun(this_df):
            return np.mean(this_df['iC2_CPEven']*np.array(this_df['iFO'])[:,it_sum])
        return NNQ_fun


    def Create_NNQ2_fun(self,it_sum,jt_sum):
        if '_mult' in self.tr_sum_type:
            scale_coeff = np.float(self.nt)
        elif '_sum' in self.tr_sum_type:
            if '_sym' in self.tr_sum_type:
                scale_coeff = np.float(jt_sum)*2+1
                if jt_sum == self.nt//2-1:
                    scale_coeff -= 1.0
            else:
                scale_coeff = np.float(jt_sum)+1
        else:
            scale_coeff = 1.0
        def NNQ2_fun(this_df):
            return  np.mean(this_df['iC2']* \
                    np.array(this_df['iFO'])[:,it_sum]* \
                    np.array(this_df['iFO'])[:,jt_sum])/scale_coeff
        return NNQ2_fun
    ## DefWipe can be 'C2' or 'Flow' to only delete that specific data
    ## You can even call self.FlowCombinePreBS(OtherOpp) if you wanted to look at other flowed operators too!
    def Bootstrap(self,WipeData = False,tflowlist = 'PreDef',show_timer=True,DoAuto=True):
        # self.C2 = [ ip , it , iconf ]
        # self.Op = [ itflow , iconf ]
        if tflowlist != 'PreDef': self.tflowlist = tflowlist
        ## Do not wipe the data, we need it for combining on the configuration level!
        # self.NN.Bootstrap(WipeData=False)
        # self.FOFull.FlowBootstrap(tflowlist=tflowlist,WipeData=False,show_timer=False)
        # ReadAll = (len(self.tflowlist) == 0)

        if DoAuto:
            NNQFulllcfg,NNQIntlcfg = [],[]
        if self.is_improved:
            Qcfg,QIntcfg = [],[]
        if self.is_P5NNQsub:
            P5NNcfg,P5NNIntcfg = [],[]
        NNQFulll,NNQFulllAvg,NNQFulllStd = [],[],[]
        NNQIntl,NNQIntlAvg,NNQIntlStd = [],[],[]
        indexintl = []
        indexl = []
        if show_timer: thistimer = Timer(   linklist=[self.pform,self.boot_nt_list,
                                                        self.FOFull.tflowlist,
                                                        self.boot_tsum_list],
                                            name='Booting '+self.name)
        self.FOFull.Op_cfgs['tsrc_list'] = self.NN.C2_cfgs['tsrc_list']
        for icp,ip in enumerate(self.pform):
            ic2pt = self.NN.pform.index(ip)
            if DoAuto:
                NNQFulllcfg.append([])
                NNQIntlcfg.append([])
            if self.is_improved:
                Qcfg.append([])
                QIntcfg.append([])
            if self.is_P5NNQsub:
                P5NNcfg.append([])
                P5NNIntcfg.append([])
            for it in self.boot_nt_list:
                if DoAuto:
                    NNQFulllcfg[-1].append([])
                    NNQIntlcfg[-1].append([])
                if self.is_improved:
                    Qcfg[-1].append([])
                    QIntcfg[-1].append([])
                if self.is_P5NNQsub:
                    P5NNcfg[-1].append([])
                    P5NNIntcfg[-1].append([])
                strit = tstr(it+1)
                def C2_ip_it_fun(x):
                    return np.array(x[ic2pt][it])
                tdata = self.NNQFull_cfgs['C2'].apply(C2_ip_it_fun).to_frame('iC2')
                tdata_xsrc_avg = tdata['iC2'].apply(lambda x : np.mean(x))
                if self.is_NNQsub:
                    tdata.loc[:,'iC2_CPEven'] = self.NNQFull_cfgs['C2_2'].apply(C2_ip_it_fun)
                    tdata_xsrc_avg2 = tdata['iC2_CPEven'].apply(lambda x : np.mean(x))

                # tdata = self.NN.C2_cfgs['C2'].apply(lambda x : np.array(x[icp][it])).to_frame('iC2')
                for ictflow,tflows in enumerate(self.FOFull.tflowlist):
                    strtflows = tflowstr(tflows)
                    op_int = self.NNQFull_cfgs['Op'].apply(lambda x : np.sum(x[ictflow]))
                    int_data = tdata_xsrc_avg * op_int
                    if '_FO2' in self.tr_sum_type:
                        int_data = int_data * op_int
                    if DoAuto:
                        NNQFulllcfg[-1][-1].append([])
                        NNQIntlcfg[-1][-1].append(int_data.values)
                    if self.is_improved:
                        Qcfg[-1][-1].append([])
                        if self.is_P5NNQsub:
                            P5NNcfg[-1][-1].append([])
                            P5NNIntcfg[-1][-1].append(tdata_xsrc_avg.values)
                            QIntcfg[-1][-1].append(op_int.values)
                        if self.is_Qsub:
                            QIntcfg[-1][-1].append(op_int.values)
                        elif self.is_Q2sub:
                            Q2list = op_int.values**2
                            QIntcfg[-1][-1].append(Q2list-np.mean(Q2list))
                        elif self.is_NNQsub:
                            NNQlist = (tdata_xsrc_avg2 * op_int).values
                            QIntcfg[-1][-1].append(NNQlist)
                        elif self.is_NNQ2sub or self.is_NNSingQ2sub:
                            NNQ2list = (tdata_xsrc_avg * op_int**2/np.float(self.nt)).values
                            QIntcfg[-1][-1].append(NNQ2list)


                    Intboot = BootStrap(self.nboot, name='G_{2}('+','.join([ip,strit,strtflows])+')',cfgvals=int_data)
                    indexintl.append((ip,strit,strtflows))
                    NNQIntl.append(Intboot)
                    NNQIntlAvg.append(Intboot.Avg)
                    NNQIntlStd.append(Intboot.Std)
                    # print self.tflowlist , strtflows
                    # print 'test',ip,strit,strtflows

                    tdata['iFO'] = self.NNQFull_cfgs[['tsrc_list','Op']].apply(self.Create_Tsrc_Fun(ictflow,it),axis=1)['iFO']
                    for it_sum in self.boot_tsum_list:
                        tsum_str = tstr(it_sum).replace('t','ts')
                        thisdata = tdata.apply(self.Create_NNQ_fun(it_sum),axis=1)
                        if DoAuto:
                            NNQFulllcfg[-1][-1][-1].append(thisdata.values)
                        if self.is_P5NNQsub:
                            P5NN_vals = tdata['iC2'].apply(lambda x: np.mean(np.array(x))).values
                            Q_vals = tdata['iFO'].apply(lambda x: np.mean(np.array(x)[:,it_sum])).values
                            # Qcfg[-1][-1][-1].append(Q_vals/np.float(self.nt))
                            P5NNcfg[-1][-1][-1].append(P5NN_vals)
                            Qcfg[-1][-1][-1].append(Q_vals)
                        if self.is_Qsub:
                            Q_vals = tdata['iFO'].apply(lambda x: np.mean(np.array(x)[:,it_sum])).values
                            # Qcfg[-1][-1][-1].append(Q_vals/np.float(self.nt))
                            Qcfg[-1][-1][-1].append(Q_vals/np.float(self.nt))
                        elif self.is_Q2sub:
                            def ThisFun(x):
                                return np.mean(np.array(x)[:,it_sum]**2)
                            Q2_vals = tdata['iFO'].apply(ThisFun).values
                            Qcfg[-1][-1][-1].append( Q2_vals - np.mean(Q2_vals))
                        elif self.is_NNQsub:
                            this_NNQ = tdata.apply(self.Create_NNQ_CPEven_fun(it_sum),axis=1)
                            Qcfg[-1][-1][-1].append(this_NNQ.values/np.float(self.nt) )
                        elif self.is_NNSingQ2sub:
                            this_NNQ2 = tdata.apply(self.Create_NNQ2_fun(it_sum,it_sum),axis=1)
                            Qcfg[-1][-1][-1].append(this_NNQ2.values )
                        elif self.is_NNQ2sub:
                            Qcfg[-1][-1][-1].append([])
                            for jt_sum in self.boot_spec_tsum_list:
                                this_NNQ2 = tdata.apply(self.Create_NNQ2_fun(it_sum,jt_sum),axis=1)
                                Qcfg[-1][-1][-1][-1].append(this_NNQ2.values )
                        indexl.append((ip,strit,strtflows,tsum_str))
                        # for icfg,icfg_data in zip(self.FOFull.flowcfglist,thisdata):

                        if hasattr(self,'IExpect') and self.IExpect:
                            raise EnvironmentError('IExpect not implemented after pandas rework. on the TODO list!')
                            # if isinstance(self.IExpect,float):
                            #     self.NNQFullboot[ip][strit][strtflows] = BootStrap(self.nboot, name='G_{2}('+','.join([ip,strit,strtflows])+')',cfgvals=tdata,
                            #                                                    thiseps=self.IExpect,QforImprov=tflowdata)
                            #     # thisfile = './temp/'+'_'([ip,strit,strtflows,epsstr(self.IExpect)])+'.py3p'
                            #     # if os.path.isfile(thisfile):
                            #     #     with open(thisfile,'rb') as f:
                            #     #         self.NNQFullboot[ip][strit][strtflows] = pickle.load(f)
                            #     # else:
                            #     #     for ieps in np.arange(0.1,self.IExpect,0.1):
                            #     #         mkdir_p('./temp')
                            #     #         with open(thisfile.replace(epsstr(self.IExpect),epsstr(ieps)),'wb') as f:
                            #     #             pickle.dump(BootStrap(self.nboot, name='G_{2}('+','.join([ip,strit,strtflows])+')',data=tdata,
                            #     #                                                thiseps=ieps,QforImprov=tflowdata),f)
                            # elif self.IExpect == True:
                            #     self.NNQFullboot[ip][strit][strtflows] = BootStrap(self.nboot, name='G_{2}('+','.join([ip,strit,strtflows])+')',cfgvals=tdata,QforImprov=tflowdata)
                            # elif self.IExpect == False:
                            #     self.NNQFullboot[ip][strit][strtflows] = BootStrap(self.nboot, name='G_{2}('+','.join([ip,strit,strtflows])+')',cfgvals=thisdata)
                        else:
                            thisboot = BootStrap(self.nboot, name='G_{2}('+','.join([ip,strit,strtflows,tsum_str])+')',cfgvals=thisdata)
                            NNQFulll.append(thisboot)
                            NNQFulllAvg.append(thisboot.Avg)
                            NNQFulllStd.append(thisboot.Std)
                        # thistimer.Lap()
                        if show_timer: thistimer.Lap((ip,strit,strtflows,tsum_str))
        if len(indexl) > 0:
            indicies = pa.MultiIndex.from_tuples(indexl,names=self.NNQFull_col_names)
            if DoAuto:
                self.NNQFull_cfgs.loc[:,'NNQFull'] = pa.Series(np.rollaxis(np.array(NNQFulllcfg),-1).tolist(),index=self.NNQFull_cfgs.index)
                if self.is_Qsub or self.is_Q2sub or self.is_NNSingQ2sub or self.is_NNQsub or self.is_P5NNQsub:
                    self.NNQFull_cfgs.loc[:,'FlowOp'] = pa.Series(np.rollaxis(np.array(Qcfg),-1).tolist(),index=self.NNQFull_cfgs.index)
                    if self.is_P5NNQsub:
                        self.NNQFull_cfgs.loc[:,'P5C2'] = pa.Series(np.rollaxis(np.array(P5NNcfg),-1).tolist(),index=self.NNQFull_cfgs.index)
                elif self.is_NNQ2sub:
                    self.NNQFull_cfgs.loc[:,'NNQ2'] = pa.Series(np.rollaxis(np.array(Qcfg),-1).tolist(),index=self.NNQFull_cfgs.index)

            self.NNQFull_Stats.loc[:,'boot'] = pa.Series(NNQFulll,index=indicies)
            self.NNQFull_Stats.loc[:,'Avg'] = pa.Series(NNQFulllAvg,index=indicies)
            self.NNQFull_Stats.loc[:,'Std'] = pa.Series(NNQFulllStd,index=indicies)
        if len(indexintl) > 0:
            indicies = pa.MultiIndex.from_tuples(indexintl,names=self.NNQInt_col_names)
            if DoAuto:
                self.NNQFull_cfgs.loc[:,'NNQInt'] = pa.Series(np.rollaxis(np.array(NNQIntlcfg),-1).tolist(),index=self.NNQFull_cfgs.index)
                if self.is_Qsub or self.is_Q2sub or self.is_NNSingQ2sub or self.is_NNQsub or self.is_P5NNQsub:
                    self.NNQFull_cfgs.loc[:,'FlowOpInt'] = pa.Series(np.rollaxis(np.array(QIntcfg),-1).tolist(),index=self.NNQFull_cfgs.index)
                    if self.is_P5NNQsub:
                        self.NNQFull_cfgs.loc[:,'C2Int'] = pa.Series(np.rollaxis(np.array(P5NNIntcfg),-1).tolist(),index=self.NNQFull_cfgs.index)
                elif self.is_NNQ2sub:
                    self.NNQFull_cfgs.loc[:,'NNQ2Int'] = pa.Series(np.rollaxis(np.array(QIntcfg),-1).tolist(),index=self.NNQFull_cfgs.index)

            self.NNQInt_Stats.loc[:,'boot'] = pa.Series(NNQIntl,index=indicies)
            self.NNQInt_Stats.loc[:,'Avg'] = pa.Series(NNQIntlAvg,index=indicies)
            self.NNQInt_Stats.loc[:,'Std'] = pa.Series(NNQIntlStd,index=indicies)
        self.RemoveVals(WipeData = WipeData)

    def ObsToName(self,thisop):
        thisobs = self.FOFull.FlowObsToName(thisop)
        return thisobs.replace('<',r'\alpha_{').replace('>',r'}')

    # ## thisInfo must contain correct information for combining the correlators (the smearings, the eigenvectors, etc...)
    # def GetCombCorr(self,thisInfo):
    #     from SetsOfTwoPtCorrs import SetOfNNQFull
    #     sinklab,jsmcoeffs,sm2ptList,InfoList = GetInfoList(thisInfo)

    #     C2Set = SetOfNNQFull(cfglist=self.NNQFull_cfgs[['configs','xsrc_list']],InfoDict = InfoList)
    #     C2Set.LoadPickle(CheckCfgs=CheckCfgs)
    #     thisname,thisleglab,jsmlab = CreateNameAndLL(C2Set.SetC2,sinklab)
    #     thisfun = GenSinkFun(jsmcoeffs)
    #     thisfile = self.filename.replace(self.jsm,sinklab)
    #     thisLL = self.LegLab.replace(self.jsm,sinklab)
    #     # self.SetCustomName(thisfile,thisLL)
    #     # self.jsm = jsmlab
    #     return C2Set.CombineCorrs(thisfun,filename=thisfile,LegLab=thisLL,jsmlab=jsmlab)


    def SetCustomName(self,string='',stringLL='',fo_for_cfgs='PreDef'):
        if string == '':
            self.name = '_'.join([  self.NN.dim_label,self.NN.kappafolder,
                                    self.NN.stream_name,self.NN.tsrc,self.NN.ism,self.jsm,
                                    self.FOFull.Observ])+self.tr_sum_name
            self.filename = '_'.join([self.NN.MesOrBar,self.NN.stream_name,self.ProjStr,self.NN.tsrc,self.NN.ism,self.jsm,self.FOFull.Observ])+self.tr_sum_name
            if 'CPOdd' != self.ProjStr:
                self.name = self.name + '_'+self.ProjStr
                self.filename = self.filename + '_'+self.ProjStr
            if hasattr(self,'IExpect'):
                if isinstance(self.IExpect,float):
                    self.filename += '_Ieps'+epsstr(self.IExpect)
                elif self.IExpect:
                    self.filename += '_Ieps'+epsstr(Qeps)
        else:
            self.filename = self.name = string
        if stringLL == '':
            self.LegLab = ('$'+'\ '.join([self.NN.dim_ll,self.NN.latparams.GetPionMassLab(),self.NN.stream_name,
                                        self.jsm,self.ObsToName(self.FOFull.Observ)])
                                        +self.tr_sum_name.replace('_','\ ')+'$') ## customise this how you like
            if 'CPOdd' != self.ProjStr:
                self.LegLab = self.LegLab + ' $'+self.ProjStr+'$'
            if hasattr(self,'IExpect'):
                if isinstance(self.IExpect,float):
                    self.LegLab += ' Ieps'+epsstr(self.IExpect)
                elif self.IExpect:
                    self.LegLab += ' Ieps'+epsstr(Qeps)
        else:
            self.LegLab = stringLL
        if isinstance(fo_for_cfgs,str):
            if fo_for_cfgs == 'PreDef':
                fo_for_cfgs = self.fo_for_cfgs
        if isinstance(fo_for_cfgs,str):
            self.NNQFulldir = outputdir + '/'+self.NN.dim_label+self.NN.kappafolder+'/alpha/'+fo_for_cfgs+'/'
        else:
            self.NNQFulldir = outputdir + '/'+self.NN.dim_label+self.NN.kappafolder+'/alpha/'
        mkdir_p(self.NNQFulldir+'/Pickle/')
        mkdir_p(self.NNQFulldir+'/Excel/')
        self.HumanFile = self.NNQFulldir+self.filename+'.xml'
        self.ExcelFile = self.NNQFulldir+'/Excel/'+self.filename+'.xlsx'
        self.PickleFile = self.NNQFulldir+'/Pickle/'+self.filename+'.py3p'


    def ImportCfgList(self,cfglist):
        if 'configs' in cfglist:
            self.FOFull.FlowImportCfgList(cfglist)
            self.NN.ImportCfgList(cfglist)
            self.NNQFull_cfgs = self.NN.C2_cfgs
            # self.NNQFull_cfgs['Op'] = pa.Series(self.FOFull.Op_cfgs.loc[:,'Op'].values,self.NNQFull_cfgs.index)
            self.Update_Stream_List()

    def Update_Stream_List(self):
        self.stream_list = self.NN.stream_list
        self.ncfg_list = self.NN.ncfg_list
        self.xsrcAvg = self.NN.xsrcAvg
        self.nMeas = self.NN.nMeas
        self.nstream = len(self.stream_list)
        if hasattr(self,'stream_name'):
            old_sn = self.stream_name
            self.stream_name = '-'+''.join([istr.replace('-','') for istr in self.stream_list])+'-'
            if hasattr(self,'filename'):
                self.filename = self.filename.replace(old_sn,self.stream_name)
            if hasattr(self,'name'):
                self.name = self.name.replace(old_sn,self.stream_name)
            if hasattr(self,'LegLab'):
                self.LegLab = self.LegLab.replace(old_sn,self.stream_name)
        else:
            self.stream_name = '-'+''.join([istr.replace('-','') for istr in self.stream_list])+'-'



    def RemoveVals(self,WipeData = True):
        if 'C2' == WipeData:
            self.NN.RemoveVals()
            if 'C2' in self.NNQFull_cfgs: del self.NNQFull_cfgs['C2']
            if 'C2_2' in self.NNQFull_cfgs: del self.NNQFull_cfgs['C2_2']
        elif 'Flow' ==  WipeData:
            self.FOFull.FlowRemoveVals()
            if 'Op' in self.NNQFull_cfgs: del self.NNQFull_cfgs['Op']
        elif WipeData:
            self.RemoveVals(WipeData='C2')
            self.RemoveVals(WipeData='Flow')

    ## self.C2 = [ ip, it , iconf ]
    def Read(self):
        self.NN.Read(Full=True)
        self.FOFull.FlowRead()
        self.NNQFull_cfgs = self.NN.C2_cfgs
        if self.is_NNQsub:
            self.NN2.Read(Full=True)
            self.NNQFull_cfgs.loc[:,'C2_2'] = pa.Series(self.NN2.C2_cfgs.loc[:,'C2'].values,index=self.NNQFull_cfgs.index)
        self.NNQFull_cfgs.loc[:,'Op'] = pa.Series(self.FOFull.Op_cfgs.loc[:,'Op'].values,index=self.NNQFull_cfgs.index)
        if not hasattr(self,'ncfg_list'):
            self.ncfg_list = [icfg.size for istream,icfg in self.NNQFull_cfgs['configs'].groupby(level='stream')]


    # def ReadAndBoot(self):
    #     self.Read()
    #     self.Bootstrap()

    # def ImportFitRanges(self,fit_range):
    #     self.fit_range = fit_range
    #     if isinstance(self.fit_range,basestring):
    #         self.fit_range = [self.fit_range]


    def CheckCol(self,thiscol):
        if 'PreDefine' == thiscol :
            if 'Not Set' == self.thiscol:
                raise IOError('Pass in color to initialize it')
        else:
            self.thiscol = thiscol

    def CheckSym(self,thissym):
        if 'PreDefine' == thissym :
            if 'Not Set' == self.thissym:
                raise IOError('Pass in symbol to initialize it')
        else:
            self.thissym = thissym

    def GetShift(self,xlims,thisshift):
        if thisshift != 'PreDefine': self.thisshift = thisshift
        # if self.thisshiftscale == 'Not Set':
        xlen = np.abs(xlims[1]-xlims[0])
        self.thisshiftscale = self.thisshift*xlen
        return self.thisshiftscale
        # else:
        #     return self.thisshiftscale

    def SetTflowFit(self,tflow):
        # if tflow not in self.tflowlist and tflow != 'Not Set':
        #     print 'Warning: ', tflow, ' not in flow list, tflowfit not added'
        # else:
        self.tflowfit = tflow

    def SetTFit(self,tsink):
        # if tsink not in ['t'+str(it) for it in xrange(1,self.nt+1)] and tsink != 'Not Set':
        #     print 'Warning: ', tsink, ' not in  list, tsinkfit not added'
        # else:
        self.tsinkfit = tsink

    def SetFunction(self,Fun,npar):
        self.Fun = Fun
        self.npar = npar



    # def Fit_TauDep(self,tau_min,tau_max,imom,itsink,itflow,tau_fit_fun='PreDef'):
    #     if tau_fit_fun != 'PreDef': self.tau_fit_fun = tau_fit_fun
    #     fitr = 'taufitr'+str(untstr(tau_min))+'-'+str(untstr(tau_max))
    #     tau_fit_name = '_'.join([fitr,imom,itsink,itflow])
    #     fit_data = self.NNQFull_Stats['Alphaboot'][imom,itsink,itflow]
    #     fit_data = fit_data[tau_min:tau_max]
    #     parse_data = pa.DataFrame(columns=['ydata','xdata'],index=[0])
    #     parse_data['ydata'][0] = fit_data.values
    #     parse_data['xdata'][0] = [map(untstr,list(fit_data.index)),[self.nt]*len(fit_data.index)]
    #     fit_class = ff.Fitting(Funs=self.tau_fit_fun,data=parse_data,name=tau_fit_name)
    #     fit_class.FitBoots()
    #     return fit_class
    #
    # def Fit_All_Tau(self,fit_range,mom_list='All',tsink_list='All',tflow_list='All',tau_fit_fun='PreDef'):
    #     flist,aflist,sflist,index_list = [],[],[],[]
    #     if mom_list == 'All': mom_list = self.pform
    #     if tsink_list == 'All': tsink_list = ['t'+str(it+1) for it in self.boot_nt_list]
    #     if tflow_list == 'All': tflow_list = self.tflowlist
    #     for ifitr,imom,itsink,itflow in itertools.product(fit_range,mom_list,tsink_list,tflow_list):
    #         if isinstance(ifitr,basestring):
    #             ifit_min,ifit_max = map(int,ifitr.replace('taufitr','').split('-'))
    #         else:
    #             ifit_min,ifit_max = ifitr
    #         flist.append(self.Fit_TauDep(ifit_min,ifit_max,imom,itsink,itflow,tau_fit_fun))
    #         aflist.append(flist[-1].fit_data['ParamsAvg'])
    #         sflist.append(flist[-1].fit_data['ParamsStd'])
    #         ikey = flist[-1].name.split('_')
    #         index_list.append(ikey[:-2]+ [ikey[-2]+'_'+ikey[-1]])
    #     indicies = pa.MultiIndex.from_tuples(index_list,names=['fit_range','momentum','source_sink_sep','flow_times'])
    #     if len(index_list) > 0:
    #         if 'boot' in self.Tau_Fit_Stats.columns:
    #             if tuple(index_list[-1]) not in self.Tau_Fit_Stats.index:
    #                 this_df = pa.DataFrame()
    #                 this_df.loc[:,'boot'] = pa.Series(flist,index=indicies)
    #                 this_df.loc[:,'Avg'] = pa.Series(aflist,index=indicies)
    #                 this_df.loc[:,'Std'] = pa.Series(sflist,index=indicies)
    #                 self.Tau_Fit_Stats = self.Tau_Fit_Stats.append(this_df)
    #         else:
    #             self.Tau_Fit_Stats.loc[:,'boot'] = pa.Series(flist,index=indicies)
    #             self.Tau_Fit_Stats.loc[:,'Avg'] = pa.Series(aflist,index=indicies)
    #             self.Tau_Fit_Stats.loc[:,'Std'] = pa.Series(sflist,index=indicies)
    #         self.Write()

    def ImportFitRanges(self,fit_range,tsum_fit_range):
        if not isinstance(fit_range,str):
            print('setsoffits implementation now only takes 1 fit range to scan over, choosing first')
            fit_range = fit_range[0]
        if not isinstance(tsum_fit_range,str):
            print('setsoffits implementation now only takes 1 fit range to scan over, choosing first')
            tsum_fit_range = tsum_fit_range[0]
        if 'PreDef' == fit_range:
            fit_range = self.PredefFits[0]
        if 'PreDef' == tsum_fit_range:
            tsum_fit_range = self.PredefFits[1]
        # self.fit_range = list(itertools.product(fit_range,tsum_fit_range))
        self.fit_range = [fit_range,tsum_fit_range]

    # def ImportFitRangesNP(self,fit_range,tsum_fit_range):
    #     if isinstance(fit_range,basestring):
    #         fit_range = [fit_range]
    #     if isinstance(tsum_fit_range,basestring):
    #         tsum_fit_range = [tsum_fit_range]
    #     if 'PreDef' == fit_range[0]:
    #         fit_range = [self.PredefFits[0]]
    #     if 'PreDef' == tsum_fit_range[0]:
    #         tsum_fit_range = [self.PredefFits[1]]
    #     self.fit_range = [[ifitr,itsumfitr] for ifitr,itsumfitr in zip(fit_range,tsum_fit_range)]



    def FitAlpha(   self,fit_range='PreDef',tsum_ranges='PreDef',iGuess='PreDef',
                    EstDir=False,thistflowlist='PreDef',WipeFit=True,show_timer=True):
        if self.is_NNQ2sub:
            print('fitting alpha not implemented (yet) for NNQ2Sub')
            return
        if self.Fun == 'None': raise IOError('Please define funciton using .SetFunction(Fun,npar) for fitting before calling .FitAlpha()')
        self.ImportFitRanges(fit_range,tsum_ranges)
        if iGuess != 'PreDef': self.AlphaiGuess = iGuess
        if thistflowlist == 'PreDef':
            if self.tflowfit == 'Not Set':
                raise IOError('tflowfit not set, please set by using .SetTflowFit()')
        else:
            if thistflowlist == 'All':
                self.tflowfit = self.tflowlist
            else:
                self.tflowfit = thistflowlist
        lFit,ilist = [],[]
        ifit,isumfit = self.fit_range
        if self.min_fit_len_tsum == 'Default':
            self.min_fit_len_tsum = unxmlfitr_int(isumfit.replace('tsumfitr','fitr'))
            self.min_fit_len_tsum = max(0,self.min_fit_len_tsum[1] - self.min_fit_len_tsum[0]-3)
        # llist = [ival[0][0] for ival in self.NNQFull_Stats['Alphaboot'].groupby(level=('momentum','flow_time'))]
        # if show_timer: thistimer = Timer( linklist=llist,name='Alpha Fits ')
        for (ip,itflow),pdata in self.NNQFull_Stats['Alphaboot'].groupby(level=('momentum','flow_time')):
            this_key = (ip,itflow)
            if ip != 'p000' or itflow not in self.tflowfit or ifit == 'None' or isumfit == 'None':
                # if show_timer: thistimer.Lap(flag=ip+' not being fit ')
                continue
            # ifit_str = '_'.join([ifit,isumfit])
            if (ip,itflow) in list(self.NNQFull_Fit_Stats.index):
                if WipeFit:
                    self.NNQFull_Fit_Stats.drop([(ip,itflow)],inplace=True)
                else:
                    # if show_timer: thistimer.Lap(flag=','.join([ip,itflow])+' already present')
                    continue

            # ifitmin,ifitmax = np.array(unxmlfitr(ifit),dtype=np.float64)
            # isumfitmin,isumfitmax = np.array(unxmlfitr(isumfit.replace('tsum','')),dtype=np.float64)

            # ##debugging
            # print
            # print 'FitData:'
            # for it,iy in zip(tdatarange[0],ydata):
            #     print it, iy
            # print
            fit_info = {}
            if EstDir:
                fit_info['Funs'] = (self.Fun,self.npar,'Estimate')
            else:
                fit_info['Funs'] = (self.Fun,self.npar)
            fit_info['iGuess'] = self.AlphaiGuess
            fit_info['name'] = self.name + ' Fits'
            # fit_info['paramlist'] = [r'\alpha']
            this_fit = sff.SetOfFitFuns(data=pdata.loc[ip,:,itflow,:])
            this_test = this_fit.ScanBox_fitr(ifit,isumfit,fit_info=fit_info,
                                            min_fit_len_1=0,min_fit_len_2=self.min_fit_len_tsum)
            if this_test:
                ilist.append(this_key)
                lFit.append(this_fit)
        if len(ilist) != 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.NNQFull_Fit_col_names)
            if 'boot' in self.NNQFull_Fit_Stats.columns:
                this_df = pa.DataFrame()
                this_df.loc[:,'boot'] = pa.Series(lFit,index=indicies)
                self.NNQFull_Fit_Stats = self.NNQFull_Fit_Stats.append(this_df)
            else:
                self.NNQFull_Fit_Stats.loc[:,'boot'] = pa.Series(lFit,index=indicies)
            def DoFit(val):
                val.DoFits(show_timer=show_timer)
                val.SortChi()
                return val
            self.NNQFull_Fit_Stats.loc[:,'boot'] = self.NNQFull_Fit_Stats['boot'].apply(DoFit)
            self.Write()

    ## implement setsoffits here
    def AlphaPlotVstiFitr(  self,plot_class,flowtime,thistfitr,fit_max,thismom='p000',
                            thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',
                            fit_min_list='All',Auto=False,WipeFit=False):
        print('Setsoffits have not implemented vstifitr yet, skipping')
        return plot_class
        # if self.is_NNQ2sub:
        #     print 'plotting anything other than PlotAlphaTsum not implemented for NNQ2Sub'
        #     return plot_class
        # self.CheckCol(thiscol)
        # self.CheckSym(thissym)
        # self.thisshiftscale = 'Not Set'
        # if isinstance(thistfitr,(list,tuple,np.ndarray)):
        #     thistfitr = 'fitr'+'-'.join(map(str,thistfitr))
        # if not isinstance(fit_max,basestring):
        #     fit_max = str(int(fit_max))
        # fitr_list = []
        # if fit_min_list != 'All':
        #     for imin in fit_min_list:
        #         if isinstance(imin,basestring):
        #             tsum_fitr = 'tsumfitr'+imin+'-'+fit_max
        #         else:
        #             tsum_fitr = 'tsumfitr'+str(imin)+'-'+fit_max
        #         fit_str = thistfitr+'_'+tsum_fitr
        #         if WipeFit or 'boot' not in self.NNQFull_Fit_Stats.columns or (thismom,flowtime,fit_str) not in self.NNQFull_Fit_Stats['boot'].index:
        #             print thismom,flowtime,fit_str , 'Not found, attempting to fit now'
        #             self.FitAlpha(fit_range=[thistfitr],thistflowlist = [flowtime],tsum_ranges=[tsum_fitr],WipeFit=WipeFit)
        #         fitr_list.append(fit_str)
        # else:
        #     if 'boot' not in self.NNQFull_Fit_Stats.columns or (thismom,flowtime) not in self.NNQFull_Fit_Stats['boot'].index:
        #         raise EnvironmentError('No fit ranges found for '+' '.join([thismom,flowtime]))
        #     else:
        #         alpha_hold = self.NNQFull_Fit_Stats.loc[(thismom,flowtime,slice(None)),'boot']
        #         alpha_hold = pa.Series(alpha_hold.values,[ival[2] for ival in list(alpha_hold.index)])
        #         fitr_list = list(alpha_hold[alpha_hold.index.str.contains(thistfitr+"_tsumfitr.*-"+fit_max)].index)
        #
        # thisdata = self.NNQFull_Fit_Stats.loc[(thismom,flowtime,fitr_list),'boot']
        # thisdata = thisdata.apply(lambda x : x.fit_data['Params'].iloc[0])
        # xlist = list(thisdata.index)
        # if isinstance(xlist[0],tuple):
        #     xlist = [ix[-1] for ix in xlist]
        # plotx = map(lambda x : x.replace(thistfitr+'_tsumfitr',''),xlist)
        # plotx = map(lambda x : int(x.replace('-'+fit_max,'')),plotx)
        # ploty,plotyerr = zip(*[(ival.Avg,ival.Std) for ival in thisdata.values])
        # plotshift = self.thisshift*np.abs(plotx[-1]-plotx[0])
        # thisleglab = self.LegLab +'$\ '+ tflowTOLeg(flowtime) + '\ t_{e}='+fit_max+'$'
        # if thismom != 'p000': thisleglab = thisleglab+'\ '+thismom.replace('p','p=')
        # # if Auto: thisleglab  = thisleglab + ' Auto'
        # hold_series = pa.Series()
        # hold_series['x_data'] = plotx
        # hold_series['y_data'] = ploty
        # hold_series['yerr_data'] = plotyerr
        # hold_series['xerr_data'] = None
        # hold_series['type'] = 'error_bar'
        # hold_series['fit_class'] = None
        # hold_series['label'] = thisleglab
        # hold_series['symbol'] = self.thissym
        # hold_series['color'] = self.thiscol
        # hold_series['shift'] = plotshift
        # hold_series['otherXvals'] = None
        # hold_series['xaxis'] = None
        # hold_series['xdatarange'] = None
        # hold_series['ShowPar'] = r'\alpha'
        # hold_series['Phys'] = None
        # plot_class.AppendData(hold_series)
        # initial_time = tstr(thistfitr.replace('fitr','').split('-')[0])
        # plot_class = self.PlotIntLine(plot_class,thismom,flowtime,initial_time,self.thiscol)
        # return plot_class


    def AlphaPlotVstfFitr(  self,plot_class,flowtime,thistfitr,fit_min,thismom='p000',
                            thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',
                            fit_max_list='All',WipeFit=False):
        print('Setsoffits have not implemented vstffitr yet, skipping')
        return plot_class
        # if self.is_NNQ2sub:
        #     print 'plotting anything other than PlotAlphaTsum not implemented for NNQ2Sub'
        #     return plot_class
        # self.CheckCol(thiscol)
        # self.CheckSym(thissym)
        # self.thisshiftscale = 'Not Set'
        # if isinstance(thistfitr,(list,tuple,np.ndarray)):
        #     thistfitr = 'fitr'+'-'.join(map(str,thistfitr))
        # if not isinstance(fit_min,basestring):
        #     fit_min = str(int(fit_min))
        # fitr_list = []
        # if fit_max_list != 'All':
        #     for imax in fit_max_list:
        #         if isinstance(imax,basestring):
        #             tsum_fitr = 'tsumfitr'+fit_min+'-'+imax
        #         else:
        #             tsum_fitr = 'tsumfitr'+fit_min+'-'+str(imax)
        #         fit_str = thistfitr+'_'+tsum_fitr
        #         if WipeFit or 'boot' not in self.NNQFull_Fit_Stats.columns or (thismom,flowtime,fit_str) not in self.NNQFull_Fit_Stats['boot'].index:
        #             print thismom,flowtime,fit_str , 'Not found, attempting to fit now'
        #             self.FitAlpha(fit_range=[thistfitr],thistflowlist = [flowtime],tsum_ranges=[tsum_fitr],WipeFit=WipeFit)
        #         fitr_list.append(fit_str)
        # else:
        #     if 'boot' not in self.NNQFull_Fit_Stats.columns or (thismom,flowtime) not in self.NNQFull_Fit_Stats['boot'].index:
        #         raise EnvironmentError('No fit ranges found for '+' '.join([thismom,flowtime]))
        #     else:
        #         alpha_hold = self.NNQFull_Fit_Stats.loc[(thismom,flowtime,slice(None)),'boot']
        #         alpha_hold = pa.Series(alpha_hold.values,[ival[2] for ival in list(alpha_hold.index)])
        #         fitr_list = list(alpha_hold[alpha_hold.index.str.contains(thistfitr+'_tsumfitr'+fit_min+'-')].index)
        #
        # thisdata = self.NNQFull_Fit_Stats.loc[(thismom,flowtime,fitr_list),'boot']
        # # thisdata = thisdata[fitr_list]
        # thisdata = thisdata.apply(lambda x : x.fit_data['Params'].iloc[0])
        #
        # xlist = list(thisdata.index)
        # if isinstance(xlist[0],tuple):
        #     xlist = [ix[-1] for ix in xlist]
        # plotx = map(lambda x : int(x.replace(thistfitr+'_tsumfitr'+fit_min+'-','')),xlist)
        # ploty,plotyerr = zip(*[(ival.Avg,ival.Std) for ival in thisdata.values])
        # plotshift = self.thisshift*np.abs(plotx[-1]-plotx[0])
        # thisleglab = self.LegLab +'$\ '+ tflowTOLeg(flowtime) + '\ t_{i}='+fit_min+'$'
        # if thismom != 'p000': thisleglab = thisleglab+'\ '+thismom.replace('p','p=')
        # # if Auto: thisleglab  = thisleglab + ' Auto'
        # hold_series = pa.Series()
        # hold_series['x_data'] = plotx
        # hold_series['y_data'] = ploty
        # hold_series['yerr_data'] = plotyerr
        # hold_series['xerr_data'] = None
        # hold_series['type'] = 'error_bar'
        # hold_series['fit_class'] = None
        # hold_series['label'] = thisleglab
        # hold_series['symbol'] = self.thissym
        # hold_series['color'] = self.thiscol
        # hold_series['shift'] = plotshift
        # hold_series['otherXvals'] = None
        # hold_series['xaxis'] = None
        # hold_series['xdatarange'] = None
        # hold_series['ShowPar'] = r'\alpha'
        # hold_series['Phys'] = None
        # plot_class.AppendData(hold_series)
        # initial_time = tstr(thistfitr.replace('fitr','').split('-')[0])
        # plot_class = self.PlotIntLine(plot_class,thismom,flowtime,initial_time,self.thiscol)
        # return plot_class



    ## make sure all xlims are the same, because this sets the window of xlim (needed for calculating shifts)
    ## thistflowlist is array of tflows to plot over (currently supports length 1 arrays).
    ## ToDo: make color and symbol and shift iterate over this. Currently handeled externally (SetsOfCorrs.py) with repeated calls
    def AlphaPlot(  self,plot_class,xlims=defxlimAlpha,momlist=['p000'],
                    thistflowlist = 'PreDef',tsum_list='PreDef',spec_tsum_list='PreDef',
                    thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',Auto=False):
        if self.is_NNQ2sub:
            print('plotting anything other than PlotAlphaTsum not implemented for NNQ2Sub')
            return plot_class
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        self.thisshiftscale = 'Not Set'
        df_col_name = 'Alphaboot'
        if Auto:
            df_col_name = 'AlphaAuto'
        if df_col_name not in self.NNQInt_Stats or df_col_name not in self.NNQFull_Stats:
            print(df_col_name,'not computed for ',self.name)
            return plot_class
        # if isinstance(xlims[0],basestring): xlims = map(untstr,xlims)
        thisshift = self.GetShift(xlims,thisshift)
        if thistflowlist == 'PreDef':
            thistflowlist = self.tflowfit
        elif not isinstance(thistflowlist,list):
            thistflowlist = [thistflowlist]
        if tsum_list == 'PreDef':
            tsum_list = self.tsum_fit_list
        if spec_tsum_list == 'PreDef':
            spec_tsum_list = self.tsum_fit_list
        elif not isinstance(tsum_list,list):
            tsum_list = [tsum_list]
        if len(momlist) == 0: momlist = self.NNQFull_Stats.index.get_level_values('momentum')
        ip,it_sum,iflow,jt_sum = momlist[0],tsum_list[0],thistflowlist[0],spec_tsum_list[0]
        # for ip,it_sum,iflow in itertools.product(momlist,tsum_list,thistflowlist):
        if 'All' in ip: ip = 'p000'
        if 'Int' in it_sum:
            thisdata = self.NNQInt_Stats[df_col_name]
            # pdata = thisdata.swaplevel(1,2)[ip,iflow]
            this_key = (ip,slice(None),iflow)
            if (ip,iflow) not in thisdata.swaplevel(1,2):
                print(ip,iflow,' not found in Op Integrated Alpha keys')
                return plot_class
        else:
            if self.is_NNQ2sub:
                thisdata = self.NNQ_subNNQ2_Stats[df_col_name]
                # pdata = thisdata.swaplevel(1,3)[ip,it_sum,iflow]
                this_key = (ip,slice(None),iflow,it_sum,jt_sum)
                if (ip,jt_sum,iflow,it_sum) not in thisdata.swaplevel(1,4):
                    print(ip, iflow,it_sum,jt_sum,' not found in Alpha keys')
                    return plot_class
            else:
                thisdata = self.NNQFull_Stats[df_col_name]
                # pdata = thisdata.swaplevel(1,3)[ip,it_sum,iflow]
                this_key = (ip,slice(None),iflow,it_sum)
                if (ip,it_sum,iflow) not in thisdata.swaplevel(1,3):
                    print(ip, it_sum,iflow,' not found in Alpha keys')
                    return plot_class
        # if len(momlist) > 1: thisleglab = thisleglab+'\ '+ip.replace('p','p=')
        thisleglab = self.LegLab
        if Auto: thisleglab  = thisleglab + ' Auto'
        dataavg = thisdata.apply(lambda x : x.Avg)
        dataerr = thisdata.apply(lambda x : x.Std)
        # dataavg,dataerr = zip(*[(idata.Avg,idata.Std) for idata in databoot])
        # tlist,databoot = pdata.keys()[xlims[0]:xlims[1]],pdata.values[xlims[0]:xlims[1]]
        hold_series = pa.Series()
        # hold_series['x_data'] = np.array(map(untstr,tlist))
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = dataavg
        hold_series['yerr_data'] = dataerr
        hold_series['key_select'] = this_key
        hold_series['xerr_data'] = None
        hold_series['type'] = 'error_bar_vary'
        hold_series['fit_class'] = None
        hold_series['label'] = thisleglab
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        hold_series['otherXvals'] = None
        hold_series['xaxis'] = None
        hold_series['xdatarange'] = None
        hold_series['ShowPar'] = r'\alpha'
        hold_series['Phys'] = None
        hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
        plot_class.AppendData(hold_series)
        return plot_class

    def AlphaPlotImproveFit(  self,plot_class,fitr,xlims=defxlimAlpha,momlist=['p000'],
                            thistflowlist = 'PreDef',WipeFit=False,
                            thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',
                            this_par='First'):
        if self.is_NNQ2sub:
            print('plotting anything other than PlotAlphaTsum not implemented for NNQ2Sub')
            return plot_class
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        self.thisshiftscale = 'Not Set'
        # if isinstance(xlims[0],basestring): xlims = map(untstr,xlims)
        thisshift = self.GetShift(xlims,thisshift)
        if thistflowlist == 'PreDef':
            thistflowlist = self.tflowfit
        elif not isinstance(thistflowlist,list):
            thistflowlist = [thistflowlist]
        if len(momlist) == 0: momlist = self.NNQFull_Stats.index.get_level_values('momentum')
        ip,iflow = momlist[0],thistflowlist[0]
        # for ip,iflow in itertools.product(momlist,thistflowlist):
        # thisleglab = self.LegLab
        # if len(momlist) > 1: thisleglab = thisleglab+'\ '+ip.replace('p','p=')
        fit_scan = 'fitr'+str(xlims[0])+'-'+str(xlims[1]+1)
        fitr_list = ['fitr'+str(ix)+'-'+str(ix+1) for ix in range(xlims[0],xlims[1]+1)]
        if 'tsumfitr' not in fitr:
            fitr = fitr.replace('fitr','tsumfitr')

        if WipeFit or 'boot' not in self.NNQFull_Fit_Stats:
            self.FitAlpha(fit_range=fit_scan,thistflowlist=[iflow],tsum_ranges=fitr,WipeFit=WipeFit)
        fit_data = self.Get_Extrapolation()
        this_fun = fit_data.index[0][3]
        if this_par == 'First':
            this_par = fit_data.index[0][6]
        elif isinstance(this_par,int):
            this_par = fit_data.index.levels[6][this_par]
        this_key_list = [(ip,iflow,this_fun,ifitr,fitr,this_par) for ifitr in fitr_list]
        this_key_select = (ip,iflow,this_fun,slice(None),fitr,this_par)
        if WipeFit or any([ikey not in fit_data.index for ikey in this_key_list]):
            self.FitAlpha(fit_range=fit_scan,thistflowlist=[iflow],tsum_ranges=fitr,WipeFit=WipeFit)
            fit_data = self.Get_Extrapolation()
        fit_data = fit_data.loc[(slice(None),slice(None),slice(None),fitr_list,slice(None),slice(None))]


        hold_series = null_series
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = fit_data.apply(lambda x : x.Avg)
        hold_series['yerr_data'] = fit_data.apply(lambda x : x.Std)
        hold_series['key_select'] = this_key_select
        hold_series['type'] = 'error_bar_vary'
        hold_series['label'] = self.LegLab
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        hold_series['ShowPar'] = r'\alpha'
        hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
        plot_class.AppendData(hold_series)
        return plot_class

        # plot_class.get_xlim(*xlims)



    def AlphaPlotTflow( self,plot_class,tflowplot = defxlimOp,
                        momlist=['p000'],thistsinklist = 'PreDef',tsum_list='PreDef',
                        thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',Auto=False):
        # def WithinBounds(tlist):
        #     output,indicies = [],[]
        #     for ic,ilist in enumerate(tlist):
        #         if untflowstr(ilist) in tflowplot:
        #             output.append(ilist)
        #             indicies.append(ic)
        #     return output,indicies
        if self.is_NNQ2sub:
            print('plotting anything other than PlotAlphaTsum not implemented for NNQ2Sub')
            return plot_class
        df_col_name = 'Alphaboot'
        if Auto:
            df_col_name = 'AlphaAuto'
            if df_col_name not in self.NNQInt_Stats or df_col_name not in self.NNQFull_Stats:
                print('Auto not computed for ',self.name)
                return plot_class
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        tflowplot = list(map(untflowstr,tflowplot))
        if tsum_list == 'PreDef': tsum_list = self.tsum_fit_list
        if thistsinklist == 'PreDef':
            thistsinklist = self.tsinkfit
        elif not isinstance(thistsinklist,list):
            thistsinklist = [thistsinklist]
        if tsum_list == 'PreDef':
            tsum_list = self.tsum_fit_list
        elif not isinstance(tsum_list,list):
            tsum_list = [tsum_list]
        if len(momlist) == 0: momlist = self.NNQFull_Stats.index.get_level_values('momentum')
        ip, it_sum,itsink = momlist[0],tsum_list[0],thistsinklist[0]
        # for ip,it_sum,itsink in itertools.product(momlist,tsum_list,thistsinklist):
        if 'All' in ip: ip = 'p000'
        if 'Int' in it_sum:
            thisdata = self.NNQInt_Stats[df_col_name]
            this_key = (ip,itsink,slice(None))
            if (ip,itsink) not in thisdata:
                print(ip,itsink,' not found in Op Integrated Alpha keys')
                return plot_class
            tflowdata = thisdata[ip,itsink]
            thisleglab = self.LegLab
        else:
            thisdata = self.NNQFull_Stats[df_col_name]
            this_key = (ip,itsink,slice(None),it_sum)
            if (ip,itsink,it_sum) not in thisdata.swaplevel(2,3):
                print(ip, itsink,it_sum,' not found in Alpha keys')
                return plot_class
            tflowdata = thisdata.swaplevel(2,3)[ip,itsink,it_sum]
            thisleglab = self.LegLab
        if Auto: thisleglab = thisleglab + ' Auto'
        tflowlist = list(tflowdata.keys())
        # tflowlist,indicies = WithinBounds(np.array(tflowlist))
        if len(tflowlist) > 0 and isinstance(tflowlist[0],(list,tuple,np.ndarray)):
            tflowlist = [itf[2] for itf in tflowlist]
        tflowlist = TflowToPhys(np.array(tflowlist),self.NN.latparams.latspace)
        plotshift = self.GetShift([tflowlist[0],tflowlist[-1]],thisshift)
        # print
        # print thisshift
        # print plotshift
        # print [tflowlist[0],tflowlist[-1]]
        # print np.array(tflowlist)+plotshift
        # databoot = np.array(databoot)[indicies]
        dataavg = thisdata.apply(lambda x : x.Avg)
        dataerr = thisdata.apply(lambda x : x.Std)
        # dataavg,dataerr = zip(*[(idata.Avg,idata.Std) for idata in databoot])
        hold_series = pa.Series()
        hold_series['x_data'] = np.array(tflowlist)
        hold_series['y_data'] = dataavg
        hold_series['yerr_data'] = dataerr
        hold_series['xerr_data'] = None
        hold_series['type'] = 'error_bar_vary'
        hold_series['key_select'] = this_key
        hold_series['fit_class'] = None
        hold_series['label'] = thisleglab
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['otherXvals'] = None
        hold_series['xaxis'] = None
        hold_series['shift'] = plotshift
        hold_series['xdatarange'] = None
        hold_series['ShowPar'] = r'\alpha'
        hold_series['Phys'] = None
        hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
        plot_class.AppendData(hold_series)
        return plot_class

    def PlotIntTsink(   self,plot_class,this_mom,this_flow,xlims=defxlimAlpha,Auto=False,
                        thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine'):
        if self.is_NNQ2sub:
            print('plotting anything other than PlotAlphaTsum not implemented for NNQ2Sub')
            return plot_class
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        if isinstance(xlims[0],str): xlims = list(map(untstr,xlims))
        # xlim_str = [tstr(ix+1) for ix in xlims]
        thisshift = self.GetShift(xlims,thisshift)
        df_col_name = 'Alphaboot'
        if Auto:
            df_col_name = 'AlphaAuto'
            if df_col_name not in self.NNQInt_Stats or df_col_name not in self.NNQFull_Stats:
                print('Auto not computed for ',self.name)
                return plot_class

        if isinstance(this_mom,(list,tuple,np.ndarray)):
            this_mom = this_mom[0]
        if (this_mom,this_flow) in self.NNQInt_Stats[df_col_name].swaplevel(1,2).index:
            plot_data = self.NNQInt_Stats.loc[:,df_col_name]
            # print 'DEBUG'
            # print plot_data
            # print list(plot_data.keys())
            # # print plot_data.apply(lambda x : x.Avg)
            # print xlim_str
            # print plot_data.loc[xlim_str[0]:xlim_str[1]]
            # plot_data = plot_data.loc[xlim_str[0]:xlim_str[1]]
            # xdata = map(untstr,list(plot_data.keys()))
            mid = plot_data.apply(lambda x : x.Avg)
            err = plot_data.apply(lambda x : x.Std)
            # up,down = mid+err,mid-err
            thisleglab = self.LegLab + r' $\alpha_{Int}$'
            if Auto: thisleglab += ' Auto'
            hold_series = pa.Series()
            hold_series['x_data'] = 'from_keys'
            hold_series['y_data'] = mid
            hold_series['yerr_data'] = err
            hold_series['key_select'] = (this_mom,slice(None),this_flow)
            hold_series['xerr_data'] = None
            hold_series['type'] = 'error_bar_vary'
            hold_series['fit_class'] = None
            hold_series['label'] = thisleglab
            hold_series['symbol'] = self.thissym
            hold_series['color'] = self.thiscol
            hold_series['shift'] = thisshift
            hold_series['otherXvals'] = None
            hold_series['xaxis'] = None
            hold_series['xdatarange'] = None
            hold_series['ShowPar'] = r'\alpha'
            hold_series['Phys'] = None
            hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
            plot_class.AppendData(hold_series)
            # pl.axhline(mid,color=thiscol,label = )
            # pl.axhspan(up,down, alpha=fillalpha, color=)
        else:
            print((this_mom,this_flow) , 'has no NNQ Integrated result')
        return plot_class


    def PlotIntLine(self,plot_class,this_mom,this_flow,this_tsink,thiscol,Auto=False):
        df_col_name = 'Alphaboot'
        if Auto:
            df_col_name = 'AlphaAuto'
        if df_col_name not in self.NNQInt_Stats:
            print('Auto not computed for ',self.name)
            return plot_class
        if (this_mom,this_tsink,this_flow) in list(self.NNQInt_Stats[df_col_name].keys()):
            mid = self.NNQInt_Stats[df_col_name].apply(lambda x : x.Avg)
            err = self.NNQInt_Stats[df_col_name].apply(lambda x : x.Std)
            # up,down = mid+err,mid-err
            linelab = r'ACTUAL_VALUE_\alpha_{Int}'
            if Auto: linelab = 'Auto ' + linelab
            hold_series = pa.Series()
            hold_series['x_data'] = None
            hold_series['y_data'] = mid
            hold_series['yerr_data'] = err
            hold_series['xerr_data'] = None
            hold_series['key_select'] = (this_mom,this_tsink,this_flow)
            hold_series['type'] = 'hline_vary'
            hold_series['fit_class'] = None
            hold_series['label'] = linelab
            hold_series['symbol'] = self.thissym
            hold_series['color'] = thiscol
            hold_series['shift'] = 0.0
            hold_series['otherXvals'] = None
            hold_series['xaxis'] = None
            hold_series['xdatarange'] = None
            hold_series['ShowPar'] = r'\alpha'
            hold_series['Phys'] = None
            hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
            plot_class.AppendData(hold_series)
            # pl.axhline(mid,color=thiscol,label = )
            # pl.axhspan(up,down, alpha=fillalpha, color=)
        else:
            print((this_mom,this_tsink,this_flow) , 'has no NNQ Integrated result')
        return plot_class



    def Plot_Tau_Fit(   self,plot_class,fitr,mom,tsink,
                        tflow,thiscol,thisshift,xlims='Data',WipeFit=False):
        if self.is_NNQ2sub:
            print('PlotAlpha or PlotAlphaTsum are only done for NNQ2Sub')
            return plot_class
        if isinstance(fitr,(list,tuple,np.ndarray)):
            fitr = '_'.join(fitr)
        fitr_list = fitr.split('_')
        if WipeFit or 'boot' not in self.NNQFull_Fit_Stats.columns :
            print(mom,tflow , 'Not found, attempting to fit now')
            self.FitAlpha(fit_range=fitr_list[0],thistflowlist = [tflow],tsum_ranges=fitr_list[1],WipeFit=WipeFit)
        fit_data = sff.PullSOFSeries(self.NNQFull_Fit_Stats['boot'])
        fit_data = Series_fix_key(fit_data,4,'fittwor','tsumfitr')
        this_fun = fit_data.index[0][2]
        this_key = (mom,tflow,this_fun,fitr_list[0],fitr_list[1])
        if this_key not in fit_data:
            print(this_key , 'Not found, attempting to fit now')
            self.FitAlpha(fit_range=fitr_list[0],thistflowlist = [tflow],tsum_ranges=fitr_list[1],WipeFit=WipeFit)
        fit_data = sff.PullSOFSeries(self.NNQFull_Fit_Stats['boot'])
        fit_data = Series_fix_key(fit_data,4,'fittwor','tsumfitr')
        this_fun = fit_data.index[0][2]
        this_key = (mom,tflow,this_fun,fitr_list[0],fitr_list[1])
        if this_key not in fit_data:
            print(list(fit_data.index))
            print(this_key)
            raise IOError('error with fitting')
        # if isinstance(fit_plot,pa.Series):
        #     fit_plot = fit_plot.iloc[0]
        hold_series = null_series
        hold_series['key_select'] = this_key
        hold_series['type'] = 'fit_vary'
        hold_series['fit_class'] = fit_data
        hold_series['label'] = 'Fit '
        hold_series['color'] = thiscol
        hold_series['shift'] = thisshift
        hold_series['xdatarange'] = xlims
        hold_series['otherXvals'] = [untstr(tsink)]
        hold_series['xaxis'] = 1
        hold_series['ShowPar'] = r'\alpha'
        hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
        plot_class.AppendData(hold_series)
        return plot_class


    def Plot_Tsink_Fit( self,plot_class,fitr,mom,tsum,
                        tflow,thiscol,thisshift,xlims='Data',WipeFit=False):
        if self.is_NNQ2sub:
            print('PlotAlpha or PlotAlphaTsum are only done for NNQ2Sub')
            return plot_class
        if isinstance(fitr,(list,tuple,np.ndarray)):
            fitr = '_'.join(fitr)
        fitr_list = fitr.split('_')
        if isinstance(mom,(list,tuple,np.ndarray)):
            mom = mom[0]
        if WipeFit or 'boot' not in self.NNQFull_Fit_Stats.columns :
            print(mom,tflow , 'Not found, attempting to fit now')
            self.FitAlpha(fit_range=fitr_list[0],thistflowlist = [tflow],tsum_ranges=fitr_list[1],WipeFit=WipeFit)
        fit_data = sff.PullSOFSeries(self.NNQFull_Fit_Stats['boot'])
        fit_data = Series_fix_key(fit_data,4,'fittwor','tsumfitr')

        this_fun = fit_data.index[0][2]
        this_key = (mom,tflow,this_fun,fitr_list[0],fitr_list[1])
        if this_key not in fit_data:
            print(this_key , 'Not found, attempting to fit now')
            self.FitAlpha(fit_range=fitr_list[0],thistflowlist = [tflow],tsum_ranges=fitr_list[1],WipeFit=WipeFit)
        fit_data = sff.PullSOFSeries(self.NNQFull_Fit_Stats['boot'])
        fit_data = Series_fix_key(fit_data,4,'fittwor','tsumfitr')

        this_fun = fit_data.index[0][2]
        this_key = (mom,tflow,this_fun,fitr_list[0],fitr_list[1])
        if this_key not in fit_data:
            print(list(fit_data.index))
            print(this_key)
            raise IOError('error with fitting')

        # if isinstance(fit_plot,pa.Series):
        #     fit_plot = fit_plot.iloc[0]
        hold_series = null_series
        hold_series['key_select'] = this_key
        hold_series['type'] = 'fit_vary'
        hold_series['fit_class'] = fit_data
        hold_series['label'] = 'Fit '
        hold_series['color'] = thiscol
        hold_series['shift'] = thisshift
        hold_series['xdatarange'] = xlims
        hold_series['otherXvals'] = [untstr(tsum)]
        hold_series['xaxis'] = 0
        hold_series['ShowPar'] = r'\alpha'
        hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
        plot_class.AppendData(hold_series)
        return plot_class



    def AlphaPlotTsum( self,plot_class,momlist=['p000'],thistsinklist = 'PreDef',
                        thistflowlist='PreDef',fit_list='None',WipeFit=False,spec_tsum_list='PreDef',
                        thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',Auto=False):
        if fit_list == 'None':
            fit_list = ['None']
        if isinstance(spec_tsum_list,str):
            spec_tsum_list = [spec_tsum_list]
        if isinstance(fit_list,str):
            fit_list = fit_list.split('_')
        if not isinstance(fit_list[0],(list,tuple,np.ndarray)) and fit_list[0] != 'None':
            fit_list = [[fit_list[0],fit_list[1]]]
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        df_col_name = 'Alphaboot'
        if Auto:
            df_col_name = 'AlphaAuto'
        this_logic = df_col_name not in self.NNQInt_Stats
        if self.is_NNQ2sub:
            this_logic = this_logic or df_col_name not in self.NNQ_subNNQ2_Stats
        else:
            this_logic = this_logic or df_col_name not in self.NNQFull_Stats
        if this_logic:
            print(df_col_name,' not computed for ',self.name)
            if self.is_NNQ2sub:
                print(self.NNQ_subNNQ2_Stats.columns)
            else:
                print(self.NNQFull_Stats.columns)
            return plot_class
        if thistsinklist == 'PreDef':
            thistsinklist = self.tsinkfit
        elif not isinstance(thistsinklist,list):
            thistsinklist = [thistsinklist]
        if thistflowlist == 'PreDef':
            thistflowlist = self.tsum_fit_list
        elif not isinstance(thistflowlist,list):
            thistflowlist = [thistflowlist]


        if len(momlist) == 0: momlist = self.NNQFull_Stats.index.get_level_values('momentum')
        # for fitr,ip,itflow,itsink in itertools.product(fit_list,momlist,thistflowlist,thistsinklist):
        fitr,ip,itflow,itsink,jt_sum  = fit_list[0],momlist[0],thistflowlist[0],thistsinklist[0],spec_tsum_list[0]
        if self.is_NNQ2sub:
            thisdata = self.NNQ_subNNQ2_Stats[df_col_name]
            # pdata = thisdata.swaplevel(1,3)[ip,it_sum,iflow]
            this_key = (ip,itsink,itflow,slice(None),jt_sum)
            if (ip,itsink,itflow,jt_sum) not in thisdata.swaplevel(3,4):
                print(ip,itsink,itflow,jt_sum,' not found in Alpha keys')
                return plot_class
            # tflowdata = thisdata[ip,itsink,itflow,slice(None),jt_sum]
        else:
            this_key = (ip,itsink,itflow,slice(None))
            thisdata = self.NNQFull_Stats[df_col_name]
            if (ip,itsink,itflow) not in thisdata:
                print(ip,itsink, itsink ,' not found in Alphaboot keys')
                return plot_class
            # tflowdata = thisdata[ip,itsink,itflow]
        thisleglab = self.LegLab
        # if len(momlist) > 1: thisleglab = thisleglab+'\ '+ip.replace('p','p=')
        if Auto: thisleglab += ' Auto'
        # this_slice = slice(self.boot_tsum_list[0],self.boot_tsum_list[-1])

        # tsum_plot = map(untstr,tflowdata.keys()[this_slice])
        # if len(tsum_plot) == 0:
        #     plotshift = thisshift
        # else:
        #     plotshift = self.GetShift([tsum_plot[0],tsum_plot[-1]],thisshift)
        # print
        # print thisshift
        # print plotshift
        # print [tsum_plot[0],tsum_plot[-1]]
        # print np.array(tsum_plot)+plotshift
        dataavg = thisdata.apply(lambda x : x.Avg)
        dataerr = thisdata.apply(lambda x : x.Std)
        # databoot = np.array(databoot)
        # dataavg,dataerr = zip(*[(idata.Avg,idata.Std) for idata in databoot])
        hold_series = pa.Series()
        hold_series['x_data'] = 'from_keys'
        hold_series['key_select'] = this_key
        hold_series['y_data'] = dataavg
        hold_series['yerr_data'] = dataerr
        hold_series['xerr_data'] = None
        hold_series['type'] = 'error_bar_vary'
        hold_series['fit_class'] = None
        hold_series['label'] = thisleglab
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['otherXvals'] = None
        hold_series['xaxis'] = None
        hold_series['shift'] = thisshift
        hold_series['xdatarange'] = None
        hold_series['ShowPar'] = r'\alpha'
        hold_series['Phys'] = None
        hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
        plot_class.AppendData(hold_series)
        plot_class = self.PlotIntLine(plot_class,ip,itflow,itsink,self.thiscol,Auto=Auto)
        if 'None' not in fitr and not Auto:
            plot_class = self.Plot_Tau_Fit(plot_class,fitr,ip,itsink,itflow,self.thiscol,thisshift,WipeFit=WipeFit)
        return plot_class
    # plot_class.get_xlim(*xlims)

    def AlphaPlotCovar( self,plot_class,momlist=['p000'],thistsinklist = 'PreDef',
                        thistflowlist='PreDef',spec_tsum_list='PreDef',
                        thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine'):
        if isinstance(spec_tsum_list,str):
            spec_tsum_list = [spec_tsum_list]
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        df_col_name = 'AlphaAuto'
        if self.is_NNQ2sub:
            this_logic = df_col_name not in self.NNQ_subNNQ2_Stats
        else:
            this_logic = df_col_name not in self.NNQFull_Stats
        if this_logic:
            print(df_col_name,' not computed for ',self.name)
            if self.is_NNQ2sub:
                print(self.NNQ_subNNQ2_Stats.columns)
            else:
                print(self.NNQFull_Stats.columns)
            return plot_class
        if thistsinklist == 'PreDef':
            thistsinklist = self.tsinkfit
        elif not isinstance(thistsinklist,list):
            thistsinklist = [thistsinklist]
        if thistflowlist == 'PreDef':
            thistflowlist = self.tsum_fit_list
        elif not isinstance(thistflowlist,list):
            thistflowlist = [thistflowlist]


        if len(momlist) == 0: momlist = self.NNQFull_Stats.index.get_level_values('momentum')
        # for fitr,ip,itflow,itsink in itertools.product(fit_list,momlist,thistflowlist,thistsinklist):
        ip,itflow,itsink,jt_sum  = momlist[0],thistflowlist[0],thistsinklist[0],spec_tsum_list[0]
        if self.is_NNQ2sub:
            thisdata = self.NNQ_subNNQ2_Stats[df_col_name]
            # pdata = thisdata.swaplevel(1,3)[ip,it_sum,iflow]
            this_key = (ip,itsink,itflow,slice(None),jt_sum,'C12')
            if (ip,itsink,itflow,jt_sum) not in thisdata.swaplevel(3,4):
                print(ip,itsink,itflow,jt_sum,' not found in Alpha keys')
                return plot_class
            # tflowdata = thisdata[ip,itsink,itflow,slice(None),jt_sum]
            covarl,indexl = [],[]
            for ikey,idata in thisdata.items():
                for ic,icovar in enumerate(idata.covar):
                    for jc,jcovar in enumerate(icovar):
                        if jc >= ic:
                            indexl.append(ikey+('C'+str(ic+1)+str(jc+1),))
                            covarl.append(jcovar)
            if len(indexl) > 0:
                indicies = pa.MultiIndex.from_tuples(indexl,names=list(self.NNQ_subNNQ2_Stats.index.names) + ['Covar_matrix'])
                ploty = pa.Series(covarl,index=indicies)
            else:
                ploty = pa.Series()
        else:
            this_key = (ip,itsink,itflow,slice(None),'C12')
            thisdata = self.NNQFull_Stats[df_col_name]
            if (ip,itsink,itflow) not in thisdata:
                print(ip,itsink, itsink ,' not found in Alphaboot keys')
                return plot_class
            covarl,indexl = [],[]
            # try:
            #     tflowdata = thisdata.loc[(ip,itsink,itflow,slice(None))]
            # except Exception as err:
            #     out_str = 'not in tflowdata:\n' + ', '.join([ip,itsink,itflow,'None'])
            #     out_str += '\ntotal_list:\n'+'\n'.join([', '.join(ival) for ival in list(thisdata.index)])
            #     raise EnvironmentError(str(err)+'\n'+out_str)
            for ikey,idata in thisdata.items():
                for ic,icovar in enumerate(idata.covar):
                    for jc,jcovar in enumerate(icovar):
                        if jc >= ic:
                            indexl.append(ikey+('C'+str(ic+1)+str(jc+1),))
                            covarl.append(jcovar)
            if len(indexl) > 0:
                indicies = pa.MultiIndex.from_tuples(indexl,names=list(self.NNQFull_Stats.index.names) + ['Covar_matrix'])
                ploty = pa.Series(covarl,index=indicies)
            else:
                ploty = pa.Series()

        thisleglab = self.LegLab + ' Covar'
        # if len(momlist) > 1: thisleglab = thisleglab+'\ '+ip.replace('p','p=')
        # this_slice = slice(self.boot_tsum_list[0],self.boot_tsum_list[-1])

        # tsum_plot = map(untstr,tflowdata.keys()[this_slice])
        # if len(tsum_plot) == 0:
        #     plotshift = thisshift
        # else:
        #     plotshift = self.GetShift([tsum_plot[0],tsum_plot[-1]],thisshift)
        # print
        # print thisshift
        # print plotshift
        # print [tsum_plot[0],tsum_plot[-1]]
        # print np.array(tsum_plot)+plotshift
        # databoot = np.array(databoot)
        # dataavg,dataerr = zip(*[(idata.Avg,idata.Std) for idata in databoot])
        hold_series = pa.Series()
        hold_series['x_data'] = 'from_keys'
        hold_series['key_select'] = this_key
        hold_series['y_data'] = ploty
        hold_series['yerr_data'] = None
        hold_series['xerr_data'] = None
        hold_series['type'] = 'scatter_vary'
        hold_series['fit_class'] = None
        hold_series['label'] = thisleglab
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['otherXvals'] = None
        hold_series['xaxis'] = None
        hold_series['shift'] = thisshift
        hold_series['xdatarange'] = None
        hold_series['ShowPar'] = r'\alpha'
        hold_series['Phys'] = None
        hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
        plot_class.AppendData(hold_series)
        return plot_class
    # plot_class.get_xlim(*xlims)


    #state is 1 or 2 state fit
    ## fitr is fit range e.g. fitr5-10
    ## fit range of 'Data' is just the original fit range
    def AlphaFitPlot(   self,plot_class,fitr,xlims = 'Data',momlist=['p000'],
                        thistflowlist = 'PreDef',thispar='First',WipeFit=False,
                        thiscol='PreDefine',thisshift='PreDefine',Auto=False):
        # autocorrelation has not fits (currenlty)
        if Auto: return plot_class
        if fitr == 'None': return plot_class
        self.CheckCol(thiscol)
        plotshift = self.GetShift(plot_class.get_xlim(),thisshift)
        if isinstance(fitr,(tuple,list,np.ndarray)):
            fitr = '_'.join(fitr)
        if len(momlist) == 0: momlist = self.NNQFull_Stats.index.get_level_values('momentum')
        if thistflowlist == 'PreDef':
            thistflowlist = self.tflowfit
        elif not isinstance(thistflowlist,list):
            thistflowlist = [thistflowlist]
        ip,iflow = momlist[0],thistflowlist[0]
        fitr_list = fitr.split('_')

        fit_data = self.Get_Extrapolation()

        this_fun = fit_data.index[0][3]
        if thispar == 'First':
            thispar = fit_data.index[0][6]
        elif isinstance(thispar,int):
            thispar = fit_data.index.levels[6][thispar]
        this_key = (ip,iflow,this_fun,fitr_list[0],fitr_list[1],thispar)
        if this_key not in fit_data:
            print(this_key , 'Not found, attempting to fit now')
            self.FitAlpha(fit_range=fitr_list[0],thistflowlist = [iflow],tsum_ranges=fitr_list[1],WipeFit=WipeFit)
        fit_data = self.Get_Extrapolation()

        this_fun = fit_data.index[0][3]
        if thispar == 'First':
            thispar = fit_data.index[0][6]
        elif isinstance(thispar,int):
            thispar = fit_data.index.levels[6][thispar]
        this_key = (ip,iflow,this_fun,fitr_list[0],fitr_list[1])
        if this_key not in fit_data:
            print(list(fit_data.index))
            print(this_key)
            raise IOError('error with fitting')
        hold_series = null_series
        hold_series['y_data'] = fit_data.apply(lambda x : x.Avg)
        hold_series['yerr_data'] = fit_data.apply(lambda x : x.Avg)
        hold_series['key_select'] = this_key
        hold_series['type'] = 'hline_vary'
        hold_series['label'] = r'ACTUAL_VALUE_\alpha_{fit}'
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = plotshift
        hold_series['ShowPar'] = r'\alpha'
        hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
        plot_class.AppendData(hold_series)
        return plot_class

    def AlphaPlotDist(self,plot_class,tflow,tsum,tsink,spec_tsum=None,mom='p000',
                        thiscol='PreDefine',Auto=False):
        # autocorrelation has not fits (currenlty)
        if Auto: return plot_class
        self.CheckCol(thiscol)
        if self.is_NNQ2sub:
            ikey = (mom,tsink,tflow,tsum,spec_tsum)
            if 'Alphaboot' not in self.NNQ_subNNQ2_Stats or ikey not in self.NNQ_subNNQ2_Stats['Alphaboot']:
                print(ikey ,'is not in NNQFull_Stats, skipping distribution plot')
                return plot_class
            this_boot = self.NNQ_subNNQ2_Stats.loc[:,'Alphaboot']
        else:
            ikey = (mom,tsink,tflow,tsum)
            if 'Alphaboot' not in self.NNQFull_Stats or ikey not in self.NNQFull_Stats['Alphaboot']:
                print(ikey ,'is not in NNQFull_Stats, skipping distribution plot')
                return plot_class
            this_boot = self.NNQFull_Stats.loc[:,'Alphaboot']
        # this_boot = self.NNQFull_Stats.loc[ikey,'Alphaboot']
        # thisleglab = self.LegLab +'$\ '+ '\ '.join([mom,tsink,tflowTOLeg(tflow),tsum]) + '$'
        thisleglab = self.LegLab #+'$\ '+ '\ '.join([mom,tsink,tflowTOLeg(tflow),tsum]) + '$'
        hold_series = pa.Series()
        hold_series['type'] = 'histogram_vary'
        hold_series['shift'] = self.thisshift
        hold_series['symbol'] = self.thissym
        hold_series['key_select'] = ikey
        hold_series['boot_data'] = this_boot
        hold_series['label'] = thisleglab
        hold_series['color'] = self.thiscol
        hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
        plot_class.AppendData(hold_series)
        return plot_class

    #state is 1 or 2 state fit
    ## fitr is fit range e.g. fitr5-10
    ## fit range of 'Data' is just the original fit range
    def AlphaPlotVsEnergy(  self,plot_class,fitr,momlist=['p000'],thistflowlist = 'PreDef',
                            thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',
                            Auto=False,PhysEnergy=True,WipeFit=False,thispar=r'\alpha'):
        # autocorrelation has not fits (currenlty)
        if self.is_NNQ2sub:
            print('PlotAlpha or PlotAlphaTsum are only done for NNQ2Sub')
            return plot_class
        if Auto: return plot_class
        if fitr == 'None': return plot_class
        if isinstance(fitr,(list,tuple,np.ndarray)):
            fitr = '_'.join(fitr)
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        thisshift = self.GetShift(plot_class.get_xlim(),thisshift)
        if len(momlist) == 0: momlist = self.NNQFull_Stats.index.get_level_values('momentum')
        if thistflowlist == 'PreDef':
            thistflowlist = self.tflowfit
        elif not isinstance(thistflowlist,list):
            thistflowlist = [thistflowlist]
        imom,iflow = momlist[0],thistflowlist[0]
        fitr_list = fitr.split('_')

        fit_data = self.Get_Extrapolation()

        this_fun = fit_data.index[0][3]
        if thispar == 'First':
            thispar = fit_data.index[0][6]
        elif isinstance(thispar,int):
            thispar = fit_data.index.levels[6][thispar]
        this_key = (imom,iflow,this_fun,fitr_list[0],fitr_list[1],thispar)
        if this_key not in fit_data:
            print(this_key , 'Not found, attempting to fit now')
            self.FitAlpha(fit_range=fitr_list[0],thistflowlist = [iflow],tsum_ranges=fitr_list[1],WipeFit=WipeFit)
        fit_data = self.Get_Extrapolation()

        this_fun = fit_data.index[0][3]
        if thispar == 'First':
            thispar = fit_data.index[0][6]
        elif isinstance(thispar,int):
            thispar = fit_data.index.levels[6][thispar]
        this_key = (imom,iflow,this_fun,fitr_list[0],fitr_list[1])
        if this_key not in fit_data:
            print(list(fit_data.index))
            print(this_key)
            raise IOError('error with fitting')
        fit_data = self.NN.latparams.Fmt_pmom_index(fit_data,0)
        this_key_select = (slice(None),) + this_key[1:]

        # ploty.append(this_fit.fit_data['ParamsAvg'][r'\alpha'])
        # plotyerr.append(this_fit.fit_data['ParamsStd'][r'\alpha'])
        # thisplotx = [self.NN.latparams.EnergyFromRecRel(kmom,self.NN.latparams.GetNucleonMass(Phys=False)) for kmom in self.NNQFull_Fit_Stats.loc[(slice(None),iflow,fitr),'boot'].iterkeys()]
        # ## into GeV
        # if PhysEnergy: thisplotx = map(lambda x : x*self.NN.latparams.hbarcdivlat,thisplotx)
        # plotshift = self.GetShift([thisplotx[0],thisplotx[-1]],thisshift)
        hold_series = null_series
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = fit_data.apply(lambda x : x.Avg)
        hold_series['yerr_data'] = fit_data.apply(lambda x : x.Std)
        hold_series['key_select'] = this_key_select
        hold_series['type'] = 'error_bar_vary'
        hold_series['label'] = self.LegLab
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
        plot_class.AppendData(hold_series)
        return plot_class

    ## Comparisons

    def __str__(self):
        return self.name


    def __iter__(self):
        return self

    def items(self):
        return list(self.NNQFull_Stats['boot'].items())

    def iteritems(self):
        return iter(self.NNQFull_Stats['boot'].items())


    def itemsAvgStd(self):
        outlist = []
        for ivals in self.NNQFull_Stats.itertuples():
            outlist.append((ivals[0],ivals[2],ivals[3]))
        return outlist


    def values(self):
        return self.NNQFull_Stats['boot'].values

    def valuesAvgStd(self):
        outlist = []
        for ivals in self.NNQFull_Stats.itertuples(index=False):
            outlist.append((ivals[1],ivals[2]))
        return outlist


    def keys(self):
        return self.NNQFull_Stats.index


    def __setitem__(self, key, value ):
        self.NNQFull_Stats.loc[key,'boot'] = value

    def __getitem__(self, key):
        return self.NNQFull_Stats.loc[key,'boot']

    def __reversed__(self):
        self.NNQFull_Stats = self.NNQFull_Stats.iiloc[::-1]
        return self

    def __next__(self):
        if self.current >= len(list(self.keys())):
            self.current = 0
            raise StopIteration
        else:
            self.current += 1
            if self.current >= len(list(self.keys()))+1:
                return 0.0
            else:
                thiskey = list(self.keys())[self.current-1]
                return self[thiskey]



    def EffMass(self):
        if 'boot' not in self.NNQFull_Stats: return
        lEff,lEffA,lEffS = [],[],[]
        ilist = []
        for (ip,itflow,it_sum),tdata in self.NNQFull_Stats['boot'].groupby(level=('momentum','flow_time','flow_op_trange')):
            next_tdata = np.roll(tdata.swaplevel(1,3)[ip,it_sum,itflow],1)
            for (it,idata),idata_shift in zip(iter(tdata.swaplevel(1,3)[ip,it_sum,itflow].items()),next_tdata):
                thisEff = (idata/idata_shift).Log()
                thisEff.Stats()
                ilist.append((ip,it,itflow,it_sum))
                lEff.append(thisEff)
                lEffA.append(thisEff.Avg)
                lEffS.append(thisEff.Std)
        indicies = pa.MultiIndex.from_tuples(ilist,names=self.NNQFull_col_names)
        self.NNQFull_Stats.loc[:,'EffM'] = pa.Series(lEff,index=indicies)
        self.NNQFull_Stats.loc[:,'EffMAvg'] = pa.Series(lEffA,index=indicies)
        self.NNQFull_Stats.loc[:,'EffMStd'] = pa.Series(lEffS,index=indicies)

    def RatFunAndDer(self):
        if self.is_Qsub or self.is_Q2sub:
            def chit(*vals):
                return vals[0]/vals[1]-vals[2]
            def chitder(*vals):
                return [1.0/vals[1],-vals[0]/(vals[1]**2.),-vals[1]/vals[1]]
        elif self.is_P5NNQsub:
            if 'OnlySub' in self.tr_sum_type:
                def chit(*vals):
                    return (vals[2]*vals[3])/vals[1]
                def chitder(*vals):
                    return [vals[0]-vals[0],(vals[2]*vals[3])/(vals[1]**2.),-vals[3]/vals[1],-vals[2]/vals[1]]
            else:
                def chit(*vals):
                    return (vals[0]-(vals[2]*vals[3]))/vals[1]
                def chitder(*vals):
                    return [1.0/vals[1],((vals[2]*vals[3])-vals[0])/(vals[1]**2.),-vals[3]/vals[1],-vals[2]/vals[1]]
        elif self.is_NNQ2sub or self.is_NNSingQ2sub:
            def chit(*vals):
                return (vals[0]+vals[2])/vals[1]
            def chitder(*vals):
                return [1.0/vals[1],-(vals[2]+vals[0])/(vals[1]**2.),1.0/vals[1]]
        # elif self.is_NNQ2sub or self.is_NNSingQ2sub:
        #     def chit(*vals):
        #         return (vals[0]-vals[2])/vals[1]
        #     def chitder(*vals):
        #         return [1.0/vals[1],(vals[2]-vals[0])/(vals[1]**2.),-1.0/vals[1]]
        elif self.is_NNQsub:
            def chit(*vals):
                return (vals[0]-vals[2])/vals[1]
            def chitder(*vals):
                return [1.0/vals[1],(vals[2]-vals[0])/(vals[1]**2.),-1.0/vals[1]]
        else:
            def chit(*vals):
                return vals[0]/vals[1]
            def chitder(*vals):
                return [1.0/vals[1],-vals[0]/(vals[1]**2.)]
        return (chit,chitder)

# add instead of subtract
    # def RatFunAndDer(self):
    #     if self.is_Qsub or self.is_Q2sub:
    #         def chit(*vals):
    #             return vals[0]/vals[1]+vals[2]
    #         def chitder(*vals):
    #             return [1.0/vals[1],-vals[0]/(vals[1]**2.),vals[1]/vals[1]]
    #     elif self.is_NNQ2sub or self.is_NNSingQ2sub or self.is_NNQsub:
    #         def chit(*vals):
    #             return (vals[0]+vals[2])/vals[1]
    #         def chitder(*vals):
    #             return [1.0/vals[1],(vals[2]+vals[0])/(vals[1]**2.),1.0/vals[1]]
    #     else:
    #         def chit(*vals):
    #             return vals[0]/vals[1]
    #         def chitder(*vals):
    #             return [1.0/vals[1],-vals[0]/(vals[1]**2.)]
    #     return (chit,chitder)

    def AlphaRatio(self,DefWipeG2=False,DoAuto=True,show_timer= False):
        cpevenInfo = deepcopy(self.Info)
        cpevenInfo['Interp'] = 'CPEven'
        cpevenInfo['pmom'] = 'All'

        if self.IsComb:
            if DoAuto:
                ## NOT FEASABLE, too much memory me thinks...
                DoAuto = False
                # raise EnvironmentError('Autocorrelation analysis not implemented for combining Alpha Fulls')

            ## TODO, autocorrelation needs to be done for combining correlators...
            from SetsOfTwoPtCorrs import SetOfTwoPt
            sinklab,jsmcoeffs,sm2ptList,InfoList = GetInfoList(cpevenInfo)
            C2Set = SetOfTwoPt(cfglist=self.NNQFull_cfgs[['configs','xsrc_list']],InfoDict = InfoList)
            C2Set.LoadPickle(DefWipe=DefWipeG2,CheckCfgs=True,CheckMom=True)
            # for iset in C2Set.SetC2.values():
            #     if len(iset.C2) == 0 and DoAuto:
            #         iset.LoadPickle(DefWipe=True,CheckCfgs=True,CheckMom=True,WipeData=not DoAuto)
            thisname,thisleglab,jsmlab = CreateNameAndLL(C2Set.SetC2,sinklab)
            thisfun = GenSinkFun(jsmcoeffs)
            thisfile = self.filename.replace(self.jsm,sinklab)
            thisLL = self.LegLab.replace(self.jsm,sinklab)
            self.SetCustomName(thisfile,thisLL)
            self.jsm = jsmlab
            self.NNCPEven = C2Set.CombineCorrs(thisfun,filename=thisfile,LegLab=thisLL,jsmlab=jsmlab)
            thisname = self.NNCPEven.name + '_ForAlpha'
            self.NNCPEven.SetCustomName(string=thisname)
            self.NNCPEven.Stats()
        else:
            self.NNCPEven = TwoPointCorr(self.nboot,self.NNQFull_cfgs[['configs','xsrc_list']],cpevenInfo,self.thissym,self.thiscol,self.thisshift)
            # self.NNCPEven = self.NNCPEven.GetCombCorr(cpevenInfo)
            thisname = self.NNCPEven.name + '_ForAlpha'
            self.NNCPEven.SetCustomName(string=thisname)
            self.NNCPEven.LoadPickle(DefWipe=DefWipeG2,CheckCfgs=True,CheckMom=True,WipeData=not DoAuto)
            # if 'C2' not in self.NNCPEven.C2_cfgs and DoAuto:
            #     self.NNCPEven.LoadPickle(DefWipe=True,CheckCfgs=True,CheckMom=True,WipeData=not DoAuto)
        print('NNCPEven complete: ', self.NNCPEven.name)


        self.CPEvenFile = self.NNCPEven.HumanFile
        # if DoAuto and 'NNQFull' in self.NNQFull_cfgs and self.IsComb:
        #     raise EnvironmentError('Autocorrelation analysis not implemented for Alpha Full')

        if self.is_NNQ2sub:
            self.RatioNNQ2sub(DefWipeG2=DefWipeG2,DoAuto=DoAuto,show_timer= show_timer)
            return


        lAl,lAlAvg,lAlStd = [],[],[]
        ilist = []
        this_gen = itertools.product(   enumerate(self.NN.pform),enumerate(self.boot_nt_list),
                                        enumerate(self.boot_tsum_list),enumerate(self.FOFull.tflowlist))
        if show_timer: thistimer = Timer(linklist=list(this_gen),name='Autcor '+self.name)
        this_gen = itertools.product(   enumerate(self.NN.pform),enumerate(self.boot_nt_list),
                                        enumerate(self.boot_tsum_list),enumerate(self.FOFull.tflowlist))
        for ((icp,ip),(ict,it),(ict_sum,it_sum),(ictflow,itflow)) in this_gen:
            it_sum = tstr(it_sum).replace('t','ts')
            it = tstr(it+1)
            if (ip,it) not in self.NNCPEven.C2_Stats['boot'].index:
                if show_timer: thistimer.Lap()
                continue
            if ip not in self.pform:
                if show_timer: thistimer.Lap()
                continue

            idata = self.NNQFull_Stats.loc[(ip,it,itflow,it_sum),'boot']
            ilist.append((ip,it,itflow,it_sum))
            idata.Stats()
            thisboot = idata/self.NNCPEven.C2_Stats.loc[(ip,it),'boot']
            if self.is_Qsub or self.is_Q2sub or self.is_NNSingQ2sub or self.is_P5NNQsub:
                q_data = self.NNQFull_cfgs['FlowOp'].apply(lambda x : x[icp][ict][ictflow][ict_sum])
                q_boot = BootStrap(self.nboot, name='q('+','.join([ip,it,itflow,it_sum])+')',cfgvals=q_data)
                if self.is_P5NNQsub:
                    P5NN_data = self.NNQFull_cfgs['P5C2'].apply(lambda x : x[icp][ict][ictflow][ict_sum])
                    P5NN_boot = BootStrap(self.nboot, name='C2('+','.join([ip,it,itflow,it_sum])+')',cfgvals=P5NN_data)
                    if 'OnlySub' in self.tr_sum_type:
                        thisboot = ((q_boot * P5NN_boot)/self.NNCPEven.C2_Stats.loc[(ip,it),'boot'])
                    else:
                        thisboot = thisboot - ((q_boot * P5NN_boot)/self.NNCPEven.C2_Stats.loc[(ip,it),'boot'])
                else:
                    thisboot = thisboot - q_boot
            thisboot.Stats()
            # print 'DEBUG'
            # print 'iscomb',self.IsComb
            # print 'nncpeven',self.NNCPEven
            # print 'nncpodd',self.NN
            # print 'flow',self.FOFull
            # print (ip,it,itflow, idata.Avg,idata.Std,
            #                     self.NNCPEven.C2_Stats['boot'][ip,it].Avg,self.NNCPEven.C2_Stats['boot'][ip,it].Std,
            #                     thisboot.Avg, thisboot.Std)
            lAl.append(thisboot)
            lAlAvg.append(thisboot.Avg)
            lAlStd.append(thisboot.Std)
            if show_timer: thistimer.Lap((ip,it,itflow,it_sum))
        indicies = pa.MultiIndex.from_tuples(ilist,names=self.NNQFull_col_names)
        self.NNQFull_Stats.loc[:,'Alphaboot'] = pa.Series(lAl,index=indicies)
        self.NNQFull_Stats.loc[:,'AlphaAvg'] = pa.Series(lAlAvg,index=indicies)
        self.NNQFull_Stats.loc[:,'AlphaStd'] = pa.Series(lAlStd,index=indicies)

        lAl,lAlAvg,lAlStd = [],[],[]
        ilist = []
        this_gen = itertools.product(enumerate(self.NN.pform),enumerate(self.boot_nt_list),enumerate(self.FOFull.tflowlist))
        if show_timer: thistimer = Timer(linklist=list(this_gen),name='Autcor Int '+self.name)
        this_gen = itertools.product(enumerate(self.NN.pform),enumerate(self.boot_nt_list),enumerate(self.FOFull.tflowlist))
        for ((icp,ip),(ict,it),(ictflow,itflow)) in this_gen:
            it = tstr(it+1)
            if (ip,it) not in self.NNCPEven.C2_Stats['boot'].index:
                if show_timer: thistimer.Lap()
                continue
            if ip not in self.pform:
                if show_timer: thistimer.Lap()
                continue
            if (ip,it,itflow) not in self.NNQInt_Stats.index:
                if show_timer: thistimer.Lap()
                continue
            idata = self.NNQInt_Stats.loc[(ip,it,itflow),'boot']
            if not isinstance(idata,BootStrap): continue
            ilist.append((ip,it,itflow))
            idata.Stats()
            thisboot = idata/self.NNCPEven.C2_Stats.loc[(ip,it),'boot']
            if self.is_Qsub or self.is_Q2sub or self.is_NNSingQ2sub or self.is_P5NNQsub:
                q_data = self.NNQFull_cfgs['FlowOpInt'].apply(lambda x : x[icp][ict][ictflow])
                q_boot = BootStrap(self.nboot, name='q('+','.join([ip,it,itflow])+')',cfgvals=q_data)
                if self.is_P5NNQsub:
                    P5NN_data = self.NNQFull_cfgs['C2Int'].apply(lambda x : x[icp][ict][ictflow])
                    P5NN_boot = BootStrap(self.nboot, name='C2('+','.join([ip,it,itflow])+')',cfgvals=P5NN_data)
                    if 'OnlySub' in self.tr_sum_type:
                        thisboot = ((q_boot * P5NN_boot)/self.NNCPEven.C2_Stats.loc[(ip,it),'boot'])
                    else:
                        thisboot = thisboot - ((q_boot * P5NN_boot)/self.NNCPEven.C2_Stats.loc[(ip,it),'boot'])
                else:
                    thisboot = thisboot - q_boot
            thisboot.Stats()
            # print 'DEBUG'
            # print 'iscomb',self.IsComb
            # print 'nncpeven',self.NNCPEven
            # print 'nncpodd',self.NN
            # print 'flow',self.FOFull
            # print (ip,it,itflow, idata.Avg,idata.Std,
            #                     self.NNCPEven.C2_Stats['boot'][ip,it].Avg,self.NNCPEven.C2_Stats['boot'][ip,it].Std,
            #                     thisboot.Avg, thisboot.Std)
            lAl.append(thisboot)
            lAlAvg.append(thisboot.Avg)
            lAlStd.append(thisboot.Std)
            if show_timer: thistimer.Lap((ip,it,itflow))
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.NNQInt_col_names)
            self.NNQInt_Stats.loc[:,'Alphaboot'] = pa.Series(lAl,index=indicies)
            self.NNQInt_Stats.loc[:,'AlphaAvg'] = pa.Series(lAlAvg,index=indicies)
            self.NNQInt_Stats.loc[:,'AlphaStd'] = pa.Series(lAlStd,index=indicies)
        else:
            print()
            print('WARNING: no NNQInt was found for:')
            print(self.name)
            print()


        if DoAuto and 'NNQFull' in self.NNQFull_cfgs and not self.IsComb:
            thisfuns = self.RatFunAndDer()
            lAl,lAlAvg,lAlStd = [],[],[]
            ilist = []
            this_gen = itertools.product(   enumerate(self.NN.pform),enumerate(self.boot_nt_list),
                                            enumerate(self.boot_tsum_list),enumerate(self.FOFull.tflowlist))
            if show_timer: thistimer = Timer(linklist=list(this_gen),name='Autcor '+self.name)
            this_gen = itertools.product(   enumerate(self.NN.pform),enumerate(self.boot_nt_list),
                                            enumerate(self.boot_tsum_list),enumerate(self.FOFull.tflowlist))
            for ((icp,ip),(ict,it),(ict_sum,it_sum),(ictflow,itflow)) in this_gen:
                tsum_str = tstr(it_sum).replace('t','ts')
                it_str = tstr(it+1)
                # cant do these checks XO.... dictionaries are better!
                # if ip not in self.NNCPEven.C2.keys(): continue
                # if it not in self.NNCPEven.C2[ip].keys(): continue
                # if ('C2' not in self.NNCPEven.C2_cfgs
                #     or icp*self.nt+ict >= self.NNCPEven.C2_cfgs['C2'].size): continue
                auto_df = self.NNQFull_cfgs['NNQFull'].apply(lambda x : x[icp][ict][ictflow][ict_sum]).to_frame(name='NNQFull')
                auto_df['CPEven'] = pa.Series(self.NNCPEven.C2_cfgs['C2'].apply(lambda x : x[icp][it]).values,index=self.NNQFull_cfgs.index)
                if self.is_improved:
                    auto_df['FlowOp'] = self.NNQFull_cfgs['FlowOp'].apply(lambda x : x[icp][ict][ictflow][ict_sum])
                    if self.is_P5NNQsub:
                        auto_df['C2'] = self.NNQFull_cfgs['P5C2'].apply(lambda x : x[icp][ict][ictflow][ict_sum])

                # print 'debug',idata[5],twoptdata[5]
                # if len(idata) != len(twoptdata):
                #     errorstr = '''
                #     NNQ and NN have different configuration lists, \n
                #     maybe force delete the pickled files for NN, \n
                #     or select a different -stats from- directory corresponding to operator Q or W
                #     '''
                #     raise EnvironmentError(errorstr)
                ilist.append((ip,it_str,itflow,tsum_str))

                thisAuto =  AutoCorrelate(  Fun=thisfuns,Sparam=self.Sparam,
                                            name=   self.name + ip+' $t='+str(it)+ \
                                                    '\ '+tflowTOLeg(itflow)+'$ '+\
                                                    tsum_str,data=auto_df,save_covar=True)
                lAl.append(thisAuto)
                lAlAvg.append(thisAuto.Avg)
                lAlStd.append(thisAuto.Std)
                if show_timer: thistimer.Lap((ip,it_str,itflow,tsum_str))
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.NNQFull_col_names)
            self.NNQFull_Stats.loc[:,'AlphaAuto'] = pa.Series(lAl,index=indicies)
            self.NNQFull_Stats.loc[:,'AlphaAutoAvg'] = pa.Series(lAlAvg,index=indicies)
            self.NNQFull_Stats.loc[:,'AlphaAutoStd'] = pa.Series(lAlStd,index=indicies)


            lAl,lAlAvg,lAlStd = [],[],[]
            ilist = []
            this_gen = itertools.product(enumerate(self.NN.pform),enumerate(self.boot_nt_list),enumerate(self.FOFull.tflowlist))
            if show_timer: thistimer = Timer(linklist=list(this_gen),name='Autcor Int '+self.name)
            this_gen = itertools.product(enumerate(self.NN.pform),enumerate(self.boot_nt_list),enumerate(self.FOFull.tflowlist))
            for ((icp,ip),(ict,it),(ictflow,itflow)) in this_gen:
                it_str = tstr(it+1)
                # cant do these checks XO.... dictionaries are better!
                # if ip not in self.NNCPEven.C2.keys(): continue
                # if it not in self.NNCPEven.C2[ip].keys(): continue
                # if ('C2' not in self.NNCPEven.C2_cfgs
                #     or icp*self.nt+ict >= self.NNCPEven.C2_cfgs['C2'].size): continue
                auto_df = self.NNQFull_cfgs['NNQInt'].apply(lambda x : x[icp][ict][ictflow]).to_frame(name='NNQInt')
                auto_df['CPEven'] = pa.Series(self.NNCPEven.C2_cfgs['C2'].apply(lambda x : x[icp][it]).values,index=self.NNQFull_cfgs.index)
                if self.is_improved:
                    auto_df['FlowOpInt'] = self.NNQFull_cfgs['FlowOpInt'].apply(lambda x : x[icp][ict][ictflow])
                    if self.is_P5NNQsub:
                        auto_df['C2Int'] = self.NNQFull_cfgs['C2Int'].apply(lambda x : x[icp][ict][ictflow])
                # print 'debug',idata[5],twoptdata[5]
                # if len(idata) != len(twoptdata):
                #     errorstr = '''
                #     NNQ and NN have different configuration lists, \n
                #     maybe force delete the pickled files for NN, \n
                #     or select a different -stats from- directory corresponding to operator Q or W
                #     '''
                #     raise EnvironmentError(errorstr)
                ilist.append((ip,it_str,itflow))

                thisAuto =  AutoCorrelate(  Fun=thisfuns,Sparam=self.Sparam,
                                            name=self.name + ip+' $t='+str(it)+'\ '+tflowTOLeg(itflow)+'$',data=auto_df,save_covar=True)
                lAl.append(thisAuto)
                lAlAvg.append(thisAuto.Avg)
                lAlStd.append(thisAuto.Std)
                if show_timer: thistimer.Lap((ip,it_str,itflow))
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.NNQInt_col_names)
            self.NNQInt_Stats.loc[:,'AlphaAuto'] = pa.Series(lAl,index=indicies)
            self.NNQInt_Stats.loc[:,'AlphaAutoAvg'] = pa.Series(lAlAvg,index=indicies)
            self.NNQInt_Stats.loc[:,'AlphaAutoStd'] = pa.Series(lAlStd,index=indicies)
        if 'NNQFull' in self.NNQFull_cfgs: del self.NNQFull_cfgs['NNQFull']
        if 'NNQInt' in self.NNQFull_cfgs: del self.NNQFull_cfgs['NNQInt']
        if 'FlowOp' in self.NNQFull_cfgs: del self.NNQFull_cfgs['FlowOp']
        if 'FlowOpInt' in self.NNQFull_cfgs: del self.NNQFull_cfgs['FlowOpInt']
        if 'P5C2' in self.NNQFull_cfgs: del self.NNQFull_cfgs['P5C2']
        if 'C2Int' in self.NNQFull_cfgs: del self.NNQFull_cfgs['C2Int']
        self.RemoveVals()

    def RatioNNQ2sub(self,DefWipeG2=False,DoAuto=True,show_timer= False):
        lAl,lAlAvg,lAlStd = [],[],[]
        ilist = []
        this_gen = itertools.product(   enumerate(self.NN.pform),enumerate(self.boot_nt_list),
                                        enumerate(self.boot_tsum_list),
                                        enumerate(self.boot_spec_tsum_list),
                                        enumerate(self.FOFull.tflowlist))
        if show_timer: thistimer = Timer(linklist=list(this_gen),name='Bootstrap '+self.name)
        this_gen = itertools.product(   enumerate(self.NN.pform),enumerate(self.boot_nt_list),
                                        enumerate(self.boot_tsum_list),
                                        enumerate(self.boot_spec_tsum_list),
                                        enumerate(self.FOFull.tflowlist))
        for ((icp,ip),(ict,it),(ict_sum,it_sum),(jct_sum,jt_sum),(ictflow,itflow)) in this_gen:
            it_sum = tstr(it_sum).replace('t','ts')
            jt_sum = tstr(jt_sum).replace('t','tss')
            it = tstr(it+1)
            if (ip,it) not in self.NNCPEven.C2_Stats['boot'].index:
                if show_timer: thistimer.Lap()
                continue
            if ip not in self.pform:
                if show_timer: thistimer.Lap()
                continue

            idata = self.NNQFull_Stats.loc[(ip,it,itflow,it_sum),'boot']
            ilist.append((ip,it,itflow,it_sum,jt_sum))
            # idata.Stats()
            NNq2_data = self.NNQFull_cfgs['NNQ2'].apply(lambda x : x[icp][ict][ictflow][ict_sum][jct_sum])
            NNq2_boot = BootStrap(self.nboot, name='NNQ2('+','.join([ip,it,itflow,it_sum,jt_sum])+')',cfgvals=NNq2_data)
            thisboot = (idata - NNq2_boot) / self.NNCPEven.C2_Stats.loc[(ip,it),'boot']
            thisboot.Stats()
            # print 'DEBUG'
            # print 'iscomb',self.IsComb
            # print 'nncpeven',self.NNCPEven
            # print 'nncpodd',self.NN
            # print 'flow',self.FOFull
            # print (ip,it,itflow, idata.Avg,idata.Std,
            #                     self.NNCPEven.C2_Stats['boot'][ip,it].Avg,self.NNCPEven.C2_Stats['boot'][ip,it].Std,
            #                     thisboot.Avg, thisboot.Std)
            lAl.append(thisboot)
            lAlAvg.append(thisboot.Avg)
            lAlStd.append(thisboot.Std)
            if show_timer: thistimer.Lap((ip,it,itflow,it_sum,jt_sum))
        indicies = pa.MultiIndex.from_tuples(ilist,names=self.NNQ_subNNQ2_col_names)
        self.NNQ_subNNQ2_Stats.loc[:,'Alphaboot'] = pa.Series(lAl,index=indicies)
        self.NNQ_subNNQ2_Stats.loc[:,'AlphaAvg'] = pa.Series(lAlAvg,index=indicies)
        self.NNQ_subNNQ2_Stats.loc[:,'AlphaStd'] = pa.Series(lAlStd,index=indicies)

        lAl,lAlAvg,lAlStd = [],[],[]
        ilist = []
        this_gen = itertools.product(   enumerate(self.NN.pform),
                                        enumerate(self.boot_nt_list),
                                        enumerate(self.FOFull.tflowlist))
        if show_timer: thistimer = Timer(linklist=list(this_gen),name='Autcor Int '+self.name)
        this_gen = itertools.product(   enumerate(self.NN.pform),
                                        enumerate(self.boot_nt_list),
                                        enumerate(self.FOFull.tflowlist))
        for ((icp,ip),(ict,it),(ictflow,itflow)) in this_gen:
            it = tstr(it+1)
            if (ip,it) not in self.NNCPEven.C2_Stats['boot'].index:
                if show_timer: thistimer.Lap()
                continue
            if ip not in self.pform:
                if show_timer: thistimer.Lap()
                continue
            idata = self.NNQInt_Stats.loc[(ip,it,itflow),'boot']
            ilist.append((ip,it,itflow))
            # idata.Stats()
            NNq2_data = self.NNQFull_cfgs['NNQ2Int'].apply(lambda x : x[icp][ict][ictflow])
            NNq2_boot = BootStrap(self.nboot, name='NNQ2('+','.join([ip,it,itflow])+')',cfgvals=NNq2_data)
            thisboot = (idata - NNq2_boot)/self.NNCPEven.C2_Stats.loc[(ip,it),'boot']
            thisboot.Stats()
            # print 'DEBUG'
            # print 'iscomb',self.IsComb
            # print 'nncpeven',self.NNCPEven
            # print 'nncpodd',self.NN
            # print 'flow',self.FOFull
            # print (ip,it,itflow, idata.Avg,idata.Std,
            #                     self.NNCPEven.C2_Stats['boot'][ip,it].Avg,self.NNCPEven.C2_Stats['boot'][ip,it].Std,
            #                     thisboot.Avg, thisboot.Std)
            lAl.append(thisboot)
            lAlAvg.append(thisboot.Avg)
            lAlStd.append(thisboot.Std)
            if show_timer: thistimer.Lap((ip,it,itflow))
        indicies = pa.MultiIndex.from_tuples(ilist,names=self.NNQInt_col_names)
        self.NNQInt_Stats.loc[:,'Alphaboot'] = pa.Series(lAl,index=indicies)
        self.NNQInt_Stats.loc[:,'AlphaAvg'] = pa.Series(lAlAvg,index=indicies)
        self.NNQInt_Stats.loc[:,'AlphaStd'] = pa.Series(lAlStd,index=indicies)


        if DoAuto and 'NNQFull' in self.NNQFull_cfgs and not self.IsComb:
            thisfuns = self.RatFunAndDer()
            lAl,lAlAvg,lAlStd = [],[],[]
            ilist = []
            this_gen = itertools.product(   enumerate(self.NN.pform),
                                            enumerate(self.boot_nt_list),
                                            enumerate(self.boot_tsum_list),
                                            enumerate(self.boot_spec_tsum_list),
                                            enumerate(self.FOFull.tflowlist))
            if show_timer: thistimer = Timer(linklist=list(this_gen),name='Autcor '+self.name)
            this_gen = itertools.product(   enumerate(self.NN.pform),
                                            enumerate(self.boot_nt_list),
                                            enumerate(self.boot_tsum_list),
                                            enumerate(self.boot_spec_tsum_list),
                                            enumerate(self.FOFull.tflowlist))
            for ((icp,ip),(ict,it),(ict_sum,it_sum),(jct_sum,jt_sum),(ictflow,itflow)) in this_gen:
                tsum_str = tstr(it_sum).replace('t','ts')
                jtsum_str = tstr(jt_sum).replace('t','tss')
                it_str = tstr(it+1)
                # cant do these checks XO.... dictionaries are better!
                # if ip not in self.NNCPEven.C2.keys(): continue
                # if it not in self.NNCPEven.C2[ip].keys(): continue
                # if ('C2' not in self.NNCPEven.C2_cfgs
                #     or icp*self.nt+ict >= self.NNCPEven.C2_cfgs['C2'].size): continue
                auto_df = self.NNQFull_cfgs['NNQFull'].apply(lambda x : x[icp][ict][ictflow][ict_sum]).to_frame(name='NNQFull')
                auto_df['CPEven'] = pa.Series(self.NNCPEven.C2_cfgs['C2'].apply(lambda x : x[icp][it]).values,index=self.NNQFull_cfgs.index)
                auto_df['NNQ2'] = self.NNQFull_cfgs['NNQ2'].apply(lambda x : x[icp][ict][ictflow][ict_sum][jct_sum])
                # print 'debug',idata[5],twoptdata[5]
                # if len(idata) != len(twoptdata):
                #     errorstr = '''
                #     NNQ and NN have different configuration lists, \n
                #     maybe force delete the pickled files for NN, \n
                #     or select a different -stats from- directory corresponding to operator Q or W
                #     '''
                #     raise EnvironmentError(errorstr)
                ilist.append((ip,it_str,itflow,tsum_str,jtsum_str))

                thisAuto =  AutoCorrelate(  Fun=thisfuns,Sparam=self.Sparam,
                                            name=   self.name + ip+' $t='+str(it)+ \
                                                    '\ '+tflowTOLeg(itflow)+'$ '+tsum_str+jtsum_str,
                                            data=auto_df,save_covar=True)
                lAl.append(thisAuto)
                lAlAvg.append(thisAuto.Avg)
                lAlStd.append(thisAuto.Std)
                if show_timer: thistimer.Lap((ip,it_str,itflow,tsum_str,jtsum_str))
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.NNQ_subNNQ2_col_names)
            self.NNQ_subNNQ2_Stats.loc[:,'AlphaAuto'] = pa.Series(lAl,index=indicies)
            self.NNQ_subNNQ2_Stats.loc[:,'AlphaAutoAvg'] = pa.Series(lAlAvg,index=indicies)
            self.NNQ_subNNQ2_Stats.loc[:,'AlphaAutoStd'] = pa.Series(lAlStd,index=indicies)


            lAl,lAlAvg,lAlStd = [],[],[]
            ilist = []
            this_gen = itertools.product(   enumerate(self.NN.pform),
                                            enumerate(self.boot_nt_list),
                                            enumerate(self.FOFull.tflowlist))
            if show_timer: thistimer = Timer(linklist=list(this_gen),name='Autcor Int '+self.name)
            this_gen = itertools.product(   enumerate(self.NN.pform),
                                            enumerate(self.boot_nt_list),
                                            enumerate(self.FOFull.tflowlist))
            for ((icp,ip),(ict,it),(ictflow,itflow)) in this_gen:
                it_str = tstr(it+1)
                # cant do these checks XO.... dictionaries are better!
                # if ip not in self.NNCPEven.C2.keys(): continue
                # if it not in self.NNCPEven.C2[ip].keys(): continue
                # if ('C2' not in self.NNCPEven.C2_cfgs
                #     or icp*self.nt+ict >= self.NNCPEven.C2_cfgs['C2'].size): continue
                auto_df = self.NNQFull_cfgs['NNQInt'].apply(lambda x : x[icp][ict][ictflow]).to_frame(name='NNQInt')
                auto_df['CPEven'] = pa.Series(self.NNCPEven.C2_cfgs['C2'].apply(lambda x : x[icp][it]).values,index=self.NNQFull_cfgs.index)
                auto_df['NNQ2Int'] = self.NNQFull_cfgs['NNQ2Int'].apply(lambda x : x[icp][ict][ictflow])
                # print 'debug',idata[5],twoptdata[5]
                # if len(idata) != len(twoptdata):
                #     errorstr = '''
                #     NNQ and NN have different configuration lists, \n
                #     maybe force delete the pickled files for NN, \n
                #     or select a different -stats from- directory corresponding to operator Q or W
                #     '''
                #     raise EnvironmentError(errorstr)
                ilist.append((ip,it_str,itflow))

                thisAuto =  AutoCorrelate(  Fun=thisfuns,Sparam=self.Sparam,
                                            name=self.name + ip+' $t='+str(it)+'\ '+tflowTOLeg(itflow)+'$',data=auto_df,save_covar=True)
                lAl.append(thisAuto)
                lAlAvg.append(thisAuto.Avg)
                lAlStd.append(thisAuto.Std)
                if show_timer: thistimer.Lap((ip,it,itflow))
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.NNQInt_col_names)
            self.NNQInt_Stats.loc[:,'AlphaAuto'] = pa.Series(lAl,index=indicies)
            self.NNQInt_Stats.loc[:,'AlphaAutoAvg'] = pa.Series(lAlAvg,index=indicies)
            self.NNQInt_Stats.loc[:,'AlphaAutoStd'] = pa.Series(lAlStd,index=indicies)
        if 'NNQFull' in self.NNQFull_cfgs: del self.NNQFull_cfgs['NNQFull']
        if 'NNQInt' in self.NNQFull_cfgs: del self.NNQFull_cfgs['NNQInt']
        if 'NNQ2' in self.NNQFull_cfgs: del self.NNQFull_cfgs['NNQ2']
        if 'NNQ2Int' in self.NNQFull_cfgs: del self.NNQFull_cfgs['NNQ2Int']
        self.RemoveVals()


    # def ReadAndBootAndEffM(self,DefWipe=False):
    #     if os.path.isfile(self.PickleFile) and not DefWipe:
    #         self.LoadPickle()
    #     else:
    #         self.ReadAndBoot()
    #         self.EffMass()

    def RemoveAllBoots(self):
        self.NNQFull_Stats = self.NNQFull_Stats.applymap(RemoveAllBoots)


    def Write(self,WipeData=True):
        def FixDictArray(thislist,ikey):
            outDict = OrderedDict()
            for ic,ilist in enumerate(thislist):
                outDict[ikey+str(ic+1)] = ilist
            return outDict

        # self.FOFull.FlowWrite()
        # if not self.IsComb:
        #     self.NN.Write()

        if not hasattr(self,'show_cfgs'):
            self.show_cfgs = False

        outDict = ODNested()
        if self.IsComb:
            outDict['TwoPtCorr_File'] = 'Not outputted'
        #     self.NN.HumanFile = self.NN.HumanFile.replace(self.NN.jsm,self.jsm)
        else:
            outDict['TwoPtCorr_File'] = self.NN.HumanFile
        outDict['FlowOp_File'] = self.FOFull.flowHumanFile
        if self.CPEvenFile != 'None':
            outDict['TwoPtCorr_CPEven_file'] = self.CPEvenFile
        excel_params = pa.Series(deepcopy(outDict))

        outDict['Fits']['Function'] = self.Fun.__name__
        if self.AlphaiGuess == 'None':
            outDict['Fits']['Initial_Guess'] = 'Default (see FitFunctions.py)'
        else:
            for iparam in range(1,self.npar+1):
                thispar = 'Par'+str(iparam)
                outDict['Fits']['Initial_Guess'][thispar] = self.AlphaiGuess[iparam]
        if 'boot' in self.NNQFull_Fit_Stats:
            for (ip,itflow),fitdict in self.NNQFull_Fit_Stats['boot'].items():
                if not isinstance(fitdict,sff.SetOfFitFuns): continue
                outDict['Fits'][ip][itflow] = fitdict.GetOutputDict()


        # if 'boot' in self.Tau_Fit_Stats:
        #     for (ifitr,imom,itsink,itflow),idata in self.Tau_Fit_Stats['boot'].iteritems():
        #         if isinstance(idata,ff.Fitting):
        #             continue
        #         outDict['Tau_Fits'][imom][itflow][itsink][ifitr] = fitdict.GetOutputDict()


        if self.is_NNQ2sub:

            for col_key,col_data in self.NNQ_subNNQ2_Stats.items():
                if 'Avg' in col_key or 'Std' in col_key: continue
                for (ip,it,itflow,it_sum,jt_sum),idata in col_data.items():
                    outDict[col_key][ip][it][itflow][it_sum][jt_sum] = AvgStdToFormat(idata.Avg,idata.Std)
        else:
            for col_key,col_data in self.NNQFull_Stats.items():
                if 'Avg' in col_key or 'Std' in col_key: continue
                for (ip,it,itflow,it_sum),idata in col_data.items():
                    outDict[col_key][ip][it][itflow][it_sum] = AvgStdToFormat(idata.Avg,idata.Std)

        for col_key,col_data in self.NNQInt_Stats.items():
            if 'Avg' in col_key or 'Std' in col_key: continue
            for (ip,it,itflow),idata in col_data.items():
                if hasattr(idata,'Avg'):
                    outDict[col_key+'_Int'][ip][it][itflow] = AvgStdToFormat(idata.Avg,idata.Std)

        if any(np.array(self.ncfg_list) > 0) and self.show_cfgs:
            for ((istream,iccfg),icfg),ixsrc_list in zip(iter(self.NNQFull_cfgs['configs'].items()),self.NNQFull_cfgs['xsrc_list'].values):
                outDict['cfglist'][istream][icfg] = 'n_xsrc_'+str(len(ixsrc_list))

        if WipeData: self.RemoveVals()

        ## Write human readable file with data
        WriteXml(self.HumanFile,{'Results':outDict})
        ## pickles rest of data for reading

        if any(np.array(self.ncfg_list) > 0):
            for ((istream,iccfg),icfg),ixsrc_list in zip(iter(self.NNQFull_cfgs['configs'].items()),self.NNQFull_cfgs['xsrc_list'].values):
                outDict['cfglist'][istream][icfg] = 'n_xsrc_'+str(len(ixsrc_list))
        ## TODO, include fitting in the excell file
        WriteExcel( self.ExcelFile,{'NNQFull_results':deepcopy(self.NNQFull_Stats),
                                    'NNQInt_results':deepcopy(self.NNQInt_Stats)},
                    cfglist=self.NNQFull_cfgs,params=excel_params)

        self.RemoveFuns()
        self.NNQFull_Stats.to_pickle(self.PickleFile+'.stats')
        self.NNQFull_Stats_name = self.PickleFile+'.stats'
        out_dict = deepcopy(self.__dict__)
        del out_dict['NNQFull_Stats']
        WritePickle(self.PickleFile,out_dict)
        self.GetFuns()
        ##################################

    def ReadAndWrite(self,WipeData=True,SaveMem=False,NoFit=False):
        self.Read()
        self.Bootstrap(WipeData=WipeData)
        self.Stats()
        if NoFit:
            self.FitAlpha(WipeFit=WipeData)
        if SaveMem:
            self.RemoveAllBoots()
        self.Write()


    def FixRead(self,readdict):
        readdict['Info']['outdir'] = self.Info['outdir']
        readdict['FOFull'].latparams.Info = self.FOFull.latparams.Info
        readdict['FOFull'].latparams.outputdir = self.FOFull.latparams.outputdir
        readdict['FOFull'].latparams.PickleFile = self.FOFull.latparams.PickleFile
        readdict['FOFull'].latparams.HumanFile = self.FOFull.latparams.HumanFile
        readdict['NN'].latparams.outputdir = self.NN.latparams.outputdir
        readdict['NN'].latparams.PickleFile = self.NN.latparams.PickleFile
        readdict['NN'].latparams.HumanFile = self.NN.latparams.HumanFile
        readdict['FOFull'].flowHumanFile = self.FOFull.flowHumanFile
        readdict['FOFull'].flowPickleFile = self.FOFull.flowPickleFile
        readdict['NN'].HumanFile = self.NN.HumanFile
        readdict['NN'].PickleFile = self.NN.PickleFile
        readdict['NNQFull_Fit_col_names'] = self.NNQFull_Fit_col_names
        ## this is to check updated function
        if hasattr(self,'Fun') and 'Fun' in list(readdict.keys()):
            if self.Fun.__name__ != readdict['Fun'].__name__:
                readdict['Fun'] = self.Fun
                readdict['npar'] = self.npar
                readdict['NNQFull_Fit_Stats'] = pa.DataFrame()
        readdict['PickleFile'] = self.PickleFile
        readdict['ExcelFile'] = self.ExcelFile
        readdict['HumanFile'] = self.HumanFile
        readdict['NNQFull_Fit_col_names'] = self.NNQFull_Fit_col_names
        readdict['pform'] = self.pform
        readdict['boot_tsum_list'] = self.boot_tsum_list
        if not os.path.isfile(readdict['NNfile']):
            readdict['NNfile'] = self.NNfile
        if not os.path.isfile(readdict['FOFullfile']):
            readdict['FOFullfile'] = self.FOFullfile
        if 'is_improved' not in list(readdict.keys()):
            readdict['is_improved'] = False
        # if not os.path.isfile(readdict['ppparams'].PickleFile):
        #     readdict['ppparams'].PickleFile = self.ppparams.PickleFile
        #     readdict['ppparams'].HumanFile = self.ppparams.HumanFile
        #     readdict['qparams'].PickleFile = self.qparams.PickleFile
        #     readdict['qparams'].HumanFile = self.qparams.HumanFile
        #     readdict['Info']['outdir'] = self.Info['outdir']
        return readdict

    ## TODO: Pickle load and dump self.NN and self.FOFull from their respective directories.
    ##       At this point, I am just dumping this whole class instead. (doubles up data)
    def LoadPickle(self,WipeData=True,DefWipe = False,CheckCfgs=False,CheckMom=True,legacy=False,SaveMem=False,NoFit=False):
        print('Loading Pickle',self.PickleFile)
        read_pickle = self.PickleFile.replace('_'.join([self.NN.MesOrBar,self.NN.stream_name])+'_','')
        if os.path.isfile(self.PickleFile) and not DefWipe and not legacy:
            loadeddict = ReadPickleWrap(self.PickleFile)
            loadeddict = self.FixRead(loadeddict)
            checklist = self.thischecklist
            if CheckCfgs: checklist = checklist + ['NNQFull_cfgs']
            if CheckMom: checklist = checklist + ['checkmomlist']
            loadeddict['checkmomlist'] = [self.NN.latparams.pformTOpstr(ip,nopref=True) for ip in loadeddict['NNQInt_Stats'].index.get_level_values('momentum')]
            self.GetFuns()
            if CheckClass(self.__dict__,loadeddict,checklist) and 'AlphaAuto' in loadeddict['NNQInt_Stats']:
                print('loading pickle:')
                print(self.PickleFile)
                print()
                self.__dict__.update(loadeddict)
                if hasattr(self,'NNQFull_Stats_name') and os.path.isfile(self.NNQFull_Stats_name):
                    self.NNQFull_Stats = pa.read_pickle(self.NNQFull_Stats_name)
                else:
                    self.Write()
                if Wipe_All_Fits:
                    self.NNQFull_Fit_Stats = pa.DataFrame()
                elif len(self.NNQFull_Fit_Stats.values)> 0 and not isinstance(self.NNQFull_Fit_Stats.values[0][0],sff.SetOfFitFuns):
                    print('Legacy file, fits need to be recomputed')
                    self.NNQFull_Fit_Stats = pa.DataFrame()
            else:
                if os.path.isfile(self.PickleFile+'.bak'):
                    print('using backupfile for '+self.name)
                    self.PickleFile = self.PickleFile+'.bak'
                    self.HumanFile = self.HumanFile+'.bak'
                    self.ExcelFile = self.ExcelFile.replace('.xlsx','.bak.xlsx')
                    self.LoadPickle(DefWipe=DefWipe,WipeData=WipeData,CheckCfgs=CheckCfgs,CheckMom=CheckMom)
                    return
                print('Warning, file ' , self.PickleFile , ' has different parameters than this instance, reading and writing over file now')
                print('Attempting legacy file read')
                print()
                self.temp = self.PickleFile,self.HumanFile,self.ExcelFile
                self.PickleFile = self.PickleFile.replace('.bak','')
                self.HumanFile = self.HumanFile.replace('.bak','')
                self.ExcelFile = self.ExcelFile.replace('.bak','')
                self.LoadPickle(DefWipe=DefWipe,WipeData=WipeData,CheckCfgs=CheckCfgs,CheckMom=CheckMom,legacy=True)
                return
                # self.ReadAndWrite(WipeData=WipeData)
            # self.Write()
        elif os.path.isfile(read_pickle) and not DefWipe:
            loadeddict = ReadPickleWrap(read_pickle)
            loadeddict = self.FixRead(loadeddict)
            checklist = self.thischecklist
            if CheckCfgs: checklist = checklist + ['NNQFull_cfgs']
            if CheckMom: checklist = checklist + ['checkmomlist']
            loadeddict['checkmomlist'] = [self.NN.latparams.pformTOpstr(ip,nopref=True) for ip in loadeddict['NNQFull_Stats'].index.get_level_values('momentum')]
            self.GetFuns()
            if CheckClass(self.__dict__,loadeddict,checklist) and 'AlphaAuto' in loadeddict['NNQInt_Stats']:
                print('Legacy file found, reading and writing to proper file')
                print(read_pickle)
                print()
                self.__dict__.update(loadeddict)
                if SaveMem:
                    self.RemoveAllBoots()
                if hasattr(self,'temp'):
                    self.PickleFile,self.HumanFile,self.ExcelFile = self.temp
                    del self.temp
                self.Write()
            else:
                if os.path.isfile(self.PickleFile+'.bak'):
                    print('using legacy backupfile for '+self.name)
                    self.PickleFile = self.PickleFile+'.bak'
                    self.HumanFile = self.HumanFile+'.bak'
                    self.ExcelFile = self.ExcelFile.replace('.xlsx','.bak.xlsx')
                    self.LoadPickle(DefWipe=DefWipe,WipeData=WipeData,CheckCfgs=CheckCfgs,CheckMom=CheckMom,legacy=True)
                    return
                print('Warning, legacy file ' , read_pickle , ' has different parameters than this instance, reading and writing over file now')
                print()
                if hasattr(self,'temp'):
                    self.PickleFile,self.HumanFile,self.ExcelFile = self.temp
                    del self.temp
                self.ReadAndWrite(WipeData=WipeData,SaveMem=SaveMem,NoFit=NoFit)
        else:
            if not DefWipe:
                print('Pickle file not found:')
                print(self.PickleFile)
                print('or legacy file')
                print(read_pickle)
                print()
            if hasattr(self,'temp'):
                self.PickleFile,self.HumanFile,self.ExcelFile = self.temp
                del self.temp
            self.ReadAndWrite(WipeData=WipeData,SaveMem=SaveMem,NoFit=NoFit)
    ## Operator overloading


    ## unitary arithmetic operations

    def __neg__(self):
        for ival in list(self.keys()):
            self[ival] = -self[ival]
        return self


    def __pos__(self):
        return self

    def __abs__(self):
        for ival in list(self.keys()):
            self[ival] = abs(self[ival])
        return self

    ## Casting, just casts the bootstrap values, ignore return value (only needed since return has to be same type).


    def __complex__(self):
        for ival in list(self.keys()):
            complex(self[ival])
        self.Stats()
        return 1+0j

    def __int__(self):
        for ival in list(self.keys()):
            int(self[ival])
        self.Stats()
        return 1

    def __float__(self):
        for ival in list(self.keys()):
            np.float64(self[ival])
        self.Stats()
        return 1.0

    def sin(self):
        result = NNQFullCorr(thisnboot=self.nboot, cfglist=self.NNQFull_cfgs[['configs','xsrc_list']],
                         Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.NNQFull_Stats['boot'] = np.sin(self.NNQFull_Stats.loc[:,'boot'])
        result.NNQInt_Stats['boot'] = np.sin(self.NNQInt_Stats.loc[:,'boot'])
        return result

    def cos(self):
        result = NNQFullCorr(thisnboot=self.nboot, cfglist=self.NNQFull_cfgs[['configs','xsrc_list']],
                         Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.NNQFull_Stats['boot'] = np.cos(self.NNQFull_Stats.loc[:,'boot'])
        result.NNQInt_Stats['boot'] = np.cos(self.NNQInt_Stats.loc[:,'boot'])
        return result


    def Overload_wrap(self,NNQFull2,this_fun):
        result = NNQFullCorr(thisnboot=self.nboot, cfglist=self.NNQFull_cfgs[['configs','xsrc_list']],
                              Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.NNQFull_Stats = pa.DataFrame(columns=['boot'],index=self.NNQFull_Stats.index)
        result.NNQInt_Stats = pa.DataFrame(columns=['boot'],index=self.NNQInt_Stats.index)
        # result.UpdateName('+',self,NNQFull2)
        if isinstance(NNQFull2,NNQFullCorr):
            # for (ip,it,itflow,it_sum),iNNQFull in self.NNQFull_Stats['boot'].iteritems():
            #     print 'DEBUG',ip,it,itflow,it_sum
            #     if (ip,it,itflow,it_sum) not in NNQFull2.NNQFull_Stats.index: continue
                result.NNQFull_Stats.loc[:,'boot'] = this_fun(self.NNQFull_Stats.loc[:,'boot'],
                                                                NNQFull2.NNQFull_Stats.loc[:,'boot'])
                if  'boot' in self.NNQInt_Stats.columns:
                    result.NNQInt_Stats.loc[:,'boot'] = this_fun(   self.NNQInt_Stats.loc[:,'boot'],
                                                                    NNQFull2.NNQInt_Stats.loc[:,'boot'])
                    result.NNQInt_Stats = result.NNQInt_Stats.dropna()
        elif isinstance(NNQFull2,NNQCorr):
            for (ip,it,itflow,it_sum),iNNQFull in self.NNQFull_Stats['boot'].items():
                if (ip,it,itflow) not in NNQFull2.NNQ_Stats.index: continue
                result[(ip,it,itflow,it_sum)] = this_fun(iNNQFull, NNQFull2[(ip,it,itflow)])
                if 'boot' in self.NNQInt_Stats.columns and it_sum == self.boot_tsum_list[0]:
                    result.NNQInt_Stats.loc[(ip,it,itflow),'boot'] = this_fun(result.NNQInt_Stats.loc[(ip,it,itflow),'boot'],
                                                                                NNQFull2[(ip,it,itflow)])
        elif isinstance(NNQFull2,TwoPointCorr):
            for (ip,it,itflow,it_sum),iNNQFull in self.NNQFull_Stats['boot'].items():
                if (ip,it) not in NNQFull2.C2_Stats.index: continue
                result[(ip,it,itflow,it_sum)] = this_fun(iNNQFull, NNQFull2[(ip,it)])
                if 'boot' in self.NNQInt_Stats.columns and it_sum == self.boot_tsum_list[0]:
                    result.NNQInt_Stats.loc[(ip,it,itflow),'boot'] = this_fun(result.NNQInt_Stats.loc[(ip,it,itflow),'boot'],
                                                                        NNQFull2[(ip,it)])

        elif isinstance(NNQFull2,FlowOp):
            for (ip,it,itflow,it_sum),iNNQFull in self.NNQFull_Stats['boot'].items():
                if itflow not in NNQFull2.Op_Stats.index: continue
                result[(ip,it,itflow,it_sum)] = this_fun(iNNQFull, NNQFull2[itflow])
                if 'boot' in self.NNQInt_Stats.columns and it_sum == self.boot_tsum_list[0]:
                    result.NNQInt_Stats.loc[(ip,it,itflow),'boot'] = this_fun(result.NNQInt_Stats.loc[(ip,it,itflow),'boot'], NNQFull2[itflow])
        else:
            try:
                result.NNQFull_Stats.loc[:,'boot'] = this_fun(self.NNQFull_Stats.loc[:,'boot'], NNQFull2)
                if 'boot' in self.NNQInt_Stats.columns :
                    result.NNQInt_Stats.loc[:,'boot'] = this_fun(result.NNQInt_Stats.loc[:,'boot'], NNQFull2)
            except Exception as e:

                print(type(NNQFull2), NNQFull2)
                print(self.NNQFull_Stats)
                print(self.NNQInt_Stats)
                print(e)
                raise EnvironmentError('Invalid value to combine with NNQFullCorr class')
        return result


    def __add__(self,NNQFull2):
        if any([ipar in str(type(NNQFull2)) for ipar in NNQFullCorr.parentlist]):
            return NNQFull2.__radd__(self)
        else:
            return self.Overload_wrap(NNQFull2,op.add)

    def __sub__(self,NNQFull2):
        if any([ipar in str(type(NNQFull2)) for ipar in NNQFullCorr.parentlist]):
            return NNQFull2.__rsub__(self)
        else:
            return self.Overload_wrap(NNQFull2,op.sub)

    def __mul__(self,NNQFull2):
        if any([ipar in str(type(NNQFull2)) for ipar in NNQFullCorr.parentlist]):
            return NNQFull2.__rmul__(self)
        else:
            return self.Overload_wrap(NNQFull2,op.mul)

    def __div__(self,NNQFull2):
        if any([ipar in str(type(NNQFull2)) for ipar in NNQFullCorr.parentlist]):
            return NNQFull2.__rdiv__(self)
        else:
            return self.Overload_wrap(NNQFull2,op.truediv)

    def __pow__(self,NNQFull2):
        if any([ipar in str(type(NNQFull2)) for ipar in NNQFullCorr.parentlist]):
            return NNQFull2.__rpow__(self)
        else:
            return self.Overload_wrap(NNQFull2,op.pow)


    ## Right multiplication functions


    def __radd__(self,NNQFull2):
        if any([ipar in str(type(NNQFull2)) for ipar in NNQFullCorr.parentlist]):
            return NNQFull2.__add__(self)
        else:
            return self.Overload_wrap(NNQFull2,flip_fun_2arg(op.add))

    def __rsub__(self,NNQFull2):
        if any([ipar in str(type(NNQFull2)) for ipar in NNQFullCorr.parentlist]):
            return NNQFull2.__sub__(self)
        else:
            return self.Overload_wrap(NNQFull2,flip_fun_2arg(op.sub))

    def __rmul__(self,NNQFull2):
        if any([ipar in str(type(NNQFull2)) for ipar in NNQFullCorr.parentlist]):
            return NNQFull2.__mul__(self)
        else:
            return self.Overload_wrap(NNQFull2,flip_fun_2arg(op.mul))

    def __rdiv__(self,NNQFull2):
        if any([ipar in str(type(NNQFull2)) for ipar in NNQFullCorr.parentlist]):
            return NNQFull2.__div__(self)
        else:
            return self.Overload_wrap(NNQFull2,flip_fun_2arg(op.truediv))


    def __rpow__(self,NNQFull2):
        if any([ipar in str(type(NNQFull2)) for ipar in NNQFullCorr.parentlist]):
            return NNQFull2.__pow__(self)
        else:
            return self.Overload_wrap(NNQFull2,flip_fun_2arg(op.pow))


class FlowedTwoPtCorr(object):
    """ C2(t_f,t) uses bootstrap class
        t_f = flow time
        t = time

    """

    ## Info is a dictionary containing information for the correlator:
    ## see comments below to see what is required.
    ## missing elements just ignores that field

    ## TODO: create variational flag and interface with this function.

    ## array of all classes that inherit this class. Append this when adding new classes please!
    parentlist = ['SetOfCorrs.SetOfTwoPt','TwoPtCorrelators.NNQCorr','TwoPtCorrelators.NNQFullCorr','VarMethod.VariationalTwoPt']

    def __init__(self, thisnboot=nboot, cfglist={},Info={},thissym='Not Set',thiscol='Not Set',thisshift=0.0,name=''):

        ## these class elemnts are checked when loading in pickled file
        self.thischecklist = ['nboot','nt','kud','ks','tflowlist']

        self.latparams = mp.LatticeParameters(Info=Info)
        self.latparams.LoadPickle()

        self.nt = self.latparams.nt
        self.nxyz = self.latparams.nxyz
        self.current = 0
        ## Initialising plotting parameters
        self.thisshift = thisshift
        self.thiscol = thiscol
        self.thissym = thissym
        ## the actual shift ammount needs to be scaled by the x axis size
        ## see GetShift
        self.thisshiftscale = 'Not Set'

        #self.C2 = [ ip , it , ,istream , iconf ]
        # self.C2     = np.array([])

        self.C2_cfgs = pa.DataFrame()
        '''
        self.C2_cfgs DataFrame:

        columns =   config_numbers , Op[itflow, it]

        rows =      stream , configuration
        '''

        self.C2_Col_Names = ['flow_time','source_sink_separation']
        self.C2_Stats = pa.DataFrame()

        '''
        self.C2_Stats DataFrame:

        columns =   boot,   Avg,        Std
                    EffM,   EffMAvg,    EffMStd

        rows =      flow_time, time separation
        '''



        self.C2_Fit_Col_Names = ['states','flow_time', 'time_fit_range']
        self.C2_Fit_Stats = pa.DataFrame()

        '''
        self.C2_Fit_Stats DataFrame:

        columns =   boot,   Avg,        Std

        rows =      states fitted to, flow_time, time fit range

        Avg/Std elements are pa.Series instances with:
            index = fit parameters
            values = Avg/Std values
        '''


        # # FitDict = [ States_Fitted_To, ip, Fit_range ] fitting class
        # self.FitDict    = ODNested()
        # # FitDictAvg/Std [ States_Fitted_To , ip ,  Fit_range , Fit_Parameter ] bs
        # self.FitDictAvg = ODNested()
        # self.FitDictStd = ODNested()


        if 'Adjoint_Flow' in list(Info.keys()): self.Adjoint = Info['Adjoint_Flow']
        else: self.Adjoint = False

        self.nboot = thisnboot


        ## toggle to show exact configurations and sources used in xml outputfiles
        if 'show_cfgs' in list(Info.keys()): self.show_cfgs = Info['show_cfgs']
        else: self.show_cfgs = False


        ## used to pick out different configuraiton/source lists.
        ## for example, if there are less flowed operators than 2 point correlators, passing 'TopCharge' will read that result
        if 'fo_for_cfgs' in list(Info.keys()): self.fo_for_cfgs = Info['fo_for_cfgs']
        else: self.fo_for_cfgs = False



        ## Info can contain fit ranges to be calculated
        ## must bein form Info['Fits'] = [ ('state#' , 'fit#-#' ) , ...]
        if 'Fits' in list(Info.keys()): self.PredefFits = Info['Fits']
        else: self.PredefFits = []



        ## Info can contain projector used (or type of meson in psibar Gamma psi)
        if 'MesOrBar' in list(Info.keys()): self.MesOrBar = Info['MesOrBar']
        else: self.MesOrBar = 'Meson' ## Default to baryon


        ## Info can contain initial guess for fits, make sure it array of correct length
        if 'iGuess' in list(Info.keys()): self.iGuess = np.array(Info['iGuess'])
        else: self.iGuess = 'None' ## Defaults to value set for function in FitFunctions.py

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


        # ## x_source list has been pushed to cfglist, this is depriciated
        # if 'xsrc' in Info.keys(): self.xsrc = ['xsrc'+str(isrc) f
        # else: self.xsrc = ['xsrc1']

        ## kappa up down (in integer form)
        if 'kud' in list(Info.keys()): self.kud = 'kud'+str(Info['kud'])
        else: self.kud = ''

        ## kappa strange (in integer form)
        if 'ks' in list(Info.keys()): self.ks = 'ks'+str(Info['ks'])
        else: self.ks = ''
        listout = []

        ## this is used display only these flowtimes, Still reads and writes all
        if 'tflowlist' in list(Info.keys()): self.tflowlist = Info['tflowlist']
        else: self.tflowlist = defxlimOp ## default to read all flow times

        for ict,itflow in enumerate(self.tflowlist):
            if 't_f' not in str(itflow):
                listout.append(tflowstr(itflow))
            else:
                listout.append(itflow)
        self.tflowlist = np.array(listout)

        self.kappafolder = 'Kud0'+self.kud.replace('kud','')+'Ks0'+self.ks.replace('ks','')
        self.dim_label = 'RC'+str(self.nxyz)+'x'+str(self.nt)
        self.dim_ll = str(self.nxyz)+r'^{3}\times ' + str(self.nt)
        self.flowfiledir = datadir+'Flow/'+self.dim_label+'_'+self.kappafolder+'-*-/'

        self.Create_Streams_Dirs()
        ## cfglist { configuration : [ x-source numbers ] }
        if isinstance(cfglist,(list,np.ndarray,pa.DataFrame)):
            self.ImportCfgList(cfglist)
        else:
            self.GetCfgList()
        self.Info = Info



        self.SetCustomName(name)
        self.LogScale = 1.0


    # ## thisInfo must contain correct information for combining the correlators (the smearings, the eigenvectors, etc...)
    # def GetCombCorr(self,thisInfo):
    #     from SetsOfTwoPtCorrs import SetOfTwoPt
    #     sinklab,jsmcoeffs,sm2ptList,InfoList = GetInfoList(thisInfo)

    #     C2Set = SetOfTwoPt(cfglist=self.NNQ_cfgs[['configs','xsrc_list']],InfoDict = InfoList)
    #     C2Set.LoadPickle(CheckCfgs=CheckCfgs)
    #     thisname,thisleglab,jsmlab = CreateNameAndLL(C2Set.SetC2,sinklab)
    #     thisfun = GenSinkFun(jsmcoeffs)
    #     thisfile = self.filename.replace(self.jsm,sinklab)
    #     thisLL = self.LegLab.replace(self.jsm,sinklab)
    #     # self.SetCustomName(thisfile,thisLL)
    #     # self.jsm = jsmlab
    #     return C2Set.CombineCorrs(thisfun,filename=thisfile,LegLab=thisLL,jsmlab=jsmlab)

    def Create_Streams_Dirs(self):
        if isinstance(self.stream_list,str) and self.stream_list == 'all':
            self.All_Streams()
        else:
            self.Check_Streams()
        self.stream_name = '-'+''.join([istr.replace('-','') for istr in self.stream_list])+'-'
        self.dir_list = {}
        for istream in self.stream_list:
            thisdir = self.flowfiledir.replace('-*-',istream)
            if os.path.isdir(thisdir):
                self.dir_list[istream] = thisdir
            else:
                self.dir_list[istream] = thisdir.replace(self.dim_label+'_','')

    def All_Streams(self):
        self.stream_list = []
        for ifolder in glob.glob(self.flowfiledir):
            self.stream_list.append(re.findall('-.*-',ifolder)[-1])
        for ifolder in glob.glob(self.flowfiledir.replace(self.dim_label+'_','')):
            this_stream = re.findall('-.*-',ifolder)[-1]
            if this_stream not in self.stream_list:
                self.stream_list.append(this_stream)
        if len(self.stream_list) == 0:
            raise IOError(self.flowfiledir+'\n No Streams found')

    def Check_Streams(self):
        for istream in self.stream_list:
            this_check = os.path.isdir(self.flowfiledir.replace('-*-',istream))
            this_check = this_check or os.path.isdir(self.flowfiledir.replace('-*-',istream).replace(self.dim_label+'_',''))
            if not this_check:
                raise IOError(  self.flowfiledir.replace('-*-',istream) + '\n stream directory not found, as well as \n'+
                                self.flowfiledir.replace('-*-',istream).replace(self.dim_label+'_',''))


    def GetCfgList(self):
        def CheckCfgFolder(thisdir):
            outbool = self.Adjoint and os.path.isfile(thisdir+'/_Adjoint_'+self.MesOrBar.lower())
            return outbool or os.path.isfile(thisdir+'/_'+self.MesOrBar.lower())

        cfglist,ilist = [],[]
        for istream in self.stream_list:
            idir = self.dir_list[istream]
            if len(glob.glob(idir+'/*')) == 1:
                idir = idir + 'Comb/'
            for icfg,ifile in enumerate(glob.glob(idir+'/*')):
                if 'cfg' in ifile and CheckCfgFolder(ifile+'/'):
                    ilist.append((istream,icfg))
                    cfglist.append(ifile.replace(idir+'cfg',''))
        if len(cfglist) == 0:
            print('No configs found for directory:')
            print(' , '.join(self.dir_list))
        indicies = pa.MultiIndex.from_tuples(ilist,names=['stream','config_number'])
        cfg_df = pa.DataFrame(cfglist,columns=['configs'],index=indicies)
        cfg_df.loc[:,'int_configs'] = pa.Series(list(map(int,cfglist)),index=indicies)
        cfg_df = cfg_df.sort_values('int_configs').sort_index(level='stream',sort_remaining=False)
        del cfg_df['int_configs']
        self.ImportCfgList(cfg_df)



    # def MakeOutputFile(self):
    #     self.name = '_'.join([self.tsrc,self.ism,self.jsm])
    #     self.HumanFile = '/'+self.kappafolder+'/G2/'+self.name+'.xml'
    #     self.PickleFile = '/'+self.kappafolder+'/G2/Pickle/'+self.name+'.py3p'


    ## Routines for Plotting ###############################################

    def ImportPlotParams(self,thiscol,thissym,thisshift):
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        self.thisshift = thisshift

    def CheckCol(self,thiscol):
        if 'PreDefine' == thiscol :
            if 'Not Set' == self.thiscol:
                raise IOError('Pass in color to initialize it')
        else:
            self.thiscol = thiscol

    def CheckSym(self,thissym):
        if 'PreDefine' == thissym :
            if 'Not Set' == self.thissym:
                raise IOError('Pass in symbol to initialize it')
        else:
            self.thissym = thissym


    def GetShift(self,xlims,thisshift):
        if thisshift != 'PreDefine': self.thisshift = thisshift
        # if self.thisshiftscale == 'Not Set':
        xlen = np.abs(xlims[1]-xlims[0])
        self.thisshiftscale = self.thisshift*xlen
        return self.thisshiftscale
        # else:
        #     return self.thisshiftscale




    ## make sure all xlims are the same, because this sets the window of xlim (needed for calculating shifts)
    def LogPlot(self,plot_class,xlims=defxlim,tflowlist=[],thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',norm=True):
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        self.thisshiftscale = 'Not Set'
        thisshift = self.GetShift(xlims,thisshift)
        if len(tflowlist) == 0:
            tflowlist = self.C2_Stats.index.get_level_values('flow_time')
        for itflow in tflowlist:
            if itflow not in self.C2_Stats.index.get_level_values('flow_times'):
                print(itflow, ' not found in effective mass list')
                continue
            pdata = self.C2_Stats['boot'].xs(itflow, level='flow_times')
            tlist,databoot = pdata.index,pdata.values
            dataavg,dataerr = [],[]
            if norm:
                scale = databoot[xlims[0]].Log()
                scale.Stats()
                scale = scale.Avg
            else:
                scale = 1.0
            for idata in databoot[xlims[0]:xlims[1]]:
                logdata = idata.Log()/scale
                logdata.Stats()
                dataavg.append(logdata.Avg)
                dataerr.append(logdata.Std)
            tlist = tlist[xlims[0]:xlims[1]]
            hold_series = pa.Series()
            hold_series['x_data'] = np.array(list(map(untstr,tlist)))
            hold_series['y_data'] = dataavg
            hold_series['yerr_data'] = dataerr
            hold_series['xerr_data'] = None
            hold_series['type'] = 'error_bar'
            hold_series['fit_class'] = None
            hold_series['label'] = self.LegLab+'$\ '+ tflowTOLeg(itflow) + '$'
            hold_series['symbol'] = self.thissym
            hold_series['color'] = self.thiscol
            hold_series['shift'] = thisshift
            hold_series['xdatarange'] = xlims
            hold_series['scale'] = self.LogScale
            hold_series['Phys'] = None
            hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
            plot_class.AppendData(hold_series)
            self.LogScale=scale
        return plot_class
        # plot_class.get_xlim(*xlims)

    ## make sure all xlims are the same, because this sets the window of xlim (needed for calculating shifts)
    def Plot(   self,plot_class,xlims=defxlim,tflowlist=[],
                thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',norm=True):
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        self.thisshiftscale = 'Not Set'
        thisshift = self.GetShift(xlims,thisshift)
        if len(tflowlist) == 0:
            tflowlist = self.C2_Stats.index.get_level_values('flow_time')
        for itflow in tflowlist:
            if itflow not in self.C2_Stats['boot']:
                print(itflow, ' not found in effective mass list')
                continue
            pdata = self.C2_Stats['boot'].xs(itflow, level='flow_time')
            tlist,databoot = pdata.index,pdata.values
            dataavg,dataerr = [],[]
            if norm:
                scale = databoot[xlims[0]]
                scale.Stats()
                scale = scale.Avg
            else:
                scale = 1.0
            for idata in databoot[xlims[0]:xlims[1]]:
                data = idata/scale
                data.Stats()
                dataavg.append(data.Avg)
                dataerr.append(data.Std)
            tlist = tlist[xlims[0]:xlims[1]]
            hold_series = pa.Series()
            hold_series['x_data'] = np.array(list(map(untstr,tlist)))
            hold_series['y_data'] = dataavg
            hold_series['yerr_data'] = dataerr
            hold_series['xerr_data'] = None
            hold_series['type'] = 'error_bar'
            hold_series['fit_class'] = None
            hold_series['label'] = self.LegLab
            hold_series['symbol'] = self.thissym
            hold_series['color'] = self.thiscol
            hold_series['shift'] = thisshift
            hold_series['xdatarange'] = xlims
            hold_series['scale'] = self.LogScale
            hold_series['Phys'] = None
            hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
            plot_class.AppendData(hold_series)
            self.Scale=scale
        return plot_class
        # plot_class.get_xlim(*xlims)


    ## make sure all xlims are the same, because this sets the window of xlim (needed for calculating shifts)
    def EffMassPlot(self,plot_class,xlims = defxlim,tflowlist=[],thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',Phys=True):
        def DoPhys(value):
            thisval = value*self.latparams.hbarcdivlat
            thisval.Stats()
            return thisval

        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        ## resets thisshiftscale so that you can plot EffMass and Log plot next to eachother,

        self.thisshiftscale = 'Not Set'
        thisshift = self.GetShift(xlims,thisshift)
        if len(tflowlist) == 0:
            tflowlist = self.C2_Stats.index.get_level_values('flow_time')
        m_getattr = lambda x, y: getattr(y,x)
        for itflow in tflowlist:
            if itflow not in self.C2_Stats['EffM']:
                print(itflow, ' not found in effective mass list')
                continue
            pdata = self.C2_Stats['EffM'].xs(itflow, level='flow_time')
            if Phys: pdata = pdata.apply(DoPhys)
            dataavg = pdata.apply(partial(m_getattr,'Avg'))
            dataerr = pdata.apply(partial(m_getattr,'Std'))
            tlist = pdata.index[xlims[0]:xlims[1]]
            dataavg = dataavg.iloc[xlims[0]:xlims[1]]
            dataerr = dataerr.iloc[xlims[0]:xlims[1]]
            ## np.abs makes negative effective masses coming from G2^-1 be positive.
            hold_series = pa.Series()
            hold_series['x_data'] = np.array(list(map(untstr,tlist)))
            hold_series['y_data'] = dataavg
            hold_series['yerr_data'] = dataerr
            hold_series['xerr_data'] = None
            hold_series['type'] = 'error_bar'
            hold_series['fit_class'] = None
            hold_series['label'] = self.LegLab+'$\ '+ tflowTOLeg(itflow) + '$'
            hold_series['symbol'] = self.thissym
            hold_series['color'] = self.thiscol
            hold_series['shift'] = thisshift
            hold_series['xdatarange'] = xlims
            hold_series['scale'] = self.LogScale
            hold_series['Phys'] = None
            hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
            plot_class.AppendData(hold_series)
        return plot_class

    ## make sure all xlims are the same, because this sets the window of xlim (needed for calculating shifts)
    def FlowRatPlot(self,plot_class,xlims = defxlim,tflowlist=[],thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine'):
        def DoStats(value):
            # thisval = value*self.latparams.hbarcdivlat
            thisvalue = value
            # thisvalue = thisvalue.Log()
            if isinstance(thisvalue,BootStrap): thisvalue.Stats()
            return thisvalue

        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        ## resets thisshiftscale so that you can plot EffMass and Log plot next to eachother,

        self.thisshiftscale = 'Not Set'
        thisshift = self.GetShift(xlims,thisshift)
        if len(tflowlist) == 0:
            tflowlist = self.C2_Stats.index.get_level_values('flow_time')
        m_getattr = lambda x, y: getattr(y,x)
        for itflow in tflowlist:
            if itflow not in self.C2_Stats['boot']:
                print(itflow, ' not found in effective mass list')
                continue
            pdata = self.C2_Stats['boot'].xs(itflow, level='flow_time')/self.C2_Stats['boot'].xs('t_f0.00', level='flow_time')
            pdata = pdata.apply(DoStats)
            dataavg = pdata.apply(partial(m_getattr,'Avg'))
            dataerr = pdata.apply(partial(m_getattr,'Std'))
            tlist = pdata.index[xlims[0]:xlims[1]]
            dataavg = dataavg.iloc[xlims[0]:xlims[1]]
            dataerr = dataerr.iloc[xlims[0]:xlims[1]]

            ## np.abs makes negative effective masses coming from G2^-1 be positive.
            hold_series = pa.Series()
            hold_series['x_data'] = np.array(list(map(untstr,tlist)))
            hold_series['y_data'] = dataavg
            hold_series['yerr_data'] = dataerr
            hold_series['xerr_data'] = None
            hold_series['type'] = 'error_bar'
            hold_series['fit_class'] = None
            hold_series['label'] = self.LegLab+'$\ '+ tflowTOLeg(itflow) + '$'
            hold_series['symbol'] = self.thissym
            hold_series['color'] = self.thiscol
            hold_series['shift'] = thisshift
            hold_series['xdatarange'] = xlims
            hold_series['scale'] = self.LogScale
            hold_series['Phys'] = None
            hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
            plot_class.AppendData(hold_series)
        return plot_class

    ## make sure all xlims are the same, because this sets the window of xlim (needed for calculating shifts)
    def AllFlowRatPlot(self,plot_class,tsink,thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine'):
        def DoStats(value):
            # thisval = value*self.latparams.hbarcdivlat
            thisvalue = value
            # thisvalue = thisvalue.Log()
            if isinstance(thisvalue,BootStrap): thisvalue.Stats()
            return thisvalue

        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        ## resets thisshiftscale so that you can plot EffMass and Log plot next to eachother,

        self.thisshiftscale = 'Not Set'
        thisshift = self.GetShift([untflowstr(self.tflowlist[0]),untflowstr(self.tflowlist[-1])],thisshift)
        dataavg,dataerr,tlist = [],[],[]
        for itflow in self.tflowlist:
            pdata = self.C2_Stats['boot'].xs(itflow, level='flow_time')[tsink]/self.C2_Stats['boot'].xs('t_f0.00', level='flow_time')[tsink]
            # pdata = pdata.Log()
            pdata.Stats()
            dataavg.append(pdata.Avg)
            dataerr.append(pdata.Std)
            tlist.append(TflowToPhys(np.array([itflow]),self.latparams.latspace)[0])
        ## np.abs makes negative effective masses coming from G2^-1 be positive.
        hold_series = pa.Series()
        hold_series['x_data'] = np.array(list(map(untflowstr,tlist)))
        hold_series['y_data'] = dataavg
        hold_series['yerr_data'] = dataerr
        hold_series['xerr_data'] = None
        hold_series['type'] = 'error_bar'
        hold_series['fit_class'] = None
        hold_series['label'] = self.LegLab+'$\ '+ tsinkTOLeg(tsink) + '$'
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        hold_series['xdatarange'] = None
        hold_series['scale'] = self.LogScale
        hold_series['Phys'] = None
        hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
        plot_class.AppendData(hold_series)
        return plot_class

    ## state is 1 or 2 state fit
    ## fitr is fit range e.g. fitr5-10
    ## fit range of 'Data' is just the original fit range
    ## make sure to plot this after LogPlot or the xlims wont work properly
    def LogFitPlot( self,plot_class,state,fitr,xlims = 'Data',
                    tflowlist=[],thiscol='PreDefine',thisshift='PreDefine'):
        self.CheckCol(thiscol)
        thisshift = self.GetShift(plot_class.get_xlim(),thisshift)
        if len(tflowlist) == 0:
            tflowlist = self.C2_Fit_Stats.index.get_level_values('flow_time')
        if isinstance(state,int): state = 'state'+str(state)
        for itflow in tflowlist:
            if (state,itflow,fitr) not in list(self.C2_Fit_Stats.index):
                print(state,itflow,fitr, ' not found in FitDict, performing fit now')
                self.Fit(int(state.replace('state','')),[fitr])
            hold_series = pa.Series()
            hold_series['x_data'] = None
            hold_series['y_data'] = None
            hold_series['yerr_data'] = None
            hold_series['xerr_data'] = None
            hold_series['type'] = 'log_fit'
            hold_series['fit_class'] = self.C2_Fit_Stats['boot'][state,itflow,fitr]
            hold_series['label'] = self.C2_Fit_Stats['boot'][state,itflow,fitr].name
            hold_series['symbol'] = self.thissym
            hold_series['color'] = self.thiscol
            hold_series['shift'] = thisshift
            hold_series['xdatarange'] = xlims
            hold_series['scale'] = self.LogScale
            hold_series['Phys'] = None
            hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
            plot_class.AppendData(hold_series)
        return plot_class

    ## state is 1 or 2 state fit
    ## fitr is fit range e.g. fitr5-10
    ## fit range of 'Data' is just the original fit range
    ## make sure to plot this after EffMassPlot or the xlims wont work properly
    def EffMassFitPlot( self,plot_class,state,fitr,xlims = 'Data',tflowlist=[],
                        thiscol='PreDefine',thisshift='PreDefine',Phys=True):
        self.CheckCol(thiscol)
        thisshift = self.GetShift(plot_class.get_xlim(),thisshift)
        if len(tflowlist) == 0:
            tflowlist = self.C2_Stats.index.get_level_values('flow_time')
        if isinstance(state,int): state = 'state'+str(state)
        for itflow in tflowlist:
            if (state,itflow,fitr) not in list(self.C2_Fit_Stats.index):
                print(state,itflow,fitr, ' not found in FitDict, performing fit now')
                self.Fit(int(state.replace('state','')),[fitr])
            hold_series = pa.Series()
            hold_series['x_data'] = None
            hold_series['y_data'] = None
            hold_series['yerr_data'] = None
            hold_series['xerr_data'] = None
            hold_series['type'] = 'effm_fit'
            hold_series['fit_class'] = self.C2_Fit_Stats['boot'][state,itflow,fitr]
            hold_series['label'] = self.C2_Fit_Stats['boot'][state,itflow,fitr].name
            hold_series['symbol'] = self.thissym
            hold_series['color'] = self.thiscol
            hold_series['shift'] = thisshift
            hold_series['xdatarange'] = xlims
            hold_series['scale'] = self.LogScale
            hold_series['Phys'] = Phys
            hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
            plot_class.AppendData(hold_series)
        return plot_class
    ## ##################### ###############################################

    ## state is 1 or 2 state fit
    ## fitr is fit range e.g. fitr5-10
    ## fit range of 'Data' is just the original fit range
    ## make sure to plot this after EffMassPlot or the xlims wont work properly
    def RecRelPlot(self,plot_class,state,fitr,thiscol='PreDefine',thisshift='PreDefine',thissym='PreDef',Phys=True):
        def DoPhys(value):
            thisval = value*self.latparams.hbarcdivlat
            thisval.Stats()
            return thisval

        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        if isinstance(state,int): state = 'state'+str(state)
        # fitr_params = self.FitDict[state][fitr]
        ploty,plotyerr,plotx,plotxerr = [],[],[],[]
        for imom in self.C2_Stats.index.get_level_values('flow_time'):
            if (state,imom,fitr) not in list(self.C2_Fit_Stats.index):
                print(state,imom,fitr , ' has not been done, performing fit now')
                self.Fit(int(state.replace('state','')),[fitr])
            this_fit = self.C2_Fit_Stats['boot'][state,imom,fitr]
            if isinstance(this_fit,pa.Series):
                this_fit = this_fit.iloc[0]
            thisboot = this_fit.fit_data['Params']['Energy']
            thiszboot = this_fit.fit_data['Params']['Energy']
            if Phys:
                thisboot = thisboot*self.latparams.hbarcdivlat
                thisboot.Stats()
            ploty.append(thisboot.Avg)
            plotyerr.append(thisboot.Std)
            xboot = self.latparams.EnergyFromRecRel(imom,thiszboot)
            if Phys:
                xboot = xboot*self.latparams.hbarcdivlat
            xboot.Stats()
            plotx.append(xboot.Avg)
            plotxerr.append(xboot.Std)
        thisshift = self.GetShift((min(plotx),max(plotx)),thisshift)
        hold_series = pa.Series()
        hold_series['x_data'] = np.array(plotx)
        hold_series['y_data'] = ploty
        hold_series['yerr_data'] = plotyerr
        hold_series['xerr_data'] = plotxerr
        hold_series['type'] = 'effm_fit'
        hold_series['fit_class'] = None
        hold_series['label'] = self.LegLab+' '+'_'.join((state,fitr))
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        hold_series['xdatarange'] = None
        hold_series['scale'] = self.LogScale
        hold_series['Phys'] = Phys
        hold_series['fmt_class'] = KeyForamtting(self.NN.latparams)
        plot_class.AppendData(hold_series)
        return plot_class



    def SetCustomName(self,string='',stringLL='',fo_for_cfgs='PreDef'):
        if not hasattr(self,'stream_name'):
            self.stream_name = '-def-'
        if string == '':
            self.name = '_'.join([self.dim_label,self.kappafolder,self.stream_name,self.MesOrBar])
            Prefix = 'Flowed'
            if self.Adjoint:
                Prefix = '_Adjoint'
            self.filename = '_'.join([Prefix,self.MesOrBar,self.stream_name])
            # self.filename = self.name
        else:
            self.filename = self.name = string
        if stringLL == '':
            if self.Adjoint:
                self.LegLab = '$'+'\ '.join(['Adjoint',self.dim_ll,self.latparams.GetPionMassLab(),self.stream_name,self.MesOrBar])+'$' ## customise this how you like
            else:
                self.LegLab = '$'+'\ '.join([self.dim_ll,self.latparams.GetPionMassLab(),self.stream_name,self.MesOrBar])+'$' ## customise this how you like
        else:
            self.LegLab = stringLL


        if isinstance(fo_for_cfgs,str):
            if fo_for_cfgs == 'PreDef':
                fo_for_cfgs = self.fo_for_cfgs
        if isinstance(fo_for_cfgs,str):
            self.G2dir = outputdir + '/'+self.dim_label+self.kappafolder+'/G2/'+fo_for_cfgs+'/'
        else:
            self.G2dir = outputdir + '/'+self.dim_label+self.kappafolder+'/G2/'
        mkdir_p(self.G2dir+'/Pickle/')
        mkdir_p(self.G2dir+'/Excel/')
        self.HumanFile = self.G2dir+self.filename+'.xml'
        self.ExcelFile = self.G2dir+'/Excel/'+self.filename+'.xlsx'
        self.PickleFile = self.G2dir+'/Pickle/'+self.filename+'.py3p'

    def ImportCfgList(self,cfglist,stream_list='PreDefine'):
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
        self.C2_cfgs = cfglist
        self.Update_Stream_List()
        self.ncfg_list = [icfg.size for istream,icfg in self.C2_cfgs['configs'].groupby(level='stream')]


    def Update_Stream_List(self):
        if not hasattr(self,'stream_name'):
            self.stream_name = '-def-'
        self.stream_list = self.C2_cfgs.index.levels[0]
        self.nstream = len(self.stream_list)
        if hasattr(self,'stream_name'):
            old_sn = self.stream_name
            self.stream_name = '-'+''.join([istr.replace('-','') for istr in self.stream_list])+'-'
            if hasattr(self,'filename'):
                self.filename = self.filename.replace(old_sn,self.stream_name)
            if hasattr(self,'name'):
                self.name = self.name.replace(old_sn,self.stream_name)
            if hasattr(self,'LegLab'):
                self.LegLab = self.LegLab.replace(old_sn,self.stream_name)
        else:
            self.stream_name = '-'+''.join([istr.replace('-','') for istr in self.stream_list])+'-'


    def Stats(self):
        if 'boot' not in self.C2_Stats: return
        lAvg,lStd = [],[]
        for (it_flow,it),iC2 in self.items():
            iC2.Stats()
            lAvg.append(iC2.Avg)
            lStd.append(iC2.Std)
        indicies = pa.MultiIndex.from_tuples(self.C2_Stats.index,names=self.C2_Col_Names)
        self.C2_Stats.loc[:,'Avg'] = pa.Series(lAvg,index=indicies)
        self.C2_Stats.loc[:,'Std'] = pa.Series(lStd,index=indicies)
        self.EffMass()


    def RemoveVals(self):
        if 'C2' in self.C2_cfgs:del self.C2_cfgs['C2']


    def Read(self,show_timer=True,Full=False):
        # def FixTime(this_data):
        #     output = []
        #     for flow_data in this_data:
        #         max_index = flow_data.index(max(flow_data))
        #         rolled = np.roll(np.array(flow_data),-max_index).tolist()
        #         output.append(rolled)
        #     return output

        def ReadMesonFile(thisfile):
            if not os.path.isfile(thisfile):
                print('Warning, file not found',thisfile)
                return [],[]

            read_the_flow = True
            read_this_flow = False
            tflow_list, meson_list = [],[]
            with open(thisfile,'r') as ifile:
                for line in ifile:
                    strip_line = line.replace('\n','')
                    if len(strip_line) == 0:
                        read_the_flow = True
                        read_this_flow = False
                    elif read_the_flow:
                        read_the_flow = False
                        tflow = tflowstr(strip_line.split(' ')[0])
                        if tflow in self.tflowlist:
                            read_this_flow = True
                            tflow_list.append(tflow)
                            meson_list.append([])
                    elif read_this_flow:
                        t_sep,val,cmplx = strip_line.split(' ')
                        meson_list[-1].append(np.float64(val))
            return tflow_list,meson_list

        ## making interpolator number more human readable
        thisdata = []
        if show_timer: thistimer = Timer(linklist=self.C2_cfgs['configs'].values,name='Read '+self.name)
        thisfile = '_'+self.MesOrBar.lower()
        if self.Adjoint:
            thisfile = '_Adjoint'+thisfile
        for (istream,iccfg),icfg in self.C2_cfgs['configs'].items():
            istream_dir = self.dir_list[istream]
            idir = 'cfg'+icfg
            try:
                this_tflow_list,file_data = ReadMesonFile(istream_dir+'/'+idir+'/'+thisfile)
                if self.CheckTsrc(file_data):
                    err_str = 'TsrcErr With: \n'+istream_dir+'/'+idir+'/'+thisfile + ' \n'
                    raise IOError(err_str)
            except Exception as err:
                err_str = 'Error with file: \n'+istream_dir+'/'+idir+'/'+thisfile + ' \n'
                err_str += str(err)
                raise IOError(err_str)
            if len(this_tflow_list) < len(self.tflowlist):
                print(', '.join(this_tflow_list))
                print(', '.join(self.tflowlist))
                raise IOError('not enough flow times in file \n'+istream_dir+'/'+idir+'/'+thisfile)
            thisdata.append(file_data)
            if show_timer: thistimer.Lap()
        self.C2_cfgs.loc[:,'C2'] = pa.Series(thisdata,index=self.C2_cfgs.index)

    def CheckTsrc(self,this_data):
        this_data = np.array(this_data)
        left_test = np.any(this_data[0,:self.latparams.nt/2-5] < this_data[0,1:self.latparams.nt/2-4])
        right_test = np.any(this_data[0,self.latparams.nt/2+5:-1] > this_data[0,self.latparams.nt/2+6:])
        return left_test or right_test


    def Bootstrap(self,WipeData=True):
        if 'C2' not in self.C2_cfgs:
            self.Read()
        lC2,lC2Avg,lC2Std = [],[],[]
        ilist = []
        rlist = None
        for ictflow,itflow in enumerate(self.tflowlist):
            for it in range(self.latparams.nt):
                strit = tstr(it+1)
                ilist.append((itflow,strit))
                # def debugthis(val):
                #     print val[ictflow][it]
                # print self.C2_cfgs['C2'].apply(debugthis)
                this_lambda = lambda x : x[ictflow][it]
                tdata = self.C2_cfgs['C2'].apply(this_lambda)

                thisboot = BootStrap(self.nboot, name='G_{2}('+itflow+','+strit+')',cfgvals=tdata,rand_list=rlist)
                rlist = thisboot.Get_Rand_List()
                lC2.append(thisboot)
                lC2Avg.append(thisboot.Avg)
                lC2Std.append(thisboot.Std)
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.C2_Col_Names)
            self.C2_Stats.loc[:,'boot'] = pa.Series(lC2,index=indicies)
            self.C2_Stats.loc[:,'Avg'] = pa.Series(lC2Avg,index=indicies)
            self.C2_Stats.loc[:,'Std'] = pa.Series(lC2Std,index=indicies)
        if WipeData: del self.C2_cfgs['C2']


    # def ReadAndBoot(self):
    #     self.Read()
    #     self.Bootstrap()

    def ImportFitRanges(self,fit_range):
        self.fit_range = fit_range
        if isinstance(self.fit_range,str):
            self.fit_range = [self.fit_range]


    ## Fitting the correlator
    ## fit_range is either list of fit ranges to fit the correlator to, or a single fit range
    ## formatted as fitr#-# where # are the min and max fit values.
    ## states correspond to how many states to fit to (currenlty only working for states = 1 or 2)
    def Fit(self,state='PreDef',fit_range='PreDef',iGuess='PreDef',EstDir=False,Ratio=False,WipeFit=True):
        print()
        print('****** Warning******')
        print('TODO: fits in ImpSOF have problem related to exp(mass) in fit, please fix before using')
        print('****** Warning******')
        print()
        def RunPredefFits():
            if state == 'PreDef' and fit_range == 'PreDef':
                for izip in self.PredefFits:
                    if len(izip) == 2:
                        self.Fit(state=izip[0],fit_range=[izip[1]],EstDir=EstDir,Ratio=Ratio,WipeFit=WipeFit)
                    else:
                        self.Fit(state=izip[0],fit_range=[izip[1]],iGuess=izip[2],EstDir=EstDir,Ratio=Ratio,WipeFit=WipeFit)
            elif state == 'PreDef':
                for izip in self.PredefFits:
                    if izip[1] == fit_range:
                        if len(izip) == 2:
                            self.Fit(state=izip[0],fit_range=[izip[1]],EstDir=EstDir,Ratio=Ratio,WipeFit=WipeFit)
                        else:
                            self.Fit(state=izip[0],fit_range=[izip[1]],iGuess=izip[2],EstDir=EstDir,Ratio=Ratio,WipeFit=WipeFit)
            elif fit_range == 'PreDef':
                for izip in self.PredefFits:
                    if izip[0] == state:
                        if len(izip) == 2:
                            self.Fit(state=izip[0],fit_range=[izip[1]],EstDir=EstDir,Ratio=Ratio,WipeFit=WipeFit)
                        else:
                            self.Fit(state=izip[0],fit_range=[izip[1]],iGuess=izip[2],EstDir=EstDir,Ratio=Ratio,WipeFit=WipeFit)

        if 'PreDef' in [state,fit_range]:
            RunPredefFits()
            return
        if iGuess != 'PreDef': self.iGuess = iGuess
        if isinstance(state,str):
            if 'state' in state:
                state = state.replace('state','')
            state = int(state)
        self.ImportFitRanges(fit_range)
        if state == 1:
            # thisfitfun = C2OneStateFitFunNoExp
            thisfitfun = C2OneStateFitFun
            npar = 2
        elif state == 2:
            # thisfitfun = C2TwoStateFitFunNoExp
            thisfitfun = C2TwoStateFitFun
            npar = 4
        elif state > 2:
            raise IOError('fitting state is only implemented up to 2 state so far. Cannot do '+str(state))
        ststr = 'state'+str(state)
        legstate = str(state)+'SF'
        lFit,lFitA,lFitS = [],[],[]
        ilist = []
        for itflow,tflowdata in self.C2_Stats['boot'].groupby(level='flow_time'):
            # print 'Fitting state',str(state) , ' ' , ifit
            for ifit in self.fit_range:
                if (ststr,itflow,ifit) in list(self.C2_Fit_Stats.index):
                    if WipeFit:
                        self.C2_Fit_Stats.drop([(ststr,itflow,ifit)],inplace=True)
                    else:
                        continue
                ilist.append((ststr,itflow,ifit))
                ifitmin,ifitmax = list(map(int,unxmlfitr(ifit)))
                # print '     Fitting ', itflow
                tdatarange = [[]]
                ydata = []
                for it,tdata in tflowdata.items():
                    tint = untstr(it[1])
                    if ifitmin <= tint <= ifitmax:
                        tdatarange[0].append(tint)
                        ydata.append(tdata)
                if EstDir:
                    thisff = ff.Fitting(Funs=[thisfitfun,npar,'Estimate'],data=[tdatarange,ydata],name=legstate+' '+ifit,iGuess=self.iGuess)
                else:
                    thisff = ff.Fitting(Funs=[thisfitfun,npar],data=[tdatarange,ydata],name=legstate+' '+ifit,iGuess=self.iGuess)
                thistest = thisff.FitBoots()
                if not np.isnan(thistest):

                    if state == 1:
                        thisff.UnExpParams(paramkeys=('Energy',))
                        thisff.ImportFunction( C2OneStateFitFunNoExp,2)
                    elif state == 2:
                        thisff.UnExpParams(paramkeys=('Energy','DE'))
                        thisff.ImportFunction( C2TwoStateFitFunNoExp,4)

                    lFit.append(thisff)
                    lFitA.append(thisff.fit_data['ParamsAvg'])
                    lFitS.append(thisff.fit_data['ParamsStd'])
        indicies = pa.MultiIndex.from_tuples(ilist,names=self.C2_Fit_Col_Names)
        if 'boot' in self.C2_Fit_Stats.columns:
            this_df = pa.DataFrame()
            this_df.loc[:,'boot'] = pa.Series(lFit,index=indicies)
            this_df.loc[:,'Avg'] = pa.Series(lFitA,index=indicies)
            this_df.loc[:,'Std'] = pa.Series(lFitS,index=indicies)
            self.C2_Fit_Stats = self.C2_Fit_Stats.append(this_df)
        else:
            self.C2_Fit_Stats.loc[:,'boot'] = pa.Series(lFit,index=indicies)
            self.C2_Fit_Stats.loc[:,'Avg'] = pa.Series(lFitA,index=indicies)
            self.C2_Fit_Stats.loc[:,'Std'] = pa.Series(lFitS,index=indicies)
        self.Write()



    ## Comparisons

    def __str__(self):
        return self.name


    def __iter__(self):
        return self

    def items(self):
        return list(self.C2_Stats['boot'].items())

    def iteritems(self):
        return iter(self.C2_Stats['boot'].items())

    def itemsAvgStd(self):
        outlist = []
        for ivals in self.C2_Stats.itertuples():
            outlist.append((ivals[0],ivals[2],ivals[3]))
        return outlist


    def values(self):
        return self.C2_Stats['boot'].values

    def valuesAvgStd(self):
        outlist = []
        for ivals in self.C2_Stats.itertuples(index=False):
            outlist.append((ivals[1],ivals[2]))
        return outlist


    def keys(self):
        return self.C2_stats.index


    def __setitem__(self, key, value ):
        self.C2_Stats.loc[key,'boot'] = value

    def __getitem__(self, key):
        return self.C2_Stats.loc[key,'boot']

    def __reversed__(self):
        ## TODO this is incorrect, fix if required
        self.C2_Stats = self.C2_Stats.iiloc[::-1]
        return self


    def __next__(self):
        if self.current >= len(list(self.keys())):
            self.current = 0
            raise StopIteration
        else:
            self.current += 1
            if self.current >= len(list(self.keys()))+1:
                return 0.0
            else:
                thiskey = list(self.keys())[self.current-1]
                return self[thiskey]


    def EffMass(self):
        lEM,lEMA,lEMS = [],[],[]
        for (itflow,it),idata in self.items():
            nextt = (untstr(it) % (self.nt-1)) + 1
            thisEmass = (idata/self.C2_Stats['boot'][itflow,tstr(nextt)]).Log()
            thisEmass.Stats()
            lEM.append(thisEmass)
            lEMA.append(thisEmass.Avg)
            lEMS.append(thisEmass.Std)
        self.C2_Stats.loc[:,'EffM'] = pa.Series(lEM,index=self.C2_Stats.index)
        self.C2_Stats.loc[:,'EffMAvg'] = pa.Series(lEMA,index=self.C2_Stats.index)
        self.C2_Stats.loc[:,'EffMStd'] = pa.Series(lEMS,index=self.C2_Stats.index)

    # def ReadAndBootAndEffM(self,DefWipe=False):
    #     if os.path.isfile(self.PickleFile) and not DefWipe:
    #         self.LoadPickle()
    #     else:
    #         self.ReadAndBoot()
    #         self.EffMass()


    def Write(self):
        def FixDictArray(thislist,ikey):
            outDict = OrderedDict()
            for ic,ilist in enumerate(thislist):
                outDict[ikey+str(ic+1)] = ilist
            return outDict

        if not hasattr(self,'show_cfgs'):
            self.show_cfgs = False

        outDict = ODNested()
        outDict['name'] = self.name
        outDict['kud'] = self.kud
        outDict['ks'] = self.ks
        outDict['nt'] = self.nt
        outDict['nxyz'] = self.nxyz

        for istream,incfg in zip(self.stream_list,self.ncfg_list):
            outDict['ncfg_'+istream] = incfg
        outDict['nboot'] = self.nboot

        excel_params = pa.Series(deepcopy(outDict))

        if 'boot' in self.C2_Fit_Stats:
            for (istate,itflow,ifitr),fitdict in self.C2_Fit_Stats['boot'].items():
                outDict['FlowFits'][istate][itflow][ifitr] = fitdict.GetOutputDict()


        # for ststr,statedict in self.FitDict.iteritems():
        #     if len(statedict.keys()) == 0: continue
        #     for ip,pdict in statedict.iteritems():
        #         for ifit,fitdict in pdict.iteritems():
        #             for ipar,pardict in fitdict.Params.iteritems():
        #                 thisfmtflag = 'f'
        #                 if ipar == 'A_Ep': thisfmtflag = 'e'
        #                 outDict['Fits'][ststr][ip][ifit][ipar] = AvgStdToFormat(pardict.Avg,pardict.Std,frmtflag=thisfmtflag)
        #             outDict['Fits'][ststr][ip][ifit]['Chi^2_pdf'] = AvgStdToFormat(fitdict.Chi2DoF.Avg,fitdict.Chi2DoF.Std)
        #             outDict['Fits'][ststr][ip][ifit]['MaxIter'] = fitdict.MaxIter
        #             outDict['Fits'][ststr][ip][ifit]['Prec'] = fitdict.Prec


        if 'boot' in self.C2_Stats.columns:
            for (itflow,it),tflowdata in self.C2_Stats['boot'].items():
                outDict['C2boot'][itflow][it] = AvgStdToFormat(tflowdata.Avg,tflowdata.Std,frmtflag='e')
        if 'EffM' in self.C2_Stats.columns:
            for (itflow,it),tflowdata in self.C2_Stats['EffM'].items():
                outDict['EffM'][itflow][it] = AvgStdToFormat(tflowdata.Avg,tflowdata.Std,frmtflag='e')

        if self.show_cfgs:
            outDict['cfglist'] = Series_TO_ODict(self.C2_cfgs['configs'])



        ## Write human readable file with data
        WriteXml(self.HumanFile,{'Results':outDict})

        ## TODO, include fitting in the excell file
        WriteExcel(self.ExcelFile,{'C2_results':deepcopy(self.C2_Stats)},cfglist=self.C2_cfgs,params=excel_params)

        outDict['cfglist'] = Series_TO_ODict(self.C2_cfgs['configs'])

        ## pickles rest of data for reading
        WritePickle(self.PickleFile,self.__dict__)
    ##################################

    def ReadAndWrite(self,WipeData=True,NoFit=False):
        self.Read()
        self.Bootstrap(WipeData=WipeData)
        self.EffMass()
        if NoFit:
            self.Fit()
        self.Write()


    def FixRead(self,readdict):
        readdict['Info']['outdir'] = self.Info['outdir']
        readdict['latparams'].outputdir = self.latparams.outputdir
        readdict['latparams'].PickleFile = self.latparams.PickleFile
        readdict['latparams'].HumanFile = self.latparams.HumanFile
        readdict['PickleFile'] = self.PickleFile
        readdict['ExcelFile'] = self.ExcelFile
        readdict['HumanFile'] = self.HumanFile
        return readdict

    def LoadPickle(self,DefWipe=False,WipeData = True,CheckCfgs=False,NoFit=False):
        # print 'filename ' , self.PickleFile
        if os.path.isfile(self.PickleFile) and not DefWipe:
            # print 'Loading Pickle for ' , self.name
            loadeddict = ReadPickleWrap(self.PickleFile)
            loadeddict = self.FixRead(loadeddict)
            checklist = self.thischecklist
            if CheckCfgs: checklist = checklist + ['C2_cfgs']
            if CheckClass(self.__dict__,loadeddict,checklist):
                self.__dict__.update(loadeddict)
                if Wipe_All_Fits:
                    self.C2_Fit_Stats = pa.DataFrame()
                elif len(self.C2_Fit_Stats.values)> 0 and not isinstance(self.C2_Fit_Stats.values[0][0],sff.SetOfFitFuns):
                    print('Legacy file, fits need to be recomputed')
                    self.C2_Fit_Stats = pa.DataFrame()
            else:
                if os.path.isfile(self.PickleFile+'.bak'):
                    print('using backupfile for '+self.name)
                    self.PickleFile = self.PickleFile+'.bak'
                    self.HumanFile = self.HumanFile+'.bak'
                    self.ExcelFile = self.ExcelFile.replace('.xlsx','.bak.xlsx')
                    self.LoadPickle(DefWipe=DefWipe,WipeData=WipeData,CheckCfgs=CheckCfgs)
                    return
                print('Warning, file ' , self.HumanFile , ' has different parameters than this instance, reading and writing over file now')
                print()
                self.ReadAndWrite(WipeData=WipeData,NoFit=NoFit)
                # self.Write()
        else:
            self.ReadAndWrite(WipeData=WipeData,NoFit=NoFit)
    ## Operator overloading

    #Too Much work, only works for one operator. Please use SetCustomName after doing all algebraic manipulation


    def UpdateName(self,opp,C2,C22):
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
        if isinstance(C22,TwoPointCorr) and isinstance(C2,TwoPointCorr):
            op1LL = C2.LegLab[1:-1] ## removes $ symbols, re adds at end
            op2LL = C22.LegLab[1:-1] ## removes $ symbols, re adds at end
            if C2.name == C22.name:
                outname = C2.name,GetOppSelf(opp)
                outLL = op1LL,' Comb ',GetLegSelf(oppLL)
            else:
                for iname,iname2 in zip(C2.name.split('_'),C22.name.split('_')):
                    if iname in iname2:
                        outname.append(iname2)
                    elif iname2 in iname:
                        outname.append(iname)
                    else:
                        outname.append(iname+opp+iname2)
                outname = '_'.join(outname)
                for iLegLab,iLegLab2 in zip(op1LL.split('\ '),op2LL.split('\ ')):
                    if iLegLab in iLegLab2:
                        outLL.append(iLegLab2)
                    elif iLegLab2 in iLegLab:
                        outLL.append(iLegLab)
                    else:
                        outLL.append(iLegLab+oppLL+iLegLab2)
                outLL = '_'.join(outLL)
        elif isinstance(C2,TwoPointCorr):
            op1LL = C2.LegLab[1:-1] ## removes $ symbols, re adds at end
            if isinstance(C22,int):
                outname  =  C2.name,opp+str(C22)
                outLL  =  op1LL,oppLL+str(C22)
            else:
                if C22 > 0.01:
                    outname = C2.name,opp+'{:.2f}'.format(C22)
                    outLL = op1LL,oppLL+'{:.2f}'.format(C22)
        else:
            op2LL = C22.LegLab[1:-1] ## removes $ symbols, re adds at end
            if isinstance(C2,int):
                outname  =  str(C2) + opp ,C22.name
                outLL  =  str(C2) + oppLL,op2LL
            else:
                if C2 > 0.01:
                    outname = '{:.2f}'.format(C2) + opp,C22.name
                    outLL = '{:.2f}'.format(C2) + oppLL,op2LL
        # print r'$'+r'\ '.join(outLL)+'$'
        outLL = list(outLL)
        for ic,iLL in enumerate(outLL):
            if '^{}' in iLL:
                outLL[ic] = iLL.replace('^{}','^{')+'}'
        self.SetCustomName('_'.join(outname),r'$'+r'\ '.join(outLL)+'$')




    ## Numerical operatoions


    def __add__(self,C22):
        if any([ipar in str(type(C22)) for ipar in TwoPointCorr.parentlist]):
            return C22.__radd__(self)
        else:
            result = TwoPointCorr(thisnboot=self.nboot, cfglist=self.C2_cfgs,
                                  Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
            result.SetCustomName(string=self.filename,stringLL=self.LegLab,fo_for_cfgs=self.fo_for_cfgs)
            result.C2_Stats = pa.DataFrame(columns=['boot'],index=self.C2_Stats.index)
            # result.UpdateName('+',self,C22)
            if isinstance(C22,TwoPointCorr):
                result.C2_Stats = self.C2_Stats + C22.C2_Stats.ix[:, C22.C2_Stats.columns != 'Auto']
            else:
                try:
                    result.C2_Stats = self.C2_Stats + C22
                except Exception as err:
                    print(type(C22))
                    raise EnvironmentError('Invalid value to combine with TwoPointCorr class')
            return result

    def __sub__(self,C22):
        if any([ipar in str(type(C22)) for ipar in TwoPointCorr.parentlist]):
            return C22.__rsub__(self)
        else:
            result = TwoPointCorr(thisnboot=self.nboot, cfglist=self.C2_cfgs,
                                  Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
            result.SetCustomName(string=self.filename,stringLL=self.LegLab,fo_for_cfgs=self.fo_for_cfgs)
            result.C2_Stats = pa.DataFrame(columns=['boot'],index=self.C2_Stats.index)
            # result.UpdateName('-',self,C22)
            if isinstance(C22,TwoPointCorr):
                result.C2_Stats = self.C2_Stats - C22.C2_Stats.ix[:, C22.C2_Stats.columns != 'Auto']
            else:
                try:
                    result.C2_Stats = self.C2_Stats - C22
                except Exception as err:
                    print(type(C22))
                    raise EnvironmentError('Invalid value to combine with TwoPointCorr class')
            return result

    def __mul__(self,C22):
        if any([ipar in str(type(C22)) for ipar in TwoPointCorr.parentlist]):
            return C22.__rmul__(self)
        else:
            result = TwoPointCorr(thisnboot=self.nboot, cfglist=self.C2_cfgs,
                                  Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
            result.SetCustomName(string=self.filename,stringLL=self.LegLab,fo_for_cfgs=self.fo_for_cfgs)
            result.C2_Stats = pa.DataFrame(columns=['boot'],index=self.C2_Stats.index)
            # result.UpdateName('x',self,C22)
            if isinstance(C22,TwoPointCorr):
                result.C2_Stats = self.C2_Stats * C22.C2_Stats.ix[:, C22.C2_Stats.columns != 'Auto']
            else:
                try:
                    result.C2_Stats = self.C2_Stats * C22
                except Exception as err:
                    print(type(C22))
                    raise EnvironmentError('Invalid value to combine with TwoPointCorr class')
            return result

    def __div__(self,C22):
        if any([ipar in str(type(C22)) for ipar in TwoPointCorr.parentlist]):
            return C22.__rdiv__(self)
        else:
            result = TwoPointCorr(thisnboot=self.nboot, cfglist=self.C2_cfgs,
                                  Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
            result.SetCustomName(string=self.filename,stringLL=self.LegLab,fo_for_cfgs=self.fo_for_cfgs)
            result.C2_Stats = pa.DataFrame(columns=['boot'],index=self.C2_Stats.index)
            # result.UpdateName('div',self,C22)
            if isinstance(C22,TwoPointCorr):
                result.C2_Stats = self.C2_Stats / C22.C2_Stats.ix[:, C22.C2_Stats.columns != 'Auto']
            else:
                try:
                    result.C2_Stats = self.C2_Stats / C22
                except Exception as err:
                    print(type(C22))
                    print(C22)
                    raise EnvironmentError('Invalid value to combine with TwoPointCorr class')
            return result


    def __pow__(self,C22):
        if any([ipar in str(type(C22)) for ipar in TwoPointCorr.parentlist]):
            return C22.__rpow__(self)
        else:
            result = TwoPointCorr(thisnboot=self.nboot, cfglist=self.C2_cfgs,
                                  Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
            result.SetCustomName(string=self.filename,stringLL=self.LegLab,fo_for_cfgs=self.fo_for_cfgs)
            result.C2_Stats = pa.DataFrame(columns=['boot'],index=self.C2_Stats.index)
            # result.UpdateName('pow',self,C22)
            if isinstance(C22,TwoPointCorr):
                result.C2_Stats = self.C2_Stats ** C22.C2_Stats.ix[:, C22.C2_Stats.columns != 'Auto']
            else:
                try:
                    result.C2_Stats = self.C2_Stats ** C22
                except Exception as err:
                    print(type(C22))
                    print(C22)
                    raise EnvironmentError('Invalid value to combine with TwoPointCorr class')
            return result


    ## Right multiplication functions


    def __radd__(self,C22):
        if isinstance(C22,TwoPointCorr):
            return C22 + self
        result = TwoPointCorr(thisnboot=self.nboot, cfglist=self.C2_cfgs,
                              Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.SetCustomName(string=self.filename,stringLL=self.LegLab,fo_for_cfgs=self.fo_for_cfgs)
        result.C2_Stats = pa.DataFrame(columns=['boot'],index=self.C2_Stats.index)
        # result.UpdateName('+',C22,self)
        try:
            result.C2_Stats =  C22 + self.C2_Stats
        except Exception as err:
            print(type(C22))
            print(C22)
            raise EnvironmentError('Invalid value to combine with TwoPointCorr class')
        return result

    def __rsub__(self,C22):
        if isinstance(C22,TwoPointCorr):
            return C22 - self
        result = TwoPointCorr(thisnboot=self.nboot, cfglist=self.C2_cfgs,
                              Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.SetCustomName(string=self.filename,stringLL=self.LegLab,fo_for_cfgs=self.fo_for_cfgs)
        result.C2_Stats = pa.DataFrame(columns=['boot'],index=self.C2_Stats.index)
        # result.UpdateName('+',C22,self)
        try:
            result.C2_Stats =  C22 - self.C2_Stats
        except Exception as err:
            print(type(C22))
            print(C22)
            raise EnvironmentError('Invalid value to combine with TwoPointCorr class')
        return result

    def __rmul__(self,C22):
        if isinstance(C22,TwoPointCorr):
            return C22 * self
        result = TwoPointCorr(thisnboot=self.nboot, cfglist=self.C2_cfgs,
                              Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.SetCustomName(string=self.filename,stringLL=self.LegLab,fo_for_cfgs=self.fo_for_cfgs)
        result.C2_Stats = pa.DataFrame(columns=['boot'],index=self.C2_Stats.index)
        # result.UpdateName('+',C22,self)
        try:
            result.C2_Stats =  C22 * self.C2_Stats
        except Exception as err:
            print(type(C22))
            print(C22)
            raise EnvironmentError('Invalid value to combine with TwoPointCorr class')
        return result

    def __rdiv__(self,C22):
        if isinstance(C22,TwoPointCorr):
            return C22 / self
        result = TwoPointCorr(thisnboot=self.nboot, cfglist=self.C2_cfgs,
                              Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.SetCustomName(string=self.filename,stringLL=self.LegLab,fo_for_cfgs=self.fo_for_cfgs)
        result.C2_Stats = pa.DataFrame(columns=['boot'],index=self.C2_Stats.index)
        # result.UpdateName('+',C22,self)
        try:
            result.C2_Stats =  C22 / self.C2_Stats
        except Exception as err:
            print(type(C22))
            print(C22)
            raise EnvironmentError('Invalid value to combine with TwoPointCorr class')
        return result


    def __rpow__(self,C22):
        if isinstance(C22,TwoPointCorr):
            return C22 ** self
        result = TwoPointCorr(thisnboot=self.nboot, cfglist=self.C2_cfgs,
                              Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.SetCustomName(string=self.filename,stringLL=self.LegLab,fo_for_cfgs=self.fo_for_cfgs)
        result.C2_Stats = pa.DataFrame(columns=['boot'],index=self.C2_Stats.index)
        # result.UpdateName('+',C22,self)
        try:
            result.C2_Stats =  C22 ** self.C2_Stats
        except Exception as err:
            print(type(C22))
            print(C22)
            raise EnvironmentError('Invalid value to combine with TwoPointCorr class')
        return result


    ## unitary arithmetic operations

    def __neg__(self):
        for ival in list(self.keys()):
            self[ival] = -self[ival]
        return self


    def __pos__(self):
        return self

    def __abs__(self):
        for ival in list(self.keys()):
            self[ival] = abs(self[ival])
        return self

    ## Casting, just casts the bootstrap values, ignore return value (only needed since return has to be same type).


    def __complex__(self):
        for ival in list(self.keys()):
            complex(self[ival])
        self.Stats()
        return 1+0j

    def __int__(self):
        for ival in list(self.keys()):
            int(self[ival])
        self.Stats()
        return 1

    def __float__(self):
        for ival in list(self.keys()):
            np.float64(self[ival])
        self.Stats()
        return 1.0







def TestTPCorr(DefWipe=False):
    from Params import defInfo
    data = TwoPointCorr(Info=defInfo)
    data.LoadPickle(DefWipe=DefWipe)

    this_info = pa.Series()
    this_info['save_file'] = this_dir+'/TestGraphs/TestFit2.pdf'
    this_info['title'] = 'Test Two Point Corr'
    import PlotData as jpl
    data_plot = jpl.Plotting(plot_info=this_info)
    data_plot = data.EffMassPlot(data_plot,xlims=(2,20),momlist=['p000'],thiscol='Blue',thissym='o',thisshift=0.1)
    data_plot = data.EffMassFitPlot(data_plot,1,'fitr9-20',momlist=['p000'])
    data_plot = data.EffMassFitPlot(data_plot,2,'fitr5-20',momlist=['p000'])
    data_plot.PrintData()
    data_plot.PlotAll()
    data.Write()
    return data

def TestNNQCorr(DefWipe=False):
    this_info = pa.Series()
    this_info['save_file'] = this_dir+'/TestGraphs/TestNNQ.pdf'
    this_info['title'] = 'Test Alpha'
    from Params import defInfoFlow
    dataflow = NNQCorr(Info=defInfoFlow)
    dataflow.LoadPickle(DefWipe=DefWipe)
    import PlotData as jpl
    data_plot = jpl.Plotting(plot_info=this_info)
    data_plot = dataflow.AlphaPlot(data_plot,xlims=[2,10],
                                   thiscol='Blue',momlist=['p000'],
                                   thissym='o',thisshift=0.0,thistflowlist = ['t_f6.01'])
    data_plot = dataflow.AlphaFitPlot(data_plot,'fitr8-15',
                                      momlist=['p000'],thistflowlist = ['t_f6.01'])
    data_plot.PrintData()
    data_plot.PlotAll()
    dataflow.Write()
    return dataflow

def TestNNQFullCorr(DefWipe=False):
    this_info = pa.Series()
    this_info['save_file'] = this_dir+'/TestGraphs/TestNNQFull.pdf'
    this_info['title'] = 'Test Alpha Full'
    from Params import defInfoFlow
    defInfoFlow['pmom'] = ['0 0 0']
    defInfoFlow['Interp'] = 'CPOdd'
    defInfoFlow['tflowlist'] = np.array([4.01,5.01,6.01])
    defInfoFlow['tsum_fit_list'] = ['ts28','ts30','ts32']
    defInfoFlow['tsum_list'] = ['ts28','ts30','ts32']
    defInfoFlow['nrandt'] = 1
    # defInfoFlow['tsum_fit_list'] = ['ts'+str(it) for it in range(1,11)]
    # defInfoFlow['tsum_list'] = ['ts'+str(it) for it in range(1,11)]
    dataflow = NNQFullCorr(Info=defInfoFlow)
    dataflow.LoadPickle(DefWipe=True)
    import PlotData as jpl
    data_plot = jpl.Plotting(plot_info=this_info)
    data_plot = dataflow.AlphaPlot( data_plot,xlims=[0,20],thiscol='Blue',tsum_list=['ts10'],
                                    momlist=['p000'],thissym='o',thisshift=0.0,thistflowlist = ['t_f6.01'])
    data_plot = dataflow.AlphaFitPlot(data_plot,['fitr8-15','tsumfitr0-10'],momlist=['p000'],thistflowlist = ['t_f6.01'])
    data_plot.PrintData()
    data_plot.PlotAll()
    dataflow.Write()
    # dataflow.Fit_All_Tau([['ts1','ts10']])
    return dataflow

def ReformatC2(this_file,out_filename):
    if 'CFGNUMB' not in out_filename:
        raise IOError('CFGNUMB not in out_filename.')
    meta,data = ReadWithMeta(this_file)
    data = data.reset_index()['C2'].values
    data = np.array([np.array(idata) for idata in data])
    data = np.swapaxes(data,0,1)
    for imom,momdata in enumerate(data):
        mom_file = out_filename.replace('MOMNUMB','p_'+meta[0][imom].replace(' ','_'))
        try:
            os.mkdir(mom_file.replace('Nucleon_CFGNUMB.txt',''))
        except Exception as err:
            pass
        for ic,idata in enumerate(momdata):
            np.savetxt(mom_file.replace('CFGNUMB',str(ic).zfill(4)),idata)
    return data

def TestFlowCorr(DefWipe=False):
    from Params import defInfoFlow
    import PlotData as jpl
    defInfoFlow['tflowlist'] = np.arange(0,10,0.4)
    defInfoFlow['Adjoint_Flow'] = False
    this_info = pa.Series()
    this_info['save_file'] = this_dir+'/TestGraphs/TestFlowMes.pdf'
    this_info['title'] = 'Test Flowed Fermions'
    this_info['x_label'] = 'source_sink_separation (lattice units)'
    this_info['y_label'] = 'effective mass'
    # this_info['ylims'] = [-5,5]
    data = FlowedTwoPtCorr(Info=defInfoFlow)
    data.LoadPickle(DefWipe=True)
    data_plot = jpl.Plotting(plot_info=this_info)
    this_xrange = (0,63)
    data_plot = data.EffMassPlot(data_plot,xlims=this_xrange,tflowlist=['t_f0.00'],thiscol='Black',thissym='o',thisshift=0)
    data_plot = data.EffMassPlot(data_plot,xlims=this_xrange,tflowlist=['t_f0.40'],thiscol='Blue',thissym='o',thisshift=0)
    # data_plot = data.EffMassFitPlot(data_plot,1,'fitr20-25',tflowlist=['t_f0.00'])
    data_plot = data.EffMassPlot(data_plot,xlims=this_xrange,tflowlist=['t_f2.00'],thiscol='Red',thissym='o',thisshift=0.001)
    # data_plot = data.EffMassFitPlot(data_plot,1,'fitr20-25',tflowlist=['t_f2.00'])
    data_plot = data.EffMassPlot(data_plot,xlims=this_xrange,tflowlist=['t_f6.00'],thiscol='Green',thissym='o',thisshift=-0.001)
    # data_plot = data.EffMassFitPlot(data_plot,1,'fitr20-25',tflowlist=['t_f6.00'])
    data_plot = data.EffMassPlot(data_plot,xlims=this_xrange,tflowlist=['t_f10.00'],thiscol='Orange',thissym='o',thisshift=0.002)
    # data_plot = data.EffMassFitPlot(data_plot,1,'fitr20-25',tflowlist=['t_f10.00'])
    data_plot.PrintData()
    data_plot.PlotAll()


    this_info = pa.Series()
    this_info['save_file'] = this_dir+'/TestGraphs/TestFlowRatMes.pdf'
    this_info['title'] = 'Ratio of flowed to non-flowed correlators'
    this_info['x_label'] = 'source_sink_separation (lattice units)'
    this_info['y_label'] = 'Ratio'
    this_info['xlims'] = [0,20]
    data_plot = jpl.Plotting(plot_info=this_info)
    this_xrange = (0,63)
    data = FlowedTwoPtCorr(Info=defInfoFlow)
    data.LoadPickle(DefWipe=DefWipe)
    data_plot = data.FlowRatPlot(data_plot,xlims=this_xrange,tflowlist=['t_f0.00'],thiscol='Black',thissym='o',thisshift=0)
    data_plot = data.FlowRatPlot(data_plot,xlims=this_xrange,tflowlist=['t_f0.40'],thiscol='Blue',thissym='o',thisshift=0)
    data_plot = data.FlowRatPlot(data_plot,xlims=this_xrange,tflowlist=['t_f2.00'],thiscol='Red',thissym='o',thisshift=0.001)
    data_plot = data.FlowRatPlot(data_plot,xlims=this_xrange,tflowlist=['t_f6.00'],thiscol='Green',thissym='o',thisshift=-0.001)
    data_plot = data.FlowRatPlot(data_plot,xlims=this_xrange,tflowlist=['t_f10.00'],thiscol='Orange',thissym='o',thisshift=0.002)
    data.Write()
    data_plot.PrintData()
    data_plot.PlotAll()

    this_info = pa.Series()
    this_info['save_file'] = this_dir+'/TestGraphs/AllFlowRatMes.pdf'
    this_info['title'] = 'Ratio of flowed to non-flowed correlators'
    this_info['x_label'] = 'flow time $\sqrt{8t_f}[fm]'
    this_info['y_label'] = 'Ratio'
    data_plot = jpl.Plotting(plot_info=this_info)
    data = FlowedTwoPtCorr(Info=defInfoFlow)
    data.LoadPickle(DefWipe=DefWipe
    )
    data.AllFlowRatPlot(data_plot,tsink='t28',thiscol='Blue',thissym='o',thisshift=0)
    data_plot.PrintData()
    data_plot.PlotAll()
    return data

def ReadBootFile(this_file):
    from FileIO import ReadXml
    data = ReadXml(this_file)
    out_data = []
    for imom,mom_data in data['Results'].items():
        for it,t_data in mom_data.items():
            out_data.append(t_data)
    return np.array(out_data).reshape(len(data['Results'].keys()),-1),list(data['Results'].keys())

def AvgBootData(this_path):
    data = []
    for ifile in glob.glob(this_path+'/*'):
        this_data,mom_list = ReadBootFile(ifile)
        data.append(this_data)
    return np.mean(np.array(data),axis=0),np.std(np.array(data),axis=0),mom_list

def AvgBootDir(this_path):
    for ifile in glob.glob(this_path+'/*'):
        if os.path.isdir(ifile):
            print('Formatting:' ,ifile)
            out_file_Avg = ifile+'_Avg.txt'
            out_file_Err = ifile+'_Err.txt'
            this_data,this_data_err,mom_list = AvgBootData(ifile)
            for imom,mom_data,mom_data_err in zip(mom_list,this_data,this_data_err):
                np.savetxt(out_file_Avg.replace('.txt','_'+imom+'.txt'),mom_data)
                np.savetxt(out_file_Err.replace('.txt','_'+imom+'.txt'),mom_data_err)

def TestLoadFile(this_file):
    return TwoPointCorr.Construct_2pt_File(this_file)

if __name__ == '__main__':
    tpdata = TestTPCorr()
    nnqdata = TestNNQCorr()
    # nnqfull = TestNNQFullCorr()
    # flowdata = TestFlowCorr()
    # this_file = '/home/jackdra/LQCD/Results/DebugResults/RC32x64Kud01375400Ks01364000/G2/Pickle/MassForFF.py3p'
    # test_object = TestLoadFile(this_file)
    print('Testing 2ptcorr complete, see tpdata and nnqdata variables')
