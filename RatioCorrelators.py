#!/usr/bin/env python


#IMPORT THIS FIRST, put in base import file
# import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
# import matplotlib.pyplot as pl
##
from Params import nboot,outputdir,Debug_Mode,Wipe_All_Fits,this_dir
from FileIO import ReadFuns,WriteFuns
# from PredefFitFuns import ConstantFitFun,RF_Prop_Fit
from PredefFitFuns import ConstantFitFun
import SetsOfFits as sff
from NullPlotData import null_series
from TwoPtCorrelators import NNQCorr
# from TwoPtCorrelators import TwoPointCorr,NNQCorr,NNQFullCorr
# from FlowOpps import FlowOp

from MiscFuns import ODNested,mkdir_p,CheckClass,CombineCfgs,op_str_list,MultiSplit,op_dict
from MiscFuns import flip_fun_2arg,GenSinkFun,CreateNameAndLL,GetInfoList,Series_fix_key
from XmlFormatting import untstr,tstr,unxmlfitr,AvgStdToFormat,AvgFormat
from XmlFormatting import rm_tsumfitr,unxmlfitr_int,KeyForamtting
# from XmlFormatting import rm_tsumfitr,tflowTOLeg,tsumTOLeg,tsinkTOLeg
from FileIO import WriteXml,WritePickle,ReadPickleWrap
from FileIO import WriteExcel,Construct_File_Object
from copy import deepcopy
# import cPickle as pickle
# import dill as pickle
import pandas as pa
# import warnings
# import itertools

import operator as op
import os
import numpy as np
from SetsOfTwoPtCorrs import SetOfTwoPt,SetOfNNQ
from ThreePtCorrelators import ThreePointCorr,NJNQCorr,NJNQFullCorr
# import FitFunctions as ff

from collections import OrderedDict

# from FlowOpps import FlowOp




"""
How to use:

Ratio(cfglist=[], Info={})
creates class using:
 configuration list cfglist
 Information Dictionary Info (See below for what can be in Info)


.LoadPickle()
Does everything basically.

"""

## function to center the x-axis around zero
def CenterT(tlist):
    return tlist - np.mean(tlist),-tlist[0]+np.mean(tlist)+1

## for some reason, CROMA ouptut of the three-point correlator is shifted by 3 in time compared to the two-point correlator
twopttshift = 1


class RatioCorr(object):
    """

    Rat(q,tau,pp,t) uses bootstrap class
    pp  = sink momentum
    t   = source sink separation time
    q   = current insersion momentum
    tau = current insersion time

    Info is a dictionary containing information for the correlator:
    see comments below to see what is required.
    missing elements just sets field to default


    """
    def Construct_Rat_File(this_file):
        return Construct_File_Object(this_file,RatioCorr)


    ## TODO: create variational flag and interface with this function.

    ## array of all classes that inherit this class. Append this when adding new classes please!
    parentlist = ['SetOfCorrs.SetsOfRatios','RatioCorrelators.RatioFOCorr','RatioCorrelators.RatioFOFullCorr']

    def __init__(self, thisnboot=nboot, cfglist={}, Info={},name='',
                 filename='',thissym='Not Set',thiscol='Not Set',thisshift=0.0,man_load_cfgs=False):


        self.current = 0
        ## Initialising plotting parameters
        self.thisshift = thisshift
        self.thiscol = thiscol
        self.thissym = thissym
        ## the actual shift ammount needs to be scaled by the x axis size
        ## see GetShift
        self.thisshiftscale = 'Not Set'



        self.thischecklist = ['C3file','Sym']


        self.C3class = ThreePointCorr(thisnboot=thisnboot, cfglist=cfglist, Info=Info,name=name,filename=filename)
        self.C3file = self.C3class.PickleFile
        self.checkmomlist = self.C3class.checkmomlist
        self.latparams = self.C3class.latparams


        if 'pmom' not in list(Info.keys()):
            Info['pmom'] = [self.latparams.GetAvgmomstr(ip) for ip in self.C3class.qmom]
        else:
            if isinstance(Info['pmom'],str):
                if Info['pmom'] == 'All':
                    Info['pmom'] = self.latparams.GetMomList(list_type='str_Avg',list_order='squared')
                else:
                    Info['pmom'] = [Info['pmom']]
            Info['pmom'] = list(Info['pmom'])
            for ip in self.C3class.qmom:
                fmt_p = self.latparams.GetAvgmomstr(ip)
                if fmt_p not in Info['pmom']:
                    Info['pmom'].append(fmt_p)
            Info['pmom'] = self.latparams.MomSqrdSort(Info['pmom'])

        if '0 0 0' not in Info['pmom']:
            Info['pmom'] = ['0 0 0'] + Info['pmom']
        if 'AllC2p' in list(Info.keys()) and Info['AllC2p']:
            Info['pmom'] = self.latparams.GetMomList(list_type='str_Avg',list_order='squared')


        self.nt = self.latparams.nt
        self.nxyz = self.latparams.nxyz
        self.dim_label = 'RC'+str(self.nxyz)+'x'+str(self.nt)
        self.dim_ll = str(self.nxyz)+r'^{3}\times ' + str(self.nt)

        # Info['pmom'] = self.latparams.momstrlistAvg

        ## both current and sink momentum are requred for the two point correlator construction when constructing the ratio function
        ## but the two point correlator has been averaged over, so we need to get the averaged momentum anaolgy
        ## bah, this was too hard to keep unique files for the two-point correlator.
        ## easier just to read and write the whole 2 point mom file.
        ## e.g. q = [-1,1,-1] corresponds to p = [1 1 1] in the two-point correlator
        # ppforp = self.latparams.GetAvgmomstr(Info['ppmom'])
        # if ppforp not in Info['pmom']: Info['pmom'].append(ppforp)

        if 'show_cfgs' in list(Info.keys()): self.show_cfgs = Info['show_cfgs']
        else: self.show_cfgs = True

        ## used to pick out different configuraiton/source lists.
        ## for example, if there are less flowed operators than 2 point correlators, passing 'TopCharge' will read that result
        if 'fo_for_cfgs' in list(Info.keys()): self.fo_for_cfgs = Info['fo_for_cfgs']
        else: self.fo_for_cfgs = False


        ## Does symmetrising over smearings
        if 'Symmetric' in list(Info.keys()): self.Sym = Info['Symmetric']
        else: self.Sym = False ## default to no symmetrizing, (usually full matrix isn't computed nowerdays.

        self.sinklab,self.jsmcoeffs,self.sm2ptList,InfoList = GetInfoList(deepcopy(Info))

        # self.C2filelist = [ikey.PickleFile for ikey in self.C2classlist.itervalues()]
        self.C2Set = SetOfTwoPt(cfglist=cfglist,InfoDict = InfoList)
        self.C2class = 'Not Set'

        self.defgammamomlist = OrderedDict()
        for igamma in self.C3class.gammalist:
            self.defgammamomlist[igamma] = self.C3class.qform

        self.Rat_Col_Names = ['current_operator','momentum','source_sink_separation']
        self.Rat_Stats = pa.DataFrame()

        '''
        self.Rat_Stats DataFrame:

        columns =   boot,   Avg,        Std

        rows =      current_operator,   momenta,    time separation
        '''

        self.Rat_Fit_Col_Names = ['current_operator','momentum']
        self.Rat_Fit_Stats = pa.DataFrame()

        '''
        self.Rat_Fit_Stats DataFrame:

        columns =   boot

        rows =      current_operator, momenta

        Avg/Std elements are pa.Series instances with:
            index = fit parameters
            values = Avg/Std values
        boot elements are SetOfFitFuns objects
        '''

        # # Ratboot [ igamma , ip , it ] bs
        # self.Ratboot = ODNested()
        # self.RatAvg  = ODNested()
        # self.RatStd  = ODNested()
        #
        # ## self.FitDict = { igamma , ip , ifit } = thisff
        # self.FitDict    = ODNested()
        # self.FitDictAvg = ODNested()
        # self.FitDictStd = ODNested()

        self.nboot = thisnboot
        self.centert = 'Not Set'
        self.cfglist = cfglist


        ## Info can contain fit ranges to be calculated
        ## must bein form Info['FitsRF'] = [ 'fit#-#' , ...]
        if 'FitsRF' in list(Info.keys()): self.PredefFits = Info['FitsRF']
        else: self.PredefFits = []

        if 'RFFitFun' in list(Info.keys()): self.Fun,self.npar = Info['RFFitFun']
        else: self.Fun,self.npar = (ConstantFitFun,1)

        ## Info can contain initial guess for fits, make sure it array of correct length
        if 'iGuessRF' in list(Info.keys()): self.iGuess = np.array(Info['iGuessRF'])
        else: self.iGuess = 'None' ## Defaults to value set for function in FitFunctions.py

        ## is needed?
        self.Info = Info
        self.SqrtFact = 'Not Set'

        self.SetCustomName(string=name,stringfn=filename)
        if not man_load_cfgs:
            self.LoadCfgs(cfglist,name,filename)

    def LoadCfgs(self,cfglist,name='',filename=''):
        if len(cfglist) == 0:
            self.CombineCfgLists()
        else:
            self.ImportCfgList(cfglist)
        self.SetCustomName(string=name,stringfn=filename)


    ## just for consistancy with other classes
    def GetCfgList(self):
        self.CombineCfgLists()

    def CombineCfgLists(self):
        if isinstance(self.C2class,str):
            thiscfglist = CombineCfgs(self.C3class.C3_cfgs,self.C2Set.cfglist)
        else:
            thiscfglist = CombineCfgs(self.C3class.C3_cfgs,self.C2class.C2_cfgs)
        self.ImportCfgList(thiscfglist)





    def SetCustomName(self,string='',stringfn='',stringLL=''):
        if string == '':
            self.name = 'Rat_'+self.C3class.name
        else:
            self.name = string
        if stringfn == '':
            self.filename = 'Rat_'+self.C3class.filename
        else:
            self.filename = stringfn
        if stringLL == '':
            # self.LegLab = '$'+'\ '.join([self.C3class.DS,self.C3class.ism,self.C3class.jsm])+'$' ## customise this how you like
            # self.LegLab = '$'+'\ '.join([self.latparams.GetPionMassLab(),self.C3class.DS,self.C3class.ism,self.C3class.jsm,self.C3class.Proj,self.C3class.gammalist[0]])+'$' ## customise this how you like
            self.LegLab = '$'+'\ '.join([   self.dim_ll,self.latparams.GetPionMassLab(),self.C3class.stream_name,
                                            self.C3class.DS,self.C3class.Proj])+'$' ## customise this how you like
        else:
            self.LegLab = stringLL
        self.SubDir = '_'.join([self.C3class.Proj,self.C3class.DS,self.C3class.ppform])
        if isinstance(self.fo_for_cfgs,str):
            self.Ratdir = outputdir + '/'+self.dim_label+self.C3class.kappafolder+'/Rat/'+self.fo_for_cfgs+'/'+self.SubDir+'/'
        else:
            self.Ratdir = outputdir + '/'+self.dim_label+self.C3class.kappafolder+'/Rat/'+self.SubDir+'/'
        mkdir_p(self.Ratdir+'/Pickle/')
        mkdir_p(self.Ratdir+'/Excel/')
        self.HumanFile = self.Ratdir+self.filename+'.xml'
        self.ExcelFile = self.Ratdir+'/Excel/'+self.filename+'.xlsx'
        self.PickleFile = self.Ratdir+'/Pickle/'+self.filename+'.py3p'


    def ImportCfgList(self,cfglist):
        self.cfglist = cfglist
        self.C2Set.ImportCfgList(cfglist)
        if self.C2class != 'Not Set': self.C2class.ImportCfgList(cfglist)
        self.C3class.ImportCfgList(cfglist)
        self.stream_name = self.C3class.stream_name


    def Stats(self):
        if hasattr(self.C2class,'Stats'):
            self.C2class.Stats()
        if hasattr(self.C3class,'Stats'):
            self.C3class.Stats()
        if 'boot' not in self.Rat_Stats: return
        lAvg,lStd = [],[]
        for (igamma,ip,it),iRat in self.items():
            iRat.Stats()
            lAvg.append(iRat.Avg)
            lStd.append(iRat.Std)
        indicies = pa.MultiIndex.from_tuples(self.Rat_Stats.index,names=self.Rat_Col_Names)
        self.Rat_Stats.loc[:,'Avg'] = pa.Series(lAvg,index=indicies)
        self.Rat_Stats.loc[:,'Std'] = pa.Series(lStd,index=indicies)

    def CombineWithAlpha(self,alpha_in,flow_time,this_fun=op.truediv):
        tsink_str = self.C3class.tsink.replace('tsink','t')
        if isinstance(this_fun,str):
            this_fun = op_dict[this_fun.replace('alpha_rat_','').replace('_corr','')]
        if not isinstance(alpha_in,NNQCorr):
            print(type(alpha_in))
            print(alpha_in)
            raise AttributeError('Please pass in NNQCorr class instance')
        if 'boot' not in self.Rat_Stats:
            raise AttributeError('please compute Rat_Stats before combining with alpha')
        if 'Alphaboot' not in alpha_in.NNQ_Stats :
            raise AttributeError('please compute NNQ_Stats in alpha before combining with alpha')
        if ('p000',tsink_str,flow_time) not in alpha_in.NNQ_Stats['Alphaboot']:
            print(alpha_in.NNQ_Stats['Alphaboot'])
            raise AttributeError('alpha ratio needs p000,'+tsink_str+','+flow_time+' result')
        alpha_p0 = alpha_in.NNQ_Stats.loc[('p000',tsink_str,flow_time),'Alphaboot']
        vall,ilist,corrl = [],[],[]

        for (igamma,ip,it),irat in self.Rat_Stats['boot'].items():
            ilist.append((igamma,ip,flow_time,it))
            output = this_fun(irat,alpha_p0[('p000',tsink_str,flow_time)])
            output.Stats()
            vall.append(output)
            corrl.append(irat.corr(alpha_p0[('p000',tsink_str,flow_time)]))
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.Rat_Col_Names)
            self.Rat_Stats.loc[:,'alpha_rat_'+this_fun.__name__] = pa.Series(vall,index=indicies)
            self.Rat_Stats.loc[:,'alpha_rat_'+this_fun.__name__+'_corr'] = pa.Series(corrl,index=indicies)
            self.Write()


    def RemoveVals(self):
        # for iC2 in self.C2classlist.itervalues():
        #     iC2.RemoveVals()
        self.C2class.RemoveVals()
        self.C3class.RemoveVals()

    ## self.Rat = [ igamma , ip, it , iconf ]
    ## Configs must be the same, so we must check for them.
    def Read(self,WipeData=True,DefWipe=False,CheckMom=True):
        if isinstance(self.C2class,str): self.C2Set.LoadPickle(DefWipe=DefWipe,CheckCfgs=True,CheckMom=CheckMom,NoFit=True)
        self.C3class.LoadPickle(WipeData=WipeData,DefWipe=DefWipe,CheckCfgs=True,CheckMom=CheckMom)

    def CreateC2Sink(self):
        if isinstance(self.C2class,str):
            thisname,thisleglab,jsmlab = CreateNameAndLL(self.C2Set.SetC2,self.sinklab)
            thisfun = GenSinkFun(self.jsmcoeffs)
            self.C2class = self.C2Set.CombineCorrs(thisfun,filename=thisname,LegLab=thisleglab,jsmlab=jsmlab)
            self.C2Set = 'Complete'

    def Bootstrap(self,WipeData=True,DefWipe=False):
        if isinstance(self.C2class,str) :
            self.Read(WipeData=WipeData,DefWipe=DefWipe)
            self.CreateC2Sink()
        # self.C2class.Bootstrap(WipeData=WipeData)
        self.C3class.Bootstrap(WipeData=WipeData,DefWipe=DefWipe)


    def CalcRatio(self):
        self.CalcSqrtFac()
        ## if DS has a divide in it, we ignore the sqrt factor as they cancel
        div_bool = 'divide' in self.C3class.DS
        ## double the square root factor if we multiply two ds types
        mult_bool = 'multiply' in self.C3class.DS
        lRat,lRatAvg,lRatStd = [],[],[]
        ilist = []
        tsinkint = int(self.C3class.tsink.replace('tsink',''))
        for (igamma,ip,it),tC3 in self.C3class.C3_Stats['boot'].items():
            ip2pt = self.latparams.GetAvgmomform(ip,pref='p')
            ict = int(untstr(it))
            if ict > tsinkint: continue
            if (ip2pt,it) not in self.SqrtFact.index: continue
            # warnings.warn('DEBUG, please check if tsink for 2-pt corr is correct via g4.')
            # print igamma, ip, it, self.C3class.tsink
            # print 'SqrtFact'
            # print self.SqrtFact
            # print list(self.SqrtFact.index), ip2pt,it
            # print type(self.SqrtFact[ip2pt,it])
            # print self.SqrtFact[ip2pt,it]
            # print self.SqrtFact[ip2pt,it].Avg

            if div_bool:
                this_boot = tC3
            elif mult_bool:
                this_boot = tC3*(self.SqrtFact[ip2pt,it]**2)
            else:
                this_boot = tC3*self.SqrtFact[ip2pt,it]
            this_boot.Stats()
            ilist.append((igamma,ip,it))
            lRat.append(this_boot)
            lRatAvg.append(this_boot.Avg)
            lRatStd.append(this_boot.Std)
            # self.Ratboot[igamma][ip][it] = tC3/self.C2class[ip.replace('q','p'),'t'+str(int(self.C3class.tsink.replace('tsink',''))+1)]
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.Rat_Col_Names)
            self.Rat_Stats.loc[:,'boot'] = pa.Series(lRat,index=indicies)
            self.Rat_Stats.loc[:,'Avg'] = pa.Series(lRatAvg,index=indicies)
            self.Rat_Stats.loc[:,'Std'] = pa.Series(lRatStd,index=indicies)



    def CalcSqrtFac(self):
        SqrtFact_hold = []
        ilist = []
        ipp = self.latparams.GetAvgmomform(self.C3class.ppform,pref='p')
        if Debug_Mode:
            print()
            print('Outputting sqrt factor debug info:')
            print()
        for iq in self.C3class.qform:
            ip = self.latparams.GetAvgmomform(iq,pref='p')
            sqrt_fac = self.CalcSqrtFacMom(self.C2class.C2_Stats['boot'][ipp],self.C2class.C2_Stats['boot'][ip],self.C3class.tsink)
            for it,tdata in enumerate(sqrt_fac):
                if (ip,tstr(it+1)) not in ilist:
                    if Debug_Mode:
                        for ierr in tdata.GetNaNIndex():
                            print(iq, tstr(it+1), ' iboot'+str(ierr).zfill(4))
                    SqrtFact_hold.append(tdata)
                    ilist.append((ip,tstr(it+1)))
            # if neg_sqrt: warnings.warn('negative square root in square root factor for '+ipp+iq)
        indicies = pa.MultiIndex.from_tuples(ilist,names=['momentum','source_sink_separation'])
        self.SqrtFact = pa.Series(SqrtFact_hold,index=indicies)

    def CheckNegForSqrt(self,thisboot,negsqrt):
        for ict in range(len(thisboot.bootvals)):
            if thisboot.bootvals[ict] < 0.0:
                negsqrt = True
                # if G2_abs:
                #     thisboot.bootvals[ict] = abs(thisboot.bootvals[ict])
                # else:
                thisboot.bootvals[ict] = float('NaN')
                # thisboot.bootvals[ict] = 0.0
        return thisboot,negsqrt


    ## taking squareroot of each element before combinding makes the numbers not too small
    def CalcSqrtFacMom(self,data2ptpp,data2ptp,thistsink):
        dictout = []
        thistsink = thistsink.replace('tsink','t')
        thistsinkp1 = 't'+str(int(thistsink.replace('t',''))+twopttshift)
        inttsink = int(thistsink.replace('t',''))
        # NegSQRTRatFac = False
        # data2ptpp[thistsinkp1],NegSQRTRatFac = self.CheckNegForSqrt(data2ptpp[thistsinkp1],NegSQRTRatFac)
        # data2ptp[thistsinkp1],NegSQRTRatFac = self.CheckNegForSqrt(data2ptp[thistsinkp1],NegSQRTRatFac)

        data2ptpp = data2ptpp.apply(lambda x : x.Sqrt())
        data2ptp = data2ptp.apply(lambda x : x.Sqrt())
        # NegSQRTRatFac = any(data2ptpp.isnull().values) or any(data2ptp.isnull().values)
        two = data2ptpp[thistsinkp1]*data2ptp[thistsinkp1]
        # two = two**(0.5)
        thisshift = twopttshift
        for it in range(thisshift,inttsink+thisshift):
            thiststr = tstr(it)
            tflip = tstr(inttsink-it+twopttshift+thisshift)
            # data2ptpp[thiststr],NegSQRTRatFac = self.CheckNegForSqrt(data2ptpp[thiststr],NegSQRTRatFac)
            # data2ptp[thiststr],NegSQRTRatFac = self.CheckNegForSqrt(data2ptp[thiststr],NegSQRTRatFac)
            # data2ptpp[tflip],NegSQRTRatFac = self.CheckNegForSqrt(data2ptpp[tflip],NegSQRTRatFac)
            # data2ptp[tflip],NegSQRTRatFac = self.CheckNegForSqrt(data2ptp[tflip],NegSQRTRatFac)
            one = data2ptpp[thiststr]/data2ptp[thiststr]
            three = data2ptp[tflip]/data2ptpp[tflip]
            inter = one*three
            ott = inter/two
            # ott.Sqrt()
            ott.Stats()
            dictout.append(ott)
        return dictout


    def ImportFitRanges(self,fit_range):
        if 'PreDef' == fit_range:
            self.fit_range = self.PredefFits
        else:
            self.fit_range = fit_range
        if not isinstance(self.fit_range,str):
            print('SetsOfFits implementation now accepts 1 fit range to scan between. Selecting first one')
            self.fit_range = self.fit_range[0]
        self.fit_range = rm_tsumfitr(self.fit_range)


    ## NB, iGuess here is not implemented yet, needs to be passed into the function
    def Fit(self,fit_range='PreDef',iGuess='PreDef',EstDir=False,WipeFit=True,show_timer=False,min_fitr=2):
        self.ImportFitRanges(fit_range)
        if iGuess != 'PreDef': self.iGuess = iGuess
        lFit,ilist = [],[]
        for (igamma,ip),pdata in self.Rat_Stats['boot'].groupby(level=('current_operator','momentum')):
            this_key = (igamma,ip)
            if this_key in self.Rat_Fit_Stats.index:
                if WipeFit:
                    self.Rat_Fit_Stats.drop([(igamma,ip)],inplace=True)
                else:
                    continue
            if isinstance(self.fit_range,str):
                ifitmin,ifitmax = np.array(unxmlfitr(self.fit_range),dtype=np.float64)
            elif isinstance(self.fit_range,(tuple,list,np.ndarray)):
                ifitmin,ifitmax = np.array(self.fit_range,dtype=np.float64)
            else:
                print(self.fit_range)
                raise IOError('incorrect format for fitr#-# ' )
            # ifitmin,ifitmax = np.array(unxmlfitr(ifit),dtype=np.float64)

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
            fit_info['iGuess'] = self.iGuess
            fit_info['name'] = '_'.join(this_key)
            thisff = sff.SetOfFitFuns(data=pdata[this_key],name=fit_info['name'])
            thisff.ScanRange_fitr(self.fit_range,fit_info=fit_info,min_fit_len=min_fitr)
            lFit.append(thisff)
            ilist.append(this_key)

        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.Rat_Fit_Col_Names)
            if 'boot' in self.Rat_Fit_Stats.columns:
                this_df = pa.DataFrame()
                this_df.loc[:,'boot'] = pa.Series(lFit,index=indicies)
                self.Rat_Fit_Stats = self.Rat_Fit_Stats.append(this_df)
            else:
                self.Rat_Fit_Stats.loc[:,'boot'] = pa.Series(lFit,index=indicies)
            def DoFits(val):
                val.DoFits(show_timer=show_timer)
                return val
            self.Rat_Fit_Stats['boot'] = self.Rat_Fit_Stats['boot'].apply(DoFits)
            self.Write()


    def GetFuns(self):
        if self.Rat_Fit_Stats is not None and 'boot' in self.Rat_Fit_Stats.columns:
            for ival in self.Rat_Fit_Stats['boot'].values:
                ival.GetFuns()
        if hasattr(self.C2class,'GetFuns'):
            self.C2class.GetFuns()
        if hasattr(self.C2Set,'GetFuns'):
            self.C2Set.GetFuns()
        if not hasattr(self,'Fun'):
            self.Fun = ReadFuns(self.Fun_name)[0]

    def RemoveFuns(self):
        if self.Rat_Fit_Stats is not None and 'boot' in self.Rat_Fit_Stats.columns:
            for ival in self.Rat_Fit_Stats['boot'].values:
                ival.RemoveFuns()
        if hasattr(self.C2class,'GetFuns'):
            self.C2class.RemoveFuns()
        if hasattr(self.C2Set,'GetFuns'):
            self.C2Set.RemoveFuns()
        if hasattr(self,'Fun'):
            self.Fun_name = self.Fun.__name__
            WriteFuns(self.Fun)
            del self.Fun

    def Write(self):
        def FixDictArray(thislist,ikey):
            outDict = OrderedDict()
            for ic,ilist in enumerate(thislist):
                outDict[ikey+str(ic+1)] = ilist
            return outDict

        # C2name = self.C2class.name.replace(self.C2class.jsm,self.C3class.jsm)
        # C2LL = self.C2class.LegLab.replace(self.C2class.jsm,self.C3class.jsm)
        # self.C2class.SetCustomName(C2name,C2LL)
        if hasattr(self.C3class,'Write'):
            self.C3class.Write()

        if not hasattr(self,'show_cfgs'):
            self.show_cfgs = True

        outDict = ODNested()
        outDict['name'] = self.name
        if hasattr(self.C2class,'Write'):
            self.C2class.Write()
            outDict['C2file'] = self.C2class.HumanFile
        outDict['C3file'] = self.C3class.HumanFile
        outDict['kud'] = self.C3class.kud
        outDict['ks'] = self.C3class.ks
        outDict['nt'] = self.nt
        outDict['nxyz'] = self.nxyz

        for istream,incfg in zip(self.C3class.stream_list,self.C3class.ncfg_list):
            outDict['ncfg_'+istream] = incfg
        outDict['nxsrc_avg'] = self.C3class.xsrcAvg
        outDict['nMeas'] = self.C3class.nMeas
        outDict['nboot'] = self.nboot
        outDict['tsrc'] = self.C3class.tsrc
        outDict['tsink'] = self.C3class.tsink
        # outDict['xsrc'] = FixDictArray(self.xsrc,'xsrc')
        outDict['ism'] = self.C3class.ism
        outDict['sink_type'] = self.C3class.jsm
        outDict['Projector'] = self.C3class.Proj
        outDict['Doub_Sing'] = self.C3class.DS
        outDict['ppmom'] = self.C3class.ppmom

        for icg,igamma in enumerate(self.C3class.gammalist):
            outDict['gamma'+str(icg+1)] = igamma
        for icq,iqmom in enumerate(self.C3class.qmom):
            outDict['qmom'+str(icq+1)] = iqmom

        excel_params = pa.Series(deepcopy(outDict))
        outDict['gamma_list'] = FixDictArray(self.C3class.gammalist,'gamma')
        outDict['qmom_list'] = FixDictArray(self.C3class.qmom,'qmom')
        for icg,igamma in enumerate(self.C3class.gammalist):
            del outDict['gamma'+str(icg+1)]
        for icq,iqmom in enumerate(self.C3class.qmom):
            del outDict['qmom'+str(icq+1)]



        for (igamma,ip),fitdata in self.Rat_Fit_Stats['boot'].items():
            outDict['Fits'][igamma][ip] = fitdata.GetOutputDict()

        for col_key,col_data in self.Rat_Stats.items():
            if 'Avg' in col_key or 'Std' in col_key: continue
            for (igamma,ip,it),idata in col_data.items():
                outDict[col_key][igamma][ip][it] = AvgStdToFormat(idata.Avg,idata.Std)

        y = getattr(self, "SqrtFact", None)
        if y is not None and hasattr(self.SqrtFact,'iteritems'):
            for (ip,it),idata in self.SqrtFact.items():
                if np.isnan(idata.Avg): idata.Stats()
                outDict['SqrtFact'][ip][it] = AvgStdToFormat(idata.Avg,idata.Std)
            for (ip,it),idata in self.C2class.C2_Stats['boot'].items():
                if np.isnan(idata.Avg): idata.Stats()
                outDict['G2'][ip][it] = AvgStdToFormat(idata.Avg,idata.Std)


        if any(np.array(self.C3class.ncfg_list) > 0) and self.show_cfgs:
            for ((istream,iccfg),icfg),ixsrc_list in zip(iter(self.C3class.C3_cfgs['configs'].items()),self.C3class.C3_cfgs['xsrc_list'].values):
                outDict['cfglist'][istream][icfg] = 'n_xsrc_'+str(len(ixsrc_list))

        ## Write human readable file with data
        WriteXml(self.HumanFile,{'Results':outDict})

        if any(np.array(self.C3class.ncfg_list) > 0):
            for ((istream,iccfg),icfg),ixsrc_list in zip(iter(self.C3class.C3_cfgs['configs'].items()),self.C3class.C3_cfgs['xsrc_list'].values):
                outDict['cfglist'][istream][icfg] = 'n_xsrc_'+str(len(ixsrc_list))

        ## TODO, include fitting in the excell file
        WriteExcel(self.ExcelFile,{'Rat_results':deepcopy(self.Rat_Stats)},cfglist=self.C3class.C3_cfgs,params=excel_params)

        ## pickles rest of data for reading
        self.RemoveFuns()
        WritePickle(self.PickleFile,self.__dict__)
        self.GetFuns()

    ##################################


    def ReadAndWrite(self,WipeData=True,DefWipe=False,CheckCfgs=False):
        self.Read(WipeData=WipeData,DefWipe=DefWipe)
        self.CreateC2Sink()
        self.Bootstrap(WipeData=WipeData,DefWipe=DefWipe)
        self.CalcRatio()
        self.Fit()
        self.Write()

    ## add things you want to update to the depreciated pickled files here.
    def FixDict(self,thisdict):
        thisdict['LegLab'] = self.LegLab
        return thisdict

    def GetTflowList(self,load_dict):
        return load_dict['Rat_Stats'].index.levels[2]

    def FixRead(self,readdict):
        rewrite = False
        if not isinstance(readdict['C2class'],str):
            if not isinstance(readdict['C2Set'],str):
                readdict['C2Set'] = 'Complete'
                rewrite = True
        return readdict,rewrite

    def LoadPickle(self,DefWipe=False,WipeData = True,CheckCfgs=False,CheckMom=True):
        # print
        # print 'NJN'
        # print self.HumanFile
        # print self.PickleFile
        if os.path.isfile(self.PickleFile) and not DefWipe:
            print('Loading Pickle file ' , self.PickleFile)
            loadeddict = ReadPickleWrap(self.PickleFile)
            loadeddict,rewrite = self.FixRead(loadeddict)
            checklist = self.thischecklist
            if CheckCfgs: checklist = checklist + ['cfglist']
            if CheckMom: checklist = checklist + ['checkmomlist']
            if CheckClass(self.__dict__,loadeddict,checklist):
                loadeddict = self.FixDict(loadeddict)
                self.__dict__.update(loadeddict)
                if rewrite: self.Write()
                if Wipe_All_Fits:
                    self.Rat_Fit_Stats = pa.DataFrame()
                elif len(self.Rat_Fit_Stats.values)> 0 and not isinstance(self.Rat_Fit_Stats.values[0][0],sff.SetOfFitFuns):
                    print('Legacy file, fits need to be recomputed')
                    self.Rat_Fit_Stats = pa.DataFrame()
            else:
                if os.path.isfile(self.PickleFile+'.bak'):
                    self.PickleFile = self.PickleFile+'.bak'
                    self.HumanFile = self.HumanFile+'.bak'
                    self.ExcelFile = self.ExcelFile.replace('.xlsx','.bak.xlsx')
                    print('using backupfile for '+self.HumanFile)
                    self.LoadPickle(DefWipe=DefWipe,WipeData=WipeData,CheckCfgs=CheckCfgs,CheckMom=CheckMom)
                    return
                print('Warning, file ' , self.HumanFile , ' has different parameters than this instance, reading and writing over file now')
                print()
                self.ReadAndWrite(WipeData=WipeData,DefWipe=DefWipe)
            # self.Write()
        else:
            if not os.path.isfile(self.PickleFile):
                print()
                print('pickle file not found: \n', self.PickleFile)
            self.ReadAndWrite(WipeData=WipeData,DefWipe=DefWipe)
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
        if self.thisshiftscale == 'Not Set':
            xlen = np.abs(xlims[1]-xlims[0])
            self.thisshiftscale = self.thisshift*xlen
            return self.thisshiftscale
        else:
            return self.thisshiftscale




    ## make sure all xlims are the same, because this sets the window of xlim (needed for calculating shifts)
    ## gammamomlist is a dictionary where the keys are the gamma matricies, and the values are lists of momenta
    ## TODO: SHIFTS NEED TO BE FIXED!! LIKE NOW. A better way needs to be implemented.
    def Plot(self,plot_class,xlims='tsink',gammamomlist={},thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine'):
        if xlims == 'tsink': xlims = [0,int(self.C3class.tsink.replace('tsink',''))+1]
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        self.thisshiftscale = 'Not Set'
        thisshift = self.GetShift(xlims,thisshift)
        if 'boot' not in self.Rat_Stats.columns:
            raise EnvironmentError('no data to plot, maybe try to read in and bootstrap before plotting')
        if len(list(gammamomlist.keys())) == 0:
            gammamomlist = self.defgammamomlist


        # for (igamma,ip),Rat_data in self.Rat_Stats['boot'].groupby(level=('current_operator','momentum')):
            # if igamma not in gammamomlist.keys(): continue
            # if ip not in gammamomlist[igamma]: continue
        igamma = list(gammamomlist.keys())[0]
        ip = self.latparams.TOpform(list(gammamomlist.values())[0][0],pref='q')
        ip = ip.replace('qq','q')
        if (igamma,ip) not in self.Rat_Stats['boot'].index:
            print(', '.join([igamma,ip]) ,' not computed')
            return plot_class
        # tlist,databoot = zip(*list(self.Rat_Stats['boot'][igamma,ip].iteritems()))
        this_key = (igamma,ip,slice(None))
        databoot = self.Rat_Stats['boot']
        dataavg = databoot.apply(lambda x : x.Avg)
        dataerr = databoot.apply(lambda x : x.Std)
        tlist = list(databoot[this_key].keys())
        # dataavg,dataerr = [],[]
        # for idata in databoot[xlims[0]:xlims[1]]:
        #     dataavg.append(idata.Avg)
        #     dataerr.append(idata.Std)
        tlist = tlist[xlims[0]:xlims[1]]
        # print 'plotting ' , igamma,ip
        # pleg,gammaleg = '',''
        # if len(gammamomlist[igamma]) > 1:
        # pleg = ' $'+ip+'$'
        # if '_cmplx' in igamma:
        #     gammaleg = ' $i'+igamma.replace('_cmplx','')+'$'
        # else:
        #     gammaleg = ' $'+igamma+'$'
        # if len(gammamomlist.keys()) > 1:
        tlistplot,self.centert = CenterT(np.array(list(map(untstr,tlist)))-twopttshift+thisshift)
        hold_series = pa.Series()
        # hold_series['x_data'] = tlistplot
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = dataavg
        hold_series['yerr_data'] = dataerr
        hold_series['xerr_data'] = None
        hold_series['key_select'] = this_key
        hold_series['type'] = 'error_bar_vary'
        hold_series['fit_class'] = None
        hold_series['label'] = self.LegLab#+pleg+gammaleg
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        hold_series['xdatarange'] = None
        hold_series['Phys'] = None
        hold_series['ShowPar'] = None
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class

    ## state is 1 or 2 state fit
    ## fitr is fit range e.g. fitr5-10
    ## fit range of 'Data' is just the original fit range
    ## tflowlist is just for consistancy with other function in NJNQ, does nothing
    def GetBestFitr(self,igamma,ip):
        if 'boot' not in self.Rat_Fit_Stats.columns or (igamma,ip) not in self.Rat_Fit_Stats.index:
            errstr = ' '.join(igamma,ip)+ ' not foun in ratio fit list, cannot find best'
            raise EnvironmentError(errstr)
        chi_str,cut_fit_result = self.Rat_Fit_Stats.at[(igamma,ip),'boot'].GetBestChi()
        return r'$\chi^{2}_{pdf}={:.3f}$'.format(chi_str),cut_fit_result

    def Get_Extrapolation(self,keep_keys=slice(None),fmted=True):
        if 'boot' in self.Rat_Fit_Stats:
            return sff.PullSOFSeries_Wpar(self.Rat_Fit_Stats['boot'],keep_keys=keep_keys,fmted=fmted)
        else:
            return pa.Series()

    def FitPlot(self,plot_class,fitr,xlims = 'Data',gammamomlist={},
                thiscol='PreDefine',thisshift='PreDefine',tflowlist=[]):
        self.CheckCol(thiscol)
        thisshift = self.GetShift(xlims,thisshift)
        if self.centert == 'Not Set':
            print('print centert not set, defaulting to 0')
            self.centert = 0.9

        if 'boot' not in self.Rat_Fit_Stats.columns:
            raise EnvironmentError('No fitting has been done, maybe try to Fit before plotting fit')
        if len(list(gammamomlist.keys())) == 0:
            gammamomlist = self.defgammamomlist
        # for igamma,momlist in gammamomlist.iteritems():
        #     for ip in momlist:
        igamma,momlist = next(iter(gammamomlist.items()))
        ip = self.latparams.TOpform(momlist[0],pref='q')
        ip = ip.replace('qq','q')
        if fitr == 'Best':
            chi_legend,fitr = self.GetBestFitr(igamma,ip)
        if 'boot' not in self.Rat_Fit_Stats.columns:
            print('no fits done, attempting fit now')
            self.Fit(fit_range=fitr)
        fit_data = sff.PullSOFSeries(self.Rat_Fit_Stats['boot'])
        this_fun = fit_data.index[0][2]
        this_key = (igamma,ip,this_fun,fitr)
        if this_key not in fit_data:
            print(this_key , 'not present in Rat Fit Stats, fitting now')
            self.Fit(fit_range=fitr)
            fit_data = sff.PullSOFSeries(self.Rat_Fit_Stats['boot'])
            this_fun = fit_data.index[0][2]
            this_key = (igamma,ip,this_fun,fitr)

        hold_series = null_series
        hold_series['key_select'] = this_key
        hold_series['type'] = 'fit_vary'
        hold_series['fit_class'] = fit_data
        # hold_series['label'] = fit_data[this_key].name +' ' + chi_legend
        hold_series['label'] = 'ACTUAL_VALUE_fit'
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        hold_series['xdatarange'] = xlims
        hold_series['ShowPar'] = r'Par1'
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class


    ## ##################### ###############################################

    ## Comparisons

    def __str__(self):
        return '_'.join([self.SubDir,self.name])


    def __iter__(self):
        return self

    def iteritems(self):
        return iter(self.Rat_Stats['boot'].items())

    def items(self):
        return list(self.Rat_Stats['boot'].items())

    def itemsAvgStd(self):
        outlist = []
        for ivals in self.Rat_Stats.itertuples():
            outlist.append((ivals[0],ivals[2],ivals[3]))
        return outlist

    def values(self):
        return self.Rat_Stats['boot'].values
    def itervalues(self):
        return iter(self.Rat_Stats['boot'].values)

    def valuesAvgStd(self):
        outlist = []
        for ivals in self.Rat_Stats.itertuples():
            outlist.append((ivals[2],ivals[3]))
        return outlist

    def keys(self):
        return self.Rat_Stats.index


    def __setitem__(self, key, value ):
        self.Rat_Stats['boot'][key] = value

    def __getitem__(self, key):
        return self.Rat_Stats['boot'][key]

    def __reversed__(self):
        self.Rat_Stats = self.Rat_Stats.iiloc[::-1]
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


## BLOCK FOR OVERLOADING ARITHMETIC OPERATIONS ####
    def ReduceIndicies(self):
        indexl = []
        for (igamma,ip,it) in list(self.Rat_Stats['boot'].keys()):
            first_gamma = MultiSplit(igamma,op_str_list+['FO','FOFull'])[0]
            this_gamma = first_gamma.replace('FO','')
            this_gamma = this_gamma.replace('FOFull','')
            this_gamma = this_gamma.replace('__','_')
            if this_gamma[-1] == '_':
                this_gamma = this_gamma[:-1]
            this_p = MultiSplit(ip,op_str_list+['FO','FOFull','_'])[0]
            indexl.append((this_gamma,this_p,it))
        indicies = pa.MultiIndex.from_tuples(indexl,names=self.Rat_Col_Names)
        this_series = pa.Series(self.Rat_Stats['boot'].values,index=indicies)
        self.Rat_Stats = pa.DataFrame(columns=['boot'],index=indicies)
        self.Rat_Stats.loc[:,'boot'] = this_series
        self.Stats()


    def sin(self):
        result = RatioCorr(thisnboot=self.nboot, cfglist=self.C3class.C3_cfgs[['configs','xsrc_list']],
                         Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift,
                         name = 'sin('+self.name+')',filename='sin_'+self.filename)
        result.Rat_Stats.loc[:,'boot'] = np.sin(self.Rat_Stats.loc[:,'boot'])
        return result

    def cos(self):
        result = RatioCorr(thisnboot=self.nboot, cfglist=self.C3class.C3_cfgs[['configs','xsrc_list']],
                         Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift,
                         name = 'cos('+self.name+')',filename='cos_'+self.filename)
        result.Rat_Stats = pa.DataFrame(columns = ['boot'],index=self.Rat_Stats.index)
        result.Rat_Stats.loc[:,'boot'] = np.cos(self.Rat_Stats.loc[:,'boot'])
        return result


    def Overload_wrap(self,Rat2,this_fun):
        ## TODO, include some name and filename formatting
        thisname,thisfilename = self.name,self.filename
        if hasattr(Rat2,'name'):
            thisname += '_'+this_fun.__name__+'_'+Rat2.name
        if hasattr(Rat2,'filename'):
            thisfilename += '_'+this_fun.__name__+'_'+Rat2.filename
        result = RatioCorr(thisnboot=self.nboot, cfglist=self.C3class.C3_cfgs[['configs','xsrc_list']],
                         Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.Rat_Stats = pa.DataFrame(columns = ['boot'])
        # result.UpdateName('+',self,Rat2)
        if 'boot' in self.Rat_Stats:
            if isinstance(Rat2,RatioCorr):
                # result.Rat_Stats.loc[:,'boot'] = this_fun(self.Rat_Stats.loc[:,'boot'],Rat2.Rat_Stats.loc[:,'boot'])
                if 'boot' in Rat2.Rat_Stats:
                    vall,ilist = [],[]
                    for (igamma,ip,it),irat in self.Rat_Stats['boot'].items():
                        if it not in Rat2.Rat_Stats['boot'].index.swaplevel(0,2):
                            out_str = 'missmatch in euclidian times when combining ratio functions \n'
                            out_str += str(self) + '\n'
                            out_str += this_fun.__name__+'\n'
                            out_str += str(Rat2)
                            raise EnvironmentError(out_str)
                        for (igamma2,ip2,itdump),irat2 in Rat2.Rat_Stats.loc[(slice(None),slice(None),it),'boot'].items():
                            igammacomb,ipcomb = igamma,ip
                            if igamma2 not in igamma:
                                igammacomb += '_'+this_fun.__name__+'_'+igamma2
                            if ip2 not in ip2:
                                ipcomb += '_'+this_fun.__name__+'_'+ip2
                            ilist.append((igammacomb,ipcomb,it))
                            vall.append(this_fun(irat,irat2))
                    if len(ilist) > 0:
                        indicies = pa.MultiIndex.from_tuples(ilist,names=self.Rat_Col_Names)
                        result.Rat_Stats.loc[:,'boot'] = pa.Series(vall,index=indicies)
                else:
                    raise EnvironmentError('Run Bootstrap on second class before combining ratio functions')
            else:
                try:
                    result.Rat_Stats.loc[:,'boot'] = this_fun(self.Rat_Stats.loc[:,'boot'], Rat2)
                except Exception as err:
                    print(type(Rat2), Rat2)
                    print(type(self))
                    print(self)
                    print(self.name)
                    print(type(self.Rat_Stats['boot'].iloc[0]))
                    print(self.Rat_Stats['boot'].iloc[0]/Rat2)
                    print(this_fun(self.Rat_Stats['boot'].iloc[0],Rat2))
                    raise EnvironmentError(str(err)+'\nInvalid value to combine with RatioCorr class')
            return result
        else:
            raise EnvironmentError('Run Bootstrap on first class before combining ratio functions')

    def __add__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioCorr.parentlist]):
            return Rat2.__radd__(self)
        else:
            return self.Overload_wrap(Rat2,op.add)

    def __sub__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioCorr.parentlist]):
            return Rat2.__rsub__(self)
        else:
            return self.Overload_wrap(Rat2,op.sub)

    def __mul__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioCorr.parentlist]):
            return Rat2.__rmul__(self)
        else:
            return self.Overload_wrap(Rat2,op.mul)

    def __div__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioCorr.parentlist]):
            return Rat2.__rdiv__(self)
        else:
            return self.Overload_wrap(Rat2,op.truediv)

    def __pow__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioCorr.parentlist]):
            return Rat2.__rpow__(self)
        else:
            return self.Overload_wrap(Rat2,op.pow)


    ## Right multiplication functions


    def __radd__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioCorr.parentlist]):
            return Rat2.__add__(self)
        else:
            return self.Overload_wrap(Rat2,flip_fun_2arg(op.add))

    def __rsub__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioCorr.parentlist]):
            return Rat2.__sub__(self)
        else:
            return self.Overload_wrap(Rat2,flip_fun_2arg(op.sub))

    def __rmul__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioCorr.parentlist]):
            return Rat2.__mul__(self)
        else:
            return self.Overload_wrap(Rat2,flip_fun_2arg(op.mul))

    def __rdiv__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioCorr.parentlist]):
            return Rat2.__div__(self)
        else:
            return self.Overload_wrap(Rat2,flip_fun_2arg(op.truediv))


    def __rpow__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioCorr.parentlist]):
            return Rat2.__pow__(self)
        else:
            return self.Overload_wrap(Rat2,flip_fun_2arg(op.pow))




####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################

##UP TO HERE TODO:

class RatioFOCorr(object):
    """

    Rat(q,tau,pp,t,tflow) uses bootstrap class
    pp  = sink momentum
    t   = source sink separation time
    q   = current insersion momentum
    tau = current insersion time
    tflow = flow time

    Info is a dictionary containing information for the correlator:
    see comments below to see what is required.
    missing elements just sets field to default


    """

    def Construct_RatFO_File(this_file):
        return Construct_File_Object(this_file,RatioFOCorr)

    ## TODO: create variational flag and interface with this function.

    ## array of all classes that inherit this class. Append this when adding new classes please!
    parentlist = ['SetOfCorrs.SetsOfRatios','RatioCorrelators.RatioFOFullCorr']

    def __init__(self, thisnboot=nboot, cfglist={}, Info={},name='',
                 filename='',thissym='Not Set',thiscol='Not Set',thisshift=0.0
                 ,man_load_cfgs=False):


        self.current = 0
        ## Initialising plotting parameters
        self.thisshift = thisshift
        self.thiscol = thiscol
        self.thissym = thissym
        ## the actual shift ammount needs to be scaled by the x axis size
        ## see GetShift
        self.thisshiftscale = 'Not Set'

        self.C3FOclass = NJNQCorr(thisnboot=thisnboot, cfglist=cfglist, Info=Info,name=name,filename=filename)
        self.C3file = self.C3FOclass.PickleFile
        if 'fo_for_cfgs' in list(Info.keys()): self.fo_for_cfgs = Info['fo_for_cfgs']
        else: self.fo_for_cfgs = False
        self.tflow_list = self.C3FOclass.FO.tflowlist
        self.thischecklist = ['C3file','Sym','fo_for_cfgs','tflow_list']
        self.checkmomlist = self.C3FOclass.NJN.checkmomlist

        self.latparams = self.C3FOclass.NJN.latparams
        Info['pmom'] = [self.latparams.GetAvgmomstr(ip) for ip in self.C3FOclass.NJN.qmom]
        if '0 0 0' not in Info['pmom']:
            Info['pmom'] = ['0 0 0'] + Info['pmom']

        self.nt = self.latparams.nt
        self.nxyz = self.latparams.nxyz
        self.dim_label = 'RC'+str(self.nxyz)+'x'+str(self.nt)
        self.dim_ll = str(self.nxyz)+r'^{3}\times ' + str(self.nt)

        ## both current and sink momentum are requred for the two point correlator construction when constructing the ratio function
        ## but the two point correlator has been averaged over, so we need to get the averaged momentum anaolgy
        ## bah, this was too hard to keep unique files for the two-point correlator.
        ## easier just to read and write the whole 2 point mom file.
        Info['pmom'] = self.latparams.GetMomList(list_type='str_Avg',list_order='squared_rev')
        ## e.g. q = [-1,1,-1] corresponds to p = [1 1 1] in the two-point correlator
        # ppforp = self.latparams.GetAvgmomstr(Info['ppmom'])
        # if ppforp not in Info['pmom']: Info['pmom'].append(ppforp)


        if 'show_cfgs' in list(Info.keys()): self.show_cfgs = Info['show_cfgs']
        else: self.show_cfgs = True

        if 'Rat_comb_alpha_fun' in list(Info.keys()):
            if Info['Rat_comb_alpha_fun'] == True:
                print('Rat_comb_alpha_fun has been set to false, turn back on if you want to use it again!.')
            self.alpha_rat = False
        else: self.alpha_rat = False



        ## Does symmetrising over smearings
        if 'Symmetric' in list(Info.keys()): self.Sym = Info['Symmetric']
        else: self.Sym = False ## default to no symmetrizing, (usually full matrix isn't computed nowerdays.

        # self.C2filelist = [ikey.PickleFile for ikey in self.C2classlist.itervalues()]

        if 'Rat_alpha' in list(Info.keys()):
            if  Info['Rat_alpha'] is True:
                print('WARNING, Rat_alpha hard coded to be off, turn this on if you decide you want to use it again!')
            # self.Rat_alpha = Info['Rat_alpha']
            self.Rat_alpha = False
        else: self.Rat_alpha = False
        if self.Rat_alpha:
            Info2pt = deepcopy(Info)
            Info2pt['Interp'] = 'CPOdd'
            self.sinklab,self.jsmcoeffs,self.sm2ptList,InfoList = GetInfoList(Info2pt)
            self.C2Set = SetOfNNQ(cfglist=cfglist,InfoDict = InfoList)
        else:
            Info2pt = deepcopy(Info)
            Info2pt['Interp'] = 'CPEven'
            self.sinklab,self.jsmcoeffs,self.sm2ptList,InfoList = GetInfoList(Info2pt)
            self.C2Set = SetOfTwoPt(cfglist=cfglist,InfoDict = InfoList)
        self.C2class = 'Not Set'

        if Debug_Mode:
            print('testing getinfolist')
            print(self.sinklab)
            print(self.jsmcoeffs)
            print(self.sm2ptList)
            print()

        if self.C3FOclass.FO.Observ == '':
            for iinfo in InfoList: iinfo['fo_for_cfgs'] = False
        else:
            for iinfo in InfoList: iinfo['fo_for_cfgs'] = self.C3FOclass.FO.Observ

        self.defgammamomlist = OrderedDict()
        for igamma in self.C3FOclass.NJN.gammalist:
            self.defgammamomlist[igamma] = self.C3FOclass.NJN.qform

        self.Rat_Col_Names = ['current_operator','momentum','flow_times','source_sink_separation']
        self.Rat_Stats = pa.DataFrame()

        '''
        self.Rat_Stats DataFrame:

        columns =   boot,   Avg,        Std

        rows =      current_operator,   momenta,    flow_times,    time separation
        '''

        self.Rat_Fit_Col_Names = ['current_operator','momentum','flow_times']
        self.Rat_Fit_Stats = pa.DataFrame()

        '''
        self.Rat_Fit_Stats DataFrame:

        columns =   boot,   Avg,        Std

        rows =      current_operator, momenta,  flow_times, time fit range

        Avg/Std elements are pa.Series instances with:
            index = fit parameters
            values = Avg/Std values
        '''



        self.nboot = thisnboot
        self.centert = 'Not Set'


        ## Info can contain fit ranges to be calculated
        ## must bein form Info['FitsRF'] = [ 'fit#-#' , ...]
        if 'FitsRF' in list(Info.keys()): self.PredefFits = Info['FitsRF']
        else: self.PredefFits = []


        if 'RFFitFun' in list(Info.keys()): self.Fun,self.npar = Info['RFFitFun']
        else: self.Fun,self.npar = (ConstantFitFun,1)

        ## Info can contain initial guess for fits, make sure it array of correct length
        if 'iGuessRF' in list(Info.keys()): self.iGuess = np.array(Info['iGuessRF'])
        else: self.iGuess = 'None' ## Defaults to value set for function in FitFunctions.py

        ## is needed?
        self.Info = Info

        self.SetCustomName(string=name,stringfn=filename)
        if not man_load_cfgs:
            self.LoadCfgs(cfglist,name,filename)

    def LoadCfgs(self,cfglist,name='',filename=''):
        if len(cfglist) == 0:
            self.CombineCfgLists()
        else:
            self.ImportCfgList(cfglist)
        self.SetCustomName(string=name,stringfn=filename)

    ## just for consistancy with other classes
    def GetCfgList(self):
        self.CombineCfgLists()

    def CombineCfgLists(self):
        if isinstance(self.C2class,str):
            thiscfglist = CombineCfgs(self.C3FOclass.NJNQ_cfgs,self.C2Set.cfglist)
            # print
            # print len(thiscfglist.keys()), len(self.C3FOclass.cfglist.keys()), len(self.C2Set.cfglist.keys())
        else:
            thiscfglist = CombineCfgs(self.C3FOclass.NJNQ_cfgs,self.C2class.C2_cfgs)
            # print
            # print len(thiscfglist.keys()), len(self.C3FOclass.cfglist.keys()), len(self.C2class.cfglist.keys())
        self.ImportCfgList(thiscfglist)




    def SetCustomName(self,string='',stringfn='',stringLL=''):
        if string == '':
            self.name = 'Rat_'+self.C3FOclass.name
        else:
            self.name = string
        if stringfn == '':
            self.filename = 'Rat_'+self.C3FOclass.filename
        else:
            self.filename = stringfn
        if stringLL == '':
            thisops = self.C3FOclass.FO.Observ
            self.LegLab = '$'+'\ '.join([   self.dim_ll,self.latparams.GetPionMassLab(),self.C3FOclass.NJN.stream_name,thisops,
                                            self.C3FOclass.DS,self.C3FOclass.NJN.Proj])+'$' ## customise this how you like
        else:
            self.LegLab = stringLL
        if self.Rat_alpha:
            self.filename += '_Rata'
            self.name += '_Rata'
            self.LegLab += r'\ $G^{Q}_{2}$'
        self.SubDir = '_'.join([self.C3FOclass.NJN.Proj,self.C3FOclass.DS,self.C3FOclass.NJN.ppform])
        self.Ratdir = outputdir + '/'+self.dim_label+self.C3FOclass.NJN.kappafolder+'/Rat/'+self.C3FOclass.FO.Observ+'/'+self.SubDir+'/'
        mkdir_p(self.Ratdir+'/Pickle/')
        mkdir_p(self.Ratdir+'/Excel/')
        self.HumanFile = self.Ratdir+self.filename+'.xml'
        self.ExcelFile = self.Ratdir+'/Excel/'+self.filename+'.xlsx'
        self.PickleFile = self.Ratdir+'/Pickle/'+self.filename+'.py3p'



    def ImportCfgList(self,cfglist):
        self.cfglist = cfglist
        self.C2Set.ImportCfgList(cfglist)
        if self.C2class != 'Not Set': self.C2class.ImportCfgList(cfglist)
        self.C3FOclass.ImportCfgList(cfglist)

    # def ImportCfgList(self,cfglist):
    #     for iC2 in self.C2classlist.itervalues():
    #         iC2.ImportCfgList(cfglist)
    #     if self.C2class != 'Not Set': self.C2class.ImportCfgList(cfglist)
    #     self.C3FOclass.ImportCfgList(cfglist)


    def Stats(self):
        if hasattr(self.C2class,'Stats'):
            self.C2class.Stats()
        if hasattr(self.C3FOclass,'Stats'):
            self.C3FOclass.Stats()
        if 'boot' not in self.Rat_Stats: return
        lAvg,lStd = [],[]
        for (igamma,ip,it,itflow),iRat in self.items():
            iRat.Stats()
            lAvg.append(iRat.Avg)
            lStd.append(iRat.Std)
        indicies = pa.MultiIndex.from_tuples(self.Rat_Stats.index,names=self.Rat_Col_Names)
        self.Rat_Stats.loc[:,'Avg'] = pa.Series(lAvg,index=indicies)
        self.Rat_Stats.loc[:,'Std'] = pa.Series(lStd,index=indicies)


    def RemoveVals(self):
        # for iC2 in self.C2classlist.itervalues():
        #     iC2.RemoveVals()
        self.C2class.RemoveVals()
        self.C3FOclass.RemoveVals()

    ## self.Rat = [ igamma , ip, it , iconf ]
    ## configs must be the same, so we check for them
    def Read(self,WipeData=True,DefWipe=False):
        if isinstance(self.C2class,str): self.C2Set.LoadPickle(DefWipe=DefWipe,CheckCfgs=True,NoFit=True)
        self.C3FOclass.LoadPickle(WipeData=WipeData,DefWipe=DefWipe,CheckCfgs=True,CheckMom=True)

    def CreateC2Sink(self):
        if isinstance(self.C2class,str):
            thisname,thisleglab,jsmlab = CreateNameAndLL(self.C2Set.SetC2,self.sinklab)
            thisfun = GenSinkFun(self.jsmcoeffs)
            if Debug_Mode:
                print('check file')
                print(thisname)
                print(thisleglab)
                print(jsmlab)
                print(thisfun.__name__)
            self.C2class = self.C2Set.CombineCorrs(thisfun,filename=thisname,LegLab=thisleglab,jsmlab=jsmlab)
            self.C2Set = 'Complete'

    def Bootstrap(self,WipeData=True,DefWipe=False,CheckCfgs=False):
        if isinstance(self.C2class,str) :
            self.Read(WipeData=WipeData,DefWipe=DefWipe)
            self.CreateC2Sink()
        # self.C2class.Bootstrap()
        self.C3FOclass.Bootstrap(WipeData=WipeData,DefWipe=DefWipe)


    def CalcRatio(self):
        self.CalcSqrtFac()
        ## if DS has a divide in it, we ignore the sqrt factor as they cancel
        div_bool = 'divide' in self.C3FOclass.DS
        ## double the square root factor if we multiply two ds types
        mult_bool = 'multiply' in self.C3FOclass.DS
        lRat,lRatAvg,lRatStd = [],[],[]
        ilist = []
        tsinkint = int(self.C3FOclass.tsink.replace('tsink',''))
        for (igamma,ip,it,itflow),tC3FO in self.C3FOclass.NJNQ_Stats['boot'].items():
            ip2pt = self.latparams.GetAvgmomform(ip,pref='p')
            ict = int(untstr(it))
            if ict > tsinkint: continue
            if self.Rat_alpha:
                if (ip2pt,it,itflow) not in self.SqrtFact.index: continue
                if div_bool:
                    this_boot = tC3FO
                elif mult_bool:
                    this_boot = tC3FO*(self.SqrtFact[ip2pt,it,itflow]**2)
                else:
                    this_boot = tC3FO*self.SqrtFact[ip2pt,it,itflow]
            else:
                if (ip2pt,it) not in self.SqrtFact.index: continue
                if div_bool:
                    this_boot = tC3FO
                elif mult_bool:
                    this_boot = tC3FO*(self.SqrtFact[ip2pt,it]**2)
                else:
                    this_boot = tC3FO*self.SqrtFact[ip2pt,it]
            # warnings.warn('DEBUG, please check if tsink for 2-pt corr is correct via g4.')
            # print igamma, ip, it, self.C3class.tsink, tC3.Avg,self.SqrtFact[ip2pt,it].Avg

            this_boot.Stats()
            ilist.append((igamma,ip,itflow,it))
            lRat.append(this_boot)
            lRatAvg.append(this_boot.Avg)
            lRatStd.append(this_boot.Std)
            # self.Ratboot[igamma][ip][it] = tC3/self.C2class[ip.replace('q','p'),'t'+str(int(self.C3class.tsink.replace('tsink',''))+1)]
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.Rat_Col_Names)
            self.Rat_Stats.loc[:,'boot'] = pa.Series(lRat,index=indicies)
            self.Rat_Stats.loc[:,'Avg'] = pa.Series(lRatAvg,index=indicies)
            self.Rat_Stats.loc[:,'Std'] = pa.Series(lRatStd,index=indicies)



    def CalcSqrtFac(self):
        SqrtFact_hold = []
        ilist = []
        ipp = self.latparams.GetAvgmomform(self.C3FOclass.NJN.ppform,pref='p')
        if Debug_Mode:
            print()
            print('Outputting sqrt factor debug info:')
            print()
        for iq in self.C3FOclass.NJN.qform:
            ip = self.latparams.GetAvgmomform(iq,pref='p')
            if self.Rat_alpha:
                for itflow in self.tflow_list:
                    sqrt_fac = self.CalcSqrtFacMom(self.C2class.NNQ_Stats['boot'][ipp,:,itflow],
                                                            self.C2class.NNQ_Stats['boot'][ip,:,itflow],
                                                            self.C3FOclass.NJN.tsink)
                    for it,tdata in enumerate(sqrt_fac):
                        if (ip,tstr(it+1),itflow) not in ilist:
                            if Debug_Mode:
                                for ierr in tdata.GetNaNIndex():
                                    print(iq, itflow, tstr(it+1), ' iboot'+str(ierr).zfill(4))
                            tdata.Stats()
                            SqrtFact_hold.append(tdata)
                            ilist.append((ip,tstr(it+1),itflow))
            else:
                sqrt_fac = self.CalcSqrtFacMom(self.C2class.C2_Stats['boot'][ipp],
                                                        self.C2class.C2_Stats['boot'][ip],
                                                        self.C3FOclass.NJN.tsink)
                for it,tdata in enumerate(sqrt_fac):
                    if (ip,tstr(it)) not in ilist:
                        if Debug_Mode:
                            for ierr in tdata.GetNaNIndex():
                                print(iq, tstr(it), ' iboot'+str(ierr).zfill(4))
                        tdata.Stats()
                        SqrtFact_hold.append(tdata)
                        ilist.append((ip,tstr(it)))
            # if neg_sqrt: warnings.warn('negative square root in square root factor for '+ipp+iq)
        if self.Rat_alpha:
            indicies = pa.MultiIndex.from_tuples(ilist,names=['momentum','source_sink_separation','flow_time'])
        else:
            indicies = pa.MultiIndex.from_tuples(ilist,names=['momentum','source_sink_separation'])
        self.SqrtFact = pa.Series(SqrtFact_hold,index=indicies)

    def CheckNegForSqrt(self,thisboot,negsqrt):
        for ict in range(len(thisboot.bootvals)):
            if thisboot.bootvals[ict] < 0.0:
                negsqrt = True
                # if G2_abs:
                #     thisboot.bootvals[ict] = abs(thisboot.bootvals[ict])
                # else:
                thisboot.bootvals[ict] = float('NaN')
                # thisboot.bootvals[ict] = 0.0
        return thisboot,negsqrt


    ## taking squareroot of each element before combinding makes the numbers not too small
    def CalcSqrtFacMom(self,data2ptpp,data2ptp,thistsink):
        dictout = []
        thistsink = thistsink.replace('tsink','t')
        thistsinkp1 = 't'+str(int(thistsink.replace('t',''))+twopttshift)
        inttsink = int(thistsink.replace('t',''))
        # data2ptpp[thistsinkp1],NegSQRTRatFac = self.CheckNegForSqrt(data2ptpp[thistsinkp1],NegSQRTRatFac)
        # data2ptp[thistsinkp1],NegSQRTRatFac = self.CheckNegForSqrt(data2ptp[thistsinkp1],NegSQRTRatFac)
        data2ptpp = data2ptpp.apply(lambda x : x.Sqrt())
        data2ptp = data2ptp.apply(lambda x : x.Sqrt())
        # NegSQRTRatFac = any(data2ptpp.isnull().values) or any(data2ptp.isnull().values)
        # if Debug_Mode:
        #     print
        #     print 'p'
        #     for ip in data2ptp:
        #         print '     ',ip.Avg,ip.Std
        #     print 'pp'
        #     for ipp in data2ptpp:
        #         print '     ',ip.Avg,ip.Std
        #     print

        two = data2ptpp[thistsinkp1]*data2ptp[thistsinkp1]
        # two = two.Sqrt()
        # two = (data2ptpp[thistsinkp1]**0.5)*(data2ptp[thistsinkp1]**0.5)
        # two = two**(0.5)
        thisshift = twopttshift
        for it in range(thisshift,inttsink+thisshift):
            thiststr = tstr(it)
            tflip = tstr(inttsink-it+twopttshift+thisshift)
            # data2ptpp[thiststr],NegSQRTRatFac = self.CheckNegForSqrt(data2ptpp[thiststr],NegSQRTRatFac)
            # data2ptp[thiststr],NegSQRTRatFac = self.CheckNegForSqrt(data2ptp[thiststr],NegSQRTRatFac)
            # data2ptpp[tflip],NegSQRTRatFac = self.CheckNegForSqrt(data2ptpp[tflip],NegSQRTRatFac)
            # data2ptp[tflip],NegSQRTRatFac = self.CheckNegForSqrt(data2ptp[tflip],NegSQRTRatFac)
            # one = (data2ptpp[thiststr]**0.5)/(data2ptp[thiststr]**0.5)
            # three = (data2ptp[tflip]**0.5)/(data2ptpp[tflip]**0.5)
            one = data2ptpp[thiststr]/data2ptp[thiststr]
            three = data2ptp[tflip]/data2ptpp[tflip]
            inter = one*three
            ott = inter/two
            # ott.Sqrt()
            ott.Stats()
            # ott,NegSQRTRatFac = self.CheckNegForSqrt(ott,NegSQRTRatFac)
            # ott.Stats()
            dictout.append(ott)
        return dictout

    # def DoAlphaCombine(self,this_fun=op.truediv):
    #     if not hasattr(self,'alpha_class'):
    #         from SetsOfTwoPtCorrs import SetOfNNQ
    #         Info2pt = deepcopy(self.Info)
    #         Info2pt['pmom'] = ['0 0 0']
    #         Info2pt['Interp'] = 'CPOdd'
    #         sinklab,jsmcoeffs,sm2ptList,InfoList = GetInfoList(Info2pt)
    #         FOSet = SetOfNNQ(cfglist=self.cfglist,InfoDict=InfoList)
    #
    #         FOSet.LoadPickle(CheckCfgs=True)
    #
    #
    #         thisname,thisleglab,jsmlab = CreateNameAndLL(FOSet.SetC2,sinklab)
    #         thisfun = GenSinkFun(jsmcoeffs)
    #         # print thisname,thisleglab
    #         # print jsmcoeffs
    #         # print sinklab
    #         # print thisfun.__name__
    #         # print FOSet.SetC2.keys()
    #         # print FOSet.SetC2[FOSet.SetC2.keys()[0]].NNQ_Stats['boot']
    #         # for ikey,ival in FOSet.SetC2.iteritems():
    #         #     print ikey
    #         #     print ival.NNQ_Stats
    #         #     print
    #         self.alpha_class = FOSet.CombineCorrs(thisfun,filename=thisname,
    #                                         LegLab=thisleglab,jsmlab =jsmlab)
    #         del FOSet
    #     self.CombineWithAlpha(self.alpha_class,this_fun=this_fun)
    #
    # def CombineWithAlpha(self,alpha_in,this_fun):
    #     tsink_str = self.C3FOclass.NJN.tsink.replace('tsink','t')
    #     if isinstance(this_fun,basestring):
    #         this_fun = op_dict[this_fun.replace('alpha_rat_','').replace('_corr','')]
    #     if not isinstance(alpha_in,NNQCorr):
    #         print type(alpha_in)
    #         print alpha_in
    #         raise AttributeError('Please pass in NNQCorr class instance')
    #     if 'boot' not in self.Rat_Stats:
    #         raise AttributeError('please compute Rat_Stats before combining with alpha')
    #     if 'Alphaboot' not in alpha_in.NNQ_Stats :
    #         raise AttributeError('please compute NNQ_Stats in alpha before combining with alpha')
    #     if ('p000',tsink_str) not in alpha_in.NNQ_Stats['Alphaboot']:
    #         print alpha_in.NNQ_Stats['Alphaboot']
    #         raise AttributeError('alpha ratio needs p000,'+tsink_str+' result')
    #     alpha_p0 = alpha_in.NNQ_Stats.loc[('p000',tsink_str,slice(None)),'Alphaboot']
    #     vall,ilist,corrl = [],[],[]
    #     for (igamma,ip,itflow,it),irat in self.Rat_Stats['boot'].iteritems():
    #         if itflow not in alpha_p0.index.swaplevel(0,2):
    #             out_str = 'missmatch in flow times when CombineWithAlpha called \n'
    #             out_str += itflow
    #             raise EnvironmentError(out_str)
    #         ilist.append((igamma,ip,itflow,it))
    #         output = this_fun(irat,alpha_p0[('p000',tsink_str,itflow)])
    #         output.Stats()
    #         vall.append(output)
    #         corrl.append(irat.corr(alpha_p0[('p000',tsink_str,itflow)]))
    #     if len(ilist) > 0:
    #         indicies = pa.MultiIndex.from_tuples(ilist,names=self.Rat_Col_Names)
    #         self.Rat_Stats.loc[:,'alpha_rat_'+this_fun.__name__] = pa.Series(vall,index=indicies)
    #         self.Rat_Stats.loc[:,'alpha_rat_'+this_fun.__name__+'_corr'] = pa.Series(corrl,index=indicies)
    #         self.Write()

    def ImportFitRanges(self,fit_range):
        if 'PreDef' == fit_range:
            self.fit_range = self.PredefFits
        else:
            self.fit_range = fit_range
        if not isinstance(self.fit_range,str):
            # print 'fit_range in setoffits implementation is a single value varied inwards, selecting first'
            self.fit_range = self.fit_range[0]
        self.fit_range = rm_tsumfitr(self.fit_range)

    def DoAlphaCombine(self,this_fun=op.truediv):
        if not hasattr(self,'alpha_class'):
            from SetsOfTwoPtCorrs import SetOfNNQ
            Info2pt = deepcopy(self.Info)
            Info2pt['pmom'] = ['0 0 0']
            Info2pt['Interp'] = 'CPOdd'
            sinklab,jsmcoeffs,sm2ptList,InfoList = GetInfoList(Info2pt)
            FOSet = SetOfNNQ(cfglist=self.cfglist,InfoDict=InfoList)

            FOSet.LoadPickle(CheckCfgs=True)


            thisname,thisleglab,jsmlab = CreateNameAndLL(FOSet.SetC2,sinklab)
            thisfun = GenSinkFun(jsmcoeffs)
            # print thisname,thisleglab
            # print jsmcoeffs
            # print sinklab
            # print thisfun.__name__
            # print FOSet.SetC2.keys()
            # print FOSet.SetC2[FOSet.SetC2.keys()[0]].NNQ_Stats['boot']
            # for ikey,ival in FOSet.SetC2.iteritems():
            #     print ikey
            #     print ival.NNQ_Stats
            #     print
            self.alpha_class = FOSet.CombineCorrs(thisfun,filename=thisname,
                                            LegLab=thisleglab,jsmlab =jsmlab)
            del FOSet
        self.CombineWithAlpha(self.alpha_class,this_fun=this_fun)

    def CombineWithAlpha(self,alpha_in,this_fun=op.truediv):
        tsink_str = self.C3FOclass.NJN.tsink.replace('tsink','t')
        if isinstance(this_fun,str):
            if this_fun == 'None': return
            this_fun = op_dict[this_fun.replace('alpha_rat_','').replace('_corr','')]
        if not isinstance(alpha_in,NNQCorr):
            print(type(alpha_in))
            print(alpha_in)
            raise AttributeError('Please pass in NNQCorr class instance')
        if 'boot' not in self.Rat_Stats:
            raise AttributeError('please compute Rat_Stats before combining with alpha')
        if 'Alphaboot' not in alpha_in.NNQ_Stats :
            raise AttributeError('please compute NNQ_Stats in alpha before combining with alpha')
        if ('p000',tsink_str) not in alpha_in.NNQ_Stats['Alphaboot']:
            print(alpha_in.NNQ_Stats['Alphaboot'])
            raise AttributeError('alpha ratio needs p000,'+tsink_str+' result')
        alpha_p0 = alpha_in.NNQ_Stats.loc[('p000',tsink_str,slice(None)),'Alphaboot']
        vall,ilist,corrl = [],[],[]
        for (igamma,ip,itflow,it),irat in self.Rat_Stats['boot'].items():
            if itflow not in alpha_p0.index.swaplevel(0,2):
                out_str = 'missmatch in flow times when CombineWithAlpha called \n'
                out_str += itflow
                raise EnvironmentError(out_str)
            ilist.append((igamma,ip,itflow,it))
            output = this_fun(irat,alpha_p0[('p000',tsink_str,itflow)])
            output.Stats()
            vall.append(output)
            corrl.append(irat.corr(alpha_p0[('p000',tsink_str,itflow)]))
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.Rat_Col_Names)
            self.Rat_Stats.loc[:,'alpha_rat_'+this_fun.__name__] = pa.Series(vall,index=indicies)
            self.Rat_Stats.loc[:,'alpha_rat_'+this_fun.__name__+'_corr'] = pa.Series(corrl,index=indicies)
            self.Write()

    def Fit(self,fit_range='PreDef',iGuess='PreDef',EstDir=False,tflowfit=[],WipeFit=True,
            show_timer=False,min_fitr=2):
        self.ImportFitRanges(fit_range)
        if iGuess != 'PreDef': self.iGuess = iGuess
        lFit,ilist = [],[]
        # print 'DEBUG'
        # print self.Rat_Stats['boot']
        for (igamma,ip,itflow),pdata in self.Rat_Stats['boot'].groupby(level=('current_operator','momentum','flow_times')):
            this_key = (igamma,ip,itflow)
            if (igamma,ip,itflow) in self.Rat_Fit_Stats.index:
                if WipeFit:
                    self.Rat_Fit_Stats.drop([(igamma,ip,itflow)],inplace=True)
                else:
                    continue

            # ##debugging
            # print
            # print 'FitData:'
            # for it,iy in zip(tdatarange[0],ydata):
            #     print it, iy
            # print
            # print 'DEBUG2'
            # print ydata
            # print np.array(ydata)
            # print [np.array(tdatarange),np.array(ydata)]
            fit_info = {}
            if EstDir:
                fit_info['Funs'] = (self.Fun,self.npar,'Estimate')
            else:
                fit_info['Funs'] = (self.Fun,self.npar)
            fit_info['iGuess'] = self.iGuess
            fit_info['name'] = '_'.join(this_key)
            thisff = sff.SetOfFitFuns(data=pdata[this_key],name=fit_info['name'])
            thisff.ScanRange_fitr(self.fit_range,fit_info=fit_info,min_fit_len=min_fitr)
            lFit.append(thisff)
            ilist.append(this_key)
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.Rat_Fit_Col_Names)
            if 'boot' in self.Rat_Fit_Stats.columns:
                this_df = pa.DataFrame()
                this_df.loc[:,'boot'] = pa.Series(lFit,index=indicies)
                self.Rat_Fit_Stats = self.Rat_Fit_Stats.append(this_df)
            else:
                self.Rat_Fit_Stats.loc[:,'boot'] = pa.Series(lFit,index=indicies)
            def DoFits(val):
                val.DoFits(show_timer=show_timer)
                return val
            self.Rat_Fit_Stats['boot'] = self.Rat_Fit_Stats['boot'].apply(DoFits)
            self.Write()


    def Write(self):
        def FixDictArray(thislist,ikey):
            outDict = OrderedDict()
            for ic,ilist in enumerate(thislist):
                outDict[ikey+str(ic+1)] = ilist
            return outDict

        if not hasattr(self,'show_cfgs'):
            self.show_cfgs = True

        # C2name = self.C2class.name.replace(self.C2class.jsm,self.C3FOclass.NJN.jsm)
        # C2LL = self.C2class.LegLab.replace(self.C2class.jsm,self.C3FOclass.NJN.jsm)
        # self.C2class.SetCustomName(C2name,C2LL)
        if hasattr(self.C3FOclass,'Write'):
            self.C3FOclass.Write()

        outDict = ODNested()
        outDict['name'] = self.name
        if hasattr(self.C2class,'Write'):
            self.C2class.Write()
            outDict['C2file'] = self.C2class.HumanFile
        outDict['C3FOfile'] = self.C3FOclass.HumanFile
        outDict['kud'] = self.C3FOclass.NJN.kud
        outDict['ks'] = self.C3FOclass.NJN.ks
        outDict['nt'] = self.nt
        outDict['nxyz'] = self.nxyz
        for istream,incfg in zip(self.C3FOclass.NJN.stream_list,self.C3FOclass.NJN.ncfg_list):
            outDict['ncfg_'+istream] = incfg
        outDict['nxsrc_avg'] = self.C3FOclass.NJN.xsrcAvg
        outDict['nMeas'] = self.C3FOclass.NJN.nMeas
        outDict['nboot'] = self.nboot
        outDict['tsrc'] = self.C3FOclass.NJN.tsrc
        outDict['tsink'] = self.C3FOclass.NJN.tsink
        # outDict['xsrc'] = FixDictArray(self.xsrc,'xsrc')
        outDict['ism'] = self.C3FOclass.NJN.ism
        outDict['sink_type'] = self.C3FOclass.NJN.jsm
        outDict['Projector'] = self.C3FOclass.NJN.Proj
        outDict['Doub_Sing'] = self.C3FOclass.DS
        outDict['ppmom'] = self.C3FOclass.NJN.ppmom

        for ict,itflow in enumerate(self.C3FOclass.FO.tflowlist):
            outDict['t_flow'+str(ict+1)] = str(itflow)
        for icg,igamma in enumerate(self.C3FOclass.NJN.gammalist):
            outDict['gamma'+str(icg+1)] = igamma
        for icq,iqmom in enumerate(self.C3FOclass.NJN.qmom):
            outDict['qmom'+str(icq+1)] = iqmom

        excel_params = pa.Series(deepcopy(outDict))
        outDict['tflow_list'] = FixDictArray(self.tflow_list,'tflow')
        outDict['gamma_list'] = FixDictArray(self.C3FOclass.NJN.gammalist,'gamma')
        outDict['qmom_list'] = FixDictArray(self.C3FOclass.NJN.qmom,'qmom')
        for ict,itflow in enumerate(self.C3FOclass.FO.tflowlist):
            del outDict['t_flow'+str(ict+1)]
        for icg,igamma in enumerate(self.C3FOclass.NJN.gammalist):
            del outDict['gamma'+str(icg+1)]
        for icq,iqmom in enumerate(self.C3FOclass.NJN.qmom):
            del outDict['qmom'+str(icq+1)]


        for (igamma,ip,itflow),fitdata in self.Rat_Fit_Stats['boot'].items():
            outDict['Fits'][igamma][ip][itflow] = fitdata.GetOutputDict()

        for col_key,col_data in self.Rat_Stats.items():
            if 'Avg' in col_key or 'Std' in col_key: continue
            for (igamma,ip,it,itflow),idata in col_data.items():
                if hasattr(idata,'Avg'):
                    outDict[col_key][igamma][ip][it][itflow] = AvgStdToFormat(idata.Avg,idata.Std)
                else:
                    outDict[col_key][igamma][ip][it][itflow] = AvgFormat(idata)

        if self.Rat_alpha:
            for (ip,it,itflow),idata in self.SqrtFact.items():
                if np.isnan(idata.Avg): idata.Stats()
                outDict['SqrtFact_G2Q'][itflow][ip][it] = AvgStdToFormat(idata.Avg,idata.Std)
            for (ip,it,itflow),idata in self.C2class.NNQ_Stats['boot'].items():
                if np.isnan(idata.Avg): idata.Stats()
                outDict['G2Q'][itflow][ip][it] = AvgStdToFormat(idata.Avg,idata.Std)
        else:
            y = getattr(self, "SqrtFact", None)
            if y is not None and hasattr(self.SqrtFact,'iteritems'):
                for (ip,it),idata in self.SqrtFact.items():
                    if np.isnan(idata.Avg): idata.Stats()
                    outDict['SqrtFact'][ip][it] = AvgStdToFormat(idata.Avg,idata.Std)
                for (ip,it),idata in self.C2class.C2_Stats['boot'].items():
                    if np.isnan(idata.Avg): idata.Stats()
                    outDict['G2'][ip][it] = AvgStdToFormat(idata.Avg,idata.Std)


        if any(np.array(self.C3FOclass.NJN.ncfg_list) > 0) and self.show_cfgs:
            for ((istream,iccfg),icfg),ixsrc_list in zip(iter(self.C3FOclass.NJNQ_cfgs['configs'].items()),self.C3FOclass.NJNQ_cfgs['xsrc_list'].values):
                outDict['cfglist'][istream][icfg] = 'n_xsrc_'+str(len(ixsrc_list))


        ## Write human readable file with data
        WriteXml(self.HumanFile,{'Results':outDict})
        if any(np.array(self.C3FOclass.NJN.ncfg_list) > 0):
            for ((istream,iccfg),icfg),ixsrc_list in zip(iter(self.C3FOclass.NJNQ_cfgs['configs'].items()),self.C3FOclass.NJNQ_cfgs['xsrc_list'].values):
                outDict['cfglist'][istream][icfg] = 'n_xsrc_'+str(len(ixsrc_list))
        ## TODO, include fitting in the excell file
        WriteExcel(self.ExcelFile,{'Rat_results':deepcopy(self.Rat_Stats)},cfglist=self.C3FOclass.NJNQ_cfgs,params=excel_params)

        ## pickles rest of data for reading
        self.RemoveFuns()
        WritePickle(self.PickleFile,self.__dict__)
        self.GetFuns()

    def GetFuns(self):
        if self.Rat_Fit_Stats is not None and 'boot' in self.Rat_Fit_Stats.columns:
            for ival in self.Rat_Fit_Stats['boot'].values:
                ival.GetFuns()
        if hasattr(self.C2class,'GetFuns'):
            self.C2class.GetFuns()
        if hasattr(self.C2Set,'GetFuns'):
            self.C2Set.GetFuns()
        self.C3FOclass.GetFuns()

    def RemoveFuns(self):
        if self.Rat_Fit_Stats is not None and 'boot' in self.Rat_Fit_Stats.columns:
            for ival in self.Rat_Fit_Stats['boot'].values:
                ival.RemoveFuns()
        if hasattr(self.C2class,'GetFuns'):
            self.C2class.RemoveFuns()
        if hasattr(self.C2Set,'GetFuns'):
            self.C2Set.RemoveFuns()
        self.C3FOclass.RemoveFuns()

    ##################################


    def ReadAndWrite(self,WipeData=True,DefWipe=False,CheckCfgs=False):
        self.Read(WipeData=WipeData,DefWipe=DefWipe)
        self.CreateC2Sink()
        self.Bootstrap(WipeData=WipeData,DefWipe=DefWipe,CheckCfgs=CheckCfgs)
        self.CalcRatio()
        self.Fit()
        if hasattr(self,'alpha_rat') and self.alpha_rat is not False:
            print('This has been hard coded to be false, remove these statements to bring it back')
            # self.DoAlphaCombine(self.alpha_rat)
        self.Write()

    def FixRead(self,readdict):
        readdict['alpha_rat'] = self.alpha_rat
        rewrite = False
        if not isinstance(readdict['C2class'],str):
            if not isinstance(readdict['C2Set'],str):
                readdict['C2Set'] = 'Complete'
                rewrite = True
        return readdict,rewrite

    def LoadPickle(self,DefWipe=False,WipeData = True,CheckCfgs=False,CheckMom=True):
        # print
        # print 'NJNQ'
        # print self.HumanFile
        # print self.PickleFile
        if os.path.isfile(self.PickleFile) and not DefWipe:
            print('Loading Pickle file ' , self.PickleFile)
            loadeddict = ReadPickleWrap(self.PickleFile)
            loadeddict,rewrite = self.FixRead(loadeddict)
            checklist = self.thischecklist
            if CheckCfgs: checklist = checklist + ['cfglist']
            if CheckMom: checklist = checklist + ['checkmomlist']
            if CheckClass(self.__dict__,loadeddict,checklist):
                self.__dict__.update(loadeddict)
                if rewrite: self.Write()
                if Wipe_All_Fits:
                    self.Rat_Fit_Stats = pa.DataFrame()
                elif len(self.Rat_Fit_Stats.values)> 0 and not isinstance(self.Rat_Fit_Stats.values[0][0],sff.SetOfFitFuns):
                    print('Legacy file, fits need to be recomputed')
                    self.Rat_Fit_Stats = pa.DataFrame()
            else:
                if os.path.isfile(self.PickleFile+'.bak'):
                    print()
                    print('using backupfile for '+self.name)
                    self.PickleFile = self.PickleFile+'.bak'
                    self.HumanFile = self.HumanFile+'.bak'
                    self.ExcelFile = self.ExcelFile.replace('.xlsx','.bak.xlsx')
                    self.LoadPickle(DefWipe=DefWipe,WipeData=WipeData,CheckCfgs=CheckCfgs,CheckMom=CheckMom)
                    return
                print('Warning, file ' , self.HumanFile , ' has different parameters than this instance, reading and writing over file now')
                print()
                self.ReadAndWrite(WipeData=WipeData,DefWipe=DefWipe)
            # self.Write()
        else:
            if not os.path.isfile(self.PickleFile):
                print()
                print('pickle file not found: \n', self.PickleFile)
            self.ReadAndWrite(WipeData=WipeData,DefWipe=DefWipe)

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
        if self.thisshiftscale == 'Not Set':
            xlen = np.abs(xlims[1]-xlims[0])
            self.thisshiftscale = self.thisshift*xlen
            return self.thisshiftscale
        else:
            return self.thisshiftscale


    def Importtflow(self,tflowlist):
        if len(tflowlist) == 0:
            self.tflowfit = self.C3FOclass.tflowfit
        else:
            self.tflowfit = tflowlist



    ## make sure all xlims are the same, because this sets the window of xlim (needed for calculating shifts)
    ## gammamomlist is a dictionary where the keys are the gamma matricies, and the values are lists of momenta
    def Plot(self,plot_class,xlims='tsink',gammamomlist={},tflowlist = [],thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine',alpha_fun='Pre_Def'):
        if xlims == 'tsink': xlims = [0,int(self.C3FOclass.NJN.tsink.replace('tsink',''))+1]
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        self.thisshiftscale = 'Not Set'
        thisshift = self.GetShift(xlims,thisshift)
        rs_col_name = 'boot'
        this_leg = self.LegLab
        if alpha_fun is not False:
            if hasattr(alpha_fun,'__name__'):
                rs_col_name = alpha_fun.__name__
            elif alpha_fun == 'Pre_Def':
                rs_col_name = self.alpha_rat
            else:
                rs_col_name = alpha_fun
            if rs_col_name is not False:
                if '_corr' not in rs_col_name:
                    plot_class = self.Plot(plot_class,xlims=xlims,gammamomlist=gammamomlist,tflowlist =tflowlist,
                                           thiscol=thiscol,thissym=thissym,thisshift=thisshift,alpha_fun=rs_col_name+'_corr')
                this_leg += ' '+rs_col_name+r' \alpha'
                rs_col_name = 'alpha_rat_'+rs_col_name
            else:
                rs_col_name = 'boot'

        if rs_col_name not in self.Rat_Stats.columns:
            if 'boot' not in rs_col_name:
                print('no data to plot for column '+rs_col_name+ ', attemping to run alpha computation (THIS MAY TAKE ALONG TIME!)')
                self.DoAlphaCombine(rs_col_name)
            else:
                print('no data to plot for ratio function, skipping')
        if rs_col_name not in self.Rat_Stats.columns:
            raise AttributeError('error trying to compute ratio alpha combintation '+rs_col_name)

        if len(list(gammamomlist.keys())) == 0:
            gammamomlist = self.defgammamomlist
        # for (igamma,ip,itflow),Rat_data in self.Rat_Stats[rs_col_name].groupby(level=('current_operator','momentum','flow_times')):
        #     if itflow not in tflowlist: continue
        #     if igamma not in gammamomlist.keys(): continue
        #     if ip not in gammamomlist[igamma]: continue

        igamma = list(gammamomlist.keys())[0]
        ip = self.latparams.TOpform(list(gammamomlist.values())[0][0],pref='q')
        ip = ip.replace('qq','q')
        itflow = tflowlist[0]
        if (igamma,ip,itflow) not in self.Rat_Stats[rs_col_name].index: return
        this_key = (igamma,ip,itflow,slice(None))
        databoot = self.Rat_Stats[rs_col_name]

        tlist = list(databoot.loc[this_key].keys())
        # tlist,databoot = zip(*list(self.Rat_Stats[rs_col_name][igamma,ip,itflow].iteritems()))
        # dataavg,dataerr = [],[]
        # for idata in databoot[xlims[0]:xlims[1]]:
        #     dataavg.append(idata.Avg)
        #     dataerr.append(idata.Std)
        tlist = tlist[xlims[0]:xlims[1]]
        # print 'plotting ' , igamma,ip
        # pleg,gammaleg = '',''
        # if len(gammamomlist[igamma]) > 1:
        #     pleg = ' $'+ip+'$'
        # if len(gammamomlist.keys()) > 1:
        #     gammaleg = ' $'+igamma+'$'
        # pleg = ' $'+ip+'$'
        # if '_cmplx' in igamma:
        #     gammaleg = ' $i'+igamma.replace('_cmplx','')+'$'
        # else:
        #     gammaleg = ' $'+igamma+'$'
        # tflowleg = ' $'+tflowTOLeg(itflow)+'$'
        tlistplot,self.centert = CenterT(np.array(list(map(untstr,tlist)))-twopttshift+thisshift)
        hold_series = pa.Series()
        hold_series['x_data'] = tlistplot
        if '_corr' in rs_col_name:
            hold_series['y_data'] = databoot
            hold_series['yerr_data'] = None
            hold_series['type'] = 'scatter_vary'
        else:
            hold_series['y_data'] = databoot.apply(lambda x : x.Avg)
            hold_series['yerr_data'] = databoot.apply(lambda x : x.Std)
            hold_series['type'] = 'error_bar_vary'
        hold_series['xerr_data'] = None
        hold_series['key_select'] = this_key
        hold_series['fit_class'] = None
        # hold_series['label'] = self.LegLab+pleg+gammaleg+tflowleg
        hold_series['label'] =  this_leg
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        hold_series['xdatarange'] = None
        hold_series['Phys'] = None
        hold_series['ShowPar'] = None
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class

    def GetBestFitr(self,igamma,ip,itflow):
        if 'boot' not in self.Rat_Fit_Stats.columns or (igamma,ip,itflow) not in self.Rat_Fit_Stats.index:
            errstr = ' '.join(igamma,ip,str(itflow))+' not foun in ratio fit list, cannot find best'
            raise EnvironmentError(errstr)
        chi_str,cut_fit_result = self.Rat_Fit_Stats.at[(igamma,ip,itflow),'boot'].GetBestChi()
        return r'$\chi^{2}_{pdf}={:.3f}$'.format(chi_str),cut_fit_result

    ## state is 1 or 2 state fit
    ## fitr is fit range e.g. fitr5-10
    ## fit range of 'Data' is just the original fit range
    def FitPlot(self,plot_class,fitr,xlims = 'Data',gammamomlist={},
                tflowlist = [],thiscol='PreDefine',thisshift='PreDefine'):
        # if self.alpha_rat is not False:
        #     print 'alpha_rat is not false, skipping fit plot'
        #     return plot_class
        self.CheckCol(thiscol)
        thisshift = self.GetShift(xlims,thisshift)
        self.Importtflow(tflowlist)
        if self.centert == 'Not Set':
            print('print centert not set, defaulting to 0')
            self.centert = 0.0
        if 'boot' not in self.Rat_Fit_Stats.columns:
            raise EnvironmentError('No fitting has been done, maybe try to Fit before plotting fit')
        if len(list(gammamomlist.keys())) == 0:
            gammamomlist = self.defgammamomlist
        # if (igamma,ip,fitr) not in self.Rat_Fit_Stats.index:
        # for igamma,momlist in gammamomlist.iteritems():
        #     for ip in momlist:
        #         for itflow in tflowlist:
        igamma,momlist = next(iter(gammamomlist.items()))
        ip = self.latparams.TOpform(list(gammamomlist.values())[0][0],pref='q')
        ip = ip.replace('qq','q')
        itflow = tflowlist[0]
        if fitr == 'Best':
            chi_legend,fitr = self.GetBestFitr(igamma,ip,itflow)
        if 'boot' not in self.Rat_Fit_Stats.columns:
            print('no fits done, attempting fit now')
            self.Fit(fit_range=fitr)
        fit_data = sff.PullSOFSeries(self.Rat_Fit_Stats['boot'])
        this_fun = fit_data.index[0][3]
        this_key = (igamma,ip,itflow,this_fun,fitr)
        if this_key not in fit_data:
            print(this_key , 'not present in Rat Fit Stats, fitting now')
            self.Fit(fit_range=fitr)
            fit_data = sff.PullSOFSeries(self.Rat_Fit_Stats['boot'])
            this_fun = fit_data.index[0][3]
            this_key = (igamma,ip,itflow,this_fun,fitr)

        hold_series = null_series
        hold_series['key_select'] = this_key
        hold_series['type'] = 'fit_vary'
        hold_series['fit_class'] = fit_data
        # hold_series['label'] = fit_data[this_key].name +' ' + chi_legend
        hold_series['label'] = 'ACTUAL_VALUE_fit'
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        hold_series['xdatarange'] = xlims
        hold_series['ShowPar'] = r'Par1'
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class






    ## ##################### ###############################################


    ## Comparisons

    def __str__(self):
        return '_'.join([self.SubDir,self.name])


    def __iter__(self):
        return self

    def iteritems(self):
        return iter(self.Rat_Stats['boot'].items())

    def items(self):
        return list(self.Rat_Stats['boot'].items())

    def itemsAvgStd(self):
        outlist = []
        for ivals in self.Rat_Stats.itertuples():
            outlist.append((ivals[0],ivals[2],ivals[3]))
        return outlist

    def values(self):
        return self.Rat_Stats['boot'].values

    def itervalues(self):
        return iter(self.Rat_Stats['boot'].values)

    def valuesAvgStd(self):
        outlist = []
        for ivals in self.Rat_Stats.itertuples():
            outlist.append((ivals[2],ivals[3]))
        return outlist

    def keys(self):
        return self.Rat_Stats.index


    def __setitem__(self, key, value ):
        self.Rat_Stats['boot'][key] = value

    def __getitem__(self, key):
        return self.Rat_Stats['boot'][key]

    def __reversed__(self):
        self.Rat_Stats = self.Rat_Stats.iiloc[::-1]
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

    def ReduceIndicies(self):
        indexl = []
        for (igamma,ip,it,itflow) in list(self.Rat_Stats['boot'].keys()):
            first_gamma = MultiSplit(igamma,op_str_list+['FO','FOFull'])[0]
            this_gamma = first_gamma.replace('FO','')
            this_gamma = this_gamma.replace('FOFull','')
            this_gamma = this_gamma.replace('__','_')
            if this_gamma[-1] == '_':
                this_gamma = this_gamma[:-1]
            this_p = MultiSplit(ip,op_str_list+['FO','FOFull','_'])[0]
            indexl.append((this_gamma,this_p,it,itflow))
        indicies = pa.MultiIndex.from_tuples(indexl,names=self.Rat_Col_Names)
        this_series = pa.Series(self.Rat_Stats['boot'].values,index=indicies)
        self.Rat_Stats = pa.DataFrame(columns=['boot'],index=indicies)
        self.Rat_Stats.loc[:,'boot'] = this_series
        self.Stats()


    def sin(self):
        result = RatioFOCorr(thisnboot=self.nboot, cfglist=self.cfglist[['configs','xsrc_list']],
                         Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift,
                         name = 'sin('+self.name+')',filename='sin_'+self.filename)
        result.Rat_Stats.loc[:,'boot'] = np.sin(self.Rat_Stats.loc[:,'boot'])
        return result

    def cos(self):
        result = RatioFOCorr(thisnboot=self.nboot, cfglist=self.cfglist[['configs','xsrc_list']],
                         Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift,
                         name = 'cos('+self.name+')',filename='cos_'+self.filename)
        result.Rat_Stats = pa.DataFrame(columns = ['boot'],index=self.Rat_Stats.index)
        result.Rat_Stats.loc[:,'boot'] = np.cos(self.Rat_Stats.loc[:,'boot'])
        return result


    def Overload_wrap(self,Rat2,this_fun):
        ## TODO, include some name and filename formatting
        thisname,thisfilename = self.name,self.filename
        if hasattr(Rat2,'name'):
            thisname += '_'+this_fun.__name__+'_'+Rat2.name
        if hasattr(Rat2,'filename'):
            thisfilename += '_'+this_fun.__name__+'_'+Rat2.filename
        result = RatioFOCorr(thisnboot=self.nboot, cfglist=self.cfglist[['configs','xsrc_list']],
                         Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.Rat_Stats = pa.DataFrame(columns = ['boot'])
        # result.UpdateName('+',self,Rat2)
        if 'boot' in self.Rat_Stats:
            if isinstance(Rat2,RatioFOCorr):
                # result.Rat_Stats.loc[:,'boot'] = this_fun(self.Rat_Stats.loc[:,'boot'],Rat2.Rat_Stats.loc[:,'boot'])
                if 'boot' in Rat2.Rat_Stats:
                    vall,ilist = [],[]
                    for (igamma,ip,itflow,it),irat in self.Rat_Stats['boot'].items():
                        if (it,itflow) not in Rat2.Rat_Stats['boot'].index.swaplevel(0,3).swaplevel(1,2):
                            out_str = 'missmatch in euclidian times / flow times when combining ratio functions \n'
                            out_str += str(self) + '\n'
                            out_str += this_fun.__name__+'\n'
                            out_str += str(Rat2)
                            raise EnvironmentError(out_str)
                        for (igamma2,ip2,itflowdump,itdump),irat2 in Rat2.Rat_Stats.loc[(slice(None),slice(None),itflow,it),'boot'].items():
                            igammacomb,ipcomb = igamma,ip
                            if igamma2 not in igamma:
                                igammacomb += '_'+this_fun.__name__+'_'+igamma2
                            if ip2 not in ip:
                                ipcomb += '_'+this_fun.__name__+'_'+ip2
                            ilist.append((igammacomb,ipcomb,itflow,it))
                            vall.append(this_fun(irat,irat2))
                    if len(ilist) > 0:
                        indicies = pa.MultiIndex.from_tuples(ilist,names=self.Rat_Col_Names)
                        result.Rat_Stats.loc[:,'boot'] = pa.Series(vall,index=indicies)
                else:
                    raise EnvironmentError('Run Bootstrap on second class before combining ratio functions')
            elif isinstance(Rat2,RatioCorr):
                # result.Rat_Stats.loc[:,'boot'] = this_fun(self.Rat_Stats.loc[:,'boot'],Rat2.Rat_Stats.loc[:,'boot'])
                if 'boot' in Rat2.Rat_Stats:
                    vall,ilist = [],[]
                    for (igamma,ip,itflow,it),irat in self.Rat_Stats['boot'].items():
                        if it not in Rat2.Rat_Stats['boot'].index.swaplevel(0,2):
                            out_str = 'missmatch in euclidian times when combining ratio functions \n'
                            out_str += str(self) + '\n'
                            out_str += this_fun.__name__+'\n'
                            out_str += str(Rat2)+'\n'
                            out_str += it
                            raise EnvironmentError(out_str)
                        for (igamma2,ip2,itdump),irat2 in Rat2.Rat_Stats.loc[(slice(None),slice(None),it),'boot'].items():
                            igammacomb,ipcomb = igamma+'_FO_',ip
                            # if igamma != igamma2:
                            igammacomb += this_fun.__name__+'_'+igamma2
                            if ip2 not in ip:
                                ipcomb += '_'+this_fun.__name__+'_'+ip2
                            ilist.append((igammacomb,ipcomb,itflow,it))
                            vall.append(this_fun(irat,irat2))
                    if len(ilist) > 0:
                        indicies = pa.MultiIndex.from_tuples(ilist,names=self.Rat_Col_Names)
                        result.Rat_Stats.loc[:,'boot'] = pa.Series(vall,index=indicies)
                else:
                    raise EnvironmentError('Run Bootstrap on second class before combining ratio functions')
            else:
                try:
                    result.Rat_Stats.loc[:,'boot'] = this_fun(self.Rat_Stats.loc[:,'boot'], Rat2)
                except Exception as err:
                    print(type(Rat2), Rat2)
                    print(self.Rat_Stats)
                    raise EnvironmentError(str(err)+'\nInvalid value to combine with RatioFOCorr class')
            return result
        else:
            raise EnvironmentError('Run Bootstrap on first class before combining ratio functions')

    def __add__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOCorr.parentlist]):
            return Rat2.__radd__(self)
        else:
            return self.Overload_wrap(Rat2,op.add)

    def __sub__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOCorr.parentlist]):
            return Rat2.__rsub__(self)
        else:
            return self.Overload_wrap(Rat2,op.sub)

    def __mul__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOCorr.parentlist]):
            return Rat2.__rmul__(self)
        else:
            return self.Overload_wrap(Rat2,op.mul)

    def __div__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOCorr.parentlist]):
            return Rat2.__rdiv__(self)
        else:
            return self.Overload_wrap(Rat2,op.truediv)

    def __pow__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOCorr.parentlist]):
            return Rat2.__rpow__(self)
        else:
            return self.Overload_wrap(Rat2,op.pow)


    ## Right multiplication functions


    def __radd__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOCorr.parentlist]):
            return Rat2.__add__(self)
        else:
            return self.Overload_wrap(Rat2,flip_fun_2arg(op.add))

    def __rsub__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOCorr.parentlist]):
            return Rat2.__sub__(self)
        else:
            return self.Overload_wrap(Rat2,flip_fun_2arg(op.sub))

    def __rmul__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOCorr.parentlist]):
            return Rat2.__mul__(self)
        else:
            return self.Overload_wrap(Rat2,flip_fun_2arg(op.mul))

    def __rdiv__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOCorr.parentlist]):
            return Rat2.__div__(self)
        else:
            return self.Overload_wrap(Rat2,flip_fun_2arg(op.truediv))


    def __rpow__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOCorr.parentlist]):
            return Rat2.__pow__(self)
        else:
            return self.Overload_wrap(Rat2,flip_fun_2arg(op.pow))



####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################

##UP TO HERE TODO:

class RatioFOFullCorr(object):
    """

    Rat(q,tau,pp,t,tflow) uses bootstrap class
    pp  = sink momentum
    t   = source sink separation time
    q   = current insersion momentum
    tau = current insersion time
    tflow = flow time
    tsum = summed current flow time.

    Info is a dictionary containing information for the correlator:
    see comments below to see what is required.
    missing elements just sets field to default


    """
    def Construct_RatFOFull_File(this_file):
        return Construct_File_Object(this_file,RatioFOFullCorr)


    ## TODO: create variational flag and interface with this function.

    ## array of all classes that inherit this class. Append this when adding new classes please!
    parentlist = ['SetOfCorrs.SetsOfRatios']

    def __init__(self, thisnboot=nboot, cfglist={}, Info={},name='',
                 filename='',thissym='Not Set',thiscol='Not Set',thisshift=0.0
                 ,man_load_cfgs=False):


        self.current = 0
        ## Initialising plotting parameters
        self.thisshift = thisshift
        self.thiscol = thiscol
        self.thissym = thissym
        ## the actual shift ammount needs to be scaled by the x axis size
        ## see GetShift
        self.thisshiftscale = 'Not Set'

        self.C3FOclass = NJNQFullCorr(thisnboot=thisnboot, cfglist=cfglist, Info=Info,name=name,filename=filename)
        self.C3file = self.C3FOclass.PickleFile
        self.tflow_list = self.C3FOclass.FO.tflowlist
        self.thischecklist = ['C3file','Sym','boot_tsum_list','tflow_list','Fun']
        self.checkmomlist = self.C3FOclass.NJN.checkmomlist
        self.boot_tsum_list = self.C3FOclass.boot_tsum_list
        self.latparams = self.C3FOclass.NJN.latparams
        Info['pmom'] = [self.latparams.GetAvgmomstr(ip) for ip in self.C3FOclass.NJN.qmom]
        if '0 0 0' not in Info['pmom']:
            Info['pmom'] = ['0 0 0'] + Info['pmom']
        self.nt = self.latparams.nt
        self.nxyz = self.latparams.nxyz
        self.dim_label = 'RC'+str(self.nxyz)+'x'+str(self.nt)
        self.dim_ll = str(self.nxyz)+r'^{3}\times ' + str(self.nt)

        ## both current and sink momentum are requred for the two point correlator construction when constructing the ratio function
        ## but the two point correlator has been averaged over, so we need to get the averaged momentum anaolgy
        ## bah, this was too hard to keep unique files for the two-point correlator.
        ## easier just to read and write the whole 2 point mom file.
        Info['pmom'] = self.latparams.GetMomList(list_type='str_Avg')
        ## e.g. q = [-1,1,-1] corresponds to p = [1 1 1] in the two-point correlator
        # ppforp = self.latparams.GetAvgmomstr(Info['ppmom'])
        # if ppforp not in Info['pmom']: Info['pmom'].append(ppforp)


        if 'show_cfgs' in list(Info.keys()): self.show_cfgs = Info['show_cfgs']
        else: self.show_cfgs = True


        ## Does symmetrising over smearings
        if 'Symmetric' in list(Info.keys()): self.Sym = Info['Symmetric']
        else: self.Sym = False ## default to no symmetrizing, (usually full matrix isn't computed nowerdays.

        self.sinklab,self.jsmcoeffs,self.sm2ptList,InfoList = GetInfoList(deepcopy(Info))




        # self.C2filelist = [ikey.PickleFile for ikey in self.C2classlist.itervalues()]
        if self.C3FOclass.FO.Observ == '':
            for iinfo in InfoList: iinfo['fo_for_cfgs'] = False
        else:
            for iinfo in InfoList: iinfo['fo_for_cfgs'] = self.C3FOclass.FO.Observ

        self.C2Set = SetOfTwoPt(cfglist=cfglist,InfoDict = InfoList)
        self.C2class = 'Not Set'

        self.defgammamomlist = OrderedDict()
        for igamma in self.C3FOclass.NJN.gammalist:
            self.defgammamomlist[igamma] = self.C3FOclass.NJN.qform

        self.Rat_Col_Names = [  'current_operator','momentum','flow_times',
                                'source_sink_separation','flow_op_trange']
        self.Rat_Stats = pa.DataFrame()

        '''
        self.Rat_Stats DataFrame:

        columns =   boot,   Avg,        Std

        rows =      current_operator,   momenta,    flow_times,
                    time separation,    flow op trange
        '''

        self.RatInt_Col_Names = [  'current_operator','momentum','flow_times',
                                'source_sink_separation']
        self.RatInt_Stats = pa.DataFrame()

        '''
        self.RatInt_Stats DataFrame:

        columns =   boot,   Avg,        Std

        rows =      current_operator,   momenta,    flow_times,
                    time separation
        '''


        self.Rat_Fit_Col_Names = [  'current_operator','momentum','flow_times']
        self.Rat_Fit_Stats = pa.DataFrame()

        '''
        self.Rat_Fit_Stats DataFrame:

        columns =   boot,   Avg,        Std

        rows =      current_operator, momenta,  flow_times, time fit range, flow op trange

        Avg/Std elements are pa.Series instances with:
            index = fit parameters
            values = Avg/Std values
        '''


        self.RatInt_Fit_Col_Names = [  'current_operator','momentum','flow_times' ]
        self.RatInt_Fit_Stats = pa.DataFrame()

        '''
        self.RatInt_Fit_Stats DataFrame:

        columns =   boot,   Avg,        Std

        rows =      current_operator, momenta,  flow_times, time fit range

        Avg/Std elements are pa.Series instances with:
            index = fit parameters
            values = Avg/Std values
        '''



        self.nboot = thisnboot
        self.centert = 'Not Set'


        ## Info can contain fit ranges to be calculated
        ## must bein form Info['FitsRF'] = [ ('taufit#-#','fit#-#') , ...]
        if 'FitsRFFull' in list(Info.keys()): self.PredefFits = Info['FitsRFFull']
        else: self.PredefFits = []


        if 'RFFitFun' in list(Info.keys()): self.Fun,self.npar = Info['RFFitFun']
        # else: self.Fun,self.npar = (RF_Prop_Fit,3)
        else: self.Fun,self.npar = (ConstantFitFun,1)

        ## Info can contain initial guess for fits, make sure it array of correct length
        if 'iGuessRF' in list(Info.keys()): self.iGuess = np.array(Info['iGuessRF'])
        else: self.iGuess = 'None' ## Defaults to value set for function in FitFunctions.py

        if 'min_fit_len_tsum' in list(Info.keys()): self.min_fit_len_tsum = Info['min_fit_len_tsum']
        else: self.min_fit_len_tsum = 'Default'


        ## is needed?
        self.Info = Info

        self.SetCustomName(string=name,stringfn='')
        if not man_load_cfgs:
            self.LoadCfgs(cfglist,name,filename)

    def LoadCfgs(self,cfglist,name='',filename=''):
        if len(cfglist) == 0:
            self.CombineCfgLists()
        else:
            self.ImportCfgList(cfglist)
        self.SetCustomName(string=name,stringfn='')

    ## just for consistancy with other classes
    def GetCfgList(self):
        self.CombineCfgLists()

    def CombineCfgLists(self):
        if isinstance(self.C2class,str):
            thiscfglist = CombineCfgs(self.C3FOclass.NJNQFull_cfgs,self.C2Set.cfglist)
            # print
            # print len(thiscfglist.keys()), len(self.C3FOclass.cfglist.keys()), len(self.C2Set.cfglist.keys())
        else:
            thiscfglist = CombineCfgs(self.C3FOclass.NJNQFull_cfgs,self.C2class.C2_cfgs)
            # print
            # print len(thiscfglist.keys()), len(self.C3FOclass.cfglist.keys()), len(self.C2class.cfglist.keys())
        self.ImportCfgList(thiscfglist)




    def SetCustomName(self,string='',stringfn='',stringLL=''):
        if string == '':
            self.name = 'Rat_'+self.C3FOclass.name
        else:
            self.name = string
        if stringfn == '':
            self.filename = 'Rat_'+self.C3FOclass.filename
        else:
            self.filename = stringfn
        if stringLL == '':
            thisops = self.C3FOclass.FO.Observ
            self.LegLab = '$'+'\ '.join([   self.dim_ll,self.latparams.GetPionMassLab(),self.C3FOclass.NJN.stream_name,thisops,
                                            self.C3FOclass.DS,self.C3FOclass.NJN.Proj]) \
                                            +self.C3FOclass.tr_sum_name.replace('_','\ ')+'$' ## customise this how you like
        else:
            self.LegLab = stringLL
        self.SubDir = '_'.join([self.C3FOclass.NJN.Proj,self.C3FOclass.DS,self.C3FOclass.NJN.ppform])
        self.Ratdir =   outputdir + '/'+self.dim_label+self.C3FOclass.NJN.kappafolder+'/Rat/'+ \
                        self.C3FOclass.FO.Observ+'/'+self.SubDir+'/'
        mkdir_p(self.Ratdir+'/Pickle/')
        mkdir_p(self.Ratdir+'/Excel/')
        self.HumanFile = self.Ratdir+self.filename+'.xml'
        self.ExcelFile = self.Ratdir+'/Excel/'+self.filename+'.xlsx'
        self.PickleFile = self.Ratdir+'/Pickle/'+self.filename+'.py3p'



    def ImportCfgList(self,cfglist):
        self.cfglist = cfglist
        self.C2Set.ImportCfgList(cfglist)
        if self.C2class != 'Not Set': self.C2class.ImportCfgList(cfglist)
        self.C3FOclass.ImportCfgList(cfglist)

    # def ImportCfgList(self,cfglist):
    #     for iC2 in self.C2classlist.itervalues():
    #         iC2.ImportCfgList(cfglist)
    #     if self.C2class != 'Not Set': self.C2class.ImportCfgList(cfglist)
    #     self.C3FOclass.ImportCfgList(cfglist)


    def Stats(self):
        if hasattr(self.C2class,'Stats'):
            self.C2class.Stats()
        if hasattr(self.C3FOclass,'Stats'):
            self.C3FOclass.Stats()
        if 'boot' not in self.Rat_Stats: return
        lAvg,lStd = [],[]
        for (igamma,ip,it,itflow,its),iRat in self.items():
            iRat.Stats()
            lAvg.append(iRat.Avg)
            lStd.append(iRat.Std)
        indicies = pa.MultiIndex.from_tuples(self.Rat_Stats.index,names=self.Rat_Col_Names)
        self.Rat_Stats.loc[:,'Avg'] = pa.Series(lAvg,index=indicies)
        self.Rat_Stats.loc[:,'Std'] = pa.Series(lStd,index=indicies)


    def RemoveVals(self):
        # for iC2 in self.C2classlist.itervalues():
        #     iC2.RemoveVals()
        self.C2class.RemoveVals()
        self.C3FOclass.RemoveVals()

    ## self.Rat = [ igamma , ip, it , iconf ]
    ## configs must be the same, so we check for them
    def Read(self,WipeData=True,DefWipe=False):
        if isinstance(self.C2class,str):
            self.C2Set.LoadPickle(DefWipe=DefWipe,CheckCfgs=True,NoFit=True)
        # else:
        #     print('WARNING: C2class (combined from C2Set) is already present, yet we are trying to read')
        self.C3FOclass.LoadPickle(WipeData=WipeData,DefWipe=DefWipe,CheckCfgs=True,CheckMom=True)

    def CreateC2Sink(self):
        if isinstance(self.C2class,str):
            thisname,thisleglab,jsmlab = CreateNameAndLL(self.C2Set.SetC2,self.sinklab)
            thisfun = GenSinkFun(self.jsmcoeffs)
            self.C2class = self.C2Set.CombineCorrs(thisfun,filename=thisname,LegLab=thisleglab,jsmlab=jsmlab)
            self.C2Set = 'Complete'

    def Bootstrap(self,WipeData=True,DefWipe=False,CheckCfgs=False):
        if isinstance(self.C2class,str) :
            self.Read(WipeData=WipeData,DefWipe=DefWipe)
            self.CreateC2Sink()
        # self.C2class.Bootstrap()
        self.C3FOclass.Bootstrap(WipeData=WipeData,DefWipe=DefWipe)


    def CalcRatio(self):
        self.CalcSqrtFac()
        ## if DS has a divide in it, we ignore the sqrt factor as they cancel
        div_bool = 'divide' in self.C3FOclass.DS
        ## double the square root factor if we multiply two ds types
        mult_bool = 'multiply' in self.C3FOclass.DS
        tsinkint = int(self.C3FOclass.tsink.replace('tsink',''))
        lRat,lRatAvg,lRatStd = [],[],[]
        ilist = []
        for (igamma,ip,it,itflow,its),tC3FO in self.C3FOclass.NJNQFull_Stats['boot'].items():
            ip2pt = self.latparams.GetAvgmomform(ip,pref='p')
            ict = int(untstr(it))
            if ict > tsinkint: continue
            if (ip2pt,it) not in self.SqrtFact.index: continue
            # warnings.warn('DEBUG, please check if tsink for 2-pt corr is correct via g4.')
            # print igamma, ip, it, self.C3class.tsink, tC3.Avg,self.SqrtFact[ip2pt,it].Avg

            if div_bool:
                this_boot = tC3FO
            elif mult_bool:
                this_boot = tC3FO*(self.SqrtFact[ip2pt,it]**2)
            else:
                this_boot = tC3FO*self.SqrtFact[ip2pt,it]
            this_boot.Stats()
            ilist.append((igamma,ip,itflow,it,its))
            lRat.append(this_boot)
            lRatAvg.append(this_boot.Avg)
            lRatStd.append(this_boot.Std)
            # self.Ratboot[igamma][ip][it] = tC3/self.C2class[ip.replace('q','p'),'t'+str(int(self.C3class.tsink.replace('tsink',''))+1)]
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.Rat_Col_Names)
            self.Rat_Stats.loc[:,'boot'] = pa.Series(lRat,index=indicies)
            self.Rat_Stats.loc[:,'Avg'] = pa.Series(lRatAvg,index=indicies)
            self.Rat_Stats.loc[:,'Std'] = pa.Series(lRatStd,index=indicies)

        lRat,lRatAvg,lRatStd = [],[],[]
        ilist = []
        for (igamma,ip,it,itflow),tC3FO in self.C3FOclass.NJNQInt_Stats['boot'].items():
            ip2pt = self.latparams.GetAvgmomform(ip,pref='p')
            ict = int(untstr(it))
            if ict > tsinkint: continue
            if (ip2pt,it) not in self.SqrtFact.index: continue
            # warnings.warn('DEBUG, please check if tsink for 2-pt corr is correct via g4.')
            # print igamma, ip, it, self.C3class.tsink, tC3.Avg,self.SqrtFact[ip2pt,it].Avg

            this_boot = tC3FO*self.SqrtFact[ip2pt,it]
            this_boot.Stats()
            ilist.append((igamma,ip,itflow,it))
            lRat.append(this_boot)
            lRatAvg.append(this_boot.Avg)
            lRatStd.append(this_boot.Std)
            # self.Ratboot[igamma][ip][it] = tC3/self.C2class[ip.replace('q','p'),'t'+str(int(self.C3class.tsink.replace('tsink',''))+1)]
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.RatInt_Col_Names)
            self.RatInt_Stats.loc[:,'boot'] = pa.Series(lRat,index=indicies)
            self.RatInt_Stats.loc[:,'Avg'] = pa.Series(lRatAvg,index=indicies)
            self.RatInt_Stats.loc[:,'Std'] = pa.Series(lRatStd,index=indicies)



    def CalcSqrtFac(self):
        SqrtFact_hold = []
        ilist = []
        ipp = self.latparams.GetAvgmomform(self.C3FOclass.NJN.ppform,pref='p')
        if Debug_Mode:
            print()
            print('Outputting sqrt factor debug info:')
            print()
        for iq in self.C3FOclass.NJN.qform:
            ip = self.latparams.GetAvgmomform(iq,pref='p')
            sqrt_fac = self.CalcSqrtFacMom(self.C2class.C2_Stats['boot'][ipp],
                                           self.C2class.C2_Stats['boot'][ip],self.C3FOclass.NJN.tsink)
            for it,tdata in enumerate(sqrt_fac):
                if (ip,tstr(it+1)) not in ilist:
                    if Debug_Mode:
                        for ierr in tdata.GetNaNIndex():
                            print(iq, tstr(it+1), ' iboot'+str(ierr).zfill(4))
                    SqrtFact_hold.append(tdata)
                    ilist.append((ip,tstr(it+1)))
            # if neg_sqrt: warnings.warn('negative square root in square root factor for '+ipp+iq)
        indicies = pa.MultiIndex.from_tuples(ilist,names=['momentum','source_sink_separation'])
        self.SqrtFact = pa.Series(SqrtFact_hold,index=indicies)


    def CheckNegForSqrt(self,thisboot,negsqrt):
        for ict in range(len(thisboot.bootvals)):
            if thisboot.bootvals[ict] < 0.0:
                negsqrt = True
                # if G2_abs:
                #     thisboot.bootvals[ict] = abs(thisboot.bootvals[ict])
                # else:
                thisboot.bootvals[ict] = float('NaN')
                # thisboot.bootvals[ict] = 0.0
        return thisboot,negsqrt


    ## taking squareroot of each element before combinding makes the numbers not too small
    def CalcSqrtFacMom(self,data2ptpp,data2ptp,thistsink):
        dictout = []
        thistsink = thistsink.replace('tsink','t')
        thistsinkp1 = 't'+str(int(thistsink.replace('t',''))+twopttshift)
        inttsink = int(thistsink.replace('t',''))
        # NegSQRTRatFac = False
        # data2ptpp[thistsinkp1],NegSQRTRatFac = self.CheckNegForSqrt(data2ptpp[thistsinkp1],NegSQRTRatFac)
        # data2ptp[thistsinkp1],NegSQRTRatFac = self.CheckNegForSqrt(data2ptp[thistsinkp1],NegSQRTRatFac)

        data2ptpp = data2ptpp.apply(lambda x : x.Sqrt())
        data2ptp = data2ptp.apply(lambda x : x.Sqrt())
        # NegSQRTRatFac = any(data2ptpp.isnull().values) or any(data2ptp.isnull().values)
        two = data2ptpp[thistsinkp1]*data2ptp[thistsinkp1]
        # two = two.Sqrt()
        thisshift = twopttshift
        for it in range(thisshift,inttsink+thisshift):
            thiststr = tstr(it)
            tflip = tstr(inttsink-it+twopttshift+thisshift)
            # data2ptpp[thiststr],NegSQRTRatFac = self.CheckNegForSqrt(data2ptpp[thiststr],NegSQRTRatFac)
            # data2ptp[thiststr],NegSQRTRatFac = self.CheckNegForSqrt(data2ptp[thiststr],NegSQRTRatFac)
            # data2ptpp[tflip],NegSQRTRatFac = self.CheckNegForSqrt(data2ptpp[tflip],NegSQRTRatFac)
            # data2ptp[tflip],NegSQRTRatFac = self.CheckNegForSqrt(data2ptp[tflip],NegSQRTRatFac)
            one = data2ptpp[thiststr]/data2ptp[thiststr]
            three = data2ptp[tflip]/data2ptpp[tflip]
            inter = one*three
            ott = inter/two
            # ott.Sqrt()
            ott.Stats()
            dictout.append(ott)
        if Debug_Mode:
            print('Printing two poitn correlators to be combine:')
            for ic,(ip,ipp,icomb) in enumerate(zip(data2ptp,data2ptpp,dictout)):
                icomb.Stats()
                print(ic,ip.Avg,ipp.Avg,icomb.Avg)
        return dictout

    def ImportFitRanges(self,tau_fit_range,tsum_fit_range):
        if not isinstance(tau_fit_range,str):
            print('tau_fit_range for setsoffits takes 1 fit range to vary over, selecting first one')
            tau_fit_range = tau_fit_range[0]
        if not isinstance(tsum_fit_range,str):
            print('tsum_fit_range for setsoffits takes 1 fit range to vary over, selecting first one')
            tsum_fit_range = tsum_fit_range[0]
        if 'PreDef' == tau_fit_range:
            tau_fit_range = self.PredefFits[0][0]
        if 'PreDef' == tsum_fit_range:
            tsum_fit_range = self.PredefFits[0][1]
        self.fit_range = [tau_fit_range,tsum_fit_range]


    def Fit(self,fit_range='PreDef',tsum_ranges='PreDef',iGuess='PreDef',
            EstDir=False,WipeFit=True,show_timer=False,min_fitr=2):
        self.ImportFitRanges(fit_range,tsum_ranges)
        if iGuess != 'PreDef': self.iGuess = iGuess
        lFit,ilist = [],[]
        itaufit,isumfit = self.fit_range
        if self.min_fit_len_tsum == 'Default':
            self.min_fit_len_tsum = unxmlfitr_int(isumfit.replace('tsumfitr','fitr'))
            self.min_fit_len_tsum = max(0,self.min_fit_len_tsum[1] - self.min_fit_len_tsum[0]-3)
        for (igamma,ip,itflow),pdata in self.Rat_Stats['boot'].groupby(level=('current_operator','momentum','flow_times')):
            this_key = (igamma,ip,itflow)
            if (igamma,ip,itflow) in self.Rat_Fit_Stats.index:
                if WipeFit:
                    self.Rat_Fit_Stats.drop([(igamma,ip,itflow)],inplace=True)
                else:
                    continue
            fit_info = {}
            if EstDir:
                fit_info['Funs'] = (self.Fun,self.npar,'Estimate')
            else:
                fit_info['Funs'] = (self.Fun,self.npar)
            fit_info['iGuess'] = self.iGuess
            fit_info['name'] = '_'.join(this_key)
            thisff = sff.SetOfFitFuns(data=pdata[this_key],name=fit_info['name'])
            thisff.ScanBox_fitr(itaufit,isumfit,fit_info=fit_info,
                                min_fit_len_1=min_fitr,min_fit_len_2=self.min_fit_len_tsum)
            lFit.append(thisff)
            ilist.append(this_key)
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.Rat_Fit_Col_Names)
            if 'boot' in self.Rat_Fit_Stats.columns:
                this_df = pa.DataFrame()
                this_df.loc[:,'boot'] = pa.Series(lFit,index=indicies)
                self.Rat_Fit_Stats = self.Rat_Fit_Stats.append(this_df)
            else:
                self.Rat_Fit_Stats.loc[:,'boot'] = pa.Series(lFit,index=indicies)
            def DoFits(val):
                val.DoFits(show_timer=show_timer)
                return val
            self.Rat_Fit_Stats['boot'] = self.Rat_Fit_Stats['boot'].apply(DoFits)
            self.Write()

        # lFit,lFitA,lFitS = [],[],[]
        # ilist = []
        # for (igamma,ip,itflow),pdata in self.RatInt_Stats['boot'].groupby(level=('current_operator','momentum','flow_times')):
        #     for ifit in self.fit_range:
        #         ifitmin,ifitmax = np.array(unxmlfitr(ifit),dtype=np.float64)
        #         tdatarange,ydata = [[]],[]
        #         for it,tdata in pdata.iteritems():
        #             tint = untstr(it[3])
        #             if ifitmin <= tint <= ifitmax:
        #                 tdatarange[0].append(tint)
        #                 ydata.append(tdata)
        #
        #         # ##debugging
        #         # print
        #         # print 'FitData:'
        #         # for it,iy in zip(tdatarange[0],ydata):
        #         #     print it, iy
        #         # print
        #
        #         if EstDir:
        #             thisff = ff.Fitting(Funs=[self.Fun,self.npar,'Estimate'],data=[np.array(tdatarange),np.array(ydata)],name=ifit,iGuess=self.iGuess)
        #         else:
        #             thisff = ff.Fitting(Funs=[self.Fun,self.npar],data=[np.array(tdatarange),np.array(ydata)],name=ifit,iGuess=self.iGuess)
        #         thisff.FitBoots()
        #         ilist.append((igamma,ip,itflow,ifit))
        #         lFit.append(thisff)
        #         lFitA.append(thisff.fit_data['ParamsAvg'])
        #         lFitS.append(thisff.fit_data['ParamsStd'])
        # indicies = pa.MultiIndex.from_tuples(ilist,names=self.RatInt_Fit_Col_Names)
        # if 'boot' in self.RatInt_Fit_Stats.columns:
        #     this_df = pa.DataFrame()
        #     this_df.loc[:,'boot'] = pa.Series(lFit,index=indicies)
        #     this_df.loc[:,'Avg'] = pa.Series(lFitA,index=indicies)
        #     this_df.loc[:,'Std'] = pa.Series(lFitS,index=indicies)
        #     self.RatInt_Fit_Stats = self.RatInt_Fit_Stats.append(this_df)
        # else:
        #     self.RatInt_Fit_Stats.loc[:,'boot'] = pa.Series(lFit,index=indicies)
        #     self.RatInt_Fit_Stats.loc[:,'Avg'] = pa.Series(lFitA,index=indicies)
        #     self.RatInt_Fit_Stats.loc[:,'Std'] = pa.Series(lFitS,index=indicies)


    def Write(self):
        def FixDictArray(thislist,ikey):
            outDict = OrderedDict()
            for ic,ilist in enumerate(thislist):
                outDict[ikey+str(ic+1)] = ilist
            return outDict

        if not hasattr(self,'show_cfgs'):
            self.show_cfgs = True

        # C2name = self.C2class.name.replace(self.C2class.jsm,self.C3FOclass.NJN.jsm)
        # C2LL = self.C2class.LegLab.replace(self.C2class.jsm,self.C3FOclass.NJN.jsm)
        # self.C2class.SetCustomName(C2name,C2LL)
        if hasattr(self.C3FOclass,'Write'):
            self.C3FOclass.Write()

        outDict = ODNested()
        outDict['name'] = self.name
        if hasattr(self.C2class,'Write'):
            self.C2class.Write()
            outDict['C2file'] = self.C2class.HumanFile
        outDict['C3FOfile'] = self.C3FOclass.HumanFile
        outDict['kud'] = self.C3FOclass.NJN.kud
        outDict['ks'] = self.C3FOclass.NJN.ks
        outDict['nt'] = self.nt
        outDict['nxyz'] = self.nxyz
        for istream,incfg in zip(self.C3FOclass.NJN.stream_list,self.C3FOclass.NJN.ncfg_list):
            outDict['ncfg_'+istream] = incfg
        outDict['nxsrc_avg'] = self.C3FOclass.NJN.xsrcAvg
        outDict['nMeas'] = self.C3FOclass.NJN.nMeas
        outDict['nboot'] = self.nboot
        outDict['tsrc'] = self.C3FOclass.NJN.tsrc
        outDict['tsink'] = self.C3FOclass.NJN.tsink
        # outDict['xsrc'] = FixDictArray(self.xsrc,'xsrc')
        outDict['ism'] = self.C3FOclass.NJN.ism
        outDict['sink_type'] = self.C3FOclass.NJN.jsm
        outDict['Projector'] = self.C3FOclass.NJN.Proj
        outDict['Doub_Sing'] = self.C3FOclass.DS
        outDict['ppmom'] = self.C3FOclass.NJN.ppmom

        for ict,itflow in enumerate(self.C3FOclass.FO.tflowlist):
            outDict['t_flow'+str(ict+1)] = str(itflow)
        for icg,igamma in enumerate(self.C3FOclass.NJN.gammalist):
            outDict['gamma'+str(icg+1)] = igamma
        for icq,iqmom in enumerate(self.C3FOclass.NJN.qmom):
            outDict['qmom'+str(icq+1)] = iqmom

        excel_params = pa.Series(deepcopy(outDict))
        outDict['tflow_list'] = FixDictArray(self.C3FOclass.FO.tflowlist,'tflow')
        outDict['gamma_list'] = FixDictArray(self.C3FOclass.NJN.gammalist,'gamma')
        outDict['qmom_list'] = FixDictArray(self.C3FOclass.NJN.qmom,'qmom')
        for ict,itflow in enumerate(self.C3FOclass.FO.tflowlist):
            del outDict['t_flow'+str(ict+1)]
        for icg,igamma in enumerate(self.C3FOclass.NJN.gammalist):
            del outDict['gamma'+str(icg+1)]
        for icq,iqmom in enumerate(self.C3FOclass.NJN.qmom):
            del outDict['qmom'+str(icq+1)]

        if 'boot' in self.Rat_Fit_Stats.columns:
            for (igamma,ip,itflow),fitdata in self.Rat_Fit_Stats['boot'].items():
                outDict['Fits'][igamma][ip][itflow] = fitdata.GetOutputDict()

        if 'boot' in self.RatInt_Fit_Stats.columns:
            for (igamma,ip,ifit,itflow),fitdata in self.RatInt_Fit_Stats['boot'].items():
                outDict['Fits'][igamma][ip][ifit][itflow] = fitdata.GetOutputDict()

        for col_key,col_data in self.Rat_Stats.items():
            if 'Avg' in col_key or 'Std' in col_key: continue
            for (igamma,ip,it,itflow,its),idata in col_data.items():
                outDict[col_key][igamma][ip][it][itflow][its] = AvgStdToFormat(idata.Avg,idata.Std)


        for col_key,col_data in self.RatInt_Stats.items():
            if 'Avg' in col_key or 'Std' in col_key: continue
            for (igamma,ip,it,itflow),idata in col_data.items():
                outDict[col_key+'_Int'][igamma][ip][it][itflow] = AvgStdToFormat(idata.Avg,idata.Std)

        y = getattr(self, "SqrtFact", None)
        if y is not None and hasattr(self.SqrtFact,'iteritems'):
            for (ip,it),idata in self.SqrtFact.items():
                outDict['SqrtFact'][ip][it] = AvgStdToFormat(idata.Avg,idata.Std)


        if any(np.array(self.C3FOclass.NJN.ncfg_list) > 0) and self.show_cfgs:
            for ((istream,iccfg),icfg),ixsrc_list in zip(   iter(self.C3FOclass.NJNQFull_cfgs['configs'].items()),
                                                            self.C3FOclass.NJNQFull_cfgs['xsrc_list'].values):
                outDict['cfglist'][istream][icfg] = 'n_xsrc_'+str(len(ixsrc_list))


        ## Write human readable file with data
        WriteXml(self.HumanFile,{'Results':outDict})
        if any(np.array(self.C3FOclass.NJN.ncfg_list) > 0):
            for ((istream,iccfg),icfg),ixsrc_list in zip(   iter(self.C3FOclass.NJNQFull_cfgs['configs'].items()),
                                                            self.C3FOclass.NJNQFull_cfgs['xsrc_list'].values):
                outDict['cfglist'][istream][icfg] = 'n_xsrc_'+str(len(ixsrc_list))
        ## TODO, include fitting in the excell file
        WriteExcel( self.ExcelFile,{'Rat_results':deepcopy(self.Rat_Stats),'RatInt_results':deepcopy(self.RatInt_Stats)},
                    cfglist=self.C3FOclass.NJNQFull_cfgs,params=excel_params)

        ## pickles rest of data for reading
        self.RemoveFuns()
        WritePickle(self.PickleFile,self.__dict__)
        self.GetFuns()

    def GetFuns(self):
        if self.Rat_Fit_Stats is not None and 'boot' in self.Rat_Fit_Stats.columns:
            for ival in self.Rat_Fit_Stats['boot'].values:
                ival.GetFuns()
        if hasattr(self.C2class,'GetFuns'):
            self.C2class.GetFuns()
        if hasattr(self.C2Set,'GetFuns'):
            self.C2Set.GetFuns()
        self.C3FOclass.GetFuns()

    def RemoveFuns(self):
        if self.Rat_Fit_Stats is not None and 'boot' in self.Rat_Fit_Stats.columns:
            for ival in self.Rat_Fit_Stats['boot'].values:
                ival.RemoveFuns()
        if hasattr(self.C2class,'GetFuns'):
            self.C2class.RemoveFuns()
        if hasattr(self.C2Set,'GetFuns'):
            self.C2Set.RemoveFuns()
        self.C3FOclass.RemoveFuns()


    ##################################


    def ReadAndWrite(self,WipeData=True,DefWipe=False,CheckCfgs=False):
        self.Read(WipeData=WipeData,DefWipe=DefWipe)
        self.CreateC2Sink()
        self.Bootstrap(WipeData=WipeData,DefWipe=DefWipe,CheckCfgs=CheckCfgs)
        self.CalcRatio()
        self.Fit()
        self.Write()

    def GetTflowList(self,load_dict):
        return load_dict['Rat_Stats'].index.levels[2]

    def FixRead(self,readdict):
        rewrite = False
        if not isinstance(readdict['C2class'],str):
            if not isinstance(readdict['C2Set'],str):
                readdict['C2Set'] = 'Complete'
                rewrite = True
        return readdict,rewrite

    def LoadPickle(self,DefWipe=False,WipeData = True,CheckCfgs=False,CheckMom=True):
        # print
        # print 'NJNQ'
        # print self.HumanFile
        # print self.PickleFile
        if os.path.isfile(self.PickleFile) and not DefWipe:
            print('Loading Pickle file ' , self.PickleFile)
            loadeddict = ReadPickleWrap(self.PickleFile)
            checklist = self.thischecklist
            loadeddict,rewrite = self.FixRead(loadeddict)
            if CheckCfgs: checklist = checklist + ['cfglist']
            if CheckMom: checklist = checklist + ['checkmomlist']
            if 'tflow_list' not in  list(loadeddict.keys()):
                loadeddict['tflow_list'] = self.GetTflowList(loadeddict)
            if CheckClass(self.__dict__,loadeddict,checklist):
                self.__dict__.update(loadeddict)
                if rewrite: self.Write()
                if Wipe_All_Fits:
                    self.Rat_Fit_Stats = pa.DataFrame()
                elif len(self.Rat_Fit_Stats.values)> 0 and not isinstance(self.Rat_Fit_Stats.values[0][0],sff.SetOfFitFuns):
                    print('Legacy file, fits need to be recomputed')
                    self.Rat_Fit_Stats = pa.DataFrame()
            else:
                if os.path.isfile(self.PickleFile+'.bak'):
                    print('using backupfile for '+self.name)
                    self.PickleFile = self.PickleFile+'.bak'
                    self.HumanFile = self.HumanFile+'.bak'
                    self.ExcelFile = self.ExcelFile.replace('.xlsx','.bak.xlsx')
                    self.LoadPickle(DefWipe=DefWipe,WipeData=WipeData,CheckCfgs=CheckCfgs,CheckMom=CheckMom)
                    return
                print('Warning, file ' , self.HumanFile , ' has different parameters than this instance, reading and writing over file now')
                print()
                self.ReadAndWrite(WipeData=WipeData,DefWipe=DefWipe)
            # self.Write()
        else:
            if not os.path.isfile(self.PickleFile):
                print()
                print('pickle file not found: \n', self.PickleFile)
            self.ReadAndWrite(WipeData=WipeData,DefWipe=DefWipe)

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
        if self.thisshiftscale == 'Not Set':
            xlen = np.abs(xlims[1]-xlims[0])
            self.thisshiftscale = self.thisshift*xlen
            return self.thisshiftscale
        else:
            return self.thisshiftscale


    def Importtflow(self,tflowlist):
        if len(tflowlist) == 0:
            self.tflowfit = self.C3FOclass.tflowfit
        else:
            self.tflowfit = tflowlist


    ## make sure all xlims are the same, because this sets the window of xlim (needed for calculating shifts)
    ## gammamomlist is a dictionary where the keys are the gamma matricies, and the values are lists of momenta
    def Plot(   self,plot_class,xlims='tsink',gammamomlist={},tflowlist = [],tsum_list=[],
                thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine'):
        if xlims == 'tsink': xlims = [0,int(self.C3FOclass.NJN.tsink.replace('tsink',''))+1]
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        self.thisshiftscale = 'Not Set'
        thisshift = self.GetShift(xlims,thisshift)
        if 'boot' not in self.Rat_Stats.columns:
            raise EnvironmentError('no data to plot, maybe try to read in and bootstrap before plotting')
        if len(list(gammamomlist.keys())) == 0:
            gammamomlist = self.defgammamomlist
        if len(tflowlist) == 0:
            itflow = self.C3FOclass.tflowlist
        if len(tsum_list) == 0:
            tsum_list = ['ts'+str(its) for its in self.C3FOclass.boot_tsum_list]
        # for (igamma,ip,itflow),Rat_data in self.Rat_Stats['boot'].groupby(level=('current_operator','momentum','flow_times')):
        #     if itflow not in tflowlist: continue
        #     if igamma not in gammamomlist.keys(): continue
        #     if ip not in gammamomlist[igamma]: continue
        its = tsum_list[0]
        igamma = list(gammamomlist.keys())[0]
        ip = self.latparams.TOpform(list(gammamomlist.values())[0][0],pref='q')
        ip = ip.replace('qq','q')

        itflow = tflowlist[0]
        if (igamma,ip,itflow,its) not in self.Rat_Stats['boot'].reorder_levels([0,1,2,4,3]).index:
            print('Warning: data not found to plot', igamma,ip,itflow,its)
            return plot_class
        databoot = self.Rat_Stats.loc[:,'boot']
        dataavg = databoot.apply(lambda x : x.Avg)
        dataerr = databoot.apply(lambda x : x.Std)
        this_key = (igamma,ip,itflow,slice(None),its)
        tlist = list(self.Rat_Stats.loc[this_key,'boot'].keys())
        if len(tlist) > 0 and isinstance(tlist[0],(list,tuple,np.ndarray)):
            tlist = [it[3] for it in tlist]

        # dataavg,dataerr = [],[]
        # for idata in databoot[xlims[0]:xlims[1]]:
        #     dataavg.append(idata.Avg)
        #     dataerr.append(idata.Std)
        tlist = tlist[xlims[0]:xlims[1]]
        # print 'plotting ' , igamma,ip
        # pleg,gammaleg = '',''
        # if len(gammamomlist[igamma]) > 1:
        #     pleg = ' $'+ip+'$'
        # if len(gammamomlist.keys()) > 1:
        #     gammaleg = ' $'+igamma+'$'
        # pleg = ' $'+ip+'$'
        # if '_cmplx' in igamma:
        #     gammaleg = ' $i'+igamma.replace('_cmplx','')+'$'
        # else:
        #     gammaleg = ' $'+igamma+'$'
        # tflowleg = ' $'+tflowTOLeg(itflow)+'$'
        # tsumleg = ' $'+tsumTOLeg(its)+'$'
        tlistplot,self.centert = CenterT(np.array(list(map(untstr,tlist)))-twopttshift+thisshift)
        hold_series = pa.Series()
        # hold_series['x_data'] = tlistplot
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = dataavg
        hold_series['yerr_data'] = dataerr
        hold_series['xerr_data'] = None
        hold_series['key_select'] = this_key
        hold_series['type'] = 'error_bar_vary'
        hold_series['fit_class'] = None
        hold_series['label'] = self.LegLab#+pleg+gammaleg+tflowleg+tsumleg
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        hold_series['xaxis'] = None
        hold_series['otherXvals'] = None
        hold_series['xdatarange'] = None
        hold_series['Phys'] = None
        hold_series['ShowPar'] = None
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class


    ## make sure all xlims are the same, because this sets the window of xlim (needed for calculating shifts)
    ## gammamomlist is a dictionary where the keys are the gamma matricies, and the values are lists of momenta
    def PlotTsum(   self,plot_class,xlims='All',gammamomlist={},tflowlist = [],tcurr=[],
                    thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine'):
        if xlims == 'All':
            if 'sym' in self.C3FOclass.tr_sum_type:
                xlims = [0,self.nt//2]
            else:
                xlims = [0,self.nt]
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        self.thisshiftscale = 'Not Set'
        thisshift = self.GetShift(xlims,thisshift)
        if 'boot' not in self.Rat_Stats.columns:
            raise EnvironmentError('no data to plot, maybe try to read in and bootstrap before plotting')
        if len(list(gammamomlist.keys())) == 0:
            gammamomlist = self.defgammamomlist
        if len(tflowlist) == 0:
            tflowlist = self.C3FOclass.tflowlist
        if len(tcurr) == 0:
            tcurr = ['t'+str(int(str(self.C3FOclass.NJN.tsink).replace('tsink',''))//2)]
        # for (igamma,ip,itflow),Rat_data in self.Rat_Stats['boot'].groupby(level=('current_operator','momentum','flow_times')):
        #     if itflow not in tflowlist: continue
        #     if igamma not in gammamomlist.keys(): continue
        #     if ip not in gammamomlist[igamma]: continue
        it = tcurr[0]
        igamma = list(gammamomlist.keys())[0]
        ip = self.latparams.TOpform(list(gammamomlist.values())[0][0],pref='q')
        ip = ip.replace('qq','q')
        itflow = tflowlist[0]
        if (igamma,ip,itflow,it) not in self.Rat_Stats['boot']:
            print('Warning: data not found to plot', igamma,ip,itflow,it)
            return plot_class
        this_key = (igamma,ip,itflow,it,slice(None))
        databoot = self.Rat_Stats.loc[:,'boot']
        dataavg = databoot.apply(lambda x : x.Avg)
        dataerr = databoot.apply(lambda x : x.Std)
        # tslist,databoot = zip(*list(self.Rat_Stats['boot'][igamma,ip,itflow,it].iteritems()))
        tslist = list(self.Rat_Stats.loc[this_key,'boot'].keys())
        # dataavg,dataerr = [],[]
        # for idata in databoot[xlims[0]:xlims[1]]:
        #     dataavg.append(idata.Avg)
        #     dataerr.append(idata.Std)
        if len(tslist) > 0 and isinstance(tslist[0],(tuple,list,np.ndarray)):
            tslist = [this_t[-1] for this_t in tslist]
        tslist = np.array(list(map(untstr,tslist[xlims[0]:xlims[1]])))
        # print 'plotting ' , igamma,ip
        # pleg,gammaleg = '',''
        # if len(gammamomlist[igamma]) > 1:
        #     pleg = ' $'+ip+'$'
        # if len(gammamomlist.keys()) > 1:
        #     gammaleg = ' $'+igamma+'$'
        # pleg = ' $'+ip+'$'
        # if '_cmplx' in igamma:
        #     gammaleg = ' $i'+igamma.replace('_cmplx','')+'$'
        # else:
        #     gammaleg = ' $'+igamma+'$'
        # tflowleg = ' $'+tflowTOLeg(itflow)+'$'
        # tleg = ' $'+tsinkTOLeg(it)+'$'
        hold_series = pa.Series()
        hold_series['x_data'] = 'from_keys'
        # hold_series['x_data'] = tslist
        hold_series['y_data'] = dataavg
        hold_series['yerr_data'] = dataerr
        hold_series['xerr_data'] = None
        hold_series['key_select'] = this_key
        hold_series['type'] = 'error_bar_vary'
        hold_series['fit_class'] = None
        hold_series['label'] = self.LegLab#+pleg+gammaleg+tflowleg+tleg
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        hold_series['xdatarange'] = None
        hold_series['xaxis'] = None
        hold_series['otherXvals'] = None
        hold_series['Phys'] = None
        hold_series['ShowPar'] = None
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class


    def GetBestFitr(self,igamma,ip,itflow):
        if 'boot' not in self.Rat_Fit_Stats.columns or (igamma,ip,itflow) not in self.Rat_Fit_Stats.index:
            errstr = ' '.join(igamma,ip,str(itflow))+' not foun in ratio fit list, cannot find best'
            raise EnvironmentError(errstr)
        chi_str,cut_fit_result = self.Rat_Fit_Stats.at[(igamma,ip,itflow),'boot'].GetBestChi()
        return r'$\chi^{2}_{pdf}={:.3f}$'.format(chi_str),cut_fit_result

    ## state is 1 or 2 state fit
    ## fitr is fit range e.g. fitr5-10
    ## fit range of 'Data' is just the original fit range
    def FitPlot(    self,plot_class,fitr,xlims = 'Data',tsum_list=[],gammamomlist={},
                    tflowlist = [],thiscol='PreDefine',thisshift='PreDefine'):
        self.CheckCol(thiscol)
        thisshift = self.GetShift(xlims,thisshift)
        self.Importtflow(tflowlist)
        if self.centert == 'Not Set':
            print('print centert not set, defaulting to 0')
            self.centert = 0.0
            # raise IOError('Run Plot() before running FitPlot() to find the center value')
        # if 'boot' not in self.Rat_Fit_Stats.columns:
        #     raise EnvironmentError('No fitting has been done, maybe try to Fit before plotting fit')
        if len(list(gammamomlist.keys())) == 0:
            gammamomlist = self.defgammamomlist
        if len(tflowlist) == 0:
            tflowlist = self.C3FOclass.tflowlist
        if len(tsum_list) == 0:
            tsum_list = [its for its in self.C3FOclass.boot_tsum_list]
        # if (igamma,ip,fitr) not in self.Rat_Fit_Stats.index:
        # for igamma,momlist in gammamomlist.iteritems():
        #     for ip in momlist:
        #         for itflow in tflowlist:
        igamma,momlist = next(iter(gammamomlist.items()))
        ip = self.latparams.TOpform(list(gammamomlist.values())[0][0],pref='q')
        ip = ip.replace('qq','q')
        itflow = tflowlist[0]
        itsum = tsum_list[0]
        if isinstance(fitr,str):
            fitr = fitr.split('_')
        if isinstance(itsum,str):
            itsum = untstr(itsum)
        if fitr == 'Best':
            chi_legend,fitr = self.GetBestFitr(igamma,ip,itflow)

        if 'boot' not in self.Rat_Fit_Stats.columns:
            print('no fits done, attempting fit now')
            self.Fit(fit_range=fitr[0],tsum_ranges=fitr[1])
        fit_data = sff.PullSOFSeries(self.Rat_Fit_Stats['boot'])
        fit_data = Series_fix_key(fit_data,5,'fittwor','tsumfitr')
        this_fun = fit_data.index[0][3]
        this_key = (igamma,ip,itflow,this_fun,fitr[0],fitr[1])
        if this_key not in fit_data:
            print(this_key , 'not present in Rat Fit Stats, fitting now')
            self.Fit(fit_range=fitr)
            fit_data = sff.PullSOFSeries(self.Rat_Fit_Stats['boot'])
            fit_data = Series_fix_key(fit_data,5,'fittwor','tsumfitr')
            this_fun = fit_data.index[0][3]
            this_key = (igamma,ip,itflow,this_fun,fitr[0],fitr[1])

        hold_series = null_series
        hold_series['key_select'] = this_key
        hold_series['type'] = 'fit_vary'
        hold_series['fit_class'] = fit_data
        # hold_series['label'] = fit_data[this_key].name +' ' + chi_legend
        hold_series['label'] = 'ACTUAL_VALUE_fit'
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        hold_series['xaxis'] = 0
        hold_series['xdatarange'] = xlims
        hold_series['otherXvals'] = [itsum]
        hold_series['ShowPar'] = r'Par1'
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class

                                    ## state is 1 or 2 state fit
    ## fitr is fit range e.g. fitr5-10
    ## fit range of 'Data' is just the original fit range
    def FitSumPlot( self,plot_class,fitr,xlims = 'Data',tcurr=[],gammamomlist={},
                    tflowlist = [],thiscol='PreDefine',thisshift='PreDefine'):
        self.CheckCol(thiscol)
        thisshift = self.GetShift(xlims,thisshift)
        self.Importtflow(tflowlist)
        # if 'boot' not in self.Rat_Fit_Stats.columns:
        #     raise EnvironmentError('No fitting has been done, maybe try to Fit before plotting fit')
        if len(list(gammamomlist.keys())) == 0:
            gammamomlist = self.defgammamomlist
        if len(tflowlist) == 0:
            tflowlist = self.C3FOclass.tflowlist
        if len(tcurr) == 0:
            tcurr = [int(str(self.C3FOclass.NJN.tsink).replace('tsink',''))//2]
        # if (igamma,ip,fitr) not in self.Rat_Fit_Stats.index:
        # for igamma,momlist in gammamomlist.iteritems():
        #     for ip in momlist:
        #         for itflow in tflowlist:
        igamma,momlist = next(iter(gammamomlist.items()))
        ip = self.latparams.TOpform(list(gammamomlist.values())[0][0],pref='q')
        ip = ip.replace('qq','q')
        itflow = tflowlist[0]
        tcurr = tcurr[0]
        if isinstance(fitr,str):
            fitr = fitr.split('_')
        if isinstance(tcurr,str):
            tcurr = untstr(tcurr)

        if 'boot' not in self.Rat_Fit_Stats.columns:
            print('no fits done, attempting fit now')
            self.Fit(fit_range=fitr[0],tsum_ranges=fitr[1])
        fit_data = sff.PullSOFSeries(self.Rat_Fit_Stats['boot'])
        fit_data = Series_fix_key(fit_data,5,'fittwor','tsumfitr')
        this_fun = fit_data.index[0][3]
        this_key = (igamma,ip,itflow,this_fun,fitr[0],fitr[1])
        if this_key not in fit_data:
            print(this_key , 'not present in Rat Fit Stats, fitting now')
            self.Fit(fit_range=fitr)
            fit_data = sff.PullSOFSeries(self.Rat_Fit_Stats['boot'])
            fit_data = Series_fix_key(fit_data,5,'fittwor','tsumfitr')
            this_fun = fit_data.index[0][3]
            this_key = (igamma,ip,itflow,this_fun,fitr[0],fitr[1])

        hold_series = null_series
        hold_series['type'] = 'fit_vary'
        hold_series['fit_class'] = fit_data
        hold_series['key_select'] = this_key
        # hold_series['label'] = fit_class[this_key].name
        hold_series['label'] = 'ACTUAL_VALUE_fit'
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        hold_series['xdatarange'] = xlims
        hold_series['xaxis'] = 1
        hold_series['otherXvals'] = [tcurr]
        hold_series['ShowPar'] = r'Par1'
        hold_series['fmt_class'] = KeyForamtting(self.latparams)
        plot_class.AppendData(hold_series)
        return plot_class




    ## ##################### ###############################################


    ## Comparisons

    def __str__(self):
        return '_'.join([self.SubDir,self.name])


    def __iter__(self):
        return self

    def iteritems(self):
        return iter(self.Rat_Stats['boot'].items())

    def items(self):
        return list(self.Rat_Stats['boot'].items())

    def itemsAvgStd(self):
        outlist = []
        for ivals in self.Rat_Stats.itertuples():
            outlist.append((ivals[0],ivals[2],ivals[3]))
        return outlist

    def values(self):
        return self.Rat_Stats['boot'].values

    def itervalues(self):
        return iter(self.Rat_Stats['boot'].values)

    def valuesAvgStd(self):
        outlist = []
        for ivals in self.Rat_Stats.itertuples():
            outlist.append((ivals[2],ivals[3]))
        return outlist

    def keys(self):
        return self.Rat_Stats.index


    def __setitem__(self, key, value ):
        self.Rat_Stats['boot'][key] = value

    def __getitem__(self, key):
        return self.Rat_Stats['boot'][key]

    def __reversed__(self):
        self.Rat_Stats = self.Rat_Stats.iiloc[::-1]
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


    def ReduceIndicies(self):
        indexl = []
        for (igamma,ip,it,itflow,its) in list(self.Rat_Stats['boot'].keys()):
            first_gamma = MultiSplit(igamma,op_str_list+['FO','FOFull'])[0]
            this_gamma = first_gamma.replace('FO','')
            this_gamma = this_gamma.replace('FOFull','')
            this_gamma = this_gamma.replace('__','_')
            if this_gamma[-1] == '_':
                this_gamma = this_gamma[:-1]
            this_p = MultiSplit(ip,op_str_list+['FO','FOFull','_'])[0]
            indexl.append((this_gamma,this_p,it,itflow,its))
        indicies = pa.MultiIndex.from_tuples(indexl,names=self.Rat_Col_Names)
        this_series = pa.Series(self.Rat_Stats['boot'].values,index=indicies)
        self.Rat_Stats = pa.DataFrame(columns=['boot'],index=indicies)
        self.Rat_Stats.loc[:,'boot'] = this_series
        self.Stats()


    def sin(self):
        result = RatioFOFullCorr(thisnboot=self.nboot, cfglist=self.cfglist[['configs','xsrc_list']],
                         Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift,
                         name = 'sin('+self.name+')',filename='sin_'+self.filename)
        result.Rat_Stats = pa.DataFrame(columns = ['boot'],index=self.Rat_Stats.index)
        result.Rat_Stats.loc[:,'boot'] = np.sin(self.Rat_Stats.loc[:,'boot'])
        return result

    def cos(self):
        result = RatioFOFullCorr(thisnboot=self.nboot, cfglist=self.cfglist[['configs','xsrc_list']],
                         Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift,
                         name = 'cos('+self.name+')',filename='cos_'+self.filename)
        result.Rat_Stats = pa.DataFrame(columns = ['boot'],index=self.Rat_Stats.index)
        result.Rat_Stats.loc[:,'boot'] = np.cos(self.Rat_Stats.loc[:,'boot'])
        return result


    def Overload_wrap(self,Rat2,this_fun):
        ## TODO, include some name and filename formatting
        thisname,thisfilename = self.name,self.filename
        if hasattr(Rat2,'name'):
            thisname += '_'+this_fun.__name__+'_'+Rat2.name
        if hasattr(Rat2,'filename'):
            thisfilename += '_'+this_fun.__name__+'_'+Rat2.filename
        result = RatioFOFullCorr(   thisnboot=self.nboot, cfglist=self.cfglist[['configs','xsrc_list']],
                                    Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.Rat_Stats = pa.DataFrame(columns = ['boot'])
        # result.UpdateName('+',self,Rat2)
        if 'boot' in self.Rat_Stats:
            if isinstance(Rat2,RatioFOFullCorr):
                # result.Rat_Stats.loc[:,'boot'] = this_fun(self.Rat_Stats.loc[:,'boot'],Rat2.Rat_Stats.loc[:,'boot'])
                if 'boot' in Rat2.Rat_Stats:
                    vall,ilist = [],[]
                    for (igamma,ip,itflow,it,its),irat in self.Rat_Stats['boot'].items():
                        if (its,it,itflow) not in Rat2.Rat_Stats['boot'].index.swaplevel(0,4).swaplevel(1,3):
                            out_str = ' missmatch in euclidian times / flow time / op times \
                                        when combining ratio functions \n'
                            out_str += str(self) + '\n'
                            out_str += this_fun.__name__+'\n'
                            out_str += str(Rat2)
                            raise EnvironmentError(out_str)
                        for (igamma2,ip2,itflowdump,itdump,itsdump),irat2 in Rat2.Rat_Stats.loc[(slice(None),slice(None),itflow,it,its),'boot'].items():
                            igammacomb,ipcomb = igamma,ip
                            if igamma2 not in igamma:
                                igammacomb += '_'+this_fun.__name__+'_'+igamma2
                            if ip2 not in ip:
                                ipcomb += '_'+this_fun.__name__+'_'+ip2
                            ilist.append((igammacomb,ipcomb,itflow,it,its))
                            vall.append(this_fun(irat,irat2))
                    if len(ilist) > 0:
                        indicies = pa.MultiIndex.from_tuples(ilist,names=self.Rat_Col_Names)
                        result.Rat_Stats.loc[:,'boot'] = pa.Series(vall,index=indicies)
                else:
                    raise EnvironmentError('Run Bootstrap on second class before combining ratio functions')
            elif isinstance(Rat2,RatioFOCorr):
                # result.Rat_Stats.loc[:,'boot'] = this_fun(self.Rat_Stats.loc[:,'boot'],Rat2.Rat_Stats.loc[:,'boot'])
                if 'boot' in Rat2.Rat_Stats:
                    vall,ilist = [],[]
                    for (igamma,ip,itflow,it,its),irat in self.Rat_Stats['boot'].items():
                        if (it,itflow) not in Rat2.Rat_Stats['boot'].index.swaplevel(0,3).swaplevel(1,2):
                            out_str = 'missmatch in euclidian times / flow times when combining ratio functions \n'
                            out_str += str(self) + '\n'
                            out_str += this_fun.__name__+'\n'
                            out_str += str(Rat2)
                            raise EnvironmentError(out_str)
                        for (igamma2,ip2,itflowdump,itdump),irat2 in Rat2.Rat_Stats.loc[(slice(None),slice(None),itflow,it),'boot'].items():
                            igammacomb,ipcomb = igamma+'_FOFull_',ip
                            igammacomb += '_'+this_fun.__name__+'_'+igamma2+'_FO'
                            if ip2 not in ip:
                                ipcomb += '_'+this_fun.__name__+'_'+ip2
                            ilist.append((igammacomb,ipcomb,itflow,it,its))
                            vall.append(this_fun(irat,irat2))
                    if len(ilist) > 0:
                        indicies = pa.MultiIndex.from_tuples(ilist,names=self.Rat_Col_Names)
                        result.Rat_Stats.loc[:,'boot'] = pa.Series(vall,index=indicies)
                else:
                    raise EnvironmentError('Run Bootstrap on second class before combining ratio functions')
            elif isinstance(Rat2,RatioCorr):
                # result.Rat_Stats.loc[:,'boot'] = this_fun(self.Rat_Stats.loc[:,'boot'],Rat2.Rat_Stats.loc[:,'boot'])
                if 'boot' in Rat2.Rat_Stats:
                    vall,ilist = [],[]
                    for (igamma,ip,itflow,it,its),irat in self.Rat_Stats['boot'].items():
                        if it not in Rat2.Rat_Stats['boot'].index.swaplevel(0,2):
                            out_str = 'missmatch in euclidian times when combining ratio functions \n'
                            out_str += str(self) + '\n'
                            out_str += this_fun.__name__+'\n'
                            out_str += str(Rat2)+'\n'
                            out_str += it
                            raise EnvironmentError(out_str)
                        for (igamma2,ip2,itdump),irat2 in Rat2.Rat_Stats.loc[(slice(None),slice(None),it),'boot'].items():
                            igammacomb,ipcomb = igamma+'_FOFull_',ip
                            # if igamma != igamma2:
                            igammacomb += this_fun.__name__+'_'+igamma2
                            if ip2 not in ip:
                                ipcomb += '_'+this_fun.__name__+'_'+ip2
                            ilist.append((igammacomb,ipcomb,itflow,it,its))
                            vall.append(this_fun(irat,irat2))
                    if len(ilist) > 0:
                        indicies = pa.MultiIndex.from_tuples(ilist,names=self.Rat_Col_Names)
                        result.Rat_Stats.loc[:,'boot'] = pa.Series(vall,index=indicies)
                else:
                    raise EnvironmentError('Run Bootstrap on second class before combining ratio functions')
            else:
                try:
                    result.Rat_Stats.loc[:,'boot'] = this_fun(self.Rat_Stats.loc[:,'boot'], Rat2)
                except Exception as err:
                    print(type(Rat2), Rat2)
                    print(self.Rat_Stats)
                    raise EnvironmentError(str(err)+'\nInvalid value to combine with RatioFOCorr class')
            return result
        else:
            raise EnvironmentError('Run Bootstrap on first class before combining ratio functions')

    def __add__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOFullCorr.parentlist]):
            return Rat2.__radd__(self)
        else:
            return self.Overload_wrap(Rat2,op.add)

    def __sub__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOFullCorr.parentlist]):
            return Rat2.__rsub__(self)
        else:
            return self.Overload_wrap(Rat2,op.sub)

    def __mul__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOFullCorr.parentlist]):
            return Rat2.__rmul__(self)
        else:
            return self.Overload_wrap(Rat2,op.mul)

    def __div__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOFullCorr.parentlist]):
            return Rat2.__rdiv__(self)
        else:
            return self.Overload_wrap(Rat2,op.truediv)

    def __pow__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOFullCorr.parentlist]):
            return Rat2.__rpow__(self)
        else:
            return self.Overload_wrap(Rat2,op.pow)


    ## Right multiplication functions


    def __radd__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOFullCorr.parentlist]):
            return Rat2.__add__(self)
        else:
            return self.Overload_wrap(Rat2,flip_fun_2arg(op.add))

    def __rsub__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOFullCorr.parentlist]):
            return Rat2.__sub__(self)
        else:
            return self.Overload_wrap(Rat2,flip_fun_2arg(op.sub))

    def __rmul__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOFullCorr.parentlist]):
            return Rat2.__mul__(self)
        else:
            return self.Overload_wrap(Rat2,flip_fun_2arg(op.mul))

    def __rdiv__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOFullCorr.parentlist]):
            return Rat2.__div__(self)
        else:
            return self.Overload_wrap(Rat2,flip_fun_2arg(op.truediv))


    def __rpow__(self,Rat2):
        if any([ipar in str(type(Rat2)) for ipar in RatioFOFullCorr.parentlist]):
            return Rat2.__pow__(self)
        else:
            return self.Overload_wrap(Rat2,flip_fun_2arg(op.pow))




    ## Operator overloading Starting to get tedious, leave to last or think of better way of implementing (maybe use predefined funtion passed in?)

def TestRat(DefWipe=False):
    from Params import defInfo
    defInfo['gammas'] = ['g4']
    defInfo['qmom'] = ['0 0 0']
    defInfo['Projector'] = 'P4'
    this_info = pa.Series()
    this_info['save_file'] = this_dir+'/TestGraphs/testRat.pdf'
    this_info['title'] = 'Test Ratio'
    data = RatioCorr(Info=defInfo)
    data.LoadPickle(DefWipe=DefWipe)
    import PlotData as jpl
    data_plot = jpl.Plotting(plot_info=this_info)
    data_plot = data.Plot(data_plot,thiscol='Blue',thissym='o',thisshift=0.0)
    data_plot = data.FitPlot(data_plot,defInfo['FitsRF'][0])
    data_plot.PrintData()
    data_plot.PlotAll()
    return data


def TestRat2(DefWipe=False):
    from Params import defInfo
    defInfo['gammas'] = ['g3g5_cmplx']
    defInfo['qmom'] = ['0 0 0']
    defInfo['Projector'] = 'P3'
    this_info = pa.Series()
    this_info['save_file'] = this_dir+'/TestGraphs/testRat2.pdf'
    this_info['title'] = 'Test Ratio2'
    data = RatioCorr(Info=defInfo)
    data.LoadPickle(DefWipe=DefWipe)
    import PlotData as jpl
    data_plot = jpl.Plotting(plot_info=this_info)
    data_plot = data.Plot(data_plot,thiscol='Blue',thissym='o',thisshift=0.0)
    data_plot = data.FitPlot(data_plot,defInfo['FitsRF'][0])
    data_plot.PrintData()
    data_plot.PlotAll()
    return data

def TestRatFO(DefWipe=False):
    from Params import defInfoFlow
    defInfoFlow['gammas'] = ['g4']
    defInfoFlow['qmom'] = ['0 0 1']
    defInfoFlow['Projector'] = 'P3'
    this_info = pa.Series()
    this_info['save_file'] = this_dir+'/TestGraphs/testRatFO.pdf'
    this_info['title'] = 'Test Ratio FO'
    dataflow = RatioFOCorr(Info=defInfoFlow)
    dataflow.LoadPickle(DefWipe=DefWipe)
    import PlotData as jpl
    data_plot = jpl.Plotting(plot_info=this_info)
    data_plot = dataflow.Plot(data_plot,thiscol='Blue',thissym='o',thisshift=0.0,tflowlist=[defInfoFlow['tflowfit'][0]])
    data_plot = dataflow.FitPlot(data_plot,defInfoFlow['FitsRF'][0],tflowlist=[defInfoFlow['tflowfit'][0]])
    data_plot.PrintData()
    data_plot.PlotAll()
    return dataflow


def TestRatFO2(DefWipe=False):
    from Params import defInfoFlow
    defInfoFlow['gammas'] = ['g4']
    defInfoFlow['qmom'] = ['0 1 -1']
    defInfoFlow['Projector'] = 'P3'
    defInfoFlow['tflowlist'] = np.array(['t_f4.01','t_f5.01','t_f6.01'])
    defInfoFlow['tflowfit'] = defInfoFlow['tflowlist']
    this_info = pa.Series()
    this_info['save_file'] = this_dir+'/TestGraphs/testRatFO2.pdf'
    this_info['title'] = 'Test Ratio FO2'
    dataflow = RatioFOCorr(Info=defInfoFlow)
    dataflow.LoadPickle(DefWipe=DefWipe)
    import PlotData as jpl
    data_plot = jpl.Plotting(plot_info=this_info)
    data_plot = dataflow.Plot(data_plot,thiscol='Blue',thissym='o',thisshift=0.0,tflowlist=[defInfoFlow['tflowfit'][0]])
    data_plot = dataflow.FitPlot(data_plot,defInfoFlow['FitsRF'][0],tflowlist=[defInfoFlow['tflowfit'][0]])
    data_plot.PrintData()
    data_plot.PlotAll()
    return dataflow

def TestRatFOFull(DefWipe=False):
    from Params import defInfoFlow
    this_info = pa.Series()
    this_info['save_file'] = this_dir+'/TestGraphs/testRatFOFull.pdf'
    this_info['title'] = 'Test Ratio FOFull'
    defInfoFlow['DoubSing'] = 'Neutron'
    defInfoFlow['gammas'] = ['g4']
    defInfoFlow['qmom'] = ['0 0 1']
    defInfoFlow['Projector'] = 'P3'
    defInfoFlow['tflowlist'] = np.array(['t_f4.01','t_f5.01','t_f6.01'])
    defInfoFlow['tflowfit'] = defInfoFlow['tflowlist']
    defInfoFlow['tsum_fit_list'] = ['ts28','ts29','ts30','ts31','ts32']
    defInfoFlow['tsum_list'] = ['ts28','ts30','ts32']
    defInfoFlow['FitsRFFull'] = [('taufitr4-8','tsumfitr28-32')]
    dataflow = RatioFOFullCorr(Info=defInfoFlow)
    dataflow.LoadPickle(DefWipe=DefWipe)
    import PlotData as jpl
    data_plot = jpl.Plotting(plot_info=this_info)
    data_plot = dataflow.Plot(data_plot,thiscol='Blue',thissym='o',thisshift=0.0,tflowlist=[defInfoFlow['tflowfit'][0]],tsum_list=['ts32'])
    data_plot = dataflow.FitPlot(data_plot,defInfoFlow['FitsRFFull'][0],tflowlist=[defInfoFlow['tflowfit'][0]],tsum_list=['ts32'])
    data_plot.PrintData()
    data_plot.PlotAll()

    # this_info = pa.Series()
    # this_info['save_file'] = this_dir+'/TestGraphs/testRatFOTsum.pdf'
    # this_info['title'] = 'Test Ratio FOFull Tsum'
    # defInfoFlow['DoubSing'] = 'Proton'
    # defInfoFlow['gammas'] = ['g4']
    # defInfoFlow['qmom'] = ['0 0 1']
    # defInfoFlow['Projector'] = 'P3'
    # # defInfoFlow['FitsRFFull'] = [('taufitr3-7','fitr1-10')]
    # dataflow = RatioFOFullCorr(Info=defInfoFlow)
    # dataflow.LoadPickle(DefWipe=False)
    # import PlotData as jpl
    # data_plot = jpl.Plotting(plot_info=this_info)
    # data_plot = dataflow.PlotTsum(data_plot,thiscol='Blue',thissym='o',thisshift=0.0)
    # dataflow.Rat_Fit_Stats = pa.DataFrame()
    # data_plot = dataflow.FitSumPlot(data_plot,defInfoFlow['FitsRFFull'][0],tcurr=['t6'],thiscol='Blue',thisshift=0.0)
    # defInfoFlow['DoubSing'] = 'Neutron'
    # dataflow = RatioFOFullCorr(Info=defInfoFlow)
    # dataflow.LoadPickle(DefWipe=False)
    # data_plot = dataflow.PlotTsum(data_plot,thiscol='Red',thissym='*',thisshift=0.01)
    # dataflow.Rat_Fit_Stats = pa.DataFrame()
    # data_plot = dataflow.FitSumPlot(data_plot,defInfoFlow['FitsRFFull'][0],tcurr=['t6'],thiscol='Red',thisshift=0.01)
    # data_plot.PrintData()
    # data_plot.PlotAll()
    return dataflow


if __name__ == '__main__':
    # dataRatFOFull = TestRatFOFull()
    dataRatFO = TestRatFO()
    dataRat = TestRat()
    print('Ratio function test complete, see variables dataRat and dataRatFO')
