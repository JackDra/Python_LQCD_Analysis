#!/usr/bin/env python

# import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
# import matplotlib.pyplot as pl
# import PlotData as jpl

# from BootStrapping import BootStrap
import numpy as np
from copy import deepcopy,copy
from MiscFuns import ODNested,mkdir_p,CreateNameAndLL,op_dict
from MiscFuns import flip_fun_2arg,GetInfoList,GenSinkFun,ReduceCfgs
import MomParams as mp
from Params import defFitDict,defmom2list,outputdir,scratchdir,Debug_Mode,Wipe_All_Fits
from Params import defFitAlpha,defMcore,nchi_threshold,nchi_fit_threshold,myeps
from MultiWrap import DoMulticore
# from Params import plot_the_ff
# import multiprocessing
# import glob,sys
import operator as op
# from FFFuns import *
from SetsOfRatios import SetOfRat
from SetsOfTwoPtCorrs import SetOfNNQ,SetOfNNQFull
from TwoPtCorrelators import TwoPointCorr,NNQCorr,NNQFullCorr
from FormFactors import FormFact,FFRenorm,FFQ2FitFun,FF_Imp
from QuantityLists import ObsListAbr
from BootStrapping import BootStrap

# from PredefFitFuns import FormFactorO1,FormFactorO2,FormFactorO3,LinearFitFun,DPfitfun
# from PredefFitFuns import FormFactorO1,LinearFitFun,DPfitfun
from PredefFitFuns import FormFactorO1
import FitFunctions as ff
from collections import OrderedDict
from XmlFormatting import QsqrdToFormat,unxmlfitr,tflowTOLeg,AvgStdToFormat,KeyForamtting
import os
import dill as pickle
import pandas as pa
# from FileIO import WriteXml,WritePickle,ReadPickleWrap,BackupFile,WriteFuns,ReadFuns
from FileIO import WriteXml,WritePickle,ReadPickleWrap,BackupFile
from GammaMatricies import CorrSpinTrace
from DSDefines import DSfundict
from TimeStuff import Timer
import itertools

## title for ratio comparison
## $G_{3}$ Fit Range Comparison $m_{\pi}=411MeV$ Neutron $F_{3}/2m_{N}$
## ylabel for CP odd form factor F3
## $\frac{F_{3}}{2m_{N}}$

class FFEquation(object):
    """

    This class uses FormFact from FormFactors.py and RatioCorr from RatioCorrelators.py

    sets up system of equations to solve:

    Af = R , with indicies:
    A_1a * f_a + a_1b * f_b + ... = R_1
    A_2a * f_a + a_2b * f_b + ... = R_2
    ...


    A = coefficients of form factors [ A_1a , A_1b , ... ]
                                     [ A_2a , A_2b , ... ]
                                     [ ...               ]

    R = extracted ratio function results (from fits to correlators or ratio functions)
    R = [ R_2, R_1, ... ]

    f = form factors to solve for [f_a, f_b,...]

    default q2list is usually [0,1,2,3,4] , (depends on what you calculate, set in MomParams.py if you like )


    """

    parentlist = []

    def reInit(self):
        self = FFEquation(self.currtype,q2list=self.q2list,pp2list=self.pp2list,cfglist=self.cfglistin, Info=self.Info,
                          fit2pt=self.fit2ptin,fitAlpha=self.fitAlphain,name=self.namein,
                          thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)


    def __init__(self,  currtype,q2list=defmom2list,pp2list = [0],cfglist={},
                        cfglist2pt={},Info={},fit2pt='From Def',fitAlpha='From Def',
                        name = '',thissym='Not Set',thiscol='Not Set',thisshift=0.0):

        self.current = 0
        self.thisshift = thisshift
        self.thiscol = thiscol
        self.thissym = thissym
        self.cfglistin = cfglist
        self.cfglist = cfglist
        self.q2list = q2list
        self.pp2list = pp2list
        self.Info = copy(Info)
        self.fit2ptin = fit2pt
        self.fitAlphain = fitAlpha
        self.nchi_threshold = nchi_threshold
        self.nchi_fit_threshold = nchi_fit_threshold
        if isinstance(fitAlpha,(list,tuple,np.ndarray)):
            self.fitAlphain = '_'.join(self.fitAlphain)
        self.currtype = currtype
        self.namein = name
        self.thisshiftscale = 'Not Set'
        self.Renormed = False
        self.Inefm = False

        if 'Rat_alpha' in list(Info.keys()): self.Rat_alpha = Info['Rat_alpha']
        else: self.Rat_alpha = False ## default to no symmetrizing, (usually full matrix isn't computed nowerdays.


        if 'Alpha_Full' in list(Info.keys()): self.doAlphafull = Info['Alpha_Full']
        else: self.doAlphafull = False

        if 'Chi_Cutoff' in list(Info.keys()): self.Do_Chi_cutoff = Info['Chi_Cutoff']
        else:  self.Do_Chi_cutoff = False

        if self.Do_Chi_cutoff:
            if 'Chi2pdf_Min' in list(Info.keys()): self.Chi2pdf_Min = Info['Chi2pdf_Min']
            else:  self.Chi2pdf_Min = 0.3
            if 'Chi2pdf_Max' in list(Info.keys()): self.Chi2pdf_Max = Info['Chi2pdf_Max']
            else:  self.Chi2pdf_Max = 0.7
        else:
            self.Chi2pdf_Min = 0
            self.Chi2pdf_Max = 0

        if 'RF_fit_range_minimum' in list(Info.keys()): self.RF_fit_range_minimum = Info['RF_fit_range_minimum']
        else:  self.RF_fit_range_minimum = 3

        if 'RF_fit_max' in list(Info.keys()): self.RF_fit_max = Info['RF_fit_max']
        else:  self.RF_fit_max = 11

        if 'RF_fit_min' in list(Info.keys()): self.RF_fit_min = Info['RF_fit_min']
        else:  self.RF_fit_min = 2


        if 'Rat_tsum_range' in list(Info.keys()): self.Rat_tsum_range = Info['Rat_tsum_range']
        else:  self.Rat_tsum_range = 'tsumfitr2-20'


        if 'Important_FF' in list(Info.keys()): self.Important_FF = Info['Important_FF']
        else:  self.Important_FF = True


        if 'Reduce_System' in list(Info.keys()): self.Reduce_System = Info['Reduce_System']
        else:  self.Reduce_System = True

        if 'Rat_comb_alpha_fun' in list(Info.keys()): self.ralpha_mult_fun = Info['Rat_comb_alpha_fun']
        else: self.ralpha_mult_fun = op.truediv

        if isinstance(self.ralpha_mult_fun,str):
            self.ralpha_mult_fun = op_dict[self.ralpha_mult_fun]


        self.doflow = ('Top' in self.currtype) or ('Wein' in self.currtype)



        if self.doflow:
            if 'Top' in self.currtype:
                self.flowtype = 'TopCharge'
                self.flowlab = 'Top'
            if 'Wein' in self.currtype:
                self.flowtype = 'Weinberg'
                self.flowlab = 'Wein'
            if 'Full_Rat' in list(Info.keys()): self.doflowfull = Info['Full_Rat']
            else: self.doflowfull = False
        else:
            self.flowtype = ''
            self.flowlab = ''
            self.doflowfull = False

        self.cfg_fromfile = False
        if isinstance(cfglist,str):
            self.cfg_fromfile = cfglist
            if os.path.isfile(cfglist):
                print('found previous cfglist, skipping search')
                with open(cfglist,'rb') as f:
                    self.cfglist = pickle.load(f)
            else:
                print('previous config file not found at')
                print(cfglist)
                self.cfglist={}
        else:
            self.cfglist = cfglist


        self.cfg_fromfile_2pt = False
        if isinstance(cfglist2pt,str):
            self.cfg_fromfile_2pt = cfglist2pt
            if os.path.isfile(cfglist2pt):
                print('found previous cfglist2pt, skipping search')
                with open(cfglist2pt,'rb') as f:
                    self.cfglist2pt = pickle.load(f)
            else:
                print('previous config2pt file not found at ')
                print(cfglist2pt)
                self.cfglist2pt={}
        else:
            self.cfglist2pt = cfglist2pt





        ## Info can contain fit ranges to be calculated
        ## must bein form Info['FlowFitFun'] = (Function Object , number of parameters )
        if 'FFiGuess' in list(Info.keys()):
            self.iGuess = Info['FFiGuess']
        else:
            self.iGuess = 'None'

        ## Info can contain fit ranges to be calculated
        ## must bein form Info['FFFits'] = [ 'fit#-#' , ...]
        if 'FFFits' in list(Info.keys()): self.PredefFits = Info['FFFits']
        else: self.PredefFits = []

        if 'tflowfit' in list(Info.keys()): self.tflows = Info['tflowfit']
        else: self.tflows = ['none']

        ## makes alpha have the same config list as everything else
        if 'AlphaFree' in list(Info.keys()): self.alphafree = Info['AlphaFree']
        else: self.alphafree = False

        # ## makes mass have same config list as everything else
        # if 'MassFree' in Info.keys(): self.MassFree = Info['MassFree']
        # else: self.MassFree = False



        ## DoubSing tells whether to read doublet (doub) or singlet (sing) contribution
        if 'DoubSing' in list(Info.keys()):
            if Info['DoubSing'] not in list(DSfundict.keys()):
                print('DSList (from DSDefines.py):')
                print(list(DSfundict.keys()))
                raise IOError("Doublet,Singlet type (Info['DoubSing']) not recognised")
            self.DS = Info['DoubSing']
        else: raise IOError('Please set DoubSing in Info for ThreePtCorrelators ')

        ## Info can contain fit ranges to be calculated
        ## must bein form Info['FlowFitFun'] = (Function Object , number of parameters )
        ## Fit Functions are now hard coded into FormFactors.py
        ## modify FFFitFuns in FormFactors.py to change
        self.SetFunction(self.currtype,self.DS)

        self.ppparams = mp.LatticeParameters(mom2list=pp2list,Info=Info)
        self.ppparams.LoadPickle()
        self.qparams = mp.LatticeParameters(mom2list=q2list,Info=Info)
        self.qparams.LoadPickle()
        self.nt = self.qparams.nt
        self.nxyz = self.qparams.nxyz
        self.dim_label = 'RC'+str(self.nxyz)+'x'+str(self.nt)
        self.dim_ll = str(self.nxyz)+r'^{3}\times ' + str(self.nt)

        self.SetCustomName(name)

        ## if flowed run
        ##      self.FitDict [itflow , iFF , ifit] fitting class instance
        ## else
        ##      self.FitDict [iff , ifitr ] fitting instance
        # self.FitDict = ODNested()
        # self.FitDictAvg = ODNested()
        # self.FitDictStd = ODNested()
        if self.doflow:
            self.FF_Col_Names = ['flow_time','Chi2pdf_number','Q_momentum']
        else:
            self.FF_Col_Names = ['Chi2pdf_number','Q_momentum']
        self.FF_Stats = pa.DataFrame()
        '''
        self.FF_Stats DataFrame:

        columns =   boot,   Avg,        Std

        rows =     'flow_time' (if flowed), Chi2pdf_number, 'Q_momentum'

        data is Fitting class instance (which includes form factor numbers as parameters)
        '''


        if self.doflow:
            self.FF_Fit_Col_Names = ['flow_time','Chi2pdf_number','form_factor_number','Q_fit_range']
        else:
            self.FF_Fit_Col_Names = ['Chi2pdf_number','form_factor_number','Q_fit_range']
        self.FF_Fit_Stats = pa.DataFrame()
        '''
        self.FF_Fit_Stats DataFrame:

        columns =   boot,   Avg,        Std

        rows =      flow_time(if flowed), Chi2pdf_number, form_factor_number,   fit_range

        data is Fitting class instance, fit parameters are for fitting form factor
        '''



    def LoadPreFF(self,DefWipe=False):
        ## set mass later
        thisfile = 'Pre'+'_'.join([ self.currtype,'q'+''.join(map(str,self.qparams.mom2list)),
                                    'pp'+''.join(map(str,self.ppparams.mom2list)),self.tflows[0]])
        self.thispreFF = FormFact(  self.currtype,self.qparams.mom2list,self.ppparams.mom2list,
                                    Info=self.Info,name=thisfile,tflowlist=self.tflows)


        ## gamma maricies for all operators in current type need to be found
        Info2pt = deepcopy(self.Info)
        # Info2pt['pmom'] = ['0 0 0']
        Info2pt['Interp'] = 'CPEven'
        self.PreC2forMass = TwoPointCorr(cfglist=self.cfglist2pt, Info=Info2pt)
        self.PreC2forMass.SetCustomName('MassForFF','Mass')
        ## we don't care what mass we use here, so if there is a file, just use it (checkcfgs=false)
        self.PreC2forMass.LoadPickle(DefWipe=DefWipe,CheckCfgs=False)

        self.PreGet2ptFit(self.fit2ptin)

        ## using above mass (with incorrect configuration list, see below!), the system of equations is set up in self.thispreFF:
        print('Finding acceptable operator+momentum combinations')
        self.thispreFF.LoadPickle(self.PreC2forMass,thistflowlist = self.tflows,DefWipe=DefWipe)

    def LoadRatios(self,DefWipe=False,Mcore=False):
        ## set up dictionary of RatioCorr instances.
        ## self.thispreFF.keysOpMoms = [ ( icurrOp , ipkey , iqkey ) ]
        self.InfoList = []
        self.namelist = []
        for iOp,ip,iq in self.thispreFF.keysOpMoms():
            thisInfo = deepcopy(self.Info)
            thisInfo['ppmom'] = self.ppparams.stripmom(ip)
            thisInfo['qmom'] = [self.qparams.stripmom(iq)]
            thisProj,thisGamma,FO = self.FormatProjGamma(iOp)
            thisInfo['Observable'] = FO
            thisInfo['gammas'] = [thisGamma]
            thisInfo['Projector'] = thisProj
            thisInfo['AllC2p'] = True
            ## namelist is required to make RatioFuns output be unique for every result.
            ## thisProj not needed, because output has projector in different folders
            self.namelist.append('_'.join(['FFRats',thisProj,thisGamma,iq,FO]))
            # if self.doflow and ('Top' in thisGamma or 'Wein' in thisGamma):
            #     self.namelist[-1] = self.namelist[-1]+'_'+self.tflow[0]
            if self.doflowfull and ('Top' in thisGamma or 'Wein' in thisGamma):
                self.namelist[-1] = self.namelist[-1]+'_RatFull'
            self.InfoList.append(thisInfo)

        # #save memory?
        # del self.thispreFF

        ## set up SetsOfRat with all momenta and operators
        ## Heads up: self.Ratios eventually gets overridden, look at last instance in __init__ to see form of dictionary
        if self.doflow and len(self.tflows) == 0:
            raise IOError("Please set Info['tflowfit'] in information dictionary if working with CP odd form factor extraction")
        print('Setting up Ratio Functions')
        self.Ratios = SetOfRat(cfglist=self.cfglist,InfoDict=self.InfoList,name=self.name,setnames = self.namelist)
        print('Obtaining ratio functions and setting up system of equations')
        self.Ratios.LoadPickleCommonC2(DefWipe=DefWipe,CheckCfgs=True,Mcore=Mcore)
        self.cfglist = self.Ratios.cfglist
        self.stream_name = list(self.Ratios.SetRat.values())[0].stream_name
        self.SetLegLab(self.LegLab,self.stream_name)
        if isinstance(self.cfg_fromfile,str):
            try:
                with open(self.cfg_fromfile,'wb') as f:
                    pickle.dump( self.cfglist , f )
            except:
                raise IOError('Error pickling ' + self.cfg_fromfile)




    def LoadProperFF(self,DefWipe=False):
        ## overwrite C2forMass with version with proper configuration list
        ## all 2point correlators in .Ratios will be identical, so just pick first one
        self.C2forMass = list(self.Ratios.SetRat.values())[0].C2class
        self.Get2ptFit(self.fit2ptin)

        if self.doflow:
            Info2pt = deepcopy(self.Info)
            if self.doAlphafull:
                Info2pt['pmom'] = ['0 0 0']
            Info2pt['Interp'] = 'CPOdd'
            if 'sm2ptRat' in list(Info2pt.keys()):
                sinklab,jsmcoeffs,sm2ptList,InfoList = GetInfoList(Info2pt)
                if self.alphafree:
                    if self.doAlphafull:
                        FOSet = SetOfNNQFull(cfglist=self.cfglist2pt,InfoDict=InfoList)
                    else:
                        FOSet = SetOfNNQ(cfglist=self.cfglist2pt,InfoDict=InfoList)
                    ## this step removes all gauge fields configurations that do not apear in 3-point correlator
                    if isinstance(self.cfglist,dict) and len(list(self.cfglist.keys())):
                        self.cfglist2pt = ReduceCfgs(FOSet.cfglist,self.cfglist)
                        FOSet.ImportCfgList(self.cfglist2pt)
                    else:
                        self.cfglist2pt = FOSet.cfglist
                    # FOSet = SetOfNNQ(cfglist=self.cfglist2pt,InfoDict=InfoList)
                else:
                    if self.doAlphafull:
                        FOSet = SetOfNNQFull(cfglist=self.cfglist,InfoDict=InfoList)
                    else:
                        FOSet = SetOfNNQ(cfglist=self.cfglist,InfoDict=InfoList)

                FOSet.LoadPickle(DefWipe=DefWipe,CheckCfgs=True)


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
                self.C2FOforAlpha = FOSet.CombineCorrs(thisfun,filename=thisname,
                                                LegLab=thisleglab,jsmlab =jsmlab)
                ## for memeory conservation?
                if self.alphafree:
                    if self.doAlphafull:
                        self.cfglist2pt = self.C2FOforAlpha.NNQFull_cfgs
                    else:
                        self.cfglist2pt = self.C2FOforAlpha.NNQ_cfgs
                    if isinstance(self.cfg_fromfile_2pt,str):
                        try:
                            with open(self.cfg_fromfile_2pt,'wb') as f:
                                pickle.dump( self.cfglist2pt , f )
                        except:
                            raise IOError('Error pickling ' + self.cfg_fromfile_2pt)

                # print self.C2FOforAlpha.NNQAvg
                # print self.C2FOforAlpha.AlphaAvg
            else:
                if self.alphafree:
                    if self.doAlphafull:
                        self.C2FOforAlpha = NNQFullCorr(cfglist=self.cfglist, Info=Info2pt)
                    else:
                        self.C2FOforAlpha = NNQCorr(cfglist=self.cfglist, Info=Info2pt)
                else:
                    if self.doAlphafull:
                        self.C2FOforAlpha = NNQFullCorr(cfglist=self.cfglist2pt, Info=Info2pt)
                    else:
                        self.C2FOforAlpha = NNQCorr(cfglist=self.cfglist2pt, Info=Info2pt)
                self.C2FOforAlpha.LoadPickle(DefWipe=DefWipe,WipeData=True,CheckCfgs=True)

            if self.alphafree:
                if self.doAlphafull:
                    self.cfglist2pt = self.C2FOforAlpha.NNQFull_cfgs
                else:
                    self.cfglist2pt = self.C2FOforAlpha.NNQ_cfgs
                if isinstance(self.cfg_fromfile_2pt,str):
                    try:
                        with open(self.cfg_fromfile_2pt,'wb') as f:
                            pickle.dump( self.cfglist2pt , f )
                    except:
                        raise IOError('Error pickling ' + self.cfg_fromfile_2pt)
            self.GetAlphaFit(self.fitAlphain)
        else:
            ## just to make sure alpha is not used, but these parameters are passed in with default values anyway.
            self.tflows = ['none']
            self.C2FOforAlpha = 'PreDef'
        print('Loading proper system of equations.')
        self.thisFF = FormFact( self.currtype,self.qparams.mom2list,self.ppparams.mom2list,
                                Info=self.Info,tflowlist=self.tflows)
        self.thisFF.LoadPickle(self.C2forMass,DefWipe=DefWipe,thistflowlist=self.tflows,thisalpha=self.C2FOforAlpha)
        print('Loading proper system of equations complete.')


    def LineupReducedRF(self,Mcore=defMcore):
        if not hasattr(self,'thisFF'):
            raise EnvironmentError('Please initialse form factors with .LoadProperFF before lining up reduced ffs')
        if not hasattr(self,'Ratios'):
            raise EnvironmentError('Please initialse ratio functions with .LoadRatios before lining up reduced ffs')
        # self.thisFF.GetRedFFParms()
        print('Reducing the system of equations.')
        self.Ratios.ReduceRatios(self.thisFF.combine_keys_forRF,self.thisFF.coeffs_forRF,show_timer=True,Mcore=Mcore)
        print('Reducing the system of equations complete.')

    def PullFits(self):
        ## Rewrite ratios with the fitted results
        ## self.Ratios = { set , States_Fitted_To , Fit_range , ip  } Fitting class instance
        # if self.fitratio == 'From Def':
        #     self.fitratio = defFitAllGamma
            ## don't forget fitratio = fit#-#
        ## output { set , igamma , ip  } Fitting class instance
        this_col,this_fun = 'boot',None
        if '_mal' in self.currtype:
            if hasattr(self,'ralpha_mult_fun'):
                this_col = 'alpha_rat_'+self.ralpha_mult_fun.__name__
                this_fun = self.ralpha_mult_fun
            else:
                this_col = 'alpha_rat_div'
                this_fun = op.truediv
        if isinstance(self.RfitIndicies,bool) and self.RfitIndicies is not False:
            if this_col is not 'boot':
                raise NotImplementedError('will do later')
            self.RatioFits,self.RatioFitsFlow = self.Ratios.GetFitsFromIndex(
                                                    self.RfitIndicies,self.RfitIndicies_flow,
                                                    self.RF_fit_min,self.RF_fit_max,
                                                    thisflowlist=self.tflows,
                                                    tsum_fit=self.Rat_tsum_range)
        elif self.Do_Chi_cutoff:
            if self.Important_FF:
                if self.currtype not in FF_Imp:
                    raise EnvironmentError(self.currtype + 'not set up for Important_FF, see FormFactors.py')
                self.RatioFits,self.RatioFitsFlow = self.Ratios.GetChi2PdfFits(self.RF_fit_min,self.RF_fit_max,
                                                                chi_min=self.Chi2pdf_Min,chi_max=self.Chi2pdf_Max,
                                                                min_fitr=self.RF_fit_range_minimum,
                                                                thisflowlist=self.tflows,
                                                                tsum_fit=self.Rat_tsum_range,
                                                                nchi_threshold=self.nchi_threshold,
                                                                ralpha_info=(self.C2FOforAlpha,this_fun,this_col),
                                                                imp_list=FF_Imp[self.currtype])
            else:
                self.RatioFits,self.RatioFitsFlow = self.Ratios.GetChi2PdfFits(
                                                self.RF_fit_min,self.RF_fit_max,
                                                chi_min=self.Chi2pdf_Min,chi_max=self.Chi2pdf_Max,
                                                min_fitr=self.RF_fit_range_minimum,
                                                thisflowlist=self.tflows,
                                                tsum_fit=self.Rat_tsum_range,
                                                ralpha_info=(self.C2FOforAlpha,this_fun,this_col),
                                                nchi_threshold=self.nchi_threshold)
        else:
            if this_col is not 'boot':
                raise NotImplementedError('will do later')
            self.RatioFits,self.RatioFitsFlow = self.Ratios.GetBestFits(self.RF_fit_min,self.RF_fit_max,
                                                            min_fitr=self.RF_fit_range_minimum,
                                                            thisflowlist=self.tflows,tsum_fit=self.Rat_tsum_range)

        self.RfitIndicies,self.RfitIndicies_flow = self.RatioFits.index,self.RatioFitsFlow.index
        ## get usefull things from self.Ratio before we overwrite it!
        self.RatioFileList = [irat.HumanFile for irat in self.Ratios.SetRat.values()]
        # self.SetLegLab(' '.join(['$'+self.currtype+'$',self.Ratios.SetRat[self.Ratios.SetRat.keys()[0]].LegLab]))

        # print 'DEBUG'
        # if 'boot' in self.RatioFits:
        #     print self.RatioFits['boot']
        # if 'boot' in self.RatioFitsFlow:
        #     print self.RatioFitsFlow['boot']
        ## abit of cleaning up, make the keys line up and get rid of redundant dictionaries
        lboot,lchi,lAvg,lStd,ilist = [],[],[],[],[]
        if self.Reduce_System: this_iter = iter(self.thisFF.keysOpMomsRed())
        else:this_iter = iter(self.thisFF.keysOpMoms())
        shown_error = False
        for iOp,ip,iq in this_iter:
            this_Op = '_'.join(iOp[1:])
            if self.doflow:
                for itflow in self.tflows:
                    if self.flowlab in this_Op:
                        # print self.RatioFitsFlow.loc[(slice(None),this_Op,iq,itflow),:]
                        try:
                            thisdata = self.RatioFitsFlow.loc[(slice(None),this_Op,iq,itflow,slice(None)),'boot']
                        except Exception as err:
                            if 'boot' in self.RatioFitsFlow:
                                if not shown_error:
                                    print(self.RatioFitsFlow['boot'])
                                    print(self.RatioFits['boot'])
                                    shown_error = True
                                print('no FO solutions found for ',this_Op,iq,itflow)
                            else:
                                print((slice(None),this_Op,iq,itflow,slice(None)))
                                print(self.RatioFitsFlow['boot'])
                                print(self.RatioFits['boot'])
                                raise EnvironmentError(str(err)+'\nNo Fits found for RatioFitsFlow')
                            continue
                        # if (slice(None),this_Op,iq,itflow,slice(None)) not in self.RatioFitsFlow:
                        # chidata =
                        for ifitr,fit_data in thisdata.groupby('time_fit_range'):
                            fit_chi_data = self.RatioFitsFlow.loc[(slice(None),this_Op,iq,itflow,ifitr),'chi2pdf']
                            # chidata[(slice(None),ifitr)]
                            if isinstance(fit_data,pa.Series):
                                for ikey,ival in fit_data.items():
                                    if iOp[0]+'_' in ikey[0]:
                                        fit_data = ival
                                        fit_chi_data = fit_chi_data[ikey]
                                        break
                            ilist.append(('_'.join(iOp),ip,iq,itflow,ifitr))
                            lboot.append(fit_data)
                            lchi.append(fit_chi_data)
                            lAvg.append(fit_data.Avg)
                            lStd.append(fit_data.Std)
                    else:
                        # if (slice(None),this_Op,iq,slice(None)) not in self.RatioFits:
                        #     print 'no solutions found for ',this_Op,iq
                        #     continue
                        try:
                            thisdata = self.RatioFits.loc[(slice(None),this_Op,iq,slice(None)),'boot']
                        except Exception as err:
                            if 'boot' in self.RatioFitsFlow:
                                if not shown_error:
                                    print(self.RatioFitsFlow['boot'])
                                    if 'boot' in self.RatioFits:
                                        print(self.RatioFits['boot'])
                                    shown_error = True
                                print('no solutions found for ',this_Op,iq)
                            else:
                                raise EnvironmentError(str(err)+'\nNo Fits found for RatioFits')
                            continue
                        # chidata = self.RatioFits.loc[(slice(None),this_Op,iq,slice(None)),'chi2pdf']
                        for ifitr,fit_data in thisdata.groupby('time_fit_range'):
                            # fit_chi_data = chidata[(slice(None),ifitr)]
                            fit_chi_data = self.RatioFits.loc[(slice(None),this_Op,iq,ifitr),'chi2pdf']
                            if isinstance(fit_data,pa.Series):
                                for ikey,ival in fit_data.items():
                                    if iOp[0]+'_' in ikey[0]:
                                        fit_data = ival
                                        fit_chi_data = fit_chi_data[ikey]
                                        break
                            ilist.append(('_'.join(iOp),ip,iq,itflow,ifitr))
                            lboot.append(fit_data)
                            lchi.append(fit_chi_data)
                            lAvg.append(fit_data.Avg)
                            lStd.append(fit_data.Std)
            else:
                # if (slice(None),this_Op,iq,slice(None)) not in self.RatioFits:
                #     print 'no solutions found for ',this_Op,iq
                #     continue
                try:
                    thisdata = self.RatioFits.loc[(slice(None),this_Op,iq,slice(None)),'boot']
                except Exception as err:
                    if 'boot' in self.RatioFits:
                        if not shown_error:
                            if 'boot' in self.RatioFitsFlow:
                                print(self.RatioFitsFlow['boot'])
                            print(self.RatioFits['boot'])
                            shown_error = True
                        print('no solutions found for ',this_Op,iq)
                    else:
                        raise EnvironmentError(str(err)+'\nNo Fits found for RatioFits')
                    continue
                for ifitr,fit_data in thisdata.groupby('time_fit_range'):
                    fit_chi_data = self.RatioFits.loc[(slice(None),this_Op,iq,ifitr),'chi2pdf']
                    if isinstance(fit_data,pa.Series):
                        for ikey,ival in fit_data.items():
                            if iOp[0]+'_' in ikey[0]:
                                fit_chi_data = fit_chi_data[ikey]
                                fit_data = ival
                                break
                    ilist.append(('_'.join(iOp),ip,iq,ifitr))
                    lboot.append(fit_data)
                    lchi.append(fit_chi_data)
                    lAvg.append(fit_data.Avg)
                    lStd.append(fit_data.Std)
            # # firstgamma = setdata.keys()[0]
            # # firstiq = setdata[firstgamma].keys()[0]
            # # ## its probably 'Par1' since constant fit function, but just to be safe.
            # # thisdata = setdata[firstgamma][firstiq]
            # if isinstance(thisdata,ff.Fitting):
            #     if self.doflow:
            #         for itflow in self.tflows:
            #             data = thisdata[itflow]
            #             rathold[iOp][ip][iq][itflow] = thisdata.Params.values()[0]
            #     else:
            #         rathold[iOp][ip][iq] = thisdata.Params.values()[0]
            # else:
            #     ## we assume all the flow times are the same for every operator
            #     ## (they should be because the info dict only specifies 1 flow time).
            #     for itflow,flowdata in thisdata.iteritems():
            #         rathold[iOp][ip][iq][itflow] = flowdata.Params.values()[0]
        self.Ratios = pa.DataFrame()
        if len(ilist) > 0:
            if self.doflow:
                these_names = ['current_operator','sink_momentum','current_momentum','flow_time','fit_range']
            else:
                these_names = ['current_operator','sink_momentum','current_momentum','fit_range']
            indicies = pa.MultiIndex.from_tuples(ilist,names=these_names)
            self.Ratios.loc[:,'boot'] = pa.Series(lboot,index=indicies)
            self.Ratios.loc[:,'chi2pdf'] = pa.Series(lboot,index=indicies)
            self.Ratios.loc[:,'Avg'] = pa.Series(lAvg,index=indicies)
            self.Ratios.loc[:,'Std'] = pa.Series(lStd,index=indicies)
        # Final form of self.Ratios = { iOp , ip , iq (, itf if has flow time) } Bootstrap
        # Final form of self.Ratios = { iOp , ip , iq } Bootstrap

        ## solves the system of equations to get the form factors out
        ## results stored in self.ffres, self.ffresAvg, self.ffresStd (see SolveRatios())

        self.thisMassboot = self.C2forMass.C2_Fit_Stats.at[(self.fit2pt[0],'p000'),'boot']
        self.thisMassboot = self.thisMassboot[self.fit2pt[1]].iloc[0].fit_data['Params']['Energy']
        self.thisMass = pa.Series([self.thisMassboot],index=['_'.join(self.fit2pt)])
        self.thisMassAvg = self.thisMassboot.Avg
        # mlist,ilist = [],[]
        # for (istate,ifitr),setdata in self.C2forMass.C2_Fit_Stats['boot'][self,'p000',self.fit2pt[1]].iteritems():
        #     # print setdata.keys(),'p000'
        #     # print setdata['p000'].keys(),self.fit2pt[1]
        #     # print setdata['p000'][self.fit2pt[1]].Params.keys(),'Energy'
        #     mlist.append(setdata.)
        #     self.thisMass['_'.join([iset,'_'.join(self.fit2pt)])] = setdata['p000'][self.fit2pt[1]].Params['Energy']
        #     self.thisMassboot = self.thisMass.values()[0]
        #     self.thisMassAvg = self.thisMassboot.Avg




    ## deletes all other data, and just leaves form factor results.
    ## recomended to do this before running analysis of form factor results
    ## to reduce memory usage
    def DelResidual(self):
        if hasattr(self,'Ratios'):
            del self.Ratios
        self.Ratios = {}
        if hasattr(self,'RatioFits'):
            del self.RatioFits
        self.RatioFits = {}
        ## mass might be usefull!
        # del self.C2forMass


    def GetQsqrd(self,iq,ipp,Phys=True):
        ## Q^2 = -q^2 = -q^mu q_mu
        if Phys:
            ## TODO? use thisMassboot instead? but then the yaxis will have error too!
            pmu,qmu,ppmu = CorrSpinTrace().Create4Mom(iq*self.qparams.qunit,ipp*self.qparams.qunit,self.thisMassAvg)
        else:
            pmu,qmu,ppmu = CorrSpinTrace().Create4Mom(iq,ipp,self.thisMassAvg)
        ## meh, dont like negative Qsqrds
        return np.abs(-np.sum(np.array(qmu)**2))
        # return -np.sum(np.array(qmu)**2)

    def MakeQsqrdList(self):
        thisqvec = self.qparams.momlist
        thisppvec = self.ppparams.momlist
        self.Qsqrdlist = np.array([])
        for ipp in thisppvec:
            for iq in thisqvec:
                thisQsqrd = self.GetQsqrd(iq,ipp)
                if all(np.abs(self.Qsqrdlist -thisQsqrd) > 0.0001):
                    self.Qsqrdlist = np.append(self.Qsqrdlist,thisQsqrd)
        self.Qsqrdlist = sorted(self.Qsqrdlist)
        # self.Qsqrdlist = sorted(list(set(Qsqrdlist)))
        self.QsqrdKeyList = ['Q^{2}'+QsqrdToFormat(iQ) for iQ in self.Qsqrdlist]

    def SetupRatiosFO(self,show_timer=True):
        self.chi_Qsqrd = {}

        for iQkey,iQsqrd in zip(self.QsqrdKeyList,self.Qsqrdlist):
            xdata,ydata = [],[]
            if self.Reduce_System:
                this_iter = iter(self.thisFF.itemsParListRed())
            else:
                this_iter = iter(self.thisFF.itemsParList())
            for iOp,ip,iq,imass,itflow,parlabs,parvals in this_iter:
                # print itflowcyc, itflow
                if itflow not in self.tflows: continue
                ## imass dissapears from here on, so we assume only 1 fit range is being used in the calculation
                qvec = self.qparams.TOpvec(self.qparams.stripmom(iq))
                ppvec = self.ppparams.TOpvec(self.qparams.stripmom(ip))
                thisQsqrd = self.GetQsqrd(np.array(qvec),np.array(ppvec))
                # print thisQsqrd,iQsqrd
                # print thisQsqrd, iQsqrd
                if abs(thisQsqrd - iQsqrd) < myeps:
                    ydata.append([])
                    xdata.append(parvals)
                    for ifitr,fitdata in self.Ratios.loc[(iOp,ip,iq,itflow),'boot'].items():
                        # print iOp,ip,iq,parlabs
                        # print [ipar.Avg for ipar in parvals]
                        # print self.Ratios[iOp][ip][iq][itflow].Avg
                        # print thisQsqrd, iQsqrd
                        # print self.Ratios.keys(),iOp
                        # print self.Ratios[iOp].keys(),ip
                        # print self.Ratios[iOp][ip].keys(),iq
                        # print self.Ratios[iOp][ip][itflow].keys(),itflow
                        ydata[-1].append(fitdata)
            # print iQkey
            # for ix,iy in zip(xdata,ydata):
            #     print 'ix',[iix.Avg for iix in ix],'iy',  iy.Avg

            thisparlab = parlabs
            thisfun = self.thisFF.fffun
            thisnpar = self.thisFF.npar
            ### TODO, look through ydata, and find which coeficients are zero. If none, raise exception.
            ## this deals with the case at Q^2=0 where only F1 can be solved in the vector case.
            ## maybe there are other Q^2's and current types where this occurs?
            OnlyGe = False
            if len(xdata) == 0:
                print('No ratio functions found for ' , itflow, iQkey)
                continue
            elif len(xdata) == 1:
                thisfun = FormFactorO1
                thisnpar = 1
                thisparlab = [thisparlab[0]]
                # print xdata
                # print ydata
                xdata = np.array([np.rollaxis(np.array(xdata),1)[0]])
                OnlyGe = True
            # elif 'Tensor' in self.currtype and len(ydata) == 2:
            #     # TODO:
            #     print 'DEBUG: 2 equations for 3 params'
            #     print xdata[0][0].Avg, xdata[0][1].Avg, xdata[0][2].Avg, ydata[0].Avg, ydata[0].Std
            #     print xdata[1][0].Avg, xdata[1][1].Avg, xdata[1][2].Avg, ydata[1].Avg, ydata[1].Std
            #     raise ArithmeticError('Found system of equations less than number of parameters, this needs to be implemented for tensor form factor solutions')
            else:
                xdata = np.rollaxis(np.array(xdata),1)
            # if len(ydata[0]) ==0 :
            #     print 'No ratio functions found for ',iQkey
            #     if show_timer: thistimer.Lap()
            #     continue

            print(iQkey, ' has ',len(list(ydata)),' different momenta')
            avg_fits = np.mean([len(iy) for iy in ydata])
            total_fits = np.sum([len(iy) for iy in ydata])
            print('Average number of fits per momenta is ', avg_fits)
            this_list = list(itertools.product(*ydata))
            print('Total number of fits is ', total_fits,' giving total number of system of equations = ',  len(this_list))
            if len(this_list) > 1000:
                for iinc in range(len(this_list)//1000+2):
                    this_dir = self.SystemDir+'/Chi_block_'+str(iinc)+'/'
                    mkdir_p(this_dir)
            else:
                this_dir = self.SystemDir+'/Chi_block_0/'
                mkdir_p(this_dir)
            # print len(list(itertools.product(*ydata))),' different fit range combinations found for ',iQkey
            # if show_timer: thistimer = Timer(   linklist=range(len(list(itertools.product(*ydata)))),
            #                                     name='Solving Form Factors for '+iQkey)
            this_file = self.SystemDir+iQkey+'_info.dat'
            with open(this_file,'wb') as f:
                pickle.dump([OnlyGe,itflow,len(this_list)],f)
            if itflow not in list(self.chi_Qsqrd.keys()):
                self.chi_Qsqrd[itflow] = {}
            self.chi_Qsqrd[itflow][iQkey] = len(this_list)

            if show_timer: thistimer = Timer(   linklist=list(range(len(this_list))),
                                                name='Setting Form Factors for '+iQkey)
            for icy,iylist in enumerate(itertools.product(*ydata)):
                thisfit = ff.Fitting(Funs=[thisfun,thisnpar],data=[xdata,np.array(iylist)],name=self.name+'_Chi'+str(icy),paramlist=thisparlab)
                mul_of_1000 = icy//1000
                this_dir = self.SystemDir+'/Chi_block_'+str(mul_of_1000)+'/'
                this_file = this_dir+iQkey+'_Chi'+str(icy)+'.dat'
                with open(this_file,'wb') as f:
                    pickle.dump(thisfit,f)
                if show_timer: thistimer.Lap()

        this_file = self.SystemDir+'_Instance.dat'
        with open(this_file,'wb') as f:
            pickle.dump(self.__dict__,f)


    def LoadSystemResultsFO(self,show_timer=True):
        this_qfile = self.SystemDir+'_info.dat'
        with open(this_qfile,'rb') as f:
            self.QsqrdKeyList,self.Qsqrdlist = pickle.load(f)
        ffboot,ffAvg,ffStd,ilist = [],[],[],[]
        for iQkey,iQsqrd in zip(self.QsqrdKeyList,self.Qsqrdlist):
            this_file = self.SystemDir+iQkey+'_info.dat'
            with open(this_file,'rb') as f:
                OnlyGe,itflow,ydata_len = pickle.load(f)

            if show_timer: thistimer = Timer(   linklist=list(range(ydata_len)),
                                                name='Loading Form Factors for '+iQkey)
            for icy in range(ydata_len):
                mul_of_1000 = icy//1000
                this_dir = self.SystemDir+'/Chi_block_'+str(mul_of_1000)+'/'
                this_file = this_dir+iQkey+'_Chi'+str(icy)+'.res'

                if os.path.isfile(this_file):
                    with open(this_file,'rb') as f:
                        thisfit = pickle.load(f)
                else:
                    print('WARNING: FNF: ',this_file)
                # else:
                #     print 'FNF:',this_file

                if OnlyGe and ('Vector' == self.currtype or 'VectorTop' in self.currtype or 'VectorWein' in self.currtype):
                    thisfit.fit_data.loc['G_{E}',:] = thisfit.fit_data.loc['F_{1}',:]
                ffboot.append(thisfit)
                ffAvg.append(thisfit.fit_data['ParamsAvg'])
                ffStd.append(thisfit.fit_data['ParamsStd'])
                ilist.append((itflow,'Chi'+str(icy),iQkey))
                if show_timer: thistimer.Lap()
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.FF_Col_Names)
            self.FF_Stats.loc[:,'boot'] = pa.Series(ffboot,index=indicies)
            self.FF_Stats.loc[:,'Avg'] = pa.Series(ffAvg,index=indicies)
            self.FF_Stats.loc[:,'Std'] = pa.Series(ffStd,index=indicies)
        self.Write()
        if 'BNL' in self.currtype:
            self.ConvertToBNL()
        if 'Top' in self.currtype or 'Wein' in self.currtype:
            self.F3div2M()


    def TOefm(self):
        self.F3div2M()
        if not self.doflow: return
        if self.Inefm: return
        F3 = 'F_{3}'
        if 'BNL' in self.currtype: F3 = 'tilde{F}_{3}'
        F3divm = r'frac{'+F3+r'}{2m_{N}}'
        if not hasattr(self.qparams,'hbarcdivlat'):
            from Params import latspace,hbarc
            self.qparams.hbarcdivlat = hbarc/latspace
        for ikey,Qsqrdres in self.FF_Stats['boot'].items():
            # if F3 in self.ffres[iflow][iQsqrd].Params.keys():
            #     self.ffres[iflow][iQsqrd].Params[F3] = (Qsqrdres.Params[F3]/self.qparams.hbarcdivlat)
            ## this is to correct legacy resuts in which Qsqrdres.Parmams[F3divm] = (Qsqrdres.Params[F3]*self.qparams.hbarc)/(2*thismass)
            ## to                                        Qsqrdres.Parmams[F3divm] = (Qsqrdres.Params[F3]*self.qparams.latspace)/(2*thismass)
            if F3divm in list(Qsqrdres.fit_data['Params'].keys()):
                self.FF_Stats.at[ikey,'boot'].fit_data.at[F3divm,'Params'] = Qsqrdres.fit_data.at[F3divm,'Params']/self.qparams.hbarcdivlat
                self.FF_Stats.at[ikey,'boot'].Stats()
                self.FF_Stats.at[ikey,'Avg'].at[F3divm] = self.FF_Stats.at[ikey,'boot'].fit_data.at[F3divm,'Params'].Avg
                self.FF_Stats.at[ikey,'Std'].at[F3divm] = self.FF_Stats.at[ikey,'boot'].fit_data.at[F3divm,'Params'].Std
        self.Inefm = True


    def F3div2M(self,show_timer=False):
        if not self.doflow: return
        F3 = 'F_{3}'
        if 'BNL' in self.currtype: F3 = 'tilde{F}_{3}'
        F3divm = r'frac{'+F3+r'}{2m_{N}}'
        # thismass = self.C2forMass.C2_Fit_Stats.at[(self.fit2pt[0],'p000',self.fit2pt[1]),'boot']
        # ## ends up being a bootstrapped quantitiy
        # thismass = thismass.fit_data.at['Energy','Params']
        # if not hasattr(self.qparams,'hbarc'):
        #     from Params import hbarc
        #     self.qparams.hbarc = hbarc
        # if not hasattr(self.qparams,'latspace'):
        #     from Params import latspace
        #     self.qparams.latspace = latspace
        if show_timer: thistimer = Timer(   linklist=len(self.FF_Stats['boot']),
                                            name='Creating F3div2M form factors')
        for ikey,Qsqrdres in self.FF_Stats['boot'].items():
            if F3 not in Qsqrdres.fit_data.index:
                if show_timer: thistimer.Lap(ikey)
                continue
            if F3divm in Qsqrdres.fit_data.index:
                if show_timer: thistimer.Lap(ikey)
                continue
            # F3hold = deepcopy(Qsqrdres.Params[F3])
            this_df = pa.DataFrame(columns=Qsqrdres.fit_data.columns)
            for icol,colval in self.FF_Stats.at[ikey,'boot'].fit_data.items():
                for jkey,keyval in colval.items():
                    this_df.at[jkey,icol] = keyval
            this_df.loc[F3divm,:] = Qsqrdres.fit_data.loc[F3,:]
            this_df.at[F3divm,'Params'] = (Qsqrdres.fit_data.at[F3,'Params']*self.qparams.latspace)/(2*self.thisMassboot)
            this_df.at[F3divm,'Params'].Stats()
            self.FF_Stats.at[ikey,'boot'].fit_data = this_df
            self.FF_Stats.at[ikey,'boot'].fit_data.at[F3divm,'Params'].Stats()
            self.FF_Stats.at[ikey,'boot'].Stats()
            self.FF_Stats.at[ikey,'Avg'].at[F3divm] = self.FF_Stats.at[ikey,'boot'].fit_data.at[F3divm,'Params'].Avg
            self.FF_Stats.at[ikey,'Std'].at[F3divm] = self.FF_Stats.at[ikey,'boot'].fit_data.at[F3divm,'Params'].Std
            if show_timer: thistimer.Lap(ikey)
        self.Inefm = True



    def SetupRatios(self,show_timer=True):
        print('Setting up system of equations')
        self.MakeQsqrdList()
        this_qfile = self.SystemDir+'_info.dat'
        with open(this_qfile,'wb') as f:
            pickle.dump([self.QsqrdKeyList,self.Qsqrdlist],f)


        if self.doflow:
            self.SetupRatiosFO(show_timer=show_timer)
            return
        self.chi_Qsqrd = {}
        for iQkey,iQsqrd in zip(self.QsqrdKeyList,self.Qsqrdlist):
            xdata,ydata = [],[]
            if self.Reduce_System:
                this_iter = iter(self.thisFF.itemsParListRed())
            else:
                this_iter = iter(self.thisFF.itemsParList())
            for iOp,ip,iq,imass,parlabs,parvals in this_iter:
                ## imass dissapears from here on, so we assume only 1 fit range is being used in the calculation
                # print 'DEBUG'
                # print iOp,ip,iq,parlabs, parvals
                # print self.Ratios.keys(), iOp
                # print self.Ratios[iOp].keys() , ip
                # print self.Ratios[iOp][ip].keys() , iq
                # print [ipar.Avg for ipar in parvals],self.Ratios[iOp][ip][iq].Avg
                # print thisQsqrd, iQsqrd
                qvec = self.qparams.TOpvec(self.qparams.stripmom(iq))
                ppvec = self.ppparams.TOpvec(self.qparams.stripmom(ip))
                thisQsqrd = self.GetQsqrd(np.array(qvec),np.array(ppvec))
                if abs(thisQsqrd - iQsqrd) < myeps:
                    xdata.append(parvals)
                    # print 'DEBUG'
                    # print self.Ratios.keys(),iOp
                    # print self.Ratios[iOp].keys(),ip
                    # print self.Ratios[iOp][ip].keys(),iq
                    # print self.Ratios[iOp][ip][iq].keys()
                    ydata.append([])
                    if (iOp,ip,iq) not in self.Ratios['boot']: continue
                    for ifitr,fitdata in self.Ratios.loc[(iOp,ip,iq),'boot'].items():
                        ydata[-1].append(fitdata)
            # print 'DEBUG'
            # print iQkey
            # for ix,iy in zip(xdata,ydata):
            #     print 'ix',ix[0].Avg,ix[1].Avg,'iy',  iy.Avg,iy.Std

            thisparlab = parlabs
            thisfun = self.thisFF.fffun
            thisnpar = self.thisFF.npar
            ### TODO, look through ydata, and find which coeficients are zero. If none, raise exception.
            ## this deals with the case at Q^2=0 where only F1 can be solved in the vector case.
            ## maybe there are other Q^2's and current types where this occurs?
            OnlyGe = False
            if len(xdata) == 0:
                print()
                print('No ratio functions found for ' , iQkey)
                continue
            elif len(xdata) == 1:
                thisfun = FormFactorO1
                thisnpar = 1
                thisparlab = [thisparlab[0]]
                # print xdata
                # print ydata
                xdata = np.array([np.rollaxis(np.array(xdata),1)[0]])
                OnlyGe = True
            elif 'Tensor' in self.currtype and len(xdata) == 2:
                # TODO:
                print('DEBUG: 2 equations for 3 params')
                print(xdata[0][0].Avg, xdata[0][1].Avg, xdata[0][2].Avg, ydata[0].Avg, ydata[0].Std)
                print(xdata[1][0].Avg, xdata[1][1].Avg, xdata[1][2].Avg, ydata[1].Avg, ydata[1].Std)
                raise ArithmeticError('Found system of equations less than number of parameters, this needs to be implemented for tensor form factor solutions')
            else:
                xdata = np.rollaxis(np.array(xdata),1)
            if len(xdata) ==0 :
                print()
                print('No ratio functions found for ',iQkey)
                continue
            # print 'debug:'
            # print xdata.shape, np.array(ydata).shape
            # for ix, iy in zip(xdata,ydata):
            #     print ix[0].Avg,iy.Avg,iy.Std
            #     print 'TypeTest',type(iy)
            # print
            print(iQkey, ' has ',len(list(ydata)),' different momenta')
            avg_fits = np.mean([len(iy) for iy in ydata])
            total_fits = np.sum([len(iy) for iy in ydata])
            print('Average number of fits per momenta is ', avg_fits)
            this_list = list(itertools.product(*ydata))
            print('Total number of fits is ', total_fits,' giving total number of system of equations = ',  len(this_list))

            # print len(list(itertools.product(*ydata))),' different fit range combinations found for ',iQkey
            this_file = self.SystemDir+iQkey+'_info.dat'
            with open(this_file,'wb') as f:
                pickle.dump([OnlyGe,len(this_list)],f)
            self.chi_Qsqrd[iQkey] = len(this_list)
            if len(this_list) > 1000:
                for iinc in range(len(this_list)//1000+2):
                    this_dir = self.SystemDir+'/Chi_block_'+str(iinc)+'/'
                    mkdir_p(this_dir)
            else:
                this_dir = self.SystemDir+'/Chi_block_0/'
                mkdir_p(this_dir)

            if show_timer: thistimer = Timer(   linklist=list(range(len(this_list))),
                                                name='Setting up Form Factors for '+iQkey)
            for icy,iylist in enumerate(itertools.product(*ydata)):
                mul_of_1000 = icy//1000
                this_dir = self.SystemDir+'/Chi_block_'+str(mul_of_1000)+'/'
                this_file = this_dir+iQkey+'_Chi'+str(icy)+'.dat'
                if not os.path.isfile(this_file):
                    thisfit = ff.Fitting(Funs=[thisfun,thisnpar],data=[np.array(xdata),np.array(iylist)],name=self.name+'_Chi'+str(icy),paramlist=thisparlab)
                    with open(this_file,'wb') as f:
                        pickle.dump(thisfit,f)
                if show_timer: thistimer.Lap()
        this_file = self.SystemDir+'_Instance.dat'
        with open(this_file,'wb') as f:
            pickle.dump(self.__dict__,f)


    def SolveSystem(self,Mcore = defMcore,show_timer=True):
        # self.MakeQsqrdList()
        def this_fun(this_in):
            ## this part requires pathos version of multiprocessing,
            ## but this downgrades python to 3.6
            # this_job_number = multiprocessing.Process()._identity[0]
            # this_file,ichi = this_in
            # if int(this_job_number)%10 == 0:
            #     dat_len = len(list(glob.glob(this_file.replace(ichi,'*'))))
            #     res_len = len(list(glob.glob(this_file.replace(ichi,'*').replace('.dat','.res'))))
            #     sys.stdout.write(str(res_len*100/dat_len)+ '% done \r')
            #     sys.stdout.flush()
            try:
                with open(this_file,'rb') as f:
                    thisfit = pickle.load(f)
                thisfit.FitBoots()
                with open(this_file.replace('.dat','.res'),'wb') as f:
                    thisfit = pickle.dump(thisfit,f)
            except EOFError:
                print('EOFError',this_file)

        this_qfile = self.SystemDir+'_info.dat'
        if not os.path.isfile(this_qfile):
            raise EnvironmentError('FNF: '+this_qfile+'\n Has FFSetup been completed?')
        with open(this_qfile,'rb') as f:
            this_QsqrdKeyList,this_Qsqrdlist = pickle.load(f)
        for iQkey,iQsqrd in zip(this_QsqrdKeyList,this_Qsqrdlist):
            this_file = self.SystemDir+iQkey+'_info.dat'
            if not os.path.isfile(this_file):
                print('FNF:',this_file)
                continue
            with open(this_file,'rb') as f:
                temp = pickle.load(f)
            file_list = []
            for icy in range(temp[-1]):
                mul_of_1000 = icy//1000
                this_dir = self.SystemDir+'/Chi_block_'+str(mul_of_1000)+'/'
                this_file = this_dir+iQkey+'_Chi'+str(icy)+'.dat'
                if os.path.isfile(this_file) and not os.path.isfile(this_file.replace('.dat','.res')):
                    file_list.append((this_file,str(icy)))
            print('Solving', len(file_list), 'equations for ' , iQkey)
            DoMulticore(this_fun,file_list)
            print()

    def DelSystemResults(self,show_timer=True):
        this_I_file = self.SystemDir+'_Instance.dat'
        if os.path.isfile(this_I_file): os.remove(this_I_file)

        this_qfile = self.SystemDir+'_info.dat'
        if os.path.isfile(this_qfile):
            with open(this_qfile,'rb') as f:
                this_QsqrdKeyList,this_Qsqrdlist = pickle.load(f)

            for iQkey,iQsqrd in zip(this_QsqrdKeyList,this_Qsqrdlist):
                this_dat_file = self.SystemDir+iQkey+'_info.dat'
                if os.path.isfile(this_dat_file):
                    with open(this_dat_file,'rb') as f:
                        temp = pickle.load(f)
                    if show_timer: thistimer = Timer(   linklist=list(range(temp[-1])),
                                                        name='Deleting Form Factors for '+iQkey)

                    for icy in range(temp[-1]):
                        mul_of_1000 = icy//1000
                        this_dir = self.SystemDir+'/Chi_block_'+str(mul_of_1000)+'/'
                        this_file = this_dir+iQkey+'_Chi'+str(icy)+'.dat'
                        if os.path.isfile(this_file):
                            os.remove(this_file)
                        this_file = this_file.replace('.dat','.res')
                        if os.path.isfile(this_file):
                            os.remove(this_file)
                        if show_timer: thistimer.Lap()
                    os.remove(this_dat_file)
            os.remove(this_qfile)

    def LoadSystemResults(self,show_timer=True):
        print('Loading system of equations')

        this_file = self.SystemDir+'_Instance.dat'
        if not os.path.isfile(this_file):
            raise EnvironmentError('FNF: '+this_file + '\n as FFSetup and FFSolve been completed?')
        with open(this_file,'rb') as f:
            self.__dict__.update(pickle.load(f))
        if self.doflow:
            self.LoadSystemResultsFO(show_timer=show_timer)
            return
        ffboot,ffAvg,ffStd,ilist = [],[],[],[]
        # self.MakeQsqrdList()

        for iQkey,iQsqrd in zip(self.QsqrdKeyList,self.Qsqrdlist):

            this_file = self.SystemDir+iQkey+'_info.dat'
            if not os.path.isfile(this_file): continue
            with open(this_file,'rb') as f:
                OnlyGe,ydata_len = pickle.load(f)
            if show_timer: thistimer = Timer(   linklist=list(range(ydata_len)),
                                                name='Loading Form Factors for '+iQkey)

            for icy in range(ydata_len):
                mul_of_1000 = icy//1000
                this_dir = self.SystemDir+'/Chi_block_'+str(mul_of_1000)+'/'
                this_file = this_dir+iQkey+'_Chi'+str(icy)+'.res'

                if os.path.isfile(this_file):
                    try:
                        with open(this_file,'rb') as f:
                            thisfit = pickle.load(f)
                    except Exception as err:
                        print()
                        print(err)
                        print('File:',this_file)
                else:
                    print('WARNING: FNF: ',this_file)

                ## at Q^2 = 0, only G_{E} can be deturmined
                if OnlyGe and ('Vector' == self.currtype):
                    thisfit.fit_data.loc['G_{E}',:] = thisfit.fit_data['F_{1}',:]

                ilist.append(('Chi'+str(icy),iQkey))
                ffboot.append(thisfit)
                ffAvg.append(thisfit.fit_data['ParamsAvg'])
                ffStd.append(thisfit.fit_data['ParamsStd'])
                if show_timer: thistimer.Lap()
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.FF_Col_Names)
            self.FF_Stats.loc[:,'boot'] = pa.Series(ffboot,index=indicies)
            self.FF_Stats.loc[:,'Avg'] = pa.Series(ffAvg,index=indicies)
            self.FF_Stats.loc[:,'Std'] = pa.Series(ffStd,index=indicies)
        self.Write()


    def Undo_Over_Renorm(self):
        ## self.Renromed is just a bool to keep track if you have renormalised the quantity or not
        if 'boot' in self.FF_Stats:
            # self.FF_Stats.loc[:,'boot'] = (self.FF_Stats.loc[:,'boot']*FFRenorm[self.currtype]).apply(DoStats)
            # self.FF_Stats.loc[:,'Avg'] = self.FF_Stats.loc[:,'boot'].apply(DoAvg)
            # self.FF_Stats.loc[:,'Std'] = self.FF_Stats.loc[:,'boot'].apply(DoStd)
            for ikey,iff in self.FF_Stats.loc[:,'boot'].items():
                hold_df = pa.DataFrame( columns=iff.fit_data.columns,
                                        index=iff.fit_data.index)
                for iparam,ival in iff.fit_data.loc[:,'Params'].items():
                    hold_df.loc[iparam,'Params'] = ival/FFRenorm[self.currtype]
                self.FF_Stats.loc[ikey,'boot'].fit_data.update(hold_df)
                self.FF_Stats.loc[ikey,'boot'].Stats()
                self.FF_Stats.loc[ikey,'Avg'].update(iff.fit_data.loc[:,'ParamsAvg'])
                self.FF_Stats.loc[ikey,'Std'].update(iff.fit_data.loc[:,'ParamsStd'])
        if 'boot' in self.FF_Fit_Stats:
            # self.FF_Fit_Stats.loc[:,'boot'] = (self.FF_Fit_Stats.loc[:,'boot']*FFRenorm[self.currtype]).apply(DoStats)
            # self.FF_Fit_Stats.loc[:,'Avg'] = self.FF_Fit_Stats.loc[:,'boot'].apply(DoAvg)
            # self.FF_Fit_Stats.loc[:,'Std'] = self.FF_Fit_Stats.loc[:,'boot'].apply(DoStd)
            for ikey,iff_fit in self.FF_Fit_Stats.loc[:,'boot'].items():
                hold_df_fit = pa.DataFrame( columns=iff_fit.fit_data.columns,
                                            index=iff_fit.fit_data.index)
                for iparam,ival in iff_fit.fit_data.loc[:,'Params'].items():
                    hold_df_fit.loc[iparam,'Params'] = ival/FFRenorm[self.currtype]
                self.FF_Fit_Stats.loc[ikey,'boot'].fit_data.update(hold_df_fit)
                self.FF_Fit_Stats.loc[ikey,'boot'].Stats()
                self.FF_Fit_Stats.loc[ikey,'Avg'].update(iff_fit.fit_data.loc[:,'ParamsAvg'])
                self.FF_Fit_Stats.loc[ikey,'Std'].update(iff_fit.fit_data.loc[:,'ParamsStd'])

    def Renorm(self):
        ## self.Renromed is just a bool to keep track if you have renormalised the quantity or not
        if not self.Renormed:
            if 'boot' in self.FF_Stats:
                # self.FF_Stats.loc[:,'boot'] = (self.FF_Stats.loc[:,'boot']*FFRenorm[self.currtype]).apply(DoStats)
                # self.FF_Stats.loc[:,'Avg'] = self.FF_Stats.loc[:,'boot'].apply(DoAvg)
                # self.FF_Stats.loc[:,'Std'] = self.FF_Stats.loc[:,'boot'].apply(DoStd)
                for ikey,iff in self.FF_Stats.loc[:,'boot'].items():
                    hold_df = pa.DataFrame( columns=iff.fit_data.columns,
                                            index=iff.fit_data.index)
                    for iparam,ival in iff.fit_data.loc[:,'Params'].items():
                        if Debug_Mode:
                            print('DEBUG')
                            print(self.currtype)
                            print(list(FFRenorm.keys()))
                            print()
                            print(iparam,'Params')
                            print(iff.fit_data)
                            print()
                        hold_df.loc[iparam,'Params'] = ival*FFRenorm[self.currtype]
                    self.FF_Stats.loc[ikey,'boot'].fit_data.update(hold_df)
                    self.FF_Stats.loc[ikey,'boot'].Stats()
                    self.FF_Stats.loc[ikey,'Avg'].update(iff.fit_data.loc[:,'ParamsAvg'])
                    self.FF_Stats.loc[ikey,'Std'].update(iff.fit_data.loc[:,'ParamsStd'])
            if 'boot' in self.FF_Fit_Stats:
                # self.FF_Fit_Stats.loc[:,'boot'] = (self.FF_Fit_Stats.loc[:,'boot']*FFRenorm[self.currtype]).apply(DoStats)
                # self.FF_Fit_Stats.loc[:,'Avg'] = self.FF_Fit_Stats.loc[:,'boot'].apply(DoAvg)
                # self.FF_Fit_Stats.loc[:,'Std'] = self.FF_Fit_Stats.loc[:,'boot'].apply(DoStd)
                for ikey,iff_fit in self.FF_Fit_Stats.loc[:,'boot'].items():
                    hold_df_fit = pa.DataFrame( columns=iff_fit.fit_data.columns,
                                                index=iff_fit.fit_data.index)
                    for iparam,ival in iff_fit.fit_data.loc[:,'Params'].items():
                        hold_df_fit.loc[iparam,'Params'] = ival*FFRenorm[self.currtype]
                    self.FF_Fit_Stats.loc[ikey,'boot'].fit_data.update(hold_df_fit)
                    self.FF_Fit_Stats.loc[ikey,'boot'].Stats()
                    self.FF_Fit_Stats.loc[ikey,'Avg'].update(iff_fit.fit_data.loc[:,'ParamsAvg'])
                    self.FF_Fit_Stats.loc[ikey,'Std'].update(iff_fit.fit_data.loc[:,'ParamsStd'])
            self.Renormed = True


    def FormatProjGamma(self,thisoplist):
        if 'P' not in thisoplist[0]:
            raise IOError('Please try to have the Projector in the first element of the oplist ')
        thisProj = thisoplist[0]
        thisGamma = '_'.join(thisoplist[1:])
        if 'Top' in thisGamma :
            thisGamma.replace('Top','')
            FO = 'TopCharge'
        elif 'Wein' in thisGamma:
            thisGamma.replace('Wein','')
            FO = 'Weinberg'
        else:
            FO = ''
        return thisProj,thisGamma,FO

    def Get2ptFit(self,fit2pt):
        if fit2pt == 'From Def':
            if self.C2forMass.name not in list(defFitDict.keys()):
                self.fit2pt = list(defFitDict.values())[0]
            else:
                self.fit2pt = defFitDict[self.C2forMass.name]
        else:
            ## don't forget fit2pt = (state# , fitr#-#)
            self.fit2pt = fit2pt

        # if 'boot' not in self.C2forMass.C2_Fit_Stats or (self.fit2pt[0],'p000',self.fit2pt[1]) not in self.C2forMass.C2_Fit_Stats['boot'].keys():
        print('Fitting two-point correlator for mass', '_'.join(self.fit2pt))
        self.C2forMass.Fit(self.fit2pt[0],[self.fit2pt[1]])
        print('Fitting complete.')
        # elif self.fit2pt[1] not in self.C2forMass.C2_Fit_Stats['boot'][self.fit2pt[0],'p000'].keys():
        #     print 'Fitting two-point correlator for mass', self.fit2pt
        #     self.C2forMass.Fit(self.fit2pt[0],[self.fit2pt[1]])

    def PreGet2ptFit(self,fit2pt):
        if fit2pt == 'From Def':
            if self.PreC2forMass.name not in list(defFitDict.keys()):
                self.fit2pt = list(defFitDict.values())[0]
            else:
                self.fit2pt = defFitDict[self.PreC2forMass.name]
        else:
            ## don't forget fit2pt = (state# , fitr#-#)
            self.fit2pt = fit2pt

        # if 'boot' not in self.C2forMass.C2_Fit_Stats or (self.fit2pt[0],'p000',self.fit2pt[1]) not in self.C2forMass.C2_Fit_Stats['boot'].keys():
        print('Fitting Pre two-point correlator for mass', '_'.join(self.fit2pt))
        self.PreC2forMass.Fit(self.fit2pt[0],[self.fit2pt[1]])
        # elif self.fit2pt[1] not in self.C2forMass.C2_Fit_Stats['boot'][self.fit2pt[0],'p000'].keys():
        #     print 'Fitting two-point correlator for mass', self.fit2pt
        #     self.C2forMass.Fit(self.fit2pt[0],[self.fit2pt[1]])


    def GetAlphaFit(self,fitAlpha):
        print('Extracting alpha fit from', self.C2FOforAlpha.name)
        if fitAlpha == 'From Def':
            self.fitAlpha = defFitAlpha[self.C2FOforAlpha.name]
        else:
            ## don't forget fitAlpha = fitr#-#
            self.fitAlpha = fitAlpha
        if self.doAlphafull:
            # if  'boot' not in self.C2FOforAlpha.NNQ_Fit_Stats or \
            #     ('p000',self.fitAlpha) not in self.C2FOforAlpha.NNQFull_Fit_Stats['boot'].swaplevel(1,2):
            print('Fitting NNQFull for Alpha', self.fitAlpha)
            fitr_list = self.fitAlpha.split('_')
            self.C2FOforAlpha.FitAlpha(fit_range=fitr_list[0],thistflowlist = self.tflows,tsum_ranges=fitr_list[1])
        else:
            # if  'boot' not in self.C2FOforAlpha.NNQ_Fit_Stats or \
            #     ('p000',self.fitAlpha) not in self.C2FOforAlpha.NNQ_Fit_Stats['boot'].swaplevel(1,2):
            print('Fitting NNQ for Alpha', self.fitAlpha)
            self.C2FOforAlpha.FitAlpha(self.fitAlpha,thistflowlist =self.tflows)

    def SetCustomName(self,name= '',LegLab = '' ):
        if name == '':
            if self.doflow:
                if self.doflowfull:
                    self.filename = '_'.join((  self.currtype,self.tflows[0],'R'+self.Rat_tsum_range,
                                                'A'+self.fitAlphain,'RFmin'+str(self.RF_fit_range_minimum)))
                    self.name = '_'.join((  self.DS,self.currtype,'qmax'+str(self.ppparams.mom2max),
                                            'ppmax'+str(self.qparams.mom2max),self.tflows[0],
                                            'A'+self.fitAlphain,'RFmin'+str(self.RF_fit_range_minimum)))
                else:
                    self.filename = '_'.join((  self.currtype,self.tflows[0],'A'+self.fitAlphain,
                                                'RFmin'+str(self.RF_fit_range_minimum)))
                    self.name = '_'.join((  self.DS,self.currtype,'qmax'+str(self.ppparams.mom2max),
                                            'ppmax'+str(self.qparams.mom2max),self.tflows[0],
                                            'A'+self.fitAlphain,'RFmin'+str(self.RF_fit_range_minimum)))

            else:
                self.filename = '_'.join((self.currtype,'RFmin'+str(self.RF_fit_range_minimum)))
                self.name = '_'.join(   (self.DS,self.currtype,'qmax'+str(self.ppparams.mom2max),
                                        'ppmax'+str(self.qparams.mom2max),'RFmin'+str(self.RF_fit_range_minimum)))
            if self.Do_Chi_cutoff:
                if self.Important_FF:
                    if self.Reduce_System:
                        self.filename += '_RIChi{:}-{:}'.format(int(100*self.Chi2pdf_Min),int(100*self.Chi2pdf_Max))
                        self.name += '_RIChi{:}-{:}'.format(int(100*self.Chi2pdf_Min),int(100*self.Chi2pdf_Max))
                    else:
                        self.filename += '_IChi{:}-{:}'.format(int(100*self.Chi2pdf_Min),int(100*self.Chi2pdf_Max))
                        self.name += '_IChi{:}-{:}'.format(int(100*self.Chi2pdf_Min),int(100*self.Chi2pdf_Max))
                else:
                    if self.Reduce_System:
                        self.filename += '_RChi{:}-{:}'.format(int(100*self.Chi2pdf_Min),int(100*self.Chi2pdf_Max))
                        self.name += '_RChi{:}-{:}'.format(int(100*self.Chi2pdf_Min),int(100*self.Chi2pdf_Max))
                    else:
                        self.filename += '_Chi{:}-{:}'.format(int(100*self.Chi2pdf_Min),int(100*self.Chi2pdf_Max))
                        self.name += '_Chi{:}-{:}'.format(int(100*self.Chi2pdf_Min),int(100*self.Chi2pdf_Max))
            else:
                if self.Reduce_System:
                    self.filename += '_BestRChi'
                    self.name += '_BestRChi'
                else:
                    self.filename += '_BestChi'
                    self.name += '_BestChi'
            if hasattr(self,'isCustomRFitR') and self.isCustomRFitR:
                self.filename += '_RfitRCust'
                self.name += '_RfitRCust'
            if self.Rat_alpha:
                self.name += '_Rata'
                self.filename += '_Rata'
        else:
            self.name = self.filename = name
        self.SetLegLab(LegLab=LegLab)
        self.subdir = '_'.join((self.DS,'qmax'+str(self.ppparams.mom2max),'ppmax'+str(self.qparams.mom2max)))
        # self.SubDir = '_'.join([self.C3class.Proj,self.C3class.DS,self.C3class.ppform])

        self.FFdir = outputdir + '/'+self.ppparams.kappafolder+'/FormFactors/'+self.subdir+'/'
        self.FFdir_scratch = scratchdir + '/'+self.ppparams.kappafolder+'/FormFactors/'+self.subdir+'/'
        mkdir_p(self.FFdir+'/Pickle/')
        mkdir_p(self.FFdir+'/Excel/')
        self.HumanFile = self.FFdir+self.filename+'.xml'
        self.ExcelFile = self.FFdir+'/Excel/'+self.filename+'.xlsx'
        self.PickleFile = self.FFdir+'/Pickle/'+self.filename+'.py3p'
        self.SystemDir = self.FFdir_scratch+'/System/'+self.filename+'/'
        mkdir_p(self.SystemDir)


    def SetLegLab(self,LegLab = '',stream_name=''):
        if LegLab == '':
            self.LegLab = '$'+'\ '.join([self.currtype,self.DS,self.dim_ll,self.qparams.GetPionMassLab(),stream_name])+'$'
        else:
            if len(stream_name) > 0:
                self.LegLab = LegLab+'$\ '+stream_name+'$'
            else:
                self.LegLab = LegLab
        if self.Rat_alpha:
            self.LegLab += r'\ R_{\alpha}'

    def ConvertToBNL(self,show_timer=True):
        if not self.doflow: return
        if hasattr(self.C2FOforAlpha,'NNQFull_Fit_Stats'):
            # thisalpha = self.C2FOforAlpha.NNQFull_Fit_Stats.loc[('p000',self.tflows[0],self.fitAlpha),'boot']
            fit_split = self.fitAlpha.split('_')
            fit_split[1] = fit_split[1].replace('tsumfitr','fittwor')
            thisalpha = list(self.thisFF.alpha.values())[0].values()[0]
            # thisalpha = self.C2FOforAlpha.NNQFull_Fit_Stats.loc[('p000',self.tflows[0]),'boot']\
            #                                                 [fit_split[0],fit_split[1]].iloc[0]
        else:
            thisalpha = list(self.thisFF.alpha.values())[0].values()[0]
            # thisalpha = self.C2FOforAlpha.NNQ_Fit_Stats.loc[('p000',self.tflows[0]),'boot'][self.fitAlpha].iloc[0]
        ## ends up being a bootstrapped quantitiy

        # thisalpha = thisalpha.fit_data['Params'].values[0]
        twocos_alpha = (2*thisalpha).cos()
        twosin_alpha = (2*thisalpha).sin()
        mintwosin_alpha = (-2*thisalpha).sin()
        if show_timer: thistimer = Timer(   linklist=len(self.FF_Stats['boot']),
                                            name='Creating BNL form factors')
        for ikey,Qsqrdres in self.FF_Stats.iterrows():
            if ('tilde{F}_{2}' not in list(Qsqrdres['boot'].fit_data['Params'].keys())) or ('tilde{F}_{3}' not in list(Qsqrdres['boot'].fit_data['Params'].keys())): continue
            # F2hold = deepcopy(Qsqrdres.fit_data.loc['tilde{F}_{2}','Params'])
            # F3hold = deepcopy(Qsqrdres.fit_data.loc['tilde{F}_{3}','Params'])
            F2hold = Qsqrdres['boot'].fit_data.at['tilde{F}_{2}','Params']
            F3hold = Qsqrdres['boot'].fit_data.at['tilde{F}_{3}','Params']
            hold_df = pa.DataFrame(columns=Qsqrdres['boot'].fit_data.columns)
            for icol,colval in Qsqrdres['boot'].fit_data.items():
                for jkey,ival in colval.items():
                    hold_df.at[jkey,icol] = ival
            # hold_df.loc['tilde{F}_{2}',:] = self.FF_Stats.loc[ikey,'boot'].fit_data.loc['tilde{F}_{2}',:]
            # hold_df.loc['tilde{F}_{3}',:] = self.FF_Stats.loc[ikey,'boot'].fit_data.loc['tilde{F}_{3}',:]
            hold_df.at['tilde{F}_{2}','Params'] = twocos_alpha * F2hold + mintwosin_alpha*F3hold
            hold_df.at['tilde{F}_{3}','Params'] = twosin_alpha * F2hold + twocos_alpha*F3hold
            Qsqrdres['boot'].fit_data = hold_df
            # self.FF_Stats.at[ikey,'boot'].fit_data.update(hold_df)
            Qsqrdres['boot'].Stats()
            Qsqrdres['Avg'].at['tilde{F}_{2}'] = Qsqrdres['boot'].fit_data.at['tilde{F}_{2}','Params'].Avg
            Qsqrdres['Avg'].at['tilde{F}_{3}'] = Qsqrdres['boot'].fit_data.at['tilde{F}_{3}','Params'].Avg
            Qsqrdres['Std'].at['tilde{F}_{2}'] = Qsqrdres['boot'].fit_data.at['tilde{F}_{2}','Params'].Std
            Qsqrdres['Std'].at['tilde{F}_{3}'] = Qsqrdres['boot'].fit_data.at['tilde{F}_{3}','Params'].Std
            if show_timer: thistimer.Lap(ikey)

    def FixRead(self,readdict):
        if 'Renormed' not in list(readdict.keys()):
            readdict['Renormed'] = False
        if 'Inefm' not in list(readdict.keys()):
            readdict['Inefm'] = False
        # if 'RenormedFix' not in readdict.keys():
        #     readdict['RenormedFix'] = False
        if not isinstance(readdict['tflows'],(list,tuple,np.ndarray)):
            readdict['tflows'] = self.tflows
        # if not os.path.isfile(str(readdict['PickleFile'])):
        #     readdict['PickleFile'] = self.PickleFile
        #     readdict['HumanFile'] = self.HumanFile
        print('THIS LINE IS HERE TO MAKE THE FILE POINT TO ITSELF!')
        readdict['PickleFile'] = self.PickleFile
        readdict['HumanFile'] = self.HumanFile
        if hasattr(readdict['ppparams'],'PickleFile') and not os.path.isfile(readdict['ppparams'].PickleFile):
            readdict['ppparams'].PickleFile = self.ppparams.PickleFile
            readdict['ppparams'].HumanFile = self.ppparams.HumanFile
            readdict['qparams'].PickleFile = self.qparams.PickleFile
            readdict['qparams'].HumanFile = self.qparams.HumanFile
            readdict['Info']['outdir'] = self.Info['outdir']
        if 'QsqrdKeylist' in list(readdict.keys()):
            readdict['QsqrdKeyList']  = readdict['QsqrdKeylist']
        if 'isCustomRFitR' not in list(readdict.keys()):
            readdict['isCustomRFitR']  = False
        if Debug_Mode:
            print('Debug mode fixes broken file for currtype, fixing now')
            readdict['currtype'] = self.currtype

        return readdict


    def SetupPickle(self,DefWipe=False,HardWipe=False,Mcore=defMcore,
                    RatioFitLists=[False,False]):
        # print 'Loading Pickle for ' , self.PickleFile
        # print 'Loading Pickle for ' , self.PickleFile

        if not hasattr(self,'isCustomRFitR') or self.isCustomRFitR:
            if self.doflow:
                self.RfitIndicies,self.RfitIndicies_flow = RatioFitLists
            else:
                self.RfitIndicies,self.RfitIndicies_flow = RatioFitLists[0],False
            self.isCustomRFitR = not isinstance(self.RfitIndicies,bool)
        self.SetCustomName()

        if os.path.isfile(self.PickleFile) and not DefWipe:
            thispickle = ReadPickleWrap(self.PickleFile)
            thispickle = self.FixRead(thispickle)

            if self.doflow and self.tflows[0] != thispickle['tflows'][0]:
                if os.path.isfile(self.PickleFile+'.bak'):
                    self.PickleFile = self.PickleFile+'.bak'
                    self.HumanFile = self.HumanFile+'.bak'
                    self.ExcelFile = self.ExcelFile.replace('.xlsx','.bak.xlsx')
                    return self.SetupPickle(DefWipe=DefWipe,HardWipe=HardWipe,Mcore=Mcore)
                print('flow times incorrect in file:')
                print(self.PickleFile)
                print('file is ' , thispickle['tflows'][0])
                print('class is ' , self.tflows[0])
                if os.path.isfile(self.SystemDir+'_Instance.dat'):
                    print('Instance file for computing system of equations found:')
                    print(self.SystemDir+'_Instance.dat')
                    print('Jumping to solving system of equations')
                    print()
                else:
                    self.LoadPreFF(DefWipe=HardWipe)
                    self.LoadRatios(DefWipe=HardWipe,Mcore=Mcore)
                    self.LoadProperFF(DefWipe=HardWipe)
                    if self.Reduce_System: self.LineupReducedRF(Mcore=Mcore)
                    self.PullFits()
                    self.SetupRatios()
                    return self.RfitIndicies,self.RfitIndicies_flow
            if 'RfitIndicies' in list(thispickle.keys()) and (not self.doflow or 'RfitIndicies_flow' in list(thispickle.keys())):
                return thispickle['RfitIndicies'],thispickle['RfitIndicies_flow']
            else:
                return None,None
        else:
            print('Cannot find file',self.PickleFile)
            if os.path.isfile(self.SystemDir+'_Instance.dat') and not DefWipe:
                print('Instance file for computing system of equations found:')
                print(self.SystemDir+'_Instance.dat')
                print('Jumping to solving system of equations')
                print()
            else:
                self.LoadPreFF(DefWipe=HardWipe)
                self.LoadRatios(DefWipe=HardWipe,Mcore=Mcore)
                self.LoadProperFF(DefWipe=HardWipe)
                if self.Reduce_System: self.LineupReducedRF(Mcore=Mcore)
                self.PullFits()
                self.SetupRatios()
            return self.RfitIndicies,self.RfitIndicies_flow

    def SolvePickle(self,DefWipe=False,HardWipe=False,Mcore=defMcore,
                    RatioFitLists=[False,False]):
        # print 'Loading Pickle for ' , self.PickleFile
        if not hasattr(self,'isCustomRFitR') or self.isCustomRFitR:
            if self.doflow:
                self.RfitIndicies,self.RfitIndicies_flow = RatioFitLists
            else:
                self.RfitIndicies,self.RfitIndicies_flow = RatioFitLists[0],False
            self.isCustomRFitR = not isinstance(self.RfitIndicies,bool)
        self.SetCustomName()

        if os.path.isfile(self.PickleFile) and not DefWipe:
            thispickle = ReadPickleWrap(self.PickleFile)
            thispickle = self.FixRead(thispickle)
            if self.doflow and self.tflows[0] != thispickle['tflows'][0]:
                if os.path.isfile(self.PickleFile+'.bak'):
                    print('using backupfile for '+self.name)
                    self.PickleFile = self.PickleFile+'.bak'
                    self.HumanFile = self.HumanFile+'.bak'
                    self.ExcelFile = self.ExcelFile.replace('.xlsx','.bak.xlsx')
                    return self.SolvePickle(DefWipe=DefWipe,HardWipe=HardWipe,Mcore=Mcore)
                print('flow times incorrect in file:')
                print(self.PickleFile)
                print('file is ' , thispickle['tflows'][0])
                print('class is ' , self.tflows[0])
                self.SolveSystem(Mcore=Mcore)
                return self.RfitIndicies,self.RfitIndicies_flow
            if 'RfitIndicies' in list(thispickle.keys()) and (not self.doflow or 'RfitIndicies_flow' in list(thispickle.keys())):
                return thispickle['RfitIndicies'],thispickle['RfitIndicies_flow']
            else:
                return None,None
        else:
            print('Cannot find file',self.PickleFile)
            self.SolveSystem(Mcore=Mcore)
            return self.RfitIndicies,self.RfitIndicies_flow

    def LoadPickle( self,DefWipe=False,HardWipe=False,Mcore=defMcore,OnlyFits=False,DoFits=True,
                    RatioFitLists=[False,False]):
        print('Loading Pickle for ' , self.PickleFile)
        if not hasattr(self,'isCustomRFitR') or self.isCustomRFitR:
            if self.doflow:
                self.RfitIndicies,self.RfitIndicies_flow = RatioFitLists
            else:
                self.RfitIndicies,self.RfitIndicies_flow = RatioFitLists[0],False
            self.isCustomRFitR = not isinstance(self.RfitIndicies,bool)
        self.SetCustomName()

        if os.path.isfile(self.PickleFile) and not DefWipe:
            thispickle = ReadPickleWrap(self.PickleFile)
            thispickle = self.FixRead(thispickle)
            if (not self.doflow) or self.tflows[0] == thispickle['tflows'][0]:
                self.__dict__.update(thispickle)
                ## if the file name is broken, just use default location.
                if 'FF_Stats_name' in list(thispickle.keys()) and thispickle['FF_Stats_name'] is None:
                    self.FF_Stats_name = self.PickleFile+'.FFstats'
                if 'FF_Fit_Stats_name' in list(thispickle.keys()) and thispickle['FF_Fit_Stats_name'] is None:
                    self.FF_Fit_Stats_name = self.PickleFile+'.FFFitstats'
                if os.path.isfile(self.FF_Stats_name) and not OnlyFits:
                    self.FF_Stats = pa.read_pickle(self.FF_Stats_name)
                    print('Reading stats ')
                    print(self.FF_Stats_name)
                if os.path.isfile(self.FF_Fit_Stats_name) and not Wipe_All_Fits:
                    print('Reading fit stats ')
                    print(self.FF_Fit_Stats_name)
                    self.FF_Fit_Stats = pa.read_pickle(self.FF_Fit_Stats_name)
                print('Loaded results from: \n',self.PickleFile)
                if not OnlyFits:
                    # if 'BNL' in self.currtype:
                    #     self.ConvertToBNL()
                    # if 'Top' in self.currtype or 'Wein' in self.currtype:
                    #     self.F3div2M()
                    self.TOefm()
                    self.Renorm()
                    if DoFits: self.FitFFs(show_timer=True)
                    self.Write()
            else:
                if os.path.isfile(self.PickleFile+'.bak'):
                    print('using backupfile for '+self.name)
                    self.PickleFile = self.PickleFile+'.bak'
                    self.HumanFile = self.HumanFile+'.bak'
                    self.ExcelFile = self.ExcelFile.replace('.xlsx','.bak.xlsx')
                    self.LoadPickle(DefWipe=DefWipe,HardWipe=HardWipe,Mcore=Mcore)
                    return
                print('flow times incorrect in file:')
                print(self.PickleFile)
                print('file is ' , thispickle['tflows'][0])
                print('class is ' , self.tflows[0])
                self.LoadSystemResults()
                self.Renorm()
                if DoFits: self.FitFFs(show_timer=True)
                self.Write()
        else:
            print('Cannot find file',self.PickleFile)
            self.LoadSystemResults()
            self.Renorm()
            if DoFits: self.FitFFs(show_timer=True)
            self.Write()

    def GetProj(self,thislist):
        for ilist in thislist:
            if 'P' == ilist[0]:
                return ilist
        return ''

    ## FFout = [ imass , ipp , ( , itflow) thisqsqrd,'_'.join([gammakey,iq]) ] = ( Coeff_1 FF_1 + ... + Coeff_n + FF_n )
    ## Ratout { set , igamma , iq  } Fitting class instance
    ## output is same format as FFout
    ## calculated in a way where set actually contains both gamma and ip, and set has igamma and ip within it
    ## so we just pull the first key for igamma and ip
    def GetEquationsOut(self,FFout,Ratout,Ratout_flow):
        output = ODNested()
        if self.doflow:
            for imass,massdata in FFout.items():
                for ipp,ppdata in massdata.items():
                    for itflow,flowdata in ppdata.items():
                        for iqsqrd,psqrdata in flowdata.items():
                            for ipop,popdata in psqrdata.items():
                                thisProj,thisGamma,FO = self.FormatProjGamma(ipop)
                                FFkey = '_'.join([self.GetProj(ipop.split('_')),self.DS,ipp,'FFRats',ipop,''])
                                if 'tilde' in FFkey: FFkey = FFkey.replace('tilde',r'\tilde')
                                if 'frac' in FFkey: FFkey = FFkey.replace('frac',r'\frac')
                                if self.flowlab in ipop: FFkey += self.flowtype
                                if self.doflowfull and ('Top' in FFkey or 'Wein' in FFkey): FFkey += '_RatFull'
                                # thisRat = Ratout[FFkey]
                                # while not isinstance(thisRat,ff.Fitting):
                                #     try:
                                #         thisRat = thisRat[thisRat.keys()[0]]
                                #     except:
                                #         print thisRat
                                #         raise EnvironmentError('Error occured when trying to get ratio function for equations out')
                                # print Ratout.keys(),FFkey
                                # raise EnvironmentError('This part needs to be implemented.')
                                # if FFkey not in Ratout.keys(): continue
                                # thisRat = thisRat[thisRat.keys()[0]]
                                if any([flowed_op in FFkey for flowed_op in ObsListAbr]):
                                    if 'boot' not in Ratout_flow:
                                        raise EnvironmentError('Ratout_flow not initialised')
                                    elif FFkey not in list(Ratout_flow['boot'].keys()):
                                        err_hold = [', '.join(ival) for ival in list(Ratout_flow['boot'].keys())]
                                        raise EnvironmentError(FFkey + " not in Ratout.keys(): \n" + '\n'.join(err_hold))
                                    thisRat = Ratout_flow['boot'][FFkey].values[0]
                                else:
                                    if 'boot' not in Ratout:
                                        raise EnvironmentError('Ratout not initialised')
                                    elif FFkey not in list(Ratout['boot'].keys()):
                                        err_hold = [', '.join(ival) for ival in list(Ratout_flow['boot'].keys())]
                                        raise EnvironmentError(FFkey + " not in Ratout.keys() \n" + '\n'.join(err_hold))
                                    thisRat = Ratout['boot'][FFkey].values[0]
                                output[imass][ipp][itflow][iqsqrd][ipop] = popdata + ' =  ('+AvgStdToFormat(thisRat.Avg,thisRat.Std,dec=5,numlen=8)+ ') '
                                # print Ratout[Ratout.keys()[0]], ipop
                                # print imass,ipp,ipop
        else:
            for imass,massdata in FFout.items():
                for ipp,ppdata in massdata.items():
                    for iqsqrd,psqrdata in ppdata.items():
                        for ipop,popdata in psqrdata.items():
                            FFkey = '_'.join([self.GetProj(ipop.split('_')),self.DS,ipp,'FFRats',ipop,''])
                            if 'boot' not in Ratout:
                                raise EnvironmentError('Ratout not initialised')
                            elif FFkey not in list(Ratout['boot'].keys()):
                                raise EnvironmentError(FFkey + " not in Ratout.keys() " + ', '.join(list(Ratout['boot'].keys())))
                            thisRat = Ratout['boot'][FFkey].iloc[0]
                            output[imass][ipp][iqsqrd][ipop] = popdata + ' =  ('+AvgStdToFormat(thisRat.Avg,thisRat.Std,dec=5,numlen=8)+ ') '
                            # print Ratout[Ratout.keys()[0]], ipop
                            # print imass,ipp,ipop
        return output

    def GetFuns(self):
        self.thisFF.GetFuns()
        if hasattr(self,'thispreFF'):
            self.thispreFF.GetFuns()
        self.C2FOforAlpha.GetFuns()
        self.PreC2forMass.GetFuns()
        self.C2forMass.GetFuns()
        # if not hasattr(self,'Fun'):
        #     self.Fun = ReadFuns(self.Fun_name)[0]

    def RemoveFuns(self):
        self.thisFF.RemoveFuns()
        if hasattr(self,'thispreFF'):
            self.thispreFF.RemoveFuns()
        self.C2FOforAlpha.RemoveFuns()
        if hasattr(self,'PreC2forMass'):
            self.PreC2forMass.RemoveFuns()
        self.C2forMass.RemoveFuns()
        # if hasattr(self,'Fun'):
        #     self.Fun_name = self.Fun.__name__
        #     WriteFuns(self.Fun)
        #     del self.Fun

    def Write(self,WipeData=True,WipeAllFits=False):
        def FixDictArray(thislist,ikey):
            outDict = OrderedDict()
            for ic,ilist in enumerate(thislist):
                outDict[ikey+str(ic+1)] = ilist
            return outDict
        try:
        # print 'DEBUG, please uncomment this'
            if not hasattr(self,'preEq'):
                if self.Reduce_System:
                    preEq = self.GetEquationsOut(self.thisFF.outDict['Equations_Reduced'],self.RatioFits,self.RatioFitsFlow)
                else:
                    preEq = self.GetEquationsOut(self.thisFF.outDict['Equations'],self.RatioFits,self.RatioFitsFlow)
        # except Exception as err:
        except:
            # print err
            if hasattr(self,'outDict') and isinstance(self.outDict,(dict,OrderedDict)):
                if 'Equations' in list(self.outDict.keys()):
                    preEq = self.outDict['Equations']
                else:
                    preEq = 'No Equations'
            else:
                preEq = 'No outDict'

        self.outDict = ODNested()

        self.outDict['nt'] = self.nt
        self.outDict['nxyz'] = self.nxyz
        self.outDict['name'] = self.name
        self.outDict['Mass_File'] = self.C2forMass.HumanFile
        if self.doflow: self.outDict['Alpha_File'] = self.C2FOforAlpha.HumanFile
        self.outDict['Coeff_File'] = self.thisFF.HumanFile
        self.outDict['Ratio_Files'] = FixDictArray(self.RatioFileList,'Rat')
        if self.Renormed:
            self.outDict['Renorm'] = FFRenorm[self.currtype]
        else:
            self.outDict['Renorm'] = 'None'
        self.outDict['F_{3}/2m_in_efm'] = self.Inefm
        self.outDict['nchi_fit_threshold'] = self.nchi_fit_threshold
        self.outDict['nchi_threshold'] = self.nchi_threshold
        if '_mal' in self.currtype:
            self.outDict['ralpha_mult_fun'] = self.ralpha_mult_fun.__name__

        if self.doflow:
            for (itflow,ichi,iQ),Qdata in self.FF_Stats['boot'].items():
                self.outDict[itflow][ichi][iQ] = Qdata.GetOutputDict()
            if 'boot' in self.FF_Fit_Stats:
                for (itflow,ichi,iFF,ifitr),fitdata in self.FF_Fit_Stats['boot'].items():
                    self.outDict[itflow][ichi]['FF_Fits'+ifitr][iFF] = fitdata.GetOutputDict()

        else:
            for (ichi,iQ),Qdata in self.FF_Stats['boot'].items():
                self.outDict[ichi][iQ] = Qdata.GetOutputDict()
            if 'boot' in self.FF_Fit_Stats:
                for (ichi,iFF,ifitr),fitdata in self.FF_Fit_Stats['boot'].items():
                    self.outDict[ichi]['FF_Fits'+ifitr][iFF] = fitdata.GetOutputDict()

        bakThis = False
        if preEq == 'No outDict':
            print('Cannot pull equations from previous result, omitting from output (this occurs when using older results)')
            bakThis = True
        elif preEq == 'No Equations':
            print('No Equations found in previous file. (this occurs when using older results)')
        else:
            self.outDict['Equations'] = preEq

        if WipeData: self.DelResidual()
        ## Write human readable file with data
        if bakThis: BackupFile(self.HumanFile)
        ## pickles rest of data for reading
        if WipeAllFits: self.WipeAllFits()
        self.RemoveFuns()
        if hasattr(self,'FF_Stats'):
            self.FF_Stats_name = self.PickleFile+'.FFstats'
            if os.path.isfile(self.FF_Stats_name):
                os.rename(self.FF_Stats_name,self.FF_Stats_name+'.bak')
            self.FF_Stats.to_pickle(self.FF_Stats_name)
        if hasattr(self,'FF_Fit_Stats'):
            self.FF_Fit_Stats_name = self.PickleFile+'.FFFitstats'
            if os.path.isfile(self.FF_Fit_Stats_name):
                os.rename(self.FF_Fit_Stats_name,self.FF_Fit_Stats_name+'.bak')
            self.FF_Fit_Stats.to_pickle(self.FF_Fit_Stats_name)
        out_dict = deepcopy(self.__dict__)
        if 'FF_Stats' in list(out_dict.keys()):
            del out_dict['FF_Stats']
        if 'FF_Fit_Stats' in list(out_dict.keys()):
            del out_dict['FF_Fit_Stats']
        WritePickle(self.PickleFile,out_dict)
        self.GetFuns()

        WriteXml(self.HumanFile,{'Results':self.outDict})
        self.DelSystemResults()


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

    def CheckFF(self,thisFF,resub=True):
        thiskey = list(self.FF_Stats['boot'].keys())[-1]
        if self.doflow:
            calcfflist = list(self.FF_Stats['boot'][thiskey].fit_data['Params'].keys())
            if thisFF not in calcfflist:
                if ('tilde{F}_{3}' in thisFF and 'tilde{F}_{3}' not in calcfflist[-1] and self.qparams.mom2max != 0) and resub:
                    self.ConvertToBNL()
                    if ('frac' in thisFF and self.qparams.mom2max != [0] ):
                        self.F3div2M()
                    self.Write()
                    return self.CheckFF(thisFF,False)
                if ('frac' in thisFF and self.qparams.mom2max != [0] ) and resub:
                    self.F3div2M()
                    self.Write()
                    return self.CheckFF(thisFF,False)
                else:
                    print('WARNING, ' , thisFF , ' not found in self.FF_Stats last key')
                    print(list(self.FF_Stats['boot'][thiskey].fit_data['Params'].keys()))
                    print()
                    return False
            else:
                return True
        else:
            if thisFF not in list(self.FF_Stats['boot'][thiskey].fit_data['Params'].keys()):
                print('WARNING, ' , thisFF , ' not found in self.FF_Stats last key')
                print(list(self.FF_Stats['boot'][thiskey].fit_data['Params'].keys()))
                print()
                return False
            else:
                return True

    def CheckChi(self,this_chi,chi_list):
        chi_min,chi_max = chi_list.split('-')
        def RemoveChi(ichi):
            return int(ichi.replace('Chi',''))
        return RemoveChi(chi_min) <= RemoveChi(this_chi) <= RemoveChi(chi_max)

    def Plot(   self,plot_class,thisFF,thiscol='PreDefine',thisshift='PreDefine',
                thissym='PreDefine',plottflow='PreDefine',this_chi_list='Chi0-Chi20'):
        self.CheckCol(thiscol)
        self.CheckSym(thissym)
        self.thisshiftscale = 'Not Set'
        self.thisshift = self.GetShift([0,max(self.Qsqrdlist)],thisshift)
        if not self.CheckFF(thisFF):
            print(thisFF, 'Not found in ',self.name)
            return plot_class
        thisFFleg = thisFF.replace('tilde',r'\tilde')
        thisFFleg = r'$'+thisFFleg.replace('frac',r'\frac')+r'$'
        if plottflow == 'PreDefine': plottflow = self.tflows
        if isinstance(plottflow,(tuple,list,np.ndarray)):
            plottflow = plottflow[0]
        if isinstance(this_chi_list,(tuple,list,np.ndarray)):
            this_chi_list = this_chi_list[0]
        if isinstance(thisFF,(tuple,list,np.ndarray)):
            thisFF = thisFF[0]
        # for iQ,Avg,Std in zip(np.array(self.Qsqrdlist)+thisshift,self.ffresAvg.values(),self.ffresStd.values()):
        #     print abs(iQ),Avg[thisFF],Std[thisFF]
        vall,errl,ilist = [],[],[]
        medl,medupl,meddownl = [],[],[]
        plotQsqrd = []
        chi_min,chi_max = this_chi_list.split('-')
        if self.doflow:
            this_key = (plottflow,chi_min,slice(None),thisFF)
        else:
            this_key = (chi_min,slice(None),thisFF)
        for ikey,fit_data in self.FF_Stats['boot'].items():
            if not self.CheckChi(ikey[-2],this_chi_list):
                continue
            iQ2key = ikey[-1]
            for ipar,par_data in fit_data.fit_data['Params'].items():
                # if ipar not in thisFF and plot_the_ff: continue
                if iQ2key in self.QsqrdKeyList:
                    iQ2 = self.Qsqrdlist[self.QsqrdKeyList.index(iQ2key)]
                    plotQsqrd.append(iQ2)
                par_data.Stats()
                vall.append(par_data.Avg)
                errl.append(par_data.Std)
                medl.append(par_data.Median)
                medupl.append(par_data.MedStd_Up)
                meddownl.append(par_data.MedStd_Down)
                ilist.append(tuple(list(ikey)[:-1] + ['Qsqrd'+str(iQ2),ipar]) )
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=list(self.FF_Stats.index.names) + ['Form_Factors'])
            ploty = pa.Series(vall,index=indicies)
            plotyerr = pa.Series(errl,index=indicies)
            ploty_median = pa.Series(medl,index=indicies)
            plotyerr_up = pa.Series(medupl,index=indicies)
            plotyerr_down = pa.Series(meddownl,index=indicies)
            # if this_key not in ploty.index:
            #     print 'Warning, key', ', '.join(this_key) ,'not in FF_Stats'
        else:
            print('No form factor to plot for ',self.name)
            return plot_class

        hold_series = pa.Series()
        hold_series['boot_data'] = None
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = ploty
        hold_series['yerr_data'] = plotyerr
        hold_series['y_data_median'] = ploty_median
        hold_series['yerr_data_up'] = plotyerr_up
        hold_series['yerr_data_down'] = plotyerr_down
        hold_series['xerr_data'] = None
        hold_series['xerr_data_up'] = None
        hold_series['xerr_data_down'] = None
        hold_series['type'] = 'error_bar_vary'
        hold_series['key_select'] = this_key
        hold_series['fit_class'] = None
        hold_series['label'] = self.LegLab #' '.join([thisFFleg,,tflowleg,ichi])
        if hasattr(self,'doflowfull') and self.doflowfull:
            hold_series['label'] += r'$\alpha_{fit}$'
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = self.thisshift
        hold_series['xdatarange'] = None
        hold_series['scale'] = None
        hold_series['Phys'] = None
        hold_series['ShowPar'] = None
        hold_series['Median'] = False
        hold_series['fmt_class'] = KeyForamtting(self.ppparams)
        plot_class.AppendData(hold_series)
        return plot_class

    def ImportFitRanges(self,fit_ranges):
        if 'PreDef' == fit_ranges:
            self.fit_ranges = self.PredefFits
        else:
            self.fit_ranges = fit_ranges
        this_bool =  isinstance(self.fit_ranges,str) and self.fit_ranges == 'All'
        this_bool_2 = isinstance(self.fit_ranges,(np.ndarray,list,tuple)) and 'All' in list(self.fit_ranges)
        # if (this_bool or this_bool_2) and not self.Do_Chi_cutoff:
        if (this_bool or this_bool_2):
            self.fit_ranges = []
            for iQmin,iQmax in itertools.product(list(range(len(self.Qsqrdlist))),list(range(len(self.Qsqrdlist)))):
                if iQmax > iQmin:
                    self.fit_ranges.append('fitr'+str(iQmin)+'-'+str(iQmax))
        if isinstance(self.fit_ranges,str):
            self.fit_ranges = [self.fit_ranges]

    def SetFunction(self,thiscurrtype,thisDS):
        if '_'.join([thisDS,thiscurrtype]) in list(FFQ2FitFun.keys()):
            self.Fun = FFQ2FitFun['_'.join([thisDS,thiscurrtype])]
        else:
            self.Fun = FFQ2FitFun[thiscurrtype]


    ## wrapper for Fit to fit all form factors
    def FitFFs(self,thisFF='All',fit_ranges='All',iGuess='PreDef',
                    EstDir=False,WipeFit=False,show_timer=False):
        if thisFF == 'All':
            imax = np.argmax([len(list(ival.fit_data['Params'].keys())) for ival in self.FF_Stats['boot'].values])
            for ipar in list(self.FF_Stats['boot'].values[imax].fit_data['Params'].keys()):
                self.Fit(   ipar,fit_ranges=fit_ranges,iGuess=iGuess,
                            EstDir=EstDir,WipeFit=WipeFit,show_timer=show_timer)
        else:
            self.Fit(thisFF,fit_ranges=fit_ranges,iGuess=iGuess,EstDir=EstDir,WipeFit=WipeFit,show_timer=show_timer)

    def WipeAllFits(self):
        print('DEBUG, wiping all fits')
        self.FF_Fit_Stats = pa.DataFrame()

    def Stats(self):
        def stat_fun(val):
            val.Stats()
            return val
        if 'boot' in self.FF_Stats:
            self.FF_Stats.loc[:,'boot'] = self.FF_Stats.loc[:,'boot'].apply(stat_fun)
        if 'boot' in self.FF_Fit_Stats:
            self.FF_Fit_Stats.loc[:,'boot'] = self.FF_Fit_Stats.loc[:,'boot'].apply(stat_fun)
    ## fit_ranges are just related to the indexing of the momenta, (since we have already converted to 2pi/L etc..
    ## thisFFfun is a dictionary of touples { iform_factor ; (the_function, n_parameters) }
    def Fit(    self,thisFF,thisFFfun='PreDef',fit_ranges='All',iGuess='PreDef',
                EstDir=False,WipeFit=False,show_timer=True):
        if thisFFfun != 'PreDef': self.Fun[thisFF] = thisFFfun
        if self.Fun == 'None' or thisFF not in list(self.Fun.keys()):
            return False
            # print thisFF
            # print self.Fun
            # raise IOError('Please define funciton using .SetFunction(Fun) for fitting before calling .Fit()')
        if not self.CheckFF(thisFF): return
        self.ImportFitRanges(fit_ranges)
        # out_str = 'Debugging \n fit_ranges\n'
        # out_str += '\n'.join(fit_ranges)
        # out_str += 'self fit_ranges\n'
        # out_str += '\n'.join(self.fit_ranges)
        # out_str += '\n'.join(self.fit_ranges)
        # raise Exception(out_str)
        if iGuess != 'PreDef': self.iGuess = iGuess
        fboot,fAvg,fStd,ilist = [],[],[],[]
        # isAllChi = isinstance(chi_list,basestring) and chi_list == 'All'
        # if 'boot' in self.FF_Fit_Stats:
        #     print self.FF_Fit_Stats['boot']
        # raise EnvironmentError('DEBUGGING')
        for ifit in self.fit_ranges:
            if ifit == 'All':
                ifitmin,ifitmax = 0,len(self.Qsqrdlist)-1
                ifit = 'fitr'+str(ifitmin)+'-'+str(ifitmax)
            else:
                ifitmin,ifitmax = np.array(unxmlfitr(ifit))

            ifitmin,ifitmax = int(ifitmin),int(ifitmax)
            if self.doflow:
                this_Qkeylist = self.QsqrdKeyList[ifitmin:ifitmax+1]
                this_Qlist = self.Qsqrdlist[ifitmin:ifitmax+1]
                # if self.isCustomRFitR:
                #     if not hasattr(self,'chi_fit_list') or not isinstance(self.chi_fit_list,(list,tuple,np.ndarray,pa.Series,pa.DataFrame)):
                #         print self.name
                #         print self.chi_fit_list
                #         raise EnvironmentError('chi_fit_list is not present even though it is needed when isCustomRFitR is set to true')
                #     this_chi_list = self.chi_fit_list
                # else:
                for itflow in self.tflows:
                    this_chi_list = []
                    for iQ in this_Qkeylist:
                        if iQ not in list(self.chi_Qsqrd[itflow].keys()):
                            continue
                        this_chi_len = self.chi_Qsqrd[itflow][iQ]
                        if this_chi_len > self.nchi_fit_threshold:
                            this_rlist = np.random.choice(this_chi_len,size=self.nchi_fit_threshold,replace=False)
                            this_chi_list.append(['Chi'+str(ichi) for ichi in this_rlist])
                        else:
                            this_chi_list.append(['Chi'+str(ichi) for ichi in range(this_chi_len)])
                # self.chi_fit_list = pa.Series(this_chi_list,index=this_Qkeylist)
                    # this_chi_list.append(list(self.FF_Stats['boot'].index.get_loc_level((slice(None),iQ))[1]))
                    # this_chi_list.append(list(self.FF_Stats['boot'].index.get_loc_level((slice(None),iQ))[1]))
                if self.Do_Chi_cutoff and show_timer:
                    print()
                    print('number of fits for ',self.tflows[0],thisFF,ifit,' is ' ,len(list(itertools.product(*this_chi_list))))
                if show_timer: thistimer = Timer(  linklist=list(range(len(list(itertools.product(*this_chi_list))))),
                                    name='Fitting Form Factors '+', '.join([self.tflows[0],thisFF,ifit]))

                for ifit_chi,ichi_list in enumerate(itertools.product(*this_chi_list)):
                    chi_key = 'Chi'+str(ifit_chi)
                    fit_chi_key = (self.tflows[0],chi_key,thisFF,ifit)
                    if 'boot' in self.FF_Fit_Stats and fit_chi_key in self.FF_Fit_Stats['boot']:
                        if WipeFit:
                            self.FF_Fit_Stats.drop([fit_chi_key],inplace=True)
                        else:
                            if show_timer: thistimer.Lap(('fit already computed: ',chi_key))
                            continue
                    tdatarange,ydata = [[]],[]
                    for ichi,iQkey,iQ in zip(ichi_list,this_Qkeylist,this_Qlist):
                        # print 'DEBUG'
                        # print (ichi,iQkey)
                        # print self.FF_Stats['boot']
                        this_key = (self.tflows[0],ichi,iQkey)
                        key_str = ','.join(this_key)
                        if this_key not in self.FF_Stats.index:
                            # if show_timer: thistimer.Lap(('key not present ',thisFF,key_str))
                            continue
                        pardata = self.FF_Stats.loc[this_key,'boot']
                        # pardata = self.FF_Stats.loc[(ichi,iQkey),'boot']
                        if isinstance(pardata,pa.Series):
                            pardata = pardata.iloc[0]
                        if thisFF not in list(pardata.fit_data['Params'].keys()):
                            # if show_timer: thistimer.Lap(('key not present ',thisFF,key_str))
                            continue
                        tdatarange[0].append(float(iQ))
                        ydata.append(pardata.fit_data['Params'][thisFF])

                    if len(ydata) < self.Fun[thisFF][1]:
                        if show_timer: thistimer.Lap(('not enough npar ',thisFF,key_str))
                        # print 'ydata length ' ,len(ydata) , ', npar len' , self.Fun[thisFF][1]
                        # print 'WARNING, not enough fit parameters to fit function ' , self.Fun[thisFF][0].__name__ , ' to ' ,itflow,thisFF , ' with fit range ' , ifit
                        continue
                    if EstDir:
                        thisff = ff.Fitting(Funs=[self.Fun[thisFF][0],self.Fun[thisFF][1],'Estimate'],data=[np.array(tdatarange),np.array(ydata)],name=ifit+'$'+tflowTOLeg(itflow)+'$',iGuess=self.iGuess)
                    else:
                        thisff = ff.Fitting(Funs=self.Fun[thisFF],data=[np.array(tdatarange),np.array(ydata)],name=ifit+'$'+tflowTOLeg(self.tflows[0])+'$',iGuess=self.iGuess)
                    testff = thisff.FitBoots()
                    if np.isnan(testff):
                        if show_timer: thistimer.Lap(chi_key)
                        continue
                    ## self.FitDict [itflow , iFF , ifit]
                    ilist.append(fit_chi_key)
                    fboot.append(thisff)
                    fAvg.append(thisff.fit_data['ParamsAvg'])
                    fStd.append(thisff.fit_data['ParamsStd'])
                    if show_timer: thistimer.Lap(chi_key)
            else:
                this_Qkeylist = self.QsqrdKeyList[ifitmin:ifitmax+1]
                this_Qlist = self.Qsqrdlist[ifitmin:ifitmax+1]
                this_chi_list = []
                for iQ in this_Qkeylist:
                    if iQ not in list(self.chi_Qsqrd.keys()):
                        continue
                    this_chi_len = self.chi_Qsqrd[iQ]
                    if this_chi_len > self.nchi_fit_threshold:
                        this_rlist = np.random.choice(this_chi_len,size=self.nchi_fit_threshold,replace=False)
                        this_chi_list.append(['Chi'+str(ichi) for ichi in this_rlist])
                    else:
                        this_chi_list.append(['Chi'+str(ichi) for ichi in range(this_chi_len)])
                    # this_chi_list.append(list(self.FF_Stats['boot'].index.get_loc_level((slice(None),iQ))[1]))
                    # this_chi_list.append(list(self.FF_Stats['boot'].index.get_loc_level((slice(None),iQ))[1]))
                if self.Do_Chi_cutoff and show_timer:
                    print()
                    print('number of fits for ',thisFF,ifit,' is ' ,len(list(itertools.product(*this_chi_list))))
                if show_timer: thistimer = Timer(  linklist=list(range(len(list(itertools.product(*this_chi_list))))),
                                    name='Fitting Form Factors '+thisFF+', '+ifit)


                for ifit_chi,ichi_list in enumerate(itertools.product(*this_chi_list)):
                    chi_key = 'Chi'+str(ifit_chi)
                    fit_chi_key = (chi_key,thisFF,ifit)
                    if 'boot' in self.FF_Fit_Stats and fit_chi_key in self.FF_Fit_Stats['boot']:
                        if WipeFit:
                            self.FF_Fit_Stats.drop([fit_chi_key],inplace=True)
                        else:
                            if show_timer: thistimer.Lap(('fit already computed: ',chi_key))
                            continue
                    tdatarange,ydata = [[]],[]
                    for ichi,iQkey,iQ in zip(ichi_list,this_Qkeylist,this_Qlist):
                        # print 'DEBUG'
                        # print (ichi,iQkey)
                        # print self.FF_Stats['boot']
                        this_key = (ichi,iQkey)
                        key_str = ','.join(this_key)
                        if this_key not in self.FF_Stats.index:
                            # if show_timer: thistimer.Lap(('key not present ',thisFF,key_str))
                            continue
                        pardata = self.FF_Stats.loc[this_key,'boot']
                        # pardata = self.FF_Stats.loc[(ichi,iQkey),'boot']
                        if isinstance(pardata,pa.Series):
                            pardata = pardata.iloc[0]
                        if thisFF not in list(pardata.fit_data['Params'].keys()):
                            # if show_timer: thistimer.Lap(('key not present ',thisFF,key_str))
                            continue
                        tdatarange[0].append(float(iQ))
                        ydata.append(pardata.fit_data['Params'][thisFF])

                    # ##debugging
                    # print
                    # print 'FitData:'
                    # for it,iy in zip(tdatarange[0],ydata):
                    #     print it, iy
                    # print
                    if len(ydata) < self.Fun[thisFF][1] and ifit_chi == 0:
                        if show_timer: thistimer.Lap(('not enough npar ',thisFF,key_str))
                        # print
                        # print 'ydata length ' ,len(ydata) , ', npar len' , self.Fun[thisFF][1]
                        # print 'WARNING, not enough fit parameters to fit function ' , self.Fun[thisFF][0].__name__ , ' to ' ,thisFF , ' with fit range ' , ifit
                        # print
                        # thistimer.Lap(fit_chi_key)
                        continue
                    if EstDir:
                        thisff = ff.Fitting(Funs=[self.Fun[thisFF][0],self.Fun[thisFF][1],'Estimate'],data=[np.array(tdatarange),np.array(ydata)],name=ifit,iGuess=self.iGuess)
                    else:
                        thisff = ff.Fitting(Funs=self.Fun[thisFF],data=[np.array(tdatarange),np.array(ydata)],name=ifit,iGuess=self.iGuess)
                    testff = thisff.FitBoots()
                    if np.isnan(testff):
                        if show_timer: thistimer.Lap(chi_key)
                        continue
                    ## self.FitDict [iff , ifitr ] fitting instance
                    ilist.append(fit_chi_key)
                    fboot.append(thisff)
                    fAvg.append(thisff.fit_data['ParamsAvg'])
                    fStd.append(thisff.fit_data['ParamsStd'])
                    if show_timer: thistimer.Lap(chi_key)
        # debug_str = ' '.join(map(str,[self.fit_ranges,WipeFit]))
        # raise Exception('Debugging, '+debug_str)
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.FF_Fit_Col_Names)
            if 'boot' in self.FF_Fit_Stats.columns:
                this_df = pa.DataFrame()
                this_df['boot'] = pa.Series(fboot,index=indicies)
                this_df['Avg'] = pa.Series(fAvg,index=indicies)
                this_df['Std'] = pa.Series(fStd,index=indicies)
                self.FF_Fit_Stats = self.FF_Fit_Stats.append(this_df)
            else:
                self.FF_Fit_Stats['boot'] = pa.Series(fboot,index=indicies)
                self.FF_Fit_Stats['Avg'] = pa.Series(fAvg,index=indicies)
                self.FF_Fit_Stats['Std'] = pa.Series(fStd,index=indicies)
            self.Write()
            return True
        else:
            # print 'No Fits were done.'
            return False

    def ValsAndFitPlot( self,plot_class,thisFF,fitr='All',xlims=[0,'ToDataMax'],
                        thiscol='PreDefine',thisshift='PreDefine',thissym='PreDefine',
                        plottflow='PreDefine',this_chi_list='Chi0-Chi20'):
        if plottflow == 'First': plottflow = [self.tflows[0]]

        plot_class = self.Plot( plot_class,thisFF,this_chi_list=this_chi_list,thiscol=thiscol,thisshift=thisshift,
                                thissym=thissym,plottflow=plottflow)
        if plottflow == 'PreDefine': plottflow = self.tflows
        if self.doflow or self.doflowfull:
            for itflow in plottflow:
                plot_class = self.FitPlot(  plot_class,thisFF,this_chi_list=this_chi_list,
                                            fitr=fitr,tflow=itflow,xlims=xlims,thiscol=thiscol,thisshift=thisshift)
        else:
            plot_class = self.FitPlot(  plot_class,thisFF,this_chi_list=this_chi_list,
                                        fitr=fitr,xlims=xlims,thiscol=thiscol,thisshift=thisshift)
        return plot_class

    ## state is 1 or 2 state fit
    ## fitr is fit range e.g. fitr5-10
    ## fit range of 'Data' is just the original fit range
    def FitPlot(self,plot_class,thisFF,fitr='All',tflow='First',
                xlims = [0,'ToDataMax'],this_chi_list='Chi0-Chi20',thiscol='PreDefine',thisshift='PreDefine'):
        self.CheckCol(thiscol)
        check = True
        # thisshift = self.GetShift(pl.xlim(),thisshift)
        thisshift = self.thisshift
        if isinstance(tflow,(list,tuple,np.ndarray)): tflow = tflow[0]
        if 'boot' not in self.FF_Fit_Stats.columns:
            print('First Fit not present, please run Fit before calling fitplot')
            return plot_class
            # raise EnvironmentError('First Fit not present, please run Fit')
        if fitr == 'All':
            ifitmin,ifitmax = 0,len(self.Qsqrdlist)-1
            fitr = 'fitr'+str(ifitmin)+'-'+str(ifitmax)
        if isinstance(this_chi_list,(list,tuple,np.ndarray)):
            this_chi_list = this_chi_list[0]
        chi_min,chi_max = this_chi_list.split('-')
        if self.doflow:
            if tflow == 'First':
                 tflow = self.tflows[0]
            this_key = (tflow,chi_min,thisFF,fitr)
        else:
            this_key = (chi_min,thisFF,fitr)
        # print 'DEBUG',this_key, this_key not in self.FF_Fit_Stats['boot']
        # print self.FF_Fit_Stats['boot'].index

        if this_key not in self.FF_Fit_Stats['boot']:
            print(', '.join(this_key) , ' has not been done, performing fit now')
            check = self.Fit(thisFF,fit_ranges=[fitr])
        if not check:
            raise AttributeError('Warning, error with fitting ',', '.join(this_key))
        vall,indexl = [],[]
        for ikey,ival in self.FF_Fit_Stats['boot'].items():
            if self.CheckChi(ikey[-3],this_chi_list):
                vall.append(ival)
                indexl.append(ikey)
        if len(indexl) > 0:
            indicies = pa.MultiIndex.from_tuples(indexl,names=self.FF_Fit_Stats['boot'].index.names)
            fitr_params = pa.Series(vall,index=indicies)
        else:
            print('Warning, no fits found that satisfy',this_chi_list)
            return plot_class
        if self.doflow:
            this_key = (tflow,chi_min,thisFF,fitr)
        else:
            this_key = (chi_min,thisFF,fitr)
        # if isinstance(fitr_params,pa.Series): fitr_params = fitr_params.iloc[-1]
        ## TODO, plot from zero to data max
        hold_series = pa.Series()
        hold_series['x_data'] = None
        hold_series['y_data'] = None
        hold_series['yerr_data'] = None
        hold_series['xerr_data'] = None
        hold_series['boot_data'] = None
        hold_series['key_select'] = this_key
        hold_series['type'] = 'fit_vary'
        hold_series['fit_class'] = fitr_params
        hold_series['label'] = 'fit '
        hold_series['symbol'] = self.thissym
        hold_series['color'] = self.thiscol
        hold_series['shift'] = thisshift
        hold_series['xdatarange'] = xlims
        hold_series['scale'] = None
        hold_series['ShowPar'] = 'Par2'
        hold_series['fmt_class'] = KeyForamtting(self.ppparams)
        plot_class.AppendData(hold_series)
        return plot_class

    def Get_Extrapolation(self):
        outl,indexl = [],[]
        if 'boot' not in self.FF_Fit_Stats.columns:
            print('First Fit not present, please run Fit')
            return pa.Series()
        bool_index = np.array(self.FF_Fit_Col_Names)!='Chi2pdf_number'
        col_names = np.array(self.FF_Fit_Col_Names)[bool_index].tolist()
        for ikey,chi_data in self.FF_Fit_Stats['boot'].groupby(col_names):
            for fit_key,fit_index in zip(['$Q^{2}F^{\'}(0)$','$F(0)$'],[0,1]):
                if self.doflow:
                    indexl.append(chi_data.index[0][0:1]+chi_data.index[0][2:]+(fit_key,))
                else:
                    indexl.append(chi_data.index[0][1:]+(fit_key,))
                outl.append([])
                for ic,(ichi,ichi_data) in enumerate(chi_data.items()):
                    if ic > self.nchi_fit_threshold**(len(self.q2list)): continue
                    plot_chi_data = ichi_data
                    if isinstance(plot_chi_data,pa.Series):
                        ## assumes first parameter is constant part of linear function
                        plot_chi_data = chi_data.iloc[0]
                    outl[-1] += list(plot_chi_data.fit_data['Params'].iloc[fit_index].bootvals.values)
                outl[-1] = BootStrap(len(outl[-1]),bootvals=outl[-1],name=', '.join(ikey))
        if len(indexl) > 0:
            indicies = pa.MultiIndex.from_tuples(indexl,names=col_names + ['fit parameter'])
            return pa.Series(outl,index=indicies).dropna()
        else:
            return pa.Series()

    def PlotFitDistribution(self,plot_class,thisFF,fitr='All',tflow='First',thiscol='PreDefine',chi_colapse=True):
        self.CheckCol(thiscol)
        # thisshift = self.GetShift(pl.xlim(),thisshift)
        if isinstance(tflow,(list,tuple,np.ndarray)): tflow = tflow[0]
        if 'boot' not in self.FF_Fit_Stats.columns:
            print('First Fit not present, please run Fit')
            return plot_class
            # raise EnvironmentError()
        check = True
        if fitr == 'All':
            ifitmin,ifitmax = 0,len(self.Qsqrdlist)-1
            fitr = 'fitr'+str(ifitmin)+'-'+str(ifitmax)
        if self.doflow:
            if tflow == 'First':
                 tflow = self.tflows[0]
            if (tflow,fitr,thisFF) not in self.FF_Fit_Stats['boot'].swaplevel(1,3):
                print(tflow,thisFF,fitr , ' has not been done, performing fit now')
                check = self.Fit(thisFF,fit_ranges=[fitr])
            if check:
                ## From Here!!
                # this_key = (tflow,thisFF,fitr)
                # fitr_label = self.FF_Fit_Stats['boot'][(tflow,'Chi0',thisFF,fitr)].name
                if chi_colapse:
                    this_key = (tflow,thisFF,fitr,'$F(0)$')
                    plot_data = self.Get_Extrapolation()
                    if len(plot_data) == 0:
                        check = False
                else:
                    ## TODO, make this like above with an extra index for F(0) or Q2F'
                    this_key = (tflow,'Chi0',thisFF,fitr)
                    plot_data = self.FF_Fit_Stats['boot'].apply(lambda x : x.fit_data['Params'].iloc[-1])
                ## TODO, plot from zero to data max
        else:
            if (fitr,thisFF) not in self.FF_Fit_Stats['boot'].swaplevel(0,2):
                print(thisFF,fitr , ' has not been done, performing fit now')
                check = self.Fit(thisFF,fit_ranges=[fitr])
            if check:
                if chi_colapse:
                    this_key = (thisFF,fitr,'$F(0)$')
                    plot_data = self.Get_Extrapolation()
                    if len(plot_data) == 0:
                        check = False
                else:
                    this_key = ('Chi0',thisFF,fitr)
                    plot_data = self.FF_Fit_Stats['boot'].apply(lambda x : x.fit_data['Params'].iloc[-1])
        if check:
            hold_series = pa.Series()
            hold_series['boot_data'] = plot_data
            hold_series['x_data'] = None
            hold_series['y_data'] = None
            hold_series['yerr_data'] = None
            hold_series['xerr_data'] = None
            hold_series['type'] = 'histogram_vary'
            hold_series['fit_class'] = None
            hold_series['key_select'] = this_key
            hold_series['label'] = r'fit_result'
            hold_series['symbol'] = None
            hold_series['color'] = self.thiscol
            hold_series['shift'] = None
            hold_series['xdatarange'] = None
            hold_series['scale'] = None
            hold_series['ShowPar'] = None
            hold_series['fmt_class'] = KeyForamtting(self.ppparams)
            plot_class.AppendData(hold_series)
        else:
            print('Fit was unable to be constructed, no fit is plotted.')
        return plot_class

    ### OVERLOADING OPERATOR FUNCTIONS ###
    def Overload_wrap(self,FF2,this_fun):
        ## TODO, include some name and filename formatting
        thisname,thisfilename = self.name,self.filename
        if hasattr(FF2,'name'):
            thisname += '_'+this_fun.__name__+'_'+FF2.name
        if hasattr(FF2,'filename'):
            thisfilename += '_'+this_fun.__name__+'_'+FF2.filename
        result = FFEquation(  self.currtype,q2list=self.q2list,pp2list=self.pp2list,
                                cfglist=self.cfglistin, Info=self.Info,
                                fit2pt=self.fit2ptin,fitAlpha=self.fitAlphain,name=self.namein,
                                thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
        result.FF_Stats = pa.DataFrame(columns = ['boot'])
        result.Qsqrdlist = self.Qsqrdlist
        result.QsqrdKeyList = self.QsqrdKeyList
        result.chi_Qsqrd = self.chi_Qsqrd
        result.C2forMass = self.C2forMass
        if self.doflow:
            result.C2FOforAlpha = self.C2FOforAlpha
        result.thisFF = self.thisFF
        result.RatioFileList = self.RatioFileList
        result.Inefm = self.Inefm
        result.nchi_fit_threshold = self.nchi_fit_threshold
        result.nchi_threshold = self.nchi_threshold

        # result.UpdateName('+',self,FF2)
        if 'boot' in self.FF_Stats:
            if isinstance(FF2,FFEquation):
                result.RatioFileList += FF2.RatioFileList
                if self.DS != FF2.DS:
                    result.DS = '_'.join([self.DS, this_fun.__name__,FF2.DS])
                if self.Inefm != FF2.Inefm:
                    print('WARNING:, combining two form factors, where one is in [efm]')
                result.SetCustomName()

                # result.FF_Stats.loc[:,'boot'] = this_fun(self.FF_Stats.loc[:,'boot'],FF2.FF_Stats.loc[:,'boot'])
                if 'boot' in FF2.FF_Stats:
                    vall,ilist = [],[]
                    if self.doflow and FF2.doflow:
                        for (itflow,ichi,iQ2),iFF in self.FF_Stats['boot'].items():
                            if (itflow,iQ2) not in FF2.FF_Stats['boot'].index.swaplevel(1,2):
                                if ichi.lower() == 'chi0':
                                    print(itflow,iQ2,'not key found in ',FF2.name , 'skipping')
                                    # print 'DEBUG'
                                    # print FF2.FF_Stats['boot'].index.swaplevel(1,2)
                                continue
                                # out_str = ' missmatch in euclidian times / flow time / op times  '
                                # out_str += itflow+' '+iQ2+'\n'
                                #
                                # out_str += ' when combining form factors \n'
                                # out_str += str(self) + '\n'
                                # out_str += this_fun.__name__+'\n'
                                # out_str += str(FF2) + '\n'
                                # raise EnvironmentError(out_str)
                            ichi_len = len(FF2.FF_Stats.loc[(itflow,slice(None),iQ2),'boot'])
                            inside_chi = int(ichi.replace('Chi',''))*ichi_len

                            for (itflow_dump,ichi2,iQ2_dump),iFF2 in FF2.FF_Stats.loc[(itflow,slice(None),iQ2),'boot'].items():
                                this_chi = int(ichi2.replace('Chi','')) + inside_chi
                                ilist.append((itflow,'Chi'+str(this_chi),iQ2))
                                vall.append(this_fun(iFF,iFF2))
                    elif self.doflow:
                        for (iflow,ichi,iQ2),iFF in self.FF_Stats['boot'].items():
                            ichi_len = len(FF2.FF_Stats.loc[(itflow,slice(None),iQ2),'boot'])
                            inside_chi = int(ichi.replace('Chi',''))*ichi_len
                            if (iQ2,) not in FF2.FF_Stats['boot'].index.swaplevel(0,1):
                                if ichi.lower() == 'chi0':
                                    print(iQ2,'not found in ',FF2.name,'skipping')
                                continue
                                # out_str = ' missmatch in euclidian times / flow time / op times \
                                #             when combining form factors \n'
                                # out_str += str(self) + '\n'
                                # out_str += this_fun.__name__+'\n'
                                # out_str += str(FF2) + '\n'
                                # raise EnvironmentError(out_str)
                            for (ichi2,iQ2_dump),iFF2 in FF2.FF_Stats.loc[(slice(None),iQ2),'boot'].items():
                                this_chi = int(ichi2.replace('Chi','')) + inside_chi
                                ilist.append((itflow,'Chi'+str(this_chi),iQ2))
                                vall.append(this_fun(iFF,iFF2))
                    elif FF2.doflow:
                        for (ichi,iQ2),iFF in self.FF_Stats['boot'].items():
                            ichi_len = len(FF2.FF_Stats.loc[(slice(None),iQ2),'boot'])
                            inside_chi = int(ichi.replace('Chi',''))*ichi_len
                            if (iQ2,) not in FF2.FF_Stats['boot'].index.swaplevel(0,2):
                                if ichi.lower() == 'chi0':
                                    print(iQ2,'not found in ',FF2.name,'skipping')
                                continue
                                # out_str = ' missmatch in euclidian times / flow time / op times \
                                #             when combining form factors \n'
                                # out_str += str(self) + '\n'
                                # out_str += this_fun.__name__+'\n'
                                # out_str += str(FF2) + '\n'
                                # raise EnvironmentError(out_str)
                            for (itflow,ichi2,iQ2_dump),iFF2 in FF2.FF_Stats.loc[(slice(None),slice(None),iQ2),'boot'].items():
                                this_chi = int(ichi2.replace('Chi','')) + inside_chi
                                ilist.append((itflow,'Chi'+str(this_chi),iQ2))
                                vall.append(this_fun(iFF,iFF2))
                    else:
                        for (ichi,iQ2),iFF in self.FF_Stats['boot'].items():
                            ichi_len = len(FF2.FF_Stats.loc[(slice(None),iQ2),'boot'])
                            inside_chi = int(ichi.replace('Chi',''))*ichi_len
                            if (iQ2,) not in FF2.FF_Stats['boot'].index.swaplevel(0,1):
                                if ichi.lower() == 'chi0':
                                    print(iQ2,'not found in ',FF2.name,'skipping')
                                continue
                                # out_str = ' missmatch in euclidian times / flow time / op times \
                                #             when combining form factors \n'
                                # out_str += str(self) + '\n'
                                # out_str += this_fun.__name__+'\n'
                                # out_str += str(FF2) + '\n'
                                # raise EnvironmentError(out_str)
                            for (ichi2,iQ2_dump),iFF2 in FF2.FF_Stats.loc[(slice(None),iQ2),'boot'].items():
                                this_chi = int(ichi2.replace('Chi','')) + inside_chi
                                ilist.append(('Chi'+str(this_chi),iQ2))
                                vall.append(this_fun(iFF,iFF2))
                    if len(ilist) > 0:
                        indicies = pa.MultiIndex.from_tuples(ilist,names=self.FF_Col_Names)
                        result.FF_Stats.loc[:,'boot'] = pa.Series(vall,index=indicies)
                else:
                    raise EnvironmentError('Run Bootstrap on second class before combining ratio functions')
            else:
                try:
                    result.FF_Stats.loc[:,'boot'] = this_fun(self.FF_Stats.loc[:,'boot'], FF2)
                except Exception as err:
                    print(type(FF2), FF2)
                    print(self.FF_Stats)
                    raise EnvironmentError(str(err)+'\nInvalid value to combine with FF_Stats class')
            return result
        else:
            raise EnvironmentError('Run Bootstrap on first class before combining FF_Stats functions')

    def __add__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in FFEquation.parentlist]):
            return FF2.__radd__(self)
        else:
            return self.Overload_wrap(FF2,op.add)

    def __sub__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in FFEquation.parentlist]):
            return FF2.__rsub__(self)
        else:
            return self.Overload_wrap(FF2,op.sub)

    def __mul__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in FFEquation.parentlist]):
            return FF2.__rmul__(self)
        else:
            return self.Overload_wrap(FF2,op.mul)

    def __div__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in FFEquation.parentlist]):
            return FF2.__rdiv__(self)
        else:
            return self.Overload_wrap(FF2,op.truediv)

    def __pow__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in FFEquation.parentlist]):
            return FF2.__rpow__(self)
        else:
            return self.Overload_wrap(FF2,op.pow)


    ## Right multiplication functions


    def __radd__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in FFEquation.parentlist]):
            return FF2.__add__(self)
        else:
            return self.Overload_wrap(FF2,flip_fun_2arg(op.add))

    def __rsub__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in FFEquation.parentlist]):
            return FF2.__sub__(self)
        else:
            return self.Overload_wrap(FF2,flip_fun_2arg(op.sub))

    def __rmul__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in FFEquation.parentlist]):
            return FF2.__mul__(self)
        else:
            return self.Overload_wrap(FF2,flip_fun_2arg(op.mul))

    def __rdiv__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in FFEquation.parentlist]):
            return FF2.__div__(self)
        else:
            return self.Overload_wrap(FF2,flip_fun_2arg(op.truediv))


    def __rpow__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in FFEquation.parentlist]):
            return FF2.__pow__(self)
        else:
            return self.Overload_wrap(FF2,flip_fun_2arg(op.pow))

    def __str__(self):
        return self.name

    def __repr__(self):
        return str(self)


def TestFFSolve(DefWipe=False,q2list=[1,4],DS='Proton'):
    from Params import defInfo
    defInfo['Chi_Cutoff'] = True
    defInfo['Chi2pdf_Min'] = 0.5
    defInfo['Chi2pdf_Max'] = 1.0
    defInfo['Important_FF'] = False
    defInfo['Reduce_System'] = True
    defInfo['DoubSing'] = DS

    this_info = pa.Series()
    this_info['save_file'] = './TestGraphs/testFF.pdf'
    this_info['title'] = 'Test Form Factor'

    data = FFEquation('GeGm',Info=defInfo,q2list=q2list)
    data.SetupPickle(DefWipe=DefWipe)
    data.SolvePickle(DefWipe=DefWipe)
    data.LoadPickle(DefWipe=DefWipe)
    import PlotData as jpl
    data_plot = jpl.Plotting(plot_info=this_info)
    data_plot = data.Plot(data_plot,'G_{E}',thiscol='Blue',thisshift=0.0,thissym='o')
    data_plot = data.FitPlot(data_plot,'G_{E}')
    data_plot.PrintData()
    data_plot.PlotAll()
    return data

def TestFFFlowSolve(DefWipe=False,q2list=[1,4]):
    from Params import defInfoFlow
    defInfoFlow['Chi_Cutoff'] = True
    defInfoFlow['Chi2pdf_Min'] = 0.5
    defInfoFlow['Chi2pdf_Max'] = 1.0
    defInfoFlow['Important_FF'] = False
    defInfoFlow['Reduce_System'] = True
    defInfoFlow['tflowlist'] = np.array(['t_f6.01'])
    defInfoFlow['tflowfit'] = defInfoFlow['tflowlist']
    this_info = pa.Series()
    this_info['save_file'] = './TestGraphs/testFF.pdf'
    this_info['title'] = 'Test Form Factor'
    data = FFEquation('VectorTop',Info=defInfoFlow,q2list=q2list)
    data.SetupPickle(DefWipe=DefWipe)
    data.SolvePickle(DefWipe=DefWipe)
    data.LoadPickle(DefWipe=DefWipe)
    import PlotData as jpl
    data_plot = jpl.Plotting(plot_info=this_info)
    data_plot = data.Plot(data_plot,'G_{E}',thiscol='Blue',thisshift=0.0,thissym='o')
    data_plot = data.FitPlot(data_plot,'G_{E}')
    data_plot.PrintData()
    data_plot.PlotAll()
    return data

if __name__ == '__main__':
    dataFF = TestFFSolve()
