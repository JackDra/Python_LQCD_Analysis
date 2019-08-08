#!/usr/bin/env python

import numpy as np

# import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
# import matplotlib.pyplot as pl
import PlotData as jpl
from NullPlotData import null_series
from copy import copy

# from collections import OrderedDict
import pandas as pa

# from MiscFuns import NumbToFileName,mkdir_p,ODNested
from MiscFuns import mkdir_p,op_dict
from Params import graphdir,defMcore,nboot
from TimeStuff import Timer
from XmlFormatting import tflowTOLeg
from PredefFitFuns import F3ChiralLimFun,SquaredFitFun,LinearFitFun,LinearFF_2D,F3LatSpaceFun,SchiffFitFun
import FitFunctions as ff
from MomParams import physcial_pion_mass
from FileIO import Construct_Empty

from FFSolve import FFEquation

## Move Somewhere?
# params = {'legend.fontsize': 7,
#           'legend.numpoints': 1,
#           'axes.labelsize' : 10,
#           'figure.autolayout': True,
#           'axes.grid': True,
#           'errorbar.capsize': 5,
#           'axes.xmargin':0.01,
#           'axes.titlesize' : 20,
#           'axes.ymargin':0.01}
# pl.style.use(plot_style)
# pl.rcParams.update(params)


colourset8 = [ '#000080','#B22222','#00B800','#8B008B', '#0000FF','#FF0000','#000000','#32CD32','#FF0066']
markerset = ['o','s','^','v','d']
shiftmax = 3
shiftper = 0.005 ##R function shift
shiftset = [0]
for ish in np.arange(1,shiftmax+1): shiftset += [-ish*shiftper,ish*shiftper]
del ish, shiftmax,shiftper

def RemoveFlowKey(this_key):
    out_key = []
    for ikey in this_key:
        if 't_f' not in ikey:
            out_key.append(ikey)
    return out_key

class SetOfFF(object):
    """

    FF(set) Inherits from FFSolve Class

    plotting parameters (Move to different file?) ####

    InfoDict is a list of Info dictionaries to be bassed to FFSolve:
    Setup to do fits taken from dictionaries within InfoDict (see FFSlove.py)
    Although Info['gammas'] can have multiple gamma matricies,
    This code is only set up to take a list, where each list only has one gamma matrix in Info['gammas'] ( e.g. = ['g4'])




    """

    def Construct_FF_Set(this_file_list):
        set_list = [FFEquation.Construct_FF_File(ifile) for ifile in this_file_list]
        out_class = Construct_Empty(SetOfFF)
        out_class.cfglist = {}
        out_class.cfglist2pt = {}
        out_class.InfoDict = []
        out_class.Plotdir = graphdir+'/FF/'
        out_class.name = 'Set1'
        out_class.SetFF = pa.DataFrame()
        out_class.SetFF.loc[:,'FF_class'] = pa.Series(set_list,index=[ilist.namein for ilist in set_list])
        return out_class

    instances = []

    parentlist = []

    def __init__(self, cfglist={}, cfglist2pt={},InfoDict=[],name = '',setnames = []):

        self.cfglist = cfglist
        self.cfglist2pt = cfglist2pt
        self.InfoDict = InfoDict
        self.Plotdir = graphdir+'/FF/'
        mkdir_p(self.Plotdir)


        if len(name) == 0:
            self.name = 'Set'+str(int(len(SetOfFF.instances))+1)
            SetOfFF.instances.append(self.name)
        else:
            self.name = name
            # SetOfFF.instances.append(self.name)

        ## dictionary containing FormFactCorrs,capsize=5
        self.SetFF = pa.DataFrame()
        coll,syml,shiftl = [],[],[]
        Infol,q2l,pp2l = [],[],[]
        fit2ptl,fitalphal = [],[]
        # fitratl = []
        cfgl,cfgl2pt = [],[]
        FFl,ilist = [],[]
        fal = []
        if len(setnames) == 0:
            setnames = ['']*len(InfoDict)
        for ic,(iname,idict) in enumerate(zip(setnames,self.InfoDict)):
            Infol.append(idict)
            coll.append(colourset8[ic%len(colourset8)])
            syml.append(markerset[ic%len(markerset)])
            shiftl.append(shiftset[ic%len(shiftset)])
            if 'CurrType' not in list(idict.keys()):
                raise IOError('CurrType needs to be in Info Dictionaries')

            # if 'Full_Rat' in idict.keys():
            #     frat.append(idict['Full_Rat'])
            #     if frat[-1]:
            #         if 'FitsRFFull' in idict.keys():
            #             fitratl.append(idict['FitsRFFull'][0])
            #         else:
            #             fitratl.append('From Def')
            #     else:
            #         if 'FitsRF' in idict.keys():
            #             fitratl.append(idict['FitsRF'][0])
            #         else:
            #             fitratl.append('From Def')
            # else:
            #     frat.append(False)
            #     if 'FitsRF' in idict.keys():
            #         fitratl.append(idict['FitsRF'][0])
            #     else:
            #         fitratl.append('From Def')

            if 'Alpha_Full' in list(idict.keys()):
                fal.append(idict['Alpha_Full'])
                if fal[-1]:
                    if 'FitsAlphaFull' in list(idict.keys()):
                        fitalphal.append(idict['FitsAlphaFull'][0])
                    else:
                        fitalphal.append('From Def')
                else:
                    if 'FitsAlpha' in list(idict.keys()):
                        fitalphal.append(idict['FitsAlpha'][0])
                    else:
                        fitalphal.append('From Def')
            else:
                fal.append(False)
                if 'FitsAlpha' in list(idict.keys()):
                    fitalphal.append(idict['FitsAlpha'][0])
                else:
                    fitalphal.append('From Def')
            if 'q2list' not in list(idict.keys()):
                q2l.append([0,1,2,3,4])
            else:
                q2l.append(idict['q2list'])
            if 'pp2' not in list(idict.keys()):
                pp2l.append(0)
            else:
                pp2l.append(idict['pp2'])
            if 'Fits' in list(idict.keys()):
                fit2ptl.append(idict['Fits'][0])
            else:
                fit2ptl.append('From Def')

            if isinstance(self.cfglist,str):
                ## kappa up down (in integer form)
                if 'kud' in list(idict.keys()): thiskud = 'kud'+str(idict['kud'])
                else: thiskud = ''
                this_nxyz = idict['nxyzt'][0]
                this_nt = idict['nxyzt'][-1]
                this_dim = 'RC'+str(this_nxyz)+'x'+str(this_nt)
                ## kappa strange (in integer form)
                if 'ks' in list(idict.keys()): thisks = 'ks'+str(idict['ks'])
                else: thisks = ''
                kappafolder = this_dim+'Kud0'+thiskud.replace('kud','')+'Ks0'+thisks.replace('ks','')
                fotype = ''
                if 'Top' in idict['CurrType']:
                    fotype = '_TopCharge'
                elif 'Wein' in idict['CurrType']:
                    fotype = '_Weinberg'
                # print 'cfglist',str(idict['outdir']+kappafolder+'/cfglist'+fotype+'.py3p')
                cfgl.append(str(idict['outdir']+kappafolder+'/cfglist'+fotype+'.py3p'))
            else:
                cfgl.append(self.cfglist)

            if isinstance(self.cfglist2pt,str):
                ## kappa up down (in integer form)
                if 'kud' in list(idict.keys()): thiskud = 'kud'+str(idict['kud'])
                else: thiskud = ''

                ## kappa strange (in integer form)
                if 'ks' in list(idict.keys()): thisks = 'ks'+str(idict['ks'])
                else: thisks = ''
                kappafolder = this_dim+'Kud0'+thiskud.replace('kud','')+'Ks0'+thisks.replace('ks','')
                fotype = ''
                if 'Top' in idict['CurrType']:
                    fotype = '_TopCharge'
                elif 'Wein' in idict['CurrType']:
                    fotype = '_Weinberg'
                # print 'cfglist2pt',str(idict['outdir']+kappafolder+'/cfglist2pt'+fotype+'.py3p')
                cfgl2pt.append(str(idict['outdir']+kappafolder+'/cfglist2pt'+fotype+'.py3p'))
            else:
                cfgl2pt.append(self.cfglist2pt)

            FFl.append(FFEquation(  idict['CurrType'],q2list=q2l[-1],pp2list=[pp2l[-1]],
                                    fit2pt=fit2ptl[-1],fitAlpha=fitalphal[-1],
                                    cfglist=cfgl[-1],cfglist2pt=cfgl2pt[-1], Info=idict,
                                    thissym=syml[-1],thiscol=coll[-1],thisshift=shiftl[-1]))
            # FFl.append(FFEquation(  idict['CurrType'],q2list=q2l[-1],pp2list=[pp2l[-1]],
            #                         fit2pt=fit2ptl[-1],fitAlpha=fitalphal[-1],fitratio=fitratl[-1],
            #                         cfglist=cfgl[-1],cfglist2pt=cfgl2pt[-1], Info=idict,
            #                         thissym=syml[-1],thiscol=coll[-1],thisshift=shiftl[-1]))
            if str(FFl[-1]) in ilist:
                FFnext = 2
                while str(FFl[-1])+'_'+str(FFnext) in ilist:
                    FFnext += 1
                ilist.append(str(FFl[-1])+'_'+str(FFnext))
            else:
                ilist.append(str(FFl[-1]))
            # print 'Setting up classes for ' , thisFF , ' complete'
        # print
            # print 'has color',thiscol , 'symbol ' , thissym , 'thisshift',thisshift
        if len(ilist) > 0:
            self.SetFF.loc[:,'FF_class'] = pa.Series(FFl,index=ilist)
            self.SetFF.loc[:,'color'] = pa.Series(coll,index=ilist)
            self.SetFF.loc[:,'symbol'] = pa.Series(syml,index=ilist)
            self.SetFF.loc[:,'shift'] = pa.Series(shiftl,index=ilist)
            self.SetFF.loc[:,'info_dictionarie'] = pa.Series(Infol,index=ilist)
            self.SetFF.loc[:,'q2_list'] = pa.Series(q2l,index=ilist)
            self.SetFF.loc[:,'pp2_list'] = pa.Series(pp2l,index=ilist)
            self.SetFF.loc[:,'2pt_fit_range'] = pa.Series(fit2ptl,index=ilist)
            self.SetFF.loc[:,'alpha_fit_range'] = pa.Series(fitalphal,index=ilist)
            # self.SetFF.loc[:,'ratio_fit_range'] = pa.Series(fitratl,index=ilist)
            self.SetFF.loc[:,'config_list'] = pa.Series(cfgl,index=ilist)
            self.SetFF.loc[:,'config_list_2pt'] = pa.Series(cfgl2pt,index=ilist)

        self.SetFFlen = len(self.SetFF)
        ##Implement Reading cfglist if len(cfglist) == 0

        self.current = 0


    def SetPlotParams(self):
        for ic,(iset,set_series) in enumerate(self.SetFF.iterrows()):
            thiscol = colourset8[ic%len(colourset8)]
            thissym = markerset[ic%len(markerset)]
            thisshift = shiftset[ic%len(shiftset)]
            if thiscol != set_series['color']:
                print('Color has changed for ',iset)
            if thissym != set_series['symbol']:
                print('symbol has changed for ',iset)
            if thisshift != set_series['shift']:
                print('shift has changed for ',iset)
            set_series['FF_class'].ImportPlotParams(thiscol,thissym,thisshift)

    def Write(self,WipeData=False):
        for setdata in self.SetFF['FF_class'].values:
            setdata.Write(WipeData=WipeData)


    def SetupPickle(self,DefWipe=False,HardWipe=True,Mcore=defMcore,CommonRatios=False):
        # doubling this should fix some memory issues with the calculation.
        FFtimer = Timer(linklist=list(self.SetFF['FF_class'].keys()),name='Form Factor Read '+self.name,rewrite=False)
        RatioFitList = [False,False]
        for setdata in self.SetFF['FF_class'].values:
            Rfit_test,Rfit_test_flow = setdata.SetupPickle( DefWipe=DefWipe,HardWipe=HardWipe,
                                                            Mcore=Mcore,RatioFitLists=RatioFitList)
            if CommonRatios:
                if Rfit_test is False:
                    print('Warning, rfitindicies not found, must be from old code')
                    setdata.SetupPickle(DefWipe=True,HardWipe=HardWipe,Mcore=Mcore,
                                        RatioFitLists=RatioFitList)
                RatioFitList = Rfit_test,Rfit_test_flow
            FFtimer.Lap()
        # for setdata in self.SetFF['FF_class'].values:
        #     setdata.LoadPickle(DefWipe=False,HardWipe=HardWipe,Mcore=False)
        # self.SetPlotParams()
        # self.FitFFs()

    def SolvePickle(self,DefWipe=False,HardWipe=True,Mcore=defMcore,CommonRatios=False):
        # doubling this should fix some memory issues with the calculation.
        FFtimer = Timer(linklist=list(self.SetFF['FF_class'].keys()),name='Form Factor Solving '+self.name,rewrite=False)
        RatioFitList = [False,False]
        for setdata in self.SetFF['FF_class'].values:
            Rfit_test,Rfit_test_flow = setdata.SolvePickle( DefWipe=DefWipe,HardWipe=HardWipe,Mcore=Mcore,
                                                            RatioFitLists=RatioFitList)
            if CommonRatios:
                if Rfit_test is False:
                    print('Warning, rfitindicies not found, must be from old code')
                    setdata.SolvePickle(DefWipe=True,HardWipe=HardWipe,Mcore=Mcore,
                                        RatioFitLists=RatioFitList)
                RatioFitList = Rfit_test,Rfit_test_flow
            FFtimer.Lap()
        # for setdata in self.SetFF['FF_class'].values:
        #     setdata.LoadPickle(DefWipe=False,HardWipe=HardWipe,Mcore=False)
        # self.SetPlotParams()
        # self.FitFFs()

    def WipeAllFits(self):
        for setdata in self.SetFF['FF_class'].values:
            setdata.WipeAllFits()


    def LoadPickle( self,DefWipe=False,HardWipe=True,Mcore=defMcore,WipeFit=False,
                    OnlyFits=False,DoFits=True,CommonRatios=False):
        # doubling this should fix some memory issues with the calculation.
        FFtimer = Timer(linklist=list(self.SetFF['FF_class'].keys()),name='Form Factor Load '+self.name,rewrite=False)
        RatioFitList = [False,False]
        for setdata in self.SetFF['FF_class'].values:
            setdata.LoadPickle( DefWipe=DefWipe,HardWipe=HardWipe,Mcore=Mcore,
                                OnlyFits=OnlyFits,DoFits=DoFits,RatioFitLists=RatioFitList)
            if CommonRatios:
                RatioFitList = setdata.RfitIndicies,setdata.RfitIndicies_flow
            FFtimer.Lap()
        # for setdata in self.SetFF['FF_class'].values:
        #     setdata.LoadPickle(DefWipe=False,HardWipe=HardWipe,Mcore=False)
        self.SetPlotParams()
        # self.FitFFs(WipeFit=WipeFit)

    def FitFFs(self,thisFFList='All',list_of_fits='PreDef',WipeFit=False):
        if thisFFList == 'All':
            if list_of_fits=='PreDef':
                for iset in  self.SetFF['FF_class'].values:
                    iset.FitFFs(WipeFit=WipeFit)
            else:
                for ifit,iset in zip(list_of_fits, self.SetFF['FF_class'].values):
                    iset.FitFFs(fit_ranges=[ifit],WipeFit=WipeFit)
        else:
            for iFF in thisFFList:
                if list_of_fits=='PreDef':
                    for iset in  self.SetFF['FF_class'].values:
                        iset.FitFFs(thisFF=iFF,WipeFit=WipeFit)
                else:
                    for ifit,iset in zip(list_of_fits, self.SetFF['FF_class'].values):
                        iset.FitFFs(thisFF=iFF,WipeFit=WipeFit)


    def GetSet(self,thisset):
        if isinstance(thisset,int):
            if thisset >= len(self.SetFF['FF_class']):
                raise EnvironmentError('set index for form factor is to large')
            return self.SetFF['FF_class'].iloc[thisset]
        else:
            return self.SetFF.loc[thisset,'FF_class']

    def CombineFFs(self,key_list,coeff_list,this_fun='add',show_timer=False):
        if isinstance(this_fun,str):
            fun_str = this_fun
            this_fun = op_dict[this_fun]
        else:
            fun_str = this_fun.__name__

        if show_timer:
            print('Combining form factors with', fun_str)
            thistimer = Timer(  linklist=key_list,name='FF comb')

        this_iter = iter(zip(key_list,coeff_list))
        first_key,first_coeff = next(this_iter)

        this_set = self.GetSet(first_key)
        if first_coeff == 1.0:
            running_set = this_set
        else:
            running_set = this_set * first_coeff
        if show_timer: thistimer.Lap(this_set.name)
        for ikey,icoeff in this_iter:
            this_set = self.GetSet(ikey)
            if icoeff == 1.0:
                running_set = this_fun(running_set,this_set)
            else:
                running_set = this_fun(running_set,this_set * icoeff)
            if show_timer: thistimer.Lap(this_set.name)
        running_set.Stats()
        set_len = len(self.SetFF.loc[:,'FF_class'])+1
        running_set.ImportPlotParams(   colourset8[set_len%len(colourset8)],
                                        markerset[set_len%len(markerset)],
                                        shiftset[set_len%len(shiftset)])
        this_row = pa.Series()
        this_row.loc['FF_class'] = running_set
        this_row.loc['color'] =   colourset8[set_len%len(colourset8)]
        this_row.loc['symbol'] = markerset[set_len%len(markerset)]
        this_row.loc['shift'] =  shiftset[set_len%len(shiftset)]
        this_row.loc['info_dictionarie'] = running_set.Info
        this_row.loc['q2_list'] = running_set.q2list
        this_row.loc['pp2_list'] = running_set.pp2list
        this_row.loc['2pt_fit_range'] = running_set.fit2ptin
        this_row.loc['alpha_fit_range'] = running_set.fitAlphain
        # this_row.loc['ratio_fit_range'] =
        this_row.loc['config_list'] = running_set.cfglistin
        this_row.loc['config_list_2pt'] = running_set.cfglist2pt
        self.SetFF.loc[running_set.name] = this_row


    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def Plot(self,thisFFlist,plot_info,plottflowlist='PreDefine',chi_list='Chi0-Chi20',WipeData=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r' Form Factor plot '
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$ Q^{2} [GeV]^{2}$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$Form Factor$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'.pdf'
            if 'xlims' in this_pi:
                this_pi['xlims'] = [0,this_pi['xlims'][1]]
            else:
                this_pi['xlims'] = [0,None]
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        if isinstance(chi_list,str):
            chi_list = [chi_list]*len(list(self.SetFF.keys()))
        if isinstance(plottflowlist,str):
            plottflowlist = [plottflowlist]*len(list(self.SetFF.keys()))
        # for thisFF,itflow,iset in zip(thisFFlist,plottflowlist,self.SetFF['FF_class'].values):
        #     data_plot = iset.ValsAndFitPlot(data_plot,thisFF,plottflow=itflow,this_chi_list=chi_list)
        for ic,iset in enumerate(self.SetFF['FF_class'].values):
            if ic >= len(thisFFlist):
                thisFF,itflow = thisFFlist[-1],plottflowlist[-1]
            else:
                thisFF,itflow = thisFFlist[ic],plottflowlist[ic]
            data_plot = iset.ValsAndFitPlot(data_plot,thisFF,plottflow=itflow,this_chi_list=chi_list)
            # iset.Plot(thisFF,plottflow=itflow)
            # if itflow == 'PreDefine':
            #     iset.FitPlot(thisFF,tflow='First')
            # else:
            #     iset.FitPlot(thisFF,tflow=itflow)
        # thisFF = thisFF.replace('tilde',r'\tilde')
        # thisFF = thisFF.replace('frac',r'\frac')
        data_plot.PlotAll()
        self.Write(WipeData=WipeData)
        return data_plot



    def __str__(self):
        return self.name


    def PlotFitVsUDMass(self,thisFFlist,plot_info,plottflowlist='PreDefine',chi_list='Chi0'):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$F_{Q^{2}=0}$ Extrapolation vs Up-Down Mass'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$ m^{MSbar}_{ud} [MeV]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$FF(Q^{2}=0)$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_F0VsUDmass.pdf'
            if 'xlims' in this_pi:
                this_pi['xlims'] = [0,this_pi['xlims'][1]]
            else:
                this_pi['xlims'] = [0,None]
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        if isinstance(chi_list,str):
            chi_list = [chi_list]*len(list(self.SetFF.keys()))
        for ic,ichi in enumerate(chi_list):
            if '-' in ichi:
                chi_list[ic] = ichi.split('-')[0]
        if plottflowlist == 'PreDefine':
            plottflowlist = ['PreDefine']*len(list(self.SetFF.keys()))
        ploty,plotyerr,plotx = [],[],[]
        the_ff,the_tflow,the_DS = '','',''
        firstset = self.SetFF['FF_class'].values[0]
        if any([firstset.DS != iset.DS for iset in self.SetFF['FF_class'].values]):
            print('warning: not all UD combinations are the same for mpi plot')
        else:
            the_DS = firstset.DS
        if any([thisFFlist[0] != iff for iff in thisFFlist] ):
            print('warning: not all form factors are the same for mpi plot')
        else:
            the_ff = thisFFlist[0]
        if firstset.doflow:
            if any([plottflowlist[0] != itflow for itflow in plottflowlist] ) :
                print('warning: not all flow times are the same for mpi plot')
            else:
                the_tflow = tflowTOLeg(plottflowlist[0])
        # for thisFF,this_chi,itflow,(ikey,iset) in zip(thisFFlist,chi_list,plottflowlist,self.SetFF['FF_class'].iteritems()):
        for ic,iset in enumerate(self.SetFF['FF_class'].values):
            if ic >= len(thisFFlist):
                thisFF,this_chi,itflow = thisFFlist[-1],chi_list[-1],plottflowlist[-1]
            else:
                thisFF,this_chi,itflow = thisFFlist[ic],chi_list[ic],plottflowlist[ic]
            ## flowed opp
            ## self.FitDict [itflow , iFF , ifitr]
            ## else
            ## self.FitDict [iFF , ifitr]

            ## we just assume the first ifitr is all the data.
            if iset.doflow:
                if 'boot' not in iset.FF_Fit_Stats or (itflow,this_chi,thisFF) not in list(iset.FF_Fit_Stats['boot'].keys()):
                    # print itflow,this_chi,thisFF ,' fit not in \n', iset.FF_Fit_Stats['boot'].keys(), '\n fitting now'
                    print(itflow,this_chi,thisFF ,' fit not in ',iset.name,' fitting now')
                    iset.Fit(thisFF)
                if 'boot' not in iset.FF_Fit_Stats or (itflow,this_chi,thisFF) not in list(iset.FF_Fit_Stats['boot'].keys()):
                    print(itflow,this_chi,thisFF, 'not found, skipping plot')
                    continue
                this_fit_instance = iset.FF_Fit_Stats['boot'][itflow,this_chi,thisFF].values[0]
                ## assume first parameter is zero intercept
                thisres = this_fit_instance.fit_data['Params'].values[0]
                ploty.append(thisres.Avg)
                plotyerr.append(thisres.Std)
                plotx.append(iset.ppparams.GetUDMass(MeV=True)) ## in MeV
            else:
                if 'boot' not in iset.FF_Fit_Stats or (this_chi,thisFF) not in list(iset.FF_Fit_Stats['boot'].keys()):
                    # print this_chi,thisFF ,' fit not in \n', iset.FF_Fit_Stats['boot'].keys(), '\n fitting now'
                    print(this_chi,thisFF ,' fit not in ',iset.name,'fitting now')
                    iset.Fit(thisFF)
                if 'boot' not in iset.FF_Fit_Stats or (this_chi,thisFF) not in list(iset.FF_Fit_Stats['boot'].keys()):
                    print(this_chi,thisFF, 'not found, skipping plot')
                    continue
                    # print this_chi,thisFF ,' fit not in \n', iset.FF_Fit_Stats['boot'].keys(), '\n fitting now'
                this_fit_instance = iset.FF_Fit_Stats['boot'][this_chi,thisFF].values[0]
                ## assume first parameter is zero intercept
                thisres = this_fit_instance.fit_data['Params'].values[0]
                ploty.append(thisres.Avg)
                plotyerr.append(thisres.Std)
                plotx.append(iset.ppparams.GetUDMass(MeV=True)) ## in MeV
        # thisFF = thisFF.replace('tilde',r'\tilde')
        # thisFF = thisFF.replace('frac',r'\frac')
        thislab = r'$'+r'\ '.join([the_ff.replace(r'frac',r'\frac'),the_tflow,the_DS])+r'$ '+this_chi
        hold_series = pa.Series()
        hold_series['x_data'] = plotx
        hold_series['y_data'] = ploty
        hold_series['yerr_data'] = plotyerr
        hold_series['xerr_data'] = None
        hold_series['type'] = 'error_bar'
        hold_series['fit_class'] = None
        hold_series['label'] = thislab
        hold_series['symbol'] = 'o'
        hold_series['color'] = None
        hold_series['shift'] = None
        hold_series['xdatarange'] = None
        hold_series['scale'] = None
        hold_series['Phys'] = None
        hold_series['ShowPar'] = None
        hold_series['ShowEval'] = None
        data_plot.AppendData(hold_series)
        data_plot.PlotAll()
        return data_plot


    def PlotDist(self,thisFFlist,plot_info,plottflowlist='PreDefine',chi_colapse=True,Vanish_Clim=False,Ord_Two=False
                                                                        ,latspace_fit=(False,False)):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$F_{Q^{2}=0}$ Distribution over All Fit Range Bootstraps'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$FF(Q^{2}=0)$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$Frequency$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Distro.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        if plottflowlist == 'PreDefine':
            plottflowlist = ['PreDefine']*len(list(self.SetFF.keys()))


        # for ichi,thisFF,itflow,(ikey,iset) in zip(  chi_colapse,thisFFlist,plottflowlist,
        #                                             self.SetFF['FF_class'].iteritems()):
        for ic,iset in enumerate(self.SetFF['FF_class'].values):
            if ic >= len(thisFFlist):
                thisFF,ichi,itflow = thisFFlist[-1],chi_colapse[-1],plottflowlist[-1]
            else:
                thisFF,ichi,itflow = thisFFlist[ic],chi_colapse[ic],plottflowlist[ic]
            data_plot = iset.PlotFitDistribution(   data_plot,thisFF,tflow=itflow,thiscol='PreDefine',
                                                    chi_colapse=ichi)

        firstset = self.SetFF['FF_class'].values[0]
        the_mpi,the_latspace = '',''
        if any([firstset.qparams.GetPionMass() != iset.qparams.GetPionMass() for iset in self.SetFF['FF_class'].values]):
            print('warning: not all mpi values are the same for latspace plot')
        else:
            the_mpi = firstset.qparams.GetPionMassLab(MeV=True)
        if any([firstset.qparams.latspace != iset.qparams.latspace for iset in self.SetFF['FF_class'].values]):
            print('warning: not all lattice spaces are the same for mpi plot')
        else:
            the_latspace = firstset.qparams.GetLatSpaceLab()
        last_col = colourset8[(len(list(self.SetFF.keys()))+1)%len(colourset8)]
        if len(the_latspace) == 0:
            if len(the_mpi) == 0:
                data_plot = self.PlotDist_Extrap(   data_plot,'continuum',this_color=last_col,
                                                    this_fun=F3ChiralLimFun[(Vanish_Clim,Ord_Two)+latspace_fit])

                # last_col2 = colourset8[(len(self.SetFF.keys())+2)%len(colourset8)]
                # data_plot = self.PlotDist_Extrap(data_plot,'latspace',this_color=last_col2)
            else:
                data_plot = self.PlotDist_Extrap(data_plot,'latspace',this_color=last_col)
        else:
            if len(the_mpi) == 0:
                data_plot = self.PlotDist_Extrap(data_plot,'chiral',this_color=last_col,this_fun=F3ChiralLimFun[(Vanish_Clim,Ord_Two,False,False)])
        data_plot.PlotAll()
        return data_plot


    def PlotFitVsMpi(self,thisFFlist,plot_info,plottflowlist='PreDefine',Vanish_Clim=False,Ord_Two=False,
                     latspace_fit=(False,False),fitSchiff=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$F_{Q^{2}=0}$ Extrapolation vs Pion Mass'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$ m_{\pi} [MeV]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$FF(Q^{2}=0)$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_F0VsMpi.pdf'
            if 'xlims' in this_pi:
                this_pi['xlims'] = [0,this_pi['xlims'][1]]
            else:
                this_pi['xlims'] = [0,None]
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        if plottflowlist == 'PreDefine':
            plottflowlist = ['PreDefine']*len(list(self.SetFF.keys()))
        the_DS = ''
        firstset = self.SetFF['FF_class'].values[0]
        if any([firstset.DS != iset.DS for iset in self.SetFF['FF_class'].values]):
            print('warning: not all UD combinations are the same for mpi plot')
        else:
            the_DS = firstset.DS

        the_latspace = ''
        firstset = self.SetFF['FF_class'].values[0]
        if any([firstset.qparams.latspace != iset.qparams.latspace for iset in self.SetFF['FF_class'].values]):
            print('warning: not all lattice spaces are the same for mpi plot')
            for iset in self.SetFF['FF_class'].values:
                print(iset.name)
                print(iset.qparams.latspace)
                print()
        else:
            the_latspace = firstset.qparams.GetLatSpaceLab()

        if any([thisFFlist[0] != iff for iff in thisFFlist] ):
            print('warning: not all form factors are the same for mpi plot')
        if firstset.doflow:
            if any([plottflowlist[0] != itflow for itflow in plottflowlist] ) :
                print('warning: not all flow times are the same for mpi plot')
        ploty,plotyerr,plotx = [],[],[]
        this_names = None
        for ikey,iset in self.SetFF['FF_class'].items():
            thisres = iset.Get_Extrapolation()
            if thisres.index.names[0] is not None:
                this_names = thisres.index.names
            for ireskey,iresval in thisres.items():
                this_key = (iset.ppparams.GetPionMassLab(MeV=True),)+RemoveFlowKey(ireskey)## in MeV
                ploty.append(iresval.Avg)
                plotyerr.append(iresval.Std)
                plotx.append(this_key)
        if len(plotx) > 0 and this_names is not None :
            this_key = [slice(None)] + list(plotx[0][1:])
            indicies = pa.MultiIndex.from_tuples(plotx,names=['$m_{pi}$']+this_names)
            ploty = pa.Series(ploty,index=indicies)
            plotyerr = pa.Series(plotyerr,index=indicies)
        else:
            print('No Sets Found for plotting vs mpi')
            return data_plot

        # thisFF = thisFF.replace('tilde',r'\tilde')
        # thisFF = thisFF.replace('frac',r'\frac')
        thislab = r'$'+r'\ '.join([the_DS,the_latspace])+r'$ '
        hold_series = copy(null_series)
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = ploty
        hold_series['yerr_data'] = plotyerr
        hold_series['key_select'] = this_key
        hold_series['type'] = 'error_bar_vary'
        hold_series['label'] = thislab
        hold_series['symbol'] = 'o'
        data_plot.AppendData(hold_series)
        if fitSchiff:
            data_plot = self.Plot_Extrap(data_plot,'continuum',this_fun=(SchiffFitFun,3))
        else:
            if latspace_fit[0]:
                data_plot = self.Plot_Extrap(data_plot,'continuum',this_fun=F3ChiralLimFun[(Vanish_Clim,Ord_Two)+latspace_fit])
            else:
                data_plot = self.Plot_Extrap(data_plot,'chiral',this_fun=F3ChiralLimFun[(Vanish_Clim,Ord_Two)+latspace_fit])
        data_plot.PlotAll()
        return data_plot


    def PlotFitVsLatSpace(self,thisFFlist,plot_info,plottflowlist='PreDefine',mpi_fit=(False,False,False),fitSchiff=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$F_{Q^{2}=0}$ Extrapolation vs Lattice Spacing'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$ a[fm]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$FF(Q^{2}=0)$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_F0VsLatSpace.pdf'
            if 'xlims' in this_pi:
                this_pi['xlims'] = [0,this_pi['xlims'][1]]
            else:
                this_pi['xlims'] = [0,None]
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        if plottflowlist == 'PreDefine':
            plottflowlist = ['PreDefine']*len(list(self.SetFF.keys()))
        the_DS = ''
        firstset = self.SetFF['FF_class'].values[0]
        if any([firstset.DS != iset.DS for iset in self.SetFF['FF_class'].values]):
            print('warning: not all UD combinations are the same for latspace plots')
        else:
            the_DS = firstset.DS
        the_mpi = ''
        firstset = self.SetFF['FF_class'].values[0]
        if any([firstset.qparams.GetPionMass() != iset.qparams.GetPionMass() for iset in self.SetFF['FF_class'].values]):
            print('warning: not all mpi values are the same for latspace plot')
            for iset in self.SetFF['FF_class'].values:
                print(iset.name)
                print(iset.qparams.GetPionMass())
                print()
        else:
            the_mpi = firstset.qparams.GetPionMassLab(MeV=True)

        if any([thisFFlist[0] != iff for iff in thisFFlist] ):
            print('warning: not all form factors are the same for mpi plot')
        if firstset.doflow:
            if any([plottflowlist[0] != itflow for itflow in plottflowlist] ) :
                print('warning: not all flow times are the same for mpi plot')
        ploty,plotyerr,plotx = [],[],[]
        this_names = None
        for ikey,iset in self.SetFF['FF_class'].items():
            thisres = iset.Get_Extrapolation()
            if thisres.index.names[0] is not None:
                this_names = thisres.index.names
            for ireskey,iresval in thisres.items():
                this_key = (iset.ppparams.GetLatSpaceLab(),)+RemoveFlowKey(ireskey)## in MeV
                ploty.append(iresval.Avg)
                plotyerr.append(iresval.Std)
                plotx.append(this_key) ## in MeV
        if len(plotx) > 0 and this_names is not None :
            this_key = [slice(None)] + list(plotx[0][1:])
            indicies = pa.MultiIndex.from_tuples(plotx,names=['$a$']+this_names)
            ploty = pa.Series(ploty,index=indicies)
            plotyerr = pa.Series(plotyerr,index=indicies)
        else:
            print('No Sets Found for plotting vs lat spacing')
            return data_plot

        # thisFF = thisFF.replace('tilde',r'\tilde')
        # thisFF = thisFF.replace('frac',r'\frac')
        thislab = r'$'+r'\ '.join([the_DS,the_mpi])+r'$ '
        hold_series = copy(null_series)
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = ploty
        hold_series['yerr_data'] = plotyerr
        hold_series['key_select'] = this_key
        hold_series['type'] = 'error_bar_vary'
        hold_series['label'] = thislab
        hold_series['symbol'] = 'o'
        data_plot.AppendData(hold_series)
        if fitSchiff:
            data_plot = self.Plot_Extrap(data_plot,'continuum',this_fun=(SchiffFitFun,3))
        else:
            if mpi_fit[0]:
                data_plot = self.Plot_Extrap(data_plot,'continuum',this_fun=F3LatSpaceFun[mpi_fit])
            else:
                data_plot = self.Plot_Extrap(data_plot,'latspace',this_fun=F3LatSpaceFun[mpi_fit])
        data_plot.PlotAll()
        return data_plot



    def PlotFitVsBoxSize(self,thisFFlist,plot_info,plottflowlist='PreDefine'):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$F_{Q^{2}=0}$ Extrapolation vs BoxSize'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$ La[fm]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$FF(Q^{2}=0)$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_F0VsBox.pdf'
            if 'xlims' in this_pi:
                this_pi['xlims'] = [0,this_pi['xlims'][1]]
            else:
                this_pi['xlims'] = [0,None]
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        if plottflowlist == 'PreDefine':
            plottflowlist = ['PreDefine']*len(list(self.SetFF.keys()))
        the_DS = ''
        firstset = self.SetFF['FF_class'].values[0]
        if any([firstset.DS != iset.DS for iset in self.SetFF['FF_class'].values]):
            print('warning: not all UD combinations are the same for box size plots')
        else:
            the_DS = firstset.DS
        the_mpi = ''
        firstset = self.SetFF['FF_class'].values[0]
        if any([firstset.qparams.GetPionMass() != iset.qparams.GetPionMass() for iset in self.SetFF['FF_class'].values]):
            print('warning: not all mpi values are the same for box size plot')
            for iset in self.SetFF['FF_class'].values:
                print(iset.name)
                print(iset.qparams.GetPionMass())
                print()
        else:
            the_mpi = firstset.qparams.GetPionMassLab(MeV=True)

        if any([thisFFlist[0] != iff for iff in thisFFlist] ):
            print('warning: not all form factors are the same for mpi plot')
        if firstset.doflow:
            if any([plottflowlist[0] != itflow for itflow in plottflowlist] ) :
                print('warning: not all flow times are the same for mpi plot')
        ploty,plotyerr,plotx = [],[],[]
        this_names = None
        for ikey,iset in self.SetFF['FF_class'].items():
            thisres = iset.Get_Extrapolation()
            if thisres.index.names[0] is not None:
                this_names = thisres.index.names
            for ireskey,iresval in thisres.items():
                this_key = (iset.ppparams.GetBoxSizeLab(),)+RemoveFlowKey(ireskey)## in MeV
                ploty.append(iresval.Avg)
                plotyerr.append(iresval.Std)
                plotx.append(this_key)
        if len(plotx) > 0 and this_names is not None :
            this_key = [slice(None)] + list(plotx[0][1:])
            indicies = pa.MultiIndex.from_tuples(plotx,names=['$La$']+this_names)
            ploty = pa.Series(ploty,index=indicies)
            plotyerr = pa.Series(plotyerr,index=indicies)
        else:
            print('No Sets Found for plotting vs box size')
            return data_plot

        # thisFF = thisFF.replace('tilde',r'\tilde')
        # thisFF = thisFF.replace('frac',r'\frac')
        thislab = r'$'+r'\ '.join([the_DS,the_mpi])+r'$ '
        hold_series = copy(null_series)
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = ploty
        hold_series['yerr_data'] = plotyerr
        hold_series['key_select'] = this_key
        hold_series['type'] = 'error_bar_vary'
        hold_series['label'] = thislab
        hold_series['symbol'] = 'o'
        data_plot.AppendData(hold_series)
        # data_plot = self.Plot_BoxSize_Extrap(data_plot)
        data_plot.PlotAll()
        return data_plot




    def PlotDist_Extrap(self,plot_data,extrap_type,this_fun='Default',this_color=None):
        extrap_type = extrap_type.lower()
        this_eval = [0]
        if 'latspace' in extrap_type:
            if isinstance(this_fun,str):
                this_fun = (SquaredFitFun,2)
            fit_data = self.LatSpace_Extrapolation(this_fun=this_fun)
            prefix_lab = r'$a\rightarrow 0 Extrap$'
        elif 'boxsize' in extrap_type:
            if isinstance(this_fun,str):
                this_fun = (SquaredFitFun,2)
            fit_data = self.BoxSize_Extrapolation(this_fun=this_fun)
            prefix_lab = r'$La\rightarrow 0 Extrap$'
        elif 'mpi' in extrap_type or 'chiral' in extrap_type:
            if isinstance(this_fun,str):
                this_fun = (LinearFitFun,2)
            fit_data = self.Chiral_Extrapolation(this_fun=this_fun)
            this_eval = [physcial_pion_mass]
            prefix_lab = r'$m_{\pi}\rightarrow '+str(int(physcial_pion_mass))+'MeV Extrap$'
        elif 'cont' in extrap_type:
            if isinstance(this_fun,str):
                this_fun = (LinearFF_2D,3)
            fit_data = self.Continuum_Extrapolation(this_fun=this_fun)
            this_eval = [0,physcial_pion_mass]
            prefix_lab = r'$m_{\pi}\rightarrow '+str(int(physcial_pion_mass))+'MeV, a\rightarrow 0 Extrap$'
        else:
            raise IOError('extrapolation type not recognised '+str(extrap_type))
        if len(fit_data) == 0:
            return plot_data
        def GetEval(val):
            return val.Eval_Function(*this_eval)
        fit_data = fit_data.apply(GetEval)
        hold_series = copy(null_series)
        hold_series['key_select'] = fit_data.index[0]
        hold_series['type'] = 'histogram_vary'
        hold_series['boot_data'] = fit_data
        hold_series['label'] = prefix_lab
        hold_series['color'] = this_color
        plot_data.AppendData(hold_series)
        return plot_data


    def Plot_Extrap(self,plot_data,extrap_type,this_fun='Default'):
        extrap_type = extrap_type.lower()
        this_eval = 0
        if 'latspace' in extrap_type:
            if isinstance(this_fun,str):
                this_fun = (SquaredFitFun,2)
            fit_data = self.LatSpace_Extrapolation(this_fun=this_fun)
        elif 'boxsize' in extrap_type:
            if isinstance(this_fun,str):
                this_fun = (SquaredFitFun,2)
            fit_data = self.BoxSize_Extrapolation(this_fun=this_fun)
        elif 'mpi' in extrap_type or 'chiral' in extrap_type:
            if isinstance(this_fun,str):
                this_fun = (LinearFitFun,2)
            fit_data = self.Chiral_Extrapolation(this_fun=this_fun)
            this_eval = physcial_pion_mass
        elif 'cont' in extrap_type:
            if isinstance(this_fun,str):
                this_fun = (LinearFF_2D,3)
            fit_data = self.Continuum_Extrapolation(this_fun=this_fun)
            this_eval = [0,physcial_pion_mass]
        else:
            raise IOError('extrapolation type not recognised '+str(extrap_type))
        if len(fit_data) == 0:
            return plot_data
        # if isinstance(fitr_params,pa.Series): fitr_params = fitr_params.iloc[-1]
        ## TODO, plot from zero to data max
        hold_series = copy(null_series)
        hold_series['key_select'] = fit_data.index[0]
        hold_series['type'] = 'fit_vary'
        hold_series['fit_class'] = fit_data
        hold_series['label'] = 'fit '
        hold_series['xdatarange'] = [0,'ToDataMax']
        hold_series['ShowEval'] = this_eval
        plot_data.AppendData(hold_series)
        return plot_data


    def Chiral_Extrapolation(self,this_fun=(LinearFitFun,2)):
        def this_fit(this_series):
            xdata,ydata = [[]],[]
            for icol,coldata in this_series.items():
                if '_mpi' in icol:
                    xdata[0].append(coldata)
                    # ##### DEBUGGING ####
                    # xdata[0].append(coldata+100)
                else:
                    ydata.append(coldata)
                    # ##### DEBUGGING ####
                    # coldata2,rand_list = coldata.ResampleBS()
                    # ydata.append(coldata2)

            if len(list(set(xdata[0]))) < this_fun[1]:
                return float('NaN')

            # max_nboot = max([ival.nboot for ival in ydata])
            fmt_ydata,fmt_xdata = [],[[]]
            for icy,(ix,iy) in enumerate(zip(xdata[0],ydata)):
                if not isinstance(iy,float):
                    if nboot != iy.nboot:
                        data,rand_list = iy.ResampleBS(nboot)
                        fmt_xdata[0].append(ix)
                        fmt_ydata.append(data)
                    else:
                        fmt_xdata[0].append(ix)
                        fmt_ydata.append(iy)

            thisff = ff.Fitting(Funs=this_fun,data=[np.array(fmt_xdata),np.array(fmt_ydata)],name='Chiral_Ext')
            fit_result = thisff.FitBoots()
            if np.isnan(fit_result):
                return fit_result
            # print thisff.fit_data['Params']
            return thisff

        this_df = pa.DataFrame()
        for ikey,iset in self.SetFF['FF_class'].items():
            this_data = iset.Get_Extrapolation()
            if len(this_data) > 0:
                this_df[ikey] = this_data
                this_df.loc[:,ikey+'_mpi'] = int(iset.qparams.GetPionMass(MeV=True))
        this_df.loc[:,'Fit'] = this_df.apply(this_fit,axis=1)
        return this_df['Fit'].dropna()


    def LatSpace_Extrapolation(self,this_fun=(SquaredFitFun,2)):
        def this_fit(this_series):
            xdata,ydata = [[]],[]
            for icol,coldata in this_series.items():
                if '_latspace' in icol:
                    xdata[0].append(coldata)
                    # ##### DEBUGGING ####
                    # xdata[0].append(coldata+100)
                else:
                    ydata.append(coldata)
                    # ##### DEBUGGING ####
                    # ydata.append(coldata)
            if len(list(set(xdata[0]))) < this_fun[1]:
                return float('NaN')
            # max_nboot = max([ival.nboot for ival in ydata])
            fmt_ydata,fmt_xdata = [],[[]]
            for icy,(ix,iy) in enumerate(zip(xdata[0],ydata)):
                if not isinstance(iy,float):
                    if nboot != iy.nboot:
                        data,rand_list = iy.ResampleBS(nboot)
                        fmt_xdata[0].append(ix)
                        fmt_ydata.append(data)
                    else:
                        fmt_xdata[0].append(ix)
                        fmt_ydata.append(iy)

            thisff = ff.Fitting(Funs=this_fun,data=[np.array(fmt_xdata),np.array(fmt_ydata)],name='LatSpace_Ext')
            fit_result = thisff.FitBoots()
            if np.isnan(fit_result):
                return fit_result
            return thisff

        this_df = pa.DataFrame()
        for ikey,iset in self.SetFF['FF_class'].items():
            this_data = iset.Get_Extrapolation()
            if len(this_data) > 0:
                this_df[ikey] = this_data
                this_df.loc[:,ikey+'_latspace'] = float(iset.ppparams.latspace)
        this_df.loc[:,'Fit'] = this_df.apply(this_fit,axis=1)
        return this_df['Fit'].dropna()


    def BoxSize_Extrapolation(self,this_fun=(SquaredFitFun,2)):
        def this_fit(this_series):
            xdata,ydata = [[]],[]
            for icol,coldata in this_series.items():
                if '_boxsize' in icol:
                    xdata[0].append(coldata)
                    # ##### DEBUGGING ####
                    # xdata[0].append(coldata+100)
                else:
                    ydata.append(coldata)
                    # ##### DEBUGGING ####
                    # ydata.append(coldata)
            if len(list(set(xdata[0]))) < this_fun[1]:
                return float('NaN')
            # max_nboot = max([ival.nboot for ival in ydata])
            fmt_ydata,fmt_xdata = [],[[]]
            for icy,(ix,iy) in enumerate(zip(xdata[0],ydata)):
                if not isinstance(iy,float):
                    if nboot != iy.nboot:
                        data,rand_list = iy.ResampleBS(nboot)
                        fmt_xdata[0].append(ix)
                        fmt_ydata.append(data)
                    else:
                        fmt_xdata[0].append(ix)
                        fmt_ydata.append(iy)

            thisff = ff.Fitting(Funs=this_fun,data=[np.array(fmt_xdata),np.array(fmt_ydata)],name='BoxSize_Ext')
            fit_result = thisff.FitBoots()
            if np.isnan(fit_result):
                return fit_result
            return thisff

        this_df = pa.DataFrame()
        for ikey,iset in self.SetFF['FF_class'].items():
            this_data = iset.Get_Extrapolation()
            if len(this_data) > 0:
                this_df[ikey] = this_data
                this_df.loc[:,ikey+'_boxsize'] = float(iset.ppparams.Lxyz)
        this_df.loc[:,'Fit'] = this_df.apply(this_fit,axis=1)
        return this_df['Fit'].dropna()




    def Continuum_Extrapolation(self,this_fun=(LinearFF_2D,3)):
        def this_fit(this_series):
            xdata,ydata = [[],[]],[]
            for icol,coldata in this_series.items():
                if '_latspace' in icol:
                    xdata[0].append(coldata)
                    ### DEBUG ###
                    # xdata[0].append(coldata+0.1)
                    # xdata[0].append(coldata)
                    # xdata[0].append(coldata+0.1)
                    # xdata[0].append(coldata)
                    # xdata[0].append(coldata+0.1)
                elif '_mpi' in icol:
                    xdata[1].append(coldata)
                    ### DEBUG ###
                    # xdata[1].append(coldata)
                    # xdata[1].append(coldata+100)
                    # xdata[1].append(coldata+100)
                    # xdata[1].append(coldata+400)
                    # xdata[1].append(coldata+400)
                else:
                    ydata.append(coldata)
                    ### DEBUG ###
                    # ydata.append(coldata)
                    # ydata.append(coldata)
                    # ydata.append(coldata)
                    # ydata.append(coldata)
                    # ydata.append(coldata)
            if len(list(set(xdata[0]))) + len(list(set(xdata[1]))) < this_fun[1] :
                return float('NaN')
            # max_nboot = max([ival.nboot for ival in ydata])
            fmt_ydata,fmt_xdata = [],[[],[]]
            for icy,(ix0,ix1,iy) in enumerate(zip(xdata[0],xdata[1],ydata)):
                if not isinstance(iy,float):
                    if nboot != iy.nboot:
                        data,rand_list = iy.ResampleBS(nboot)
                        fmt_xdata[0].append(ix0)
                        fmt_xdata[1].append(ix1)
                        fmt_ydata.append(data)
                    else:
                        fmt_xdata[0].append(ix0)
                        fmt_xdata[1].append(ix1)
                        fmt_ydata.append(iy)
            thisff = ff.Fitting(Funs=this_fun,data=[np.array(fmt_xdata),np.array(fmt_ydata)],name='Continuum_Ext')
            fit_result = thisff.FitBoots()
            if np.isnan(fit_result):
                return fit_result
            return thisff

        this_df = pa.DataFrame()
        for ikey,iset in self.SetFF['FF_class'].items():
            this_data = iset.Get_Extrapolation()
            if len(this_data) > 0:
                this_df[ikey] = this_data
                this_df.loc[:,ikey+'_latspace'] = float(iset.Info['a'])
                this_df.loc[:,ikey+'_mpi'] = int(iset.ppparams.GetPionMass(MeV=True))
        this_df.loc[:,'Fit'] = this_df.apply(this_fit,axis=1)
        return this_df['Fit'].dropna()




def TestSetFF(this_file_list):
    return SetOfFF.Construct_FF_Set(this_file_list)


if __name__ == '__main__':
    this_file_list = []
    this_file_list.append('/home/jackdra/LQCD/Results/temp/VectorTop_t_f7.47_Rtsumfitr4-32_Afitr10-20_RFmin3_RChi50-100.p')
    this_data = TestSetFF(this_file_list)
    this_data.SetFF['FF_class'].values[0].FF
