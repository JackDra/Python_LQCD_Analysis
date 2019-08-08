#!/usr/bin/env python

import numpy as np

from collections import OrderedDict
from MiscFuns import NumbToFileName
from Params import graphdir,defFitDict
# from Params import defxlimAlpha,defxlimAlphaTflow,,defxlim
from MiscFuns import mkdir_p,CombineListCfgs
from copy import deepcopy,copy
from XmlFormatting import tflowstr#,tflowTOLeg
# import AhmedRes as ara
from PredefFitFuns import Chipt_Chit
import FitFunctions as ff
import PlotData as jpl
from NullPlotData import null_series
import pandas as pa
from FileIO import Construct_Empty

# import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
# import matplotlib.pyplot as pl
from QuantityLists import top_susc_0710_1130_plot,top_susc_1406_5363_plot
from QuantityLists import top_susc_1312_5161_plot,top_susc_1312_5161_211_plot
from QuantityLists import top_susc_1512_06746_plot,top_susc_1705_10906_plot

## Move Somewhere?
# params = {'legend.fontsize': 7,
#           'legend.numpoints': 1,
#           'axes.labelsize' : 20,
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




class SetOfFO(object):
    """

    FO(set) Inherits from FlowOp in FlowOpps.py


    plotting parameters (Move to different file?) ####


    InfoDict is a list of Info dictionaries to be bassed to TwoPtCorr:


    """
    def Construct_FO_Set(this_file_list):
        set_list = [FlowOp.Construct_FO_File(ifile) for ifile in this_file_list]
        out_class = Construct_Empty(SetOfFO)
        out_class.cfglist = {}
        out_class.InfoDict = []
        out_class.Plotdir = graphdir+'/FlowOp/'
        out_class.name = 'Set1'
        out_class.SetFO = {}
        for iset in set_list:
            out_class.SetFO[iset.name] = iset
        return out_class

    instances = []

    parentlist = []

    def __init__(self, cfglist={}, InfoDict=[],name = '',corr_cfglists=True,CompAhmed=False):
        from FlowOpps import FlowOp

        self.cfglist = cfglist
        self.InfoDict = InfoDict
        self.Plotdir = graphdir+'/FlowOp/'
        mkdir_p(self.Plotdir)
        self.AComp = CompAhmed

        # self.FullFO = 'Not Set'
        if len(name) == 0:
            self.name = 'Set'+str(int(len(SetOfFO.instances))+1)
            SetOfFO.instances.append(self.name)
        else:
            self.name = name
            # SetOfFO.instances.append(self.name)
        self.AutoList = []
        # self.RandTList = []
        self.IOverList = []

        self.SetFO = OrderedDict()
        for ic,idict in enumerate(self.InfoDict):
            thiscol = colourset8[ic%len(colourset8)]
            thissym = markerset[ic%len(markerset)]
            thisshift = shiftset[ic%len(shiftset)]

            if 'autocorr' in list(idict.keys()):
                self.AutoList.append(idict['autocorr'])
            else:
                self.AutoList.append(False)

            if 'IOver' in list(idict.keys()):
                self.IOverList.append(idict['IOver'])
            else:
                self.IOverList.append(False)


            if 'ObsComb' in list(idict.keys()) and len(idict['ObsComb']) > 0:
                thislist = []
                for iObs in idict['ObsComb']:
                    thisdict = deepcopy(idict)
                    thisdict['Observable'] = iObs
                    thislist.append(thisdict)
                thisFO = [FlowOp(cfglist=self.cfglist, Info=iterdict,thissym=thissym,thiscol=thiscol,thisshift=thisshift) for iterdict in thislist]
            else:
                thisFO = FlowOp(cfglist=self.cfglist, Info=idict,thissym=thissym,thiscol=thiscol,thisshift=thisshift)

            if str(thisFO) in list(self.SetFO.keys()):
                FOnext = 2
                while str(thisFO)+'_'+str(FOnext) in list(self.SetFO.keys()):
                    FOnext += 1
                self.SetFO[str(thisFO)+'_'+str(FOnext)] = thisFO
            else:
                self.SetFO[str(thisFO)] = thisFO
            # print 'Setting up classes for ' , thisFO , ' complete'
            # print 'has color',thiscol , 'symbol ' , thissym , 'thisshift',thisshift

        self.SetFOlen = len(self.SetFO)
        if len(list(self.cfglist.keys())) == 0:
            if corr_cfglists:
                self.GetCfgList()
            else:
                self.GetUnCorrCfgList()
        ##Implement Reading cfglist if len(cfglist) == 0

        self.current = 0


    # ## thisfun(corr1,corr2,corr3,...)
    # ## Strongly encoraged to manually set the filename somewhere
    # ## e.g. sink smear combine we should get the first name and replace the sink smear with the sink smear label
    # def CombineCorrs(self,thisfun,filename='',LegLab='',jsmlab = ''):
    #     # for items in self.SetFO.values():
    #     #     print ''
    #     #     print items
    #     #     for iitems in items.itemsAvgStd():
    #     #         print iitems
    #     output2pt = thisfun(*self.SetFO.values())
    #     output2pt.jsm = jsmlab
    #     output2pt.SetCustomName(string=filename,stringLL=LegLab)
    #     output2pt.Bootstrap()
    #     output2pt.Stats()
    #     output2pt.Write()
    #     return output2pt

    def GetUnCorrCfgList(self):
        for iFO in list(self.SetFO.values()):
            if isinstance(iFO,list):
                for jFO in iFO:
                    jFO.FlowGetCfgList()
            else:
                iFO.FlowGetCfgList()


    def GetCfgList(self):
        self.GetUnCorrCfgList()
        input_cfglist = []
        for iset in self.SetFO.values():
            if isinstance(iset,(list,tuple,np.ndarray)):
                for jset in iset:
                    input_cfglist.append(jset)
            else:
                input_cfglist.append(iset)
        cfglist = CombineListCfgs(input_cfglist,'Op_cfgs')
        self.ImportCfgList(cfglist)

    def ImportCfgList(self,cfglist):
        self.cfglist = cfglist
        for iFO in self.SetFO.values():
            if isinstance(iFO,list):
                for jFO in iFO:
                    jFO.FlowImportCfgList(cfglist)
            else:
                iFO.FlowImportCfgList(cfglist)

    def Read(self):
        for setdata in self.SetFO.values():
            if isinstance(setdata,list):
                for jFO in setdata:
                    print('Importing 2pt corr: ' , str(jFO))
                    jFO.FlowRead()
            else:
                print('Importing 2pt corr: ' , str(setdata))
                setdata.FlowRead()


    def Write(self):
        for iFO in self.SetFO.values():
            if isinstance(iFO,list):
                for jFO in iFO:
                    print('Writing 2pt corr: ' , str(jFO))
                    jFO.FlowWrite()
            else:
                print('Writing 2pt corr: ' , str(iFO))
                iFO.FlowWrite()

    def LoadPickle(self,DefWipe=False,CheckCfgs=True):
        for ic,(setkey,setdata) in enumerate(self.SetFO.items()):
            if isinstance(setdata,list):
                # print 'LoadPickle flowed operators: ' +', '.join(map( str,setdata))
                for icount,jFO in enumerate(setdata):
                    ## first statement gets checked before second, so this is ok!.
                    if icount > 0 and jFO.Observ == setdata[icount-1].Observ:
                        self.SetFO[setkey][icount] = deepcopy(setdata[icount-1])
                    else:
                        jFO.FlowLoadPickle(DefWipe=DefWipe,CheckCfgs=CheckCfgs)
                    self.SetFO[setkey][icount].FlowImportPlotParams(colourset8[ic%len(colourset8)],markerset[ic%len(markerset)]
                                                                    ,shiftset[ic%len(shiftset)])
            else:
                # print 'LoadPickle flowed operator: ' , str(setdata)
                setdata.FlowLoadPickle(DefWipe=DefWipe,CheckCfgs=CheckCfgs)
                self.SetFO[setkey].FlowImportPlotParams(colourset8[ic%len(colourset8)],markerset[ic%len(markerset)]
                                                        ,shiftset[ic%len(shiftset)])
        # self.MakeRandT(DefWipe=DefWipe)
        print()

    # def MakeRandT(self,DefWipe=False):
    #     self.FullFO = deepcopy(self.SetFO)
    #     for ic,(setkey,setdata) in enumerate(self.FullFO.iteritems()):
    #         if self.RandTList[ic]:
    #             if isinstance(setdata,list):
    #                 print 'Doing RandT for flowed operators: ' +', '.join(map( str,setdata))
    #                 for icount,jFO in enumerate(setdata):
    #                     ## first statement gets checked before second, so this is ok!.
    #                     ## No Input into PickTime corresponds to last time
    #                     self.SetFO[setkey][icount] = jFO.PickTime()
    #             else:
    #                 print 'Doing RandT for flowed operator: ' , str(setdata)
    #                 # self.SetFO[setkey][icount] = jFO.PickTime(DefWipe=DefWipe)[0]
    #                 self.SetFO[setkey] = setdata.PickTime()
    #     print


    def FlowMakeChit(self,DefWipe=False,Improve=False):
        for ic,(ikey,setdata) in enumerate(self.SetFO.items()):
            if isinstance(setdata,list):
                if len(setdata) != 2 :
                    raise EnvironmentError('Chit only set up for 2 operators at the moment')
                thisdata = setdata[0].FlowMakeChit(setdata[1],DefWipe=DefWipe,Improve=self.IOverList[ic])
                self.SetFO[ikey] = thisdata


    def FlowMakeOppSquared(self):
        for setdata in self.SetFO.values():
            setdata.FlowMakeOppSquared()


    def ReadAndWrite(self,DefWipe=False):
        for ic,(setkey,setdata) in enumerate(self.SetFO.items()):
            if isinstance(setdata,list):
                for icount,jFO in enumerate(setdata):
                    print('Importing and Writing 2pt corr: ' , str(jFO))
                    jFO.FlowReadAndWrite(DefWipe=DefWipe)
            else:
                print('Importing and Writing 2pt corr: ' , str(setdata))
                setdata.FlowReadAndWrite(DefWipe=DefWipe)

    ## returns a TwoPtCorr class object (see TwoPtCorrelators.py)
    def GetSet(self,thisset):
        return self.SetFO[thisset]


    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def PlotTauInt(self,tflowlist,plot_info):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Susceptibility'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$\sqrt{8t_{f}}[fm]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\chi$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_TauInt_'+tflowlist[0]+'.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        for ic,(itflow,(ikey,iFO)) in enumerate(zip(tflowlist,iter(self.SetFO.items()))):
            tflow = tflowstr(itflow)
            if isinstance(iFO,list):
                raise EnvironmentError('Please run combine opperators before running plot')
            if self.AutoList[ic]:
                print('plotting Tau Int ', iFO)
                if tflow not in list(iFO.Op_Stats['Auto'].keys()):
                    print('Warning , flow time ',tflow,' not found ')
                else:
                    data_plot = iFO.Op_Stats['Auto'][tflow].PlotTauInt(data_plot,
                                            thiscol=iFO.thiscol,thissym=iFO.thissym,thisshift=iFO.thisshift)
            else:
                print('Warning, ' , iFO , ' is not autocorrelated analized')
        data_plot.PlotAll()
        return data_plot

    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def PlotWopt(self,tflowlist,plot_info):
        # if ftitle=='Def': ftitle=r'Deturmination of Optimal parameter $W$ at flow time $'+tflowTOLeg(tflow)+'$'
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$\rho $ for corresponding $W$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$\rho$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$W$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Wopt_'+tflowlist[0]+'.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)

        for ic,(itflow,(ikey,iFO)) in enumerate(zip(tflowlist,iter(self.SetFO.items()))):
            tflow = tflowstr(itflow)
            if isinstance(iFO,list):
                raise EnvironmentError('Please run combine opperators before running plot')
            if self.AutoList[ic]:
                print('plotting W optimal ', iFO)
                if tflow not in list(iFO.Op_Stats['Auto'].keys()):
                    print('Warning , flow time ',tflow,' not found ')
                else:
                    data_plot = iFO.Op_Stats['Auto'][tflow].PlotWopt(data_plot,thiscol=iFO.thiscol,thissym=iFO.thissym,thisshift=iFO.thisshift)
            else:
                print('Warning, ' , iFO , ' is not autocorrelated analized')
        data_plot.PlotAll()
        return data_plot



    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def PlotMonteTime(self,tflowlist,plot_info):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Monte Carlo Time of Flowed Opp'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'Monte Carlo Time'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'Flowed Operator'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Monte.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)

        icount = 0
        for ic,(itflow,(ikey,iFO)) in enumerate(zip(tflowlist,iter(self.SetFO.items()))):
            if isinstance(iFO,list):
                raise EnvironmentError('Please run combine opperators before running plot')
            print('plotting ', iFO)
            for istream in iFO.stream_list:
                thiscol = colourset8[icount%len(colourset8)]
                thisshift = shiftset[icount%len(shiftset)]
                data_plot = iFO.PlotMonteTime(data_plot,itflow,stream=istream,thiscol=thiscol,thisshift=thisshift)
                icount += 1
        data_plot.PlotAll()
        return data_plot


    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def PlotMonteTimeDeltaEps(self,tflowlist,plot_info,epslist,MultQ=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Monte Carlo Time of Delta Function'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'Monte Carlo Time'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\delta_{FO< \epsilon }$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_DeltaMonte.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)

        icount = 0
        for ic,(ieps,itflow,(ikey,iFO)) in enumerate(zip(epslist,tflowlist,iter(self.SetFO.items()))):
            if isinstance(iFO,list):
                raise EnvironmentError('Please run combine opperators before running plot')
            print('plotting ', iFO)
            for istream in iFO.stream_list:
                thiscol = colourset8[icount%len(colourset8)]
                thisshift = shiftset[icount%len(shiftset)]
                data_plot = iFO.PlotMonteTimeDeltaEps(data_plot,itflow,stream=istream,thiscol=thiscol,thisshift=thisshift,thiseps=ieps,MultQ=MultQ)
                icount += 1
        data_plot.PlotAll()
        return data_plot




    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def Plot(self,plot_info,fitDict=defFitDict,AComp='PreDef'):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Susceptibility'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$\sqrt{8t_{f}}[fm]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\chi$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Chit.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        if AComp != 'PreDef': self.AComp = AComp

        fit_is_dict = isinstance(fitDict,(dict,OrderedDict))
        for ic,(ikey,iFO) in enumerate(self.SetFO.items()):
            if isinstance(iFO,list):
                raise EnvironmentError('Please run combine opperators before running plot')
            print('plotting ', iFO)
            if self.AutoList[ic]:
                data_plot = iFO.FlowPlotAuto(data_plot)
            else:
                data_plot = iFO.FlowPlot(data_plot)
                if fit_is_dict:
                    if ikey in list(fitDict.keys()):
                        data_plot = iFO.FlowFitPlot(data_plot,fitDict[ikey])
                else:
                    data_plot = iFO.FlowFitPlot(data_plot,fitDict[ic])

        if self.AComp:
            raise EnvironmentError('Ahmed comparison not implemented with new plotting routines')
            # thiscol = colourset8[ic+1%len(colourset8)]
            # thissym = markerset[ic+1%len(markerset)]
            # thisshift = shiftset[ic+1%len(shiftset)]
            # Adata = ara.AResChi2(Info=self.InfoDict[0],thissym=thissym,thiscol=thiscol,thisshift=thisshift)
            # Adata.ReadChi2()
            # Adata.PlotChi2(xlims=self.SetFO[self.SetFO.keys()[0]].tflowlist)
        data_plot.PlotAll()
        return data_plot




    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def PlotTauIntOverFlow(self,plot_info):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$\tau_{int}(W_{opt})$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$\sqrt{8t_{f}}[fm]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\tau_{int}$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_TauInt_WOpt.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)

        for ic,(ikey,iFO) in enumerate(self.SetFO.items()):
            if isinstance(iFO,list):
                raise EnvironmentError('Please run combine opperators before running plot')
            if self.AutoList[ic]:
                print('plotting Tau Int W optimal ', iFO)
                data_plot = iFO.PlotTauIntWOpt(data_plot)
            else:
                print('Warning, ' , iFO , ' is not autocorrelator analized')
        data_plot.PlotAll()
        return data_plot




    # ## fitlist is list of touples { set : ('state#','fitr#-#')}
    # def PlotRandTimes(self,tflowlist,randt_fitr='None',tsize='Def',xsize='Def',ysize='Def',
    #                   legloc='upper right',ftitle=r'Random Time Dependance',fxlab=r'$nt_{random}$',fylab=r'$\chi$'):
    #     if ftitle=='Def': ftitle=r'Random Time Dependance'
    #     if fxlab=='Def': fxlab=r'$nt_{random}$'
    #     if fylab=='Def': fylab=r'$\chi$'
    #
    #     if self.FullFO == 'Not Set':
    #         raise EnvironmentError('self.FullFO not set, maybe not loaded pickle yet?')
    #     for ic,(itflow,(ikey,iFO)) in enumerate(zip(tflowlist,self.FullFO.iteritems())):
    #         thiscol = colourset8[ic%len(colourset8)]
    #         thissym = markerset[ic%len(markerset)]
    #         thisshift = shiftset[ic%len(shiftset)]
    #         if not self.RandTList[ic]: continue
    #
    #         if isinstance(iFO,list):
    #             raise EnvironmentError('Please run combine opperators before running plot')
    #         print 'plotting ', iFO
    #         if self.AutoList[ic]:
    #             iFO.PlotRandTimesAuto(itflow,thiscol=thiscol,thissym=thissym,thisshift=thisshift)
    #         else:
    #             thisshift = iFO.PlotRandTimes(itflow,thiscol=thiscol,thissym=thissym,thisshift=thisshift)
    #             if isinstance(randt_fitr,basestring):
    #                 if randt_fitr != 'None':
    #                     iFO.PlotFitRT(itflow,randt_fitr,thisshift,thiscol=thiscol)
    #             else:
    #                 iFO.PlotFitRT(itflow,randt_fitr[ic],thisshift,thiscol=thiscol)
    #
    #     pl.legend(loc=legloc)
    #     if tsize == 'Def':
    #         pl.title(ftitle)
    #     else:
    #         pl.title(ftitle,fontsize=tsize)
    #     if xsize == 'Def':
    #         pl.xlabel(fxlab)
    #     else:
    #         pl.xlabel(fxlab,fontsize=xsize)
    #     if ysize == 'Def':
    #         pl.ylabel(fylab)
    #     else:
    #         pl.ylabel(fylab,fontsize=ysize)
    #     # pl.xlim(*defxlim)
    #     pl.savefig(self.Plotdir+self.name+'_RandTimes.pdf')
    #     pl.clf()


    def PlotVsPionMass(self,fit_list,plot_info):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Flowed Ops vs Pion Mass'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$m_{\pi} [MeV]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\chi$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_MpiDep.pdf'
            this_pi['xlims'] = (0,None)
            this_pi['ylims'] = (0,None)
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        if isinstance(fit_list,(tuple,list,np.ndarray)):
            fit_list = fit_list[-1]
        if str(fit_list) == 'None': return
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        vall,errl,indexl = [],[],[]
        for ic,iFO in enumerate(self.SetFO.values()):
            # shiftlist.append(shiftset[ic%len(shiftset)])
            this_pion_mass = (iFO.latparams.GetPionMassLab(Phys=True,MeV=True),)
            if isinstance(iFO,list):
                raise EnvironmentError('Please run combine opperators before running plot')
            this_loop = iFO.Get_Extrapolation(fmted=True)
            if this_loop is None:
                print('Cannot Fit: '+str(iFO.flowname)+', skipping')
                continue
            if ic == 0:
                index_names = ('Pion Mass',) + this_loop.index.names
            for ikey,ival in this_loop.items():
                this_key = (slice(None),ikey[0],fit_list,ikey[-1])
                vall.append(ival.Avg)
                errl.append(ival.Std)
                indexl.append(this_pion_mass+ikey)
        # shiftlist = np.array(shiftlist)*(max(50,max(plotx)-min(plotx)))
        if len(indexl)>0:
            indicies = pa.MultiIndex.from_tuples(indexl,names=index_names)
            ploty = pa.Series(vall,index=indicies)
            plotyerr = pa.Series(errl,index=indicies)
            hold_series = pa.Series()
            hold_series['x_data'] = 'from_keys'
            hold_series['key_select'] = this_key
            hold_series['y_data'] = ploty
            hold_series['yerr_data'] = plotyerr
            hold_series['label'] = 'This Calculation'
            hold_series['shift'] = shiftset[0%len(shiftset)]
            hold_series['color'] = colourset8[0%len(colourset8)]
            hold_series['symbol'] = markerset[0%len(markerset)]
            hold_series['type'] = 'error_bar_vary'
            data_plot.AppendData(hold_series)
        top_susc_0710_1130_plot['color'] = colourset8[1%len(colourset8)]
        top_susc_0710_1130_plot['symbol'] = markerset[1%len(markerset)]
        top_susc_0710_1130_plot['shift'] = shiftset[1%len(shiftset)]
        data_plot.AppendData(top_susc_0710_1130_plot)
        top_susc_1406_5363_plot['color'] = colourset8[2%len(colourset8)]
        top_susc_1406_5363_plot['symbol'] = markerset[2%len(markerset)]
        top_susc_1406_5363_plot['shift'] = shiftset[2%len(shiftset)]
        data_plot.AppendData(top_susc_1406_5363_plot)
        top_susc_1312_5161_plot['color'] = colourset8[3%len(colourset8)]
        top_susc_1312_5161_plot['symbol'] = markerset[3%len(markerset)]
        top_susc_1312_5161_plot['shift'] = shiftset[3%len(shiftset)]
        data_plot.AppendData(top_susc_1312_5161_plot)
        top_susc_1312_5161_211_plot['color'] = colourset8[4%len(colourset8)]
        top_susc_1312_5161_211_plot['symbol'] = markerset[4%len(markerset)]
        top_susc_1312_5161_211_plot['shift'] = shiftset[4%len(shiftset)]
        data_plot.AppendData(top_susc_1312_5161_211_plot)
        top_susc_1512_06746_plot['color'] = colourset8[5%len(colourset8)]
        top_susc_1512_06746_plot['symbol'] = markerset[5%len(markerset)]
        top_susc_1512_06746_plot['shift'] = shiftset[5%len(shiftset)]
        data_plot.AppendData(top_susc_1512_06746_plot)
        top_susc_1705_10906_plot['color'] = colourset8[6%len(colourset8)]
        top_susc_1705_10906_plot['symbol'] = markerset[6%len(markerset)]
        top_susc_1705_10906_plot['shift'] = shiftset[6%len(shiftset)]
        data_plot.AppendData(top_susc_1705_10906_plot)

        data_plot.PlotAll()
        return data_plot

    def PlotVsQuarkMass(self,tflowlist,plot_info):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Flowed Ops vs Up Down Mass'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$m_{ud}^{MSbar} [MeV]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\chi$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_UDmassDep.pdf'
            this_pi['xlims'] = (0,None)
            this_pi['ylims'] = (0,None)
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        ploty,plotyerr,plotx,yboot = [],[],[],[]
        shiftlist,smass = [],[]
        mpi_list = []
        for ic,(itflow,(ikey,iFO)) in enumerate(zip(tflowlist,iter(self.SetFO.items()))):
            shiftlist.append(shiftset[ic%len(shiftset)])
            if isinstance(iFO,list):
                raise EnvironmentError('Please run combine opperators before running plot')
            print('plotting ', iFO)
            plotx.append(iFO.latparams.GetUDMass(Phys=True,MeV=True))
            smass.append(iFO.latparams.GetSMass(Phys=True,MeV=True))
            if iFO.latparams.GetPionMass() not in mpi_list:
                mpi_list.append(iFO.latparams.GetPionMass())
            if self.AutoList[ic]:
                if itflow in list(iFO.Op_Stats['Auto'].keys()):
                    yboot.append('Auto')
                    ploty.append(iFO.Op_Stats['Auto'][itflow].Avg)
                    plotyerr.append(iFO.Op_Stats['Auto'][itflow].Std)
            else:
                if itflow in list(iFO.Op_Stats['boot'].keys()):
                    yboot.append(iFO.Op_Stats['boot'][itflow])
                    ploty.append(iFO.Op_Stats['boot'][itflow].Avg)
                    plotyerr.append(iFO.Op_Stats['boot'][itflow].Std)
        if len(ploty) > 0:
            shiftlist = np.array(shiftlist)*(max(1,max(plotx)-min(plotx)))
            hold_series = pa.Series()
            hold_series['x_data'] = np.array(plotx)
            hold_series['y_data'] = ploty
            hold_series['yerr_data'] = plotyerr
            hold_series['shift'] = shiftlist[0]
            hold_series['symbol'] = 'o'
            hold_series['color'] = 'blue'
            hold_series['label'] = '_'.join(['{:.2f}'.format(ism) for ism in smass])
            hold_series['type'] = 'error_bar'
            hold_series['xdatarange'] = None
            hold_series['otherXvals'] = None
            hold_series['ShowPar'] = None

            data_plot.AppendData(hold_series)
            if len(mpi_list) > 1 and all([iyboot != 'Auto' for iyboot in yboot]):
                fitfun = ff.Fitting(Funs=[Chipt_Chit,1,'Estimate'],data=[np.array([plotx,smass]),yboot],name='$\chi_{pt}$ fit ')
                fitfun.FitBoots()
                hold_series = pa.Series()
                hold_series['xdatarange'] = [0,max(plotx)]
                hold_series['otherXvals'] = [smass[0]]
                hold_series['shift'] = shiftlist[0]
                hold_series['ShowPar'] =r'Sig'
                hold_series['symbol'] = 'o'
                hold_series['color'] = 'blue'
                hold_series['type'] = 'fit'
                hold_series['fit_class'] = fitfun
                data_plot.AppendData(hold_series)
            data_plot.PlotAll()
        return data_plot



    def __setitem__(self, key, value ):
        self.SetFO[key[0]][key[1]] = value

    def __getitem__(self, key):
        return self.SetFO[key[0]][key[1]]

    def Stats(self):
        for iset in self.SetFO.values():
            iset.Stats()

    def __str__(self):
        return self.name

    def __iter__(self):
        return self

    def items(self):
        outdata = []
        for iset,setdata in self.SetFO.items():
            for itflow,tdata in list(setdata.items()):
                outdata.append((iset,itflow,tdata))
        return outdata

    def itemsAvgStd(self):
        outdata = []
        for iset,setdata in self.SetFO.items():
            for itflow,tdata in list(setdata.items()):
                outdata.append((iset,itflow,tdata.Avg,tdata.Std))
        return outdata

    def values(self):
        outdata = []
        for setdata in self.SetFO.values():
            for tdata in list(setdata.values()):
                outdata.append(tdata)
        return outdata



    def valuesAvgStd(self):
        outdata = []
        for setdata in self.SetFO.values():
            for tdata in list(setdata.values()):
                outdata.append((tdata.Avg,tdata.Std))
        return outdata

    def keys(self):
        outdata = []
        for iset,setdata in self.SetFO.items():
            for ip,it in list(setdata.keys()):
                outdata.append((iset,ip,it))
        return outdata

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
        for ikey in list(self.keys()):
            self[ikey] = -self[ikey]
        return self


    def __pos__(self):
        return self

    def __abs__(self):
        for ikey in list(self.keys()):
            self[ikey] = abs(self[ikey])
        return self

    ## Casting, just casts the bootstrap values, ignore return value (only needed since return has to be same type).


    def __complex__(self):
        for ikey in list(self.keys()):
            complex(self[ikey])
            self[ikey].Stats()
        return 1+0j

    def __int__(self):
        for ikey in list(self.keys()):
            int(self[ikey])
            self[ikey].Stats()
        return 1

    def __float__(self):
        for ikey in list(self.keys()):
            float(self[ikey])
            self[ikey].Stats()
        return 1.0

    ## Numerical operatoions

    def __add__(self,FO2):
        if type(FO2) in SetOfFO.parentlist:
            return FO2.__radd__(self)
        else:
            result = SetOfFO(cfglist=self.cfglist, InfoDict=self.InfoDict)
            if isinstance(FO2,SetOfFO):
                result.name = '+'.join([self.name,FO2.name])
                for (iset,itflow,iFO) in list(self.items()):
                    if iset in list(FO2.SetFO.keys()):
                        if (iset,itflow) not in list(FO2.keys()): continue
                        result[(iset,itflow)] = iFO + FO2[(iset,itflow)]
            # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
            #     result.name = '+'.join([self.name,FO2.name])
            #     for (iset,itflow,iFO) in self.items():
            #         if (itflow) not in FO2.keys(): continue
            #         result[(iset,itflow)] = iFO + FO2[(itflow)]
            else:
                try:
                    result.name = '+'.join([self.name,NumbToFileName(FO2)])
                    for iset,itflow,iFO in list(self.items()):
                        result[(iset,itflow)] = iFO + FO2
                except Exception as err:
                    raise EnvironmentError('Invalid value to combine with BootStrap class')
            return result

    def __sub__(self,FO2):
        if type(FO2) in SetOfFO.parentlist:
            return FO2.__rsub__(self)
        else:
            result = SetOfFO( cfglist=self.cfglist, InfoDict=self.InfoDict)
            if isinstance(FO2,SetOfFO):
                result.name = '-'.join([self.name,FO2.name])
                for (iset,itflow,iFO) in list(self.items()):
                    if iset in list(FO2.SetFO.keys()):
                        if (iset,itflow) not in list(FO2.keys()): continue
                        result[(iset,itflow)] = iFO - FO2[(iset,itflow)]
            # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
            #     for (iset,itflow,iFO) in self.items():
            #         result.name = '-'.join([self.name,FO2.name])
            #         if (itflow) not in FO2.keys(): continue
            #         result[(iset,itflow)] = iFO - FO2[(itflow)]
            else:
                try:
                    for iset,itflow,iFO in list(self.items()):
                        result.name = '-'.join([self.name,NumbToFileName(FO2)])
                        result[(iset,itflow)] = iFO - FO2
                except Exception as err:
                    raise EnvironmentError('Invalid value to combine with BootStrap class')
            return result


    def __mul__(self,FO2):
        if type(FO2) in SetOfFO.parentlist:
            return FO2.__rmul__(self)
        else:
            result = SetOfFO( cfglist=self.cfglist, InfoDict=self.InfoDict)
            if isinstance(FO2,SetOfFO):
                result.name = 'x'.join([self.name,FO2.name])
                for (iset,itflow,iFO) in list(self.items()):
                    if iset in list(FO2.SetFO.keys()):
                        if (iset,itflow) not in list(FO2.keys()): continue
                        result[(iset,itflow)] = iFO * FO2[(iset,itflow)]
            # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
            #     for (iset,itflow,iFO) in self.items():
            #         result.name = 'x'.join([self.name,FO2.name])
            #         if (itflow) not in FO2.keys(): continue
            #         result[(iset,itflow)] = iFO * FO2[(itflow)]
            else:
                try:
                    for iset,itflow,iFO in list(self.items()):
                        result.name = 'x'.join([self.name,NumbToFileName(FO2)])
                        result[(iset,itflow)] = iFO * FO2
                except Exception as err:
                    raise EnvironmentError('Invalid value to combine with BootStrap class')
            return result


    def __div__(self,FO2):
        if type(FO2) in SetOfFO.parentlist:
            return FO2.__rdiv__(self)
        else:
            result = SetOfFO( cfglist=self.cfglist, InfoDict=self.InfoDict)
            if isinstance(FO2,SetOfFO):
                result.name = 'div'.join([self.name,FO2.name])
                for (iset,itflow,iFO) in list(self.items()):
                    if iset in list(FO2.SetFO.keys()):
                        try:
                            if (iset,itflow) not in list(FO2.keys()): continue
                            result[(iset,itflow)] = iFO / FO2[(iset,itflow)]
                        except Exception as err:
                            if any([ibs == 0 for ibs in FO2[(iset,itflow)].bootvals]):
                                raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+FO2.name + ' '+FO2[(iset,itflow)].name)
            # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
            #     for (iset,itflow,iFO) in self.items():
            #         result.name = 'div'.join([self.name,FO2.name])
            #         try:
            #             if (itflow) not in FO2.keys(): continue
            #             result[(iset,itflow)] = iFO / FO2[(itflow)]
            #         except Exception as err:
            #             if any([ibs == 0 for ibs in FO2[(itflow)].bootvals]):
            #                 raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+FO2.name + ' '+FO2[(itflow)].name)
            else:
                if FO2 == 0:
                    raise ZeroDivisionError('Dividing by zero found in bootstrap division')
                try:
                    for iset,itflow,iFO in list(self.items()):
                        result.name = 'div'.join([self.name,NumbToFileName(FO2)])
                        result[(iset,itflow)] = iFO / FO2
                except Exception as err:
                    raise EnvironmentError('Invalid value to combine with BootStrap class')
            return result


    def __pow__(self,FO2):
        if type(FO2) in SetOfFO.parentlist:
            return FO2.__rpow__(self)
        else:
            result = SetOfFO( cfglist=self.cfglist, InfoDict=self.InfoDict)
            if isinstance(FO2,SetOfFO):
                result.name = 'pow'.join([self.name,FO2.name])
                for (iset,itflow,iFO) in list(self.items()):
                    if iset in list(FO2.SetFO.keys()):
                        if (iset,itflow) not in list(FO2.keys()): continue
                        result[(iset,itflow)] = iFO ** FO2[(iset,itflow)]
            # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
            #     for (iset,itflow,iFO) in self.items():
            #         result.name = 'pow'.join([self.name,FO2.name])
            #         if (itflow) not in FO2.keys(): continue
            #         result[(iset,itflow)] = iFO ** FO2[(itflow)]
            else:
                try:
                    for iset,itflow,iFO in list(self.items()):
                        result.name = 'pow'.join([self.name,NumbToFileName(FO2)])
                        result[(iset,itflow)] = iFO ** FO2
                except Exception as err:
                    raise EnvironmentError('Invalid value to combine with BootStrap class')
            return result


    ## Right multiplication functions


    def __radd__(self,FO2):
        result = SetOfFO( cfglist=self.cfglist, InfoDict=self.InfoDict)
        if isinstance(FO2,SetOfFO):
            result.name = '+'.join([FO2.name,self.name])
            for (iset,itflow,iFO) in list(self.items()):
                if iset in list(FO2.SetFO.keys()):
                    if (iset,itflow) not in list(FO2.keys()): continue
                    result[(iset,itflow)] = FO2[(iset,itflow)] + iFO
        # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
        #     result.name = '+'.join([FO2.name,self.name])
        #     for (iset,itflow,iFO) in self.items():
        #         if (itflow) not in FO2.keys(): continue
        #         result[(iset,itflow)] = FO2[(itflow)] + iFO
        else:
            try:
                result.name = '+'.join([NumbToFileName(FO2),self.name])
                for iset,itflow,iFO in list(self.items()):
                    result[(iset,itflow)] = FO2 + iFO
            except Exception as err:
                raise EnvironmentError('Invalid value to combine with BootStrap class')
        return result


    def __rsub__(self,FO2):
        result = SetOfFO( cfglist=self.cfglist, InfoDict=self.InfoDict)
        if isinstance(FO2,SetOfFO):
            result.name = '-'.join([FO2.name,self.name])
            for (iset,itflow,iFO) in list(self.items()):
                if iset in list(FO2.SetFO.keys()):
                    if (iset,itflow) not in list(FO2.keys()): continue
                    result[(iset,itflow)] = FO2[(iset,itflow)] - iFO
        # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
        #     result.name = '-'.join([FO2.name,self.name])
        #     for (iset,itflow,iFO) in self.items():
        #         if (itflow) not in FO2.keys(): continue
        #         result[(iset,itflow)] = FO2[(itflow)] - iFO
        else:
            try:
                result.name = '-'.join([NumbToFileName(FO2),self.name])
                for iset,itflow,iFO in list(self.items()):
                    result[(iset,itflow)] = FO2 - iFO
            except Exception as err:
                raise EnvironmentError('Invalid value to combine with BootStrap class')
        return result

    def __rmul__(self,FO2):
        result = SetOfFO( cfglist=self.cfglist, InfoDict=self.InfoDict)
        if isinstance(FO2,SetOfFO):
            result.name = 'x'.join([FO2.name,self.name])
            for (iset,itflow,iFO) in list(self.items()):
                if iset in list(FO2.SetFO.keys()):
                    if (iset,itflow) not in list(FO2.keys()): continue
                    result[(iset,itflow)] = FO2[(iset,itflow)] * iFO
        # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
        #     result.name = 'x'.join([FO2.name,self.name])
        #     for (iset,itflow,iFO) in self.items():
        #         if (itflow) not in FO2.keys(): continue
        #         result[(iset,itflow)] = FO2[(itflow)] * iFO
        else:
            try:
                result.name = 'x'.join([NumbToFileName(FO2),self.name])
                for iset,itflow,iFO in list(self.items()):
                    result[(iset,itflow)] = FO2 * iFO
            except Exception as err:
                raise EnvironmentError('Invalid value to combine with BootStrap class')
        return result


    def __rdiv__(self,FO2):
        result = SetOfFO( cfglist=self.cfglist, InfoDict=self.InfoDict)
        if isinstance(FO2,SetOfFO):
            result.name = 'div'.join([FO2.name,self.name])
            for (iset,itflow,iFO) in list(self.items()):
                if iset in list(FO2.SetFO.keys()):
                    try:
                        if (iset,itflow) not in list(FO2.keys()): continue
                        result[(iset,itflow)] = FO2[(iset,itflow)] / iFO
                    except Exception as err:
                        if any([ibs == 0 for ibs in iFO.bootvals]):
                            raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+self.name + ' '+iFO.name)
        # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
        #     result.name = 'div'.join([FO2.name,self.name])
        #     for (iset,itflow,iFO) in self.items():
        #         try:
        #             if (itflow) not in FO2.keys(): continue
        #             result[(iset,itflow)] = FO2[(itflow)] / iFO
        #         except Exception as err:
        #             if any([ibs == 0 for ibs in iFO.bootvals]):
        #                 raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+FO2.name + ' '+FO2[(itflow)].name)
        else:
            result.name = 'div'.join([NumbToFileName(FO2),self.name])
            try:
                for iset,itflow,iFO in list(self.items()):
                    result[(iset,itflow)] = FO2 / iFO
            except Exception as err:
                raise EnvironmentError('Invalid value to combine with BootStrap class')
        return result


    def __rpow__(self,FO2):
        result = SetOfFO( cfglist=self.cfglist, InfoDict=self.InfoDict)
        if isinstance(FO2,SetOfFO):
            result.name = 'pow'.join([FO2.name,self.name])
            for (iset,itflow,iFO) in list(self.items()):
                if iset in list(FO2.SetFO.keys()):
                    if (iset,itflow) not in list(FO2.keys()): continue
                    result[(iset,itflow)] = FO2[(iset,itflow)] ** iFO
        # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
        #     result.name = 'pow'.join([FO2.name,self.name])
        #     for (iset,itflow,iFO) in self.items():
        #         if (itflow) not in FO2.keys(): continue
        #         result[(iset,itflow)] = FO2[(itflow)] ** iFO
        else:
            try:
                result.name = 'pow'.join([NumbToFileName(FO2),self.name])
                for iset,itflow,iFO in list(self.items()):
                    result[(iset,itflow)] = FO2 ** iFO
            except Exception as err:
                raise EnvironmentError('Invalid value to combine with BootStrap class')
        return result




class SetOfFOFull(object):
    """

    FOFull(set) Inherits from FlowOpFullSquared in FlowOpps.py


    plotting parameters (Move to different file?) ####


    InfoDict is a list of Info dictionaries to be bassed to TwoPtCorr:


    """
    def Construct_FOFull_Set(this_file_list):
        set_list = [FlowOpFullSquared.Construct_FOFull_File(ifile) for ifile in this_file_list]
        out_class = Construct_Empty(SetOfFOFull)
        out_class.cfglist = {}
        out_class.InfoDict = []
        out_class.Plotdir = graphdir+'/FlowOpFull/'
        out_class.name = 'Set1'
        out_class.SetFO = {}
        for iset in set_list:
            out_class.SetFO[iset.name] = iset
        return out_class


    instances = []

    parentlist = []

    def __init__(self, cfglist={}, InfoDict=[],name = '',corr_cfglists=True):
        from FlowOpps import FlowOpFullSquared

        self.cfglist = cfglist
        self.InfoDict = InfoDict
        self.Plotdir = graphdir+'/FlowOpFull/'
        mkdir_p(self.Plotdir)

        # self.FullFO = 'Not Set'
        if len(name) == 0:
            self.name = 'Set'+str(int(len(SetOfFO.instances))+1)
            SetOfFO.instances.append(self.name)
        else:
            self.name = name
            # SetOfFO.instances.append(self.name)
        self.AutoList = []

        self.SetFO = OrderedDict()
        for ic,idict in enumerate(self.InfoDict):
            thiscol = colourset8[ic%len(colourset8)]
            thissym = markerset[ic%len(markerset)]
            thisshift = shiftset[ic%len(shiftset)]

            if 'autocorr' in list(idict.keys()):
                self.AutoList.append(idict['autocorr'])
            else:
                self.AutoList.append(False)

            thisFO = FlowOpFullSquared(cfglist=self.cfglist, Info=idict,thissym=thissym,thiscol=thiscol,thisshift=thisshift)

            if str(thisFO) in list(self.SetFO.keys()):
                FOnext = 2
                while str(thisFO)+'_'+str(FOnext) in list(self.SetFO.keys()):
                    FOnext += 1
                self.SetFO[str(thisFO)+'_'+str(FOnext)] = thisFO
            else:
                self.SetFO[str(thisFO)] = thisFO
            # print 'Setting up classes for ' , thisFO , ' complete'
            # print 'has color',thiscol , 'symbol ' , thissym , 'thisshift',thisshift

        self.SetFOlen = len(self.SetFO)
        if len(list(self.cfglist.keys())) == 0:
            if corr_cfglists:
                self.GetCfgList()
            else:
                self.GetUnCorrCfgList()
        else:
            self.ImportCfgList(self.cfglist)
        ##Implement Reading cfglist if len(cfglist) == 0

        self.current = 0



    def GetUnCorrCfgList(self):
        for iFO in list(self.SetFO.values()):
            iFO.FlowGetCfgList()


    def GetCfgList(self):
        self.GetUnCorrCfgList()
        input_cfglist = []
        for iset in self.SetFO.values():
            if isinstance(iset,(list,tuple,np.ndarray)):
                for jset in iset:
                    input_cfglist.append(jset)
            else:
                input_cfglist.append(iset)
        cfglist = CombineListCfgs(input_cfglist,'Op_cfgs')
        self.ImportCfgList(cfglist)

    def ImportCfgList(self,cfglist):
        self.cfglist = cfglist
        for iFO in self.SetFO.values():
            iFO.FlowImportCfgList(cfglist)

    def Read(self):
        for setdata in self.SetFO.values():
            print('Importing 2pt corr: ' , str(setdata))
            setdata.FlowRead()


    def Write(self):
        for iFO in self.SetFO.values():
            print('Writing 2pt corr: ' , str(iFO))
            iFO.FlowWrite()

    def LoadPickle(self,DefWipe=False,CheckCfgs=False,OnlyQ=False):
        for ic,(setkey,setdata) in enumerate(self.SetFO.items()):
            # print 'LoadPickle flowed operator: ' , str(setdata)
            setdata.FlowLoadPickle(DefWipe=DefWipe,CheckCfgs=CheckCfgs,OnlyQ=OnlyQ)
            self.SetFO[setkey].FlowImportPlotParams(colourset8[ic%len(colourset8)],markerset[ic%len(markerset)]
                                                    ,shiftset[ic%len(shiftset)])
        # self.MakeRandT(DefWipe=DefWipe)
        print()



    def ReadAndWrite(self,DefWipe=False):
        for ic,(setkey,setdata) in enumerate(self.SetFO.items()):
            print('Importing and Writing 2pt corr: ' , str(setdata))
            setdata.FlowReadAndWrite(DefWipe=DefWipe)

    ## returns a TwoPtCorr class object (see TwoPtCorrelators.py)
    def GetSet(self,thisset):
        return self.SetFO[thisset]


    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def PlotTauInt(self,tflowlist,tsinklist,tsumlist,plot_info):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$\tau_{int}$ at Optimal $W$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$W$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\tau_{int}$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_TauInt.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)

        for ic,(itsum,itsink,itflow,(ikey,iFO)) in enumerate(zip(tsumlist,tsinklist,tflowlist,iter(self.SetFO.items()))):
            if self.AutoList[ic]:
                if (itflow,itsink,itsum) not in iFO.Op2_avg_Stats['Auto'].index:
                    print('Warning , flow and sink time ',itflow,itsink,itsum,' not found ')
                else:
                    print('plotting Tau Int ', iFO)
                    data_plot = iFO.Op2_avg_Stats['Auto'][itflow,itsink,itsum].PlotTauInt(data_plot,
                                            thiscol=iFO.thiscol,thissym=iFO.thissym,thisshift=iFO.thisshift)
            else:
                print('Warning, ' , iFO , ' is not autocorrelated analized')
        data_plot.PlotAll()
        return data_plot


    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def PlotWopt(self,tflowlist,tsinklist,tsumlist,plot_info):
        # if ftitle=='Def': ftitle=r'Deturmination of Optimal parameter $W$ at flow time $'+tflowTOLeg(tflow)+'$'
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$\rho $ for corresponding $W$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$W$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\rho$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Wopt.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)

        for ic,(itflow,itsink,itsum,(ikey,iFO)) in enumerate(zip(tflowlist,tsinklist,tsumlist,iter(self.SetFO.items()))):
            if self.AutoList[ic]:
                print('plotting W optimal ', iFO)
                if (itflow,itsink,itsum) not in iFO.Op2_avg_Stats['Auto'].index:
                    print('Warning , flow time and sink time ',itflow,itsink,itsum,' not found ')
                else:
                    data_plot = iFO.Op2_avg_Stats['Auto'][itflow,itsink,itsum].PlotWopt(
                    data_plot,thiscol=iFO.thiscol,thissym=iFO.thissym,thisshift=iFO.thisshift)
            else:
                print('Warning, ' , iFO , ' is not autocorrelated analized')
        data_plot.PlotAll()
        return data_plot




    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def PlotTauIntOverTsink(self,tflowlist,tsumlist,plot_info):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$\tau_{int}(W_{opt})$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'source sink separation $t$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\tau_{int}$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_TauInt_WOpt_Tsink.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)

        for ic,(itflow,itsum,(ikey,iFO)) in enumerate(zip(tflowlist,tsumlist,iter(self.SetFO.items()))):
            if self.AutoList[ic]:
                print('plotting Tau Int W optimal ', iFO, ' over tsink')
                data_plot = iFO.PlotTsinkAuto_TauIntWOpt(data_plot,itflow,itsum)
            else:
                print('Warning, ' , iFO , ' is not autocorrelator analized')
        data_plot.PlotAll()
        return data_plot


    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def PlotOverTsink(self,tflowlist,tsumlist,plot_info,fitrlist='None',shifted=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$Operator Squared$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'source sink separation $t$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'Op Sqrd'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_OpSqrd_Tsink.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        if fitrlist == 'None':
            fitrlist = ['None']*len(tflowlist)

        for ic,(ifitr,itflow,itsum,(ikey,iFO)) in enumerate(zip(fitrlist,tflowlist,tsumlist,iter(self.SetFO.items()))):
            if self.AutoList[ic]:
                print('plotting Autocorrelated ', iFO, ' over tsink')
                data_plot = iFO.PlotOp2CorrTsinkAuto(data_plot,itflow,itsum)
            elif not isinstance(ifitr,str):
                print('plotting Bootstrapped ', iFO, ' over tsink')
                data_plot = iFO.PlotOp2CorrTsink(data_plot,itflow,itsum,fitr=ifitr[0],shifted=shifted)
                if not shifted:
                    iFO.Fit_All_Props(ifitr,WipeFit=False,tflow_list=[itflow],sum_list=[itsum],show_timer=True)
            else:
                print('plotting Bootstrapped ', iFO, ' over tsink')
                data_plot = iFO.PlotOp2CorrTsink(data_plot,itflow,itsum,fitr=ifitr,shifted=shifted)
        data_plot.PlotAll()
        return data_plot


    def PlotOverT(self,tflowlist,plot_info):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$Operator$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'time $t$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'Op'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Op_T.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)

        for ic,(itflow,(ikey,iFO)) in enumerate(zip(tflowlist,iter(self.SetFO.items()))):
            if self.AutoList[ic]:
                print('plotting Autocorrelated ', iFO, ' over t')
                data_plot = iFO.PlotOpTAuto(data_plot,itflow)
            else:
                print('plotting Bootstrapped ', iFO, ' over t')
                data_plot = iFO.PlotOpT(data_plot,itflow)
        data_plot.PlotAll()
        return data_plot


    def PlotVstiFitr(self,tflowlist,tsumlist,fit_max_list,plot_info,fit_min_list='All',parlist='First'):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$Operator Squared$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'min of fit range $t_{i}$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'Op Sqrd'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_OpSqrd_Fit_Min.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)

        if fit_min_list == 'All':
            fit_min_list = ['All']*len(tflowlist)
        if parlist == 'First':
            parlist = ['First']*len(tflowlist)

        for ic,(ifit_max,ifit_min,ipar,itflow,itsum,(ikey,iFO)) in enumerate(zip(fit_max_list,fit_min_list,parlist,
                                                                            tflowlist,tsumlist,iter(self.SetFO.items()))):
            print('plotting Bootstrapped ', iFO, ' over fit min ranges')
            thiscol = colourset8[ic%len(colourset8)]
            thissym = markerset[ic%len(markerset)]
            thisshift = shiftset[ic%len(shiftset)]
            data_plot = iFO.PlotVstiFitr(data_plot,itflow,itsum,ifit_max,
                                        thiscol=thiscol,thissym=thissym,thisshift=thisshift,
                                        fit_min_list=ifit_min,thisparam=ipar)
        data_plot.PlotAll()
        return data_plot


    def PlotVstfFitr(self,tflowlist,tsumlist,fit_min_list,plot_info,fit_max_list='All',parlist='First'):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$Operator Squared$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'max of fit range $t_{e}$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'Op Sqrd'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_OpSqrd_Fit_Max.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        if fit_max_list == 'All':
            fit_max_list = ['All']*len(tflowlist)
        if parlist == 'First':
            parlist = ['First']*len(tflowlist)

        for ic,(ifit_max,ifit_min,ipar,itflow,itsum,(ikey,iFO)) in enumerate(zip(fit_max_list,fit_min_list,parlist,
                                                                            tflowlist,tsumlist,iter(self.SetFO.items()))):
            thiscol = colourset8[ic%len(colourset8)]
            thissym = markerset[ic%len(markerset)]
            thisshift = shiftset[ic%len(shiftset)]
            print('plotting Bootstrapped ', iFO, ' over fit max ranges')
            data_plot = iFO.PlotVstfFitr(data_plot,itflow,itsum,ifit_min,
                                        thiscol=thiscol,thissym=thissym,thisshift=thisshift,
                                        fit_max_list=ifit_max,thisparam=ipar)
        data_plot.PlotAll()
        return data_plot


    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def PlotTauIntOverTsum(self,tflowlist,tsinklist,plot_info):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$\tau_{int}(W_{opt})$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'summed source $t_{sum}$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\tau_{int}$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_TauInt_WOpt_Tsum.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)

        for ic,(itflow,itsink,(ikey,iFO)) in enumerate(zip(tflowlist,tsinklist,iter(self.SetFO.items()))):
            if self.AutoList[ic]:
                print('plotting Tau Int W optimal ', iFO, ' over tsum')
                data_plot = iFO.PlotTsumAuto_TauIntWOpt(data_plot,itflow,itsink)
            else:
                print('Warning, ' , iFO , ' is not autocorrelator analized')
        data_plot.PlotAll()
        return data_plot


    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def PlotOverTsum(self,tflowlist,tsinklist,plot_info):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$Operator Squared$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'summed source $t_{sum}$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'Op Sqrd'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_OpSqrd_Tsum.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)

        for ic,(itflow,itsink,(ikey,iFO)) in enumerate(zip(tflowlist,tsinklist,iter(self.SetFO.items()))):
            if self.AutoList[ic]:
                print('plotting Autocorrelated ', iFO, ' over tsum')
                data_plot = iFO.PlotOp2CorrTsumAuto(data_plot,itflow,itsink)
            else:
                print('plotting Bootstrapped ', iFO, ' over tsum')
                data_plot = iFO.PlotOp2CorrTsum(data_plot,itflow,itsink)
        data_plot.PlotAll()
        return data_plot

    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def PlotTauIntOverTflow(self,tsumlist,tsinklist,plot_info):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$\tau_{int}(W_{opt})$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$\sqrt{8t_{f}}[fm]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\tau_{int}$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_TauInt_WOpt_Tflow.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)

        for ic,(itsum,itsink,(ikey,iFO)) in enumerate(zip(tsumlist,tsinklist,iter(self.SetFO.items()))):
            if self.AutoList[ic]:
                print('plotting Tau Int W optimal ', iFO, ' over tflow')
                data_plot = iFO.PlotTflowAuto_TauIntWOpt(data_plot,itsink,itsum)
            else:
                print('Warning, ' , iFO , ' is not autocorrelator analized')
        data_plot.PlotAll()
        return data_plot


    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def PlotOverTflow(self,tsumlist,tsinklist,plot_info):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Operator Squared'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$\sqrt{8t_{f}}[fm]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'Op Sqrd'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_OpSqrd_Tflow.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)

        for ic,(itsum,itsink,(ikey,iFO)) in enumerate(zip(tsumlist,tsinklist,iter(self.SetFO.items()))):
            if self.AutoList[ic]:
                print('plotting Autocorrelated ', iFO , ' over tflow')
                data_plot = iFO.PlotOp2CorrTflowAuto(data_plot,itsink,itsum)
            else:
                print('plotting Bootstrapped ', iFO, ' over tflow')
                data_plot = iFO.PlotOp2CorrTflow(data_plot,itsink,itsum)
        data_plot.PlotAll()
        return data_plot

    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def PlotFitOverTflow(self,tfitlist,tsumlist,parlist,plot_info):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Operator Squared Fitted'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$\sqrt{8t_{f}}[fm]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'Op Sqrd'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_OpSqrd_Tflow.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)

        counter = 0
        for ic,(ipar,itfit,itsum,(ikey,iFO)) in enumerate(zip(   parlist,tfitlist,tsumlist
                                                                ,iter(self.SetFO.items()))):
            thiscol = colourset8[counter%len(colourset8)]
            thissym = markerset[counter%len(markerset)]
            thisshift = shiftset[counter%len(shiftset)]
            print('plotting fitted results for ', iFO, ' over tflow')
            data_plot = iFO.Plot_Prop_Tflow(data_plot,thiscol=thiscol,thisshift=thisshift,
                                            thissym=thissym,tfitr=itfit,tsum=itsum,thisparam=ipar)
            counter += 1
            if ipar == 'A' or ipar == 'First':
                thiscol = colourset8[counter%len(colourset8)]
                thissym = markerset[counter%len(markerset)]
                thisshift = shiftset[counter%len(shiftset)]
                print('plotting chi for ',iFO)
                data_plot = iFO.PlotChiTflow(   data_plot,thiscol=thiscol,thisshift=thisshift,
                                                thissym=thissym)
                counter += 1
        data_plot.PlotAll()
        return data_plot


    def PlotVsPionMass(self,tflowlist,tfitlist,tsumlist,parlist,plot_info):
        from SetsOfFits import PullSOFSeries_Wpar
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Flowed Ops vs Pion Mass'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$m_{\pi} [MeV]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\chi$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_MpiDep.pdf'
            this_pi['xlims'] = (0,None)
            this_pi['ylims'] = (0,None)
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        vall,indexl = [],[]
        chi_vall,chi_indexl = [],[]
        # shiftlist = []
        for ic,(itflow,ipar,itfit,itsum,(ikey,iFO)) in enumerate(zip(   tflowlist,parlist,
                                                            tfitlist,tsumlist
                                                            ,iter(self.SetFO.items()))):
            if isinstance(itfit,(list,tuple,np.ndarray)):
                itfit = itfit[0]
            if isinstance(iFO,list):
                raise EnvironmentError('Please run combine opperators before running plot')
            print('plotting ', iFO)
            this_pion_mass = (iFO.latparams.GetPionMassLab(Phys=True,MeV=True),)
            if itflow in list(iFO.Chi_Stats['Auto'].keys()):
                chi_indexl.append(this_pion_mass)
                chi_vall.append(iFO.Chi_Stats['Auto'][itflow])
            if not self.AutoList[ic]:
                fit_data = PullSOFSeries_Wpar(iFO.Prop_Fit_Stats.loc[:,'boot'])
                if ipar == 'First':
                    ipar = fit_data.index[0][-1]
                elif isinstance(ipar,int):
                    ipar = fit_data.index.levels[ipar]
                if ic == 0:
                    index_names = ('Pion Mass',) + fit_data.index.names
                if (itsum,itflow,iFO.prop_fit_fun[0].__name__,itfit,ipar) not in fit_data:
                    print(itsum,itflow,iFO.prop_fit_fun.__name__,itfit,ipar, 'Not found, attempting to fit')
                    iFO.Fit_All_Props(itfit,sum_list=[itsum],tflow_list=[itflow],show_timer=True)
                    fit_data = PullSOFSeries_Wpar(iFO.Prop_Fit_Stats.loc[:,'boot'])
                for ikey,ival in fit_data.items():
                    vall.append(ival)
                    indexl.append((this_pion_mass,)+ikey)
        # shiftlist = np.array(shiftlist)*(max(50,max(plotx)-min(plotx)))
        if len(indexl)>0:
            indicies = pa.MultiIndex.from_tuples(indexl,names=index_names)
            plot_data = pa.Series(vall,index=indicies)
            # if all(np.array(parlist) == parlist[0]):
            #     Improv_leg += parlist[0:1]
            hold_series = null_series
            hold_series['x_data'] = 'from_keys'
            hold_series['key_select'] = (slice(None),)+tuple(plot_data.index[0][1:])
            hold_series['y_data'] = plot_data.apply(lambda x : x.Avg)
            hold_series['yerr_data'] = plot_data.apply(lambda x : x.Std)
            hold_series['label'] = 'Improved'
            hold_series['shift'] = shiftset[0%len(shiftset)]
            hold_series['color'] = colourset8[0%len(colourset8)]
            hold_series['symbol'] = markerset[0%len(markerset)]
            hold_series['type'] = 'error_bar_vary'
            data_plot.AppendData(hold_series)
        else:
            print('No improved results to plot vs mpi')
        if len(chi_indexl) > 0:
            indicies = pa.MultiIndex.from_tuples(chi_indexl,names=('Pion Mass',))
            plot_data = pa.Series(chi_vall,index=indicies)
            hold_series = null_series
            hold_series['x_data'] = 'from_keys'
            hold_series['key_select'] = (slice(None),)
            hold_series['y_data'] = plot_data.apply(lambda x : x.Avg)
            hold_series['yerr_data'] = plot_data.apply(lambda x : x.Std)
            hold_series['label'] =  'Integrated'
            hold_series['shift'] = shiftset[1%len(shiftset)]
            hold_series['color'] = colourset8[1%len(colourset8)]
            hold_series['symbol'] = markerset[1%len(markerset)]
            hold_series['type'] = 'error_bar_vary'
            data_plot.AppendData(hold_series)
        else:
            print('No integrated result to plot vs mpi')
        if all(np.array(parlist) == 'First'):
            top_susc_0710_1130_plot['color'] = colourset8[2%len(colourset8)]
            top_susc_0710_1130_plot['symbol'] = markerset[2%len(markerset)]
            top_susc_0710_1130_plot['shift'] = shiftset[2%len(shiftset)]
            data_plot.AppendData(top_susc_0710_1130_plot)
            top_susc_1406_5363_plot['color'] = colourset8[3%len(colourset8)]
            top_susc_1406_5363_plot['symbol'] = markerset[3%len(markerset)]
            top_susc_1406_5363_plot['shift'] = shiftset[3%len(shiftset)]
            data_plot.AppendData(top_susc_1406_5363_plot)
            top_susc_1312_5161_plot['color'] = colourset8[4%len(colourset8)]
            top_susc_1312_5161_plot['symbol'] = markerset[4%len(markerset)]
            top_susc_1312_5161_plot['shift'] = shiftset[4%len(shiftset)]
            data_plot.AppendData(top_susc_1312_5161_plot)
            top_susc_1312_5161_211_plot['color'] = colourset8[5%len(colourset8)]
            top_susc_1312_5161_211_plot['symbol'] = markerset[5%len(markerset)]
            top_susc_1312_5161_211_plot['shift'] = shiftset[5%len(shiftset)]
            data_plot.AppendData(top_susc_1312_5161_211_plot)
            top_susc_1512_06746_plot['color'] = colourset8[6%len(colourset8)]
            top_susc_1512_06746_plot['symbol'] = markerset[6%len(markerset)]
            top_susc_1512_06746_plot['shift'] = shiftset[6%len(shiftset)]
            data_plot.AppendData(top_susc_1512_06746_plot)
            top_susc_1705_10906_plot['color'] = colourset8[7%len(colourset8)]
            top_susc_1705_10906_plot['symbol'] = markerset[7%len(markerset)]
            top_susc_1705_10906_plot['shift'] = shiftset[7%len(shiftset)]
            data_plot.AppendData(top_susc_1705_10906_plot)
        data_plot.PlotAll()
        return data_plot


    def __setitem__(self, key, value ):
        self.SetFO[key[0]][key[1]] = value

    def __getitem__(self, key):
        return self.SetFO[key[0]][key[1]]

    def Stats(self):
        for iset in self.SetFO.values():
            iset.Stats()

    def __str__(self):
        return self.name

    def __iter__(self):
        return self

    def items(self):
        outdata = []
        for iset,setdata in self.SetFO.items():
            for itflow,tdata in list(setdata.items()):
                outdata.append((iset,itflow,tdata))
        return outdata

    def itemsAvgStd(self):
        outdata = []
        for iset,setdata in self.SetFO.items():
            for itflow,tdata in list(setdata.items()):
                outdata.append((iset,itflow,tdata.Avg,tdata.Std))
        return outdata

    def values(self):
        outdata = []
        for setdata in self.SetFO.values():
            for tdata in list(setdata.values()):
                outdata.append(tdata)
        return outdata



    def valuesAvgStd(self):
        outdata = []
        for setdata in self.SetFO.values():
            for tdata in list(setdata.values()):
                outdata.append((tdata.Avg,tdata.Std))
        return outdata

    def keys(self):
        outdata = []
        for iset,setdata in self.SetFO.items():
            for ip,it in list(setdata.keys()):
                outdata.append((iset,ip,it))
        return outdata

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
        for ikey in list(self.keys()):
            self[ikey] = -self[ikey]
        return self


    def __pos__(self):
        return self

    def __abs__(self):
        for ikey in list(self.keys()):
            self[ikey] = abs(self[ikey])
        return self

    ## Casting, just casts the bootstrap values, ignore return value (only needed since return has to be same type).


    def __complex__(self):
        for ikey in list(self.keys()):
            complex(self[ikey])
            self[ikey].Stats()
        return 1+0j

    def __int__(self):
        for ikey in list(self.keys()):
            int(self[ikey])
            self[ikey].Stats()
        return 1

    def __float__(self):
        for ikey in list(self.keys()):
            float(self[ikey])
            self[ikey].Stats()
        return 1.0

    ## Numerical operatoions

    def __add__(self,FO2):
        if type(FO2) in SetOfFO.parentlist:
            return FO2.__radd__(self)
        else:
            result = SetOfFO(cfglist=self.cfglist, InfoDict=self.InfoDict)
            if isinstance(FO2,SetOfFO):
                result.name = '+'.join([self.name,FO2.name])
                for (iset,itflow,iFO) in list(self.items()):
                    if iset in list(FO2.SetFO.keys()):
                        if (iset,itflow) not in list(FO2.keys()): continue
                        result[(iset,itflow)] = iFO + FO2[(iset,itflow)]
            # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
            #     result.name = '+'.join([self.name,FO2.name])
            #     for (iset,itflow,iFO) in self.items():
            #         if (itflow) not in FO2.keys(): continue
            #         result[(iset,itflow)] = iFO + FO2[(itflow)]
            else:
                try:
                    result.name = '+'.join([self.name,NumbToFileName(FO2)])
                    for iset,itflow,iFO in list(self.items()):
                        result[(iset,itflow)] = iFO + FO2
                except Exception as err:
                    raise EnvironmentError('Invalid value to combine with BootStrap class')
            return result

    def __sub__(self,FO2):
        if type(FO2) in SetOfFO.parentlist:
            return FO2.__rsub__(self)
        else:
            result = SetOfFO( cfglist=self.cfglist, InfoDict=self.InfoDict)
            if isinstance(FO2,SetOfFO):
                result.name = '-'.join([self.name,FO2.name])
                for (iset,itflow,iFO) in list(self.items()):
                    if iset in list(FO2.SetFO.keys()):
                        if (iset,itflow) not in list(FO2.keys()): continue
                        result[(iset,itflow)] = iFO - FO2[(iset,itflow)]
            # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
            #     for (iset,itflow,iFO) in self.items():
            #         result.name = '-'.join([self.name,FO2.name])
            #         if (itflow) not in FO2.keys(): continue
            #         result[(iset,itflow)] = iFO - FO2[(itflow)]
            else:
                try:
                    for iset,itflow,iFO in list(self.items()):
                        result.name = '-'.join([self.name,NumbToFileName(FO2)])
                        result[(iset,itflow)] = iFO - FO2
                except Exception as err:
                    raise EnvironmentError('Invalid value to combine with BootStrap class')
            return result


    def __mul__(self,FO2):
        if type(FO2) in SetOfFO.parentlist:
            return FO2.__rmul__(self)
        else:
            result = SetOfFO( cfglist=self.cfglist, InfoDict=self.InfoDict)
            if isinstance(FO2,SetOfFO):
                result.name = 'x'.join([self.name,FO2.name])
                for (iset,itflow,iFO) in list(self.items()):
                    if iset in list(FO2.SetFO.keys()):
                        if (iset,itflow) not in list(FO2.keys()): continue
                        result[(iset,itflow)] = iFO * FO2[(iset,itflow)]
            # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
            #     for (iset,itflow,iFO) in self.items():
            #         result.name = 'x'.join([self.name,FO2.name])
            #         if (itflow) not in FO2.keys(): continue
            #         result[(iset,itflow)] = iFO * FO2[(itflow)]
            else:
                try:
                    for iset,itflow,iFO in list(self.items()):
                        result.name = 'x'.join([self.name,NumbToFileName(FO2)])
                        result[(iset,itflow)] = iFO * FO2
                except Exception as err:
                    raise EnvironmentError('Invalid value to combine with BootStrap class')
            return result


    def __div__(self,FO2):
        if type(FO2) in SetOfFO.parentlist:
            return FO2.__rdiv__(self)
        else:
            result = SetOfFO( cfglist=self.cfglist, InfoDict=self.InfoDict)
            if isinstance(FO2,SetOfFO):
                result.name = 'div'.join([self.name,FO2.name])
                for (iset,itflow,iFO) in list(self.items()):
                    if iset in list(FO2.SetFO.keys()):
                        try:
                            if (iset,itflow) not in list(FO2.keys()): continue
                            result[(iset,itflow)] = iFO / FO2[(iset,itflow)]
                        except Exception as err:
                            if any([ibs == 0 for ibs in FO2[(iset,itflow)].bootvals]):
                                raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+FO2.name + ' '+FO2[(iset,itflow)].name)
            # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
            #     for (iset,itflow,iFO) in self.items():
            #         result.name = 'div'.join([self.name,FO2.name])
            #         try:
            #             if (itflow) not in FO2.keys(): continue
            #             result[(iset,itflow)] = iFO / FO2[(itflow)]
            #         except Exception as err:
            #             if any([ibs == 0 for ibs in FO2[(itflow)].bootvals]):
            #                 raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+FO2.name + ' '+FO2[(itflow)].name)
            else:
                if FO2 == 0:
                    raise ZeroDivisionError('Dividing by zero found in bootstrap division')
                try:
                    for iset,itflow,iFO in list(self.items()):
                        result.name = 'div'.join([self.name,NumbToFileName(FO2)])
                        result[(iset,itflow)] = iFO / FO2
                except Exception as err:
                    raise EnvironmentError('Invalid value to combine with BootStrap class')
            return result


    def __pow__(self,FO2):
        if type(FO2) in SetOfFO.parentlist:
            return FO2.__rpow__(self)
        else:
            result = SetOfFO( cfglist=self.cfglist, InfoDict=self.InfoDict)
            if isinstance(FO2,SetOfFO):
                result.name = 'pow'.join([self.name,FO2.name])
                for (iset,itflow,iFO) in list(self.items()):
                    if iset in list(FO2.SetFO.keys()):
                        if (iset,itflow) not in list(FO2.keys()): continue
                        result[(iset,itflow)] = iFO ** FO2[(iset,itflow)]
            # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
            #     for (iset,itflow,iFO) in self.items():
            #         result.name = 'pow'.join([self.name,FO2.name])
            #         if (itflow) not in FO2.keys(): continue
            #         result[(iset,itflow)] = iFO ** FO2[(itflow)]
            else:
                try:
                    for iset,itflow,iFO in list(self.items()):
                        result.name = 'pow'.join([self.name,NumbToFileName(FO2)])
                        result[(iset,itflow)] = iFO ** FO2
                except Exception as err:
                    raise EnvironmentError('Invalid value to combine with BootStrap class')
            return result


    ## Right multiplication functions


    def __radd__(self,FO2):
        result = SetOfFO( cfglist=self.cfglist, InfoDict=self.InfoDict)
        if isinstance(FO2,SetOfFO):
            result.name = '+'.join([FO2.name,self.name])
            for (iset,itflow,iFO) in list(self.items()):
                if iset in list(FO2.SetFO.keys()):
                    if (iset,itflow) not in list(FO2.keys()): continue
                    result[(iset,itflow)] = FO2[(iset,itflow)] + iFO
        # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
        #     result.name = '+'.join([FO2.name,self.name])
        #     for (iset,itflow,iFO) in self.items():
        #         if (itflow) not in FO2.keys(): continue
        #         result[(iset,itflow)] = FO2[(itflow)] + iFO
        else:
            try:
                result.name = '+'.join([NumbToFileName(FO2),self.name])
                for iset,itflow,iFO in list(self.items()):
                    result[(iset,itflow)] = FO2 + iFO
            except Exception as err:
                raise EnvironmentError('Invalid value to combine with BootStrap class')
        return result


    def __rsub__(self,FO2):
        result = SetOfFO( cfglist=self.cfglist, InfoDict=self.InfoDict)
        if isinstance(FO2,SetOfFO):
            result.name = '-'.join([FO2.name,self.name])
            for (iset,itflow,iFO) in list(self.items()):
                if iset in list(FO2.SetFO.keys()):
                    if (iset,itflow) not in list(FO2.keys()): continue
                    result[(iset,itflow)] = FO2[(iset,itflow)] - iFO
        # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
        #     result.name = '-'.join([FO2.name,self.name])
        #     for (iset,itflow,iFO) in self.items():
        #         if (itflow) not in FO2.keys(): continue
        #         result[(iset,itflow)] = FO2[(itflow)] - iFO
        else:
            try:
                result.name = '-'.join([NumbToFileName(FO2),self.name])
                for iset,itflow,iFO in list(self.items()):
                    result[(iset,itflow)] = FO2 - iFO
            except Exception as err:
                raise EnvironmentError('Invalid value to combine with BootStrap class')
        return result

    def __rmul__(self,FO2):
        result = SetOfFO( cfglist=self.cfglist, InfoDict=self.InfoDict)
        if isinstance(FO2,SetOfFO):
            result.name = 'x'.join([FO2.name,self.name])
            for (iset,itflow,iFO) in list(self.items()):
                if iset in list(FO2.SetFO.keys()):
                    if (iset,itflow) not in list(FO2.keys()): continue
                    result[(iset,itflow)] = FO2[(iset,itflow)] * iFO
        # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
        #     result.name = 'x'.join([FO2.name,self.name])
        #     for (iset,itflow,iFO) in self.items():
        #         if (itflow) not in FO2.keys(): continue
        #         result[(iset,itflow)] = FO2[(itflow)] * iFO
        else:
            try:
                result.name = 'x'.join([NumbToFileName(FO2),self.name])
                for iset,itflow,iFO in list(self.items()):
                    result[(iset,itflow)] = FO2 * iFO
            except Exception as err:
                raise EnvironmentError('Invalid value to combine with BootStrap class')
        return result


    def __rdiv__(self,FO2):
        result = SetOfFO( cfglist=self.cfglist, InfoDict=self.InfoDict)
        if isinstance(FO2,SetOfFO):
            result.name = 'div'.join([FO2.name,self.name])
            for (iset,itflow,iFO) in list(self.items()):
                if iset in list(FO2.SetFO.keys()):
                    try:
                        if (iset,itflow) not in list(FO2.keys()): continue
                        result[(iset,itflow)] = FO2[(iset,itflow)] / iFO
                    except Exception as err:
                        if any([ibs == 0 for ibs in iFO.bootvals]):
                            raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+self.name + ' '+iFO.name)
        # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
        #     result.name = 'div'.join([FO2.name,self.name])
        #     for (iset,itflow,iFO) in self.items():
        #         try:
        #             if (itflow) not in FO2.keys(): continue
        #             result[(iset,itflow)] = FO2[(itflow)] / iFO
        #         except Exception as err:
        #             if any([ibs == 0 for ibs in iFO.bootvals]):
        #                 raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+FO2.name + ' '+FO2[(itflow)].name)
        else:
            result.name = 'div'.join([NumbToFileName(FO2),self.name])
            try:
                for iset,itflow,iFO in list(self.items()):
                    result[(iset,itflow)] = FO2 / iFO
            except Exception as err:
                raise EnvironmentError('Invalid value to combine with BootStrap class')
        return result


    def __rpow__(self,FO2):
        result = SetOfFO( cfglist=self.cfglist, InfoDict=self.InfoDict)
        if isinstance(FO2,SetOfFO):
            result.name = 'pow'.join([FO2.name,self.name])
            for (iset,itflow,iFO) in list(self.items()):
                if iset in list(FO2.SetFO.keys()):
                    if (iset,itflow) not in list(FO2.keys()): continue
                    result[(iset,itflow)] = FO2[(iset,itflow)] ** iFO
        # elif type(FO2) in 'TwoPtCorrelators.TwoPointCorr':
        #     result.name = 'pow'.join([FO2.name,self.name])
        #     for (iset,itflow,iFO) in self.items():
        #         if (itflow) not in FO2.keys(): continue
        #         result[(iset,itflow)] = FO2[(itflow)] ** iFO
        else:
            try:
                result.name = 'pow'.join([NumbToFileName(FO2),self.name])
                for iset,itflow,iFO in list(self.items()):
                    result[(iset,itflow)] = FO2 ** iFO
            except Exception as err:
                raise EnvironmentError('Invalid value to combine with BootStrap class')
        return result
