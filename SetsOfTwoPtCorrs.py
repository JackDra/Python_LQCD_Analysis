#!/usr/bin/env python

import numpy as np

from collections import OrderedDict
# from MiscFuns import NumbToFileName,CreateNameAndLL,GetInfoList,GenSinkFun
from MiscFuns import NumbToFileName
from Params import graphdir,defFitDict,Debug_Mode
from Params import defxlimAlphaTflow
from MomParams import physcial_pion_mass
from MiscFuns import mkdir_p,CombineListCfgs
from PredefFitFuns import SquaredFitFun,LinearFitFun
import FitFunctions as ff
from SetsOfFits import PullSOFSeries_Wpar
from XmlFormatting import MakeValAndErr
# import AhmedRes as ara
import pandas as pa
from TimeStuff import Timer
from copy import copy
import PlotData as jpl


# import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
# import matplotlib.pyplot as pl


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




class SetOfTwoPt(object):
    """

    C2(set) Inherits from TwoPtCorrelators Class

    THIS IS WHERE THE VARIATIONAL METHOD AND SINK/SOURCE SMEARING COMBINATIONS SHOULD GO!!!!

    plotting parameters (Move to different file?) ####


    InfoDict is a list of Info dictionaries to be bassed to TwoPtCorr:
    Setup to do fits taken from dictionaries within InfoDict (see TwoPtCorrelators.py)


    """
    instances = []

    parentlist = []

    def __init__(self, cfglist={}, InfoDict=[],corr_cfglists=True,name = ''):
        from TwoPtCorrelators import TwoPointCorr

        self.cfglist = cfglist
        self.InfoDict = InfoDict
        self.Plotdir = graphdir+'/G2/'
        mkdir_p(self.Plotdir)


        if len(name) == 0:
            self.name = 'Set'+str(int(len(SetOfTwoPt.instances))+1)
            SetOfTwoPt.instances.append(self.name)
        else:
            self.name = name
            # SetOfTwoPt.instances.append(self.name)

        self.VarKeys = []
        self.powlist,self.coefflist = [],[]
        self.SetC2 = OrderedDict()
        for ic,idict in enumerate(self.InfoDict):
            thiscol = colourset8[ic%len(colourset8)]
            thissym = markerset[ic%len(markerset)]
            thisshift = shiftset[ic%len(shiftset)]
            if 'Popp' in list(idict.keys()):
                if idict['Popp'] == 'None':
                    self.powlist.append(1.0)
                else:
                    self.powlist.append(np.float64(idict['Popp']))
            else:
                self.powlist.append(1.0)
            if 'coeffopp' in list(idict.keys()):
                self.coefflist.append(np.float64(idict['coeffopp']))
            else:
                self.powlist.append(1.0)
            thisC2 = TwoPointCorr(cfglist=self.cfglist, Info=idict,thissym=thissym,thiscol=thiscol,thisshift=thisshift)

            if str(thisC2) in list(self.SetC2.keys()):
                C2next = 2
                while str(thisC2)+'_'+str(C2next) in list(self.SetC2.keys()):
                    C2next += 1
                self.SetC2[str(thisC2)+'_'+str(C2next)] = thisC2
            else:
                self.SetC2[str(thisC2)] = thisC2
            # print 'Setting up classes for ' , thisC2 , ' complete'
            # print 'has color',thiscol , 'symbol ' , thissym , 'thisshift',thisshift

        self.SetC2len = len(self.SetC2)
        if not isinstance(self.cfglist,pa.DataFrame):
            if corr_cfglists:
                self.GetCfgList()
            else:
                self.GetUnCorrCfgList()
            ##Implement Reading cfglist if len(cfglist) == 0

        self.current = 0


    ## thisfun(corr1,corr2,corr3,...)
    ## Strongly encoraged to manually set the filename somewhere
    ## e.g. sink smear combine we should get the first name and replace the sink smear with the sink smear label
    ## indicies specifies which sets are used in the combination
    def CombineCorrs(self,thisfun,filename='',LegLab='',jsmlab = '',indicies='All'):
        # for items in self.SetC2.values():
        #     print ''
        #     print items
        #     for iitems in items.itemsAvgStd():
        #         print iitems
        if indicies == 'All':
            indicies = list(range(len(list(self.SetC2.keys()))))
        output2pt = thisfun(*np.array(list(self.SetC2.values()))[indicies])
        if jsmlab != 'Previous': output2pt.jsm = jsmlab
        output2pt.IsComb = True
        output2pt.SetCustomName(string=filename,stringLL=LegLab)
        # output2pt.Bootstrap()
        # print 'DEBUG 2pt'
        # for ival in np.array(self.SetC2.values())[indicies]:
        #     print ival
        #     print ival.C2_Stats['boot']['p000','t13'].Avg,ival.C2_Stats['boot']['p000','t13'].Std
        # print output2pt
        # output2pt.C2_Stats['boot']['p000','t13'].Stats()
        # print 'output',output2pt.C2_Stats['boot']['p000','t13'].Avg,output2pt.C2_Stats['boot']['p000','t13'].Std
        output2pt.Stats()
        output2pt.Write()
        return output2pt

    def AppendSomething(self,newcorr):
        thiscol = colourset8[self.SetC2len%len(colourset8)]
        thissym = markerset[self.SetC2len%len(markerset)]
        thisshift = shiftset[self.SetC2len%len(shiftset)]
        newcorr.ImportPlotParams(thiscol,thissym,thisshift)
        if str(newcorr) in list(self.SetC2.keys()):
            C2next = 2
            while str(newcorr)+'_'+str(C2next) in list(self.SetC2.keys()):
                C2next += 1
            self.SetC2[str(newcorr)+'_'+str(C2next)] = newcorr
        else:
            self.SetC2[str(newcorr)] = newcorr
        self.VarKeys.append(list(self.SetC2.keys())[-1])
        self.SetC2len = len(self.SetC2)


    def AppendCCorrs(self,thisfun,filename='',LegLab='',jsm = 'Previous',indicies='All'):
        thisC2 = self.CombineCorrs(thisfun,filename=filename,LegLab=LegLab,jsmlab = jsm,indicies=indicies)
        self.AppendSomething(thisC2)
        # if str(thisC2) in self.SetC2.keys():
        #     C2next = 2
        #     while str(thisC2)+'_'+str(C2next) in self.SetC2.keys():
        #         C2next += 1
        #     self.SetC2[str(thisC2)+'_'+str(C2next)] = thisC2
        # else:
        #     self.SetC2[str(thisC2)] = thisC2
        # self.VarKeys.append(self.SetC2.keys()[-1])
        # self.SetC2len = len(self.SetC2)

    def GetUnCorrCfgList(self):
        for iC2 in self.SetC2.values():
            iC2.GetCfgList()


    def RemoveFuns(self):
        for iC2 in self.SetC2.values():
            iC2.RemoveFuns()
    def GetFuns(self):
        for iC2 in self.SetC2.values():
            iC2.GetFuns()

    def GetCfgList(self):
        self.GetUnCorrCfgList()
        cfglist = CombineListCfgs(np.array(list(self.SetC2.values())).flatten(),'C2_cfgs')
        self.ImportCfgList(cfglist)


    def ImportCfgList(self,cfglist):
        self.cfglist = cfglist
        for iC2 in self.SetC2.values():
            iC2.ImportCfgList(self.cfglist)

    def Read(self):
        for setdata in self.SetC2.values():
            print('Importing 2pt corr: ' , str(setdata))
            setdata.Read()


    def Write(self):
        for iC2 in self.SetC2.values():
            print('Writing 2pt corr: ' , str(iC2))
            iC2.Write()

    def LoadPickle(self,DefWipe=False,CheckCfgs=False,CheckMom=True,WipeData=True,NoFit=False):
        for ic,setdata in enumerate(self.SetC2.values()):
            print('LoadPickle 2pt corr: ' , str(setdata))
            setdata.LoadPickle(DefWipe=DefWipe,CheckCfgs=CheckCfgs,CheckMom=CheckMom,WipeData=WipeData,NoFit=NoFit)
            setdata.ImportPlotParams(colourset8[ic%len(colourset8)],markerset[ic%len(markerset)]
                                     ,shiftset[ic%len(shiftset)])
        # print

    def ReadAndWrite(self,DefWipe=False):
        for setdata in self.SetC2.values():
            print('Importing and Writing 2pt corr: ' , str(setdata))
            setdata.ReadAndWrite(DefWipe=DefWipe)

    ## returns a TwoPtCorr class object (see TwoPtCorrelators.py)
    def GetSet(self,thisset):
        return self.SetC2[thisset]



    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    ## momlist is list of list, first dimension is for the set, second is for the different momentum (usually plotted at 1 momentum)
    def EffMassPlot(self,plot_info,fitDict=defFitDict,momlist=[],OnlyVar=False,Phys=True):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = 'Effective Mass'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$t$'
            if 'ylabel' not in this_pi:
                if Phys:
                    this_pi['ylabel'] = r'$E_{p} [GeV]$'
                else:
                    this_pi['ylabel'] = r'$aE_{p}$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_EffMass.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        fit_is_dict = isinstance(fitDict,(OrderedDict,dict))
        varcount = 0
        for icount,(ikey,iC2) in enumerate(self.SetC2.items()):
            if (not OnlyVar) or ikey in  self.VarKeys:
                print('Effective mass plotting ', iC2)
                thiscol = colourset8[varcount%len(colourset8)]
                thissym = markerset[varcount%len(markerset)]
                thisshift = shiftset[varcount%len(shiftset)]
                varcount += 1
                if fit_is_dict:
                    thisxlim = [1,int(fitDict[ikey][1].split('-')[-1])]
                    data_plot = iC2.EffMassPlot(data_plot,xlims=thisxlim,momlist=momlist[icount],
                                                thiscol=thiscol,thissym=thissym,
                                                thisshift=thisshift,Phys=Phys)
                    if ikey in list(fitDict.keys()):
                        data_plot = iC2.EffMassFitPlot( data_plot,*fitDict[ikey][:2],momlist=momlist[icount],
                                            thiscol=thiscol,thisshift = thisshift,Phys=Phys)
                else:
                    thisxlim = [1,int(fitDict[icount][1].split('-')[-1])]
                    data_plot = iC2.EffMassPlot(data_plot,xlims=thisxlim,momlist=momlist[icount],
                                                thiscol=thiscol,thissym=thissym,
                                                thisshift=thisshift,Phys=Phys)
                    data_plot = iC2.EffMassFitPlot( data_plot,*fitDict[icount][:2],
                                                    momlist=momlist[icount],thiscol=thiscol,
                                                    thisshift = thisshift,Phys=Phys)
                if iC2.MesOrBar == 'Baryon':
                    hold_series = jpl.null_series
                    # hold_series['x_data'] = iC2.latparams.GetPionMass(MeV=True)
                    # hold_series['xerr_data'] = iC2.latparams.GetPionMassErr(MeV=True)
                    hold_series['y_data'] = iC2.latparams.GetNucleonMass()
                    hold_series['yerr_data'] = iC2.latparams.GetNucleonMassErr()
                    # hold_series['shift'] = shiftset[(icount+1)%len(shiftset)]
                    # hold_series['symbol'] = 'o'
                    hold_series['type'] = 'hline'
                    hold_series['color'] = colourset8[(icount+1)%len(colourset8)]
                    hold_series['label'] = 'Paper m_{\pi}='+MakeValAndErr(iC2.latparams.GetPionMass(MeV=True),
                                                                          iC2.latparams.GetPionMassErr(MeV=True))+'MeV'
                    data_plot.AppendData(hold_series)
        data_plot.PlotAll()
        return data_plot

    ## assumes first correlator is the master correlator
    ## and second correlator is the excited state subtractor
    def ExcitedSubtraction(self,mass_ratio_guess,fit_params):
        key1,key2 = list(self.SetC2.keys())[0],list(self.SetC2.keys())[1]
        # self.SetC2[key1].Fit()
        # self.SetC2[key2].Fit()

        ## TODO, pull norms out from fit within self.SetC2[key].FitDict.Params[0] or something...
        # print 'debug'
        # print self.SetC2[key1].FitDict.keys()
        # print self.SetC2[key1].FitDict['state2']['p000'].keys()
        # massfitr1 =  self.SetC2[key1].FitDict['state2']['p000'].keys()[0]
        massfitr2 =  list(self.SetC2[key2].FitDict['state2']['p000'].keys())[0]
        #bootstrapped?
        # massfit2 = self.SetC2[key2].FitDict['state2']['p000'][massfitr2].Params['Energy']
        # Zp = self.SetC2[key1].FitDict['state2']['p000'][massfitr2].Params['A_Ep']
        # Zpi = self.SetC2[key2].FitDict['state2']['p000'][massfitr1].Params['A_Ep']
        massfit2 = self.SetC2[key2].FitDict['state2']['p000'][massfitr2].fit_data['Params']['Energy']
        # Zp = self.SetC2[key1].FitDict['state2']['p000'][massfitr2].Params['A_Ep']
        # Zpi = self.SetC2[key2].FitDict['state2']['p000'][massfitr1].Params['A_Ep']
        # print 'debug'
        # print massfit2.Avg
        # print Zp.Avg
        # print Zpi.Avg
        # Gpnorm = self.SetC2[key1]/Zp
        # Gpinorm = self.SetC2[key2]/Zpi
        Gpnorm = self.SetC2[key1]
        Gpinorm = self.SetC2[key2]
        # SubCorr = Gpnorm - (Gpnorm*(Gpinorm**(-1))) + 1
        # coefflist = np.array([np.exp(-massfit2*it) for it in range(len(self.SetC2[key1].keys()))])
        coefflist = np.array([(-massfit2*it*mass_ratio_guess).Exp() for it in range(len(list(self.SetC2[key1].keys())))])
        # SubCorr = Gpnorm/2. - Gpinorm/2. + np.array(coefflist)
        SubCorr = (Gpnorm*coefflist)/(Gpinorm**mass_ratio_guess)
        thisfilename = '_'.join([self.SetC2[key1].filename,'div',self.SetC2[key2].filename])
        thisLL = '_'.join([self.SetC2[key1].LegLab,'div',self.SetC2[key2].LegLab])
        # thisLL = 'Ratio'
        SubCorr.SetCustomName(thisfilename,thisLL)
        SubCorr.Stats()
        if len(fit_params) == 0: fit_params.append('state1')
        if len(fit_params) == 1: fit_params.append('fitr4-20')
        if len(fit_params) == 2: fit_params.append([1/77000.,50,0.5,0.03])
        SubCorr.Fit(state=fit_params[0],fit_range=fit_params[1],iGuess=fit_params[2])
        # print SubCorr.FitDict.keys(),'state2'
        # print SubCorr.FitDict['state2'].keys(),'p000'
        # print SubCorr.FitDict['state2']['p000'].keys(), 'fitr5-20'
        self.AppendSomething(SubCorr)

    def Reformat2pt(self,thisdir):
        print('reformatting 2 point correlation functions into ' , thisdir)
        for iC2 in self.SetC2.values():
            print('reformatting ',str(iC2))
            iC2.Read()
            iC2.SetCustomName()
            iC2.WriteCfgsToFile(thisdir)


    def WriteBoot2pt(self,thisdir):
        print('writing boot of 2 point correlation functions into ' , thisdir)
        for iC2 in self.SetC2.values():
            print('writing boot ',str(iC2))
            iC2.Read()
            iC2.Bootstrap()
            iC2.SetCustomName()
            iC2.WriteBootToFile(thisdir)

    def DoPowers(self):
        for ipow,(ikey,iC2) in zip(self.powlist,iter(self.SetC2.items())):
            self.SetC2[ikey] = iC2**ipow
            if float(ipow) != 1.0:
                print('Doing power of',ipow)
                thisfilename =  self.SetC2[ikey].filename +'_pow'+'{:.1f}'.format(ipow)
                thisLegLab =  '$('+self.SetC2[ikey].LegLab.replace('$','') +')^{'+str(int(ipow))+'}$'
                self.SetC2[ikey].SetCustomName(string=thisfilename,stringLL=thisLegLab)
            self.SetC2[ikey].Stats()

    def DoCoeff(self):
        for icoeff,(ikey,iC2) in zip(self.coefflist,iter(self.SetC2.items())):
            self.SetC2[ikey] = iC2*icoeff
            if float(icoeff) != 1.0:
                print('Doing coefficient of',icoeff)
                # ## REMOVE THIS, TODO, DEBUG !!
                # if '32' in self.SetC2[ikey].filename:
                #     print 'debug test'
                #     print self.SetC2[ikey].filename
                #     self.SetC2[ikey] = self.SetC2[ikey] + 1
                thisfilename =  '_{:.1f}x'.format(icoeff)+self.SetC2[ikey].filename
                thisLegLab =  '$'+'{:.2f}'.format(icoeff)+'('+self.SetC2[ikey].LegLab.replace('$','') +')$'
                self.SetC2[ikey].SetCustomName(string=thisfilename,stringLL=thisLegLab)
            self.SetC2[ikey].Stats()


    def LogPlot(self,plot_info,momlist=[],norm=True,OnlyVar=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Log of $G_{2}$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$t$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$log(G_{2})$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Log.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        varcount = 0
        for icount,(ikey,iC2) in enumerate(self.SetC2.items()):
            if (not OnlyVar) or ikey in  self.VarKeys:
                print('Log plotting ', iC2)
                thiscol = colourset8[varcount%len(colourset8)]
                thissym = markerset[varcount%len(markerset)]
                thisshift = shiftset[varcount%len(shiftset)]
                varcount += 1
                data_plot = iC2.LogPlot(data_plot,momlist=momlist[icount],norm=norm,
                                        thiscol=thiscol,thissym=thissym,
                                        thisshift=thisshift)
            ## TODO, implement fit plots, just do not have time for it now
        data_plot.PlotAll()
        return data_plot


    def Plot(self,plot_info,momlist=[],norm=True,OnlyVar=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$G_{2}$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$t$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$G_{2}$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Corr.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        varcount = 0
        for icount,(ikey,iC2) in enumerate(self.SetC2.items()):
            if (not OnlyVar) or ikey in  self.VarKeys:
                print('Correlator plotting ', iC2)
                thiscol = colourset8[varcount%len(colourset8)]
                thissym = markerset[varcount%len(markerset)]
                thisshift = shiftset[varcount%len(shiftset)]
                varcount += 1
                data_plot = iC2.Plot(   data_plot,momlist=momlist[icount],norm=norm,
                                        thiscol=thiscol,thissym=thissym,
                                        thisshift=thisshift)
            ## TODO, implement fit plots, just do not have time for it now
        data_plot.PlotAll()
        return data_plot



    ## fitparams = [ (statefit1,imom1,fitr1) , ....]
    def PlotMassVsPionMass(self,fitparams,plot_info,OnlyVar=False,Phys=True,par_type='Energy'):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = 'Nucleon vs Pion Mass'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$m_{\pi} [MeV]$'
            if 'ylabel' not in this_pi:
                if Phys:
                    this_pi['ylabel'] = r'$E_{p} [GeV]$'
                else:
                    this_pi['ylabel'] = r'$aE_{p}$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_MpiDep.pdf'
            if 'xlims' in this_pi:
                this_pi['xlims'] = [0,this_pi['xlims'][1]]
            else:
                this_pi['xlims'] = [0,None]
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        varcount = 0
        if isinstance(par_type,str):
            par_type = [par_type]*len(list(self.SetC2.keys()))
        vall,indexl = [],[]
        xdataPap,ydataPap = [],[]
        xdataPaperr,ydataPaperr = [],[]
        shiftlist = []
        for icount,((ikey,iC2),(ifitstate,imom,ifitr)) in enumerate(zip(iter(self.SetC2.items()),fitparams)):
            # print 'DEBUG:', ikey, ifitstate, imom, ifitr
            if (not OnlyVar) or ikey in  self.VarKeys:
                print('pion mass dependance plotting ', iC2)
                shiftlist.append(shiftset[varcount%len(shiftset)])
                varcount += 1
                this_pion_mass = iC2.latparams.GetPionMass(Phys=True,MeV=True)
                fit_data = PullSOFSeries_Wpar(iC2.C2_Fit_Stats.loc[:,'boot'],fmted=True)
                ipar = par_type[icount]
                if ipar == 'First':
                    ipar = fit_data.index[0][-1]
                elif isinstance(ipar,int):
                    ipar = fit_data.index.levels[ipar]
                fit_fun = fit_data.index[0][2]
                if icount == 0:
                    index_names = ('Pion Mass',) + fit_data.index.names
                imom = iC2.latparams.TOpform(imom,'p')
                if (ifitstate,imom,fit_fun,ifitr,ipar) not in fit_data:
                    print(ifitstate,imom,fit_fun,ifitr,ipar, 'Not found, attempting to fit')
                    iC2.Fit(state=ifitstate,fit_range=ifitr)
                    fit_data = PullSOFSeries_Wpar(iC2.C2_Fit_Stats.loc[:,'boot'],fmted=True)
                for ikey,ival in fit_data.items():
                    vall.append(ival)
                    indexl.append((this_pion_mass,)+ikey)
                xdataPap.append(iC2.latparams.GetPionMass(MeV=True))
                xdataPaperr.append(iC2.latparams.GetPionMassErr(MeV=True))
                ydataPap.append(iC2.latparams.GetNucleonMass())
                ydataPaperr.append(iC2.latparams.GetNucleonMassErr())
        if len(indexl) > 0:
            indicies = pa.MultiIndex.from_tuples(indexl,names=index_names)
            plot_data = pa.Series(vall,index=indicies)
            hold_series = jpl.null_series
            hold_series['x_data'] = 'from_keys'
            hold_series['color'] = colourset8[0%len(colourset8)]
            hold_series['y_data'] = plot_data.apply(lambda x : x.Avg)
            hold_series['yerr_data'] = plot_data.apply(lambda x : x.Std)
            hold_series['shift'] = shiftlist[0]
            hold_series['symbol'] = 'o'
            hold_series['key_select'] = (slice(None),) + plot_data.index[0][1:]
            hold_series['type'] = 'error_bar_vary'
            hold_series['label'] = 'This Work'
            data_plot.AppendData(hold_series)
            hold_series = jpl.null_series
            hold_series['x_data'] = np.array(xdataPap)
            hold_series['xerr_data'] = np.array(xdataPaperr)
            hold_series['y_data'] = ydataPap
            hold_series['yerr_data'] = ydataPaperr
            hold_series['shift'] = shiftlist[0]*1.05
            hold_series['symbol'] = 'o'
            hold_series['type'] = 'error_bar'
            hold_series['color'] = colourset8[1%len(colourset8)]
            hold_series['label'] = 'Paper'
            data_plot.AppendData(hold_series)
            data_plot.PlotAll()
        return data_plot



    ## TODO make Nucleon Energy vs Recurrence relation energy along with fit.
    ## fitparams = [ (statefit1,imom1,fitr1) , ....]
    def PlotEnergyVSRecRel(self,fitparams,plot_info,tsize='Def',xsize='Def',ysize='Def',
                           legloc='upper right',ftitle='Lattice Energy vs Recurrence Relation',fxlab=r'$E_{p,Rec} [GeV]$',
                           fylab=r'$aM_{N}$',OnlyVar=False,Phys=True,xmore=[0,0],ymore=[0,0]):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = 'Lattice Energy vs Recurrence Relation'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$E_{p,Rec} [GeV]$'
            if 'ylabel' not in this_pi:
                if Phys:
                    this_pi['ylabel'] = r'$E_{p} [GeV]$'
                else:
                    this_pi['ylabel'] = r'$aE_{p}$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_RecRel.pdf'
            # if 'xlims' in this_pi:
            #     this_pi['xlims'] = [0,this_pi['xlims'][1]]
            # else:
            #     this_pi['xlims'] = [0,None]
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        varcount = 0
        for icount,((istate,imom,ifitr),(ikey,iC2)) in enumerate(zip(fitparams,iter(self.SetC2.items()))):
            if (not OnlyVar) or ikey in  self.VarKeys:
                print('Recurrence relation plotting ', iC2)
                thiscol = colourset8[varcount%len(colourset8)]
                thissym = markerset[varcount%len(markerset)]
                thisshift = shiftset[varcount%len(shiftset)]
                varcount += 1
                data_plot = iC2.RecRelPlot(data_plot,istate,ifitr,
                                                        thiscol=thiscol,thissym=thissym,
                                                        thisshift=thisshift,Phys=Phys)
        hold_series = pa.Series()
        hold_series['x_data'] = [0,5]
        hold_series['y_data'] = [0,5]
        hold_series['yerr_data'] = None
        hold_series['shift'] = None
        hold_series['symbol'] = None
        hold_series['type'] = 'plot'
        hold_series['label'] = 'Linear'
        data_plot.AppendData(hold_series)
        data_plot.PlotAll()
        return data_plot




    def __setitem__(self, key, value ):
        self.SetC2[key[0]][(key[1],key[2])] = value

    def __getitem__(self, key):
        return self.SetC2[key[0]][(key[1],key[2])]

    def Stats(self):
        for iset in self.SetC2.values():
            iset.Stats()

    def __str__(self):
        return self.name

    def __iter__(self):
        return self

    def items(self):
        outdata = []
        for iset,setdata in self.SetC2.items():
            for ip,it,tdata in list(setdata.items()):
                outdata.append((iset,ip,it,tdata))
        return outdata

    def itemsAvgStd(self):
        outdata = []
        for iset,setdata in self.SetC2.items():
            for ip,it,tdata in list(setdata.items()):
                outdata.append((iset,ip,it,tdata.Avg,tdata.Std))
        return outdata

    def values(self):
        outdata = []
        for setdata in self.SetC2.values():
            for tdata in list(setdata.values()):
                outdata.append(tdata)
        return outdata



    def valuesAvgStd(self):
        outdata = []
        for setdata in self.SetC2.values():
            for tdata in list(setdata.values()):
                outdata.append((tdata.Avg,tdata.Std))
        return outdata

    def keys(self):
        outdata = []
        for iset,setdata in self.SetC2.items():
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

    def __add__(self,C22):
        if type(C22) in SetOfTwoPt.parentlist:
            return C22.__radd__(self)
        else:
            result = SetOfTwoPt(cfglist=self.cfglist, InfoDict=self.InfoDict)
            if isinstance(C22,SetOfTwoPt):
                result.name = '+'.join([self.name,C22.name])
                for (iset,ip,it,iC2) in list(self.items()):
                    if iset in list(C22.SetC2.keys()):
                        if (iset,ip,it) not in list(C22.keys()): continue
                        result[(iset,ip,it)] = iC2 + C22[(iset,ip,it)]
            elif type(C22) in 'TwoPtCorrelators.TwoPointCorr':
                result.name = '+'.join([self.name,C22.name])
                for (iset,ip,it,iC2) in list(self.items()):
                    if (ip,it) not in list(C22.keys()): continue
                    result[(iset,ip,it)] = iC2 + C22[(ip,it)]
            else:
                try:
                    result.name = '+'.join([self.name,NumbToFileName(C22)])
                    for iset,ip,it,iC2 in list(self.items()):
                        result[(iset,ip,it)] = iC2 + C22
                except:
                    raise EnvironmentError('Invalid value to combine with BootStrap class')
            return result

    def __sub__(self,C22):
        if type(C22) in SetOfTwoPt.parentlist:
            return C22.__rsub__(self)
        else:
            result = SetOfTwoPt( cfglist=self.cfglist, InfoDict=self.InfoDict)
            if isinstance(C22,SetOfTwoPt):
                result.name = '-'.join([self.name,C22.name])
                for (iset,ip,it,iC2) in list(self.items()):
                    if iset in list(C22.SetC2.keys()):
                        if (iset,ip,it) not in list(C22.keys()): continue
                        result[(iset,ip,it)] = iC2 - C22[(iset,ip,it)]
            elif type(C22) in 'TwoPtCorrelators.TwoPointCorr':
                for (iset,ip,it,iC2) in list(self.items()):
                    result.name = '-'.join([self.name,C22.name])
                    if (ip,it) not in list(C22.keys()): continue
                    result[(iset,ip,it)] = iC2 - C22[(ip,it)]
            else:
                try:
                    for iset,ip,it,iC2 in list(self.items()):
                        result.name = '-'.join([self.name,NumbToFileName(C22)])
                        result[(iset,ip,it)] = iC2 - C22
                except:
                    raise EnvironmentError('Invalid value to combine with BootStrap class')
            return result


    def __mul__(self,C22):
        if type(C22) in SetOfTwoPt.parentlist:
            return C22.__rmul__(self)
        else:
            result = SetOfTwoPt( cfglist=self.cfglist, InfoDict=self.InfoDict)
            if isinstance(C22,SetOfTwoPt):
                result.name = 'x'.join([self.name,C22.name])
                for (iset,ip,it,iC2) in list(self.items()):
                    if iset in list(C22.SetC2.keys()):
                        if (iset,ip,it) not in list(C22.keys()): continue
                        result[(iset,ip,it)] = iC2 * C22[(iset,ip,it)]
            elif type(C22) in 'TwoPtCorrelators.TwoPointCorr':
                for (iset,ip,it,iC2) in list(self.items()):
                    result.name = 'x'.join([self.name,C22.name])
                    if (ip,it) not in list(C22.keys()): continue
                    result[(iset,ip,it)] = iC2 * C22[(ip,it)]
            else:
                try:
                    for iset,ip,it,iC2 in list(self.items()):
                        result.name = 'x'.join([self.name,NumbToFileName(C22)])
                        result[(iset,ip,it)] = iC2 * C22
                except:
                    raise EnvironmentError('Invalid value to combine with BootStrap class')
            return result


    def __div__(self,C22):
        if type(C22) in SetOfTwoPt.parentlist:
            return C22.__rdiv__(self)
        else:
            result = SetOfTwoPt( cfglist=self.cfglist, InfoDict=self.InfoDict)
            if isinstance(C22,SetOfTwoPt):
                result.name = 'div'.join([self.name,C22.name])
                for (iset,ip,it,iC2) in list(self.items()):
                    if iset in list(C22.SetC2.keys()):
                        try:
                            if (iset,ip,it) not in list(C22.keys()): continue
                            result[(iset,ip,it)] = iC2 / C22[(iset,ip,it)]
                        except:
                            if any([ibs == 0 for ibs in C22[(iset,ip,it)].bootvals]):
                                raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+C22.name + ' '+C22[(iset,ip,it)].name)
            elif type(C22) in 'TwoPtCorrelators.TwoPointCorr':
                for (iset,ip,it,iC2) in list(self.items()):
                    result.name = 'div'.join([self.name,C22.name])
                    try:
                        if (ip,it) not in list(C22.keys()): continue
                        result[(iset,ip,it)] = iC2 / C22[(ip,it)]
                    except:
                        if any([ibs == 0 for ibs in C22[(ip,it)].bootvals]):
                            raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+C22.name + ' '+C22[(ip,it)].name)
            else:
                if C22 == 0:
                    raise ZeroDivisionError('Dividing by zero found in bootstrap division')
                try:
                    for iset,ip,it,iC2 in list(self.items()):
                        result.name = 'div'.join([self.name,NumbToFileName(C22)])
                        result[(iset,ip,it)] = iC2 / C22
                except:
                    raise EnvironmentError('Invalid value to combine with BootStrap class')
            return result


    def __pow__(self,C22):
        if type(C22) in SetOfTwoPt.parentlist:
            return C22.__rpow__(self)
        else:
            result = SetOfTwoPt( cfglist=self.cfglist, InfoDict=self.InfoDict)
            if isinstance(C22,SetOfTwoPt):
                result.name = 'pow'.join([self.name,C22.name])
                for (iset,ip,it,iC2) in list(self.items()):
                    if iset in list(C22.SetC2.keys()):
                        if (iset,ip,it) not in list(C22.keys()): continue
                        result[(iset,ip,it)] = iC2 ** C22[(iset,ip,it)]
            elif type(C22) in 'TwoPtCorrelators.TwoPointCorr':
                for (iset,ip,it,iC2) in list(self.items()):
                    result.name = 'pow'.join([self.name,C22.name])
                    if (ip,it) not in list(C22.keys()): continue
                    result[(iset,ip,it)] = iC2 ** C22[(ip,it)]
            else:
                try:
                    for iset,ip,it,iC2 in list(self.items()):
                        result.name = 'pow'.join([self.name,NumbToFileName(C22)])
                        result[(iset,ip,it)] = iC2 ** C22
                except:
                    raise EnvironmentError('Invalid value to combine with BootStrap class')
            return result


    ## Right multiplication functions


    def __radd__(self,C22):
        result = SetOfTwoPt( cfglist=self.cfglist, InfoDict=self.InfoDict)
        if isinstance(C22,SetOfTwoPt):
            result.name = '+'.join([C22.name,self.name])
            for (iset,ip,it,iC2) in list(self.items()):
                if iset in list(C22.SetC2.keys()):
                    if (iset,ip,it) not in list(C22.keys()): continue
                    result[(iset,ip,it)] = C22[(iset,ip,it)] + iC2
        elif type(C22) in 'TwoPtCorrelators.TwoPointCorr':
            result.name = '+'.join([C22.name,self.name])
            for (iset,ip,it,iC2) in list(self.items()):
                if (ip,it) not in list(C22.keys()): continue
                result[(iset,ip,it)] = C22[(ip,it)] + iC2
        else:
            try:
                result.name = '+'.join([NumbToFileName(C22),self.name])
                for iset,ip,it,iC2 in list(self.items()):
                    result[(iset,ip,it)] = C22 + iC2
            except:
                raise EnvironmentError('Invalid value to combine with BootStrap class')
        return result


    def __rsub__(self,C22):
        result = SetOfTwoPt( cfglist=self.cfglist, InfoDict=self.InfoDict)
        if isinstance(C22,SetOfTwoPt):
            result.name = '-'.join([C22.name,self.name])
            for (iset,ip,it,iC2) in list(self.items()):
                if iset in list(C22.SetC2.keys()):
                    if (iset,ip,it) not in list(C22.keys()): continue
                    result[(iset,ip,it)] = C22[(iset,ip,it)] - iC2
        elif type(C22) in 'TwoPtCorrelators.TwoPointCorr':
            result.name = '-'.join([C22.name,self.name])
            for (iset,ip,it,iC2) in list(self.items()):
                if (ip,it) not in list(C22.keys()): continue
                result[(iset,ip,it)] = C22[(ip,it)] - iC2
        else:
            try:
                result.name = '-'.join([NumbToFileName(C22),self.name])
                for iset,ip,it,iC2 in list(self.items()):
                    result[(iset,ip,it)] = C22 - iC2
            except:
                raise EnvironmentError('Invalid value to combine with BootStrap class')
        return result

    def __rmul__(self,C22):
        result = SetOfTwoPt( cfglist=self.cfglist, InfoDict=self.InfoDict)
        if isinstance(C22,SetOfTwoPt):
            result.name = 'x'.join([C22.name,self.name])
            for (iset,ip,it,iC2) in list(self.items()):
                if iset in list(C22.SetC2.keys()):
                    if (iset,ip,it) not in list(C22.keys()): continue
                    result[(iset,ip,it)] = C22[(iset,ip,it)] * iC2
        elif type(C22) in 'TwoPtCorrelators.TwoPointCorr':
            result.name = 'x'.join([C22.name,self.name])
            for (iset,ip,it,iC2) in list(self.items()):
                if (ip,it) not in list(C22.keys()): continue
                result[(iset,ip,it)] = C22[(ip,it)] * iC2
        else:
            try:
                result.name = 'x'.join([NumbToFileName(C22),self.name])
                for iset,ip,it,iC2 in list(self.items()):
                    result[(iset,ip,it)] = C22 * iC2
            except:
                raise EnvironmentError('Invalid value to combine with BootStrap class')
        return result


    def __rdiv__(self,C22):
        result = SetOfTwoPt( cfglist=self.cfglist, InfoDict=self.InfoDict)
        if isinstance(C22,SetOfTwoPt):
            result.name = 'div'.join([C22.name,self.name])
            for (iset,ip,it,iC2) in list(self.items()):
                if iset in list(C22.SetC2.keys()):
                    try:
                        if (iset,ip,it) not in list(C22.keys()): continue
                        result[(iset,ip,it)] = C22[(iset,ip,it)] / iC2
                    except:
                        if any([ibs == 0 for ibs in iC2.bootvals]):
                            raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+self.name + ' '+iC2.name)
        elif type(C22) in 'TwoPtCorrelators.TwoPointCorr':
            result.name = 'div'.join([C22.name,self.name])
            for (iset,ip,it,iC2) in list(self.items()):
                try:
                    if (ip,it) not in list(C22.keys()): continue
                    result[(iset,ip,it)] = C22[(ip,it)] / iC2
                except:
                    if any([ibs == 0 for ibs in iC2.bootvals]):
                        raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+C22.name + ' '+C22[(ip,it)].name)
        else:
            result.name = 'div'.join([NumbToFileName(C22),self.name])
            try:
                for iset,ip,it,iC2 in list(self.items()):
                    result[(iset,ip,it)] = C22 / iC2
            except:
                raise EnvironmentError('Invalid value to combine with BootStrap class')
        return result


    def __rpow__(self,C22):
        result = SetOfTwoPt( cfglist=self.cfglist, InfoDict=self.InfoDict)
        if isinstance(C22,SetOfTwoPt):
            result.name = 'pow'.join([C22.name,self.name])
            for (iset,ip,it,iC2) in list(self.items()):
                if iset in list(C22.SetC2.keys()):
                    if (iset,ip,it) not in list(C22.keys()): continue
                    result[(iset,ip,it)] = C22[(iset,ip,it)] ** iC2
        elif type(C22) in 'TwoPtCorrelators.TwoPointCorr':
            result.name = 'pow'.join([C22.name,self.name])
            for (iset,ip,it,iC2) in list(self.items()):
                if (ip,it) not in list(C22.keys()): continue
                result[(iset,ip,it)] = C22[(ip,it)] ** iC2
        else:
            try:
                result.name = 'pow'.join([NumbToFileName(C22),self.name])
                for iset,ip,it,iC2 in list(self.items()):
                    result[(iset,ip,it)] = C22 ** iC2
            except:
                raise EnvironmentError('Invalid value to combine with BootStrap class')
        return result







##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################



class SetOfNNQ(object):

    """

    NNQ(set) Inherits from TwoPtCorrelators Class

    """
    instances = []

    parentlist = []


    ## InfoDict is a list of Info dictionaries to be bassed to TwoPtCorr:
    ## Setup to do fits taken from dictionaries within InfoDict (see TwoPtCorrelators.py)
    def __init__(self, cfglist={}, InfoDict=[],corr_cfglists=True,name = '',CompAhmed=False):
        from TwoPtCorrelators import NNQCorr

        self.cfglist = cfglist
        self.InfoDict = InfoDict
        self.Plotdir = graphdir+'/Alpha/'
        mkdir_p(self.Plotdir)

        self.AComp = CompAhmed

        if len(name) == 0:
            self.name = 'Set'+str(int(len(SetOfTwoPt.instances))+1)
            SetOfTwoPt.instances.append(self.name)
        else:
            self.name = name
            # SetOfTwoPt.instances.append(self.name)

        self.SetC2 = OrderedDict()
        self.AutoList = []
        for ic,idict in enumerate(self.InfoDict):
            thiscol = colourset8[ic%len(colourset8)]
            thissym = markerset[ic%len(markerset)]
            thisshift = shiftset[ic%len(shiftset)]
            if 'autocorr' in list(idict.keys()):
                self.AutoList.append(idict['autocorr'])
            else:
                self.AutoList.append(False)
            thisC2 = NNQCorr(cfglist=self.cfglist, Info=idict,thissym=thissym,thiscol=thiscol,thisshift=thisshift)

            if str(thisC2) in list(self.SetC2.keys()):
                C2next = 2
                while str(thisC2)+'_'+str(C2next) in list(self.SetC2.keys()):
                    C2next += 1
                self.SetC2[str(thisC2)+'_'+str(C2next)] = thisC2
            else:
                self.SetC2[str(thisC2)] = thisC2
            # print 'Setting up classes for ' , thisC2 , ' complete'
            # print 'has color',thiscol , 'symbol ' , thissym , 'thisshift',thisshift

        self.SetC2len = len(self.SetC2)
        if not isinstance(self.cfglist,pa.DataFrame):
            if corr_cfglists:
                self.GetCfgList()
            else:
                self.GetUnCorrCfgList()
        else:
            self.ImportCfgList(self.cfglist)
            ##Implement Reading cfglist if len(cfglist) == 0

        self.VarKeys = []
        self.current = 0

    def GetUnCorrCfgList(self):
        for iC2 in self.SetC2.values():
            iC2.GetCfgList()

    def RemoveFuns(self):
        for iC2 in self.SetC2.values():
            iC2.RemoveFuns()
    def GetFuns(self):
        for iC2 in self.SetC2.values():
            iC2.GetFuns()


    def GetCfgList(self):
        self.GetUnCorrCfgList()
        cfglist = CombineListCfgs(np.array(list(self.SetC2.values())).flatten(),'NNQ_cfgs')
        self.ImportCfgList(cfglist)


    def ImportCfgList(self,cfglist):
        self.cfglist = cfglist
        for iC2 in self.SetC2.values():
            iC2.ImportCfgList(self.cfglist)



    # def Read(self):
    #     for setdata in self.SetC2.itervalues():
    #         print 'Importing 2pt corr: ' , str(setdata)
    #         setdata.()


    def Write(self):
        for iC2 in self.SetC2.values():
            print('Writing 2pt corr: ' , str(iC2))
            iC2.Write()

    ## I like LoadPickle the best
    def LoadPickle(self,DefWipe=False,CheckCfgs=True,CheckMom=True,NoFit=False):
        for ic,setdata in enumerate(self.SetC2.values()):
            print('LoadPickle NNQ corr: ' , str(setdata))
            setdata.LoadPickle(DefWipe=DefWipe,CheckCfgs=CheckCfgs,CheckMom=CheckMom,NoFit=NoFit)
            setdata.ImportPlotParams(colourset8[ic%len(colourset8)],markerset[ic%len(markerset)]
                                     ,shiftset[ic%len(shiftset)])
    # def ReadAndWrite(self,DefWipe=False):
    #     for setdata in self.SetC2.itervalues():
    #         print 'Importing and Writing 2pt corr: ' , str(setdata)
    #         setdata.ReadAndWrite(DefWipe=DefWipe)

    ## returns a TwoPtCorr class object (see TwoPtCorrelators.py)
    def GetSet(self,thisset='First'):
        if thisset == 'First':
            return self.SetC2[list(self.SetC2.keys())[0]]
        else:
            return self.SetC2[thisset]


    def AppendNoiseSub(self,indicies=[0,1],filename='',LegLab=''):
        if len(indicies) != 2:
            raise EnvironmentError('please pass in 2 indicies for noise subtraction')
        for ikey,iset in self.SetC2.items():
            if 'boot' not in iset.NNQ_Stats:
                print('WARNING, ', iset , ' was not read and boostrapped, attempting to initialize now')
                iset.LoadPickle(CheckCfgs=True)
        lset,rset = np.array(list(self.SetC2.values()))[indicies]
        if 'Alphaboot' not in lset.NNQ_Stats:
            lset.AlphaRatio()
        if 'Alphaboot' not in rset.NNQ_Stats:
            rset.AlphaRatio()
        if lset.ProjStr == 'g5':
            thisC2 = lset/2.-rset
        else:
            thisC2 = lset-rset
        thisC2.ProjStr = lset.ProjStr + '-'+rset.ProjStr
        thisC2.SetCustomName(string=filename,stringLL=LegLab)
        thisC2.IsSubComb = thisC2.IsComb = True
        thisC2.Stats()
        thisC2.Write()
        thiscol = colourset8[self.SetC2len%len(colourset8)]
        thissym = markerset[self.SetC2len%len(markerset)]
        thisshift = shiftset[self.SetC2len%len(shiftset)]
        thisC2.ImportPlotParams(thiscol,thissym,thisshift)

        if str(thisC2) in list(self.SetC2.keys()):
            C2next = 2
            while str(thisC2)+'_'+str(C2next) in list(self.SetC2.keys()):
                C2next += 1
            self.SetC2[str(thisC2)+'_'+str(C2next)] = thisC2
        else:
            self.SetC2[str(thisC2)] = thisC2
        self.VarKeys.append(list(self.SetC2.keys())[-1])
        self.SetC2len = len(self.SetC2)
        self.AutoList.append(self.AutoList[-1])

    ## thisfun(corr1,corr2,corr3,...)
    ## Strongly encoraged to manually set the filename somewhere
    ## e.g. sink smear combine we should get the first name and replace the sink smear with the sink smear label
    def CombineCorrs(self,thisfun,jsmcoeffs=None,filename='',LegLab='',jsmlab = '',indicies='All'):
        if indicies == 'All':
            indicies = list(range(len(list(self.SetC2.keys()))))


        for ikey,iset in self.SetC2.items():
            if 'boot' not in iset.NNQ_Stats:
                print('WARNING, ', iset , ' was not read and boostrapped, attempting to initialize now')
                iset.LoadPickle(CheckCfgs=True)

        if Debug_Mode:
            print('debugsetsoftwopt')
            print()
            print(thisfun.__name__)
            print()
            for ival in np.array(list(self.SetC2.values()))[indicies]:
                print(ival)

        output2pt = thisfun(*np.array(list(self.SetC2.values()))[indicies])

        # output2pt = thisfun(*self.SetC2.values())
        output2pt.SetCustomName(string=filename,stringLL=LegLab)
        output2pt.IsComb = True ## this needs to be done for alpha (maybe TODO move it into __add__ etc..?)
        output2pt.jsm = jsmlab
        if jsmcoeffs is not None: output2pt.Info['jsmCoeffs'] = jsmcoeffs
        # print 'DEBUG'
        # for ijsm,ival in zip(jsmcoeffs,np.array(self.SetC2.values())[indicies]):
        #     print ival
        #     print ijsm , '*' , ival.NNQ_Stats['boot']['p000','t13','t_f6.01'].Avg,ival.NNQ_Stats['boot']['p000','t13','t_f6.01'].Std
        # print output2pt
        # output2pt.NNQ_Stats['boot']['p000','t13','t_f6.01'].Stats()
        # print 'output',output2pt.NNQ_Stats['boot']['p000','t13','t_f6.01'].Avg,output2pt.NNQ_Stats['boot']['p000','t13','t_f6.01'].Std
        output2pt.Stats()
        output2pt.Write()
        return output2pt

    def AppendCCorrs(self,thisfun,jsmcoeffs=None,filename='',LegLab='',jsm = '',indicies='All'):
        thiscol = colourset8[self.SetC2len%len(colourset8)]
        thissym = markerset[self.SetC2len%len(markerset)]
        thisshift = shiftset[self.SetC2len%len(shiftset)]
        thisC2 = self.CombineCorrs(thisfun,filename=filename,jsmcoeffs=jsmcoeffs,LegLab=LegLab,jsmlab = jsm,indicies=indicies)
        thisC2.ImportPlotParams(thiscol,thissym,thisshift)

        if str(thisC2) in list(self.SetC2.keys()):
            C2next = 2
            while str(thisC2)+'_'+str(C2next) in list(self.SetC2.keys()):
                C2next += 1
            self.SetC2[str(thisC2)+'_'+str(C2next)] = thisC2
        else:
            self.SetC2[str(thisC2)] = thisC2
        self.VarKeys.append(list(self.SetC2.keys())[-1])
        self.SetC2len = len(self.SetC2)
        self.AutoList.append(self.AutoList[-1])


    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def AlphaPlotTauIntWopt(self,plot_info,momlist=[],tflowlist='PreDef',OnlyVar=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$\tau_{int}(W_{opt})$ for $\alpha$ over $t$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$t$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\tau_{int}(W_{opt})$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Alpha_t_tauintWopt.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        varcount = 0
        for icount,(ikey,iC2) in enumerate(self.SetC2.items()):
            if ((not OnlyVar) or ikey in  self.VarKeys) and self.AutoList[icount]:
                thiscol = colourset8[varcount%len(colourset8)]
                thissym = markerset[varcount%len(markerset)]
                thisshift = shiftset[varcount%len(shiftset)]
                varcount += 1
                print('Alpha tau int over time ', iC2)
                if tflowlist == 'PreDef':
                    data_plot = iC2.AlphaPlotTauIntWopt(data_plot,momlist=momlist[icount],
                                                        thiscol=thiscol,thissym=thissym,thisshift=thisshift)
                else:
                    itflow = tflowlist[icount]
                    # thiscol,thisshift,thissym = self.GetColShiftSym(iC2,mult=ictflow)
                    data_plot = iC2.AlphaPlotTauIntWopt(data_plot,thistflowlist = [itflow],
                                                        momlist=momlist[icount],
                                                        thiscol=thiscol,thissym=thissym,
                                                        thisshift=thisshift)
        data_plot.PlotAll()
        return data_plot

    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def AlphaPlot(  self,plot_info,fitDict=defFitDict,momlist=[],tflowlist='PreDef',
                    NNQfit=False,OnlyVar=False,AComp='PreDef'):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Nucleon Mixing Angle $\alpha$ over $t$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$t$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\alpha$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Alpha_t.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        if AComp != 'PreDef': self.AComp = AComp
        fit_is_dict = isinstance(fitDict,(OrderedDict,dict))
        varcount = 0
        for icount,(ikey,iC2) in enumerate(self.SetC2.items()):
            if (not OnlyVar) or ikey in  self.VarKeys:
                thiscol = colourset8[varcount%len(colourset8)]
                thissym = markerset[varcount%len(markerset)]
                thisshift = shiftset[varcount%len(shiftset)]
                varcount += 1
                print('Alpha plotting over time ', iC2)
                if tflowlist == 'PreDef':
                    if fit_is_dict:
                        thisxlim = [0,int(fitDict[ikey].split('-')[-1])]
                        data_plot = iC2.AlphaPlot(  data_plot,xlims=thisxlim,momlist=momlist[icount],
                                                thiscol=thiscol,thissym=thissym,thisshift=thisshift,
                                                Auto=self.AutoList[icount])
                        if ikey in list(fitDict.keys()):
                            if NNQfit:
                                data_plot = iC2.AlphaNNQFitPlot(data_plot,fitDict[ikey],
                                                momlist=momlist[icount],thiscol=thiscol,thisshift=thisshift,
                                                Auto=self.AutoList[icount])
                            else:
                                data_plot = iC2.AlphaFitPlot(data_plot,fitDict[ikey],
                                                momlist=momlist[icount],thiscol=thiscol,thisshift=thisshift,
                                                Auto=self.AutoList[icount])
                    else:
                        thisxlim = [0,int(fitDict[icount].split('-')[-1])]
                        data_plot = iC2.AlphaPlot(data_plot,xlims=thisxlim,momlist=momlist[icount],
                                                thiscol=thiscol,thissym=thissym,thisshift=thisshift,
                                                Auto=self.AutoList[icount])
                        if NNQfit:
                            data_plot = iC2.AlphaNNQFitPlot(data_plot,fitDict[icount],
                                                momlist=momlist[icount],Auto=self.AutoList[icount])
                        else:
                            data_plot = iC2.AlphaFitPlot(data_plot,fitDict[icount],
                                                momlist=momlist[icount],Auto=self.AutoList[icount])
                else:
                    itflow = tflowlist[icount]
                    # thiscol,thisshift,thissym = self.GetColShiftSym(iC2,mult=ictflow)
                    if fit_is_dict:
                        thisxlim = [0,int(fitDict[ikey].split('-')[-1])]
                        data_plot = iC2.AlphaPlot(data_plot,xlims=thisxlim,thistflowlist=[itflow],
                                                momlist=momlist[icount],thiscol=thiscol,thissym=thissym,
                                                thisshift=thisshift,Auto=self.AutoList[icount])
                        if ikey in list(fitDict.keys()):
                            if NNQfit:
                                data_plot = iC2.AlphaNNQFitPlot(data_plot,fitDict[ikey],thistflowlist=[itflow],
                                                momlist=momlist[icount],thiscol=thiscol,thisshift=thisshift,
                                                Auto=self.AutoList[icount])
                            else:
                                data_plot = iC2.AlphaFitPlot(data_plot,fitDict[ikey],thistflowlist=[itflow],
                                                momlist=momlist[icount],thiscol=thiscol,thisshift=thisshift,
                                                Auto=self.AutoList[icount])
                    else:
                        thisxlim = [0,int(fitDict[icount].split('-')[-1])]
                        data_plot = iC2.AlphaPlot(data_plot,xlims=thisxlim,thistflowlist = [itflow],
                                                momlist=momlist[icount],thiscol=thiscol,thissym=thissym,
                                                thisshift=thisshift,Auto=self.AutoList[icount])
                        if NNQfit:
                            data_plot = iC2.AlphaNNQFitPlot(data_plot,fitDict[icount],thistflowlist=[itflow],
                                                momlist=momlist[icount],thiscol=thiscol,thisshift=thisshift,
                                                Auto=self.AutoList[icount])
                        else:
                            data_plot = iC2.AlphaFitPlot(data_plot,fitDict[icount],thistflowlist=[itflow],
                                                momlist=momlist[icount],thiscol=thiscol,thisshift=thisshift,
                                                Auto=self.AutoList[icount])

        # sinklab,jsmcoeffs,sm2ptList,InfoList = GetInfoList(self.InfoDict[0])
        # InfoList = self.InfoDict
        # thisname,thisleglab,jsmlab = CreateNameAndLL(self.SetC2,sinklab)
        # thisfun = GenSinkFun(jsmcoeffs)
        # # print thisname,thisleglab
        # # print jsmcoeffs
        # # print sinklab
        # # print thisfun.__name__
        # # print self.FOSet.SetC2.keys()
        # # print self.FOSet.SetC2[self.FOSet.SetC2.keys()[0]].NNQboot
        # combcorr = self.CombineCorrs(thisfun,thisname,thisleglab,jsmlab)

        # combcorr.AlphaPlot(momlist=momlist[icount])

        ## compares it to what ever the information dictornarys first element parameters are
        if self.AComp:
            raise EnvironmentError('Ahmed comparison not implemented for updated plotting version')
            # thiscol = colourset8[varcount+1%len(colourset8)]
            # thissym = markerset[varcount+1%len(markerset)]
            # thisshift = shiftset[varcount+1%len(shiftset)]
            # Adata = ara.AResAlpha(Info=self.InfoDict[0],thissym=thissym,thiscol=thiscol,thisshift=thisshift,Auto=self.AutoList[0])
            # Adata.ReadAlpha()
            # Adata.PlotAlpha(xlims=thisxlim)
        data_plot.PlotAll()
        return data_plot




    ## this is used for plotting multiple versions of the same results (over tflow or tsink if you like)
    def GetColShiftSym(self,iC2,mult=0):
        thiscol = colourset8[(colourset8.index(iC2.thiscol)+self.SetC2len*mult)%len(colourset8)]
        thissym = markerset[(markerset.index(iC2.thissym)+self.SetC2len*mult)%len(markerset)]
        thisshift = shiftset[(shiftset.index(iC2.thisshift)+self.SetC2len*mult)%len(shiftset)]
        return thiscol,thisshift,thissym

    def AlphaPlotTflowTauIntWopt(self,plot_info,tflowlist=defxlimAlphaTflow,tlist='PreDef',momlist=[],OnlyVar=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$\tau_{int}(W_{opt})$ for $\alpha$ over $t_{f}$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$\sqrt{8t_{f}}[fm]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\tau_{int}(W_{opt})$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Alpha_tflow_tauintWopt.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        tflow_is_list  = isinstance(tflowlist[0],(list,tuple,np.ndarray))
        varcount = 0
        for icount,(ikey,iC2) in enumerate(self.SetC2.items()):
            if ((not OnlyVar) or ikey in  self.VarKeys) and self.AutoList[icount]:
                thiscol = colourset8[varcount%len(colourset8)]
                thissym = markerset[varcount%len(markerset)]
                thisshift = shiftset[varcount%len(shiftset)]
                varcount += 1
                print('Alpha tau int over flow time ', iC2)
                if tlist == 'PreDef':
                    if tflow_is_list:
                        data_plot = iC2.AlphaPlotTflowTauIntWopt(data_plot,tflowplot=tflowlist[icount],
                                            momlist=momlist[icount],thiscol=thiscol,thissym=thissym,
                                            thisshift=thisshift)
                    else:
                        data_plot = iC2.AlphaPlotTflowTauIntWopt(data_plot,tflowplot=tflowlist,
                                            momlist=momlist[icount],thiscol=thiscol,thissym=thissym,
                                            thisshift=thisshift)
                else:
                    thist = tlist[icount]
                    if tflow_is_list:
                        data_plot = iC2.AlphaPlotTflowTauIntWopt(data_plot,tflowplot=tflowlist[icount],
                                            thistsinklist = [thist],momlist=momlist[icount],
                                            thiscol=thiscol,thissym=thissym,thisshift=thisshift)
                    else:
                        data_plot = iC2.AlphaPlotTflowTauIntWopt(data_plot,tflowplot=tflowlist,
                                            thistsinklist = [thist],momlist=momlist[icount],
                                            thiscol=thiscol,thissym=thissym,thisshift=thisshift)
                # if ikey in fitDict.keys():
                #     iC2.AlphaFitPlot(fitDict[ikey])
        data_plot.PlotAll()
        return data_plot


    def AlphaPlotTflow(self,plot_info,tflowlist=defxlimAlphaTflow,tlist='PreDef',momlist=[],OnlyVar=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Nucleon Mixing Angle $\alpha$ over $t_{f}$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$\sqrt{8t_{f}}[fm]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\alpha$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Alpha_tflow.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        tflow_is_list  = isinstance(tflowlist[0],(list,tuple,np.ndarray))
        varcount = 0
        for icount,(ikey,iC2) in enumerate(self.SetC2.items()):
            if (not OnlyVar) or ikey in  self.VarKeys:
                thiscol = colourset8[varcount%len(colourset8)]
                thissym = markerset[varcount%len(markerset)]
                thisshift = shiftset[varcount%len(shiftset)]
                varcount += 1
                print('Alpha plotting over flow time ', iC2)
                if tlist == 'PreDef':
                    if tflow_is_list:
                        data_plot = iC2.AlphaPlotTflow(data_plot,tflowplot=tflowlist[icount],
                                            momlist=momlist[icount],thiscol=thiscol,thissym=thissym,
                                            thisshift=thisshift,Auto=self.AutoList[icount])
                    else:
                        data_plot = iC2.AlphaPlotTflow(data_plot,tflowplot=tflowlist,momlist=momlist[icount],
                                            thiscol=thiscol,thissym=thissym,thisshift=thisshift,
                                            Auto=self.AutoList[icount])
                else:
                    thist = tlist[icount]
                    if tflow_is_list:
                        data_plot = iC2.AlphaPlotTflow(data_plot,tflowplot=tflowlist[icount],
                                            thistsinklist = [thist],momlist=momlist[icount],
                                            thiscol=thiscol,thissym=thissym,thisshift=thisshift,
                                            Auto=self.AutoList[icount])
                    else:
                        data_plot = iC2.AlphaPlotTflow(data_plot,tflowplot=tflowlist,thistsinklist=[thist],
                                            momlist=momlist[icount],thiscol=thiscol,thissym=thissym,
                                            thisshift=thisshift,Auto=self.AutoList[icount])
                # if ikey in fitDict.keys():
                #     iC2.AlphaFitPlot(fitDict[ikey])
        data_plot.PlotAll()
        return data_plot


    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    ## momlist is now all of them. THIS IS DIFFERENT TO ABOVE
    def AlphaPlotVsEnergy(self,plot_info,fitDict=defFitDict,momlist=[],tflowlist='PreDef',OnlyVar=False,WipeFit=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Nucleon Mixing Angle $\alpha$ vs energy extraction from'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$E_{p}$ [GeV]'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\alpha$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Alpha_VsEnergy.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        fit_is_dict = isinstance(fitDict,(OrderedDict,dict))
        varcount = 0
        for icount,(ikey,iC2) in enumerate(self.SetC2.items()):
            if (not OnlyVar) or ikey in  self.VarKeys:
                thiscol = colourset8[varcount%len(colourset8)]
                thissym = markerset[varcount%len(markerset)]
                thisshift = shiftset[varcount%len(shiftset)]
                varcount += 1
                print('Alpha plotting over time ', iC2)
                if tflowlist == 'PreDef':
                    if fit_is_dict and ikey in list(fitDict.keys()):
                        data_plot = iC2.AlphaPlotVsEnergy(data_plot,fitDict[ikey],momlist=momlist[icount],
                                                thiscol=thiscol,thissym=thissym,thisshift=thisshift,
                                                Auto=self.AutoList[icount],WipeFit=WipeFit)
                    else:
                        data_plot = iC2.AlphaPlotVsEnergy(data_plot,fitDict[icount],momlist=momlist[icount],
                                                thiscol=thiscol,thissym=thissym,thisshift=thisshift,
                                                Auto=self.AutoList[icount],WipeFit=WipeFit)
                else:
                    itflow = tflowlist[icount]
                    # thiscol,thisshift,thissym = self.GetColShiftSym(iC2,mult=ictflow)
                    if fit_is_dict:
                        data_plot = iC2.AlphaPlotVsEnergy(data_plot,fitDict[ikey],thistflowlist=[itflow],
                                                momlist=momlist[icount],thiscol=thiscol,thissym=thissym,
                                                thisshift=thisshift,Auto=self.AutoList[icount],WipeFit=WipeFit)
                    else:
                        data_plot = iC2.AlphaPlotVsEnergy(data_plot,fitDict[icount],thistflowlist=[itflow],
                                                momlist=momlist[icount],thiscol=thiscol,thissym=thissym,
                                                thisshift=thisshift,Auto=self.AutoList[icount],WipeFit=WipeFit)
        data_plot.PlotAll()
        return data_plot


    def AlphaPlotVsLatSpace(self,plot_info,this_fit,plottflowlist='PreDefine',momlist='p000',OnlyVar=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$\alpha_{N}$ Lattice spacing $a$ dependance'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$ a[fm]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\alpha_{N}$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_VsLatSpace.pdf'
            if 'xlims' in this_pi:
                this_pi['xlims'] = [0,this_pi['xlims'][1]]
            else:
                this_pi['xlims'] = [0,None]
            return this_pi

        if isinstance(this_fit,str):
            this_fit = [this_fit]*len(list(self.SetC2.keys()))
        if isinstance(momlist,str):
            momlist = [momlist]*len(list(self.SetC2.keys()))
        if isinstance(plottflowlist,str):
            plottflowlist = [plottflowlist]*len(list(self.SetC2.keys()))
        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        for ikey,iset in self.SetC2.items():
            if (not OnlyVar) or ikey in  self.VarKeys:
                firstset = iset
                break
        the_mpi = ''
        this_bool = False
        for ikey,iset in self.SetC2.items():
            if (not OnlyVar) or ikey in  self.VarKeys:
                if firstset.NN.latparams.GetPionMass() != iset.NN.latparams.GetPionMass():
                    this_bool = True
                    print('warning: not all mpi values are the same for latspace plot')
                    break
        if this_bool:
            for ikey,iset in self.SetC2.items():
                if (not OnlyVar) or ikey in  self.VarKeys:
                    print(iset.name)
                    print(iset.NN.latparams.GetPionMass())
                    print()
        else:
            the_mpi = firstset.NN.latparams.GetPionMassLab(MeV=True)
        print('Plotting fit results over lattice spacing')
        ploty,plotx = [],[]
        for ic,(ikey,iset) in enumerate(self.SetC2.items()):
            if (not OnlyVar) or ikey in  self.VarKeys:
                fit_data = iset.Get_Extrapolation()
                if len(fit_data.index) == 0: continue
                this_fun = fit_data.index[0][2]
                this_par = fit_data.index[0][4]
                if not isinstance(momlist[ic],str):
                    this_key = momlist[ic][0],plottflowlist[ic],this_fun,this_fit[ic],this_par
                else:
                    this_key = momlist[ic],plottflowlist[ic],this_fun,this_fit[ic],this_par
                if 'all' in this_key[0].lower():
                    this_key = list(this_key)
                    this_key[0] = fit_data.index[0][0]
                    this_key = tuple(this_key)
                if this_key not in fit_data:
                    print(this_key, ' fit has not been computed for alpha, attempting now')
                    iset.FitAlpha(fit_range=this_fit[ic])
                    fit_data = iset.Get_Extrapolation()
                    if len(fit_data.index) == 0: continue
                    this_fun = fit_data.index[0][2]
                    this_par = fit_data.index[0][4]
                    if not isinstance(momlist[ic],str):
                        this_key = momlist[ic][0],plottflowlist[ic],this_fun,this_fit[ic],this_par
                    else:
                        this_key = momlist[ic],plottflowlist[ic],this_fun,this_fit[ic],this_par
                    if 'all' in this_key[0].lower():
                        this_key = list(this_key)
                        this_key[0] = fit_data.index[0][0]
                        this_key = tuple(this_key)
                    if this_key not in fit_data:
                        print(' ERROR FITTING:' ,  this_key)
                        print(fit_data)
                        continue
                this_names = ['$a$']+list(fit_data.index.names)
                for ireskey,ires in fit_data.items():
                    ploty.append(ires)
                    plotx.append((iset.NN.latparams.GetLatSpaceLab(),)+ireskey) ## in MeV
        if len(plotx) > 0:
            plot_key = [slice(None)] + list(plotx[0][1:])
            indicies = pa.MultiIndex.from_tuples(plotx,names=this_names)
            ploty = pa.Series(ploty,index=indicies)
        else:
            print('No Sets Found for plotting vs lat spacing')
            return data_plot

        # thisFF = thisFF.replace('tilde',r'\tilde')
        # thisFF = thisFF.replace('frac',r'\frac')
        thislab = r'$'+r'\ '.join([the_mpi])+r'$ '
        hold_series = copy(jpl.null_series)
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = ploty.apply(lambda x : x.Avg)
        hold_series['yerr_data'] = ploty.apply(lambda x : x.Std)
        hold_series['key_select'] = plot_key
        hold_series['type'] = 'error_bar_vary'
        hold_series['label'] = thislab
        hold_series['symbol'] = 'o'
        data_plot.AppendData(hold_series)
        data_plot = self.Plot_Extrap(data_plot,'latspace')
        data_plot.PlotAll()
        return data_plot


    def AlphaPlotVsMpi(self,plot_info,this_fit,plottflowlist='PreDefine',momlist='p000',OnlyVar=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$\alpha_{N}$ $m_{\pi}$ dependance'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$ m_{\pi}[MeV]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\alpha_{N}$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_VsMpi.pdf'
            if 'xlims' in this_pi:
                this_pi['xlims'] = [0,this_pi['xlims'][1]]
            else:
                this_pi['xlims'] = [0,None]
            return this_pi

        if isinstance(this_fit,str):
            this_fit = [this_fit]*len(list(self.SetC2.keys()))
        if isinstance(momlist,str):
            momlist = [momlist]*len(list(self.SetC2.keys()))
        if isinstance(plottflowlist,str):
            plottflowlist = [plottflowlist]*len(list(self.SetC2.keys()))
        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        for ikey,iset in self.SetC2.items():
            if (not OnlyVar) or ikey in  self.VarKeys:
                firstset = iset
                break
        the_latspace = ''
        this_bool = False
        for ikey,iset in self.SetC2.items():
            if (not OnlyVar) or ikey in  self.VarKeys:
                if firstset.NN.latparams.latspace != iset.NN.latparams.latspace:
                    this_bool = True
                    print('warning: not all lattice spacing values are the same for mpi plot')
                    break
        if this_bool:
            for ikey,iset in self.SetC2.items():
                if (not OnlyVar) or ikey in  self.VarKeys:
                    print(iset.name)
                    print(iset.NN.latparams.latspace)
                    print()
        else:
            the_latspace = firstset.NN.latparams.GetLatSpaceLab()
        print('Plotting fit results over pion mass')

        ploty,plotx = [],[]
        for ic,(ikey,iset) in enumerate(self.SetC2.items()):
            if (not OnlyVar) or ikey in  self.VarKeys:
                fit_data = iset.Get_Extrapolation()
                if len(fit_data.index) == 0: continue
                this_fun = fit_data.index[0][2]
                this_par = fit_data.index[0][4]
                if not isinstance(momlist[ic],str):
                    this_key = momlist[ic][0],plottflowlist[ic],this_fun,this_fit[ic],this_par
                else:
                    this_key = momlist[ic],plottflowlist[ic],this_fun,this_fit[ic],this_par
                if 'all' in this_key[0].lower():
                    this_key = list(this_key)
                    this_key[0] = fit_data.index[0][0]
                    this_key = tuple(this_key)
                if this_key not in fit_data:
                    print(this_key, ' fit has not been computed for alpha, attempting now')
                    iset.FitAlpha(fit_range=this_fit[ic])
                    fit_data = iset.Get_Extrapolation()
                    this_fun = fit_data.index[0][2]
                    this_par = fit_data.index[0][4]
                    if not isinstance(momlist[ic],str):
                        this_key = momlist[ic][0],plottflowlist[ic],this_fun,this_fit[ic],this_par
                    else:
                        this_key = momlist[ic],plottflowlist[ic],this_fun,this_fit[ic],this_par
                    if 'all' in this_key[0].lower():
                        this_key = list(this_key)
                        this_key[0] = fit_data.index[0][0]
                        this_key = tuple(this_key)
                    if this_key not in fit_data:
                        print(' ERROR FITTING:' ,  this_key)
                        print(fit_data)
                        continue
                this_names = ['$m_{\pi}$']+list(fit_data.index.names)
                for ireskey,ires in iset.Get_Extrapolation().items():
                    ploty.append(ires)
                    plotx.append((iset.NN.latparams.GetPionMassLab(MeV=True),)+ireskey) ## in MeV
        if len(plotx) > 0:
            plot_key = [slice(None)] + list(plotx[0][1:])
            indicies = pa.MultiIndex.from_tuples(plotx,names=this_names)
            ploty = pa.Series(ploty,index=indicies)
        else:
            print('No Sets Found for plotting vs mpi')
            return data_plot

        # thisFF = thisFF.replace('tilde',r'\tilde')
        # thisFF = thisFF.replace('frac',r'\frac')
        thislab = r'$'+r'\ '.join([the_latspace])+r'$ '
        hold_series = copy(jpl.null_series)
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = ploty.apply(lambda x : x.Avg)
        hold_series['yerr_data'] = ploty.apply(lambda x : x.Std)
        hold_series['key_select'] = plot_key
        hold_series['type'] = 'error_bar_vary'
        hold_series['label'] = thislab
        hold_series['symbol'] = 'o'
        data_plot.AppendData(hold_series)
        data_plot = self.Plot_Extrap(data_plot,'chiral')
        data_plot.PlotAll()
        return data_plot


    def Plot_Extrap(self,plot_data,extrap_type,this_fun='Default'):
        extrap_type = extrap_type.lower()
        this_eval = 0
        if 'latspace' in extrap_type:
            if isinstance(this_fun,str):
                this_fun = (SquaredFitFun,2)
            fit_data = self.LatSpace_Extrapolation(this_fun=this_fun)
        # elif 'boxsize' in extrap_type:
        #     if isinstance(this_fun,str):
        #         this_fun = (SquaredFitFun,2)
        #     fit_data = self.BoxSize_Extrapolation(this_fun=this_fun)
        elif 'mpi' in extrap_type or 'chiral' in extrap_type:
            if isinstance(this_fun,str):
                this_fun = (LinearFitFun,2)
            fit_data = self.Chiral_Extrapolation(this_fun=this_fun)
            this_eval = physcial_pion_mass
        else:
            raise IOError('extrapolation type not recognised '+str(extrap_type))
        if len(fit_data) == 0:
            return plot_data
        # if isinstance(fitr_params,pa.Series): fitr_params = fitr_params.iloc[-1]
        ## TODO, plot from zero to data max
        hold_series = copy(jpl.null_series)
        hold_series['key_select'] = fit_data.index[0]
        hold_series['type'] = 'fit_vary'
        hold_series['fit_class'] = fit_data
        hold_series['label'] = 'fit '
        hold_series['xdatarange'] = [0,'ToDataMax']
        hold_series['ShowEval'] = this_eval
        plot_data.AppendData(hold_series)
        return plot_data


    def LatSpace_Extrapolation(self,this_fun=(SquaredFitFun,2)):
        # print 'DEBUGGING LatSpace_Extrapolation'
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
                    # ydata.append(coldata+200)
                    # ydata[-1].Stats()
            thisff = ff.Fitting(Funs=this_fun,data=[np.array(xdata),np.array(ydata)],name='LatSpace_Ext')
            fit_result = thisff.FitBoots()
            if np.isnan(fit_result):
                return fit_result
            return thisff

        this_df = pa.DataFrame()
        for ikey,iset in self.SetC2.items():
            this_data = iset.Get_Extrapolation()
            if len(this_data) > 0:
                this_df[ikey] = this_data
                this_df.loc[:,ikey+'_latspace'] = float(iset.NN.latparams.latspace)
        this_df.loc[:,'Fit'] = this_df.apply(this_fit,axis=1)
        return this_df['Fit'].dropna()

    def Chiral_Extrapolation(self,this_fun=(LinearFitFun,2)):
        # print 'DEBUGGING Chiral_Extrapolation'
        def this_fit(this_series):
            xdata,ydata = [[]],[]
            for icol,coldata in this_series.items():
                if '_chiral' in icol:
                    xdata[0].append(coldata)
                    # ##### DEBUGGING ####
                    # xdata[0].append(coldata+100)
                else:
                    ydata.append(coldata)
                    # ##### DEBUGGING ####
                    # ydata.append(coldata+200)
                    # ydata[-1].Stats()
            thisff = ff.Fitting(Funs=this_fun,data=[np.array(xdata),np.array(ydata)],name='Chiral_Ext')
            fit_result = thisff.FitBoots()
            if np.isnan(fit_result):
                return fit_result
            return thisff

        this_df = pa.DataFrame()
        for ikey,iset in self.SetC2.items():
            this_data = iset.Get_Extrapolation()
            if len(this_data) > 0:
                this_df[ikey] = this_data
                this_df.loc[:,ikey+'_chiral'] = float(iset.NN.latparams.GetPionMass(MeV=True))
        this_df.loc[:,'Fit'] = this_df.apply(this_fit,axis=1)
        return this_df['Fit'].dropna()


    def __setitem__(self, key, value ):
        self.SetC2[key[0]][key[1:]] = value

    def __getitem__(self, key):
        return self.SetC2[key[0]][key[1:]]

    def Stats(self):
        for iset in self.SetC2.values():
            iset.Stats()

    def __str__(self):
        return self.name

    def __iter__(self):
        return self

    def items(self):
        outdata = []
        for iset,setdata in self.SetC2.items():
            for ikey,tdata in list(setdata.items()):
                outdata.append((iset,)+ikey+(tdata,))
        return outdata

    def itemsAvgStd(self):
        outdata = []
        for iset,setdata in self.SetC2.items():
            for ikey,tdata in list(setdata.items()):
                outdata.append((iset,)+ikey+(tdata.Avg,tdata.Std))
        return outdata

    def values(self):
        outdata = []
        for setdata in self.SetC2.values():
            for tdata in list(setdata.values()):
                outdata.append(tdata)
        return outdata



    def valuesAvgStd(self):
        outdata = []
        for setdata in self.SetC2.values():
            for tdata in list(setdata.values()):
                outdata.append((tdata.Avg,tdata.Std))
        return outdata

    def keys(self):
        outdata = []
        for iset,setdata in self.SetC2.items():
            for ikey in list(setdata.keys()):
                outdata.append((iset,)+ikey)
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
            np.complex(self[ikey])
            self[ikey].Stats()
        return 1+0j

    def __int__(self):
        for ikey in list(self.keys()):
            int(self[ikey])
            self[ikey].Stats()
        return 1

    def __float__(self):
        for ikey in list(self.keys()):
            np.float(self[ikey])
            self[ikey].Stats()
        return 1.0









##############################################################################################
##############################################################################################
##############################################################################################
##############################################################################################



class SetOfNNQFull(object):

    """

    NNQFull(set) Inherits from TwoPtCorrelators Class

    """
    instances = []

    parentlist = []


    ## InfoDict is a list of Info dictionaries to be bassed to TwoPtCorr:
    ## Setup to do fits taken from dictionaries within InfoDict (see TwoPtCorrelators.py)
    def __init__(self, cfglist={}, InfoDict=[],corr_cfglists=True,name = ''):
        from TwoPtCorrelators import NNQFullCorr

        self.cfglist = cfglist
        self.InfoDict = InfoDict
        self.Plotdir = graphdir+'/AlphaFull/'
        mkdir_p(self.Plotdir)


        if len(name) == 0:
            self.name = 'Set'+str(int(len(SetOfTwoPt.instances))+1)
            SetOfTwoPt.instances.append(self.name)
        else:
            self.name = name
            # SetOfTwoPt.instances.append(self.name)

        self.SetC2 = OrderedDict()
        self.AutoList = []
        for ic,idict in enumerate(self.InfoDict):
            if 'autocorr' in list(idict.keys()):
                self.AutoList.append(idict['autocorr'])
            else:
                self.AutoList.append(False)
            thiscol = colourset8[ic%len(colourset8)]
            thissym = markerset[ic%len(markerset)]
            thisshift = shiftset[ic%len(shiftset)]
            thisC2 = NNQFullCorr(cfglist=self.cfglist, Info=idict,thissym=thissym,thiscol=thiscol,thisshift=thisshift)

            if str(thisC2) in list(self.SetC2.keys()):
                C2next = 2
                while str(thisC2)+'_'+str(C2next) in list(self.SetC2.keys()):
                    C2next += 1
                self.SetC2[str(thisC2)+'_'+str(C2next)] = thisC2
            else:
                self.SetC2[str(thisC2)] = thisC2
            # print 'Setting up classes for ' , thisC2 , ' complete'
            # print 'has color',thiscol , 'symbol ' , thissym , 'thisshift',thisshift

        self.SetC2len = len(self.SetC2)
        if not isinstance(self.cfglist,pa.DataFrame):
            if corr_cfglists:
                self.GetCfgList()
            else:
                self.GetUnCorrCfgList()
        else:
            self.ImportCfgList(self.cfglist)
            ##Implement Reading cfglist if len(cfglist) == 0

        self.VarKeys = []
        self.current = 0

    def GetUnCorrCfgList(self):
        for iC2 in self.SetC2.values():
            iC2.GetCfgList()

    def RemoveFuns(self):
        for iC2 in self.SetC2.values():
            iC2.RemoveFuns()
    def GetFuns(self):
        for iC2 in self.SetC2.values():
            iC2.GetFuns()


    def GetCfgList(self):
        self.GetUnCorrCfgList()
        cfglist = CombineListCfgs(np.array(list(self.SetC2.values())).flatten(),'NNQFull_cfgs')
        self.ImportCfgList(cfglist)

    def ImportCfgList(self,cfglist):
        self.cfglist = cfglist
        for iC2 in self.SetC2.values():
            iC2.ImportCfgList(self.cfglist)



    # def Read(self):
    #     for setdata in self.SetC2.itervalues():
    #         print 'Importing 2pt corr: ' , str(setdata)
    #         setdata.()


    def Write(self):
        for iC2 in self.SetC2.values():
            print('Writing 2pt corr: ' , str(iC2))
            iC2.Write()

    ## I like LoadPickle the best
    def LoadPickle(self,DefWipe=False,CheckCfgs=False,CheckMom=True,show_timer=True,SaveMem=False,NoFit=False):
        if show_timer: thistimer = Timer(linklist=list(map(str,list(self.SetC2.values()))),name='LP NNQFull corr: ')
        for ic,setdata in enumerate(self.SetC2.values()):
            setdata.LoadPickle(DefWipe=DefWipe,CheckCfgs=CheckCfgs,CheckMom=CheckMom,SaveMem=SaveMem,NoFit=NoFit)
            setdata.ImportPlotParams(colourset8[ic%len(colourset8)],markerset[ic%len(markerset)]
                                     ,shiftset[ic%len(shiftset)])
            if show_timer: thistimer.Lap()
    # def ReadAndWrite(self,DefWipe=False):
    #     for setdata in self.SetC2.itervalues():
    #         print 'Importing and Writing 2pt corr: ' , str(setdata)
    #         setdata.ReadAndWrite(DefWipe=DefWipe)

    ## returns a TwoPtCorr class object (see TwoPtCorrelators.py)
    def GetSet(self,thisset='First'):
        if thisset == 'First':
            return self.SetC2[list(self.SetC2.keys())[0]]
        else:
            return self.SetC2[thisset]


    ## thisfun(corr1,corr2,corr3,...)
    ## Strongly encoraged to manually set the filename somewhere
    ## e.g. sink smear combine we should get the first name and replace the sink smear with the sink smear label
    def CombineCorrs(self,thisfun,jsmcoeffs=None,filename='',LegLab='',jsmlab = '',indicies='All'):
        if indicies == 'All':
            indicies = list(range(len(list(self.SetC2.keys()))))


        output2pt = thisfun(*np.array(list(self.SetC2.values()))[indicies])

        # output2pt = thisfun(*self.SetC2.values())
        output2pt.SetCustomName(string=filename,stringLL=LegLab)
        output2pt.IsComb = True ## this needs to be done for alpha (maybe TODO move it into __add__ etc..?)
        output2pt.jsm = jsmlab
        if jsmcoeffs is not None: output2pt.Info['jsmCoeffs'] = jsmcoeffs
        # print 'DEBUG'
        # for ijsm,ival in zip(jsmcoeffs,np.array(self.SetC2.values())[indicies]):
        #     print ival
        #     print ijsm , '*' , ival.NNQFull_Stats['boot']['p000','t13','t_f6.01'].Avg,ival.NNQFull_Stats['boot']['p000','t13','t_f6.01'].Std
        # print output2pt
        # output2pt.NNQFull_Stats['boot']['p000','t13','t_f6.01'].Stats()
        # print 'output',output2pt.NNQFull_Stats['boot']['p000','t13','t_f6.01'].Avg,output2pt.NNQFull_Stats['boot']['p000','t13','t_f6.01'].Std
        output2pt.Stats()
        output2pt.Write()
        return output2pt

    def AppendCCorrs(self,thisfun,jsmcoeffs=None,filename='',LegLab='',jsm = '',indicies='All'):
        thiscol = colourset8[self.SetC2len%len(colourset8)]
        thissym = markerset[self.SetC2len%len(markerset)]
        thisshift = shiftset[self.SetC2len%len(shiftset)]
        thisC2 = self.CombineCorrs(thisfun,filename=filename,jsmcoeffs=jsmcoeffs,LegLab=LegLab,jsmlab = jsm,indicies=indicies)
        thisC2.ImportPlotParams(thiscol,thissym,thisshift)

        if str(thisC2) in list(self.SetC2.keys()):
            C2next = 2
            while str(thisC2)+'_'+str(C2next) in list(self.SetC2.keys()):
                C2next += 1
            self.SetC2[str(thisC2)+'_'+str(C2next)] = thisC2
        else:
            self.SetC2[str(thisC2)] = thisC2
        self.VarKeys.append(list(self.SetC2.keys())[-1])
        self.SetC2len = len(self.SetC2)
        self.AutoList.append(False)

    def AlphaPlotVstfFitr(  self,tflowlist,tfitlist,fit_min_list,plot_info,
                            momlist='p000',fit_max_list='All',OnlyVar=False,WipeFit=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Nucleon Mixing Angle $\alpha$ over fit ranges'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'max of fit range $t_{e}$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\alpha$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Alpha_Fit_Max.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        varcount = 0
        if isinstance(momlist,str): momlist = [momlist]*len(list(self.SetC2.keys()))
        if isinstance(fit_max_list,str): fit_max_list = [fit_max_list]*len(list(self.SetC2.keys()))
        for icount,(ikey,iC2) in enumerate(self.SetC2.items()):
            if (not OnlyVar) or ikey in  self.VarKeys:
                thiscol = colourset8[varcount%len(colourset8)]
                thissym = markerset[varcount%len(markerset)]
                thisshift = shiftset[varcount%len(shiftset)]
                varcount += 1
                print('Alpha plotting over max fit ranges ', iC2)
                # thiscol,thisshift,thissym = self.GetColShiftSym(iC2,mult=ictflow)
                thismom = momlist[icount]
                if not isinstance(thismom,str):
                    thismom = thismom[0]
                data_plot = iC2.AlphaPlotVstfFitr(data_plot,tflowlist[icount],tfitlist[icount],
                            fit_min_list[icount],thismom=thismom,WipeFit=WipeFit,
                            thiscol=thiscol,thissym=thissym,thisshift=thisshift,
                            fit_max_list=fit_max_list[icount])
        data_plot.PlotAll()
        return data_plot


    def AlphaPlotVstiFitr(  self,tflowlist,tfitlist,fit_max_list,plot_info,
                            momlist='p000',fit_min_list='All',OnlyVar=False,WipeFit=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Nucleon Mixing Angle $\alpha$ over fit ranges'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'min of fit range $t_{i}$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\alpha$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Alpha_Fit_Min.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        varcount = 0
        if isinstance(momlist,str): momlist = [momlist]*len(list(self.SetC2.keys()))
        if isinstance(fit_min_list,str): fit_min_list = [fit_min_list]*len(list(self.SetC2.keys()))
        for icount,(ikey,iC2) in enumerate(self.SetC2.items()):
            if (not OnlyVar) or ikey in  self.VarKeys:
                thiscol = colourset8[varcount%len(colourset8)]
                thissym = markerset[varcount%len(markerset)]
                thisshift = shiftset[varcount%len(shiftset)]
                varcount += 1
                print('Alpha plotting over min fit ranges ', iC2)
                # thiscol,thisshift,thissym = self.GetColShiftSym(iC2,mult=ictflow)
                thismom = momlist[icount]
                if not isinstance(thismom,str):
                    thismom = thismom[0]
                data_plot = iC2.AlphaPlotVstiFitr(data_plot,tflowlist[icount],tfitlist[icount],
                            fit_max_list[icount],thismom=thismom,WipeFit=WipeFit,
                            thiscol=thiscol,thissym=thissym,thisshift=thisshift,
                            fit_min_list=fit_min_list[icount])
        data_plot.PlotAll()
        return data_plot

    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def AlphaPlot(  self,plot_info,fitrlist='None',momlist='p000',tflowlist='PreDef',
                    tsum_list='PreDef',spec_tsum_list='PreDef',
                    OnlyVar=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Nucleon Mixing Angle $\alpha$ over $t$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$t$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\alpha$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Alpha_t.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        varcount = 0
        if fitrlist == 'PreDef': fitrlist = [fitrlist]*len(list(self.SetC2.keys()))
        if tflowlist == 'PreDef': tflowlist = [tflowlist]*len(list(self.SetC2.keys()))
        if tsum_list == 'PreDef': tsum_list = [tsum_list]*len(list(self.SetC2.keys()))
        if spec_tsum_list == 'PreDef': spec_tsum_list = [spec_tsum_list]*len(list(self.SetC2.keys()))
        for icount,(ikey,iC2) in enumerate(self.SetC2.items()):
            if (not OnlyVar) or ikey in  self.VarKeys:
                fitr = fitrlist[icount]
                itflow = tflowlist[icount]
                itsum = tsum_list[icount]
                itsum_spec = spec_tsum_list[icount]
                thiscol = colourset8[varcount%len(colourset8)]
                thissym = markerset[varcount%len(markerset)]
                thisshift = shiftset[varcount%len(shiftset)]
                varcount += 1
                print('Alpha plotting over time ', iC2)
                # thiscol,thisshift,thissym = self.GetColShiftSym(iC2,mult=ictflow)
                data_plot = iC2.AlphaPlot(data_plot,
                                thistflowlist = itflow,momlist=momlist[icount],
                                tsum_list=itsum,spec_tsum_list=itsum_spec,
                                thiscol=thiscol,thissym=thissym,thisshift=thisshift,Auto=self.AutoList[icount])
                if 'none' not in fitr.lower() and not self.AutoList[icount]:
                    data_plot = iC2.Plot_Tsink_Fit(data_plot,fitr,momlist[icount],itsum,itflow,thiscol,thisshift)
                thiscol = colourset8[varcount%len(colourset8)]
                thissym = markerset[varcount%len(markerset)]
                thisshift = shiftset[varcount%len(shiftset)]
                varcount += 1
                data_plot = iC2.PlotIntTsink(data_plot,momlist[icount],itflow,Auto=self.AutoList[icount],
                                thiscol=thiscol,thissym=thissym,thisshift=thisshift)
        data_plot.PlotAll()
        return data_plot



    # def AlphaPlot(self,plot_info,fitDict=defFitDict,momlist=[],tflowlist='PreDef',
    #                 parlist = 'First',OnlyVar=False,WipeFit=False,No_Fit=False):
    #     def DefPlotInfo(pi):
    #         this_pi = copy(pi)
    #         if 'title' not in this_pi:
    #             this_pi['title'] = r'Nucleon Mixing Angle $\alpha$ over $t$'
    #         if 'xlabel' not in this_pi:
    #             this_pi['xlabel'] = r'$t$'
    #         if 'ylabel' not in this_pi:
    #             this_pi['ylabel'] = r'$\alpha$'
    #         if 'save_file' not in this_pi:
    #             this_pi['save_file'] = self.Plotdir+self.name+'_Alpha_t.pdf'
    #         return this_pi
    #
    #     this_plot_info = DefPlotInfo(plot_info)
    #     data_plot = jpl.Plotting(plot_info=this_plot_info)
    #     fit_is_dict = isinstance(fitDict,(OrderedDict,dict))
    #     varcount = 0
    #     if tflowlist == 'PreDef': tflowlist = [tflowlist]*len(self.SetC2.keys())
    #     if parlist == 'First': parlist = [parlist]*len(self.SetC2.keys())
    #     for icount,(ikey,iC2) in enumerate(self.SetC2.iteritems()):
    #         if (not OnlyVar) or ikey in  self.VarKeys:
    #             itflow = tflowlist[icount]
    #             ipar = parlist[icount]
    #             if fit_is_dict: iindex = ikey
    #             else: iindex = icount
    #             thisxlim = [0,int(fitDict[iindex][0].split('-')[-1])]
    #             if not No_Fit:
    #                 thiscol = colourset8[varcount%len(colourset8)]
    #                 thissym = markerset[varcount%len(markerset)]
    #                 thisshift = shiftset[varcount%len(shiftset)]
    #                 varcount += 1
    #                 print 'Alpha plotting over time ', iC2
    #                 # thiscol,thisshift,thissym = self.GetColShiftSym(iC2,mult=ictflow)
    #                 data_plot = iC2.AlphaPlotImproveFit(data_plot,fitDict[iindex][1],xlims=thisxlim,
    #                                 thistflowlist = itflow,momlist=momlist[icount],WipeFit=WipeFit,
    #                                 thiscol=thiscol,thissym=thissym,thisshift=thisshift,thispar=ipar)
    #                 if (not fit_is_dict) or ikey in fitDict.keys():
    #                     data_plot = iC2.AlphaFitPlot(data_plot,fitDict[iindex],thistflowlist=itflow,
    #                                 momlist=momlist[icount],WipeFit=WipeFit,
    #                                 thiscol=thiscol,thisshift=thisshift,thispar=ipar)
    #                 if ipar == 'First':
    #                     thiscol = colourset8[varcount%len(colourset8)]
    #                     thissym = markerset[varcount%len(markerset)]
    #                     thisshift = shiftset[varcount%len(shiftset)]
    #                     varcount += 1
    #                     if isinstance(itflow,(list,tuple,np.ndarray)):
    #                         int_itflow = itflow[0]
    #                     else:
    #                         int_itflow = itflow
    #                     if isinstance(momlist[icount],(list,tuple,np.ndarray)):
    #                         this_mom = momlist[icount][0]
    #                     data_plot = iC2.PlotIntTsink(data_plot,this_mom,int_itflow,xlims=thisxlim,
    #                                     thiscol=thiscol,thissym=thissym,thisshift=thisshift)
    #     data_plot.PlotAll()
    #     return data_plot


    ## this is used for plotting multiple versions of the same results (over tflow or tsink if you like)
    def GetColShiftSym(self,iC2,mult=0):
        thiscol = colourset8[(colourset8.index(iC2.thiscol)+self.SetC2len*mult)%len(colourset8)]
        thissym = markerset[(markerset.index(iC2.thissym)+self.SetC2len*mult)%len(markerset)]
        thisshift = shiftset[(shiftset.index(iC2.thisshift)+self.SetC2len*mult)%len(shiftset)]
        return thiscol,thisshift,thissym


    def AlphaPlotTflow( self,plot_info,tflowlist=defxlimAlphaTflow,tlist='PreDef',
                        tsum_list='PreDef',momlist=[],OnlyVar=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Nucleon Mixing Angle $\alpha$ over $t_{f}$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$\sqrt{8t_{f}}[fm]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\alpha$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Alpha_tflow.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        varcount = 0
        if tflowlist == 'PreDef': tflowlist = [tflowlist]*len(list(self.SetC2.keys()))
        if tsum_list == 'PreDef': tsum_list = [tsum_list]*len(list(self.SetC2.keys()))
        for icount,(ikey,iC2) in enumerate(self.SetC2.items()):
            if (not OnlyVar) or ikey in  self.VarKeys:
                thiscol = colourset8[varcount%len(colourset8)]
                thissym = markerset[varcount%len(markerset)]
                thisshift = shiftset[varcount%len(shiftset)]
                varcount += 1
                print('Alpha plotting over flow time ', iC2)
                it_sum = tsum_list[icount]
                thist = tlist[icount]
                if isinstance(tflowlist[icount],(list,np.ndarray)):
                    data_plot = iC2.AlphaPlotTflow(data_plot,tflowplot=tflowlist[icount],
                                        thistsinklist = thist,tsum_list=it_sum,momlist=momlist[icount],
                                        thiscol=thiscol,thissym=thissym,thisshift=thisshift,Auto=self.AutoList[icount])
                else:
                    data_plot = iC2.AlphaPlotTflow(data_plot,tflowplot=tflowlist,thistsinklist = thist,
                                        tsum_list=it_sum,momlist=momlist[icount],
                                        thiscol=thiscol,thissym=thissym,thisshift=thisshift,Auto=self.AutoList[icount])
        data_plot.PlotAll()
        return data_plot


    def AlphaPlotTsum(  self,plot_info,tlist='PreDef',tflowlist='PreDef',momlist=[],
                        spec_tsum_list='PreDef',fit_list='None',
                        OnlyVar=False,WipeFit=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Nucleon Mixing Angle $\alpha$ over $t_{sum}$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] =r'$t_{s}$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\alpha$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Alpha_tsum.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        varcount = 0
        if tlist == 'PreDef': tlist = [tlist]*len(list(self.SetC2.keys()))
        if spec_tsum_list == 'PreDef': spec_tsum_list = [spec_tsum_list]*len(list(self.SetC2.keys()))
        if tflowlist == 'PreDef': tflowlist = [tflowlist]*len(list(self.SetC2.keys()))
        if fit_list == 'None':fit_list = ['None']*len(list(self.SetC2.keys()))
        for icount,(ikey,iC2) in enumerate(self.SetC2.items()):
            if (not OnlyVar) or ikey in  self.VarKeys:
                thiscol = colourset8[varcount%len(colourset8)]
                thissym = markerset[varcount%len(markerset)]
                thisshift = shiftset[varcount%len(shiftset)]
                varcount += 1
                print('Alpha plotting over tsum time ', iC2)
                jt_sum = spec_tsum_list[icount]
                itflow = tflowlist[icount]
                thist = tlist[icount]
                thisfitr = fit_list[icount]
                if isinstance(thisfitr[0],str):
                    data_plot = iC2.AlphaPlotTsum(data_plot,
                                        spec_tsum_list=jt_sum,
                                        thistsinklist = thist,thistflowlist=itflow,
                                        fit_list=thisfitr,momlist=momlist[icount],
                                        Auto=self.AutoList[icount],
                                        thiscol=thiscol,thissym=thissym,
                                        thisshift=thisshift,WipeFit=WipeFit)
                else:
                    data_plot = iC2.AlphaPlotTsum(data_plot,
                                        spec_tsum_list=jt_sum,
                                        thistsinklist = thist,thistflowlist=itflow,
                                        fit_list=thisfitr[0],momlist=momlist[icount],
                                        Auto=self.AutoList[icount],
                                        thiscol=thiscol,thissym=thissym,
                                        thisshift=thisshift,WipeFit=WipeFit)
                    if not self.AutoList[icount]:
                        ssfitr = [ifit[0] for ifit in thisfitr]
                        tsumfitr = [ifit[1] for ifit in thisfitr]
                        iC2.FitAlpha(   fit_ranges=ssfitr,tsum_ranges=tsumfitr,
                                        thistflowlist = itflow,NoProd=True,
                                        show_timer=True,WipeFit=WipeFit)
        data_plot.PlotAll()
        return data_plot


    def AlphaPlotCovar(  self,plot_info,tlist='PreDef',tflowlist='PreDef',momlist=[],
                        spec_tsum_list='PreDef',OnlyVar=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Covariance of Nucleon Mixing Angle $\alpha$ over $t_{sum}$'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] =r'$t_{s}$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$Covar$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Alpha_covar.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        varcount = 0
        if tlist == 'PreDef': tlist = [tlist]*len(list(self.SetC2.keys()))
        if spec_tsum_list == 'PreDef': spec_tsum_list = [spec_tsum_list]*len(list(self.SetC2.keys()))
        if tflowlist == 'PreDef': tflowlist = [tflowlist]*len(list(self.SetC2.keys()))
        for icount,(ikey,iC2) in enumerate(self.SetC2.items()):
            if (not OnlyVar) or ikey in  self.VarKeys:
                thiscol = colourset8[varcount%len(colourset8)]
                thissym = markerset[varcount%len(markerset)]
                thisshift = shiftset[varcount%len(shiftset)]
                varcount += 1
                jt_sum = spec_tsum_list[icount]
                itflow = tflowlist[icount]
                thist = tlist[icount]
                if self.AutoList[icount]:
                    print('Alpha plotting Covariance over tsum time ', iC2)
                    data_plot = iC2.AlphaPlotCovar(data_plot,
                                        spec_tsum_list=jt_sum,
                                        thistsinklist = thist,thistflowlist=itflow,
                                        momlist=momlist[icount],
                                        thiscol=thiscol,thissym=thissym,
                                        thisshift=thisshift)
        data_plot.PlotAll()
        return data_plot


    def AlphaPlotDist(  self,plot_info,tlist,tflowlist,tsum_list,
                        tsum_spec_list=None,momlist='p000',
                        OnlyVar=False,Auto=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Nucleon Mixing Angle $\alpha$ Bootstrap Distribution'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] =r'$\alpha$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$frequency$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Alpha_dist.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        if not isinstance(tsum_list,(list,tuple,np.ndarray)):
            tsum_list = [tsum_list for _ in self.SetC2.keys()]
        if not isinstance(tflowlist,(list,tuple,np.ndarray)):
            tflowlist = [tflowlist for _ in self.SetC2.keys()]
        if not isinstance(tlist,(list,tuple,np.ndarray)):
            tlist = [tlist for _ in self.SetC2.keys()]
        if not isinstance(momlist,(list,tuple,np.ndarray)):
            momlist = [momlist for _ in self.SetC2.keys()]
        if not isinstance(tsum_spec_list,(list,tuple,np.ndarray)):
            tsum_spec_list = [tsum_spec_list for _ in self.SetC2.keys()]

        varcount = 0
        for icount,(ikey,iC2) in enumerate(self.SetC2.items()):
            if (not OnlyVar) or ikey in  self.VarKeys:
                print('Alpha plotting distribution ', iC2)
                itflow = tflowlist[icount]
                itsum = tsum_list[icount]
                jtsum = tsum_spec_list[icount]
                it = tlist[icount]
                imom = momlist[icount]
                if isinstance(imom,(list,tuple,np.ndarray)):
                    imom = imom[0]

                if iC2.is_NNQ2sub:
                    NNQkey = (imom,it,itflow,itsum,jtsum)
                    if NNQkey not in iC2.NNQ_subNNQ2_Stats.index:
                        print(iC2.NNQ_subNNQ2_Stats.index)
                        print('Warning: ', NNQkey ,'not in NNQ_subNNQ2_Stats, skipping plot')
                    else:
                        this_col = colourset8[varcount%len(colourset8)]
                        varcount += 1
                        iC2.AlphaPlotDist(data_plot,itflow,itsum,it,spec_tsum=jtsum,mom=imom,
                                thiscol=this_col,Auto=self.AutoList[icount])
                else:
                    NNQkey = (imom,it,itflow,itsum)
                    if NNQkey not in iC2.NNQFull_Stats.index:
                        print(iC2.NNQFull_Stats.index)
                        print('Warning: ', NNQkey ,'not in NNQFull_Stats, skipping plot')
                    else:
                        this_col = colourset8[varcount%len(colourset8)]
                        varcount += 1
                        iC2.AlphaPlotDist(data_plot,itflow,itsum,it,mom=imom,
                                thiscol=this_col,Auto=self.AutoList[icount])
        data_plot.PlotAll()
        return data_plot

    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    ## momlist is now all of them. THIS IS DIFFERENT TO ABOVE
    def AlphaPlotVsEnergy(  self,plot_info,fitDict=defFitDict,momlist=[],
                            tflowlist='PreDef',OnlyVar=False,WipeFit=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'Nucleon Mixing Angle $\alpha$ vs energy extraction from'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] =r'$E_{p}$ [GeV]'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\alpha$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_Alpha_VsEnergy.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        fit_is_dict = isinstance(fitDict,(OrderedDict,dict))
        varcount = 0
        if tflowlist == 'PreDef': tflowlist = [tflowlist]*len(list(self.SetC2.keys()))
        for icount,(ikey,iC2) in enumerate(self.SetC2.items()):
            if (not OnlyVar) or ikey in  self.VarKeys:
                thiscol = colourset8[varcount%len(colourset8)]
                thissym = markerset[varcount%len(markerset)]
                thisshift = shiftset[varcount%len(shiftset)]
                varcount += 1
                print('Alpha plotting over time ', iC2)
                itflow = tflowlist[icount]
                # thiscol,thisshift,thissym = self.GetColShiftSym(iC2,mult=ictflow)
                if fit_is_dict: iindex = ikey
                else: iindex = icount
                data_plot = iC2.AlphaPlotVsEnergy(data_plot,fitDict[iindex],thistflowlist=itflow,
                                        momlist=momlist[icount],thiscol=thiscol,thissym=thissym,
                                        thisshift=thisshift,WipeFit=WipeFit)
        data_plot.PlotAll()
        return data_plot


    def AlphaPlotVsLatSpace(self,plot_info,this_fit,plottflowlist='PreDefine',momlist='p000',OnlyVar=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$\alpha_{N}$ Lattice spacing $a$ dependance'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$ a[fm]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\alpha_{N}$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_VsLatSpace.pdf'
            if 'xlims' in this_pi:
                this_pi['xlims'] = [0,this_pi['xlims'][1]]
            else:
                this_pi['xlims'] = [0,None]
            return this_pi

        if isinstance(this_fit,str):
            this_fit = [this_fit]*len(list(self.SetC2.keys()))
        if isinstance(momlist,str):
            momlist = [momlist]*len(list(self.SetC2.keys()))
        if not isinstance(momlist[0],str):
            momlist = [imom[0] for imom in momlist]
        if isinstance(plottflowlist,str):
            plottflowlist = [plottflowlist]*len(list(self.SetC2.keys()))
        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        for ikey,iset in self.SetC2.items():
            if (not OnlyVar) or ikey in  self.VarKeys:
                firstset = iset
                break
        the_mpi = ''
        this_bool = False
        for ikey,iset in self.SetC2.items():
            if (not OnlyVar) or ikey in  self.VarKeys:
                if firstset.NN.latparams.GetPionMass() != iset.NN.latparams.GetPionMass():
                    this_bool = True
                    print('warning: not all mpi values are the same for latspace plot')
                    break
        if this_bool:
            for ikey,iset in self.SetC2.items():
                if (not OnlyVar) or ikey in  self.VarKeys:
                    print(iset.name)
                    print(iset.NN.latparams.GetPionMass())
                    print()
        else:
            the_mpi = firstset.NN.latparams.GetPionMassLab(MeV=True)
        print('Plotting fit results over lattice spacing')
        ploty,plotx = [],[]
        for ic,(ikey,iset) in enumerate(self.SetC2.items()):
            if (not OnlyVar) or ikey in  self.VarKeys:
                this_fitr,this_tsum_fitr = this_fit[ic].split('_')
                fit_data = iset.Get_Extrapolation()
                if len(fit_data.index) == 0: continue
                this_fun = fit_data.index[0][2]
                this_par = fit_data.index[0][5]
                this_key = momlist[ic],plottflowlist[ic],this_fun,this_fitr,this_tsum_fitr,this_par
                if 'all' in this_key[0].lower():
                    this_key = list(this_key)
                    this_key[0] = fit_data.index[0][0]
                    this_key = tuple(this_key)
                if this_key not in fit_data:
                    print(this_key, ' fit has not been computed for alpha, attempting now')
                    iset.FitAlpha(  fit_range=this_fitr,thistflowlist=[plottflowlist[ic]],
                                    tsum_ranges=this_tsum_fitr)
                    fit_data = iset.Get_Extrapolation()
                    this_fun = fit_data.index[0][2]
                    this_par = fit_data.index[0][5]
                    this_key = momlist[ic],plottflowlist[ic],this_fun,this_fit[ic],this_par
                    if 'all' in this_key[0].lower():
                        this_key = list(this_key)
                        this_key[0] = fit_data.index[0][0]
                        this_key = tuple(this_key)
                    if this_key not in fit_data:
                        print(' ERROR FITTING:' ,  this_key)
                        print(fit_data)
                        continue
                this_names = ['$a$']+list(fit_data.index.names)
                for ireskey,ires in fit_data.items():
                    ploty.append(ires)
                    plotx.append((iset.NN.latparams.GetLatSpaceLab(),)+ireskey) ## in MeV
        if len(plotx) > 0:
            plot_key = [slice(None)] + list(plotx[0][1:])
            indicies = pa.MultiIndex.from_tuples(plotx,names=this_names)
            ploty = pa.Series(ploty,index=indicies)
        else:
            print('No Sets Found for plotting vs lat spacing')
            return data_plot

        # thisFF = thisFF.replace('tilde',r'\tilde')
        # thisFF = thisFF.replace('frac',r'\frac')
        thislab = r'$'+r'\ '.join([the_mpi])+r'$ '
        hold_series = copy(jpl.null_series)
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = ploty.apply(lambda x : x.Avg)
        hold_series['yerr_data'] = ploty.apply(lambda x : x.Avg)
        hold_series['key_select'] = plot_key
        hold_series['type'] = 'error_bar_vary'
        hold_series['label'] = thislab
        hold_series['symbol'] = 'o'
        data_plot.AppendData(hold_series)
        data_plot = self.Plot_Extrap(data_plot,'latspace')
        data_plot.PlotAll()
        return data_plot


    def AlphaPlotVsMpi(self,plot_info,this_fit,plottflowlist='PreDefine',momlist='p000',OnlyVar=False):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = r'$\alpha_{N}$ $m_{\pi}$ dependance'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$ m_{\pi}[MeV]$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel'] = r'$\alpha_{N}$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_VsMpi.pdf'
            if 'xlims' in this_pi:
                this_pi['xlims'] = [0,this_pi['xlims'][1]]
            else:
                this_pi['xlims'] = [0,None]
            return this_pi

        if isinstance(this_fit,str):
            this_fit = [this_fit]*len(list(self.SetC2.keys()))
        if isinstance(momlist,str):
            momlist = [momlist]*len(list(self.SetC2.keys()))
        if not isinstance(momlist[0],str):
            momlist = [imom[0] for imom in momlist]
        if isinstance(plottflowlist,str):
            plottflowlist = [plottflowlist]*len(list(self.SetC2.keys()))
        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        for ikey,iset in self.SetC2.items():
            if (not OnlyVar) or ikey in  self.VarKeys:
                firstset = iset
                break
        the_latspace = ''
        this_bool = False
        for ikey,iset in self.SetC2.items():
            if (not OnlyVar) or ikey in  self.VarKeys:
                if firstset.NN.latparams.latspace != iset.NN.latparams.latspace:
                    this_bool = True
                    print('warning: not all lattice spacing values are the same for mpi plot')
                    break
        if this_bool:
            for ikey,iset in self.SetC2.items():
                if (not OnlyVar) or ikey in  self.VarKeys:
                    print(iset.name)
                    print(iset.NN.latparams.latspace)
                    print()
        else:
            the_latspace = firstset.NN.latparams.GetLatSpaceLab()
        print('Plotting fit results over pion mass')

        ploty,plotx = [],[]
        for ic,(ikey,iset) in enumerate(self.SetC2.items()):
            if (not OnlyVar) or ikey in  self.VarKeys:
                this_fitr,this_tsum_fitr = this_fit[ic].split('_')
                fit_data = iset.Get_Extrapolation()
                if len(fit_data.index) == 0: continue
                this_fun = fit_data.index[0][2]
                this_par = fit_data.index[0][5]
                this_key = momlist[ic],plottflowlist[ic],this_fun,this_fitr,this_tsum_fitr,this_par
                if 'all' in this_key[0].lower():
                    this_key = list(this_key)
                    this_key[0] = fit_data.index[0][0]
                    this_key = tuple(this_key)
                if this_key not in fit_data:
                    print(this_key, ' fit has not been computed for alpha, attempting now')
                    iset.FitAlpha(  fit_range=this_fitr,thistflowlist=[plottflowlist[ic]],
                                    tsum_ranges=this_tsum_fitr)
                    fit_data = iset.Get_Extrapolation()
                    if len(fit_data.index) == 0: continue
                    this_fun = fit_data.index[0][2]
                    this_par = fit_data.index[0][5]
                    this_key = momlist[ic],plottflowlist[ic],this_fun,this_fit[ic],this_par
                    if 'all' in this_key[0].lower():
                        this_key = list(this_key)
                        this_key[0] = fit_data.index[0][0]
                        this_key = tuple(this_key)
                    if this_key not in fit_data:
                        print(' ERROR FITTING:' ,  this_key)
                        print(fit_data)
                        continue
                this_names = ['$m_{\pi}$']+list(fit_data.index.names)
                for ireskey,ires in fit_data.items():
                    ploty.append(ires)
                    plotx.append((iset.NN.latparams.GetPionMassLab(MeV=True),)+ireskey) ## in MeV
        if len(plotx) > 0:
            plot_key = [slice(None)] + list(plotx[0][1:])
            indicies = pa.MultiIndex.from_tuples(plotx,names=this_names)
            ploty = pa.Series(ploty,index=indicies)
        else:
            print('No Sets Found for plotting vs mpi')
            return data_plot

        # thisFF = thisFF.replace('tilde',r'\tilde')
        # thisFF = thisFF.replace('frac',r'\frac')
        thislab = r'$'+r'\ '.join([the_latspace])+r'$ '
        hold_series = copy(jpl.null_series)
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = ploty.apply(lambda x : x.Avg)
        hold_series['yerr_data'] = ploty.apply(lambda x : x.Avg)
        hold_series['key_select'] = plot_key
        hold_series['type'] = 'error_bar_vary'
        hold_series['label'] = thislab
        hold_series['symbol'] = 'o'
        data_plot.AppendData(hold_series)
        data_plot = self.Plot_Extrap(data_plot,'chiral')
        data_plot.PlotAll()
        return data_plot


    def Plot_Extrap(self,plot_data,extrap_type,this_fun='Default'):
        extrap_type = extrap_type.lower()
        this_eval = 0
        if 'latspace' in extrap_type:
            if isinstance(this_fun,str):
                this_fun = (SquaredFitFun,2)
            fit_data = self.LatSpace_Extrapolation(this_fun=this_fun)
        # elif 'boxsize' in extrap_type:
        #     if isinstance(this_fun,str):
        #         this_fun = (SquaredFitFun,2)
        #     fit_data = self.BoxSize_Extrapolation(this_fun=this_fun)
        elif 'mpi' in extrap_type or 'chiral' in extrap_type:
            if isinstance(this_fun,str):
                this_fun = (LinearFitFun,2)
            fit_data = self.Chiral_Extrapolation(this_fun=this_fun)
            this_eval = physcial_pion_mass
        else:
            raise IOError('extrapolation type not recognised '+str(extrap_type))
        if len(fit_data) == 0:
            return plot_data
        # if isinstance(fitr_params,pa.Series): fitr_params = fitr_params.iloc[-1]
        ## TODO, plot from zero to data max
        hold_series = copy(jpl.null_series)
        hold_series['key_select'] = fit_data.index[0]
        hold_series['type'] = 'fit_vary'
        hold_series['fit_class'] = fit_data
        hold_series['label'] = 'fit '
        hold_series['xdatarange'] = [0,'ToDataMax']
        hold_series['ShowEval'] = this_eval
        plot_data.AppendData(hold_series)
        return plot_data


    def LatSpace_Extrapolation(self,this_fun=(SquaredFitFun,2)):
        # print 'DEBUGGING LatSpace_Extrapolation'
        def this_fit(this_series):
            xdata,ydata = [[]],[]
            for icol,coldata in this_series.items():
                if '_latspace' in icol:
                    xdata[0].append(coldata)
                    ##### DEBUGGING ####
                    # xdata[0].append(coldata+100)
                else:
                    ydata.append(coldata)
                    ##### DEBUGGING ####
                    # ydata.append(coldata+200)
                    # ydata[-1].Stats()
            thisff = ff.Fitting(Funs=this_fun,data=[np.array(xdata),np.array(ydata)],name='LatSpace_Ext')
            fit_result = thisff.FitBoots()
            if np.isnan(fit_result):
                return fit_result
            return thisff

        this_df = pa.DataFrame()
        for ikey,iset in self.SetC2.items():
            this_data = iset.Get_Extrapolation()
            if len(this_data) > 0:
                this_df[ikey] = this_data
                this_df.loc[:,ikey+'_latspace'] = float(iset.NN.latparams.latspace)
        this_df.loc[:,'Fit'] = this_df.apply(this_fit,axis=1)
        return this_df['Fit'].dropna()

    def Chiral_Extrapolation(self,this_fun=(LinearFitFun,2)):
        # print 'DEBUGGING Chiral_Extrapolation'
        def this_fit(this_series):
            xdata,ydata = [[]],[]
            for icol,coldata in this_series.items():
                if '_chiral' in icol:
                    xdata[0].append(coldata)
                    # ##### DEBUGGING ####
                    # xdata[0].append(coldata+100)
                else:
                    ydata.append(coldata)
                    # ##### DEBUGGING ####
                    # ydata.append(coldata+200)
                    # ydata[-1].Stats()
            thisff = ff.Fitting(Funs=this_fun,data=[np.array(xdata),np.array(ydata)],name='Chiral_Ext')
            fit_result = thisff.FitBoots()
            if np.isnan(fit_result):
                return fit_result
            return thisff

        this_df = pa.DataFrame()
        for ikey,iset in self.SetC2.items():
            this_data = iset.Get_Extrapolation()
            if len(this_data) > 0:
                this_df[ikey] = this_data
                this_df.loc[:,ikey+'_chiral'] = float(iset.NN.latparams.GetPionMass(MeV=True))
        this_df.loc[:,'Fit'] = this_df.apply(this_fit,axis=1)
        return this_df['Fit'].dropna()


    def __setitem__(self, key, value ):
        self.SetC2[key[0]][key[1:]] = value

    def __getitem__(self, key):
        return self.SetC2[key[0]][key[1:]]

    def Stats(self):
        for iset in self.SetC2.values():
            iset.Stats()

    def __str__(self):
        return self.name

    def __iter__(self):
        return self

    def items(self):
        outdata = []
        for iset,setdata in self.SetC2.items():
            for ikey,tdata in list(setdata.items()):
                outdata.append((iset,)+ikey+(tdata,))
        return outdata

    def itemsAvgStd(self):
        outdata = []
        for iset,setdata in self.SetC2.items():
            for ikey,tdata in list(setdata.items()):
                outdata.append((iset,)+ikey+(tdata.Avg,tdata.Std))
        return outdata

    def values(self):
        outdata = []
        for setdata in self.SetC2.values():
            for tdata in list(setdata.values()):
                outdata.append(tdata)
        return outdata



    def valuesAvgStd(self):
        outdata = []
        for setdata in self.SetC2.values():
            for tdata in list(setdata.values()):
                outdata.append((tdata.Avg,tdata.Std))
        return outdata

    def keys(self):
        outdata = []
        for iset,setdata in self.SetC2.items():
            for ikey in list(setdata.keys()):
                outdata.append((iset,)+ikey)
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
            np.complex(self[ikey])
            self[ikey].Stats()
        return 1+0j

    def __int__(self):
        for ikey in list(self.keys()):
            int(self[ikey])
            self[ikey].Stats()
        return 1

    def __float__(self):
        for ikey in list(self.keys()):
            np.float(self[ikey])
            self[ikey].Stats()
        return 1.0
