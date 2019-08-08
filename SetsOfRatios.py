#!/usr/bin/env python

import numpy as np

# import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
# import matplotlib.pyplot as pl

from collections import OrderedDict
# from MiscFuns import NumbToFileName,mkdir_p,ODNested
from MiscFuns import mkdir_p,CombineListCfgs
from Params import graphdir,defFitDict,defMcore,nchi_threshold
import pandas as pa
# from XmlFormatting import rm_tsumfitr

from RatioCorrelators import RatioCorr,RatioFOCorr,RatioFOFullCorr
from TwoPtCorrelators import TwoPointCorr
from MultiWrap import DoMulticore
from TimeStuff import Timer
from copy import copy
from FileIO import Construct_Empty
import PlotData as jpl
import operator as op

def Mcore_CombRats(RatSet):
    if len(RatSet) == 0:
        return '','None Found'
    this_key = RatSet[0][0]
    this_set = RatSet[0][1]
    this_name = this_set.name
    for iset_key,iset_var,coeff_scale in RatSet[1:]:
        this_set = this_set + (iset_var*coeff_scale)
    this_set.ReduceIndicies()
    this_set.Stats()
    this_set.name = this_name
    return this_set,this_key


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



class SetOfRat(object):
    """

    Rat(set) Inherits from RatioCorrelators Class

    plotting parameters (Move to different file?) ####

    InfoDict is a list of Info dictionaries to be bassed to RatioCorrelators:
    Setup to do fits taken from dictionaries within InfoDict (see RatioCorrelators.py)
    Although Info['gammas'] can have multiple gamma matricies,
    This code is only set up to take a list, where each list only has one gamma matrix in Info['gammas'] ( e.g. = ['g4'])




    """

    def Construct_Rat_Set(this_file_list,fo_list='Def'):
        if fo_list == 'Def':
            fo_list = ['' for i in this_file_list]
        set_list = []
        for ifull,ifile in zip(full_list,this_file_list):
            if 'FO' in ifull:
                if 'Full' in ifull:
                    set_list.append(RatioFOFullCorr.Construct_RatFOFull_File(ifile))
                else:
                    set_list.append(RatioFOCorr.Construct_RatFO_File(ifile))
            else:
                set_list.append(RatioCorr.Construct_Rat_File(ifile))
        out_class = Construct_Empty(SetOfRat)
        out_class.cfglist = {}
        out_class.InfoDict = []
        out_class.Plotdir = graphdir+'/Rat/'
        out_class.name = 'Set1'
        out_class.SetRat = {}
        for iset in set_list:
            out_class.SetRat[iset.name] = iset
        return out_class

    instances = []

    parentlist = []

    def __init__(self, cfglist={}, InfoDict=[],corr_cfglists=True,name = '',setnames = [],thisFO='FromList'):

        self.cfglist = cfglist
        self.InfoDict = InfoDict
        self.Plotdir = graphdir+'/Rat/'
        mkdir_p(self.Plotdir)


        if len(name) == 0:
            self.name = 'Set'+str(int(len(SetOfRat.instances))+1)
            SetOfRat.instances.append(self.name)
        else:
            self.name = name
            # SetOfRat.instances.append(self.name)

        if thisFO == 'FromList':
            self.thisFO = False
            for idict in self.InfoDict:
                if 'Observable' in list(idict.keys()):
                    if idict['Observable'] != '':
                        self.thisFO = idict['Observable']
                        break
        else:
            self.thisFO = thisFO

        for idict in self.InfoDict: idict['fo_for_cfgs'] = self.thisFO

        ## dictionary containing RatioCorrs
        self.SetRat = OrderedDict()
        if len(setnames) == 0:
            setnames = ['']*len(InfoDict)
        thistimer = Timer(linklist=setnames,name='Ratio Set '+ self.name)
        for ic,(iname,idict) in enumerate(zip(setnames,self.InfoDict)):
            thiscol = colourset8[ic%len(colourset8)]
            thissym = markerset[ic%len(markerset)]
            thisshift = shiftset[ic%len(shiftset)]
            ## each flowed observable dictionary must contain the flow op in Info['Observable']
            if 'Observable' in list(idict.keys()):
                if idict['Observable'] == '':
                    thisRat = RatioCorr(    cfglist=self.cfglist, Info=idict,thissym=thissym,
                                            thiscol=thiscol,thisshift=thisshift,name = iname,filename=iname)
                elif 'Full_Rat' in idict and idict['Full_Rat']:
                    thisRat = RatioFOFullCorr(  cfglist=self.cfglist, Info=idict,thissym=thissym,
                                                thiscol=thiscol,thisshift=thisshift,name = iname,filename=iname)
                else:
                    thisRat = RatioFOCorr(  cfglist=self.cfglist, Info=idict,thissym=thissym,
                                            thiscol=thiscol,thisshift=thisshift,name = iname,filename=iname)
            else:
                thisRat = RatioCorr(cfglist=self.cfglist, Info=idict,thissym=thissym,
                                    thiscol=thiscol,thisshift=thisshift,name = iname,filename=iname)
            # if isinstance(self.cfglist,dict) and len(self.cfglist.keys()) == 0: self.cfglist = thisRat.cfglist

            if str(thisRat) in list(self.SetRat.keys()):
                Ratnext = 2
                while str(thisRat)+'_'+str(Ratnext) in list(self.SetRat.keys()):
                    Ratnext += 1
                self.SetRat[str(thisRat)+'_'+str(Ratnext)] = thisRat
            else:
                self.SetRat[str(thisRat)] = thisRat
            thistimer.Lap()
            # print 'Setting up classes for ' , thisRat , ' complete'
        # print
            # print 'has color',thiscol , 'symbol ' , thissym , 'thisshift',thisshift

        self.SetRatlen = len(self.SetRat)
        if len(self.cfglist) == 0:
            if corr_cfglists:
                self.GetCfgList()
            else:
                self.GetUnCorrCfgList()

        ##Implement Reading cfglist if len(cfglist) == 0

        self.current = 0

    def GetUnCorrCfgList(self):
        for iset,iRat in self.SetRat.items():
            iRat.GetCfgList()

    def GetCfgList(self):
        self.GetUnCorrCfgList()
        cfglist = CombineListCfgs(np.array(list(self.SetRat.values())).flatten(),'cfglist')
        self.ImportCfgList(cfglist)


    def ImportCfgList(self,cfglist):
        self.cfglist = cfglist
        for iRat in self.SetRat.values():
            iRat.ImportCfgList(cfglist)


    def CombineWithAlpha(self,alpha_in,tflow_list,this_fun=op.truediv):
        for iRat in self.SetRat.values():
            if self.CheckFO(iRat):
                iRat.CombineWithAlpha(alpha_in,this_fun=this_fun)
            else:
                iRat.CombineWithAlpha(alpha_in,tflow_list,this_fun=this_fun)


    def DoAlphaCombine(self,this_fun=op.truediv):
        for iRat in self.SetRat.values():
            if self.CheckFO(iRat):
                iRat.DoAlphaCombine(this_fun=this_fun)

    def Get_CombRats(self,this_keys,this_coeffs):
        this_len = float(len(this_coeffs))
        set_list,coeff_list = [],[]
        for ikey,icoeff in zip(this_keys,this_coeffs):
            coeff_scale = icoeff/this_len
            coeff_list.append(coeff_scale)
            for iset_key,iset_var in self.SetRat.items():
                if ikey in iset_key:
                    set_list.append((iset_key,iset_var,coeff_scale))
        return set_list

    def CombRats(self,this_keys,this_coeffs,show_timer=None):
        first = True
        this_name,this_key = '',''
        this_set = 'None Found'
        this_len = float(len(this_coeffs))
        for ikey,icoeff in zip(this_keys,this_coeffs):
            coeff_scale = icoeff/this_len
            for iset_key,iset_var in self.SetRat.items():
                if ikey in iset_key:
                    if first:
                        this_set = iset_var*coeff_scale
                        first = False
                        this_name = iset_var.name
                        this_key = iset_key
                    else:
                        this_set = this_set + (iset_var*coeff_scale)
            if show_timer is not None:
                show_timer.Lap(ikey)
                    # print 'DEBUG'
                    # print iset_var.name
                    # print type(this_set.Rat_Stats['boot'].iloc[0])
                    # print type(this_set.Rat_Stats['boot'].iloc[0].bootvals[0])
                    # print this_set.Rat_Stats['boot'].iloc[0].bootvals[0]
                    # print
        if this_set != 'None Found':
            this_set.ReduceIndicies()
            this_set.Stats()
            this_set.name = this_name
        return this_set,this_key

    def ReduceRatios(self,key_list,coeff_list,show_timer=False,Mcore=defMcore):
        hold_ratios = OrderedDict()
        thistimer = None
        if Mcore:
            this_key_llist = [ ['_'.join((jkey[0],jkey[2])) for jkey in ikey] for ikey in key_list]
            vallist = [ (ikey,icoeff) for ikey,icoeff in zip(this_key_llist,coeff_list)]
            if show_timer:
                print('Getting Data for Mcore')
            vallist = [self.Get_CombRats(*ival) for ival in vallist]
            # for ival in vallist:
            #     print ival
            if show_timer:
                # print 'Reducing Ratios: '
                thistimer = Timer(name=r'Reducing Ratios Mcore: ')
            output = DoMulticore(Mcore_CombRats,vallist)

            for ic,(this_rat,this_key) in enumerate(output):
                if isinstance(this_rat,str) and this_rat == 'None Found':
                    print()
                    print('Warning, no ratio in set for:')
                    print('\n'.join(this_key_llist[ic]))
                    print()
                    this_key = 'Error'
                else:
                    hold_ratios[this_key] = this_rat
            if show_timer:
                thistimer.Stop()
        else:
            fmt_key_list = [['_'.join((jkey[0],jkey[2])) for jkey in ikey] for ikey in key_list]
            if show_timer:
                # print 'Reducing Ratios: '
                thistimer = Timer(linklist=fmt_key_list,name=r'Reducing Ratios: ',prod_or_sum='sum')
            for ikey,icoeff in zip(key_list,coeff_list):
                this_key_list = ['_'.join((jkey[0],jkey[2])) for jkey in ikey]
                this_rat,this_key = self.CombRats(this_key_list,icoeff,show_timer=thistimer)
                if isinstance(this_rat,str) and this_rat == 'None Found':
                    print()
                    print('Warning, no ratio in set for:')
                    print('\n'.join(this_key_list))
                    print()
                    this_key = 'Error'
                else:
                    hold_ratios[this_key] = this_rat
                # if show_timer: print 'Ratio ',this_key,' complete.'
        self.SetRat = hold_ratios
        self.SetRatlen = len(self.SetRat)

    def Read(self):
        for setdata in self.SetRat.values():
            print('Importing Ratio of corr: ' , str(setdata))
            setdata.Read()
        print()

    def Write(self):
        for iRat in self.SetRat.values():
            print('Writing Ratio of corr: ' , str(iRat))
            iRat.Write()
        print()

    def Reformat3pt(self,thisdir):
        print('reformatting 3 point correlation functions into ' , thisdir)
        for iRat in self.SetRat.values():
            print('reformatting ',str(iRat))
            iRat.C3class.Read(show_timer = True)
            iRat.C3class.SetCustomName()
            if not iRat.__class__.__name__ == 'RatioFOFullCorr':
                iRat.C3class.WriteCfgsToFile(thisdir)

    def LoadPickleCommonC2(self,DefWipe=False,WipeData=True,CheckCfgs=False,Mcore = defMcore):
        if Mcore:
            print('multicore has loadpicklecommonc2 not implemented yet, running over single thread')
            # raise NotImplementedError('multicore has loadpicklecommonc2 not implemented yet')
        # else:
        if len(list(self.SetRat.keys())) > 0:
            thistimer = Timer(linklist=list(self.SetRat.keys()),name=self.name)
            thisiter = iter(self.SetRat.items())
            first_key,first_Rat = next(thisiter)
            first_Rat.LoadPickle(DefWipe=DefWipe,WipeData=WipeData,CheckCfgs=CheckCfgs)
            thistimer.Lap(first_Rat.name)
            if isinstance(first_Rat.C2class,TwoPointCorr):
                common_NN = first_Rat.C2class
                common_NNQ = None
            else:
                common_NN = None
                common_NNQ = first_Rat.C2class
            for ikey,iRat in thisiter:
                isFO = self.CheckFO(iRat)
                isFOFull = self.CheckFOFull(iRat)
                if isFO and not isFOFull:
                    if common_NNQ is not None:
                        iRat.C2class = common_NNQ
                else:
                    if common_NN is not None:
                        iRat.C2class = common_NN
                iRat.LoadPickle(DefWipe=DefWipe,WipeData=WipeData,CheckCfgs=CheckCfgs)
                if isFO and not isFOFull:
                    if common_NNQ is None:
                        common_NNQ = iRat.C2class
                else:
                    if common_NN is None:
                        common_NN = iRat.C2class
                thistimer.Lap(iRat.name)

    def LoadPickle(self,DefWipe=False,WipeData=True,CheckCfgs=False,Mcore = defMcore):
        def LPMulti(thisinput):
            iset,setdata = thisinput
            # print 'LoadPickle Ratio of corr: ' , str(setdata)
            thistimer = Timer(name=str(setdata))
            setdata.LoadPickle(DefWipe=DefWipe,WipeData=WipeData,CheckCfgs=CheckCfgs)
            setdata.ImportPlotParams(colourset8[iset%len(colourset8)],markerset[iset%len(markerset)],
                                     shiftset[iset%len(shiftset)])
            thistimer.Stop(newl=False)
            # print 'Ratio of corr: ' , str(setdata) , ' Complete! '
            return iset,setdata


        if Mcore:
            # print LPMulti
            vallist = [(iset,setdata) for iset,setdata in enumerate(self.SetRat.values())]
            # for ival in vallist:
            #     print ival
            output = DoMulticore(LPMulti,vallist)
            for icount,ikey in enumerate(self.SetRat.keys()):
                for icout,iout in output:
                    if icout == icount:
                        self.SetRat[ikey] = iout
        else:
            # for iset,setdata in enumerate(self.SetRat.itervalues()):
            [ LPMulti((iset,setdata)) for iset,setdata in enumerate(self.SetRat.values())]
        # print

    def ReadAndWrite(self,DefWipe=False):
        for setdata in self.SetRat.values():
            print('Importing and Writing Ratio of corr: ' , str(setdata))
            setdata.ReadAndWrite(DefWipe=DefWipe)
        print()

    def GetSet(self,thisset):
        return self.SetRat[thisset]

    def DoFitRanges(self,fit_min,fit_max,tsum_fit='None',min_fitr=3):
        fit_min,fit_max = str(int(fit_min)),str(int(fit_max))
        fitr = 'fitr'+fit_min+'-'+fit_max
        if isinstance(tsum_fit,(list,tuple,np.ndarray)):
            tsum_fit = 'tsumfitr'+str(int(tsum_fit[0]))+'-'+str(int(tsum_fit[1]))
        thistimer = Timer(linklist=list(self.SetRat.keys()),name='Fitting '+ self.name)
        for iset,setdata in self.SetRat.items():
            isFOFull = self.CheckFOFull(setdata)
            if isFOFull and tsum_fit == 'None':
                raise EnvironmentError('tsum_fit_range must be passed in if using FOFull')
            if isFOFull:
                setdata.Fit(fit_range=fitr,tsum_ranges=tsum_fit,WipeFit=False,min_fitr=min_fitr-1)
            else:
                setdata.Fit(fit_range=fitr,WipeFit=False,min_fitr=min_fitr-1)
            thistimer.Lap()
        # return fit_list

    def GetBestFits(self,fit_min,fit_max,thisflowlist='All',tsum_fit='None',min_fitr=3):
        # print 'DEBUG, uncomment for full run'
        print('Undergoing Fits of Ratio Functions')
        self.DoFitRanges(fit_min,fit_max,tsum_fit=tsum_fit,min_fitr=min_fitr)
        print('Selecting Best Fits')
        lfout,lfchi,iflist = [],[],[]
        lout,lchi,ilist = [],[],[]
        def Get_Chi(val):
            return val.fit_data['Chi2DoF'].iloc[0].Avg
        def Get_Res(val):
            return val.fit_data['Params'].iloc[0]
        for iset,setdata in self.SetRat.items():
            isFO = self.CheckFO(setdata)
            isFOFull = self.CheckFOFull(setdata)
            if isFO or isFOFull:
                for ikey,idata in setdata.Rat_Fit_Stats['boot'].items():
                    fit_result,chi_data = idata.GetBestChi()
                    ifit_key = idata.GetBestChi_key()
                    if isFOFull:
                        ifit_key = '_'.join(ifit_key[1:])
                    else:
                        ifit_key = ifit_key[-1]
                    iflist.append((iset,) + ikey + (ifit_key,))
                    lfout.append(fit_result.fit_data['Params'].iloc[0])
                    lfchi.append(chi_data)
            else:
                for ikey,idata in setdata.Rat_Fit_Stats['boot'].items():
                    fit_result,chi_data = idata.GetBestChi()
                    ifit_key = idata.GetBestChi_key()
                    ilist.append((iset,) + ikey + (ifit_key[-1],))
                    lout.append(fit_result.fit_data['Params'].iloc[0])
                    lchi.append(chi_data)
        output,outputflow = pa.DataFrame(),pa.DataFrame()
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,
                            names=[  'Rat_Set','current_operator','momentum','time_fit_range'])
            output.loc[:,'boot'] = pa.Series(lout,index=indicies)
            output.loc[:,'chi2pdf'] = pa.Series(lchi,index=indicies)
            # output.loc[:,'Avg'] = pa.Series(lAvg,index=indicies)
            # output.loc[:,'Std'] = pa.Series(lStd,index=indicies)
        else:
            print('WARNING:no Ratio function fits found')
        if len(iflist) > 0:
            indicies = pa.MultiIndex.from_tuples(iflist,
                        names=[ 'Rat_Set','current_operator','momentum',
                                'flow_times','time_fit_range'])
            outputflow.loc[:,'boot'] = pa.Series(lfout,index=indicies)
            outputflow.loc[:,'chi2pdf'] = pa.Series(lfchi,index=indicies)
            # outputflow.loc[:,'Avg'] = pa.Series(lfAvg,index=indicies)
            # outputflow.loc[:,'Std'] = pa.Series(lfStd,index=indicies)
        elif isFO or isFOFull:
            print('WARNING:, no flow time Ratio function fits found')
        return output,outputflow


    def GetFitsFromIndex( self,input_indicies,input_flowed_indicies,
                            fit_min,fit_max,thisflowlist='All',tsum_fit='None',
                            chi_min=0.3,chi_max=0.7,min_fitr=3):
        print('Undergoing Fits of Ratio Functions')
        def Get_Res(val):
            return val.fit_data['Params'].iloc[0]
        self.DoFitRanges(fit_min,fit_max,tsum_fit=tsum_fit,min_fitr=min_fitr)
        lfout,iflist = [],[]
        lout,ilist = [],[]
        for iset,setdata in self.SetRat.items():
            isFO = self.CheckFO(setdata)
            isFOFull = self.CheckFOFull(setdata)
            if isFO or isFOFull:
                for ikey,idata in setdata.Rat_Fit_Stats['boot'].items():
                    this_index = input_flowed_indicies[ikey+(slice(None),)]
                    fit_data = idata.Fit_Stats_fmt.iloc[this_index]
                    for ifit,ival in fit_data.items():
                        if isFOFull:
                            ifit = '_'.join(ifit[1:])
                        else:
                            ifit = ifit[-1]
                        iflist.append(ikey + (ifit,))
                        lfout.append(ival.fit_data['Params'].iloc[0])
            else:
                isFO = False
                isFOFull = False
                for ikey,idata in setdata.Rat_Fit_Stats['boot'].items():
                    this_index = input_flowed_indicies[ikey+(slice(None),)]
                    fit_data = idata.Fit_Stats_fmt.iloc[this_index]
                    for ifit,ival in fit_data.items():
                        ilist.append(ikey + (ifit[-1],))
                        lout.append(ival.fit_data['Params'].iloc[0])
        output,outputflow = pa.DataFrame(),pa.DataFrame()
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,
                            names=[  'Rat_Set','current_operator','momentum','time_fit_range'])
            output.loc[:,'boot'] = pa.Series(lout,index=indicies)
        else:
            print('WARNING: no Ratio function fits found')
        if len(iflist) > 0:
            indicies = pa.MultiIndex.from_tuples(iflist,
                        names=[ 'Rat_Set','current_operator','momentum','flow_times',
                                'time_fit_range'])
            outputflow.loc[:,'boot'] = pa.Series(lfout,index=indicies)
        elif isFO or isFOFull:
            print('WARNING:, no flow time Ratio function fits found')
        return output,outputflow



    def GetChi2PdfFits( self,fit_min,fit_max,thisflowlist='All',tsum_fit='None',ralpha_info=(None,None,'boot'),
                        chi_min=0.3,chi_max=0.7,min_fitr=3,nchi_threshold=nchi_threshold,imp_list='All'):
        this_col = ralpha_info[2]
        if this_col != 'boot':
            print('Constructing alpha combine')
            self.CombineWithAlpha(ralpha_info[0],thisflowlist[0],this_fun=ralpha_info[1])
        print('Undergoing Fits of Ratio Functions')
        self.DoFitRanges(fit_min,fit_max,tsum_fit=tsum_fit,min_fitr=min_fitr)
        print('Selecting Chi2pdf in range', chi_min,'-',chi_max)
        lfout,lfchi,iflist = [],[],[]
        lout,lchi,ilist = [],[],[]
        imp_check = isinstance(imp_list,str)
        for iset,setdata in self.SetRat.items():
            isFO = self.CheckFO(setdata)
            isFOFull = self.CheckFOFull(setdata)
            if isFO or isFOFull:
                for ikey,idata in setdata.Rat_Fit_Stats[this_col].items():
                    this_nchi = 1
                    if imp_check or setdata.C3FOclass.NJN.Proj+'_'+ikey[0] in imp_list:
                        this_nchi = nchi_threshold
                    fit_data,chi_test = idata.GetChiRange(chi_min=chi_min,chi_max=chi_max,nchi_cutoff = this_nchi)
                    if isinstance(chi_test,float):
                        print(' '.join(['Warning, not fit ranges found for ',iset,
                                'that satisfies chi range, using minimum with chi of ','{:.5f}'.format(chi_test)]))
                    for ifit,ival in fit_data.iterrows():
                        if isFOFull:
                            ifit = '_'.join(ifit[1:])
                        else:
                            ifit = ifit[-1]
                        iflist.append((iset,) + ikey + (ifit,))
                        lfout.append(ival['Fit'].fit_data['Params'].iloc[0])
                        lfchi.append(ival['Chi2pdfAvg'])
            else:
                isFO = False
                isFOFull = False
                for ikey,idata in setdata.Rat_Fit_Stats[this_col].items():
                    this_nchi = 1
                    if imp_check or setdata.C3FOclass.NJN.Proj+'_'+ikey[0] in imp_list:
                        this_nchi = nchi_threshold
                    fit_data,chi_test = idata.GetChiRange(chi_min=chi_min,chi_max=chi_max,nchi_cutoff = this_nchi)
                    if isinstance(chi_test,float):
                        print(' '.join(['Warning, not fit ranges found for ',iset,
                                'that satisfies chi range, using minimum with chi of ','{:.5f}'.format(chi_test)]))
                    for ifit,ival in fit_data.iterrows():
                        ilist.append((iset,) + ikey + (ifit[-1],))
                        lout.append(ival['Fit'].fit_data['Params'].iloc[0])
                        lchi.append(ival['Chi2pdfAvg'])
        output,outputflow = pa.DataFrame(),pa.DataFrame()
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=[  'Rat_Set','current_operator','momentum',
                                                                'time_fit_range'])
            output.loc[:,'boot'] = pa.Series(lout,index=indicies)
            output.loc[:,'chi2pdf'] = pa.Series(lchi,index=indicies)
        else:
            print('WARNING: no Ratio function fits found')
        if len(iflist) > 0:
            indicies = pa.MultiIndex.from_tuples(iflist,names=[ 'Rat_Set','current_operator','momentum','flow_times',
                                                                'time_fit_range'])
            outputflow.loc[:,'boot'] = pa.Series(lfout,index=indicies)
            outputflow.loc[:,'chi2pdf'] = pa.Series(lfchi,index=indicies)
        elif isFO or isFOFull:
            print('WARNING:, no flow time Ratio function fits found')
        return output,outputflow


    ## depreciated
    # def GetFits(self,fitr = 'All',thisflowlist = 'All'):
    #     lfout,iflist = [],[]
    #     lout,ilist = [],[]
    #     ## output { set , igamma , ip  } Fitting class instance
    #     for iset,setdata in self.SetRat.iteritems():
    #         isFO = self.CheckFO(setdata)
    #         isFOFull = self.CheckFOFull(setdata)
    #         if isFOFull:
    #             if len(setdata.Rat_Fit_Stats.index) == 0 or fitr \
    #             not in setdata.Rat_Fit_Stats.index.get_level_values('tau_tsum_fit_ranges'):
    #                 fitr_list = fitr.split('_')
    #                 setdata.Fit(fit_ranges=fitr_list[0],tsum_ranges=fitr_list[1])
    #         else:
    #             if len(setdata.Rat_Fit_Stats.index) == 0 or fitr \
    #             not in setdata.Rat_Fit_Stats.index.get_level_values('time_fit_range'):
    #                 setdata.Fit(fit_ranges=fitr)
    #         if isFO or isFOFull:
    #             for (igamma,ip,itflow,ifitr),idata in setdata.Rat_Fit_Stats['boot'].iteritems():
    #                 fit_check = fitr == 'All' or ifitr == fitr
    #                 flow_check = thisflowlist == 'All' or itflow in thisflowlist
    #                 if fit_check and flow_check:
    #                     fitdata = idata.fit_data['Params'].iloc[0]
    #                     iflist.append((iset,igamma,ip,itflow,ifitr))
    #                     lfout.append(fitdata)
    #                     lfAvg.append(fitdata.Avg)
    #                     lfStd.append(fitdata.Std)
    #         else:
    #             for (igamma,ip,ifitr),idata in setdata.Rat_Fit_Stats['boot'].iteritems():
    #                 fit_check = fitr == 'All' or ifitr == rm_tsumfitr(fitr)
    #                 if fit_check:
    #                     fitdata = idata.fit_data['Params'].iloc[0]
    #                     ilist.append((iset,igamma,ip,ifitr))
    #                     lout.append(fitdata)
    #                     lAvg.append(fitdata.Avg)
    #                     lStd.append(fitdata.Std)
    #     output,outputflow = pa.DataFrame(),pa.DataFrame()
    #     if len(ilist) > 0:
    #         indicies = pa.MultiIndex.from_tuples(ilist,names=[  'Rat_Set','current_operator','momentum',
    #                                                             'time_fit_range'])
    #         output.loc[:,'boot'] = pa.Series(lout,index=indicies)
    #         output.loc[:,'Avg'] = pa.Series(lAvg,index=indicies)
    #         output.loc[:,'Std'] = pa.Series(lStd,index=indicies)
    #     if len(iflist) > 0:
    #         indicies = pa.MultiIndex.from_tuples(iflist,names=[ 'Rat_Set','current_operator','momentum','flow_times',
    #                                                             'time_fit_range'])
    #         outputflow.loc[:,'boot'] = pa.Series(lfout,index=indicies)
    #         outputflow.loc[:,'Avg'] = pa.Series(lfAvg,index=indicies)
    #         outputflow.loc[:,'Std'] = pa.Series(lfStd,index=indicies)
    #     return output,outputflow

    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def Plot(   self,plot_info,fitDict=defFitDict,fitDictFull=defFitDict,gammamomlist={},
                tflowlist = [],tsum_list=[]):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = 'Ratio Funciton plot'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$\tau - t/2$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel']  =r'$Ratio Function$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        for icount,(ikey,iRat) in enumerate(self.SetRat.items()):
            isFO = self.CheckFO(iRat)
            isFOFull = self.CheckFOFull(iRat)
            print('plotting ', iRat, end=' ')
            fit_is_list = isinstance(fitDict,(list,tuple,np.ndarray))
            if isFOFull:
                data_plot = iRat.Plot(data_plot,gammamomlist=gammamomlist[icount],tflowlist=tflowlist,tsum_list=tsum_list)
            elif isFO:
                data_plot = iRat.Plot(data_plot,gammamomlist=gammamomlist[icount],tflowlist=tflowlist)
            else:
                data_plot = iRat.Plot(data_plot,gammamomlist=gammamomlist[icount])
            if fit_is_list:
                print('fit as well ')
                if isFOFull:
                    data_plot = iRat.FitPlot(data_plot,fitDictFull[icount],gammamomlist=gammamomlist[icount],tflowlist=tflowlist,tsum_list=tsum_list)
                elif isFO:
                    data_plot = iRat.FitPlot(data_plot,fitDict[icount],gammamomlist=gammamomlist[icount],tflowlist=tflowlist)
                else:
                    data_plot = iRat.FitPlot(data_plot,fitDict[icount],gammamomlist=gammamomlist[icount])
            else:
                if ikey in list(fitDict.keys()):
                    print('fit as well ')
                    if isFOFull:
                        data_plot = iRat.FitPlot(data_plot,fitDictFull[ikey],gammamomlist=gammamomlist[icount],tflowlist=tflowlist,tsum_list=tsum_list)
                    elif isFO:
                        data_plot = iRat.FitPlot(data_plot,fitDict[ikey],gammamomlist=gammamomlist[icount],tflowlist=tflowlist)
                    else:
                        data_plot = iRat.FitPlot(data_plot,fitDict[ikey],gammamomlist=gammamomlist[icount])
                else:
                    print()
        data_plot.PlotAll()
        return data_plot


    ## fitlist is list of touples { set : ('state#','fitr#-#')}
    def PlotTsum(self,plot_info,fitDict=defFitDict,gammamomlist={},tflowlist = [],tcurr_list=[]):
        def DefPlotInfo(pi):
            this_pi = copy(pi)
            if 'title' not in this_pi:
                this_pi['title'] = 'Ratio Funciton plot'
            if 'xlabel' not in this_pi:
                this_pi['xlabel'] = r'$t_{s}$'
            if 'ylabel' not in this_pi:
                this_pi['ylabel']  =r'$Ratio Function$'
            if 'save_file' not in this_pi:
                this_pi['save_file'] = self.Plotdir+self.name+'_tsum.pdf'
            return this_pi

        this_plot_info = DefPlotInfo(plot_info)
        data_plot = jpl.Plotting(plot_info=this_plot_info)
        for icount,(ikey,iRat) in enumerate(self.SetRat.items()):
            if not self.CheckFOFull(iRat): continue
            print('plotting tsum', iRat, end=' ')
            fit_is_list = isinstance(fitDict,(list,tuple,np.ndarray))
            data_plot = iRat.PlotTsum(data_plot,gammamomlist=gammamomlist[icount],tflowlist=tflowlist,tcurr=tcurr_list)
            if fit_is_list:
                print('fit as well ')
                data_plot = iRat.FitSumPlot(data_plot,fitDict[icount],gammamomlist=gammamomlist[icount],tflowlist=tflowlist,tcurr=tcurr_list)
            else:
                if ikey in list(fitDict.keys()):
                    print('fit as well ')
                    data_plot = iRat.FitSumPlot(data_plot,fitDict[ikey],gammamomlist=gammamomlist[icount],tflowlist=tflowlist,tcurr=tcurr_list)
        data_plot.PlotAll()
        return data_plot


    def CheckFO(self,iRat):
        return iRat.__class__.__name__ == 'RatioFOCorr'

    def CheckFOFull(self,iRat):
        return iRat.__class__.__name__ == 'RatioFOFullCorr'


    def __setitem__(self, key, value ):
        self.SetRat[key[0]][key[1:]] = value

    def __getitem__(self, key):
        return self.SetRat[key[0]][key[1:]]

    def Stats(self):
        for iset in self.SetRat.values():
            iset.Stats()

    def __str__(self):
        return self.name

    def __iter__(self):
        return self

    def items(self):
        outdata = []
        for iset,setdata in self.SetRat.items():
            #values = igamma,ip,(maybe itflow,) it, data
            for values in list(setdata.items()):
                outdata.append((iset,)+values)
        return outdata

    def itemsAvgStd(self):
        outdata = []
        for iset,setdata in self.SetRat.items():
            for values in setdata.itemsAvgStd():
                outdata.append((iset,) + values)
        return outdata

    def values(self):
        outdata = []
        for setdata in self.SetRat.values():
            for tdata in list(setdata.values()):
                outdata.append(tdata)
        return outdata



    def valuesAvgStd(self):
        outdata = []
        for setdata in self.SetRat.values():
            for tdata in list(setdata.values()):
                outdata.append((tdata.Avg,tdata.Std))
        return outdata

    def keys(self):
        outdata = []
        for iset,setdata in self.SetRat.items():
            for values in list(setdata.keys()):
                outdata.append((iset,)+values)
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

    # def __add__(self,Rat2):
    #     if type(Rat2) in SetOfRat.parentlist:
    #         return Rat2.__radd__(self)
    #     else:
    #         result = SetOfRat(cfglist=self.cfglist, InfoDict=self.InfoDict)
    #         if isinstance(Rat2,SetOfRat):
    #             result.name = '+'.join([self.name,Rat2.name])
    #             for (iset,igamma,ip,it,iRat) in self.items():
    #                 if iset in Rat2.SetRat.keys():
    #                     result[(iset,igamma,ip,it)] = iRat + Rat2[(iset,igamma,ip,it)]
    #         elif isinstance(Rat2,TwoPointCorr):
    #             result.name = '+'.join([self.name,Rat2.name])
    #             for (iset,igamma,ip,it,iRat) in self.items():
    #                 result[(iset,igamma,ip,it)] = iRat + Rat2[(igamma,ip,it)]
    #         else:
    #             try:
    #                 result.name = '+'.join([self.name,NumbToFileName(Rat2)])
    #                 for iset,igamma,ip,it,iRat in self.items():
    #                     result[(iset,igamma,ip,it)] = iRat + Rat2
    #             except Exception as err:
    #                 raise EnvironmentError('Invalid value to combine with BootStrap class')
    #         return result

    # def __sub__(self,Rat2):
    #     if type(Rat2) in SetOfRat.parentlist:
    #         return Rat2.__rsub__(self)
    #     else:
    #         result = SetOfRat( cfglist=self.cfglist, InfoDict=self.InfoDict)
    #         if isinstance(Rat2,SetOfRat):
    #             result.name = '-'.join([self.name,Rat2.name])
    #             for (iset,igamma,ip,it,iRat) in self.items():
    #                 if iset in Rat2.SetRat.keys():
    #                     result[(iset,igamma,ip,it)] = iRat - Rat2[(iset,igamma,ip,it)]
    #         elif isinstance(Rat2,TwoPointCorr):
    #             for (iset,igamma,ip,it,iRat) in self.items():
    #                 result.name = '-'.join([self.name,Rat2.name])
    #                 result[(iset,igamma,ip,it)] = iRat - Rat2[(igamma,ip,it)]
    #         else:
    #             try:
    #                 for iset,igamma,ip,it,iRat in self.items():
    #                     result.name = '-'.join([self.name,NumbToFileName(Rat2)])
    #                     result[(iset,igamma,ip,it)] = iRat - Rat2
    #             except Exception as err:
    #                 raise EnvironmentError('Invalid value to combine with BootStrap class')
    #         return result


    # def __mul__(self,Rat2):
    #     if type(Rat2) in SetOfRat.parentlist:
    #         return Rat2.__rmul__(self)
    #     else:
    #         result = SetOfRat( cfglist=self.cfglist, InfoDict=self.InfoDict)
    #         if isinstance(Rat2,SetOfRat):
    #             result.name = 'x'.join([self.name,Rat2.name])
    #             for (iset,igamma,ip,it,iRat) in self.items():
    #                 if iset in Rat2.SetRat.keys():
    #                     result[(iset,igamma,ip,it)] = iRat * Rat2[(iset,igamma,ip,it)]
    #         elif isinstance(Rat2,TwoPointCorr):
    #             for (iset,igamma,ip,it,iRat) in self.items():
    #                 result.name = 'x'.join([self.name,Rat2.name])
    #                 result[(iset,igamma,ip,it)] = iRat * Rat2[(igamma,ip,it)]
    #         else:
    #             try:
    #                 for iset,igamma,ip,it,iRat in self.items():
    #                     result.name = 'x'.join([self.name,NumbToFileName(Rat2)])
    #                     result[(iset,igamma,ip,it)] = iRat * Rat2
    #             except Exception as err:
    #                 raise EnvironmentError('Invalid value to combine with BootStrap class')
    #         return result


    # def __div__(self,Rat2):
    #     if type(Rat2) in SetOfRat.parentlist:
    #         return Rat2.__rdiv__(self)
    #     else:
    #         result = SetOfRat( cfglist=self.cfglist, InfoDict=self.InfoDict)
    #         if isinstance(Rat2,SetOfRat):
    #             result.name = 'div'.join([self.name,Rat2.name])
    #             for (iset,igamma,ip,it,iRat) in self.items():
    #                 if iset in Rat2.SetRat.keys():
    #                     try:
    #                         result[(iset,igamma,ip,it)] = iRat / Rat2[(iset,igamma,ip,it)]
    #                     except Exception as err:
    #                         if any([ibs == 0 for ibs in Rat2[(iset,igamma,ip,it)].bootvals]):
    #                             raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+Rat2.name + ' '+Rat2[(iset,igamma,ip,it)].name)
    #         elif isinstance(Rat2,TwoPointCorr):
    #             for (iset,igamma,ip,it,iRat) in self.items():
    #                 result.name = 'div'.join([self.name,Rat2.name])
    #                 try:
    #                     result[(iset,igamma,ip,it)] = iRat / Rat2[(igamma,ip,it)]
    #                 except Exception as err:
    #                     if any([ibs == 0 for ibs in Rat2[(igamma,ip,it)].bootvals]):
    #                         raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+Rat2.name + ' '+Rat2[(igamma,ip,it)].name)
    #         else:
    #             if Rat2 == 0:
    #                 raise ZeroDivisionError('Dividing by zero found in bootstrap division')
    #             try:
    #                 for iset,igamma,ip,it,iRat in self.items():
    #                     result.name = 'div'.join([self.name,NumbToFileName(Rat2)])
    #                     result[(iset,igamma,ip,it)] = iRat / Rat2
    #             except Exception as err:
    #                 raise EnvironmentError('Invalid value to combine with BootStrap class')
    #         return result


    # def __pow__(self,Rat2):
    #     if type(Rat2) in SetOfRat.parentlist:
    #         return Rat2.__rpow__(self)
    #     else:
    #         result = SetOfRat( cfglist=self.cfglist, InfoDict=self.InfoDict)
    #         if isinstance(Rat2,SetOfRat):
    #             result.name = 'pow'.join([self.name,Rat2.name])
    #             for (iset,igamma,ip,it,iRat) in self.items():
    #                 if iset in Rat2.SetRat.keys():
    #                     result[(iset,igamma,ip,it)] = iRat ** Rat2[(iset,igamma,ip,it)]
    #         elif isinstance(Rat2,TwoPointCorr):
    #             for (iset,igamma,ip,it,iRat) in self.items():
    #                 result.name = 'pow'.join([self.name,Rat2.name])
    #                 result[(iset,igamma,ip,it)] = iRat ** Rat2[(igamma,ip,it)]
    #         else:
    #             try:
    #                 for iset,igamma,ip,it,iRat in self.items():
    #                     result.name = 'pow'.join([self.name,NumbToFileName(Rat2)])
    #                     result[(iset,igamma,ip,it)] = iRat ** Rat2
    #             except Exception as err:
    #                 raise EnvironmentError('Invalid value to combine with BootStrap class')
    #         return result


    # ## Right multiplication functions


    # def __radd__(self,Rat2):
    #     result = SetOfRat( cfglist=self.cfglist, InfoDict=self.InfoDict)
    #     if isinstance(Rat2,SetOfRat):
    #         result.name = '+'.join([Rat2.name,self.name])
    #         for (iset,igamma,ip,it,iRat) in self.items():
    #             if iset in Rat2.SetRat.keys():
    #                 result[(iset,igamma,ip,it)] = Rat2[(iset,igamma,ip,it)] + iRat
    #     elif isinstance(Rat2,TwoPointCorr):
    #         result.name = '+'.join([Rat2.name,self.name])
    #         for (iset,igamma,ip,it,iRat) in self.items():
    #             result[(iset,igamma,ip,it)] = Rat2[(igamma,ip,it)] + iRat
    #     else:
    #         try:
    #             result.name = '+'.join([NumbToFileName(Rat2),self.name])
    #             for iset,igamma,ip,it,iRat in self.items():
    #                 result[(iset,igamma,ip,it)] = Rat2 + iRat
    #         except Exception as err:
    #             raise EnvironmentError('Invalid value to combine with BootStrap class')
    #     return result


    # def __rsub__(self,Rat2):
    #     result = SetOfRat( cfglist=self.cfglist, InfoDict=self.InfoDict)
    #     if isinstance(Rat2,SetOfRat):
    #         result.name = '-'.join([Rat2.name,self.name])
    #         for (iset,igamma,ip,it,iRat) in self.items():
    #             if iset in Rat2.SetRat.keys():
    #                 result[(iset,igamma,ip,it)] = Rat2[(iset,igamma,ip,it)] - iRat
    #     elif isinstance(Rat2,TwoPointCorr):
    #         result.name = '-'.join([Rat2.name,self.name])
    #         for (iset,igamma,ip,it,iRat) in self.items():
    #             result[(iset,igamma,ip,it)] = Rat2[(igamma,ip,it)] - iRat
    #     else:
    #         try:
    #             result.name = '-'.join([NumbToFileName(Rat2),self.name])
    #             for iset,igamma,ip,it,iRat in self.items():
    #                 result[(iset,igamma,ip,it)] = Rat2 - iRat
    #         except Exception as err:
    #             raise EnvironmentError('Invalid value to combine with BootStrap class')
    #     return result

    # def __rmul__(self,Rat2):
    #     result = SetOfRat( cfglist=self.cfglist, InfoDict=self.InfoDict)
    #     if isinstance(Rat2,SetOfRat):
    #         result.name = 'x'.join([Rat2.name,self.name])
    #         for (iset,igamma,ip,it,iRat) in self.items():
    #             if iset in Rat2.SetRat.keys():
    #                 result[(iset,igamma,ip,it)] = Rat2[(iset,igamma,ip,it)] * iRat
    #     elif isinstance(Rat2,TwoPointCorr):
    #         result.name = 'x'.join([Rat2.name,self.name])
    #         for (iset,igamma,ip,it,iRat) in self.items():
    #             result[(iset,igamma,ip,it)] = Rat2[(igamma,ip,it)] * iRat
    #     else:
    #         try:
    #             result.name = 'x'.join([NumbToFileName(Rat2),self.name])
    #             for iset,igamma,ip,it,iRat in self.items():
    #                 result[(iset,igamma,ip,it)] = Rat2 * iRat
    #         except Exception as err:
    #             raise EnvironmentError('Invalid value to combine with BootStrap class')
    #     return result


    # def __rdiv__(self,Rat2):
    #     result = SetOfRat( cfglist=self.cfglist, InfoDict=self.InfoDict)
    #     if isinstance(Rat2,SetOfRat):
    #         result.name = 'div'.join([Rat2.name,self.name])
    #         for (iset,igamma,ip,it,iRat) in self.items():
    #             if iset in Rat2.SetRat.keys():
    #                 try:
    #                     result[(iset,igamma,ip,it)] = Rat2[(iset,igamma,ip,it)] / iRat
    #                 except Exception as err:
    #                     if any([ibs == 0 for ibs in iRat.bootvals]):
    #                         raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+self.name + ' '+iRat.name)
    #     elif isinstance(Rat2,TwoPointCorr):
    #         result.name = 'div'.join([Rat2.name,self.name])
    #         for (iset,igamma,ip,it,iRat) in self.items():
    #             try:
    #                 result[(iset,igamma,ip,it)] = Rat2[(igamma,ip,it)] / iRat
    #             except Exception as err:
    #                 if any([ibs == 0 for ibs in iRat.bootvals]):
    #                     raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+Rat2.name + ' '+Rat2[(igamma,ip,it)].name)
    #     else:
    #         result.name = 'div'.join([NumbToFileName(Rat2),self.name])
    #         try:
    #             for iset,igamma,ip,it,iRat in self.items():
    #                 result[(iset,igamma,ip,it)] = Rat2 / iRat
    #         except Exception as err:
    #             raise EnvironmentError('Invalid value to combine with BootStrap class')
    #     return result


    # def __rpow__(self,Rat2):
    #     result = SetOfRat( cfglist=self.cfglist, InfoDict=self.InfoDict)
    #     if isinstance(Rat2,SetOfRat):
    #         result.name = 'pow'.join([Rat2.name,self.name])
    #         for (iset,igamma,ip,it,iRat) in self.items():
    #             if iset in Rat2.SetRat.keys():
    #                 result[(iset,igamma,ip,it)] = Rat2[(iset,igamma,ip,it)] ** iRat
    #     elif isinstance(Rat2,TwoPointCorr):
    #         result.name = 'pow'.join([Rat2.name,self.name])
    #         for (iset,igamma,ip,it,iRat) in self.items():
    #             result[(iset,igamma,ip,it)] = Rat2[(igamma,ip,it)] ** iRat
    #     else:
    #         try:
    #             result.name = 'pow'.join([NumbToFileName(Rat2),self.name])
    #             for iset,igamma,ip,it,iRat in self.items():
    #                 result[(iset,igamma,ip,it)] = Rat2 ** iRat
    #         except Exception as err:
    #             raise EnvironmentError('Invalid value to combine with BootStrap class')
    #     return result
