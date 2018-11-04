#!/usr/bin/env python

import numpy as np

from MiscFuns import ODNested

from BootStrapping import BootStrap
import FitFunctions as ff
# from FitFunctions import Fitting
import pandas as pa
from Params import defPrec,defMaxIter,def_min_fitr,Debug_Mode,nboot
from copy import copy
from MiscFuns import get_val_float
from XmlFormatting import xmlfitr
from TimeStuff import Timer
from XmlFormatting import unxmlfitr_int
from copy import deepcopy

## keepkeys is for keeping track of "fit funtions, fit ranges etc..."
def PullSOFSeries_Wpar(this_series,keep_keys=slice(None),fmted=True):
    if len(this_series) == 0:
        return this_series
    outl,indexl = [],[]
    for ikey,SOF_data in this_series.items():
        for fit_key,fit_data in SOF_data.Get_Extrapolation(fmted=fmted).items():
            this_key = ikey
            if not isinstance(ikey,(list,tuple,np.ndarray)):
                this_key = (ikey,)
            indexl.append(tuple(this_key)+tuple(np.array(fit_key)[keep_keys]))
            outl.append(fit_data)
    if len(indexl) > 0:
        if this_series.index.names is None:
            print(this_series)
            raise EnvironmentError('No labels set up for series')
        if fmted:
            these_key_labs = list(this_series.index.names) + list(np.array(this_series.iloc[0].Fit_Col_Names_fmt)[keep_keys])+['parameters']
        else:
            these_key_labs = list(this_series.index.names) + list(np.array(this_series.iloc[0].Fit_Col_Names)[keep_keys])+['parameters']
        indicies = pa.MultiIndex.from_tuples(indexl,names=these_key_labs)
        return pa.Series(outl,index=indicies).dropna()
    else:
        return pa.Series()


def PullSOFSeries(this_series,keep_keys=slice(None),fmted=True):
    if len(this_series) == 0:
        return this_series
    outl,indexl = [],[]
    if fmted:
        for ikey,SOF_data in this_series.items():
            for fit_key,fit_data in SOF_data.Fit_Stats_fmt['Fit'].items():
                this_key = ikey
                if not isinstance(ikey,(list,tuple,np.ndarray)):
                    this_key = (ikey,)
                indexl.append(tuple(this_key)+tuple(np.array(fit_key)[keep_keys]))
                outl.append(fit_data)
    else:
        for ikey,SOF_data in this_series.items():
            for fit_key,fit_data in SOF_data.Fit_Stats['Fit'].items():
                this_key = ikey
                if not isinstance(ikey,(list,tuple,np.ndarray)):
                    this_key = (ikey,)
                indexl.append(tuple(this_key)+tuple(np.array(fit_key)[keep_keys]))
                outl.append(fit_data)
    if len(indexl) > 0:
        if this_series.index.names is None:
            print(this_series)
            raise EnvironmentError('No labels set up for series')
        if fmted:
            these_key_labs = list(this_series.index.names) + list(np.array(this_series.iloc[0].Fit_Col_Names_fmt)[keep_keys])
        else:
            these_key_labs = list(this_series.index.names) + list(np.array(this_series.iloc[0].Fit_Col_Names)[keep_keys])
        indicies = pa.MultiIndex.from_tuples(indexl,names=these_key_labs)
        return pa.Series(outl,index=indicies).dropna()
    else:
        return pa.Series()


class SetOfFitFuns(object):
    """

    FFun(set) "Inherits" from Fitting in FitFunctions.py



    """
    instances = []

    parentlist = []

    def __init__(   self, thisPrec=defPrec,thisMI=defMaxIter,data = 'None',name='None'):

        self.name = name
        self.thisPrec = thisPrec
        self.thisMI = thisMI
        self.ChiSorted = False
        self.ChiSorted_fmt = False
        self.VaryFunctions = False
        self.par_list = 'Invalid'
        self.Fit_Stats = pa.DataFrame()

        '''
        self.Fit_Col_Names DataFrame:

        columns =   Fitting, Chi2pdf,

        ## 1D case
        rows =      Function , iGuess , fit_min , fit_max

        ## 2D case "shape"
        rows =      Function , iGuess , fit_min , fit_max, shape

        ## 2D case "box"
        rows =      Function , iGuess , fit_min_1 , fit_max_1, fit_min_2, fit_max_2

        '''


        self.Fit_Stats_fmt = pa.DataFrame()
        '''
        self.Fit_Col_Names DataFrame:

        columns =   Fitting, Chi2pdf,

        ## 1D case
        rows =      Function , fitr#-#

        ## 2D case "shape"
        rows =      Function , fitr#-#, shape

        ## 2D case "box"
        rows =      Function , fitr#-#, fittwor#-#
        '''


        if not isinstance(data,str):
            self.ImportData(data)

    ## Pull out a specific row for multiindex please!!
    ## xdata  shape [ npar_fit , fit_index ]
    ## ydata  shape [ fit_index_1 , fit_index_2 , ...]
    ## currently only implemented for 2d fits
    def ImportData_Series(self,read_data):
        if isinstance(read_data.index,pa.MultiIndex):
            xdata = list(read_data.index)
            if read_data.index.nlevels == 1:
                self.xdata = np.array([get_val_float(ix[0]) for ix in xdata])
            if read_data.index.nlevels == 2:
                self.xdata = np.array(list(zip(*[list(map(get_val_float,ix)) for ix in xdata])))
                self.xdata_split = [sorted(list(set(self.xdata[0])))]
                self.xdata_split.append(sorted(list(set(self.xdata[1]))))
                self.xdata_split = np.array(self.xdata_split)
            else:
                raise EnvironmentError('Sets of fits not implemented for more than 2 fit parameters')
        else:
            self.xdata = np.array(list(map(get_val_float,list(read_data.index))))
        self.ydata = read_data.values
    '''
    Imports data for fitting:
    acceptable formats:

    read_data = pa.DataFrame()
    columns = xdata,ydata
    N.B. for some reason they have to be inside a 1-D array, i.e. pass in [xdata]

    '''
    def ImportData(self,read_data):
        if isinstance(read_data,pa.DataFrame):
            self.xdata = np.array(read_data['xdata'][0])
            self.ydata = np.array(read_data['ydata'][0])
        elif isinstance(read_data,pa.Series):
            self.ImportData_Series(read_data)
        elif read_data == 'None':
            ## also initialises :)
            self.RemoveData()
            self.xbooted = False
        elif isinstance(read_data,(tuple,list,np.ndarray)):
            self.ImportData(*read_data)
        else:
            raise IOError('unrecognised type for fitting'+str(type(read_data)))
        ## convention is always to have dimensions xdata [ iinput, ival]
        if self.xdata.ndim == 1:
            self.xdata = np.array([self.xdata])
        self.is_bootstrapped = isinstance(self.ydata.flatten()[0],BootStrap)
        self.ndim_fit = self.xdata.shape[0]


    def Get_Extrapolation(self,fmted=True):
        if fmted:
            if 'Fit' in self.Fit_Stats_fmt:
                return ff.PullFitSeries(self.Fit_Stats_fmt['Fit'])
            else:
                return pa.Series()
        else:
            if 'Fit' in self.Fit_Stats:
                return ff.PullFitSeries(self.Fit_Stats['Fit'])
            else:
                return pa.Series()





    ## acceptable formats for function are:
    ## ['Funs'] = (fun, npar, funder)

    ## ['Funs'] = (fun, npar)
    ## ['FunDer'] = funder
    ## ['FunDer'] not initialised

    ## ['Funs'] = (fun, funder)
    ## ['npar'] = npar

    ## ['Funs'] = [fun]
    ## ['npar'] = npar
    ## ['FunDer'] = funder
    ## ['FunDer'] not initialised

    ## ['Funs'] = fun
    ## ['npar'] = npar
    ## ['FunDer'] = funder
    ## ['FunDer'] not initialised

    def AddDefault(self,this_dict):
        if 'Funs' not in list(this_dict.keys()):
            this_dict['Funs'] = 'None'
        elif isinstance(this_dict['Funs'],(list,tuple,np.ndarray)) and len(this_dict['Funs']) > 1:
            if 'FunDer' in list(this_dict.keys()) and isinstance(this_dict['Funs'][1],int):
                this_dict['Funs'] = list(this_dict['Funs']) + [this_dict['FunDer']]
            elif isinstance(this_dict['Funs'][1],int):
                this_dict['Funs'] = list(this_dict['Funs'])
            else:
                if 'npar' in list(this_dict.keys()):
                    this_dict['Funs'] = [this_dict['Funs'][0],this_dict['npar'],this_dict['Funs'][1]]
                else:
                    raise IOError('npar is needed for fit functions')
        else:
            if len(this_dict['Funs']) == 1:
                this_dict['Funs'] = this_dict['Funs'][0]
            if 'npar' in list(this_dict.keys()):
                if 'FunDer' in list(this_dict.keys()):
                    this_dict['Funs'] = [this_dict['Funs'],this_dict['npar'],this_dict['FunDer']]
                else:
                    this_dict['Funs'] = [this_dict['Funs'],this_dict['npar']]
            else:
                raise IOError('npar is needed for fit funs')
        if 'FunDer' in this_dict:
            del this_dict['FunDer']
        if 'iGuess' not in list(this_dict.keys()):
            this_dict['iGuess'] = 'None'
        if 'name' not in list(this_dict.keys()):
            this_dict['name'] = self.name
        if 'paramlist' not in list(this_dict.keys()):
            this_dict['paramlist'] = 'PreDef'
        if 'paramunits' not in list(this_dict.keys()):
            this_dict['paramunits'] = 'None'
        if self.ndim_fit == 2:
            if 'fit_max_1' not in list(this_dict.keys()):
                # print 'Warning, default fit_max_1 not present, using end point'
                this_dict['fit_max_1'] = max(self.xdata[0])
            if 'fit_min_1' not in list(this_dict.keys()):
                # print 'Warning, default fit_min_1 not present, using start point'
                this_dict['fit_min_1'] = min(self.xdata[0])
            if 'fit_max_2' not in list(this_dict.keys()):
                # print 'Warning, default fit_max_2 not present, using end point'
                this_dict['fit_max_2'] = max(self.xdata[1])
            if 'fit_min_2' not in list(this_dict.keys()):
                # print 'Warning, default fit_min_2 not present, using start point'
                this_dict['fit_min_2'] =  min(self.xdata[1])
        else:
            if 'fit_max' not in list(this_dict.keys()):
                # print 'Warning, default fit_max not present, using end point'
                this_dict['fit_max'] = max(self.xdata[0])
            if 'fit_min' not in list(this_dict.keys()):
                # print 'Warning, default fit_min not present, using start point'
                this_dict['fit_min'] = min(self.xdata[0])
        return this_dict

    ## Also Appends
    def ImportFits(self,fit_types):
        if not hasattr(self,'xdata'):
            raise EnvironmentError('pease import data before import fit information')
        fit_list,ilist = [],[]
        # print 'DEBUG'
        # print self.xdata
        # print self.ydata
        # print pa.Series(fit_types)
        # print
        for ifit in fit_types:
            # print 'DEBUG'
            # print ifit['fit_min_1']
            # print ifit['fit_max_1']
            # print ifit['fit_min_2']
            # print ifit['fit_max_2']
            ifit = self.AddDefault(ifit)
            if self.ndim_fit == 1:
                this_key = (ifit['Funs'][0].__name__,str(ifit['iGuess']),
                            str(ifit['fit_min']),str(ifit['fit_max']))
                this_xdata = self.xdata[:,ifit['fit_min']:ifit['fit_max']+1]
                this_ydata = self.ydata[ifit['fit_min']:ifit['fit_max']+1]
            elif self.ndim_fit == 2:
                if 'fit_min_1' not in list(ifit.keys()):
                    this_key = (ifit['Funs'][0].__name__,str(ifit['iGuess']),
                                str(ifit['fit_min']),str(ifit['fit_max']))
                    this_xdata = self.xdata[:,ifit['fit_min']:ifit['fit_max']+1]
                    this_ydata = self.ydata[ifit['fit_min']:ifit['fit_max']+1,ifit['fit_min']:ifit['fit_max']+1]
                elif isinstance(ifit['fit_min_2'],(list,tuple,np.ndarray)):
                    ## TODO This is not working, to fix!!
                    this_key = (ifit['Funs'][0].__name__,str(ifit['iGuess']),
                                str(ifit['fit_min_1']),str(ifit['fit_max_1']),'shape1')
                    counter = 2
                    while this_key in ilist:
                        this_key = (ifit['Funs'][0].__name__,str(ifit['iGuess']),
                                    str(ifit['fit_min_1']),str(ifit['fit_max_1']),'shape'+str(counter))
                        counter += 1
                    this_xdata,this_ydata = [],[]
                    counter = 0
                    for ic1,ix1 in enumerate(self.xdata[0]):
                        if ifit['fit_min_1'] <= ic1 <= ifit['fit_max_1']+1:
                            for ic2,ix2 in enumerate(self.xdata[1]):
                                if ifit['fit_min_2'][counter] <= ic2 <= ifit['fit_max_2'][counter]+1:
                                    this_xdata.append([ix1,ix2])
                                    this_ydata.append(self.ydata[ic1,ic2])
                                counter += 1
                    if len(this_xdata) > 0:
                        this_xdata = np.swapaxes(np.array(this_xdata),0,1)
                        this_ydata = np.array(this_ydata)
                    else:
                        raise EnvironmentError('Fit data not found in fitting...')
                else:
                    this_key = (ifit['Funs'][0].__name__,str(ifit['iGuess']),
                                str(ifit['fit_min_1']),str(ifit['fit_max_1']),
                                str(ifit['fit_min_2']),str(ifit['fit_max_2']))
                    this_xdata,this_ydata = [],[]
                    # print('DEBUG')
                    # print('fit_1 ',self.xdata_split[0][ifit['fit_min_1']],' ',self.xdata_split[0][ifit['fit_max_1']])
                    # print('fit_2 ',self.xdata_split[1][ifit['fit_min_2']],' ',self.xdata_split[1][ifit['fit_max_2']])
                    for ic,(ix1,ix2) in enumerate(zip(*self.xdata)):
                        testbool = self.xdata_split[0][ifit['fit_min_1']] <= ix1 <= self.xdata_split[0][ifit['fit_max_1']]
                        testbool = testbool and self.xdata_split[1][ifit['fit_min_2']] <= ix2 <= self.xdata_split[1][ifit['fit_max_2']]
                        if testbool:
                            this_xdata.append([ix1,ix2])
                            this_ydata.append(self.ydata[ic])
                    # this_xdata = np.array(this_xdata)
                    # if Debug_Mode:
                    #     print 'DEBUG'
                    #     print this_xdata
                    if len(this_xdata) > 0:
                        this_xdata = np.swapaxes(np.array(this_xdata),0,1)
                        this_ydata = np.array(this_ydata)
                    else:
                        print('DEBUG')
                        print('fit_1_index ',ifit['fit_min_1'],' ',ifit['fit_max_1'])
                        print('fit_2_index ',ifit['fit_min_2'],' ',ifit['fit_max_2'])
                        print('fit_1 ',self.xdata_split[0][ifit['fit_min_1']],' ',self.xdata_split[0][ifit['fit_max_1']])
                        print('fit_2 ',self.xdata_split[1][ifit['fit_min_2']],' ',self.xdata_split[1][ifit['fit_max_2']])
                        print('xdata_split_0')
                        print('\n'.join(map(str,self.xdata_split[0])))
                        print('xdata_split_1')
                        print('\n'.join(map(str,self.xdata_split[1])))
                        for ic,(ix1,ix2) in enumerate(zip(*self.xdata)):
                            testbool = self.xdata_split[0][ifit['fit_min_1']] <= ix1 <= self.xdata_split[0][ifit['fit_max_1']]
                            testbool = testbool and self.xdata_split[1][ifit['fit_min_2']] <= ix2 <= self.xdata_split[1][ifit['fit_max_2']]
                            print(ix1,ix2,testbool)
                        raise EnvironmentError('Fit data not found in fitting...')
            # print('DEBUG ndim'+str(self.ndim_fit))
            # for ikey,ival in ifit.items():
            #     print('  ',ikey,ival)
            # print('Values:')
            # for ix_0,ix_1,iy in zip(this_xdata[0],this_xdata[1],this_ydata):
            #     print('  ',ix_0,ix_1,iy.Avg,iy.Std)
            # print()

            this_fit = ff.Fitting( thisPrec=self.thisPrec,thisMI=self.thisMI,
                                Funs=ifit['Funs'],data = [this_xdata,this_ydata],
                                iGuess=ifit['iGuess'],name=ifit['name'],paramlist=ifit['paramlist'],
                                paramunits=ifit['paramunits'])
            fit_list.append(this_fit)
            ilist.append(this_key)
        if len(ilist) > 0:
            if len(ilist[0]) == 4:
                this_Fit_Col_Names = ['Function','iGuess','fit_min','fit_max']
            if len(ilist[0]) == 5:
                this_Fit_Col_Names = ['Function','iGuess','fit_min','fit_max','shape']
            if len(ilist[0]) == 6:
                this_Fit_Col_Names = ['Function','iGuess','fit_min_1','fit_max_1','fit_min_2','fit_max_2']
            if hasattr(self,'Fit_Col_Names') and len(this_Fit_Col_Names) != len(self.Fit_Col_Names):
                print('WARNING:, different shapes of fit has been calculated, ignoring result')
            else:
                self.Fit_Col_Names = this_Fit_Col_Names
                indicies = pa.MultiIndex.from_tuples(ilist,names=self.Fit_Col_Names)
                if 'Fit' in self.Fit_Stats.columns:
                    this_df = pa.DataFrame()
                    this_df.loc[:,'Fit'] = pa.Series(fit_list,index=indicies)
                    self.Fit_Stats = self.Fit_Stats.append(this_df)
                else:
                    self.Fit_Stats.loc[:,'Fit'] = pa.Series(fit_list,index=indicies)
        else:
            out_str = 'ERROR, no fit ranges were found for fit types \n'
            out_str += '\n'.join(map(str,fit_types)) + '\n'
            raise EnvironmentError(out_str)



    def DoFits(self,PullParams=True,PullChi=True,show_timer=True):
        if 'Fit' in self.Fit_Stats:
            if show_timer: this_timer = Timer(linklist=list(range(len(self.Fit_Stats.index))),name='Fitting '+self.name)
            for ikey,ifit in self.Fit_Stats.loc[:,'Fit'].items():
                successfull = ifit.FitBoots()
                if successfull:
                    if show_timer: this_timer.Lap(ikey)
                else:
                    if show_timer: this_timer.Lap(ikey + ('Alread Computed',))
            if PullParams: self.PullIdenticalParams()
            if PullChi: self.PullChi()
            self.MakeFormatted()
        else:
            print('Fit not in Fit_Stats')
        return self
        # else:
            # raise EnvironmentError('Please import data and fits before running DoFits')

    def UnExpParams(self,paramkeys='All'):
        if 'Fit' in self.Fit_Stats:
            self.Fit_Stats['Fit'].apply(lambda x : x.UnExpParams(paramkeys=paramkeys))
        if 'Fit' in self.Fit_Stats_fmt:
            self.Fit_Stats_fmt['Fit'].apply(lambda x : x.UnExpParams(paramkeys=paramkeys))

    def ImportFunction(self,Function,thisparlen,FunDer='None',thisparlab = 'PreDef',thisparunits='None'):
        vall,indexl,vall_fmt,indexl_fmt = [],[],[],[]
        if 'Fit' in self.Fit_Stats:
            self.Fit_Stats['Fit'].apply(lambda x : x.ImportFunction(Function,thisparlen,FunDer=FunDer,thisparlab = thisparlab,thisparunits=thisparunits))
            for ikey,ival in self.Fit_Stats['Fit'].items():
                new_key = list(ikey)
                new_key[0] = Function.__name__
                vall.append(ival)
                indexl.append(tuple(new_key))
            if len(indexl) > 0:
                indicies = pa.MultiIndex.from_tuples(indexl,names=self.Fit_Stats.index.names)
                self.Fit_Stats = pa.Series(vall,index=indicies).to_frame(name='Fit')
        if 'Fit' in self.Fit_Stats_fmt:
            self.Fit_Stats_fmt['Fit'].apply(lambda x : x.ImportFunction(Function,thisparlen,FunDer=FunDer,thisparlab = thisparlab,thisparunits=thisparunits))
            for ikey_fmt,ival_fmt in self.Fit_Stats_fmt['Fit'].items():
                new_key_fmt = list(ikey_fmt)
                new_key_fmt[0] = Function.__name__
                vall_fmt.append(ival_fmt)
                indexl_fmt.append(tuple(new_key_fmt))
            if len(indexl_fmt) > 0:
                indicies_fmt = pa.MultiIndex.from_tuples(indexl_fmt,names=self.Fit_Stats_fmt.index.names)
                self.Fit_Stats_fmt = pa.Series(vall_fmt,index=indicies_fmt).to_frame(name='Fit')


    def PullChi(self,failed_fit=False):
        def GetChi(x):
            if 'Chi2DoF' in x.fit_data:
                return x.fit_data.loc[:,'Chi2DoF'].iloc[0]
            else:
                return BootStrap(bootvals=[float('NaN')]*nboot)
        if 'Fit' in self.Fit_Stats:
            self.Fit_Stats.loc[:,'Chi2pdf'] = self.Fit_Stats.loc[:,'Fit'].apply(GetChi)
            self.Fit_Stats.loc[:,'Chi2pdfAvg'] =self.Fit_Stats.loc[:,'Chi2pdf'].apply(lambda x : x.Avg)
            self.Fit_Stats.loc[:,'Chi2pdfStd'] =self.Fit_Stats.loc[:,'Chi2pdf'].apply(lambda x : x.Std)
            if len(self.Fit_Stats_fmt) > 0:
                self.Fit_Stats_fmt.loc[:,'Chi2pdf'] = self.Fit_Stats_fmt.loc[:,'Fit'].apply(GetChi)
                self.Fit_Stats_fmt.loc[:,'Chi2pdfAvg'] =self.Fit_Stats_fmt.loc[:,'Chi2pdf'].apply(lambda x : x.Avg)
                self.Fit_Stats_fmt.loc[:,'Chi2pdfStd'] =self.Fit_Stats_fmt.loc[:,'Chi2pdf'].apply(lambda x : x.Std)
        elif failed_fit:
            raise EnvironmentError('Error fitting, something went wrong when fitting')
        else:
            print('Warning, called PullChi before DoFits, running fits now')
            self.DoFits()
            self.PullChi(failed_fit=True)

    def MakeFitRanges_shape(self,fit_min_1,fit_max_1,fit_min_2,fit_max_2,npar):
        out_list = []
        out_2_list = []
        for ifit_min_1 in range(fit_min_1,fit_max_1+1):
            for ifit_max_1 in range(ifit_min_1+npar,fit_max_1+1):
                out_list.append((ifit_min_1,ifit_max_1))
                for ifit in range(ifit_min_1,ifit_max_1):
                    out_2_list.append(self.MakeFitRanges(fit_min_2[ifit],fit_min_2[ifit],npar))
        return out_list,out_2_list

    def GetXdataValue(self,fit_min,fit_max,index):
        sort_xvals = np.array(self.xdata_split[index])
        if get_val_float(fit_min) not in sort_xvals:
            print('fit_min ',fit_min,' not in xdata, finding closest values above')
            fit_min = sort_xvals[sort_xvals > fit_min]
            if len(fit_min) == 0:
                out_str = 'fit_min is larger than all data for index:'+str(index)+'\n'
                out_str += 'fit_min '+fit_min +'\n'
                out_str += 'xdata: \n'
                out_str += '\n  '.join(map(str,sort_xvals))
                raise EnvironmentError(out_str)
            else:
                fit_min = fit_min[0]
        if get_val_float(fit_max) not in sort_xvals:
            print('fit_max ',fit_max,' not in xdata, finding closest values below')
            fit_max = sort_xvals[sort_xvals < fit_max]
            if len(fit_max) == 0:
                out_str = 'fit_max is smaller than all data for index:'+str(index)+'\n'
                out_str += 'fit_max '+fit_max +'\n'
                out_str += 'xdata: \n'
                out_str += '\n  '.join(map(str,sort_xvals))
                raise EnvironmentError(out_str)
            else:
                fit_max = fit_max[-1]
        return fit_min,fit_max


    def MakeFitRanges(self,fit_min,fit_max,npar,from_xdata=False,index=0):
        if from_xdata:
            fit_min,fit_max = self.GetXdataValue(fit_min,fit_max,index)
            fit_min = self.xdata_split[index].index(get_val_float(fit_min))
            fit_max = self.xdata_split[index].index(get_val_float(fit_max))
        out_list = []
        for ifit_min in range(fit_min,fit_max+1):
            for ifit_max in range(ifit_min+npar-1,fit_max+1):
                out_list.append((ifit_min,ifit_max))
        return out_list

    def MakeFitRanges_FixEnd(self,fit_min,fit_max,npar,from_xdata=False,index=0):
        if from_xdata:
            fit_min,fit_max = self.GetXdataValue(fit_min,fit_max,index)
            fit_min_i = self.xdata_split[index].index(get_val_float(fit_min))
            fit_max_i = self.xdata_split[index].index(get_val_float(fit_max))
            # print('DEBUG fit_min_i',fit_min_i)
            # print('DEBUG fit_max_i',fit_max_i)
        out_list = []
        for ifit_min in range(fit_min_i,fit_max_i-npar+1):
            out_list.append((ifit_min,fit_max_i))
        return out_list,fit_min,fit_max

    def MakeFitRanges_Sym(self,fit_min,fit_max,npar):
        out_list = []
        if (fit_max-fit_min-npar)>0:
            for icut in range((fit_max-fit_min-npar)//2):
                out_list.append((fit_min+icut,fit_max-icut))
        return out_list

    def CheckRange_fitr(self,this_fit,fit_info={},min_fit_len=def_min_fitr):
        if isinstance(this_fit,str):
            fit_min,fit_max = unxmlfitr_int(this_fit)
        elif isinstance(this_fit,(tuple,list,np.ndarray)):
            fit_min,fit_max = list(map(int,this_fit))
        else:
            print()
            raise IOError('incorrect format for fitr#-# : '+str(this_fit) )
        return self.CheckRange(fit_min,fit_max,fit_info=fit_info,min_fit_len=min_fit_len)

    def CheckRange(self,fit_min,fit_max,fit_info={},min_fit_len=def_min_fitr):
        if self.ndim_fit != 1:
            print('CheckRange must be used for 1d fitting')
            return False
        this_npar = self.AddDefault(fit_info)['Funs'][1]
        fit_ranges = self.MakeFitRanges(fit_min,fit_max,this_npar+min_fit_len)
        if len(fit_ranges) == 0:
            print('Warning, no acceptable fit ranges to check for:',fit_min,fit_max,this_npar+min_fit_len)
            return False
        if (fit_info['Funs'][0].__name__,str(fit_info['iGuess'])) not in self.Fit_Stats['Fit']:
            return False
        for ifit_min,ifit_max in fit_ranges:
            this_tuple = (fit_info['Funs'][0].__name__,str(fit_info['iGuess']),str(ifit_min),str(ifit_max))
            if this_tuple not in self.Fit_Stats['Fit'].index:
                print(this_tuple, 'not fitted')
                return False
        return True


    def ScanRange_fitr(self,this_fit,fit_info={},min_fit_len=def_min_fitr):
        if isinstance(this_fit,str):
            fit_min,fit_max = unxmlfitr_int(this_fit)
        elif isinstance(this_fit,(tuple,list,np.ndarray)):
            fit_min,fit_max = list(map(int,this_fit))
        else:
            print()
            raise IOError('incorrect format for fitr#-# : '+str(this_fit) )
        return self.ScanRange(fit_min,fit_max,fit_info=fit_info,min_fit_len=min_fit_len)

    def ScanRange_fitr_Sym(self,this_fit,fit_info={},min_fit_len=def_min_fitr):
        if isinstance(this_fit,str):
            fit_min,fit_max = unxmlfitr_int(this_fit)
        elif isinstance(this_fit,(tuple,list,np.ndarray)):
            fit_min,fit_max = list(map(int,this_fit))
        else:
            print()
            raise IOError('incorrect format for fitr#-# : '+str(this_fit) )
        return self.ScanRange_Sym(fit_min,fit_max,fit_info=fit_info,min_fit_len=min_fit_len)


    ## min_fit_len is added into npar, i.e. only accepts fit ranges past npar+min_fit_len
    ## zero is npar.
    def ScanRange(self,fit_min,fit_max,fit_info={},min_fit_len=def_min_fitr,from_xdata=False):
        if self.ndim_fit != 1:
            print('ScanRange must be used for 1d fitting')
            return False
        this_npar = self.AddDefault(fit_info)['Funs'][1]
        fit_ranges = self.MakeFitRanges(fit_min,fit_max,this_npar+min_fit_len,from_xdata=from_xdata)
        if len(fit_ranges) == 0:
            print('Warning, no acceptable fit ranges for:',fit_min,fit_max,this_npar+min_fit_len)
            return False
        parse_info = []
        for ifit_min,ifit_max in fit_ranges:
            this_info = copy(fit_info)
            this_info['fit_min'] = ifit_min
            this_info['fit_max'] = ifit_max
            parse_info.append(this_info)
        self.ImportFits(parse_info)
        return True

            ## min_fit_len is added into npar, i.e. only accepts fit ranges past npar+min_fit_len
    ## zero is npar.
    def ScanRange_Sym(self,fit_min,fit_max,fit_info={},min_fit_len=def_min_fitr):
        if self.ndim_fit != 1:
            print('ScanRange must be used for 1d fitting')
            return False
        this_npar = self.AddDefault(fit_info)['Funs'][1]
        fit_ranges = self.MakeFitRanges_Sym(fit_min,fit_max,this_npar+min_fit_len)
        if len(fit_ranges) == 0:
            print('Warning, no acceptable fit ranges for:',fit_min,fit_max,this_npar+min_fit_len)
            return False
        parse_info = []
        for ifit_min,ifit_max in fit_ranges:
            this_info = copy(fit_info)
            this_info['fit_min'] = ifit_min
            this_info['fit_max'] = ifit_max
            parse_info.append(this_info)
        self.ImportFits(parse_info)
        return True


    def ScanShape(self,fit_min_1,fit_max_1,fit_min_2,fit_max_2,fit_info={},min_fit_len=def_min_fitr):
        if self.ndim_fit != 2:
            print('ScanShape must be used for 2d fitting')
            return False
        this_npar = self.AddDefault(fit_info)['Funs'][1]
        fit_ranges,fit_ranges_2 = self.MakeFitRanges_shape(fit_min_1,fit_max_1,fit_min_2,fit_max_2,this_npar+min_fit_len)
        if len(fit_ranges) == 0:
            print('Warning, no acceptable fit ranges for shape:')
            print(fit_min_1,fit_max_1,fit_min_2,fit_max_2,this_npar+min_fit_len)
            return False
        parse_info = []
        for ifit_min_1,ifit_max_1 in fit_ranges:
            for ifit_min_2,ifit_max_2 in fit_ranges_2:
                this_info = copy(fit_info)
                this_info['fit_min_1'] = ifit_min_1
                this_info['fit_max_1'] = ifit_max_1
                this_info['fit_min_2'] = ifit_min_2
                this_info['fit_max_2'] = ifit_max_2
                parse_info.append(this_info)
        self.ImportFits(parse_info)
        return True

    def ScanBox_fitr(self,fit_1,fit_2,fit_info={},min_fit_len_1=def_min_fitr,min_fit_len_2=def_min_fitr):
        if isinstance(fit_1,str):
            fit_min_1,fit_max_1 = unxmlfitr_int(fit_1)
        elif isinstance(fit_1,(tuple,list,np.ndarray)):
            fit_min_1,fit_max_1 = list(map(int,fit_1))
        else:
            print()
            raise IOError('incorrect format for fitr#-# : '+str(fit_1) )
        if isinstance(fit_2,str):
            fit_min_2,fit_max_2 = unxmlfitr_int(fit_2)
        elif isinstance(fit_2,(tuple,list,np.ndarray)):
            fit_min_2,fit_max_2 = list(map(int,fit_2))
        else:
            print()
            raise IOError('incorrect format for fittwor#-# : '+str(fit_2) )
        return self.ScanBox(fit_min_1,fit_max_1,fit_min_2,fit_max_2,fit_info=fit_info,
                            min_fit_len_1=min_fit_len_1,min_fit_len_2=min_fit_len_2)

    def ScanBox(self,fit_min_1,fit_max_1,fit_min_2,fit_max_2,fit_info={},
                min_fit_len_1=def_min_fitr,min_fit_len_2=def_min_fitr,from_xdata=False):
        if self.ndim_fit != 2:
            print('ScanBox must be used for 2d fitting')
            return False
        this_npar = self.AddDefault(fit_info)['Funs'][1]
        fit_ranges = self.MakeFitRanges(fit_min_1,fit_max_1,this_npar+min_fit_len_1,from_xdata=from_xdata,index=0)
        fit_ranges_2 = self.MakeFitRanges(fit_min_2,fit_max_2,this_npar+min_fit_len_2,from_xdata=from_xdata,index=1)
        if len(fit_ranges) == 0 or len(fit_ranges_2) == 0:
            print('Warning, no acceptable fit ranges for:')
            print(fit_min_1,fit_max_1,fit_min_2,fit_max_2, end=' ')
            print('min fit ranges')
            print(this_npar+min_fit_len_1,this_npar+min_fit_len_2)
            return False
        parse_info = []
        for ifit_min_1,ifit_max_1 in fit_ranges:
            for ifit_min_2,ifit_max_2 in fit_ranges_2:
                # print('Adding fitr'+str(ifit_min_1)+'-'+str(ifit_max_1)+' fittwor'+str(ifit_min_2)+'-'+str(ifit_max_2))
                this_info = copy(fit_info)
                this_info['fit_min_1'] = ifit_min_1
                this_info['fit_max_1'] = ifit_max_1
                this_info['fit_min_2'] = ifit_min_2
                this_info['fit_max_2'] = ifit_max_2
                parse_info.append(this_info)
        self.ImportFits(parse_info)
        return True


    def ScanBox_FixEnd(self,fit_min_1,fit_max_1,fit_min_2,fit_max_2,fit_info={},
                min_fit_len_1=def_min_fitr,min_fit_len_2=def_min_fitr,from_xdata=False):
        if self.ndim_fit != 2:
            print('ScanBox must be used for 2d fitting')
            return False
        this_npar = self.AddDefault(fit_info)['Funs'][1]
        fit_ranges,fit_min_1,fit_max_1 = self.MakeFitRanges_FixEnd(fit_min_1,fit_max_1,this_npar+min_fit_len_1,from_xdata=from_xdata,index=0)
        fit_ranges_2,fit_min_2,fit_max_2 = self.MakeFitRanges_FixEnd(fit_min_2,fit_max_2,this_npar+min_fit_len_2,from_xdata=from_xdata,index=1)
        if len(fit_ranges) == 0 or len(fit_ranges_2) == 0:
            print('Warning, no acceptable fit ranges for:')
            print(fit_min_1,fit_max_1,fit_min_2,fit_max_2, end=' ')
            print('min fit ranges')
            print(this_npar+min_fit_len_1,this_npar+min_fit_len_2)
            print('xdata dim1')
            print('\n'.join(map(str,self.xdata_split[0])))
            print('xdata dim2')
            print('\n'.join(map(str,self.xdata_split[1])))
            return False
        parse_info = []
        for ifit_min_1,ifit_max_1 in fit_ranges:
            for ifit_min_2,ifit_max_2 in fit_ranges_2:
                # print('Adding fitr'+str(ifit_min_1)+'-'+str(ifit_max_1)+' fittwor'+str(ifit_min_2)+'-'+str(ifit_max_2))
                this_info = copy(fit_info)
                this_info['fit_min_1'] = ifit_min_1
                this_info['fit_max_1'] = ifit_max_1
                this_info['fit_min_2'] = ifit_min_2
                this_info['fit_max_2'] = ifit_max_2
                parse_info.append(this_info)
        self.ImportFits(parse_info)
        return True

    def SingleFitBox(self,fit_min_1,fit_max_1,fit_min_2,fit_max_2,fit_info={},min_fit_len=def_min_fitr):
        if self.ndim_fit != 2:
            print('ScanBox must be used for 2d fitting')
            return False
        fit_info['fit_min_1'] = fit_min_1
        fit_info['fit_max_1'] = fit_max_1
        fit_info['fit_min_2'] = fit_min_2
        fit_info['fit_max_2'] = fit_max_2
        self.ImportFits(fit_info)
        return True




    def SortChi(self):
        if not self.ChiSorted:
            if 'Chi2pdf' not in self.Fit_Stats:
                self.PullChi()
            self.Fit_Stats = self.Fit_Stats.sort_values('Chi2pdfAvg')
            if len(self.Fit_Stats_fmt) > 0:
                if 'Chi2pdf' not in self.Fit_Stats_fmt:
                    self.PullChi()
                self.Fit_Stats_fmt = self.Fit_Stats_fmt.sort_values('Chi2pdfAvg')
                self.ChiSorted_fmt = True
            self.ChiSorted = True

    def PlotFunction(self,ikey,otherXvals=[],xaxis=0,xdatarange='Data',thislineres=100,
                    color=None,shift=0.0,flowphys=False,ShowPar='None',scale=1.0,
                    plot_plane=None,hairline=False):
        if isinstance(ikey,int):
            fit_data = self.Fit_Stats[:,'Fit'].iloc[ikey]
        else:
            if ikey not in self.Fit_Stats.index:
                if ikey not in self.Fit_Stats_fmt:
                    print('WARNING: ',ikey,'not found in Fit_Stats(_fmt), skipping plot')
                else:
                    fit_data = self.Fit_Stats_fmt.loc[ikey,'Fit']
            else:
                fit_data = self.Fit_Stats.loc[ikey,'Fit']
        plot_plane = fit_data.PlotFunction(otherXvals=otherXvals,xaxis=xaxis,
                    xdatarange=xdatarange,thislineres=thislineres,
                    color=color,shift=shift,flowphys=flowphys,ShowPar=ShowPar,scale=scale,
                    plot_plane=plot_plane,hairline=hairline)
        return plot_plane


    '''
    Generic function for varying one thing, and keeping the rest FixRead
    this_list is what you are varying, flag is its index
    and fit_info is everything else


    '''
    def ImportVarySomething(self,this_list,fit_info,flag):
        parse_info = []
        for ival in this_list:
            this_info = fit_info
            this_info[flag] = ival
            parse_info.append(this_info)
        if 'Fun' in flag:
            self.VaryFunctions = True
        self.ImportFits(parse_info)

    # def PullAllParams(self):
    #     if 'Fit' in self.Fit_Stats:
    #         vall,avgl,stdl = [],[],[]
    #         for ifit,fit_data in self.Fit_Stats.loc[:,'Fit'].iteritems():
    #             vall.append([])
    #             avgl.append([])
    #             stdl.append([])
    #             for ipar in list(fit_data.fit_data.loc[:,'Params'].index):
    #                 if isinstance(ipar,int):
    #                     ipar = fit_data.fit_data.loc[:,'Params'].index[ipar]
    #                 # self.Fit_Stats.loc[ifit,'par_'+ipar] = fit_data.fit_data.apply(debug_fun)
    #                 this_val = fit_data.fit_data['Params'][ipar]
    #                 vall[-1].append(this_val)
    #
    #                 self.Fit_Stats['par_'+ipar][ifit] =
    #                 self.Fit_Stats['par_'+ipar+'_Avg'][ifit] = self.Fit_Stats['par_'+ipar][ifit].Avg
    #                 self.Fit_Stats['par_'+ipar+'_Std'][ifit] = self.Fit_Stats['par_'+ipar][ifit].Std
    #     else:
    #         print 'Warning, called PullParam before DoFits, running fits now'
    #         self.DoFits()
    #         self.PullAllParam()

    def GetFuns(self):
        if self.Fit_Stats is not None and 'Fit' in self.Fit_Stats.columns:
            self.Fit_Stats.loc[:,'Fit'].apply(lambda x : x.GetFuns())
        if self.Fit_Stats_fmt is not None and 'Fit' in self.Fit_Stats_fmt.columns:
            self.Fit_Stats_fmt.loc[:,'Fit'].apply(lambda x : x.GetFuns())

    def RemoveFuns(self):
        if self.Fit_Stats is not None and 'Fit' in self.Fit_Stats.columns:
            self.Fit_Stats.loc[:,'Fit'].apply(lambda x : x.RemoveFuns())
        if self.Fit_Stats_fmt is not None and 'Fit' in self.Fit_Stats_fmt.columns:
            self.Fit_Stats_fmt.loc[:,'Fit'].apply(lambda x : x.RemoveFuns())


    def GetBestChi(self):
        self.SortChi()
        return self.Fit_Stats.loc[:,'Fit'].iloc[0],self.Fit_Stats.loc[:,'Chi2pdfAvg'].iloc[0]


    def GetBestChi_key(self):
        self.SortChi()
        return self.Fit_Stats.index[0]

    def ReduceViaChi(self,chi_min=0,chi_max=100,nchi_cutoff = 100):
        self.SortChi()
        outputChi = deepcopy(self)
        outputChi.Fit_Stats = self.Fit_Stats[self.Fit_Stats['Chi2pdfAvg'] > chi_min]
        outputChi.Fit_Stats = outputChi.Fit_Stats[outputChi.Fit_Stats['Chi2pdfAvg'] < chi_max]
        outputChi.Fit_Stats_fmt = self.Fit_Stats_fmt[self.Fit_Stats_fmt['Chi2pdfAvg'] > chi_min]
        outputChi.Fit_Stats_fmt = outputChi.Fit_Stats_fmt[outputChi.Fit_Stats_fmt['Chi2pdfAvg'] < chi_max]
        if len(outputChi.Fit_Stats.index) > nchi_cutoff:
            rand_list = np.random.choice(len(outputChi.Fit_Stats.index),size=nchi_cutoff,replace=False)
            outputChi.Fit_Stats = outputChi.Fit_Stats.iloc[rand_list,:]
        elif len(outputChi.Fit_Stats.index) == 0:
            rand_list = np.array([min(self.Fit_Stats['Chi2pdfAvg'])])
            outputChi.Fit_Stats = self.Fit_Stats[self.Fit_Stats['Chi2pdfAvg'] == min(self.Fit_Stats['Chi2pdfAvg'])]
        else:
            rand_list = 'All'
        if len(outputChi.Fit_Stats_fmt.index) > nchi_cutoff:
            rand_list = np.random.choice(len(outputChi.Fit_Stats_fmt.index),size=nchi_cutoff,replace=False)
            outputChi.Fit_Stats_fmt = outputChi.Fit_Stats_fmt.iloc[rand_list,:]
        elif len(outputChi.Fit_Stats_fmt.index) == 0:
            rand_list = np.array([min(self.Fit_Stats_fmt['Chi2pdfAvg'])])
            outputChi.Fit_Stats_fmt = self.Fit_Stats_fmt[self.Fit_Stats_fmt['Chi2pdfAvg'] == min(self.Fit_Stats_fmt['Chi2pdfAvg'])]
        else:
            rand_list = 'All'
        return outputChi,rand_list


    def GetChiRange(self,chi_min=0,chi_max=100,nchi_cutoff = 100,fmted=True):
        self.SortChi()
        if fmted:
            this_fit = self.Fit_Stats_fmt
        else:
            this_fit = self.Fit_Stats
        temp_df = this_fit[this_fit['Chi2pdfAvg'] > chi_min]
        temp_df = temp_df[this_fit['Chi2pdfAvg'] < chi_max]
        if len(temp_df.index) > nchi_cutoff:
            rand_list = np.random.choice(len(temp_df.index),size=nchi_cutoff,replace=False)
            temp_df = temp_df.iloc[rand_list,:]
            return temp_df,rand_list
        elif len(temp_df.index) == 0:
            temp_df = this_fit[this_fit['Chi2pdfAvg'] == min(this_fit['Chi2pdfAvg'])]
            return temp_df,min(this_fit['Chi2pdfAvg'])
        else:
            return temp_df,'All'


    # def GetChiRange_fitr(self,chi_min=0,chi_max=100,nchi_cutoff = 100,fmted=True):
    #     this_data,rand_list = self.GetChiRange(chi_min=chi_min,chi_max=chi_max,nchi_cutoff = nchi_cutoff,fmted=fmted)
    #     index_list = list(this_data.index)
    #     index_list = [ival[:2] + ('fitr'+ival[2]+'-'+ival[3],) for ival in index_list]
    #     indicies = pa.MultiIndex.from_tuples(index_list,names=this_data.index.names[:2]+['fit_ranges'])
    #     output_data = pa.DataFrame(index=indicies)
    #     for icol,col_data in this_data.iteritems():
    #         output_data.loc[:,icol] = pa.Series(col_data.values,index=indicies)
    #     return output_data,rand_list

    def PullIdenticalParams(self):
        if self.VaryFunctions: return
        count = 0
        while 'Params' not in self.Fit_Stats.loc[:,'Fit'].iloc[count].fit_data.columns:
            if count+1 >= len(self.Fit_Stats.loc[:,'Fit']):
                print('Cannot pull identical parameters, none found')
                return
            count += 1
            if Debug_Mode:
                print(self.Fit_Stats.loc[:,'Fit'].iloc[count].fit_data.columns, len(self.Fit_Stats.loc[:,'Fit']))
                print(count,'Params' not in self.Fit_Stats.loc[:,'Fit'].iloc[count].fit_data.columns)
        self.par_list = list(self.Fit_Stats.loc[:,'Fit'].iloc[count].fit_data.loc[:,'Params'].index)
        for ipar in self.par_list:
            self.PullParam(ipar)

    def PullParam(self,ipar):
        if self.VaryFunctions:
            print('PullParams not implemented for different functions')
            return
        if 'Fit' in self.Fit_Stats:
            vall = []
            for ifit,fit_data in self.Fit_Stats.loc[:,'Fit'].items():
                if 'Params' not in fit_data.fit_data.columns:
                    vall.append(BootStrap(bootvals=[float('NaN')]*nboot))
                else:
                    if isinstance(ipar,int):
                        ipar = fit_data.fit_data.loc[:,'Params'].index[ipar]
                    vall.append(fit_data.fit_data['Params'][ipar])

            self.Fit_Stats.loc[:,'par_'+ipar] = pa.Series(vall,index=self.Fit_Stats.index)
            self.Fit_Stats.loc[:,'par_'+ipar+'_Avg'] = self.Fit_Stats.loc[:,'par_'+ipar].apply(lambda x : x.Avg)
            self.Fit_Stats.loc[:,'par_'+ipar+'_Std'] = self.Fit_Stats.loc[:,'par_'+ipar].apply(lambda x : x.Std)
            if 'Fit' in self.Fit_Stats_fmt:
                vall = []
                for ifit,fit_data in self.Fit_Stats_fmt.loc[:,'Fit'].items():
                    if 'Params' not in fit_data.fit_data.columns:
                        vall.append(BootStrap(bootvals=[float('NaN')]*nboot))
                    else:
                        if isinstance(ipar,int):
                            ipar = fit_data.fit_data.loc[:,'Params'].index[ipar]
                        vall.append(fit_data.fit_data['Params'][ipar])

                self.Fit_Stats_fmt.loc[:,'par_'+ipar] = pa.Series(vall,index=self.Fit_Stats_fmt.index)
                self.Fit_Stats_fmt.loc[:,'par_'+ipar+'_Avg'] = self.Fit_Stats_fmt.loc[:,'par_'+ipar].apply(lambda x : x.Avg)
                self.Fit_Stats_fmt.loc[:,'par_'+ipar+'_Std'] = self.Fit_Stats_fmt.loc[:,'par_'+ipar].apply(lambda x : x.Std)
        else:
            print('Warning, called PullParam before DoFits, running fits now')
            self.DoFits()
            self.PullParam(ipar)


    def Stats(self):
        self.Fit_Stats.loc[:,'Fit'].apply(lambda x : x.Stats())
        if 'Fit' in self.Fit_Stats_fmt:
            self.Fit_Stats_fmt.loc[:,'Fit'].apply(lambda x : x.Stats())

    def Eval_Function(self,*vals):
        return self.Fit_Stats.loc[:,'Fit'].apply(lambda x : x.Eval_Function(*vals))

    def Eval_Fun_Avg(self,*vals):
        return self.Fit_Stats.loc[:,'Fit'].apply(lambda x : x.Eval_Fun_Avg(*vals))

    def Eval_Fun_Std(self,*vals):
        return self.Fit_Stats.loc[:,'Fit'].apply(lambda x : x.Eval_Fun_Std(*vals))

    def Eval_Fun_AvgStd(self,*vals):
        return self.Fit_Stats.loc[:,'Fit'].apply(lambda x : x.Eval_Fun_AvgStd(*vals))


    def Eval_Function_fmt(self,*vals):
        return self.Fit_Stats_fmt.loc[:,'Fit'].apply(lambda x : x.Eval_Function(*vals))

    def Eval_Fun_Avg_fmt(self,*vals):
        return self.Fit_Stats_fmt.loc[:,'Fit'].apply(lambda x : x.Eval_Fun_Avg(*vals))

    def Eval_Fun_Std_fmt(self,*vals):
        return self.Fit_Stats_fmt.loc[:,'Fit'].apply(lambda x : x.Eval_Fun_Std(*vals))

    def Eval_Fun_AvgStd_fmt(self,*vals):
        return self.Fit_Stats_fmt.loc[:,'Fit'].apply(lambda x : x.Eval_Fun_AvgStd(*vals))

    def FormatIndex(self,this_index):
        if len(this_index) == 4:
            return this_index[0]+' fitr'+'-'.join(map(str,this_index[2:]))
        elif len(this_index) == 5:
            return this_index[0]+' fitr'+'-'.join(map(str,this_index[2:-1])) + ' '+this_index[-1]
        elif len(this_index) == 6:
            return this_index[0]+' fitr'+'-'.join(map(str,this_index[2:-2]))+ ' fittwor'+'-'.join(map(str,this_index[-2:]))

    def MakeFormatted(self):
        indexl = []
        if self.ndim_fit == 1:
            self.Fit_Col_Names_fmt = ['Fit Function','Fit Ranges']
            for ifun,iguess,ifit_min,ifit_max in self.Fit_Stats.index:
                indexl.append((ifun,'fit'+xmlfitr([ifit_min,ifit_max])))
        elif self.ndim_fit == 2 and len(self.Fit_Col_Names) == 6:
            self.Fit_Col_Names_fmt = ['Fit Function','Fit Ranges One','Fit Ranges Two']
            for ifun,iguess,ifit_min,ifit_max,ifit_min_2,ifit_max_2 in self.Fit_Stats.index:
                indexl.append((ifun,'fit'+xmlfitr([ifit_min,ifit_max]),'fittwo'+xmlfitr([ifit_min_2,ifit_max_2])))
        elif self.ndim_fit == 2 and len(self.Fit_Col_Names) == 5:
            self.Fit_Col_Names_fmt = ['Fit Function','Fit Ranges One','Shape']
            for ifun,iguess,ifit_min,ifit_max,ifit_min_2,ifit_max_2 in self.Fit_Stats.index:
                indexl.append((ifun,'fit'+xmlfitr([ifit_min,ifit_max]),'fittwo'+xmlfitr([ifit_min_2,ifit_max_2])))
        else:
            raise EnvironmentError('Error with making formatted fits')
        if len(indexl) > 0:
            indicies = pa.MultiIndex.from_tuples(indexl,names=self.Fit_Col_Names_fmt)
            self.Fit_Stats_fmt = pa.DataFrame(index=indicies)
            for icol,col_data in self.Fit_Stats.items():
                self.Fit_Stats_fmt.loc[:,icol] = pa.Series(col_data.values,index=indicies)

    def GetMaxFitr(self,fmted=True):
        fit_reset = self.Fit_Stats.reset_index()
        if 'fit_min_2' in self.Fit_Col_Names:
            pick_col = self.Fit_Col_Names[2:]
            fit_reset.loc[:,pick_col[0]] = fit_reset[pick_col[0]].apply(lambda x : int(x))
            fit_reset.loc[:,pick_col[1]] = fit_reset[pick_col[1]].apply(lambda x : int(x))
            fit_reset.loc[:,pick_col[2]] = fit_reset[pick_col[2]].apply(lambda x : int(x))
            fit_reset.loc[:,pick_col[3]] = fit_reset[pick_col[3]].apply(lambda x : int(x))
            fit_reset = fit_reset.loc[fit_reset[pick_col[0]] == min(fit_reset[pick_col[0]])]
            fit_reset = fit_reset.loc[fit_reset[pick_col[1]] == max(fit_reset[pick_col[1]])]
            fit_reset = fit_reset.loc[fit_reset[pick_col[2]] == min(fit_reset[pick_col[2]])]
            fit_reset = fit_reset.loc[fit_reset[pick_col[3]] == max(fit_reset[pick_col[3]])]
            fit_reset.set_index(self.Fit_Col_Names,inplace=True)
            if fmted:
                fitr = 'fitr'+str(fit_reset.index[0][2]) + '-'+str(fit_reset.index[0][3])
                fitr += '_fittwor'+str(fit_reset.index[0][4]) + '-'+str(fit_reset.index[0][5])
            else:
                fitr = fit_reset.index[0][2],fit_reset.index[0][3],fit_reset.index[0][4],fit_reset.index[0][5]
        else:
            pick_col = self.Fit_Col_Names[2:5]
            fit_reset.loc[:,pick_col[0]] = fit_reset[pick_col[0]].apply(lambda x : int(x))
            fit_reset.loc[:,pick_col[1]] = fit_reset[pick_col[1]].apply(lambda x : int(x))
            fit_reset = fit_reset.loc[fit_reset[pick_col[0]] == min(fit_reset[pick_col[0]])]
            fit_reset = fit_reset.loc[fit_reset[pick_col[1]] == max(fit_reset[pick_col[1]])]
            fit_reset.set_index(self.Fit_Col_Names,inplace=True)
            if fmted:
                fitr = 'fitr'+str(fit_reset.index[0][2]) + '-'+str(fit_reset.index[0][3])
            else:
                fitr = fit_reset.index[0][2],fit_reset.index[0][3]
        return fitr,fit_reset['Fit']



    def GetOutputDict(self):
        output_dict = ODNested()
        for ifit,fit_data in self.Fit_Stats.loc[:,'Fit'].items():
            this_fit = self.FormatIndex(ifit)
            if this_fit in list(output_dict.keys()):
                this_fit += '_2'
                while this_fit in list(output_dict.keys()):
                    this_fit = this_fit[:-1] + '_'+str(int(this_fit[-1])+1)
            output_dict[this_fit] = fit_data.GetOutputDict()
        return output_dict


    ## returns a series with all the chis combined together
    def GetChiVariation(self,chi_min=0,chi_max=100,nchi_cutoff = 100):
        chi_range,rand_list = self.GetChiRange(chi_min=0,chi_max=100,nchi_cutoff = 100)
        chi_range = chi_range['Fit']
        vall,indexl = [],[]
        for ipar in self.par_list:
            def this_fun(val):
                return val.fit_data.loc[ipar,'Params']
            this_data = chi_range.apply(this_fun)
            this_flat_data = np.array([ival.bootvals for ival in this_data.values]).flatten()
            this_nboot = len(this_flat_data)
            this_data = BootStrap(this_nboot,name='Combine Fits',bootvals=this_flat_data)
            vall.append(this_data)
            indexl.append(ipar)
        if len(indexl) > 0:
            out_df = pa.Series(vall,index=indexl).to_frame('Combine_Fits')
            return out_df,rand_list
        else:
            raise EnvironmentError('No parameters found ot do variation over')


    '''
    keying:

    dim1:
    will find fit_min and fit_max from string
    return dataframe of fitting classes

    dim2:
    fit_min , fit_max
    return dataframe of fitting classes

    dim3:
    function, fit_min , fit_max
    return dataframe of fitting classes

    dim4:
    function, iguess, fit_min , fit_max
    return dataframe of fitting classes

    dim5:
    function, iguess, fit_min , fit_max, parameter
    return bootstrap class of result
    '''

    def Get1D(self,key,fmted=True):
        if fmted:
            if len(key) == 1:
                return self.Fit_Stats_fmt.loc[[slice(None),key[0]],'Fit']
            elif len(key) == 2:
                return self.Fit_Stats_fmt.loc[key,'Fit']
            elif len(key) == 3:
                return self.Fit_Stats.loc[key[:1],'Fit'].fit_data.loc[key[1],'Params']
        else:
            if len(key) == 2:
                return self.Fit_Stats.loc[[slice(None),slice(None)] + list(key),'Fit']
            elif len(key) == 3:
                return self.Fit_Stats.loc[[key[0],slice(None)] + list(key)[1:],'Fit']
            elif len(key) == 4:
                return self.Fit_Stats.loc[key,'Fit']
            elif len(key) == 5:
                return self.Fit_Stats.loc[key[:4],'Fit'].fit_data.loc[key[4],'Params']


    def Get2D_shape(self,key,fmted=True):
        if fmted:
            if len(key) == 1:
                return self.Fit_Stats_fmt.loc[[slice(None),key[0],slice(None)],'Fit']
            elif len(key) == 2:
                return self.Fit_Stats_fmt.loc[[slice(None)] + key,'Fit']
            elif len(key) == 3:
                return self.Fit_Stats.loc[key,'Fit']
            elif len(key) == 4:
                return self.Fit_Stats.loc[key[:3],'Fit'].fit_data.loc[key[3],'Params']
        else:
            if len(key) == 2:
                return self.Fit_Stats.loc[[slice(None),slice(None)] + list(key)+[slice(None)],'Fit']
            elif len(key) == 3:
                return self.Fit_Stats.loc[[key[0],slice(None)] + list(key)[1:]+[slice(None)],'Fit']
            elif len(key) == 4:
                return self.Fit_Stats.loc[[key[0],slice(None)] + list(key)[1:],'Fit']
            elif len(key) == 5:
                return self.Fit_Stats.loc[key,'Fit']
            elif len(key) == 6:
                return self.Fit_Stats.loc[key[:5],'Fit'].fit_data.loc[key[5],'Params']

    def Get2D_box(self,key,fmted=True):
        if fmted:
            if len(key) == 1:
                return self.Fit_Stats_fmt.loc[[slice(None),key[0],slice(None)],'Fit']
            elif len(key) == 2:
                return self.Fit_Stats_fmt.loc[[slice(None)] + list(key),'Fit']
            elif len(key) == 3:
                return self.Fit_Stats.loc[key,'Fit']
            elif len(key) == 4:
                return self.Fit_Stats.loc[key[:3],'Fit'].fit_data.loc[key[3],'Params']
        else:
            if len(key) == 2:
                return self.Fit_Stats.loc[[slice(None),slice(None)] + list(key)+[slice(None),slice(None)],'Fit']
            elif len(key) == 3:
                return self.Fit_Stats.loc[[key[0],slice(None)] + list(key)[1:]+[slice(None),slice(None)],'Fit']
            elif len(key) == 5:
                return self.Fit_Stats.loc[[key[0],slice(None)] + list(key)[1:],'Fit']
            elif len(key) == 6:
                return self.Fit_Stats.loc[key,'Fit']
            elif len(key) == 7:
                return self.Fit_Stats.loc[key[:6],'Fit'].fit_data.loc[key[6],'Params']


    def __getitem__(self,key):
        if isinstance(key,str):
            key = key.split('_')
        if self.ndim_fit == 1:
            return self.Get1D(key,any(['fitr' in ikey for ikey in key]))
        elif self.ndim_fit == 2 and len(self.Fit_Col_Names) == 6:
            return self.Get2D_box(key,any(['fitr' in ikey for ikey in key]))
        elif self.ndim_fit == 2 and len(self.Fit_Col_Names) == 5:
            return self.Get2D_shape(key,any(['fitr' in ikey for ikey in key]))

    def Set1D(self,key,value,fmted=True):
        if fmted:
            if len(key) == 1:
                self.Fit_Stats_fmt.loc[[slice(None),key[0]],'Fit'] = value
            elif len(key) == 2:
                self.Fit_Stats_fmt.loc[key,'Fit'] = value
            elif len(key) == 3:
                self.Fit_Stats.loc[key[:1],'Fit'].fit_data.loc[key[1],'Params'] = value
        else:
            if len(key) == 2:
                self.Fit_Stats.loc[[slice(None),slice(None)] + list(key),'Fit'] = value
            elif len(key) == 3:
                self.Fit_Stats.loc[[key[0],slice(None)] + list(key)[1:],'Fit'] = value
            elif len(key) == 4:
                self.Fit_Stats.loc[key,'Fit'] = value
            elif len(key) == 5:
                self.Fit_Stats.loc[key[:4],'Fit'].fit_data.loc[key[4],'Params'] = value


    def Set2D_shape(self,key,value,fmted=True):
        if fmted:
            if len(key) == 1:
                self.Fit_Stats_fmt.loc[[slice(None),key[0],slice(None)],'Fit'] = value
            elif len(key) == 2:
                self.Fit_Stats_fmt.loc[[slice(None)] + key,'Fit'] = value
            elif len(key) == 3:
                self.Fit_Stats.loc[key,'Fit'] = value
            elif len(key) == 4:
                self.Fit_Stats.loc[key[:3],'Fit'].fit_data.loc[key[3],'Params'] = value
        else:
            if len(key) == 2:
                self.Fit_Stats.loc[[slice(None),slice(None)] + list(key)+[slice(None)],'Fit'] = value
            elif len(key) == 3:
                self.Fit_Stats.loc[[key[0],slice(None)] + list(key)[1:]+[slice(None)],'Fit'] = value
            elif len(key) == 4:
                self.Fit_Stats.loc[[key[0],slice(None)] + list(key)[1:],'Fit'] = value
            elif len(key) == 5:
                self.Fit_Stats.loc[key,'Fit'] = value
            elif len(key) == 6:
                self.Fit_Stats.loc[key[:5],'Fit'].fit_data.loc[key[5],'Params'] = value

    def Set2D_box(self,key,value,fmted=True):
        if fmted:
            if len(key) == 1:
                self.Fit_Stats_fmt.loc[[slice(None),key[0],slice(None)],'Fit'] = value
            elif len(key) == 2:
                self.Fit_Stats_fmt.loc[[slice(None)] + key,'Fit'] = value
            elif len(key) == 3:
                self.Fit_Stats.loc[key,'Fit'] = value
            elif len(key) == 4:
                self.Fit_Stats.loc[key[:3],'Fit'].fit_data.loc[key[3],'Params'] = value
        else:
            if len(key) == 2:
                self.Fit_Stats.loc[[slice(None),slice(None)] + list(key)+[slice(None),slice(None)],'Fit'] = value
            elif len(key) == 3:
                self.Fit_Stats.loc[[key[0],slice(None)] + list(key)[1:]+[slice(None),slice(None)],'Fit'] = value
            elif len(key) == 5:
                self.Fit_Stats.loc[[key[0],slice(None)] + list(key)[1:],'Fit'] = value
            elif len(key) == 6:
                self.Fit_Stats.loc[key,'Fit'] = value
            elif len(key) == 7:
                self.Fit_Stats.loc[key[:6],'Fit'].fit_data.loc[key[6],'Params'] = value

    def __setitem__(self,key,value):
        if isinstance(key,str):
            key = key.split('_')
        if self.ndim_fit == 1:
            self.Set1D(key,any(['fitr' in ikey for ikey in key]))
        elif self.ndim_fit == 2 and len(self.Fit_Col_Names) == 6:
            self.Set2D_box(key,any(['fitr' in ikey for ikey in key]))
        elif self.ndim_fit == 2 and len(self.Fit_Col_Names) == 5:
            self.Set2D_shape(key,any(['fitr' in ikey for ikey in key]))



    def __reversed__(self):
        self.Fit_Stats = self.Fit_Stats.iiloc[::-1]
        self.Fit_Stats_fmt = self.Fit_Stats_fmt.iiloc[::-1]
        return self

    def keys(self):
        return list(self.Fit_Stats['Fit'].keys())
    def values(self):
        return self.Fit_Stats['Fit'].values
    def items(self):
        return list(self.Fit_Stats['Fit'].items())
    def iterkeys(self):
        return iter(self.Fit_Stats['Fit'].keys())
    def itervalues(self):
        return self.Fit_Stats['Fit'].values

    def iteritems(self):
        return iter(self.Fit_Stats['Fit'].items())

    def iteritems_fmt(self):
        return iter(self.Fit_Stats_fmt['Fit'].items())


def TestSetOfFits():
    from BootStrapping import BootStrap
    from PredefFitFuns import LinearFitFun#,ConstantFitFun
    from PlotData import colourset8
    xvalssize = 5
    scale = xvalssize/3.
    nboot = 10
    values = np.random.normal(scale=scale,size=(xvalssize,nboot))
    values = np.array([ival + icv+100 for icv,ival in enumerate(values)])
    testdata = [BootStrap(nboot,name='test_bootstrap',cfgvals=ival,thisDelVal=False) for ival in values]

    testdf = pa.DataFrame()
    testdf.loc[:,'xdata'] = pa.Series([[np.arange(xvalssize)]])
    testdf.loc[:,'ydata'] = pa.Series([testdata])
    testdf.loc[:,'ydata_Avg'] = pa.Series([idata.Avg for idata in testdata])
    testdf.loc[:,'ydata_Std'] = pa.Series([idata.Std for idata in testdata])

    testdf_plot = pa.DataFrame()
    testdf_plot.loc[:,'xdata'] = pa.Series(np.arange(xvalssize))
    testdf_plot.loc[:,'ydata'] = pa.Series(testdata)
    testdf_plot.loc[:,'ydata_Avg'] = pa.Series([idata.Avg for idata in testdata])
    testdf_plot.loc[:,'ydata_Std'] = pa.Series([idata.Std for idata in testdata])

    testfit = SetOfFitFuns(data = testdf,name='Test')

    fit_info = {}
    fit_info['Funs'] = [LinearFitFun,2]
    # fit_info['Funs'] = [ConstantFitFun,1]
    # fit_info['iGuess'] = [0,0]
    testfit.ScanRange(0,len(values),fit_info,min_fit_len=2)
    testfit.DoFits()
    testfit.SortChi()

    # pl.clf()
    # pl.errorbar(testdf_plot['xdata'],testdf_plot['ydata_Avg'],testdf_plot['ydata_Std'],label='data')
    # testdf_plot.plot(x='xdata',y='ydata_Avg',label='data')
    this_data = pa.DataFrame()
    hold_series = pa.Series()
    hold_series['x_data'] = testdf_plot['xdata']
    hold_series['y_data'] = testdf_plot['ydata_Avg']
    hold_series['yerr_data'] = testdf_plot['ydata_Std']
    hold_series['fit_class'] = None
    hold_series['type'] = 'error_bar'
    hold_series['label'] = 'test data'
    hold_series['color'] = colourset8[0]
    hold_series['key_select'] = None
    hold_series['ShowPar'] = r'Par2'
    this_data.loc[:,'test_data'] = hold_series

    hold_series2 = pa.Series()
    hold_series2['type'] = 'fit_vary'
    fit_ranges = testfit.GetChiRange(fmted=True)[0].loc[:,'Fit']
    hold_series2['key_select'] = fit_ranges.index[0]
    hold_series2['fit_class'] = fit_ranges
    hold_series2['label'] = 'fit '
    hold_series2['ShowPar'] = r'Par2'
    this_data.loc[:,'test_fit_chi'] = hold_series2
    #
    #
    # for ic,(icol,(ifit,fit_data)) in enumerate(zip(colourset8[1:],testfit.GetChiRange(nchi_cutoff=3).iteritems())):
    #     hold_series = pa.Series()
    #     hold_series['type'] = 'fit'
    #     hold_series['fit_class'] = fit_data
    #     hold_series['label'] = testfit.FormatIndex(ifit)
    #     hold_series['color'] = icol
    #     this_data['test_fit_chi'+str(ic)] = hold_series



    this_info = pa.Series()
    this_info['save_file'] = './TestGraphs/TestFitChiPlot.pdf'
    this_info['title'] = 'Test Fitting'
    this_info['xlabel'] = 'Test x'
    this_info['ylabel'] = 'Test y'
    # this_info['xlims'] = [0,10]
    # this_info['ylims'] = [0,15]

    import PlotData as jpl
    data_plot = jpl.Plotting(plot_data=this_data,plot_info=this_info)
    data_plot.LoadPickle(DefWipe=True)
    data_plot.PrintData()
    data_plot.PlotAll()

    return testdata,testfit


def TestSetOfFits_2D():
    from BootStrapping import BootStrap
    from PredefFitFuns import LinearFF_2D#,ConstantFitFun
    from PlotData import colourset8
    xvalssize = (10,10)
    scale = np.prod(xvalssize)/3.
    nboot = 10
    values = np.random.normal(scale=scale,size= (np.prod(xvalssize),nboot))
    values = np.array([ival + icv+10 for icv,ival in enumerate(values)])
    testdata = [BootStrap(nboot,name='test_bootstrap',cfgvals=ival,thisDelVal=False) for ival in values]
    testdata = np.array(testdata).reshape(xvalssize)
    xdata_1 = np.array([list(range(xvalssize[1])) for _ in range(xvalssize[0])])
    xdata_2 = list(map(str,xdata_1.swapaxes(0,1).flatten()))
    xdata_1 = list(map(str,xdata_1.flatten()))
    xdata = [xdata_1,xdata_2]
    indicies = pa.MultiIndex.from_tuples(list(zip(*xdata)),names=['dim1','dim2'])
    this_plot_data = pa.Series(testdata.flatten(),index=indicies)
    # print 'DEBUG'
    # print this_data

    # testdf = pa.DataFrame()
    # testdf.loc[:,'xdata'] = pa.Series([[np.arange(xvalssize)]])
    # testdf.loc[:,'ydata'] = pa.Series([testdata])
    # testdf.loc[:,'ydata_Avg'] = pa.Series([idata.Avg for idata in testdata])
    # testdf.loc[:,'ydata_Std'] = pa.Series([idata.Std for idata in testdata])


    testfit = SetOfFitFuns(data = this_plot_data,name='Test')

    fit_info = {}
    fit_info['Funs'] = [LinearFF_2D,3]
    # fit_info['Funs'] = [ConstantFitFun,1]
    # fit_info['iGuess'] = [0,0]
    testfit.ScanBox(0,xvalssize[0],0,xvalssize[1],fit_info,min_fit_len_1=2,min_fit_len_2=2)
    testfit.DoFits()
    testfit.SortChi()

    # pl.clf()
    # pl.errorbar(testdf_plot['xdata'],testdf_plot['ydata_Avg'],testdf_plot['ydata_Std'],label='data')
    # testdf_plot.plot(x='xdata',y='ydata_Avg',label='data')
    this_data = pa.DataFrame()
    hold_series = pa.Series()
    hold_series['x_data'] = 'from_keys'
    hold_series['key_select'] = ('1',slice(None))
    hold_series['y_data'] = this_plot_data.apply(lambda x : x.Avg)
    hold_series['yerr_data'] = this_plot_data.apply(lambda x : x.Std)
    hold_series['fit_class'] = None
    hold_series['type'] = 'error_bar_vary'
    hold_series['label'] = 'test data'
    hold_series['color'] = colourset8[0]
    hold_series['ShowPar'] = r'Par2'
    hold_series['xaxis'] = None
    hold_series['otherXvals'] = None
    this_data.loc[:,'test_data'] = hold_series

    hold_series2 = pa.Series()
    hold_series2['type'] = 'fit_vary'
    hold_series2['key_select'] = 'First'
    hold_series2['fit_class'] = testfit
    hold_series2['label'] = 'fit '
    hold_series2['ShowPar'] = r'Par0'
    hold_series2['xaxis'] = 0
    hold_series2['otherXvals'] = [1]
    this_data.loc[:,'test_fit_chi'] = hold_series2
    #
    #
    # for ic,(icol,(ifit,fit_data)) in enumerate(zip(colourset8[1:],testfit.GetChiRange(nchi_cutoff=3).iteritems())):
    #     hold_series = pa.Series()
    #     hold_series['type'] = 'fit'
    #     hold_series['fit_class'] = fit_data
    #     hold_series['label'] = testfit.FormatIndex(ifit)
    #     hold_series['color'] = icol
    #     this_data['test_fit_chi'+str(ic)] = hold_series



    this_info = pa.Series()
    this_info['save_file'] = './TestGraphs/TestFitChiPlot_2D.pdf'
    this_info['title'] = 'Test Fitting 2D'
    this_info['xlabel'] = 'Test x[0]'
    this_info['ylabel'] = 'Test y'
    # this_info['xlims'] = [0,10]
    # this_info['ylims'] = [0,15]

    import PlotData as jpl
    data_plot = jpl.Plotting(plot_data=this_data,plot_info=this_info)
    data_plot.LoadPickle(DefWipe=True)
    data_plot.PrintData()
    data_plot.PlotAll()

    return testdata,testfit

    # testfit.PlotFunction()
    # pl.legend()
    # pl.savefig('./TestGraphs/TestFuns.pdf')

    #


if __name__ == '__main__':
    testdata,testfit = TestSetOfFits()
    testdata_2D,testfit_2D = TestSetOfFits_2D()
