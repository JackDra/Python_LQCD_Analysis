#!/usr/bin/env python

from copy import copy
from Params import defPrec,defMaxIter,fillalpha,TflowToPhys
from FileIO import ReadFuns,WriteFuns
from MiscFuns import Pullflag,ODNested
from scipy.optimize import leastsq
from BootStrapping import BootStrap
from XmlFormatting import AvgStdToFormat,MakeValAndErr,LegendFmt
from PredefFitFuns import LinearFitFun,SchiffFFDer,C2OSFAntiper

# import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
# import matplotlib.pyplot as pl
from PredefFitFuns import ConstFFDer,LinFFDer,C2OSFFDer,C2OSFFNoExpDer,C3OSFFDer
from PredefFitFuns import C2TSFFDer,C2TSFFNoExpDer,C3MomTSFFDer,C2OSFAntiperDer
from PredefFitFuns import TestTwoVarFFDer,FormFactorO2Der,FormFactorO1Der,FormFactorO3Der
from PredefFitFuns import DPfitfunOneParDer,ParmDivXDer,ChitFitFunDer,ParmDivXPDer,LinFFDer_x0
from PredefFitFuns import c3PFDer,c3PFDer2,c3PFDer3,c3PFDer_skip1,c3PFDer2_skip1,c3PFDer3_skip1
from PredefFitFuns import c3EFFDer,c3EFFDer2,c3EFFDer_nosqrt,c3EFFDer2_nosqrt,c3EFFDer3_nosqrt
from PredefFitFuns import DPfitfun2Der,DPfitfunDer,OORNFFDer,SquaredFFDer,SquaredFFDer_x0
from PredefFitFuns import Chi_Prop_Der,Chi_Prop_Der_2exp,LinFFDer_2D,LogFFDer,LogFFDer_x0,c3FFDer
from PredefFitFuns import Alpha_Prop_Der,Alpha_Prop_Der_2exp,RF_Prop_Der,Alpha_Prop_Der_plin,c3FFDer_NoA
from PredefFitFuns import SqPlLogFFDer,SqPlLogFFDer_x0,c3ContFFDer
from PredefFitFuns import ContFFDer,ContFFDer_x0,ContPlLogFFDer,ContPlLogFFDer_x0
from PredefFitFuns import ContFFDer_mix,ContFFDer_x0_mix,ContPlLogFFDer_mix,ContPlLogFFDer_x0_mix
from MomParams import LatticeParameters

import numpy as np
import pandas as pa
import operator as op
from MiscFuns import flip_fun_2arg
from warnings import warn
err_thresh = 0.00001
"""

This class has methods for doing least squares fitting on regular and bootstrapped (see BootStrap.py) quantities.

Fitting Contains:
 FitFunction which performs a least squares fit.
 FitBoot which performs a least squares fit to the bootstrapped data.
 Lots of functions which are predefied for fitting.
 Initial conditions for all fit functions


FitFuns is a class containing all pre-defined fit functions.
Add to this to add predefined functions for fitting
Make sure all functions used for fitting must have form:
 f(x,p), where:
 x is an array of independant variables
 p is an array of fit parameters to extract
 (look at my defined functions to see how).

"""


def PullFitSeries(this_series):
    if len(this_series) == 0:
        return this_series
    outl,indexl = [],[]
    for ikey,NNQ_data in this_series.items():
        if 'Params' not in NNQ_data.fit_data: continue
        for fit_key,fit_data in NNQ_data.fit_data['Params'].items():
            indexl.append(ikey+(fit_key,))
            outl.append(fit_data)
    if len(indexl) > 0:
        if this_series.index.names is None:
            print(this_series)
            raise EnvironmentError('No labels set up for series')
        indicies = pa.MultiIndex.from_tuples(indexl,names=this_series.index.names + ['fit parameter'])
        return pa.Series(outl,index=indicies).dropna()
    else:
        return pa.Series()



class Fitting(object):

    ## pass length 2 or 3 list data = [ xdata , ydata , yerr ] or [ xdata , ydata ]
    ## pass length 1 or 2 array Funs = [ Function , parlen ]  or [ Function , parlen , Derivative of Function ]
    ## pass n dimensional inintial guess iGuess MAKE SURE IT IS SAME DIMENSNION AS PARAMETERS IN FUNCTION

    ## pass in 'Estimate' into Derivative of Function if you explicitly want to calculate the derivative numerically
    ## or else the code will attempt to find the corresponding derivative function for the passed in Function.
    # linestyles = ['-',':','-.','--']
    linestyles = ['-']
    currlinestyle = 0
    hairline_dropoff = 20
    hairline_alpha = 1.0

    parentlist = []

    def __init__(self,thisPrec=defPrec,thisMI=defMaxIter,data = 'None', Funs='None',iGuess ='None',
                 name='None',paramlist='PreDef',paramunits='None'):

        self.FPUnits = ''
        self.FPName = ''
        self.linestyle = Fitting.linestyles[Fitting.currlinestyle]
        Fitting.currlinestyle = (Fitting.currlinestyle + 1) % len(Fitting.linestyles)
        self.hairline_dropoff = Fitting.hairline_dropoff
        self.hairline_alpha = Fitting.hairline_alpha
        self.MaxIter = thisMI
        self.Prec = thisPrec
        self.iGuess = iGuess
        self.Fun = 'None'
        self.name = name
        self.fit_data = pa.DataFrame()
        '''
        self.fit_data = DataFrame

        columns =   Params,     CoVar_thispar,      Chi2Dof
                    ParamsAvg,  CoVarAvg_thispar,   Chi2DofAvg
                    ParamsStd,  CoVarStd_thispar,   Chi2DofStd
                    Params_units
        rows = Fit Parameters
        '''


        if Funs == 'None':
            self.Fun = 'None'
            self.FunDer = 'None'
            self.parlen = 1
        else:
            self.ImportFunction(*Funs,thisparlab=paramlist,thisparunits=paramunits)


        if isinstance(data,pa.DataFrame):
            self.ImportData_DF(data)
        elif data == 'None':
            ## also initialises :)
            self.RemoveData()
            self.xbooted = False
        else:
            self.ImportData(*data)


        # print 'TypeTest2 ',type(self.ydata[0])

        self.SetParamLabels(paramlist,paramunits)

        self.ChiFun = self.ChiOfFun()
        self.ChiDer = self.ChiOfDer()



        # ## Params { fit_parameter_label , BootStrap }
        # self.Params = 'Not Computed'
        # self.Covar = 'Not Computed'
        # self.Chi2DoF = 'Not Computed'
        #
        # ## ParamsAvg = { fit_parameter_label , average }
        # self.ParamsAvg = 'Not Computed'
        # self.CovarAvg = 'Not Computed'
        # self.Chi2DoFAvg = 'Not Computed'
        #
        #
        # ## ParamsAvg = { fit_parameter_label , std }
        # self.ParamsStd = 'Not Computed'
        # self.CovarStd = 'Not Computed'
        # self.Chi2DoFStd = 'Not Computed'

    # def linestyleNext(self):
    #     output = copy(self.thislinestyle)
    #     nextindex = Fitting.linestyles.index(self.thislinestyle)+1
    #     if nextindex >= len(Fitting.linestyles): nextindex = 0
    #     self.thislinestyle = Fitting.linestyles[nextindex]
    #     return output

    def GetFuns(self):
        if hasattr(self,'Fun_name'):
            self.Fun = ReadFuns(self.Fun_name)[0]
            self.ChiFun = self.ChiOfFun()
            del self.Fun_name
        if hasattr(self,'FunDer_name'):
            self.FunDer = ReadFuns(self.FunDer_name)[0]
            self.ChiDer = self.ChiOfDer()
            del self.FunDer_name

    def RemoveFuns(self):
        if hasattr(self,'Fun') and hasattr(self,'FunDer'):
            if isinstance(self.Fun,str):
                self.Fun_name = str(self.Fun)
                self.FunDer_name = str(self.FunDer)
                self.GetFuns()
            elif hasattr(self.Fun,'file_name'):
                self.Fun_name = self.Fun.file_name
            else:
                self.Fun_name = self.Fun.__name__
            if isinstance(self.FunDer,str):
                self.FunDer_name = str(self.FunDer)
            elif hasattr(self.FunDer,'file_name'):
                self.FunDer_name = self.FunDer.file_name
            else:
                self.FunDer_name = self.FunDer.__name__
            WriteFuns(self.Fun,self.FunDer)
            del self.Fun
            del self.FunDer
            del self.ChiFun
            del self.ChiDer

    def CheckParams(self):
        if 'Params' not in self.fit_data:
            raise EnvironmentError(self.name + ' has does not Params')


    def Stats(self):
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
        self.CheckParams()
        if len(self.fit_data.index) > 0:
            hold_dataframe = pa.DataFrame(columns=self.fit_data.columns,index=self.fit_data.index)
            for ival,fitval in self.fit_data['Params'].items():
                hold_dataframe.at[ival,'Params'] = DoStats(fitval)
                hold_dataframe.at[ival,'ParamsAvg'] = DoAvg(fitval)
                hold_dataframe.at[ival,'ParamsStd'] = DoStd(fitval)
            self.fit_data.update(hold_dataframe)

    def SetGuess(self,iGuess='None'):
        if iGuess == 'None':
            self.iGuess = self.GuessFromFun()
            if self.iGuess == 'None':
                self.iGuess = [1.0]*self.parlen
        else:
            self.iGuess = iGuess

    def SetParamLabels(self,paramlist='PreDef',paramunits='None'):
        if isinstance(paramlist,str) and paramlist == 'PreDef':
            self.paramlist = self.ParamLabFromFun()
        else:
            self.paramlist = paramlist
        if isinstance(self.paramlist,str) and self.paramlist == 'Def':
            self.paramlist = ['Par'+str(i) for i in range(1,self.parlen+1)]
        if len(self.paramlist) != self.parlen:
            raise IOError('parameter list of labels '+','.join(self.paramlist)+ ' does not match parameter length '+str(self.parlen))
        if paramunits == 'None':
            self.paramunits = pa.Series(['' for _ in range(len(self.paramlist))] , index=self.paramlist)
        elif not isinstance(paramunits,pa.Series):
            self.paramunits = pa.Series(paramunits,index=self.paramlist)
        self.fit_data.loc[:,'Params_units'] = self.paramunits.copy()


    def ImportData_DF(self,this_df):
        self.ImportData(this_df['xdata'][0],this_df['ydata'][0])

    ## xdata [ param , idata ] BootStrap (param size can be larger than 1 to do multi dimensional fitting)
    ## ydata [ idata ] BootStrap
    def ImportData(self,xdata,ydata):
        self.xdata = np.array(xdata)
        ## if someone accidently forgets to add a dimension for multi dimentional fitting
        if self.xdata.ndim == 1:
            self.xdata = np.array([self.xdata])
        self.xbooted = all([isinstance(ixval,BootStrap) for ixval in self.xdata[0]])

        self.ydata = np.array(ydata)
        if all([isinstance(iyval,BootStrap) for iyval in ydata]):
            self.yerr,self.yAvg = [],[]
            for iyval in ydata:
                self.yerr.append(iyval.Std)
                self.yAvg.append(iyval.Avg)
            self.yAvg,self.yerr = np.array(self.yAvg),np.array(self.yerr)
        else:
            ## if no error is provided, error is just set to 1 (no weightings).
            self.yAvg = self.ydata
            self.yerr = [1.0]*len(ydata)
        if any(self.yerr/np.abs(self.yAvg) < err_thresh):
            warn('some errors were under the error theshold, removing weighted fit')
            self.yerr = np.ones(self.yerr.shape)

    ## setting FunDer to 'None' implies estimating the derivative numerically
    def ImportFunction(self,Function,thisparlen,FunDer='None',thisparlab = 'PreDef',thisparunits='None'):
        self.Fun = Function
        self.parlen = thisparlen
        if FunDer == 'Estimate':
            self.FunDer = 'Estimate'
        elif FunDer == 'None':
            self.FunDer = self.DerOfFun()
        else:
            self.FunDer = FunDer
        self.ChiFun = self.ChiOfFun()
        self.ChiDer = self.ChiOfDer()
        if self.name == 'None':
            self.name = self.Fun.__name__
        self.SetParamLabels(paramlist = thisparlab,paramunits=thisparunits)


    def RemoveData(self):
        self.xdata = np.array([])
        self.ydata = np.array([])
        self.yAvg = np.array([])
        self.yerr = np.array([])


    def ChiOfFun(self):
        if self.Fun == 'None':
            return 'None'
        else:
            def LSFun(par,val):
                xval = np.array(val[:-2])
                # DEBUG THIS THING IS FINIKY, WITH SOME OTHER FITTING STUFF
                # if len(xval) == 1:
                #     xval = xval[0]
                yval = np.array(val[-2])
                errval = np.array(val[-1])
                return (np.array(self.Fun(xval,par))-yval)/errval
            return LSFun



    def ChiOfDer(self):
        if self.FunDer == 'Estimate':
            return None
        else:
            def LSDerFun(par,val):
                xval = np.array(val[:-2])
                ## DEBUG THIS THING IS FINIKY, WITH SOME OTHER FITTING STUFF
                # if len(xval) == 1:
                #     xval = xval[0]
                # yval = val[-2]
                errval = val[-1]
                return np.transpose(np.array(self.FunDer(xval,par))/errval)
            return LSDerFun


    def UnExpParams(self,paramkeys='All'):
        if paramkeys == 'All':
            logparam,logavg,logstd = [],[],[]
            for parlab,ifit_data in self.fit_data.Params.iterrows():
                iparam = ifit_data['Params'].Exp()
                iparam.Stats()
                logparam.append(iparam)
                logavg.append(iparam.Avg)
                logstd.append(iparam.Std)
            self.fit_data.loc[:,'Params'] = pa.Series(logparam,index=self.fit_data.index)
            self.fit_data.loc[:,'ParamsAvg'] = pa.Series(logavg,index=self.fit_data.index)
            self.fit_data.loc[:,'ParamsStd'] = pa.Series(logstd,index=self.fit_data.index)
        else:
            # self.Params,self.ParamsAvg,self.ParamsStd = OrderedDict(),OrderedDict(),OrderedDict()
            pser = self.fit_data.loc[:,'Params'].copy()
            pavg = self.fit_data.loc[:,'ParamsAvg'].copy()
            pstd = self.fit_data.loc[:,'ParamsStd'].copy()
            for ival in paramkeys:
                if ival in pser.index:
                    expval = pser[ival].Exp()
                    expval.Stats()
                    pser[ival] = expval
                    pavg[ival] = expval.Avg
                    pstd[ival] = expval.Std
            self.fit_data.loc[:,'Params'] = pser
            self.fit_data.loc[:,'ParamsAvg'] = pavg
            self.fit_data.loc[:,'ParamsStd'] = pstd

    def LSFit(self,thisxdata,thisydata,thisyerr):

        data = np.append(thisxdata,[thisydata,thisyerr],axis=0)
        self.SetGuess(self.iGuess)
        if len(self.iGuess) != self.parlen:
            # raise IOError('Initial guess was not initiated properly.')
            print('Warning, Initial guess was not initiated properly, attempting to run with predefined initial conditions')
            self.SetGuess()

        # Debugging:
        # print()
        # print('chifun' , self.ChiFun)
        # print('iGuess' , self.iGuess)
        # print('data' , data)
        # print('ChiDer' , self.ChiDer)
        # print('MaxIter' , self.MaxIter)
        # print('Prec' , self.Prec)
        # print()
        if len(thisydata) < self.parlen:
            return float('NaN'),float('NaN'),float('NaN')

        x,covar, infodict, mesg, ier = leastsq(self.ChiFun,self.iGuess,args=data, Dfun=self.ChiDer, maxfev=self.MaxIter , xtol=self.Prec, full_output=1)
        DoF = len(thisydata) - len(thisxdata)
        if DoF == 0:
            chisqpdf = float('NaN')
        else:
            chisqpdf=sum(infodict["fvec"]*infodict["fvec"])/float(DoF)
        # if ier != 1:
        #     print x,covar
        #     raise ValueError, "Optimal parameters not found: " + mesg
        return x,covar,chisqpdf




    ##Note, xdatain must be in the form
    ## xdata [ [x1values] , [x2values] , ..... [xnvalues] ]
    ## for fitting functions over variables F( x1, x2, ... xn )
    ## xdata [ [ [x1values] , [x2values] , ..... [xnvalues] ]_Average , [ [x1values] , [x2values] , ..... [xnvalues] ]_boot1 ....  ]
    ## for bootstrapped x values
    def FitBoots(self,DefWipe=True):
        if 'Params' in self.fit_data:
            return False
        self.GetFuns()
        def GetBootIteration(data,index=1):
            old_shape = list(data.shape)
            return np.rollaxis(np.array([idata.bootvals for idata in data.flatten()]).reshape(*(old_shape+[self.nboot])),index)

        def GetXavg():
            if self.xbooted:
                return np.array(Pullflag(self.xdata,'Avg'))
            else:
                return self.xdata
        # print 'TypeTest3 ',type(self.ydata[0])

        # print GetXavg()
        # print Pullflag(self.ydata,'Avg')
        # print self.yerr
        [thisavg,thisCVavg,thisCDFavg] = self.LSFit(GetXavg(),Pullflag(self.ydata,'Avg'),self.yerr)
        if isinstance(thisavg,float) and np.isnan(thisavg):
            return float('NaN')
        self.fit_data.loc[:,'ParamsAvg'] = pa.Series(thisavg,index=self.paramlist)
        self.nboot = self.ydata.flatten()[0].nboot
        for icp,ipar in enumerate(self.paramlist):
            if isinstance(thisCVavg,np.ndarray):
                self.fit_data.loc[:,'CovarAvg_'+ipar] = pa.Series(thisCVavg[icp],index=self.paramlist)
            else:
                self.fit_data.loc[:,'CovarAvg_'+ipar] = pa.Series(np.nan,index=self.paramlist)
        if not isinstance(thisCDFavg,np.ndarray):
            self.fit_data.loc[:,'Chi2DofAvg'] = pa.Series(np.nan,index=self.paramlist)
        else:
            self.fit_data.loc[:,'Chi2DofAvg'] = pa.Series(thisCDFavg,index=self.paramlist)
        tempPar,tempCovar,tempChi2 = [],[],[]
        thisnboot = 0
        CovarError = False
        for iboot,bootdata in enumerate(GetBootIteration(self.ydata)):
            thisnboot += 1
            if self.xbooted:
                tempboot = self.LSFit(GetBootIteration(self.xdata,index=2)[iboot],bootdata,self.yerr)
            else:
                tempboot = self.LSFit(self.xdata,bootdata,self.yerr)
            tempChi2.append(tempboot[2])
            tempPar.append([])
            if isinstance(tempboot[1],np.ndarray):
                tempCovar.append([])
                for ip,ic in zip(*tempboot[:2]):
                    tempPar[-1].append(ip)
                    tempCovar[-1].append(ic)
            else:
                CovarError = True
                for ip in tempboot[0]:
                    tempPar[-1].append(ip)
                tempCovar.append(tempboot[1])

        self.fit_data.loc[:,'Chi2DoF'] = pa.Series(BootStrap(thisnboot,bootvals=tempChi2),index=self.paramlist)
        self.fit_data.loc[:,'Chi2DoFStd'] = pa.Series(self.fit_data.Chi2DoF[0].Std,index=self.paramlist)

        if CovarError:
            thisPar,thisParStd = [],[]
            for ip in np.rollaxis(np.array(tempPar),1):
                thisPar.append(BootStrap(thisnboot,bootvals=ip))
                thisParStd.append(thisPar[-1].Std)
            self.fit_data.loc[:,'Covar'] = pa.Series(np.nan,index=self.paramlist)
            self.fit_data.loc[:,'CovarStd'] = pa.Series(np.nan,index=self.paramlist)
            self.fit_data.loc[:,'Params'] = pa.Series(thisPar,index=self.paramlist)
            self.fit_data.loc[:,'ParamsStd'] = pa.Series(thisParStd,index=self.paramlist)
        else:
            thisPar,thisParStd = [],[]
            for icp,(ip,ic,ipar) in enumerate(zip(np.rollaxis(np.array(tempPar),1),np.rollaxis(np.array(tempCovar),1),self.paramlist)):
                thisPar.append(BootStrap(thisnboot,bootvals=ip))
                thisParStd.append(thisPar[-1].Std)
                thisCo,thisCoStd = [],[]
                for jcp,jc in enumerate(np.rollaxis(ic,1)):
                    thisCo.append(BootStrap(thisnboot,bootvals=jc))
                    thisCoStd.append(thisCo[-1].Std)
                self.fit_data.loc[:,'Covar_'+ipar] = pa.Series(thisCo,index=self.paramlist)
                self.fit_data.loc[:,'CovarStd_'+ipar] = pa.Series(thisCoStd,index=self.paramlist)
            self.fit_data.loc[:,'Params'] = pa.Series(thisPar,index=self.paramlist)
            self.fit_data.loc[:,'ParamsStd'] = pa.Series(thisParStd,index=self.paramlist)
        # if DefWipe: self.RemoveData()
        self.RemoveFuns()
        return True

    def Get_Chi2DoF(self):
        return self.fit_data['Chi2DoF'].iloc[0]

    def Eval_Function(self,*vals):
        self.GetFuns()
        outboot = self.Fun(vals,np.array(self.fit_data['Params'].values))
        if hasattr(self,'scale_y') and self.scale_y != 1.0:
            outboot = outboot * self.scale_y
        outboot.Stats()
        # self.RemoveFuns()
        return outboot

    def Eval_Fun_Avg(self,*vals):
        return self.Eval_Function(*vals).Avg


    def Eval_Fun_Std(self,*vals):
        return self.Eval_Function(*vals).Std

    def Eval_Fun_AvgStd(self,*vals):
        output = self.Eval_Function(*vals)
        return output.Avg,output.Std

    def PlotHairlineFun(self,otherXvals=[],xaxis=0,xdatarange='Data',thislineres=100,
                        color=None,shift=0.0,flowphys=False,ShowPar='None',ShowEval=None,
                        x_scale=1.0,y_scale=1.0,plot_plane=None,suppress_key=False,
                        x_fun=None,supress_legend=False,extrap_fade=1.5):
        self.GetFuns()
        if y_scale != 1.0:
            self.scale_y = y_scale
        self.scale_x = x_scale
        if self.scale_x == 0:
            self.scale_x = 1

        if self.xdata.shape[0] == 1:
            thisxlist = self.xdata[0]
        else:
            if not isinstance(otherXvals,(tuple,list,np.ndarray)):
                print(otherXvals)
                raise EnvironmentError('otherXvals must be list')
            if len(otherXvals) == 0:
                print(otherXvals,' is invalid for otherXvals, defaulting to first value')
                first_Xvals = copy(np.swapaxes(self.xdata,0,1)[0])
                first_Xvals = np.delete(first_Xvals,xaxis)
                return self.PlotHairlineFun(otherXvals=first_Xvals,xaxis=xaxis,xdatarange=xdatarange,
                        thislineres=thislineres,x_fun=x_fun,x_scale=x_scale,y_scale=y_scale,
                        color=color,shift=shift,flowphys=flowphys,ShowPar=ShowPar,ShowEval=ShowEval,
                        plot_plane=plot_plane,suppress_key=suppress_key,extrap_fade=extrap_fade)
            otherXvals = list(map(float,otherXvals))
            thisxlist = []
            for icx,ix in enumerate(np.swapaxes(self.xdata,0,1)):
                thisix = copy(ix)
                thisix = np.delete(thisix,xaxis)
                if ix[xaxis] not in thisxlist:
                    thisxlist.append(ix[xaxis])
        pre_extrap,post_extrap = [],[]
        if xdatarange == 'Data':
            xdatarange = min(thisxlist),max(thisxlist)
        elif xdatarange[-1] == 'ToDataMax':
            xdatarange = float(xdatarange[0]),max(thisxlist)
            if xdatarange[0] < min(thisxlist):
                pre_extrap = np.linspace(min(thisxlist),xdatarange[0],thislineres)
        elif xdatarange[0] == 'FromDataMin':
            xdatarange = min(thisxlist),float(xdatarange[-1])
            if xdatarange[-1] > max(thisxlist):
                post_extrap = np.linspace(xdatarange[-1],max(thisxlist),thislineres)
        else:
            xdatarange = float(xdatarange[0]),float(xdatarange[-1])
            if xdatarange[-1] > max(thisxlist):
                post_extrap = np.linspace(xdatarange[-1],max(thisxlist),thislineres)
            if xdatarange[0] < min(thisxlist):
                pre_extrap = np.linspace(min(thisxlist),xdatarange[0],thislineres)

        xplot = np.linspace(min(thisxlist),max(thisxlist),thislineres)
        yplotBoot,yplotUp,yplotDown,yplotAvg,yplotStd = [],[],[],[],[]
        yplotBoot_pre,yplotUp_pre,yplotDown_pre,yplotAvg_pre,yplotStd_pre = [],[],[],[],[]
        yplotBoot_post,yplotUp_post,yplotDown_post,yplotAvg_post,yplotStd_post = [],[],[],[],[]
        for ix in xplot:
            ndxval = np.insert(otherXvals,xaxis,ix)
            evalBoot = self.Eval_Function(*ndxval)
            yplotAvg.append(evalBoot.Avg)
            yplotStd.append(evalBoot.Std)
            yplotUp.append(yplotAvg[-1] + yplotStd[-1])
            yplotDown.append(yplotAvg[-1] - yplotStd[-1])
            yplotBoot.append(evalBoot.bootvals.values)
        for ix in pre_extrap:
            ndxval = np.insert(otherXvals,xaxis,ix)
            evalBoot = self.Eval_Function(*ndxval)
            yplotAvg_pre.append(evalBoot.Avg)
            yplotStd_pre.append(evalBoot.Std)
            yplotUp_pre.append(yplotAvg_pre[-1] + yplotStd_pre[-1])
            yplotDown_pre.append(yplotAvg_pre[-1] - yplotStd_pre[-1])
            yplotBoot_pre.append(evalBoot.bootvals.values)
        for ix in post_extrap:
            ndxval = np.insert(otherXvals,xaxis,ix)
            evalBoot = self.Eval_Function(*ndxval)
            yplotAvg_post.append(evalBoot.Avg)
            yplotStd_post.append(evalBoot.Std)
            yplotUp_post.append(yplotAvg_post[-1] + yplotStd_post[-1])
            yplotDown_post.append(yplotAvg_post[-1] - yplotStd_post[-1])
            yplotBoot_post.append(evalBoot.bootvals.values)
        yplotBoot = np.swapaxes(np.array(yplotBoot),0,1)
        yplotBoot_pre = np.swapaxes(np.array(yplotBoot_pre),0,1)
        yplotBoot_post = np.swapaxes(np.array(yplotBoot_post),0,1)
        # if scale is not None or scale != 1.0:
        #     xplot = scale*xplot
        if flowphys != False:
            xplot = TflowToPhys(xplot,flowphys)
            pre_extrap = TflowToPhys(pre_extrap,flowphys)
            post_extrap = TflowToPhys(post_extrap,flowphys)
        if suppress_key:
            thislab = ''
        else:
            thislab = self.name + ' '
        if ShowPar in self.fit_data.index:
            thisunits = ''
            if len(self.FPName) > 0:
                this_SP = self.FPName
            else:
                this_SP = ShowPar
            if len(self.FPUnits) > 0:
                thisunits = self.FPUnits
            else:
                if ShowPar in self.paramunits.index:
                    thisunits = self.paramunits[ShowPar]

            thislab += ''.join([this_SP,'=',MakeValAndErr(self.fit_data.Params[ShowPar].Avg,
                                                          self.fit_data.Params[ShowPar].Std,
                                                          Dec=2,latex=True),thisunits])
            # if ShowPar in self.paramunits.index:
            #     thisunits = self.paramunits[ShowPar]
            # thislab += ''.join([ShowPar,'=',MakeValAndErr(self.fit_data.Params[ShowPar].Avg,self.fit_data.Params[ShowPar].Std,Dec=2,latex=True),thisunits])
            # print 'REMOVE THIS AFTER'
            # scale_param = self.fit_data.Params[ShowPar]*scale
            # scale_param.Stats()
            # thislab += ' $'+' '.join([ShowPar,'=',MakeValAndErr(scale_param.Avg,scale_param.Std,Dec=3),thisunits])+'$'
        if ShowEval is not None and ShowEval is not False:
            if isinstance(ShowEval,(list,tuple,np.ndarray)):
                indep_var = ','.join(['{:.2f}'.format(ival) for ival in ShowEval])
            else:
                indep_var = '{:.2f}'.format(ShowEval)
                ShowEval = [ShowEval]
            this_eval = self.Eval_Fun_AvgStd(*ShowEval)
            thislab += ''.join(['f(',indep_var,')=',MakeValAndErr(this_eval[0],this_eval[1],Dec=2,latex=True)])

        if supress_legend:
            thislab = None
        else:
            thislab = str(LegendFmt(thislab))
        if plot_plane is None:
            raise EnvironmentError('Please use plot_plane implementataion.')
            # pl.plot(np.array(xplot)+shift,yplotAvg,label=thislab,color=color,linestyle=self.linestyle)
            # pl.fill_between(np.array(xplot)+shift,yplotUp,yplotDown,alpha=fillalpha,edgecolor='none',color=color)
        else:
            if self.scale_x != 1:
                if x_fun is not None and x_fun is not False:
                    def this_xfun(val):
                        return x_fun(val)*self.scale_x
                else:
                    def this_xfun(val):
                        return val*self.scale_x
            else:
                this_xfun = x_fun
            if this_xfun is not None and this_xfun is not False:
                xplot = np.array([this_xfun(ix) for ix in xplot])
                pre_extrap = np.array([this_xfun(ix) for ix in pre_extrap])
                post_extrap = np.array([this_xfun(ix) for ix in post_extrap])

            if ShowEval is not None and ShowEval is not False:
                if this_xfun is not None and this_xfun is not False:
                    plot_plane.errorbar( this_xfun(ShowEval[xaxis]),[self.Eval_Fun_Avg(*ShowEval)],[self.Eval_Fun_Std(*ShowEval)],label=thislab,color=color,fmt='x',markersize='10')
                else:
                    plot_plane.errorbar( ShowEval[xaxis],[self.Eval_Fun_Avg(*ShowEval)],[self.Eval_Fun_Std(*ShowEval)],label=thislab,color=color,fmt='x',markersize='10')
                plot_plane.plot(    np.array(xplot)+shift,yplotAvg,color=color,linestyle=self.linestyle)
                if len(pre_extrap) > 0:
                    plot_plane.plot(    np.array(pre_extrap)+shift,yplotAvg_pre,color=color,linestyle='--')
                if len(post_extrap) > 0:
                    plot_plane.plot(    np.array(post_extrap)+shift,yplotAvg_post,color=color,linestyle='--')
            else:
                plot_plane.plot(    np.array(xplot)+shift,yplotAvg,
                                label=thislab,color=color,linestyle=self.linestyle)
                if len(pre_extrap) > 0:
                    plot_plane.plot(    np.array(pre_extrap)+shift,yplotAvg_pre,
                                    color=color,linestyle='--')
                if len(post_extrap) > 0:
                    plot_plane.plot(    np.array(post_extrap)+shift,yplotAvg_post,
                                    color=color,linestyle='--')
            plot_plane.plot(np.array(xplot)+shift,yplotUp,
                                color=color,linestyle='--')
            plot_plane.plot(np.array(xplot)+shift,yplotDown,
                                color=color,linestyle='--')
            if len(pre_extrap) > 0:
                plot_plane.plot(np.array(pre_extrap)+shift,yplotUp_pre,
                                    color=color,linestyle='--')
                plot_plane.plot(np.array(pre_extrap)+shift,yplotDown_pre,
                                    color=color,linestyle='--')
            if len(post_extrap) > 0:
                plot_plane.plot(np.array(post_extrap)+shift,yplotUp_post,
                                    color=color,linestyle='--')
                plot_plane.plot(np.array(post_extrap)+shift,yplotDown_post,
                                    color=color,linestyle='--')
            yplotAvg,yplotStd = np.array(yplotAvg),np.array(yplotStd)
            for iyboot in yplotBoot:
                id = np.nanmean(np.abs((iyboot-yplotAvg)/yplotStd))
                this_alpha = self.hairline_alpha/(1+id*self.hairline_dropoff)
                # print 'DEBUG',id,this_alpha
                # for iyb,iyavg in zip(iyboot,yplotAvg):
                #     print iyb,iyavg
                plot_plane.plot(np.array(xplot)+shift,iyboot,alpha=this_alpha,
                                    color=color,linestyle=self.linestyle)

            if len(pre_extrap) > 0:
                yplotAvg_pre,yplotStd_pre = np.array(yplotAvg_pre),np.array(yplotStd_pre)
                for iyboot in yplotBoot_pre:
                    id = np.nanmean(np.abs((iyboot-yplotAvg_pre)/yplotStd_pre))
                    this_alpha = self.hairline_alpha/(1+id*self.hairline_dropoff)
                    # print 'DEBUG',id,this_alpha
                    # for iyb,iyavg in zip(iyboot,yplotAvg):
                    #     print iyb,iyavg
                    plot_plane.plot(np.array(pre_extrap)+shift,iyboot,alpha=this_alpha/extrap_fade,
                                        color=color,linestyle=self.linestyle)

                yplotAvg,yplotStd = np.array(yplotAvg),np.array(yplotStd)
            if len(post_extrap) > 0:
                for iyboot in yplotBoot_post:
                    id = np.nanmean(np.abs((iyboot-yplotAvg_post)/yplotStd_post))
                    this_alpha = self.hairline_alpha/(1+id*self.hairline_dropoff)
                    # print 'DEBUG',id,this_alpha
                    # for iyb,iyavg in zip(iyboot,yplotAvg):
                    #     print iyb,iyavg
                    plot_plane.plot(np.array(post_extrap)+shift,iyboot,alpha=this_alpha/extrap_fade,
                                        color=color,linestyle=self.linestyle)
            # self.RemoveFuns()
            return plot_plane

    ## 1d plots require a dimension to plot over for x axis
    ## todo color and shift iterators?
    # def PlotFunction(self,otherXvals,xaxis=0,xdatarange='Data',thislineres=lineres):
    def PlotFunction(   self,otherXvals=[],xaxis=0,xdatarange='Data',thislineres=100,
                        color=None,shift=0.0,flowphys=False,ShowPar='None',ShowEval=None,
                        x_scale=1.0,y_scale=1.0,plot_plane=None,suppress_key=False,x_fun=None,supress_legend=False,hairline=False,extrap_fade=1.5):
        if hairline:
            return self.PlotHairlineFun(otherXvals=otherXvals,xaxis=xaxis,xdatarange=xdatarange,thislineres=thislineres,
                        color=color,shift=shift,flowphys=flowphys,ShowPar=ShowPar,ShowEval=ShowEval,
                        y_scale=y_scale,x_scale=x_scale,plot_plane=plot_plane,suppress_key=suppress_key,x_fun=x_fun,supress_legend=supress_legend,extrap_fade=extrap_fade)
        self.GetFuns()
        if y_scale != 1.0:
            self.scale_y = y_scale
        self.scale_x = x_scale
        if self.scale_x == 0:
            self.scale_x = 1
        if self.xdata.shape[0] == 1:
            thisxlist = self.xdata[0]
        else:
            if not isinstance(otherXvals,(tuple,list,np.ndarray)):
                print(otherXvals)
                raise EnvironmentError('otherXvals must be list')
            if len(otherXvals) == 0:
                print(otherXvals,' is invalid for otherXvals, defaulting to first value')
                first_Xvals = copy(np.swapaxes(self.xdata,0,1)[0])
                first_Xvals = np.delete(first_Xvals,xaxis)
                return self.PlotFunction(otherXvals=first_Xvals,xaxis=xaxis,xdatarange=xdatarange,
                                        thislineres=thislineres,x_fun=x_fun,color=color,
                                        shift=shift,flowphys=flowphys,ShowPar=ShowPar,
                                        ShowEval=ShowEval,plot_plane=plot_plane,
                                        suppress_key=suppress_key,x_scale=x_scale,
                                        y_scale=y_scale,supress_legend=supress_legend,
                                        hairline=hairline,extrap_fade=extrap_fade)
            otherXvals = list(map(float,otherXvals))
            thisxlist = []
            for icx,ix in enumerate(np.swapaxes(self.xdata,0,1)):
                thisix = copy(ix)
                thisix = np.delete(thisix,xaxis)
                if ix[xaxis] not in thisxlist:
                    thisxlist.append(ix[xaxis])
        pre_extrap,post_extrap = [],[]
        if xdatarange == 'Data':
            xdatarange = min(thisxlist),max(thisxlist)
        elif xdatarange[-1] == 'ToDataMax':
            xdatarange = float(xdatarange[0]),max(thisxlist)
            if xdatarange[0] < min(thisxlist):
                pre_extrap = np.linspace(min(thisxlist),xdatarange[0],thislineres)
        elif xdatarange[0] == 'FromDataMin':
            xdatarange = min(thisxlist),float(xdatarange[-1])
            if xdatarange[-1] > max(thisxlist):
                post_extrap = np.linspace(xdatarange[-1],max(thisxlist),thislineres)
        else:
            xdatarange = float(xdatarange[0]),float(xdatarange[-1])
            if xdatarange[-1] > max(thisxlist):
                post_extrap = np.linspace(xdatarange[-1],max(thisxlist),thislineres)
            if xdatarange[0] < min(thisxlist):
                pre_extrap = np.linspace(min(thisxlist),xdatarange[0],thislineres)
        xplot = np.linspace(min(thisxlist),max(thisxlist),thislineres)
        yplotUp,yplotDown,yplotAvg,yplotStd = [],[],[],[]
        yplotUp_pre,yplotDown_pre,yplotAvg_pre,yplotStd_pre = [],[],[],[]
        yplotUp_post,yplotDown_post,yplotAvg_post,yplotStd_post = [],[],[],[]
        for ix in xplot:
            ndxval = np.insert(otherXvals,xaxis,ix)
            evalBoot = self.Eval_Function(*ndxval)
            yplotAvg.append(evalBoot.Avg)
            yplotStd.append(evalBoot.Std)
            yplotUp.append(yplotAvg[-1] + yplotStd[-1])
            yplotDown.append(yplotAvg[-1] - yplotStd[-1])
        for ix in pre_extrap:
            ndxval = np.insert(otherXvals,xaxis,ix)
            evalBoot = self.Eval_Function(*ndxval)
            yplotAvg_pre.append(evalBoot.Avg)
            yplotStd_pre.append(evalBoot.Std)
            yplotUp_pre.append(yplotAvg_pre[-1] + yplotStd_pre[-1])
            yplotDown_pre.append(yplotAvg_pre[-1] - yplotStd_pre[-1])
        for ix in post_extrap:
            ndxval = np.insert(otherXvals,xaxis,ix)
            evalBoot = self.Eval_Function(*ndxval)
            yplotAvg_post.append(evalBoot.Avg)
            yplotStd_post.append(evalBoot.Std)
            yplotUp_post.append(yplotAvg_post[-1] + yplotStd_post[-1])
            yplotDown_post.append(yplotAvg_post[-1] - yplotStd_post[-1])
        # if scale is not None or scale != 1.0:
        #     xplot = scale*xplot
        if flowphys != False:
            xplot = TflowToPhys(xplot,flowphys)
            pre_extrap = TflowToPhys(pre_extrap,flowphys)
            post_extrap = TflowToPhys(post_extrap,flowphys)
        if suppress_key:
            thislab = ''
        else:
            thislab = self.name + ' '
        if ShowPar in self.fit_data.index:
            thisunits = ''
            if ShowPar in self.paramunits.index:
                thisunits = self.paramunits[ShowPar]
            if len(self.FPName) > 0:
                this_SP = self.FPName
            else:
                this_SP = ShowPar
            if len(self.FPUnits) > 0:
                thisunits = self.FPUnits
            else:
                if ShowPar in self.paramunits.index:
                    thisunits = self.paramunits[ShowPar]

            thislab += ''.join([this_SP,'=',MakeValAndErr(self.fit_data.Params[ShowPar].Avg,
                                                          self.fit_data.Params[ShowPar].Std,
                                                          Dec=2,latex=True),thisunits])

            # print 'REMOVE THIS AFTER'
            # scale_param = self.fit_data.Params[ShowPar]*scale
            # scale_param.Stats()
            # thislab += ' $'+' '.join([ShowPar,'=',MakeValAndErr(scale_param.Avg,scale_param.Std,Dec=3),thisunits])+'$'
        if ShowEval is not None and ShowEval is not False:
            if isinstance(ShowEval,(list,tuple,np.ndarray)):
                indep_var = ','.join(['{:.2f}'.format(ival) for ival in ShowEval])
            else:
                indep_var = '{:.2f}'.format(ShowEval)
                ShowEval = [ShowEval]
            this_eval = self.Eval_Fun_AvgStd(*ShowEval)
            if len(self.FPName) > 0:
                this_SP = self.FPName
            else:
                this_SP = 'f('+indep_var+')='
            if len(self.FPUnits) > 0:
                thisunits = self.FPUnits
            else:
                thisunits = ''
            thislab += ''.join([this_SP,MakeValAndErr(this_eval[0],this_eval[1],Dec=2,latex=True),thisunits])


        if supress_legend:
            thislab = None
        else:
            thislab = str(LegendFmt(thislab))
        if plot_plane is None:
            raise EnvironmentError('Please use plot_plane implementataion.')
            # pl.plot(np.array(xplot)+shift,yplotAvg,label=thislab,color=color,linestyle=self.linestyle)
            # pl.fill_between(np.array(xplot)+shift,yplotUp,yplotDown,alpha=fillalpha,edgecolor='none',color=color)
        else:
            if self.scale_x != 1:
                if x_fun is not None and x_fun is not False:
                    def this_xfun(val):
                        return x_fun(val)*self.scale_x
                else:
                    def this_xfun(val):
                        return val*self.scale_x
            else:
                this_xfun = x_fun
            if this_xfun is not None and this_xfun is not False:
                xplot = np.array([this_xfun(ix) for ix in xplot])
                pre_extrap = np.array([this_xfun(ix) for ix in pre_extrap])
                post_extrap = np.array([this_xfun(ix) for ix in post_extrap])
            if ShowEval is not None and ShowEval is not False:
                if this_xfun is not None and this_xfun is not False:
                    plot_plane.errorbar( this_xfun(ShowEval[xaxis]),[self.Eval_Fun_Avg(*ShowEval)],[self.Eval_Fun_Std(*ShowEval)],label=thislab,color=color,fmt='x',markersize='10')
                else:
                    plot_plane.errorbar( ShowEval[xaxis],[self.Eval_Fun_Avg(*ShowEval)],[self.Eval_Fun_Std(*ShowEval)],label=thislab,color=color,fmt='x',markersize='10')
                plot_plane.plot(    np.array(xplot)+shift,yplotAvg,color=color,linestyle=self.linestyle)
                if len(pre_extrap) > 0:
                    plot_plane.plot(    np.array(pre_extrap)+shift,yplotAvg_pre,color=color,linestyle='--')
                if len(post_extrap) > 0:
                    plot_plane.plot(    np.array(post_extrap)+shift,yplotAvg_post,color=color,linestyle='--')
            else:
                plot_plane.plot(    np.array(xplot)+shift,yplotAvg,
                                label=thislab,color=color,linestyle=self.linestyle)
                if len(pre_extrap) > 0:
                    plot_plane.plot(    np.array(pre_extrap)+shift,yplotAvg_pre,
                                    color=color,linestyle='--')
                if len(post_extrap) > 0:
                    plot_plane.plot(    np.array(post_extrap)+shift,yplotAvg_post,
                                    color=color,linestyle='--')
            plot_plane.fill_between(np.array(xplot)+shift,yplotUp,yplotDown,alpha=fillalpha,edgecolor='none',color=color)
            if len(pre_extrap) > 0:
                plot_plane.fill_between(np.array(pre_extrap)+shift,yplotUp_pre,yplotDown_pre,alpha=fillalpha/extrap_fade,edgecolor='none',color=color)
            if len(post_extrap) > 0:
                plot_plane.fill_between(np.array(post_extrap)+shift,yplotUp_post,yplotDown_post,alpha=fillalpha/extrap_fade,edgecolor='none',color=color)
            # self.RemoveFuns()
            return plot_plane

    ## 1d plots require a dimension to plot over for x axis
    ## todo color and shift iterators?
    # def PlotFunction(self,otherXvals,xaxis=0,xdatarange='Data',thislineres=lineres):
    def PlotLog(self,xdatarange='Data',thislineres=100,scale=1,color=None,shift=0.0,plot_plane=None):
        self.GetFuns()
        if xdatarange == 'Data':
            xdatarange = self.xdata[0][0],self.xdata[0][-1]
        elif xdatarange[-1] == 'ToDataMax':
            xdatarange = xdatarange[0],self.xdata[0][-1]
        elif xdatarange[0] == 'FromDataMin':
            xdatarange = self.xdata[0][0],xdatarange[-1]
        xplot = np.linspace(xdatarange[0],xdatarange[1],thislineres)
        yplotUp,yplotDown,yplotAvg,yplotStd = [],[],[],[]
        for ix in xplot:
            yplotBoot = self.Fun([ix],np.array(self.fit_data.Params.values)).Log()/scale
            yplotBoot.Stats()
            yplotAvg.append(yplotBoot.Avg)
            yplotStd.append(yplotBoot.Std)
            yplotUp.append(yplotAvg[-1] + yplotStd[-1])
            yplotDown.append(yplotAvg[-1] - yplotStd[-1])
        # self.RemoveFuns()
        if plot_plane is None:
            raise EnvironmentError('Please use plot_plane implementataion.')
            # pl.plot(np.array(xplot)+shift,yplotAvg,label=self.name,color=color,linestyle=self.linestyle)
            # pl.fill_between(np.array(xplot)+shift,yplotUp,yplotDown,alpha=fillalpha,edgecolor='none',color=color)
        else:
            plot_plane.plot(np.array(xplot)+shift,yplotAvg,label=self.name,color=color,linestyle=self.linestyle)
            plot_plane.fill_between(np.array(xplot)+shift,yplotUp,yplotDown,alpha=fillalpha,edgecolor='none',color=color)
            return plot_plane


    ## 1d plots require a dimension to plot over for x axis
    ## todo color and shift iterators?
    # def PlotFunction(self,otherXvals,xaxis=0,xdatarange='Data',thislineres=lineres):
    def PlotEffMass(self,xdatarange='Data',thislineres=100,color=None,shift=0.0,y_scale=None,x_scale=None,Phys=True,plot_plane=None):
        self.GetFuns()
        if xdatarange == 'Data':
            xdatarange = self.xdata[0][0],self.xdata[0][-1]
        elif xdatarange[-1] == 'ToDataMax':
            xdatarange = xdatarange[0],self.xdata[0][-1]
        elif xdatarange[0] == 'FromDataMin':
            xdatarange = self.xdata[0][0],xdatarange[-1]
        xplot = np.linspace(xdatarange[0],xdatarange[1],thislineres)
        # if scale is not None:
        #     xplot = scale(xplot)
        yplotUp,yplotDown,yplotAvg,yplotStd = [],[],[],[]
        for ix in xplot:
            yplotBoot = (self.Eval_Function(ix)/self.Eval_Function(ix+1)).Log()
            # yplotBoot = (self.Fun([ix],np.array(self.fit_data.Params.values))/self.Fun([ix+1],np.array(self.fit_data.Params.values))).Log()
            # if Phys:
            #     yplotBoot = yplotBoot*latparams.hbarcdivlat
            if y_scale is not None and y_scale != 0:
                yplotBoot = yplotBoot*y_scale
            yplotBoot.Stats()
            yplotAvg.append(yplotBoot.Avg)
            yplotStd.append(yplotBoot.Std)
            yplotUp.append(yplotAvg[-1] + yplotStd[-1])
            yplotDown.append(yplotAvg[-1] - yplotStd[-1])
        if y_scale is not None and y_scale != 0:
            Energyval = self.fit_data.Params['Energy']*y_scale
            Energyval.Stats()
            thislab = self.name + ' $'+''.join(['E_{p}','=',Energyval.MakeValAndErr(latex=True),'GeV'])+'$'
        else:
            thislab = self.name + ' $'+''.join(['aE_{p}','=',self.fit_data.Params['Energy'].MakeValAndErr(latex=True)])+'$'
        # self.RemoveFuns()

        if plot_plane is None:
            pass
            # pl.plot(np.array(xplot)+shift,yplotAvg,label=thislab,color=color,linestyle=self.linestyle)
            # pl.fill_between(np.array(xplot)+shift,yplotUp,yplotDown,alpha=fillalpha,edgecolor='none',color=color)
        else:
            xplot = np.array(xplot)
            if x_scale is not None and x_scale != 0:
                xplot = xplot*x_scale
            plot_plane.plot(xplot+shift,yplotAvg,label=str(LegendFmt(thislab)),color=color,linestyle=self.linestyle)
            plot_plane.fill_between(xplot+shift,yplotUp,yplotDown,alpha=fillalpha,edgecolor='none',color=color)
            return plot_plane

    def GetOutputDict(self,show_maxiter=False,show_prec=False):
        outDict = ODNested()
        if 'Params' not in self.fit_data:
            return {}
        for thispar,idata in self.fit_data['Params'].items():
            # if 'tilde' in thispar: thispar = thispar.replace('tilde',r'tilde')
            # if 'frac' in thispar: thispar = thispar.replace('frac',r'frac')
            outDict[thispar] = AvgStdToFormat(idata.Avg,idata.Std)
        if hasattr(self.fit_data,'Chi2DoF'):
            outDict['Chi^2_pdf'] = AvgStdToFormat(self.fit_data.Chi2DoF[0].Avg,self.fit_data.Chi2DoF[0].Std)
        self.GetFuns()
        if hasattr(self.Fun,'file_name'):
            outDict['Fun'] = self.Fun.file_name
        if hasattr(self.Fun,'__name__'):
            outDict['Fun'] = self.Fun.__name__
        # self.RemoveFuns()
        if isinstance(self.iGuess,(list,tuple,np.ndarray)):
            for iiGuess,(thispar,par_data) in zip(self.iGuess,iter(self.fit_data['Params'].items())):
                outDict['iGuess_'+thispar] = iiGuess
        else:
            outDict['iGuess'] = 'Default, see FitFunctions.py'
        if show_maxiter: outDict['MaxIter'] = self.MaxIter
        if show_prec: outDict['Prec'] = self.Prec
        # print 'DEBUG, getting output'
        # for ikey,ival in outDict.iteritems():
        #     print ikey,ival
        return outDict

    def GuessFromFun(self):
        if not hasattr(self.Fun,'__name__'):
            return 'None'
        if self.Fun.__name__ == 'DPfitfun':
            return [2.7,1]
            # return [-1,-1.6]
        if self.Fun.__name__ == 'DPfitfunOnePar':
            return [1.6]
            # return [-1,-1.6]
        if self.Fun.__name__ == 'A_x_exp_minus_m_x':
            return [1.0,1.0]
        if self.Fun.__name__ == 'ParmDivX':
            return [1.0]
            # return [-1,-1.6]
        if self.Fun.__name__ == 'ChitFitFun':
            return [-.01,-1.0,np.log(100)]
            # return [-1,-1.6]
        if self.Fun.__name__ == 'RandTFitFun':
            return [0.1,0.15,1.0]
            # return [-1,-1.6]
        if self.Fun.__name__ == 'Chipt_Chit':
            return [1.0]
            # return [-1,-1.6]
        if self.Fun.__name__ == 'ParmDivXP':
            return [-1.0,1.5]
            # return [-1,-1.6]
        if self.Fun.__name__ == 'DPfitfun2':
            return [2.7,0.6]
        if self.Fun.__name__ == 'FormFactorO1':
            return [1]
        if self.Fun.__name__ == 'FormFactorO2':
            return [1,1]
        if self.Fun.__name__ == 'FormFactorO3':
            # return [1,1,1]
            return [1.0,0.5,2.0]
        if self.Fun.__name__ == 'FormFactorO':
            # return [1]*Len
            return 'None'
        if self.Fun.__name__ == 'ConstantFitFun':
            return [1]
        elif self.Fun.__name__ == 'c3ContFitFun':
            return [1,1,1]
        elif self.Fun.__name__ == 'SchiffFitFun':
            return [1,1,1]
        elif self.Fun.__name__ == 'LinearFitFun':
            return [1,1]
        elif self.Fun.__name__ == 'LinearFitFun_x0':
            return [1]
        elif self.Fun.__name__ == 'LogFitFun':
            return [1,1]
        elif self.Fun.__name__ == 'LogFitFun_x0':
            return [1]
        elif self.Fun.__name__ == 'SquaredFitFun':
            return [1,1]
        elif self.Fun.__name__ == 'SquaredFitFun_x0':
            return [1]
        elif self.Fun.__name__ == 'SqPlLogFitFun':
            return [1,1,1]
        elif self.Fun.__name__ == 'SqPlLogFitFun_x0':
            return [1,1]
        elif self.Fun.__name__ == 'ContFitFun':
            return [1,1,1]
        elif self.Fun.__name__ == 'ContFitFun_x0':
            return [1,1]
        elif self.Fun.__name__ == 'ContPlLogFitFun':
            return [1,1,1,1]
        elif self.Fun.__name__ == 'ContPlLogFitFun_x0':
            return [1,1,1]
        elif self.Fun.__name__ == 'ContFitFun_mix':
            return [1,1,1,1]
        elif self.Fun.__name__ == 'ContFitFun_x0_mix':
            return [1,1,1]
        elif self.Fun.__name__ == 'ContPlLogFitFun_mix':
            return [1,1,1,1,1]
        elif self.Fun.__name__ == 'ContPlLogFitFun_x0_mix':
            return [1,1,1,1]
        elif self.Fun.__name__ == 'LinearFF_2D':
            return [1,1,1]
        elif self.Fun.__name__ == 'C2TwoStateFitFun':
            # return [3.1861104305e-15,100,np.log(.6),np.log(0.3385347963)]
            return [3.1861104305e-7,10,np.log(.5),np.log(0.2)]
        elif self.Fun.__name__ == 'C2TwoStateFitFunNoExp':
            # return [3.1861104305e-15,100,.6,0.3385347963]
            return [3.1861104305e-7,10,.5,0.2]
        elif self.Fun.__name__ == 'C2OSFAntiper':
            return [1,0.3]
        elif 'C2atT' in self.Fun.__name__:
            return [1,0.3]
        elif self.Fun.__name__ == 'C2OneStateFitFun':
            return [1.3422466805e-10,np.log(.6)]
            # return [1.3422466805e-3,np.log(.3)]
        elif self.Fun.__name__ == 'C2OneStateFitFunNoExp':
            # return [1.3422466805e-3,.3]
            return [2.7174048910e-07,.6]
        elif self.Fun.__name__ == 'C3OneStateFitFun':
            return [0]
        elif self.Fun.__name__ == 'C2TwoStateFitFunCM':
            return 'None'
            # output = []
            # for i in xrange(Len):
            #     output.append(3.1861104305e-06)
            #     output.append(4)
            # return output+[np.log(.45),np.log(0.3385347963)]
        elif self.Fun.__name__ == 'C3MomTwoStateFitFun':
            return [0,0,0,0]
        elif self.Fun.__name__ == 'C3TwoStateFitFun2par':
            return [0,0]
        elif self.Fun.__name__ == 'C3TwoStateFitFun3par':
            return [0,0,0]
        elif self.Fun.__name__ == 'TestTwoVarFitFun':
            return [1,1]
        elif self.Fun.__name__ == 'OneOnRootNFitFun':
            return [1]
        elif self.Fun.__name__ == 'c3FitFun_NoA':
            return [np.exp(2.8),np.exp(-0.11)]
        elif self.Fun.__name__ == 'c3PolyFun':
            return [-1]
        elif self.Fun.__name__ == 'c3PolyFun2':
            return [-1,-1]
        elif self.Fun.__name__ == 'c3PolyFun3':
            return [-1,-1,-1]
        elif self.Fun.__name__ == 'c3PolyFun_skip1':
            return [-1]
        elif self.Fun.__name__ == 'c3PolyFun2_skip1':
            return [-1,-1]
        elif self.Fun.__name__ == 'c3PolyFun3_skip1':
            return [-1,-1,-1]
        elif self.Fun.__name__ == 'c3FitFun':
            return [1,1,1]
        # elif self.Fun.__name__ == 'c3FitFun_V2':
        #     return [2.8,-0.11,-1]
        elif self.Fun.__name__ == 'c3ExpFitFun':
            return [-5]
        elif self.Fun.__name__ == 'c3ExpFitFun2':
            return [-5,-10]
        elif self.Fun.__name__ == 'c3ExpFitFun_nosqrt':
            return [-5]
        elif self.Fun.__name__ == 'c3ExpFitFun2_nosqrt':
            return [-5,-10]
        elif self.Fun.__name__ == 'c3ExpFitFun3_nosqrt':
            return [-5,-10,-15]
        elif self.Fun.__name__ == 'Chi_Prop_Fit':
            return [0.20,1,1.]
        elif self.Fun.__name__ == 'Chi_Prop_Fit_2exp':
            return [0.20,1.,0.75,-1.05,0.8]
            # return [0.15,0,0,0,0]
        elif self.Fun.__name__ == 'Alpha_Prop_Fit':
            return [-0.1,4,1]
        elif self.Fun.__name__ == 'Alpha_Prop_Fit_plin':
            return [-0.1,4,1,0]
        elif self.Fun.__name__ == 'Alpha_Prop_Fit_2exp':
            return [-0.1,-4,7,0,1]
        else:
            return 'None'



    def DerOfFun(self):
        if not hasattr(self.Fun,'__name__'):
            return 'Estimate'
        if self.Fun.__name__ == 'ConstantFitFun':
            return ConstFFDer
        elif self.Fun.__name__ == 'c3ContFitFun':
            return c3ContFFDer
        elif self.Fun.__name__ == 'SchiffFitFun':
            return SchiffFFDer
        elif self.Fun.__name__ == 'LinearFitFun':
            return LinFFDer
        elif self.Fun.__name__ == 'LinearFitFun_x0':
            return LinFFDer_x0
        elif self.Fun.__name__ == 'LogFitFun':
            return LogFFDer
        elif self.Fun.__name__ == 'LogFitFun_x0':
            return LogFFDer_x0
        elif self.Fun.__name__ == 'SquaredFitFun':
            return SquaredFFDer
        elif self.Fun.__name__ == 'SquaredFitFun_x0':
            return SquaredFFDer_x0
        elif self.Fun.__name__ == 'SqPlLogFitFun':
            return SqPlLogFFDer
        elif self.Fun.__name__ == 'SqPlLogFitFun_x0':
            return SqPlLogFFDer_x0
        elif self.Fun.__name__ == 'ContFitFun':
            return ContFFDer
        elif self.Fun.__name__ == 'ContFitFun_x0':
            return ContFFDer_x0
        elif self.Fun.__name__ == 'ContPlLogFitFun':
            return ContPlLogFFDer
        elif self.Fun.__name__ == 'ContPlLogFitFun_x0':
            return ContPlLogFFDer_x0
        elif self.Fun.__name__ == 'ContFitFun_mix':
            return ContFFDer_mix
        elif self.Fun.__name__ == 'ContFitFun_x0_mix':
            return ContFFDer_x0_mix
        elif self.Fun.__name__ == 'ContPlLogFitFun_mix':
            return ContPlLogFFDer_mix
        elif self.Fun.__name__ == 'ContPlLogFitFun_x0_mix':
            return ContPlLogFFDer_x0_mix
        elif self.Fun.__name__ == 'LinearFF_2D':
            return LinFFDer_2D
        elif self.Fun.__name__ == 'C2OSFAntiper':
            return C2OSFAntiperDer
        elif self.Fun.__name__ == 'C2OneStateFitFun':
            return C2OSFFDer
        elif self.Fun.__name__ == 'C2OneStateFitFunNoExp':
            return C2OSFFNoExpDer
        elif self.Fun.__name__ == 'C3OneStateFitFun':
            return C3OSFFDer
        elif self.Fun.__name__ == 'C2TwoStateFitFun':
            return C2TSFFDer
        elif self.Fun.__name__ == 'C2TwoStateFitFunNoExp':
            return C2TSFFNoExpDer
        # elif self.Fun.__name__ == 'C2TwoStateFitFunCM':
        #     return C2TSFFCMDer
        elif self.Fun.__name__ == 'C3MomTwoStateFitFun':
            return C3MomTSFFDer
        elif self.Fun.__name__ == 'TestTwoVarFitFun':
            return TestTwoVarFFDer
        elif self.Fun.__name__ == 'FormFactorO1':
            return FormFactorO1Der
        elif self.Fun.__name__ == 'FormFactorO2':
            return FormFactorO2Der
        elif self.Fun.__name__ == 'FormFactorO3':
            return FormFactorO3Der
        elif self.Fun.__name__ == 'DPfitfunOnePar':
            return DPfitfunOneParDer
        elif self.Fun.__name__ == 'A_x_exp_minus_m_x':
            return 'Estimate'
        elif self.Fun.__name__ == 'ParmDivX':
            return ParmDivXDer
        elif self.Fun.__name__ == 'ChitFitFun':
            return ChitFitFunDer
        elif self.Fun.__name__ == 'RandTFitFun':
            return 'Estimate'
        elif self.Fun.__name__ == 'Chipt_Chit':
            return 'Estimate'
        elif self.Fun.__name__ == 'ParmDivXP':
            return ParmDivXPDer
        elif self.Fun.__name__ == 'DPfitfun2':
            return DPfitfun2Der
        elif self.Fun.__name__ == 'DPfitfun':
            return DPfitfunDer
        elif self.Fun.__name__ == 'OneOnRootNFitFun':
            return OORNFFDer
        elif self.Fun.__name__ == 'c3FitFun':
            # return 'Estimate'
            return c3FFDer
        elif self.Fun.__name__ == 'c3FitFun_NoA':
            return c3FFDer_NoA
        elif self.Fun.__name__ == 'c3PolyFun':
            return c3PFDer
        elif self.Fun.__name__ == 'c3PolyFun2':
            return c3PFDer2
        elif self.Fun.__name__ == 'c3PolyFun3':
            return c3PFDer3
        elif self.Fun.__name__ == 'c3PolyFun_skip1':
            return c3PFDer_skip1
        elif self.Fun.__name__ == 'c3PolyFun2_skip1':
            return c3PFDer2_skip1
        elif self.Fun.__name__ == 'c3PolyFun3_skip1':
            return c3PFDer3_skip1
        elif self.Fun.__name__ == 'c3ExpFitFun':
            return c3EFFDer
        elif self.Fun.__name__ == 'c3ExpFitFun2':
            return c3EFFDer2
        elif self.Fun.__name__ == 'c3ExpFitFun_nosqrt':
            return c3EFFDer_nosqrt
        elif self.Fun.__name__ == 'c3ExpFitFun2_nosqrt':
            return c3EFFDer2_nosqrt
        elif self.Fun.__name__ == 'c3ExpFitFun3_nosqrt':
            return c3EFFDer3_nosqrt
        elif self.Fun.__name__ == 'Chi_Prop_Fit_2exp':
            return Chi_Prop_Der
        elif self.Fun.__name__ == 'Chi_Prop_Fit_2exp':
            return Chi_Prop_Der_2exp
        elif self.Fun.__name__ == 'Alpha_Prop_Fit':
            return Alpha_Prop_Der
        elif self.Fun.__name__ == 'Alpha_Prop_Fit_plin':
            return Alpha_Prop_Der_plin
        elif self.Fun.__name__ == 'Alpha_Prop_Fit_2exp':
            return Alpha_Prop_Der_2exp
        elif self.Fun.__name__ == 'RF_Prop_Fit':
            return RF_Prop_Der
        else:
            return 'Estimate'


    def ParamLabFromFun(self):
        if not hasattr(self.Fun,'__name__'):
            return 'Def'
        if self.Fun.__name__ == 'ConstantFitFun':
            return 'Def'
        elif self.Fun.__name__ == 'SchiffFitFun':
            return 'Def'
        elif self.Fun.__name__ == 'c3ContFitFun':
            return 'c','mpi2_coeff','a_coeff'
        elif self.Fun.__name__ == 'LinearFitFun':
            return 'Def'
        elif self.Fun.__name__ == 'LinearFitFun_x0':
            return 'Def'
        elif self.Fun.__name__ == 'LogFitFun':
            return 'Def'
        elif self.Fun.__name__ == 'LogFitFun_x0':
            return 'Def'
        elif self.Fun.__name__ == 'SquaredFitFun':
            return 'Def'
        elif self.Fun.__name__ == 'SquaredFitFun_x0':
            return 'Def'
        elif self.Fun.__name__ == 'SqPlLogFitFun':
            return 'Def'
        elif self.Fun.__name__ == 'SqPlLogFitFun_x0':
            return 'Def'
        elif self.Fun.__name__ == 'ContFitFun':
            return 'Def'
        elif self.Fun.__name__ == 'ContFitFun_x0':
            return 'Def'
        elif self.Fun.__name__ == 'ContPlLogFitFun':
            return 'Def'
        elif self.Fun.__name__ == 'ContPlLogFitFun_x0':
            return 'Def'
        elif self.Fun.__name__ == 'ContFitFun_mix':
            return 'Def'
        elif self.Fun.__name__ == 'ContFitFun_x0_mix':
            return 'Def'
        elif self.Fun.__name__ == 'ContPlLogFitFun_mix':
            return 'Def'
        elif self.Fun.__name__ == 'ContPlLogFitFun_x0_mix':
            return 'Def'
        elif self.Fun.__name__ == 'LinearFF_2D':
            return 'Def'
        elif 'C2atT' in self.Fun.__name__:
            return 'A_Ep','Energy'
        elif self.Fun.__name__ == 'C2OSFAntiper':
            return 'A_Ep','Energy'
        elif self.Fun.__name__ == 'C2OneStateFitFun':
            return 'A_Ep','Energy'
        elif self.Fun.__name__ == 'C2OneStateFitFunNoExp':
            return 'A_Ep','Energy'
        elif self.Fun.__name__ == 'C3OneStateFitFun':
            return 'B_00',
        elif self.Fun.__name__ == 'C2TwoStateFitFun':
            return 'A_Ep','A_Epp','Energy','DE'
        elif self.Fun.__name__ == 'C2TwoStateFitFunNoExp':
            return 'A_Ep','A_Epp','Energy','DE'
        ## Not Correct, dont use
        elif self.Fun.__name__ == 'C2TwoStateFitFunCM':
            return 'Def'
            # return 'A_Ep','A_Epp','Energy','DE'
        elif self.Fun.__name__ == 'C3MomTwoStateFitFun':
            return 'B_00','B_10','B_01','B_11'
        elif self.Fun.__name__ == 'TestTwoVarFitFun':
            return 'Def'
        elif self.Fun.__name__ == 'FormFactorO1':
            return 'FF1',
        elif self.Fun.__name__ == 'FormFactorO2':
            return 'FF1','FF2'
        elif self.Fun.__name__ == 'FormFactorO3':
            return 'FF1','FF2','FF3'
        elif self.Fun.__name__ == 'DPfitfunOnePar':
            return ['F(0)']
        elif self.Fun.__name__ == 'A_x_exp_minus_m_x':
            return ['A','m']
        elif self.Fun.__name__ == 'ParmDivX':
            return 'Def'
        elif self.Fun.__name__ == 'ChitFitFun':
            return 'Def'
        elif self.Fun.__name__ == 'RandTFitFun':
            return 'y_0','y_1','m'
        elif self.Fun.__name__ == 'Chipt_Chit':
            return ['Sig']
        elif self.Fun.__name__ == 'ParmDivXP':
            return 'Def'
        elif self.Fun.__name__ == 'DPfitfun2':
            return 'Def'
        elif self.Fun.__name__ == 'DPfitfun':
            return 'F(0)','m_{EM}'
        elif self.Fun.__name__ == 'OneOnRootNFitFun':
            return 'Def'
        elif self.Fun.__name__ == 'c3PolyFun':
            return ['b']
        elif self.Fun.__name__ == 'c3PolyFun2':
            return ['b','c']
        elif self.Fun.__name__ == 'c3PolyFun3':
            return ['b','c','d']
        elif self.Fun.__name__ == 'c3PolyFun_skip1':
            return ['b']
        elif self.Fun.__name__ == 'c3PolyFun2_skip1':
            return ['b','c']
        elif self.Fun.__name__ == 'c3PolyFun3_skip1':
            return ['b','c','d']
        elif self.Fun.__name__ == 'c3FitFun_NoA':
            return ['b','c']
        elif self.Fun.__name__ == 'c3FitFun':
            return ['a','b','c']
        # elif self.Fun.__name__ == 'c3FitFun_V2':
        #     return ['b','c','d']
        elif self.Fun.__name__ == 'c3ExpFitFun':
            return ['a']
        elif self.Fun.__name__ == 'c3ExpFitFun2':
            return ['a','b']
        elif self.Fun.__name__ == 'c3ExpFitFun_nosqrt':
            return ['a']
        elif self.Fun.__name__ == 'c3ExpFitFun2_nosqrt':
            return ['a','b']
        elif self.Fun.__name__ == 'c3ExpFitFun3_nosqrt':
            return ['a','b','c']
        elif self.Fun.__name__ == 'c3FitFun_mine':
            return ['a','b','c','d']
        elif self.Fun.__name__ == 'Chi_Prop_Fit':
            return ['A','B','E']
        elif self.Fun.__name__ == 'Chi_Prop_Fit_2exp':
            return ['A','B_1','E_1','B_2','E_2']
        elif self.Fun.__name__ == 'Alpha_Prop_Fit':
            return [r'\alpha','B','E']
        elif self.Fun.__name__ == 'Alpha_Prop_Fit_plin':
            return [r'\alpha','B','E','C']
        elif self.Fun.__name__ == 'Alpha_Prop_Fit_2exp':
            return [r'\alpha','B_1','E_1','B_2','E_2']
        elif self.Fun.__name__ == 'RF_Prop_Fit':
            return 'Def'
        else:
            return 'Def'

    def Overload_wrap(self,FF2,this_fun):
        ## TODO, include some name and filename formatting
        thisname = self.name
        if hasattr(FF2,'name'):
            thisname += '_'+this_fun.__name__+'_'+FF2.name



        result = Fitting( thisPrec=self.Prec,thisMI=self.MaxIter,iGuess =self.iGuess,name=thisname)
        result.fit_data = pa.DataFrame()
        result.fit_data.loc[:,'Params_units'] = self.fit_data.loc[:,'Params_units']
        self.GetFuns()
        result.ImportFunction(self.Fun,self.parlen)
        result.xdata = self.xdata
        # result.UpdateName('+',self,FF2)
        if len(self.fit_data) != 0:
            if isinstance(FF2,Fitting):
                # result.fit_data.loc[:,'boot'] = this_fun(self.fit_data.loc[:,'boot'],FF2.fit_data.loc[:,'boot'])
                if len(FF2.fit_data) != 0:
                    try:
                        result.fit_data.loc[:,'Params'] = this_fun(self.fit_data.loc[:,'Params'],FF2.fit_data.loc[:,'Params'])
                        ## TODO, make this proper
                        FF2.GetFuns()
                        if self.Fun.__name__ != FF2.Fun.__name__:
                            print('Warning, fit classes being combine have different fit functions')
                            print(self.Fun.__name__)
                            print(FF2.Fun.__name__)
                        if self.parlen < FF2.parlen:
                            result.ImportFunction(FF2.Fun,FF2.parlen)
                        result.xdata = []
                        for ix,ix2 in zip(self.xdata,FF2.xdata):
                            result.xdata.append(sorted(set(ix).union(set(ix2))))
                        result.xdata = np.array(result.xdata)
                    except Exception as err:
                        print(self.fit_data.index)
                        print(FF2.fit_data.index)
                        raise EnvironmentError(str(err)+'\nError trying to combine fit_data, most likely different parameters\n')
                else:
                    raise EnvironmentError('Run fit on second class before combining fits')
            else:
                try:
                    result.fit_data.loc[:,'Params'] = this_fun(self.fit_data.loc[:,'Params'], FF2)
                except Exception as err:
                    print(type(FF2), FF2)
                    print(self.fit_data)
                    raise EnvironmentError(str(err)+'\nInvalid value to combine with fit_data class')
            return result
        else:
            raise EnvironmentError('Run fit on first class before combining fits')

    def __add__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in Fitting.parentlist]):
            return FF2.__radd__(self)
        else:
            return self.Overload_wrap(FF2,op.add)

    def __sub__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in Fitting.parentlist]):
            return FF2.__rsub__(self)
        else:
            return self.Overload_wrap(FF2,op.sub)

    def __mul__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in Fitting.parentlist]):
            return FF2.__rmul__(self)
        else:
            return self.Overload_wrap(FF2,op.mul)

    def __div__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in Fitting.parentlist]):
            return FF2.__rdiv__(self)
        else:
            return self.Overload_wrap(FF2,op.truediv)

    def __pow__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in Fitting.parentlist]):
            return FF2.__rpow__(self)
        else:
            return self.Overload_wrap(FF2,op.pow)


    ## Right multiplication functions


    def __radd__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in Fitting.parentlist]):
            return FF2.__add__(self)
        else:
            return self.Overload_wrap(FF2,flip_fun_2arg(op.add))

    def __rsub__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in Fitting.parentlist]):
            return FF2.__sub__(self)
        else:
            return self.Overload_wrap(FF2,flip_fun_2arg(op.sub))

    def __rmul__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in Fitting.parentlist]):
            return FF2.__mul__(self)
        else:
            return self.Overload_wrap(FF2,flip_fun_2arg(op.mul))

    def __rdiv__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in Fitting.parentlist]):
            return FF2.__div__(self)
        else:
            return self.Overload_wrap(FF2,flip_fun_2arg(op.truediv))


    def __rpow__(self,FF2):
        if any([ipar in str(type(FF2)) for ipar in Fitting.parentlist]):
            return FF2.__pow__(self)
        else:
            return self.Overload_wrap(FF2,flip_fun_2arg(op.pow))

def TestComb():
    low=-5
    high=5
    npoints=10
    nboot=200
    values = np.random.uniform(low=low,high=high,size=(npoints,nboot))
    values = np.array([ival + icv for icv,ival in enumerate(values)])
    testdata = [BootStrap(nboot,name='test_bootstrap',cfgvals=ival,thisDelVal=False) for ival in values]

    testdf = pa.DataFrame()
    testdf.loc[:,'xdata'] = pa.Series([[np.arange(npoints)]])
    testdf.loc[:,'ydata'] = pa.Series([testdata])
    testdf.loc[:,'ydata_Avg'] = pa.Series([idata.Avg for idata in testdata])
    testdf.loc[:,'ydata_Std'] = pa.Series([idata.Std for idata in testdata])

    testdf_plot = pa.DataFrame()
    testdf_plot.loc[:,'xdata'] = pa.Series(np.arange(npoints))
    testdf_plot.loc[:,'ydata'] = pa.Series(testdata)
    testdf_plot.loc[:,'ydata_Avg'] = pa.Series([idata.Avg for idata in testdata])
    testdf_plot.loc[:,'ydata_Std'] = pa.Series([idata.Std for idata in testdata])

    testfit = Fitting(data = testdf, Funs=[LinearFitFun,2],name='Test')

    testfit.FitBoots()

    low=-10
    high=20
    npoints=10
    nboot=200
    values2 = np.random.uniform(low=low,high=high,size=(npoints,nboot))
    values2 = np.array([ival + icv for icv,ival in enumerate(values)])
    testdata2 = [BootStrap(nboot,name='test_bootstrap',cfgvals=ival,thisDelVal=False) for ival in values2]


    testdf2 = pa.DataFrame()
    testdf2.loc[:,'xdata'] = pa.Series([[np.arange(npoints)]])
    testdf2.loc[:,'ydata'] = pa.Series([testdata2])
    testdf2.loc[:,'ydata_Avg'] = pa.Series([idata.Avg for idata in testdata2])
    testdf2.loc[:,'ydata_Std'] = pa.Series([idata.Std for idata in testdata2])

    testdf2_plot = pa.DataFrame()
    testdf2_plot.loc[:,'xdata'] = pa.Series(np.arange(npoints))
    testdf2_plot.loc[:,'ydata'] = pa.Series(testdata2)
    testdf2_plot.loc[:,'ydata_Avg'] = pa.Series([idata.Avg for idata in testdata2])
    testdf2_plot.loc[:,'ydata_Std'] = pa.Series([idata.Std for idata in testdata2])

    testfit2 = Fitting(data = testdf2, Funs=[LinearFitFun,2],name='Test2')

    testfit2.FitBoots()

    testcomb = np.array(testdata) + np.array(testdata2)
    for itest in testcomb:
        itest.Stats()

    testdf_comb_plot = pa.DataFrame()
    testdf_comb_plot.loc[:,'xdata'] = pa.Series(np.arange(npoints))
    testdf_comb_plot.loc[:,'ydata'] = pa.Series(testcomb)
    testdf_comb_plot.loc[:,'ydata_Avg'] = pa.Series([idata.Avg for idata in testcomb])
    testdf_comb_plot.loc[:,'ydata_Std'] = pa.Series([idata.Std for idata in testcomb])


    testfit_comb = testfit + testfit2

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
    hold_series['color'] = 'Blue'
    this_data['test_data'] = hold_series


    hold_series = pa.Series()
    hold_series['type'] = 'fit'
    hold_series['fit_class'] = testfit
    hold_series['color'] = 'previous'
    this_data['test_fit'] = hold_series


    hold_series = pa.Series()
    hold_series['x_data'] = testdf2_plot['xdata']
    hold_series['y_data'] = testdf2_plot['ydata_Avg']
    hold_series['yerr_data'] = testdf2_plot['ydata_Std']
    hold_series['fit_class'] = None
    hold_series['type'] = 'error_bar'
    hold_series['label'] = 'test data 2'
    hold_series['color'] = 'Red'
    this_data['test_data2'] = hold_series


    hold_series = pa.Series()
    hold_series['type'] = 'fit'
    hold_series['fit_class'] = testfit2
    hold_series['color'] = 'previous'
    this_data['test_fit2'] = hold_series


    hold_series = pa.Series()
    hold_series['x_data'] = testdf_comb_plot['xdata']
    hold_series['y_data'] = testdf_comb_plot['ydata_Avg']
    hold_series['yerr_data'] = testdf_comb_plot['ydata_Std']
    hold_series['fit_class'] = None
    hold_series['type'] = 'error_bar'
    hold_series['color'] = 'Green'
    this_data['test_data_comb'] = hold_series


    hold_series = pa.Series()
    hold_series['type'] = 'fit'
    hold_series['fit_class'] = testfit_comb
    hold_series['color'] = 'previous'
    this_data['test_fit_comb'] = hold_series


    this_info = pa.Series()
    this_info['save_file'] = './TestGraphs/TestFitPlot.pdf'
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

    return testfit,testfit2,testfit_comb

def TestFit(low=-5,high=5,npoints=10,nboot=200):
    values = np.random.uniform(low=low,high=high,size=(npoints,nboot))
    values = np.array([ival + icv for icv,ival in enumerate(values)])
    testdata = [BootStrap(nboot,name='test_bootstrap',cfgvals=ival,thisDelVal=False) for ival in values]

    testdf = pa.DataFrame()
    testdf.loc[:,'xdata'] = pa.Series([[np.arange(npoints)]])
    testdf.loc[:,'ydata'] = pa.Series([testdata])
    testdf.loc[:,'ydata_Avg'] = pa.Series([idata.Avg for idata in testdata])
    testdf.loc[:,'ydata_Std'] = pa.Series([idata.Std for idata in testdata])

    testdf_plot = pa.DataFrame()
    testdf_plot.loc[:,'xdata'] = pa.Series(np.arange(npoints))
    testdf_plot.loc[:,'ydata'] = pa.Series(testdata)
    testdf_plot.loc[:,'ydata_Avg'] = pa.Series([idata.Avg for idata in testdata])
    testdf_plot.loc[:,'ydata_Std'] = pa.Series([idata.Std for idata in testdata])

    testfit = Fitting(data = testdf, Funs=[LinearFitFun,2],name='Test')

    testfit.FitBoots()

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
    hold_series['color'] = 'Blue'
    this_data['test_data'] = hold_series


    hold_series = pa.Series()
    hold_series['type'] = 'fit'
    hold_series['fit_class'] = testfit
    hold_series['label'] = 'test fit'
    hold_series['color'] = 'previous'
    this_data['test_fit'] = hold_series


    this_info = pa.Series()
    this_info['save_file'] = './TestGraphs/TestFitPlot.pdf'
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


    # testfit.PlotFunction()
    # pl.legend()
    # pl.savefig('./TestGraphs/TestFuns.pdf')

    #
    return testfit

def WriteAllFuns():
    import types
    import PredefFitFuns as fitfuns
    functions = [fitfuns.__dict__.get(a) for a in dir(fitfuns) if isinstance(fitfuns.__dict__.get(a), types.FunctionType)]
    WriteFuns(*functions)

if __name__ == '__main__':
    print('Writing all Functions that have been predefined')
    WriteAllFuns()
    # data = TestFit()
##Class with Fitting to group all Functions.
