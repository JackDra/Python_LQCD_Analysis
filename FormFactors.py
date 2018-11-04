#!/usr/bin/env python

from GammaMatricies import CorrSpinTrace
from MiscFuns import ODNested,mkdir_p
import MomParams as mp
import numpy as np
import pandas as pa

from BootStrapping import BootStrap
from PredefFitFuns import FormFactorO1,FormFactorO2,FormFactorO3
from TwoPtCorrelators import TwoPointCorr,NNQCorr,NNQFullCorr
from FitFunctions import Fitting
from Params import outputdir,defInfo,alpha_guess,myeps
from FileIO import WriteExcel
from XmlFormatting import AvgStdToFormat
from FileIO import WriteXml,WritePickle,ReadPickleWrap,ReadFuns,WriteFuns
from PredefFitFuns import DPfitfun,LinearFitFun,A_x_exp_minus_m_x
# import dill as pickle
import os
from collections import OrderedDict
from TimeStuff import Timer
from copy import deepcopy


## form factor parameters defined up here ########################################################################

veclist = [1,2,3,4]

gi = [('g'+str(i),) for i in veclist]
giTop = [('g'+str(i),'Top') for i in veclist]
# giTop = ['g'+str(i)+'Top' for i in [4]]
giWein = [('g'+str(i),'Wein') for i in veclist]
gig5 = [('g'+str(i),'g5') for i in veclist]
g5gi = [('g5','g'+str(i)) for i in veclist]
Pi = [('P'+str(i),) for i in [4,3]]

gigj = []
for i,gii in enumerate(gi):
    for j,gij in enumerate(gi):
        if j>i:
            gigj.append((gii[0],gij[0]))
gjgi = []
for i,gii in enumerate(gi):
    for j,gij in enumerate(gi):
        if i>j:
            gjgi.append((gii[0],gij[0]))

# gic = [('g'+str(i),'cmplx') for i in veclist]
# giTopc = [('g'+str(i),'Top','cmplx') for i in veclist]
# # giTop = ['g'+str(i)+'Top' for i in [4]]
# giWeinc = [('g'+str(i),'Wein','cmplx') for i in veclist]
# gig5c = [('g'+str(i),'g5','cmplx') for i in veclist]
# g5gic = [('g5','g'+str(i),'cmplx') for i in veclist]

# gigjc = []
# for i,gii in enumerate(gi):
#     for j,gij in enumerate(gi):
#         if j>i:
#             gigj.append((gii[0],gij[0],'cmplx'))
# gjgic = []
# for i,gii in enumerate(gi):
#     for j,gij in enumerate(gi):
#         if i>j:
#             gjgi.append((gii[0],gij[0],'cmplx'))



PiI    = []
Pig5   = []
Pigi   = []
Pigig5 = []
Pigigj = []
PigiTop   = []
PigiWein   = []
# PiI,PiIc    = [],[]
# Pig5,Pig5c   = [],[]
# Pigi,Pigic   = [],[]
# Pigig5,Pigig5c = [],[]
# Pigigj,Pigigjc = [],[]
# PigiTop,PigiTopc   = [],[]
# PigiWein,PigiWeinc   = [],[]
for ip in Pi:
    PiI    += list(zip(ip,['I']))
    Pig5   += list(zip(ip,['g5']))
    Pigi   += [ip+igi for igi in gi]
    Pigig5 += [ip+igig5 for igig5 in gig5]
    Pigigj += [ip+igigj for igigj in gigj]
    PigiTop   += [ip+igit for igit in giTop]
    PigiWein   += [ip+igiw for igiw in giWein]
    # PiIc    += ((ip + ('I','cmplx')),)
    # Pig5c   += ((ip + ('g5','cmplx')),)
    # Pigic   += [ip+igi for igi in gic]
    # Pigig5c += [ip+igig5 for igig5 in gig5c]
    # Pigigjc += [ip+igigj for igigj in gigjc]
    # PigiTopc   += [ip+igit for igit in giTopc]
    # PigiWeinc   += [ip+igiw for igiw in giWeinc]

# CurrTypes = ['Scalar','Vector','VectorTop','VectorWein','VectorTopNoV','PsScalar','PsVector','Tensor']
##TODO implement Tensor
CurrTypes = ['Scalar','Vector','VectorTop','VectorTopBNL','VectorTopBNL_mal','VectorWein','VectorWeinBNL','PsScalar','PsVector','Tensor','GeGm']
DerCurrTypes = ['giDi']
AllCurrTypes = CurrTypes + DerCurrTypes

# CurrOpps = {'Scalar'   : PiI,
#             'Vector'   : Pigi+Pigic,
#             'VectorTop': Pigi+PigiTop+Pigic+PigiTopc,
#             'VectorWein': Pigi+PigiWein+Pigic+PigiWeinc,
#             'VectorTopNoV': PigiTop+PigiTopc,
#             'PsScalar' : Pig5+Pig5c,
#             'PsVector' : Pigig5+Pigig5c}
#             # 'Tensor'   : Pigigj+Pigigjc}


# CurrOppsNoProj = {'Scalar'   : [('I',)],
#                   'Vector'   : gi+gic,
#                   'VectorTop': gi+giTop+gic+giTopc,
#                   'VectorWein': gi+giWein+gic+giWeinc,
#                   'VectorTopNoV': giTop+giTopc,
#                   'PsScalar' : [('g5',),('g5','cmplx')],
#                   'PsVector' : gig5+gig5c}
#                   # 'Tensor'   : gigj+gigjc}

# ## changes order of gamma matricies for  for gig5 and gigj
# CurrOppsNoProjSigBack = {'Scalar'   : [('I',)],
#                          'Vector'   : gi+gic,
#                          'VectorTop': gi+giTop+gic+giTopc,
#                          'VectorWein': gi+giWein+gic+giWeinc,
#                          'VectorTopNoV': giTop+giTopc,
#                          'PsScalar' : [('g5',),('g5','cmplx')],
#                          'PsVector' : g5gi+g5gic}
#                          # 'Tensor'   : gjgi+gjgic}

CurrOpps = {'Scalar'   : PiI,
            'Vector'   : Pigi,
            'GeGm'   : Pigi,
            'VectorTop': Pigi+PigiTop,
            'VectorTopBNL': Pigi+PigiTop,
            'VectorTopBNL_mal': Pigi+PigiTop,
            'VectorWein': Pigi+PigiWein,
            'VectorWeinBNL': Pigi+PigiWein,
            'VectorTopNoV': PigiTop,
            'PsScalar' : Pig5,
            'PsVector' : Pigig5,
            'Tensor'   : Pigigj}


CurrOppsNoProj = {'Scalar'   : [('I',)],
                  'Vector'   : gi,
                  'GeGm'   : gi,
                  'VectorTop': gi+giTop,
                  'VectorTopBNL': gi+giTop,
                  'VectorTopBNL_mal': gi+giTop,
                  'VectorWein': gi+giWein,
                  'VectorWeinBNL': gi+giWein,
                  'VectorTopNoV': giTop,
                  'PsScalar' : [('g5',)],
                  'PsVector' : gig5,
                  'Tensor'   : gigj}

## changes order of gamma matricies for  for gig5 and gigj
CurrOppsNoProjSigBack = {'Scalar'   : [('I',)],
                         'Vector'   : gi,
                         'GeGm'   : gi,
                         'VectorTop': gi+giTop,
                         'VectorTopBNL': gi+giTop,
                         'VectorTopBNL_mal': gi+giTop,
                         'VectorWein': gi+giWein,
                         'VectorWeinBNL': gi+giWein,
                         'VectorTopNoV': giTop,
                         'PsScalar' : [('g5',)],
                         'PsVector' : g5gi,
                         'Tensor'   : gjgi}


FFFitFuns = {'Scalar'      : FormFactorO1,
             'Vector'      : FormFactorO2,
             'GeGm'        : FormFactorO2,
             'VectorTop'   : FormFactorO3,
             'VectorTopBNL'   : FormFactorO3,
             'VectorTopBNL_mal'   : FormFactorO3,
             'VectorWein'  : FormFactorO3,
             'VectorWeinBNL'  : FormFactorO3,
             'VectorTopNoV': FormFactorO3,
             'PsScalar'    : FormFactorO1,
             'PsVector'    : FormFactorO2,
             'Tensor'      : FormFactorO3}



NoFFPars = {'Scalar'   : 1,
            'Vector'   : 2,
            'GeGm'   : 2,
            'VectorTop': 3,
            'VectorTopBNL': 3,
            'VectorTopBNL_mal': 3,
            'VectorWein': 3,
            'VectorWeinBNL': 3,
            'VectorTopNoV': 3,
            'PsScalar' : 1,
            'PsVector' : 2,
            'Tensor'   : 3}

FFParams = {'Scalar'      : ['G_{S}'],
            'Vector'      : ['F_{1}','F_{2}'],
            'VectorTop'   : ['F_{1}','F_{2}','F_{3}'],
            'VectorTopBNL'   : ['F_{1}','tilde{F}_{2}','tilde{F}_{3}'],
            'VectorTopBNL_mal'   : ['alpha F_{1}','alpha tilde{F}_{2}','alpha tilde{F}_{3}'],
            'VectorWein'  : ['F_{1}','F_{2}','F_{3}'],
            'VectorWeinBNL'   : ['F_{1}','tilde{F}_{2}','tilde{F}_{3}'],
            'VectorTopNoV': ['F_{1}','F_{2}','F_{3}'],
            'GeGm'        : ['G_{E}','G_{M}'],
            'PsScalar'    : ['G_{P}'],
            'PsVector'    : ['G_{A}','G_{P}'],
            'Tensor'      : ['H_{T}','E_{T}','tilde{H}_{T}']}

FFQ2FitFun = {}
FFQ2FitFun['Scalar'] = {         'G_{S}'   : (DPfitfun,2)}

FFQ2FitFun['Vector'] = {         'F_{1}'   : (DPfitfun,2),
                                 'F_{2}'   : (DPfitfun,2)}

FFQ2FitFun['GeGm'] = {           'G_{E}'   : (DPfitfun,2),
                                 'G_{M}'   : (DPfitfun,2)}

FFQ2FitFun['Nucleon_GeGm'] = {   'G_{E}'  : (A_x_exp_minus_m_x,2),
                                 'G_{M}'  : (DPfitfun,2)}

FFQ2FitFun['Nucleon_Vector'] = { 'F_{1}'  : (A_x_exp_minus_m_x,2),
                                 'F_{2}'  : (DPfitfun,2)}

FFQ2FitFun['PsScalar'] = {       'G_{P}'   : (LinearFitFun,2)}

FFQ2FitFun['PsVector'] = {       'G_{A}'   : (DPfitfun,2),
                                 'G_{P}'   : (DPfitfun,2)}

FFQ2FitFun['Tensor'] = {         'H_{T}'   : (LinearFitFun,2),
                                 'E_{T}'   : (A_x_exp_minus_m_x,2),
                           'tilde{H}_{T}'  : (A_x_exp_minus_m_x,2)}

FFQ2FitFun['VectorTop'] = {      'F_{1}'   : (DPfitfun,2),
                                 'F_{2}'   : (DPfitfun,2),
                                 'F_{3}'   : (LinearFitFun,2),
                            'frac{F_{3}}{2m_{N}}'   : (LinearFitFun,2)}

FFQ2FitFun['VectorTopBNL'] = {   'F_{1}'   : (DPfitfun,2),
                                 'tilde{F}_{2}'   : (DPfitfun,2),
                                 'tilde{F}_{3}'   : (LinearFitFun,2),
                            'frac{tilde{F}_{3}}{2m_{N}}'   : (LinearFitFun,2)}

FFQ2FitFun['VectorTopBNL_mal'] = {   'alpha F_{1}'   : (DPfitfun,2),
                                 'alpha tilde{F}_{2}'   : (DPfitfun,2),
                                 'alpha tilde{F}_{3}'   : (LinearFitFun,2),
                            'alpha frac{tilde{F}_{3}}{2m_{N}}'   : (LinearFitFun,2)}

FFQ2FitFun['VectorWein'] = {     'F_{1}'   : (DPfitfun,2),
                                 'F_{2}'   : (DPfitfun,2),
                                 'F_{3}'   : (LinearFitFun,2),
                            'frac{F_{3}}{2m_{N}}'   : (LinearFitFun,2)}

FFQ2FitFun['VectorWeinBNL'] = {  'F_{1}'   : (DPfitfun,2),
                                 'tilde{F}_{2}'   : (DPfitfun,2),
                                 'tilde{F}_{3}'   : (LinearFitFun,2),
                            'frac{tilde{F}_{3}}{2m_{N}}'   : (LinearFitFun,2)}

FFQ2FitFun['VectorTopNoV'] = {   'F_{1}'   : (DPfitfun,2),
                                 'F_{2}'   : (DPfitfun,2),
                                 'F_{3}'   : (LinearFitFun,2),
                            'frac{F_{3}}{2m_{N}}'   : (LinearFitFun,2)}


FF_Imp = {}

FF_Imp['Vector'] = []
FF_Imp['Vector'].append('P4_g4')
FF_Imp['Vector'].append('P3_g1')
FF_Imp['Vector'].append('P3_g2')

FF_Imp['VectorTop'] = FF_Imp['Vector']
FF_Imp['VectorTop'].append('P3_g4_Top')

FF_Imp['VectorWein'] = FF_Imp['Vector']
FF_Imp['VectorWein'].append('P3_g4_Wein')

FF_Imp['VectorWeinBNL'] = FF_Imp['VectorWein']
FF_Imp['VectorTopBNL'] = FF_Imp['VectorTop']
FF_Imp['VectorTopBNL_mal'] = FF_Imp['VectorTop']



# normalisation for Vector Current from arXiv:1006.1164v2 theta = 0 result from table IX beta = 1.9, Iwasaki gauge action, clover term Csw=1.715
FFRenorm = {'Scalar'      : 1,
            'Vector'      : 0.7354,
            'VectorTop'   : 0.7354,
            'VectorTopBNL': 0.7354,
            'VectorTopBNL_mal': 0.7354,
            'VectorWein'  : 0.7354,
            'VectorWeinBNL': 0.7354,
            'VectorTopNoV': 0.7354,
            'GeGm'        : 0.7354,
            'PsScalar'    :1,
            'PsVector'    : 1,
            'Tensor'      : 1}



##############################################################################################



class FormFact(object):


    """

    This class should produce all form factor combinations for a particular from factor type
    (e.g. vector, pseudo-vector, tensor etc...)

    """

    def __init__(self,currtype,q2list,pp2list = [0],Info=defInfo,tflowlist=[],name=''):

        if currtype not in CurrTypes:
            raise IOError(currtype+' not recognised as a current type in FormFactors.py')
        self.currtype = currtype
        self.doflow = ('Top' in self.currtype) or ('Wein' in self.currtype)
        self.fffun = FFFitFuns[currtype]
        self.currOp = CurrOpps[currtype]
        self.currOpNP = CurrOppsNoProj[currtype]
        self.currOpNPFSig = CurrOppsNoProjSigBack[currtype]
        self.npar = NoFFPars[currtype]
        self.ffpars = FFParams[currtype]
        self.kappafolder = ''

        self.ppparams = mp.LatticeParameters(mom2list=pp2list,Info=Info)
        self.ppparams.LoadPickle()
        self.qparams = mp.LatticeParameters(mom2list=q2list,Info=Info)
        self.qparams.LoadPickle()
        self.nt = self.qparams.nt
        self.nxyz = self.qparams.nxyz
        self.dim_label = 'RC'+str(self.nxyz)+'x'+str(self.nt)
        self.dim_ll = str(self.nxyz)+r'^{3}\times ' + str(self.nt)

        tracespin = CorrSpinTrace()
        self.FFFun = tracespin[self.currtype]
        self.tflowlist = tflowlist
        # if not isinstance(mass,basestring):
        #     self.SetMass(mass)
        # else:
        #     self.mass = 'Not Set'
        # if not isinstance(alpha,basestring):
        #     self.SetAlpha(alpha,tflowlist)
        # else:
        #     self.alpha = 'Not Set'
        self.mass = 'Not Set'
        self.alpha = 'Not Set'
        self.SetCustomName(name)

        # ## self.paramDict { icurrOp , ipkey , iqkey , imass (, tflow ) , ipar } bs
        # ## (, tflow) if cp odd form factor with flow times
        # self.paramDict = 'Not Set'

        if self.doflow:
            self.param_Col_Names = ['current_operator','sink_momentum','current_momentum',
                                    'mass_and_alpha','flow_time','form_factor']
            self.paramRed_Col_Names = [ 'current_operator','sink_momentum','current_momentum',
                                        'mass_and_alpha','flow_time','form_factor']
        else:
            self.param_Col_Names = ['current_operator','sink_momentum','current_momentum',
                                    'mass','form_factor']
            self.paramRed_Col_Names = ['current_operator','sink_momentum','current_momentum',
                                    'mass','form_factor']
        self.param_Stats = pa.DataFrame()
        '''
        columns = FF FFAvg FFStd
        rows are param_Col_Names
        '''

        self.paramRed_Stats = pa.DataFrame()
        '''
        columns = FF FFAvg FFStd
        rows are param_Col_Names
        '''



    def AppendName(self):
        this_str = self.name
        if isinstance(self.mass,(dict,OrderedDict)):
            this_str += '_'+list(self.mass.keys())[0]
        if isinstance(self.alpha,(dict,OrderedDict)):
            this_str += '_'+list(self.alpha.keys())[0]+'_'+list(self.alpha.values())[0].keys()[0]
        self.SetCustomName(string=this_str)

    def SetCustomName(self,string=''):
        if string == '':
            # if len(self.tflowlist) > 0:
            #     if self.tflowlist[0] == 'none':
            #         self.name = '_'.join([  self.currtype,'q'+''.join(map(str,self.qparams.mom2list)),
            #                                 'pp'+''.join(map(str,self.ppparams.mom2list))])
            #     else:
            #         self.name = '_'.join([  self.currtype,'q'+''.join(map(str,self.qparams.mom2list)),
            #                                 'pp'+''.join(map(str,self.ppparams.mom2list)),self.tflowlist[0]])
            # else:
            self.name = '_'.join([  self.currtype,'q'+''.join(map(str,self.qparams.mom2list)),
                                    'pp'+''.join(map(str,self.ppparams.mom2list))])
        else:
            self.name = string

        self.FFcoeffDir = outputdir +'/'+self.dim_label+self.kappafolder+ 'FFcoeffs/'
        mkdir_p(self.FFcoeffDir+'/Pickle/')
        mkdir_p(self.FFcoeffDir+'/Excel/')
        self.HumanFile = self.FFcoeffDir+self.name+'.xml'
        self.ExcelFile = self.FFcoeffDir+'/Excel/'+self.name+'.xlsx'
        self.PickleFile = self.FFcoeffDir+'/Pickle/'+self.name+'.py3p'



    def GetRedFFParms(self):
        if len(self.param_Stats) == 0:
            self.GetFFParams()
        vall,keyl = [],[]
        vall_nd,qsqrdlist = [],[]
        Avgl,Stdl = [],[]
        keycombl,coeffl = [],[]
        for ikey,iFF in self.param_Stats['FF'].groupby(self.paramRed_Col_Names[:-1]):
            if isinstance(ikey,str):
                this_key = list(iFF.index[0])
                iq = this_key[2]
            else:
                this_key = ikey
                iq = ikey[2]
            thisqsqrd = self.qparams.Getpsqrdform(iq,pre='q')
            checkbool = True
            FFvals = iFF.values
            FFvalsAvg = np.array([hold_val.Avg for hold_val in FFvals])
            for ic,iffavg in enumerate(FFvalsAvg):
                if abs(iffavg) > myeps:
                    break
                elif ic == len(FFvalsAvg)-1:
                    print('Warning, form factor is all zero for:')
                    print(ikey)
            for ics,(pref_ff,prev_qsqrd) in enumerate(zip(vall_nd,qsqrdlist)):
                prev_ffAvg = np.array([hold_val.Avg for hold_val in pref_ff])
                comp_rat = prev_ffAvg[ic]/FFvalsAvg[ic]
                rat_list = []
                for denomen,numer in zip(FFvalsAvg,prev_ffAvg):
                    if denomen == 0:
                        if numer != 0:
                            rat_list.append(False)
                    else:
                        rat_list.append(numer/denomen)
                if all(rat_list == comp_rat) and prev_qsqrd == thisqsqrd:
                    checkbool = False
                    keycombl[ics].append(this_key)
                    coeffl[ics].append(comp_rat)
                    break
            if checkbool:
                keycombl.append([this_key])
                coeffl.append([1.0])
                vall += list(FFvals)
                vall_nd.append(list(FFvals))
                qsqrdlist.append(thisqsqrd)
                Avgl += [ival.Avg for ival in FFvals]
                Stdl += [ival.Std for ival in FFvals]
                keyl += list(iFF.keys())
        if len(keyl) > 0:
            indicies = pa.MultiIndex.from_tuples(keyl,names=self.paramRed_Col_Names)
            self.paramRed_Stats.loc[:,'FF'] = pa.Series(vall,index=indicies)
            self.paramRed_Stats.loc[:,'FFAvg'] = pa.Series(Avgl,index=indicies)
            self.paramRed_Stats.loc[:,'FFStd'] = pa.Series(Stdl,index=indicies)
            self.combine_keys_forRF = keycombl
            self.coeffs_forRF = coeffl
        else:
            print('Warning, error with reducing form factor list')
            print(self.param_Stats['FF'])
    ## imass is input parameter since we can use bootstrapped values
    ## wrapped by MassFFParams()
    def GetFFParms(self):
        # if thismass != 'PreDef': self.SetMass(thismass)
        # if thisalpha != 'PreDef' and self.doflow: self.SetAlpha(thisalpha,thistflowlist)

        if self.mass == 'Not Set':
            raise IOError('mass not set, either pass into MassFFParams() or use SetMass()')
        if self.alpha == 'Not Set':
            self.alpha = {}
            for itflow in self.tflowlist:
                if itflow == 'none': continue
                self.alpha[itflow] = {'Test':alpha_guess}
        plist,pAvg,pStd = [],[],[]
        ilist = []
        if self.doflow and len(self.tflowlist) == 0:
            raise IOError('no tflowlist present with flowed operators')
        looplist = [self.currOp,self.ppparams.momformlist,self.qparams.momformlist]
        thistimer = Timer(linklist = looplist,name = 'FF '+self.name)
        for icurrOp in self.currOp:
            for ip in self.ppparams.momformlist:
                for iq in self.qparams.momformlist:
                    ipvec,iqvec = list(self.ppparams.TOpvec(ip,actual=True)),list(self.qparams.TOpvec(iq,actual=True))
                    ipkey,iqkey = self.ppparams.TOpform(ip,'pp'),self.qparams.TOpform(iq,'q')
                    hold,rcheck,ccheck = self.FFFun(icurrOp,iqvec,ipvec,list(self.mass.values())[0].Avg)
                    if ccheck == True and rcheck == True:
                        print(hold)
                        raise ArithmeticError('Operator '+' '.join(icurrOp+('pp'+ip,'q'+iq)) +' has both real and imaginary components, check Form Factor coeff construction')
                    elif ccheck == True:
                        Opout = icurrOp + ('cmplx',)
                    else:
                        Opout = icurrOp
                    for imass,bootmass in self.mass.items():
                        if self.doflow:
                            for itflow in self.tflowlist:
                                for ialpha,bootalpha in self.alpha[itflow].items():
                                    holdlist = [[] for i in self.ffpars]
                                    # print bootmass.Avg,iqvec,ipvec,icurrOp
                                    # print rcheck, ccheck
                                    # print hold
                                    thisalpha = [bootalpha]*len(bootmass.bootvals)
                                    if isinstance(bootalpha,BootStrap):
                                        thisalpha = bootalpha.bootvals
                                    for bmass,balpha in zip(bootmass.bootvals,thisalpha):
                                        hold,rdump,cdump = self.FFFun(icurrOp,iqvec,ipvec,bmass,alpha=balpha)
                                        if rcheck:
                                            for icpar,ipar in enumerate(self.ffpars):
                                                holdlist[icpar].append(hold[icpar].real)
                                        if ccheck:
                                            for icpar,ipar in enumerate(self.ffpars):
                                                holdlist[icpar].append(hold[icpar].imag)
                                    for icpar,ipar in enumerate(self.ffpars):
                                        if len(holdlist[icpar]) != bootmass.nboot: continue
                                        ilist.append(('_'.join(Opout),ipkey,iqkey,'_'.join([imass,ialpha]),itflow,ipar))
                                        this_param = BootStrap(bootmass.nboot,name=ipar,bootvals=holdlist[icpar])
                                        plist.append(this_param)
                                        pAvg.append(this_param.Avg)
                                        pStd.append(this_param.Std)
                        else:
                            holdlist = [[] for i in self.ffpars]
                            for bmass in bootmass.bootvals:
                                hold,rdump,cdump = self.FFFun(icurrOp,iqvec,ipvec,bmass,alpha=self.alpha)
                                if rcheck:
                                    for icpar,ipar in enumerate(self.ffpars):
                                        holdlist[icpar].append(hold[icpar].real)
                                if ccheck:
                                    for icpar,ipar in enumerate(self.ffpars):
                                        holdlist[icpar].append(hold[icpar].imag)
                            for icpar,ipar in enumerate(self.ffpars):
                                if len(holdlist[icpar]) != bootmass.nboot: continue
                                ilist.append(('_'.join(Opout),ipkey,iqkey,imass,ipar))
                                this_param = BootStrap(bootmass.nboot,name=ipar,bootvals=holdlist[icpar])
                                plist.append(this_param)
                                pAvg.append(this_param.Avg)
                                pStd.append(this_param.Std)
                    thistimer.Lap()
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.param_Col_Names)
            self.param_Stats.loc[:,'FF'] = pa.Series(plist,index=indicies)
            self.param_Stats.loc[:,'FFAvg'] = pa.Series(pAvg,index=indicies)
            self.param_Stats.loc[:,'FFStd'] = pa.Series(pStd,index=indicies)




    ## self.mass ends up as BootStrap class
    def SetMass(self,mass):
        self.mass = ODNested()
        if isinstance(mass,TwoPointCorr):
            if 'boot' not in mass.C2_Fit_Stats:
                raise IOError('C2_Fit_Stats in TwoPointCorr not calculated, please run Fit() on TwoPointCorr instance')
            for (istate,ip),fitdata in mass.C2_Fit_Stats['boot'].items():
                if ip != 'p000': continue
                ifitr,this_data = fitdata.GetMaxFitr()
                self.mass['C2_'+'_'.join([istate,ifitr])] = this_data.iloc[0].fit_data['Params']['Energy']
        elif isinstance(mass,Fitting):
            if 'Energy' not in list(mass.fit_data['Params'].keys()):
                print(list(mass.fit_data['Params'].keys()))
                raise IOError('Fitting class instance has no parameter Energy')
            self.mass['C2_'+mass.name] = mass.fit_data['Params']['Energy']
        elif isinstance(mass,BootStrap):
            self.mass['C2_'+mass.name] = mass
        elif mass == 'Not Set':
            self.mass = 'Not Set'
        else:
            raise IOError('mass input type not recognisied (see FormFactors.py, SetMass())')

    def GetFuns(self):
        if not hasattr(self,'FFFun'):
            self.FFFun,self.fffun = ReadFuns(self.FFFun_name,self.fffun_name)

    def RemoveFuns(self):
        if hasattr(self,'FFFun') and hasattr(self.FFFun,'__name__'):
            self.FFFun_name = self.FFFun.__name__
            self.fffun_name = self.fffun.__name__
            WriteFuns(self.FFFun,self.fffun)
            del self.FFFun
            del self.fffun


    ## self.Alpha ends up as BootStrap class
    def SetAlpha(self,alpha,tflowlist):
        self.alpha = ODNested()
        self.tflowlist = tflowlist
        if isinstance(alpha,NNQFullCorr):
            # if self.kappafolder not in  ['',alpha.kappafolder+'/']:
            #     raise IOError('kappa folders used in FormFact (from FormFactors.py) have changed, check alpha and mass used.')
            # self.kappafolder = alpha.kappafolder + '/'
            # self.SetCustomName(string=self.name)
            if 'boot' not in alpha.NNQFull_Fit_Stats or 'p000' not in list(alpha.NNQFull_Fit_Stats['boot'].keys()):
                raise IOError('No zero momentum found in NNQFull_Fit_Stats in TwoPointCorr, try running Fit()')
            for itflow in tflowlist:
                if itflow not in list(alpha.NNQFull_Fit_Stats['boot']['p000'].keys()):
                    print(list(alpha.NNQFull_Fit_Stats['boot']['p000'].keys()), itflow)
                    raise IOError('No tflow found in Fit Dict in NNQFullCorr, try funning Fit() with correct flow time list')
                ifitr,this_data = alpha.NNQFull_Fit_Stats.at[('p000',itflow),'boot'].GetMaxFitr()
                ifitr = ifitr.replace('fittwor','tsumfitr')
                self.alpha[itflow]['alpha_'+ifitr] = this_data.iloc[0].fit_data['Params'].iloc[0]
            return
        elif isinstance(alpha,NNQCorr):
            # if self.kappafolder not in  ['',alpha.kappafolder+'/']:
            #     raise IOError('kappa folders used in FormFact (from FormFactors.py) have changed, check alpha and mass used.')
            # self.kappafolder = alpha.kappafolder + '/'
            # self.SetCustomName(string=self.name)
            if 'boot' not in alpha.NNQ_Fit_Stats or 'p000' not in list(alpha.NNQ_Fit_Stats['boot'].keys()):
                raise IOError('No zero momentum found in NNQ_Fit_Stats in TwoPointCorr, try running Fit()')
            for itflow in tflowlist:
                if itflow not in list(alpha.NNQ_Fit_Stats['boot']['p000'].keys()):
                    print(list(alpha.NNQ_Fit_Stats['boot']['p000'].keys()), itflow)
                    raise IOError('No tflow found in Fit Dict in NNQCorr, try funning Fit() with correct flow time list')
                ifitr,this_data = alpha.NNQ_Fit_Stats.at[('p000',itflow),'boot'].GetMaxFitr()
                self.alpha[itflow]['alpha_'+ifitr] = this_data.iloc[0].fit_data['Params'].iloc[0]
            return
        try:
            firstflow = list(alpha.keys())[0]
            if isinstance(alpha[firstflow],Fitting):
                for itflow in tflowlist:
                    if itflow not in list(alpha.keys()):
                        raise IOError('No tflow found in alpha.keys(), alpha can be dictionary with tflows')
                    if 'Energy' not in list(alpha.fit_data['Params'].keys()):
                        print(list(alpha[itflow].fit_data['Params'].keys()))
                        raise IOError('Fitting class instance has no parameter Energy')
                    self.alpha[itflow]['alpha_'+alpha.name] = alpha.fit_data['Params'].iloc[0]
            elif isinstance(alpha[firstflow],BootStrap):
                for itflow in tflowlist:
                    if itflow not in list(alpha.keys()):
                        raise IOError('No tflow found in alpha.keys(), alpha can be dictionary with tflows')
                    self.alpha[itflow]['alpha_'+alpha.name] = alpha
        except:
            if alpha == 'No Alpha' or alpha == 'PreDef':
                self.alpha = 'Not Set'
            elif self.alpha != 'Not Set':
                print(type(alpha))
                if isinstance(alpha,str):
                    print(alpha)
                raise IOError('alpha input type not recognisied (see FormFactors.py, SetAlpha())')





    def __getitem__(self,ikey):
        return self.param_Stats['FF'][ikey]

    ## if data is instance of FormFact (e.g. data = FormFact(...)
    ## for iOp,ipkey,iqkey,imass,ipar,data in data.items():
    ##    etc.....
    def items(self):
        return list(self.param_Stats['FF'].items())

    def iteritems(self):
        return iter(self.param_Stats['FF'].items())

    def iteritemsRed(self):
        return iter(self.paramRed_Stats['FF'].items())

    def itemsMass(self):
        ## self.paramDict { icurrOp , ipkey , iqkey , imass } bs
        output = []
        for icurrOp,currdata in self.paramDict.items():
            for ipkey,pdata in currdata.items():
                for iqkey,qdata in pdata.items():
                    output.append((icurrOp,ipkey,iqkey,qdata))
        return output

    def itemsParList(self):
        ## self.paramDict { icurrOp , ipkey , iqkey , imass } bs
        output = []
        if self.doflow:
            this_levels = ('current_operator','sink_momentum','current_momentum','mass_and_alpha','flow_time')
            for (icurrOp,ipkey,iqkey,imass,itflow),flow_data in self.param_Stats['FF'].groupby(level=this_levels):
                this_fdata = flow_data[icurrOp,ipkey,iqkey,imass,itflow]
                output.append((icurrOp,ipkey,iqkey,imass,itflow,list(this_fdata.keys()),this_fdata.values))
        else:
            this_levels = ('current_operator','sink_momentum','current_momentum','mass')
            for (icurrOp,ipkey,iqkey,imass),flow_data in self.param_Stats['FF'].groupby(level=this_levels):
                this_fdata = flow_data[icurrOp,ipkey,iqkey,imass]
                output.append((icurrOp,ipkey,iqkey,imass,list(this_fdata.keys()),this_fdata.values))
        return output


    def itemsParListRed(self):
        ## self.paramDict { icurrOp , ipkey , iqkey , imass } bs
        output = []
        if self.doflow:
            this_levels = ('current_operator','sink_momentum','current_momentum','mass_and_alpha','flow_time')
            for (icurrOp,ipkey,iqkey,imass,itflow),flow_data in self.paramRed_Stats['FF'].groupby(level=this_levels):
                this_fdata = flow_data[icurrOp,ipkey,iqkey,imass,itflow]
                output.append((icurrOp,ipkey,iqkey,imass,itflow,list(this_fdata.keys()),this_fdata.values))
        else:
            this_levels = ('current_operator','sink_momentum','current_momentum','mass')
            for (icurrOp,ipkey,iqkey,imass),flow_data in self.paramRed_Stats['FF'].groupby(level=this_levels):
                this_fdata = flow_data[icurrOp,ipkey,iqkey,imass]
                output.append((icurrOp,ipkey,iqkey,imass,list(this_fdata.keys()),this_fdata.values))
        return output


    def itemsAvgStd(self):
        outlist = []
        for ivals in self.param_Stats.itertuples():
            outlist.append((ivals[0],ivals[2],ivals[3]))
        return outlist

    def values(self):
        ## self.paramDict { icurrOp , ipkey , iqkey , imass , ipar } bs
        return self.param_Stats['FF'].values

    def valuesAvgStd(self):
        ## self.paramDict { icurrOp , ipkey , iqkey , imass , ipar } bs
        outlist = []
        for ivals in self.param_Stats.itertuples():
            outlist.append((ivals[2],ivals[3]))
        return outlist

    def keys(self):
        return list(self.param_Stats.keys())

    def keysOpMoms(self):
        ## self.paramDict { icurrOp , ipkey , iqkey , imass , ipar } bs
        output = []
        for (icurrOp,ipkey,iqkey),idata in self.param_Stats['FF'].groupby(level=('current_operator','sink_momentum','current_momentum')):
            output.append((icurrOp.split('_'),ipkey,iqkey))
        return output

    def keysOpMomsRed(self):
        ## self.paramDict { icurrOp , ipkey , iqkey , imass , ipar } bs
        output = []
        for (icurrOp,ipkey,iqkey),idata in self.paramRed_Stats['FF'].groupby(level=('current_operator','sink_momentum','current_momentum')):
            output.append((icurrOp.split('_'),ipkey,iqkey))
        return output

    def Write(self):
        def FixDictArray(thislist,ikey):
            outDict = OrderedDict()
            for ic,ilist in enumerate(thislist):
                outDict[ikey+str(ic+1)] = ilist
            return outDict

        self.outDict = ODNested()
        self.outDict['Type'] = self.currtype
        self.outDict['FFparams'] = FixDictArray(self.ffpars,'FF')
        self.outDict['tflowlist'] = FixDictArray(self.tflowlist,'flow')
        ## TODO, output more information about momenta if you would like
        self.outDict['nt'] = self.nt
        self.outDict['nxyz'] = self.nxyz
        self.outDict['q2max'] = self.qparams.mom2max
        self.outDict['q2min'] = self.qparams.mom2min
        self.outDict['pp2max'] = self.ppparams.mom2max
        self.outDict['pp2min'] = self.ppparams.mom2min
        excel_params = pa.Series(deepcopy(self.outDict))
        firstkey = list(self.mass.keys())[0]
        excel_params['massAvg_'+firstkey] = list(self.mass.values())[0].Avg
        excel_params['massStd_'+firstkey] = list(self.mass.values())[0].Std
        for icFF,iFF in excel_params['FFparams'].items():
            excel_params[icFF] = iFF
        for ictflow,itflow in excel_params['tflowlist'].items():
            excel_params[ictflow] = itflow
        del excel_params['FFparams']
        del excel_params['tflowlist']

        if self.doflow:
            for itflow in self.tflowlist:
                firstdata = list(self.alpha[itflow].values())[0]
                if isinstance(firstdata,BootStrap):
                    excel_params['AlphaAvg_'+firstkey+'_'+itflow] = firstdata.Avg
                else:
                    excel_params['AlphaAvg_'+firstkey+'_'+itflow] = firstdata
            for itflow in self.tflowlist:
                firstdata = list(self.alpha[itflow].values())[0]
                if isinstance(firstdata,BootStrap):
                    excel_params['AlphaStd_'+firstkey+'_'+itflow] = firstdata.Std

        for ikey,massdata in self.mass.items():
            if isinstance(massdata,BootStrap):
                self.outDict['mass'][ikey] = AvgStdToFormat(massdata.Avg,massdata.Std)
            else:
                self.outDict['mass'][ikey] = massdata
        if self.doflow:
            for itflow,alphatflow in self.alpha.items():
                firstalpha = list(alphatflow.values())[0]
                alpha_key = list(alphatflow.keys())[0]
                if isinstance(firstalpha,BootStrap):
                    self.outDict['alpha'][alpha_key][itflow] = AvgStdToFormat(firstalpha.Avg,firstalpha.Std)
                else:
                    self.outDict['alpha'][alpha_key][itflow] = firstalpha
        if self.doflow:
            for (gammakey,ipp,iq,imass,itflow,ipar),paramdata in self.items():
                thisqsqrd = self.qparams.Getpsqrdform(iq,pre='q')
                ## TODO, this is incorrectly order with respect to qsqrd, need to write function to fix this somewhere
                ## also, you can play around with the ordering of keys to get the formating you like !!
                # self.outDict[imass][ipp][itflow][self.qparams.Getpsqrdform(iq,pre='q')]['Coeff_of_'+ipar]['_'.join([gammakey,iq])] = AvgStdToFormat(paramdata.Avg,paramdata.Std)
                try:
                    self.outDict['Equations'][imass][ipp][itflow][thisqsqrd]['_'.join([gammakey,iq])] += ' ('+AvgStdToFormat(paramdata.Avg,paramdata.Std,dec=5,numlen=8)+ ') '+ipar + ' +'
                except:
                    self.outDict['Equations'][imass][ipp][itflow][thisqsqrd]['_'.join([gammakey,iq])] = ' ('+AvgStdToFormat(paramdata.Avg,paramdata.Std,dec=5,numlen=8)+ ') '+ipar + ' +'
                if ipar == self.ffpars[-1]:
                    self.outDict['Equations'][imass][ipp][itflow][thisqsqrd]['_'.join([gammakey,iq])] = self.outDict['Equations'][imass][ipp][itflow][thisqsqrd]['_'.join([gammakey,iq])][:-1]
        else:
            for (gammakey,ipp,iq,imass,ipar),paramdata in self.items():
                thisqsqrd = self.qparams.Getpsqrdform(iq,pre='q')
                ## TODO, this is incorrectly order with respect to qsqrd, need to write function to fix this somewhere
                ## also, you can play around with the ordering of keys to get the formating you like !!
                try:
                    self.outDict['Equations'][imass][ipp][thisqsqrd]['_'.join([gammakey,iq])] += ' ('+AvgStdToFormat(paramdata.Avg,paramdata.Std,dec=5,numlen=8)+ ') '+ipar + ' +'
                except:
                    self.outDict['Equations'][imass][ipp][thisqsqrd]['_'.join([gammakey,iq])] = ' ('+AvgStdToFormat(paramdata.Avg,paramdata.Std,dec=5,numlen=8)+ ') '+ipar + ' +'
                if ipar == self.ffpars[-1]:
                    self.outDict['Equations'][imass][ipp][thisqsqrd]['_'.join([gammakey,iq])] = self.outDict['Equations'][imass][ipp][thisqsqrd]['_'.join([gammakey,iq])][:-1]



        if self.doflow:
            for (gammakey,ipp,iq,imass,itflow,ipar),paramdata in self.iteritemsRed():
                thisqsqrd = self.qparams.Getpsqrdform(iq,pre='q')
                ## TODO, this is incorrectly order with respect to qsqrd, need to write function to fix this somewhere
                ## also, you can play around with the ordering of keys to get the formating you like !!
                # self.outDict[imass][ipp][itflow][self.qparams.Getpsqrdform(iq,pre='q')]['Coeff_of_'+ipar]['_'.join([gammakey,iq])] = AvgStdToFormat(paramdata.Avg,paramdata.Std)
                try:
                    self.outDict['Equations_Reduced'][imass][ipp][itflow][thisqsqrd]['_'.join([gammakey,iq])] += ' ('+AvgStdToFormat(paramdata.Avg,paramdata.Std,dec=5,numlen=8)+ ') '+ipar + ' +'
                except:
                    self.outDict['Equations_Reduced'][imass][ipp][itflow][thisqsqrd]['_'.join([gammakey,iq])] = ' ('+AvgStdToFormat(paramdata.Avg,paramdata.Std,dec=5,numlen=8)+ ') '+ipar + ' +'
                if ipar == self.ffpars[-1]:
                    self.outDict['Equations_Reduced'][imass][ipp][itflow][thisqsqrd]['_'.join([gammakey,iq])] = self.outDict['Equations_Reduced'][imass][ipp][itflow][thisqsqrd]['_'.join([gammakey,iq])][:-1]
        else:
            for (gammakey,ipp,iq,imass,ipar),paramdata in self.iteritemsRed():
                thisqsqrd = self.qparams.Getpsqrdform(iq,pre='q')
                ## TODO, this is incorrectly order with respect to qsqrd, need to write function to fix this somewhere
                ## also, you can play around with the ordering of keys to get the formating you like !!
                try:
                    self.outDict['Equations_Reduced'][imass][ipp][thisqsqrd]['_'.join([gammakey,iq])] += ' ('+AvgStdToFormat(paramdata.Avg,paramdata.Std,dec=5,numlen=8)+ ') '+ipar + ' +'
                except:
                    self.outDict['Equations_Reduced'][imass][ipp][thisqsqrd]['_'.join([gammakey,iq])] = ' ('+AvgStdToFormat(paramdata.Avg,paramdata.Std,dec=5,numlen=8)+ ') '+ipar + ' +'
                if ipar == self.ffpars[-1]:
                    self.outDict['Equations_Reduced'][imass][ipp][thisqsqrd]['_'.join([gammakey,iq])] = self.outDict['Equations_Reduced'][imass][ipp][thisqsqrd]['_'.join([gammakey,iq])][:-1]


        ## Write human readable file with data
        WriteXml(self.HumanFile,{'Results':self.outDict})

        WriteExcel(self.ExcelFile,{'FF_results':deepcopy(self.param_Stats)},params=excel_params)

        ## pickles rest of data for reading
        self.RemoveFuns()
        WritePickle(self.PickleFile,self.__dict__)
        self.GetFuns()

    def LoadKappaDir(self,mass,alpha='None'):
        if isinstance(mass,TwoPointCorr):
            if self.kappafolder.replace('/','') not in  ['',mass.kappafolder]:
                raise IOError('kappa folders used in FormFact (from FormFactors.py) have changed, check alpha and mass used.')
            self.kappafolder = mass.kappafolder + '/'
        if isinstance(alpha,NNQFullCorr) or isinstance(alpha,NNQCorr):
            if self.kappafolder.replace('/','') not in  ['',alpha.NN.kappafolder]:
                raise IOError('kappa folders used in FormFact (from FormFactors.py) have changed, check alpha and mass used.')
            self.kappafolder = alpha.NN.kappafolder + '/'
        self.SetCustomName(string=self.name)


    def LoadPickle(self,thismass,thisalpha='No Alpha',thistflowlist=[],DefWipe=False):
        # if self.doflow and thisalpha == 'No Alpha':
        #     raise EnvironmentError('alpha is needed for FormFact when using FO form factor.')
        self.LoadKappaDir(thismass,alpha=thisalpha)
        self.SetMass(thismass)
        self.SetAlpha(thisalpha,thistflowlist)
        self.AppendName()

        # print 'loading ', self.PickleFile, 'with DefWipe ' , DefWipe
        if os.path.isfile(self.PickleFile) and not DefWipe:
            # print 'Loading Pickle for ' , self.name

            thispickle = ReadPickleWrap(self.PickleFile)
            # print thistflowlist, thispickle['tflowlist'],self.doflow
            if len(thistflowlist) == 0 or thistflowlist == thispickle['tflowlist'] or (not self.doflow):
                self.__dict__.update(thispickle)
                self.Write()
            else:
                if os.path.isfile(self.PickleFile+'.bak'):
                    print('using backupfile for '+self.name)
                    self.PickleFile = self.PickleFile+'.bak'
                    self.HumanFile = self.HumanFile+'.bak'
                    self.ExcelFile = self.ExcelFile.replace('.xlsx','.bak.xlsx')
                    self.LoadPickle(thismass,thisalpha=thisalpha,thistflowlist=thistflowlist,DefWipe=DefWipe)
                    return
                self.GetFFParms()
                self.GetRedFFParms()
                self.Write()
        else:
            self.GetFFParms()
            self.GetRedFFParms()
            self.Write()



def TestFF(DefWipe=False,this_curr='VectorTop',qsqrdlist=[0,1]):
    import TwoPtCorrelators as co
    from Params import defInfo,defInfoFlow
    # import numpy as np
    data = co.TwoPointCorr(Info=defInfo)
    data.LoadPickle()
    data.Fit()
    # dataff = FormFact('Vector',[0,1])
    # dataff.LoadPickle(data)
    defInfoFlow['Interp'] = 'CPOdd'
    dataalpha = co.NNQCorr(Info=defInfoFlow)
    dataalpha.LoadPickle(DefWipe=DefWipe)
    dataalpha.FitAlpha()
    dataffalpha = FormFact( this_curr,qsqrdlist,tflowlist=['t_f6.01'])
    dataffalpha.LoadPickle(data,thisalpha=dataalpha,thistflowlist=['t_f6.01'])
    return dataffalpha


if __name__ == '__main__':
    data = TestFF()
