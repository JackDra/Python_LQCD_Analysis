#!/usr/bin/env python

from MiscFuns import ODNested
import MomParams as mp
from Params import outputdir,defxlimAlpha,chitcoeff,TflowToPhys
import numpy as np
from XmlFormatting import untstr,tflowstr,untflowstr

"""

Takes Ahmeds results and makes class out of them

"""



class AResAlpha(object):

    """

    His results for alpha ratio

    """

    def __init__(self, Info={}, thissym='Not Set',thiscol='Not Set',thisshift=0.0,name='',LLname='',Auto=True):


        ## A infront of variable is for Ahmed file convention

        ## Initialising plotting parameters
        self.thisshift = thisshift
        self.thiscol = thiscol
        self.thissym = thissym

        # Opboot { it } bs
        self.alphaAvg = ODNested()
        self.alphaStd = ODNested()

        ## kappa up down (in integer form)
        if 'kud' in list(Info.keys()): self.kud = 'kud'+str(Info['kud'])
        else: self.kud = ''

        ## kappa strange (in integer form)
        if 'ks' in list(Info.keys()): self.ks = 'ks'+str(Info['ks'])
        else: self.ks = ''

        ## Current obersvables are TopCharge and Weinberg
        if 'Observable' in list(Info.keys()): self.Observ = Info['Observable']
        else: self.Observ = 'TopCharge' ## default to TopCharge


        self.latparams = mp.LatticeParameters(Info=Info)
        self.latparams.LoadPickle()


        if self.Observ == 'TopCharge':
            self.thisconv = chitcoeff(self.latparams.nt,self.latparams.nxyz,self.latparams.latspace,'QQ')
        else:
            self.thisconv = chitcoeff(self.latparams.nt,self.latparams.nxyz,self.latparams.latspace,'WW')

        ## Ahmed only gave me a single flow time result, this will be changed when he gives me all flow times.
        if 'flowtime' in list(Info.keys()): self.flowtime = Info['flowtime']
        else: self.flowtime = 't_f6.00' ## default to read all flow times

        self.kappafolder = 'Kud0'+self.kud.replace('kud','')+'Ks0'+self.ks.replace('ks','')
        self.AObserv = self.Observ.lower().replace('weinberg','weinopp')
        self.Akud = self.kud.replace('kud','K')[:-2]
        self.Aflowtime = self.flowtime.replace('t_f','').replace('.','')
        self.Auto = Auto

        self.SetCustomName(name,LLname)
        thisdir = outputdir+'/AhmedResults/'
        if self.Auto:
            self.readdir = thisdir + 'ratio_'+self.AObserv+'_'+self.Akud+'_at_flow_time_'+self.Aflowtime+'/'
        else:
            self.readdir = thisdir + 'ratio_boot_'+self.AObserv+'_'+self.Akud+'_at_flow_time_'+self.Aflowtime+'/'
        self.thisfile = 'Tauint_at_flow_time_'+self.Aflowtime+'time_value_dvalue.txt'




    def ReadAlpha(self):
        print('reading ',self.readdir+self.thisfile)
        with open(self.readdir+self.thisfile,'r') as f:
            ## first line has info
            f.readline()
            for line in f:
                fmtline = line.split()
                self.alphaAvg['t'+str(fmtline[0])] = np.float(fmtline[1])
                self.alphaStd['t'+str(fmtline[0])] = np.float(fmtline[2])

    def PlotAlpha(self,xlims=defxlimAlpha,thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine'):
        if thisshift != 'PreDefine':
            xlen = np.abs(xlims[1]-xlims[0])
            self.thisshift = thisshift*xlen
        if thissym != 'PreDefine': self.thissym = thissym
        if thiscol != 'PreDefine': self.thiscol = thiscol
        xvals,yvals,yerr = [],[],[]
        for (it,iy),iyerr in zip(iter(self.alphaAvg.items()),iter(self.alphaStd.values())):
            itplot = untstr(it)
            if xlims[0]-1 < itplot < xlims[1]+1:
                xvals.append(itplot)
                yvals.append(iy/self.thisconv)
                yerr.append(iyerr/self.thisconv)
        if len(list(self.alphaAvg.keys())) == 0:
            raise EnvironmentError('Please read in Ahmend results before plotting with self.ReadAlpha()')
        if self.thiscol == 'Not Set' or self.thissym == 'Not Set':
            raise EnvironmentError('plotting color or symbol not set, please initialise class or set when plotting')
        # pl.errorbar(np.array(xvals)+self.thisshift,yvals,yerr,label=self.LegLab,fmt=self.thissym,color=self.thiscol)


    def SetCustomName(self,string='',stringLL=''):
        if string == '':
            self.name = self.filename = '_'.join([self.kappafolder,self.Observ])
        else:
            self.filename = self.name = string
        if stringLL == '':
            if self.Auto:
                self.LegLab = '$'+'\ '.join([self.latparams.GetPionMassLab(),self.Observ,'Ahmed'])+'$ Auto' ## customise this how you like
            else:
                self.LegLab = '$'+'\ '.join([self.latparams.GetPionMassLab(),self.Observ,'Ahmed'])+'$' ## customise this how you like
        else:
            self.LegLab = stringLL




################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################
################################################################################################################################################




class AResChi2(object):

    """

    His results for <Q^{2}> or <W^{2}> (WITHOUT COEFFICIENT)

    """

    def __init__(self, Info={}, thissym='Not Set',thiscol='Not Set',thisshift=0.0,name='',LLname=''):


        ## A infront of variable is for Ahmed file convention

        ## Initialising plotting parameters
        self.thisshift = thisshift
        self.thiscol = thiscol
        self.thissym = thissym

        # Opboot { it } bs
        self.Chi2Avg = ODNested()
        self.Chi2Std = ODNested()

        ## kappa up down (in integer form)
        if 'kud' in list(Info.keys()): self.kud = 'kud'+str(Info['kud'])
        else: self.kud = ''

        ## kappa strange (in integer form)
        if 'ks' in list(Info.keys()): self.ks = 'ks'+str(Info['ks'])
        else: self.ks = ''

        ## Current obersvables are TopCharge and Weinberg
        if 'Observable' in list(Info.keys()): self.Observ = Info['Observable']
        else: self.Observ = 'TopCharge' ## default to TopCharge

        self.latparams = mp.LatticeParameters(Info=Info)
        self.latparams.LoadPickle()


        if 'TopCharge' in self.Observ:
            self.thisconv = chitcoeff(self.latparams.nt,self.latparams.nxyz,self.latparams.latspace,'QQ')
        else:
            self.thisconv = chitcoeff(self.latparams.nt,self.latparams.nxyz,self.latparams.latspace,'WW')


        # ## this is used display only these flowtimes, Still reads and writes all
        # if 'tflowlist' in Info.keys(): self.tflowlist = Info['tflowlist']
        # else: self.tflowlist = defxlimOp ## default to read all flow times


        # listout = []
        # for ict,itflow in enumerate(self.tflowlist):
        #     if 't_f' not in str(itflow):
        #         listout.append(tflowstr(itflow))
        #     else:
        #         listout.append(itflow)
        # self.tflowfloat = [untflowstr(itflow) for itflow in listout]
        # self.tflowlist = np.array(listout)


        self.kappafolder = 'Kud0'+self.kud.replace('kud','')+'Ks0'+self.ks.replace('ks','')
        if 'topcharge' in self.Observ.lower().replace('weinberg','weinopp'):
            self.AObserv = 'topcharge'
        elif 'weinopp' in self.Observ.lower().replace('weinberg','weinopp'):
            self.AObserv = 'weinopp'
        self.Akud = self.kud.replace('kud','')[:-2]

        self.SetCustomName(name,LLname)
        thisdir = outputdir+'/AhmedResults/'
        self.readdir = thisdir + self.AObserv+self.Akud+'/'
        self.thisfile = self.AObserv+self.Akud+'.txt'




    def ReadChi2(self):
        print('reading ',self.readdir+self.thisfile)
        with open(self.readdir+self.thisfile,'r') as f:
            ## first line has info
            f.readline()
            for line in f:
                fmtline = line.split()
                tflowkey = tflowstr(float(fmtline[0])/100)
                self.Chi2Avg[tflowkey] = np.float(fmtline[1])
                self.Chi2Std[tflowkey] = np.float(fmtline[2])

    def CorrectAhmed(self,iy,iyerr):
        if 'TopCharge' in self.Observ or 'Q' in self.Observ:
            return iy**(1/4.),iy**(-3/4.)*iyerr/4.
        else:
            return iy**(1/8.),iy**(-7/8.)*iyerr/8.

    # def CorrectAhmed(self,iy,iyerr):
    #     return iy,iyerr


    def PlotChi2(self,xlims='data',thiscol='PreDefine',thissym='PreDefine',thisshift='PreDefine'):
        if len(list(self.Chi2Avg.keys())) == 0:
            raise EnvironmentError('Please read in Ahmend results before plotting with self.ReadChi2()')
        if isinstance(xlims,str) and xlims == 'data':
            xlims = list(self.Chi2Avg.keys())
        if (not isinstance(xlims[0],str)) or 't_f' not in xlims[0]:
            xlims = list(map(tflowstr,xlims))
        if thisshift != 'PreDefine':
            xlen = np.abs(xlims[1]-xlims[0])
            self.thisshift = thisshift*xlen
        if thissym != 'PreDefine': self.thissym = thissym
        if thiscol != 'PreDefine': self.thiscol = thiscol
        if self.thiscol == 'Not Set' or self.thissym == 'Not Set':
            raise EnvironmentError('plotting color or symbol not set, please initialise class or set when plotting')
        xvals,yvals,yerr = [],[],[]
        for (it,iy),iyerr in zip(iter(self.Chi2Avg.items()),iter(self.Chi2Std.values())):
            itplot = TflowToPhys(np.array([untflowstr(it)]),self.latparams.latspace)[0]
            if it in xlims:
                iycorr,iyerrcorr = self.CorrectAhmed(iy,iyerr)
                xvals.append(itplot)
                yvals.append(iycorr*self.thisconv)
                yerr.append(iyerrcorr*self.thisconv)
        # pl.errorbar(np.array(xvals)+self.thisshift,yvals,yerr,label=self.LegLab,fmt=self.thissym,color=self.thiscol)

    def SetCustomName(self,string='',stringLL=''):
        if string == '':
            self.name = self.filename = '_'.join([self.kappafolder,self.Observ])
        else:
            self.filename = self.name = string
        if stringLL == '':
            self.LegLab = '$'+'\ '.join([self.latparams.GetPionMassLab(),self.Observ,'Ahmed'])+'$' ## customise this how you like
        else:
            self.LegLab = stringLL
