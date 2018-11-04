#!/usr/bin/env python


# IMPORT THIS FIRST, put in base import file
# import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
# import matplotlib.pyplot as pl
from XmlFormatting import MakeValAndErr
from FileIO import WriteFuns,ReadFuns
from Params import defSparam
from MiscFuns import logNA
from BootStrapping import BlockCfgs,BootStrap
from Params import nboot
import numpy as np
import pandas as pa
from copy import deepcopy
# from MiscFuns import MIndex_to_numpy

##

"""
How to use:

"""


startseed = 1234


class AutoCorrelate(object):

    parentlist = ['SetOfCorrs.SetOfTwoPt','TwoPtCorrelators.TwoPointCorr','TwoPtCorrelators.TwoPointCorr','TwoPtCorrelators.NNQCorr']

    def __init__(self,  name='',Fun='NotDef', data='None',Sparam=defSparam,save_covar=True,WipeData=True):

        ## is either self.values [ ivar , monte_time ]
        ## or with replicas, self.values [ ival , ireplica , monte_time ]
        self.values = 'Not Present' ## contains regular values
        self.nval = 0.0 ## contains length of values
        ## TODO, implement numerical estimated of defivative
        ## Fun is (function,Derivate)
        if not isinstance(data,str): self.values = data
        self.Fun,self.FunDer = 'NotDef','NotDef'
        if not isinstance(Fun,str): self.Fun,self.FunDer = Fun

        self.has_reps = False

        self.Sparam = Sparam

        self.name = name ## Give it a label if you want, might help with plotting?
        self.Avg = np.float64('NaN')
        self.Std = np.float64('NaN')
        self.current = 0

        self.tau = np.float64('NaN')
        self.tauerr = np.float64('NaN')

        self.StdW = []
        self.tauW = []
        self.tauerrW = []
        self.Wopt = np.float64('NaN')

        if isinstance(data,(list,tuple,np.ndarray,pa.DataFrame)):
            self.ImportData(data)
            if self.Fun != 'NotDef' or self.FunDer != 'NotDef':
                self.uWerrMine(WipeData=WipeData,save_covar=save_covar)

    ## TODO overload operators
    def AddCoeff(self,coeff):
        self.Avg = coeff+self.Avg

    def MultCoeff(self,coeff):
        self.Avg = coeff*self.Avg
        self.Std = coeff*self.Std

    def Power(self,thispow):
        self.Avg = self.Avg**thispow
        self.Std = self.Std * self.Avg**(thispow-1)*thispow

    def RemoveVals(self):
        self.values = None
        self.nval = 0.0

    ## is either self.values [ ivar , monte_time ]
    ## or with replicas, self.values [ ival , ireplica , monte_time ]
    def ImportData(self,values,this_index=None):
        if not isinstance(values,pa.DataFrame):
            if isinstance(values,(list,np.ndarray)):
                self.values = pa.DataFrame()
                for icv,ival in enumerate(values):
                    self.values['data'+str(icv)] = pa.Series(ival,index=this_index)
            else:
                raise EnvironmentError('Importing data into BootStrap class requires numpy array or list')
        else:
            self.values = values
        self.has_reps = isinstance(self.values.index,pa.MultiIndex) and (len(self.values.index.levels[0]) > 1)
        self.nval = len(self.values)
        # if self.has_reps: self.WarnRepLen()

    # def WarnRepLen(self):
    #     for ival in self.values:
    #         if any([len(ival[0]) != len(jval) for jval in ival]):
    #             print 'Warning, replicas are not of same length'
    #             print 'Be weary of configuration list'

    def ImportFun(self,thisfun):
        self.Fun = thisfun
        self.RemoveFuns()

    def BlockCfgs(self,order):
        out_df = pa.DataFrame(columns=self.values.columns)
        ## assumes replicas are first index in pandas.MultiIndex
        if self.has_reps:
            for icol,idata in self.values.items():
                for irep,repdata in idata.groupby(level=[self.values.index.names[0]]):
                    vals_cfgs = BlockCfgs(repdata.values,order)
                    tuple_index = [(irep,i) for i in range(len(vals_cfgs))]
                    mindex = pa.MultiIndex.from_tuples(tuple_index,names=self.values.index.names)
                    out_df[icol] = pa.Series(vals_cfgs,index=mindex)
        else:
            for ikey,idata in self.values.items():
                vals_cfgs = BlockCfgs(idata.values,order)
                out_df[ikey] = pa.Series(vals_cfgs)
        return out_df



    def PlotWopt(self,plot_class,thiscol=False,thisshift=0.0,thissym=False):
        if self.Wopt == 'Not Set':
            raise IOError('Wopt not set, please import and do werr')
        xmax = int(self.Wopt*3)
        step = int(np.ceil(self.Wopt/20)) or 1
        thisshift = xmax*0.002
        GFtplot = self.GFt[:xmax:step]/self.GFt[0]
        # GFtplot = GFt[:xmax:step]
        thisshift = len(GFtplot)*thisshift
        plot_series = pa.Series()
        plot_series['type'] = 'scatter'
        plot_series['x_data'] = np.array(list(range(len(GFtplot))))
        plot_series['y_data'] = GFtplot
        plot_series['shift'] = thisshift
        plot_series['symbol'] = thissym
        plot_series['color'] = thiscol
        if self.name is not None: plot_series['label'] = self.name
        plot_class.AppendData(plot_series)
        plot_series = pa.Series()
        plot_series['type'] = 'vline'
        plot_series['x_data'] = GFtplot
        plot_series['shift'] = 'previous'
        plot_series['color'] = 'previous'
        if self.name is not None: plot_series['label'] = self.name + '_Wopt'
        plot_class.AppendData(plot_series)
        return plot_class

    def PlotTauInt(self,plot_class,thiscol=False,thisshift=0.0,thissym=False):
        if self.Wopt == 'Not Set':
            raise IOError('Wopt not set, please import and do werr')
        xmax = int(self.Wopt*3)
        step = int(np.ceil(self.Wopt/20)) or 1
        tauintplot = self.tauW[:xmax:step]
        # tauintplot = CFW[:xmax:step]
        len(tauintplot)*thisshift
        plot_series = pa.Series()
        plot_series['type'] = 'error_bar'
        plot_series['x_data'] = np.array(list(range(len(tauintplot))))
        plot_series['y_data'] = tauintplot
        plot_series['yerr_data'] = np.array(self.tauerrW)[:xmax:step]
        plot_series['shift'] = thisshift
        plot_series['symbol'] = thissym
        plot_series['color'] = thiscol
        if self.name is not None: plot_series['label'] = self.name
        plot_class.AppendData(plot_series)
        plot_series = pa.Series()
        plot_series['type'] = 'vline'
        plot_series['x_data'] = self.Wopt
        plot_series['shift'] = 'previous'
        plot_series['color'] = 'previous'
        if self.name is not None: plot_series['label'] = self.name + '_Wopt'
        plot_class.AppendData(plot_series)
        return plot_class





        # xmax = int(Wopt*3)
        # step = int(np.ceil(Wopt/20)) or 1
        # thisshift = xmax*0.002
        # fig = pl.figure(1)
        # Gplt= fig.add_subplot(211)
        # Gplt.set_ylabel(r'$\Gamma$')
        # Gplt.set_xlabel('$W$')
        # GFtplot = GFt[:xmax:step]/GFt[0]
        # # GFtplot = GFt[:xmax:step]
        # pl.errorbar(range(len(GFtplot)), GFtplot,fmt="o", color='b')
        # pl.axvline(Wopt+thisshift, color='r')
        # tplt = fig.add_subplot(212)
        # tplt.set_ylabel(r'$\tau_{\mathrm{int}}$')
        # tplt.set_xlabel(r'$W$')
        # tauintplot = tauint[:xmax:step]
        # # tauintplot = CFW[:xmax:step]
        # pl.errorbar(range(len(tauintplot)), tauintplot,
        #              dtauint[:xmax:step], fmt="o", color='b')
        # pl.axvline(Wopt+thisshift, color='r')
        # if plot == True:
        #     pl.show()
        # else:
        #     pl.savefig(plot+'.pdf')
        #     pl.clf()





    # def MyCorrelate(self,x,y,Norm=True,MinAvg=True):
    #     if self.has_reps:
    #         return self.MyCorrelateWithRep(x,y,Norm=Norm,MinAvg=MinAvg)
    #     if MinAvg:
    #         x = x-x.mean()
    #         y = y-y.mean()
    #     listout = []
    #     for it in range(x.size):
    #         itval = 0.0
    #         for index in range(x.size):
    #             if it+index < y.size:
    #                 itval += x[index]*y[index+it]
    #         if Norm: itval = itval/(x.size-it)
    #         listout.append(itval)
    #     return np.array(listout)
    #
    # def MyCorrelateWithRep(self,x,y,Norm=True,MinAvg=True):
    #     if MinAvg:
    #         for ix,iy in zip(x,y):
    #             ix = ix-np.mean(ix)
    #             iy = iy-np.mean(iy)
    #     ## if replicas are of different length, they will revert to smallest list
    #     ## and cut the bigger replicas to the smaller length size.
    #     rlen = len(x)
    #     tlen = max(map(len,x))
    #     totlen = sum(map(len,x))
    #     listout = []
    #     for it in xrange(tlen):
    #         itval = 0.0
    #         for ix,iy in zip(x,y):
    #             for index in xrange(tlen):
    #                 if it+index < iy.size and index < ix.size:
    #                     itval = itval + ix[index]*iy[index+it]
    #         if Norm: itval = itval/(totlen-(rlen*it))
    #         listout.append(itval)
    #     return np.array(listout)
    #


    ### autocorrelation work taken from https://arxiv.org/pdf/hep-lat/0306017.pdf
    def autocorr(self,x,y):
        if self.has_reps:
            return self.autocorr_Reps(x,y)
        """
        http://stackoverflow.com/q/14297012/190597
        http://en.wikipedia.org/wiki/Autocorrelation#Estimation
        """
        n = len(x)
        # variance = x.var()
        x = x-x.mean()
        y = y-y.mean()
        r = np.correlate(x, y, mode = 'full')[-n:]
        # assert np.allclose(r, np.array([(x[:n-k]*x[-(n-k):]).sum() for k in xrange(n)]))
        result = r/np.arange(n, 0, -1)
        return result


    def autocorr_Reps(self,x,y):
        totn = sum(map(len,x))
        n_max = max(list(map(len,x)))
        result,Rlist = 0.,0.
        for ix,iy in zip(x,y):
            n = len(ix)
            # variance = ix.var()
            ix = ix-ix.mean()
            iy = iy-iy.mean()
            r = np.correlate(ix, iy, mode = 'full')[-n:]
            # assert np.allclose(r, np.array([(ix[:n-k]*ix[-(n-k):]).sum() for k in xrange(n)]))
            result += np.append(r,np.zeros(n_max-n))
            Rlist += np.append(np.ones(n),np.zeros(n_max-n))
        result = result / (totn-Rlist*np.arange(1,n_max+1))
        return result

    def gW(self,tauW,thisN):
        ## using auto fitting window method used from (52)
        for iW,it in enumerate(tauW):
            if iW == 0: continue
            val = np.exp(-iW/it)-it/np.sqrt(iW*thisN)
            if val < 0.0:
                return iW
        return float('NaN')

    def VarTau(self,tau,N):
        ## Using aproximate formula (42) from paper
        return [np.sqrt(4/float(N) * (iW + 0.5 - itau) * itau**2) for iW,itau in enumerate(tau)]

    def BiasCorrect(self,CfW,N):
        ## Bias corrections using (49)
        W = np.arange(len(CfW))
        return CfW*(1+((2*W+1)/float(N)))
        # ## Testing
        # return CfW


    #data = [ variable , monte time ]
    #fun(variables) output is value
    #funder(variables) output is partial derivatives (w.r.t variables)
    ## AllOut = False
    ## average, error[Woptimal], tau[Woptimal], tauerr[Woptimal]
    ## AllOut = True
    ## average, error(W), tau(W), tauerr(W), GFt(W), Wopt
    def uWerrMine(self,data='PreDef',fun='PreDef',Sparam='PreDef',WipeData=True,save_covar=True):
        self.GetFuns()
        if data == 'PreDef':
            if isinstance(self.values,str):
                print('data not imported for ',self.name)
            else:
                data = self.values
        if fun == 'PreDef':
            if isinstance(self.Fun,str):
                print('Function and derivative have not been defined for ',self.name)
            else:
                fun,funder = self.Fun,self.FunDer
        if Sparam=='PreDef': Sparam = self.Sparam
        avgdata = [data[icol].mean() for icol in data]


        glen = 0
        if isinstance(data.index,pa.MultiIndex):
            if self.has_reps:
                data_numpy = []
                for icol in data:
                    data_numpy.append([])
                    for istream,stream_data in data[icol].groupby(level='stream'):
                        # glen = max(glen,stream_data.size)
                        glen += stream_data.size
                        data_numpy[-1].append(stream_data[istream].values)
                data_numpy = np.array(data_numpy)
            else:
                data_numpy = np.swapaxes(data.values,0,1)
                glen = data_numpy.shape[1]
        else:
            data_numpy = np.swapaxes(data.values,0,1)
            glen = data_numpy.shape[1]


        ## (31) matrix of autocorrelations w.r.t ab= variables
        G_ab_t = []
        for adat in data_numpy:
            G_ab_t.append([ self.autocorr(adat,bdat) for bdat in data_numpy])
            # G_ab_t.append([ self.MyCorrelate(adat,bdat) for bdat in data_numpy])
        G_ab_t = np.array(G_ab_t)
        if save_covar:
            self.covar = []
            G_ab_0 = G_ab_t[:,:,0]
            for ia,Ga in enumerate(G_ab_t[:,:,0]):
                for ib,Gab in enumerate(Ga):
                    self.covar.append(Gab/np.sqrt(G_ab_0[ia,ia]*G_ab_0[ib,ib]))
            self.covar = np.array(self.covar).reshape(G_ab_0.shape)
            self.covar_NoNorm = G_ab_0/glen
        ## (33) alpha function derivates (w.r.t variables)
        f_a = np.array(funder(*avgdata))
        f_ab = []
        for if_a in f_a:
            f_ab.append([if_b*if_a for if_b in f_a])
        f_ab = np.array(f_ab)

        # print 'debug' , data.shape, G_ab_t.shape, f_ab.shape, (f_ab * G_ab_t[:,:,0]).shape,np.sum(f_ab * G_ab_t[:,:,0]).shape

        # print 'debug',self.name,
        # print f_ab * G_ab_t[:,:,5]
        # print
        ## (33)
        GFt = [np.sum(f_ab * G_ab) for G_ab in np.rollaxis(G_ab_t,-1) ]

        ## (35)
        CFW = np.array([GFt[0]]+[GFt[0] + 2*np.sum(GFt[1:W+1]) for W in range(1,len(GFt))])

        # Bias corrections (49)
        CFW = self.BiasCorrect(CFW,glen)


        ## (34)
        nuF = CFW[0]


        ## equation (41)
        tauint = CFW / (2*nuF)



        ## From Paper
        tau = []
        for it,itauint in enumerate(tauint):
            if itauint <= 0.5:
                tau.append(0.00000001)
            else:
                ## (51)
                tau.append(Sparam/logNA((2*itauint+1)/(2*itauint-1)))

        # ## From Matlab
        # tau = []
        # for iGFt in GFt:
        #     if iGFt <= 0.0:
        #         tau.append(0.0000001)
        #     else:
        #         tau.append(Sparam/logNA((iGFt+1)/iGFt))

        ## (52)
        Wopt = self.gW(np.array(tau),glen)

        ## (42)
        dtauint = self.VarTau(tauint,glen)


        ## average, error(W), tau(W), tauerr(W), GFt(W), Wopt
        ## average, error[Woptimal], tau[Woptimal], tauerr[Woptimal]
        self.Avg = fun(*avgdata)
        self.Avg_paras = avgdata
        # self.Avg = np.mean([fun(*ival) for ival in np.swapaxes(data,0,1)])
        try:
            self.Std = np.sqrt(np.abs(CFW[Wopt])/float(glen))
            self.Std_corr = np.sqrt(np.abs(CFW[0])/float(glen))
            self.Wopt = Wopt
            self.GFt = GFt
            self.tau = tauint[Wopt]
            self.tauerr = dtauint[Wopt]
        except:
            print('NaN for Autocorr')
            Wopt = -1
            self.Std = np.sqrt(np.abs(CFW[Wopt])/float(glen))
            self.Std_corr = np.sqrt(np.abs(CFW[0])/float(glen))
            self.Wopt = Wopt
            self.GFt = GFt
            self.tau = tauint[Wopt]
            self.tauerr = dtauint[Wopt]

        self.StdW = np.sqrt(np.abs(CFW)/float(glen))
        self.tauW = tauint
        self.tauerrW = dtauint
        if WipeData:
            self.RemoveVals()
        self.RemoveFuns()

    def GetFuns(self):
        if not hasattr(self,'Fun') or not hasattr(self,'FunDer'):
            self.Fun,self.FunDer = ReadFuns(self.Fun_name,self.FunDer_name)
            del self.Fun_name
            del self.FunDer_name


    def RemoveFuns(self):
        self.Fun_name = self.Fun.__name__
        self.FunDer_name = self.FunDer.__name__
        WriteFuns(self.Fun,self.FunDer_name)
        del self.Fun
        del self.FunDer

    def __str__(self):
        return 'Auto_'+self.name+': '+MakeValAndErr(self.Avg,self.Std)

    def __repr__(self):
        return '<Auto_'+self.name+'_'+MakeValAndErr(self.Avg,self.Std,Dec=1)+'>'





class BlockedAutoCorr(object):
    '''
        Uses AutoCorrelate class above and repeats it for different blocked amounts until
        it hits a tau int threshold
    '''
    def __init__(self,  name='',Fun='NotDef', data='None',Sparam=defSparam,save_covar=False):
        self.cfg_list = data
        self.name = name
        self.Fun = Fun
        self.Sparam = Sparam
        self.save_covar = save_covar
        self.noblock_auto = AutoCorrelate(name=name,Fun=Fun,data=data,Sparam=Sparam,save_covar=save_covar,WipeData=False)
        self.auto_object = deepcopy(self.noblock_auto)


    ## block_cutoff is a percent of the total length!
    def RemoveAuto(self,tau_cutoff = 1,block_cutoff= 10):
        ## TODO use tauerr?
        n_block = 1
        # this_cfglist = self.noblock_auto.values
        this_cfglist = self.cfg_list
        blockc_per = int(block_cutoff*len(self.cfg_list)/100)
        # print('DEBUG')
        # import matplotlib.pyplot as pl
        while self.auto_object.tau-self.auto_object.tauerr > tau_cutoff and blockc_per > n_block:
            n_block += 1
            this_cfglist = self.noblock_auto.BlockCfgs(n_block)
            # if n_block//30 == n_block/30:
            #     print(n_block)
            #     print(this_cfglist.values[:,0])
            #     print()
            #     pl.plot(this_cfglist.values[:,0])
            #     pl.title(str(n_block))
            #     pl.savefig('/home/jackdra/LQCD/Scripts/Python_Analysis/TestGraphs/testsin'+str(n_block)+'.pdf')
            #     pl.clf()
            non_dup_cfg = pa.DataFrame(np.abs(this_cfglist.values[::n_block,:]),columns=this_cfglist.columns)
            self.auto_object = AutoCorrelate(name=self.name,Fun=self.Fun,data=non_dup_cfg,
                                             Sparam=self.Sparam,save_covar=self.save_covar)
            # print(tau_cutoff)
            # print(MakeValAndErr(self.auto_object.tau,self.auto_object.tauerr))
            # print(self.auto_object.tau > tau_cutoff)
        # print(MakeValAndErr(self.auto_object.tau-self.auto_object.tauerr,self.auto_object.tauerr))
        # print()
        if n_block >= block_cutoff:
            print()
            print('Warning, blocking to reduce autocorrelaiton hit cutoff of n_block percent='+str(block_cutoff)+'%')
            print('class: '+self.name)
            print('tau_int was '+MakeValAndErr(self.auto_object.tau,self.auto_object.tauerr))
            print()
        return this_cfglist,n_block

    ## block_cutoff is a percent of the total length!
    def GetBlockedBS(self,thisnboot=nboot,rand_list=None,tau_cutoff=0.5,block_cutoff=10,WipeData=True,noblock_comp=False):
        this_cfg,n_block = self.RemoveAuto(tau_cutoff=tau_cutoff,block_cutoff=block_cutoff)
        this_name = self.name + '_nblock'+str(n_block)
        self.n_block = n_block
        self.blocked_cfglist = this_cfg
        if 'sin' in self.name:
            print('DEBUG')
            print('plotting blocked distribution')
            non_dup_cfg = pa.DataFrame(this_cfg.values[::n_block,:],columns=this_cfg.columns)
            import matplotlib.pyplot as pl
            pl.plot(non_dup_cfg.values[:,0])
            pl.savefig('/home/jackdra/LQCD/Scripts/Python_Analysis/TestGraphs/testsin_blocked.pdf')
            pl.clf()
        out_bs = BootStrap(thisnboot=thisnboot,name=this_name,cfgvals=this_cfg,thisDelVal=WipeData,rand_list=rand_list)
        if WipeData:
            self.auto_object.RemoveVals()
            self.noblock_auto.RemoveVals()
        if noblock_comp:
            out_bs_nb = BootStrap(thisnboot=thisnboot,name=self.name,cfgvals=self.cfg_list,thisDelVal=WipeData,rand_list=rand_list)
            return out_bs,out_bs_nb
        else:
            return out_bs


def TestBlocking(tau_cutoff=0.5,block_cutoff=10):
    def thisFun(*x):
        return x[0]
    def thisDer(*x):
        return [1]
    this_size = 10000
    this_nboot = 1000
    nperiod = 100
    values = np.random.uniform(size=this_size)
    values3 = np.random.normal(loc=1,scale=1,size=this_size)
    values2 = np.sin(np.arange(this_size)*(nperiod*2*np.pi)/this_size)+values*20
    import matplotlib.pyplot as pl
    pl.plot(values2)
    pl.savefig('/home/jackdra/LQCD/Scripts/Python_Analysis/TestGraphs/testsin.pdf')
    pl.clf()
    val_df = pa.DataFrame()
    val_df['one'] = pa.Series(values)
    val_df['two'] = pa.Series(values2)
    val_df['three'] = pa.Series(values3)

    testdata = BlockedAutoCorr(Fun=[thisFun,thisDer],name='test_bootstrap_uniform',data=val_df[['one']])
    testdata2 = BlockedAutoCorr(Fun=[thisFun,thisDer],name='test_bootstrap_sin_normal',data=val_df[['two']])
    testdata3 = BlockedAutoCorr(Fun=[thisFun,thisDer],name='test_bootstrap_normal',data=val_df[['three']])
    bootdata,nbb = testdata.GetBlockedBS(thisnboot=this_nboot,noblock_comp=True,tau_cutoff=tau_cutoff,block_cutoff=block_cutoff)
    bootdata2,nbb2 = testdata2.GetBlockedBS(thisnboot=this_nboot,noblock_comp=True,tau_cutoff=tau_cutoff,block_cutoff=block_cutoff)
    bootdata3,nbb3 = testdata3.GetBlockedBS(thisnboot=this_nboot,noblock_comp=True,tau_cutoff=tau_cutoff,block_cutoff=block_cutoff)
    return testdata,testdata2,testdata3,bootdata,bootdata2,bootdata3,nbb,nbb2,nbb3

def TestAuto():
    def thisFun(*x):
        return x[0]
    def thisDer(*x):
        return [1]

    const = 100
    this_size = 20000
    values = np.random.uniform(size=this_size)

    values2 = np.arange(this_size)/this_size

    values3 = np.random.normal(loc=0.5,scale=0.25,size=this_size)
    val_df = pa.DataFrame()

    # tuple_list = []
    # for ii in range(100):
    #     tuple_list.append(('-1-',ii))
    # for ii in range(400):
    #     tuple_list.append(('-2-',ii))
    # for ii in range(1000):
    #     tuple_list.append(('-3-',ii))
    # for ii in range(500):
    #     tuple_list.append(('-4-',ii))

    tuple_list = []
    for ii in range(this_size//2):
        tuple_list.append(('-1-',ii))
    for ii in range(this_size//2):
        tuple_list.append(('-2-',ii))
    #
    # tuple_list = []
    # for ii in range(this_size):
    #     tuple_list.append(('-1-',ii))

    indicies = pa.MultiIndex.from_tuples(tuple_list,names=['stream','configs'])
    # indicies = range(this_size)
    # val_df = pa.DataFrame()
    # val_df['one'] = pa.Series(values,index=indicies)
    # val_df['two'] = pa.Series(values2,index=indicies)
    # val_df['three'] = pa.Series(values3,index=indicies)
    # def RatFun(one,two,three):
    #     return const*one/(two*three)
    #
    # def RatFunDer(one,two,three):
    #     return [const/(two*three),-const*one/(three*two**2),-const*one/(two*three**2)]

    val_df = pa.DataFrame()
    val_df['one'] = pa.Series(values,index=indicies)
    val_df['two'] = pa.Series(values2,index=indicies)
    val_df['three'] = pa.Series(values3,index=indicies)
    def RatFun(one,two):
        return const*one*two

    def RatFunDer(one,two):
        return [const*two,const*one]

    testdata = AutoCorrelate(Fun=[thisFun,thisDer],name='test_bootstrap_uniform',data=val_df[['one']])
    testdata2 = AutoCorrelate(Fun=[thisFun,thisDer],name='test_bootstrap_arange',data=val_df[['two']])
    testdata3 = AutoCorrelate(Fun=[thisFun,thisDer],name='test_bootstrap_normal',data=val_df[['three']])
    testdatarat = AutoCorrelate(Fun=[RatFun,RatFunDer],name='test_auto_ratio',data=val_df[['one','two']])

    this_info = pa.Series()
    this_info['save_file'] = './TestGraphs/test_Wopt.pdf'
    this_info['title'] = 'Test Auto Graph'
    # this_info['xlims'] = [0,10]
    # this_info['ylims'] = [0,15]
    import PlotData as jpl
    data_plot = jpl.Plotting(plot_info=this_info)
    data_plot = testdata.PlotWopt(data_plot)
    data_plot = testdata2.PlotWopt(data_plot)
    data_plot = testdata3.PlotWopt(data_plot)
    data_plot = testdatarat.PlotWopt(data_plot)
    # data_plot.LoadPickle(DefWipe=False)
    data_plot.PrintData()
    data_plot.PlotAll()

    this_info = pa.Series()
    this_info['save_file'] = './TestGraphs/test_Auto.pdf'
    this_info['title'] = 'Test Auto Graph'
    # this_info['xlims'] = [0,10]
    # this_info['ylims'] = [0,15]
    import PlotData as jpl
    data_plot = jpl.Plotting(plot_info=this_info)
    data_plot = testdata.PlotTauInt(data_plot)
    data_plot = testdata2.PlotTauInt(data_plot)
    data_plot = testdata3.PlotTauInt(data_plot)
    data_plot = testdatarat.PlotTauInt(data_plot)
    # data_plot.LoadPickle(DefWipe=False)
    data_plot.PrintData()
    data_plot.PlotAll()

    return testdata,testdata2,testdata3,testdatarat


if __name__ == '__main__':
    testdata1,testdata2,testdata3,testdatarat = TestAuto()
    # %matplotlib inline
    print('Result is in testdata1,2,3')
