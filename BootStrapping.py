#!/usr/bin/env python


# IMPORT THIS FIRST, put in base import file
import numpy as np
import pandas as pa
import warnings
from XmlFormatting import MakeValAndErr
from MiscFuns import human_format,logNA
from Params import nboot,normed_histograms,percentile
from Params import Debug_Mode
# import matplotlib
# matplotlib.use('Agg')  # Must be before importing matplotlib.pyplot or pylab!

"""
How to use:

BootStrap(nboot)
creates class with nbootstraps


.ImportData(data)
Imports data into bootstrap

.MakeBS()
Creates bootstrapped data

.Stats()
Creates .Avg and .Std

.BootHistogram()
Creates histogram (only plot command, do label and savefig after).


One thing to note, Casting to differnet type will perminatly cast the bootstrap to that type (even if passed into function and combined wiht other bs).
"""

def BlockCfgs(cfg_list,n_block):
    if not isinstance(cfg_list,np.ndarray):
        cfg_list = np.array(cfg_list)
    extra_len = (len(cfg_list)//n_block)*n_block
    if extra_len == 0:
        b_vals = np.mean([cfg_list[iblock::n_block] for iblock in range(n_block)],axis=0)
        bc_vals = np.empty((cfg_list.size,),dtype=b_vals.dtype)
        for i in range(n_block):
            bc_vals[i::n_block] = b_vals
    else:
        b_vals = np.mean([cfg_list[iblock:extra_len:n_block] for iblock in range(n_block)],axis=0)
        bc_vals = np.empty((cfg_list.size,),dtype=b_vals.dtype)
        for i in range(n_block):
            bc_vals[i:extra_len:n_block] = b_vals
        bc_vals[extra_len:] = np.mean(cfg_list[extra_len:])
    return bc_vals

def CombineBS(this_series,name='',resample=False):
    if not isinstance(this_series,pa.Series):
        raise IOError('please pass in Series to be combined\n',type(this_series))
    out_bs_list = []
    for ival in this_series.values:
        out_bs_list += list(ival.GetBootVals())
    out_boot = BootStrap(thisnboot=len(out_bs_list), name=name,bootvals=out_bs_list)
    if resample is not False:
        out_boot,dump = out_boot.ResampleBS(new_nboot=resample)
    return out_boot

def FlattenBootstrapDF(this_df,this_nboot=nboot,reset_index=True):
    if not isinstance(this_df,pa.DataFrame):
        raise IOError('please pass in DataFrame to be flattened')
    if reset_index:
        this_df.reset_index(inplace=True)
    data_list,boot_index = [[] for ival in this_df.columns],[]
    boot_range = list(range(this_nboot))
    for irow,row_data in this_df.iterrows():
        if all([not isinstance(icol,BootStrap) for icol in row_data.values]):
            boot_index.append(float('NaN'))
            for icol,col_data in enumerate(row_data.values):
                data_list[icol].append(col_data)
        else:
            boot_index += boot_range
            for icol,col_data in enumerate(row_data.values):
                if isinstance(col_data,BootStrap):
                    data_list[icol] = data_list[icol] + list(col_data.GetBootVals())
                else:
                    data_list[icol] = data_list[icol] + [col_data for _ in boot_range]
    if len(data_list) == 0:
        return this_df
    elif len(data_list[0]) == 0:
        print('No Data Fround to flatten in FlattenBootstrapDF')
        return this_df
    else:
        out_df = pa.DataFrame(columns=this_df.columns)
        out_df['boot_index'] = pa.Series(boot_index)
        for icol,col_data in zip(this_df.columns,data_list):
            out_df[icol] = pa.Series(col_data)
        return out_df



startseed = 1234


# def MeanOfBoots(bootlist):
#     return np.mean(bootlist)
#     # return np.average(bootlist)
#     # return np.sum(np.array(bootlist)/float(len(bootlist)))


class BootStrap(object):
    """ Class that a single boostrapped quantity Xi for i cfgvals
        And Constructs boostrapped cfgvals Xb for b bootstraped values

    """

    parentlist = ['SetOfCorrs.SetOfTwoPt','TwoPtCorrelators.TwoPointCorr','TwoPtCorrelators.TwoPointCorr','TwoPtCorrelators.NNQCorr']

    def __init__(self, thisnboot=nboot, name='',cfgvals='None',bootvals='None',thisDelVal=True,rand_list=None,
                 n_block = 1):
        self.nboot = thisnboot  # number of bootstrap samples

        self.nval = 0.0  # contains length of cfgvalsn
        # These are for legasy versions, that did not use pandas.
        # will eventually depreciate
        self.bootvals = pa.Series(name='bootvals')
        self.cfgvals = pa.Series(name='cfgvals')
        if n_block > 1:
            self.blocked_bootvals = {}
            self.blocked_cfgvals = {}
            self.n_block = n_block
        self.block_sel = 'block1'

        self.name = name  # Give it a label if you want, might help with plotting?
        self.Avg = np.float64(np.nan)
        self.Std = np.float64(np.nan)
        self.Median = np.float64(np.nan)
        self.MedStd = np.float64(np.nan)
        self.MedStd_Low = np.float64(np.nan)
        self.MedStd_High = np.float64(np.nan)
        self.current = 0
        self.rand_list = rand_list

        if isinstance(bootvals,(list,tuple,np.ndarray,pa.Series)):
            self.ImportBS(bootvals)
        else:
            if isinstance(cfgvals,pa.DataFrame):
                cfgvals = cfgvals.iloc[:,0]
            if isinstance(cfgvals,(list,tuple,np.ndarray,pa.Series)):
                self.ImportData(cfgvals)
                self.MakeBS(DelVal=thisDelVal)

    def DivByLast(self):
        ## hack function for if you append the bootstrap samples with
        ## the tree level result
        this_bootvals = self.bootvals[:-1]/self.bootvals[-1]
        if not hasattr(self,'n_block'):
            self.n_block = 1
        return BootStrap(thisnboot=self.nboot-1,
                         name=self.name+'_div_tree',
                         bootvals=this_bootvals,rand_list=self.rand_list,
                         n_block=self.n_block)


    def GetBootVals(self):
        return self.bootvals.values

    def GetCfgVals(self):
        return self.cfgvals.values

    def CheckBoot(self):
        if self.bootvals.size == 0:
            raise EnvironmentError(self.name + ' has does not contain boostrapped values')
        if float('NaN') in self.bootvals:
            warnings.warn('Bootstrap has NaN present, ignorning NaNs')


    def CheckBoot_Blocked(self):
        if hasattr(self,'n_block') and self.n_block > 1:
            if len(self.blocked_bootvals) == 0:
                raise EnvironmentError(self.name + ' has does not contain blocked boostrapped values')
            for iblock in self.blocked_bootvals.values():
                if float('NaN') in iblock:
                    warnings.warn('Bootstrap has NaN present, ignorning NaNs')

    def CheckCfgVals(self):
        if self.cfgvals.size == 0:
            raise EnvironmentError(self.name + ' has does not contain configuration values')

    def CheckCfgVals_Blocked(self):
        if hasattr(self,'n_block') and self.n_block > 1:
            if len(self.blocked_cfgvals) == 0:
                raise EnvironmentError(self.name + ' has does not contain blocked configuration values')

    def GetNaNIndex(self):
        return pa.isnull(self.bootvals).nonzero()[0]

    def GetNanIndex_cfgs(self):
        return pa.isnull(self.cfgvals).nonzero()[0]


    def Stats(self):
        # print self.bootvals
        self.CheckBoot()
        self.Avg = self.bootvals.mean()
        self.Std = self.bootvals.std()
        self.Median = self.bootvals.median()
        self.MedStd_Low = self.bootvals.quantile(percentile)
        self.MedStd_High = self.bootvals.quantile(1-percentile)
        self.MedStd = (self.MedStd_High-self.MedStd_Low )/2
        self.MedStd_Up = self.MedStd_High-self.Median
        self.MedStd_Down = self.Median - self.MedStd_Low
        self.Min = self.bootvals.min()
        self.Max = self.bootvals.max()
        return self

    def RemoveVals(self):
        self.cfgvals = pa.Series(name='cfgvals')
        # self.cfgvals = None
        self.nval = 0.0
        if hasattr(self,'n_block') and self.n_block > 1:
            self.blocked_cfgvals = {}

    def RemoveBoots(self):
        self.bootvals = pa.Series(name='bootvals')
        if hasattr(self,'n_block') and self.n_block > 1:
            self.blocked_bootvals = {}
        # self.bootvals = None

    def BlockCfgs(self):
        if hasattr(self,'n_block') and self.n_block > 1:
            self.blocked_cfgvals['block1'] = self.cfgvals
            for iblock in range(2,self.n_block+1):
                self.blocked_cfgvals['block'+str(iblock)] = pa.Series(self.cfgvals.values[::iblock])
                for iib in range(1,iblock):
                    self.blocked_cfgvals['block'+str(iblock)] += pa.Series(self.cfgvals.values[iib::iblock])
                self.blocked_cfgvals['block'+str(iblock)] = self.blocked_cfgvals['block'+str(iblock)]/iblock

    def ImportData(self,cfgvals):
        if isinstance(cfgvals,pa.Series):
            self.cfgvals = cfgvals
        else:
            self.cfgvals = pa.Series(cfgvals,name='cfgvals')
        self.nval = self.cfgvals.size
        if hasattr(self,'n_block') and self.n_block > 1:
            self.BlockCfgs()


    def ImportBS_Blocked(self,inputboot,this_block):
        if isinstance(this_block,(int,float)):
            this_block = 'block'+str(int(this_block))
        if isinstance(inputboot,pa.Series):
            if not isinstance(inputboot.dtype,(int,float,complex,np.dtype)) :
                print(inputboot.dtype)
                raise EnvironmentError('Class got passed to bootvals, or multiple keys exist')
            if len(inputboot) == self.nboot:
                self.blocked_bootvals[this_block] = inputboot
            else:
                raise IOError('nboots dont match up')
        else:
            if not isinstance(inputboot[0],(int,float,complex)) :
                print(type(inputboot[0]))
                raise EnvironmentError('Class got passed to bootvals'+str(inputboot[0]))
            if len(inputboot) == self.nboot:
                self.blocked_bootvals[this_block] = pa.Series(inputboot,name='bootvals')
            else:
                raise IOError('nboots dont match up')


    def ImportBS(self,inputboot,DelVal=True):
        if DelVal: self.RemoveVals()
        if isinstance(inputboot,pa.Series):
            if not isinstance(inputboot.dtype,(int,float,complex,np.dtype)) :
                print(inputboot.dtype)
                raise EnvironmentError('Class got passed to bootvals, or multiple keys exist')
            if len(inputboot) == self.nboot:
                self.bootvals = inputboot
                self.bootvals.name = 'bootvals'
            else:
                raise IOError('nboots dont match up')
        else:
            if not isinstance(inputboot[0],(int,float,complex)) :
                print(type(inputboot[0]))
                raise EnvironmentError('Class got passed to bootvals'+str(inputboot[0]))
            if len(inputboot) == self.nboot:
                self.bootvals = pa.Series(inputboot,name='bootvals')
            else:
                raise IOError('nboots dont match up')
        self.Stats()


    def Get_Rand_List(self,wipe_rlist=True):
        this_rlist = self.rand_list
        if wipe_rlist:
            self.rand_list = None
        return this_rlist

    def MakeBS_Blocked(self):
        self.CheckCfgVals_Blocked()
        for iblock,block_cfgs in self.blocked_cfgvals.items():
            this_nval = len(block_cfgs)
            myseed=startseed*this_nval//self.nboot
            # self.Avg = np.average(self.values)
            # if not isinstance(self.rand_list,np.ndarray):
            ### TODO pass random list to blocked
            np.random.seed(myseed)
            this_rand_list = np.random.randint(0,high=this_nval,size=(self.nboot,this_nval))
            tempbootvals = [block_cfgs.values[iboot].mean() for iboot in this_rand_list]
            ## above is faster than below
            # tempbootvals = []
            # for iboot in xrange(self.nboot):
            #     # print 'DEBUG',iboot
            #     # print self.data.cfgvals.sample(self.nval, random_state=myseed+iboot).head()
            #     # print self.data.cfgvals.sample(self.nval, random_state=myseed+iboot).mean()
            #     # print
            #     ## TODO, check if differnet streams need to be separated out.
            #     tempbootvals.append(self.cfgvals.sample(self.nval, replace=True,random_state=myseed+iboot).mean())
            self.ImportBS_Blocked(tempbootvals,iblock)
        # self.Stats()

    def Select_Block(self,this_block):
        if hasattr(self,'n_block') and self.n_block > 1:
            if isinstance(this_block,(int,float)):
                this_block = 'block'+str(int(this_block))
            cfg_changed,boot_changed = False,False
            if this_block in self.blocked_cfgvals.keys():
                self.cfgvals = self.blocked_cfgvals[this_block]
                cfg_changed = True
            if this_block in self.blocked_bootvals.keys():
                self.bootvals = self.blocked_bootvals[this_block]
                boot_changed = True
            if boot_changed or cfg_changed:
                self.block_sel = this_block
            else:
                print('select block failed for '+this_block)

    def GetSingleBlock(self,this_block):
        self.Select_Block(this_block)
        ##default n_block = 1
        return BootStrap(thisnboot=self.nboot, name=self.name,bootvals=self.bootvals)

    def MakeBS(self,DelVal=True):
        self.CheckCfgVals()
        myseed=startseed*self.nval//self.nboot
        # self.Avg = np.average(self.values)
        if not isinstance(self.rand_list,np.ndarray):
            np.random.seed(myseed)
            self.rand_list = np.random.randint(0,high=self.nval,size=(self.nboot,self.nval))
        tempbootvals = [self.cfgvals.values[iboot].mean() for iboot in self.rand_list]
        ## above is faster than below
        # tempbootvals = []
        # for iboot in xrange(self.nboot):
        #     # print 'DEBUG',iboot
        #     # print self.data.cfgvals.sample(self.nval, random_state=myseed+iboot).head()
        #     # print self.data.cfgvals.sample(self.nval, random_state=myseed+iboot).mean()
        #     # print
        #     ## TODO, check if differnet streams need to be separated out.
        #     tempbootvals.append(self.cfgvals.sample(self.nval, replace=True,random_state=myseed+iboot).mean())
        if hasattr(self,'n_block') and self.n_block > 1:
            self.MakeBS_Blocked()
        self.ImportBS(tempbootvals,DelVal=DelVal)
        self.Stats()




    # Operator overloading

    # Comparisons
    def MakeValAndErr(self,Dec=1,MagThresh=3,latex=False):
        self.Stats()
        return MakeValAndErr(self.Avg,self.Std,Dec=Dec,MagThresh=MagThresh,latex=latex)

    def __str__(self):
        return 'BS_'+self.name+': '+MakeValAndErr(self.Avg,self.Std)

    def __repr__(self):
        return '<BS_'+self.name+'_'+MakeValAndErr(self.Avg,self.Std,Dec=1)+'>'

    def __lt__(self,bs2):
        # self.Stats()
        if isinstance(bs2,BootStrap):
            return self.Avg+self.Std < bs2.Avg-bs2.Std
        else:
            return self.Avg+self.Std < bs2

    def __le__(self,bs2):
        # self.Stats()
        if isinstance(bs2,BootStrap):
            return self.Avg-self.Std <= bs2.Avg+bs2.Std
        else:
            return self.Avg-self.Std <= bs2


    def __gt__(self,bs2):
        # self.Stats()
        if isinstance(bs2,BootStrap):
            return self.Avg-self.Std > bs2.Avg+bs2.Std
        else:
            return self.Avg-self.Std > bs2

    def __ge__(self,bs2):
        # self.Stats()
        if isinstance(bs2,BootStrap):
            return self.Avg+self.Std >= bs2.Avg-bs2.Std
        else:
            return self.Avg+self.Std >= bs2


    def __eq__(self,bs2):
        return (self <= bs2) and (self >= bs2)


    def __ne__(self,bs2):
        return (self < bs2) or (self > bs2)

    ## If you uncomment this, it breaks the program,
    ## interestingly, it lets you cast a bootstrap class to a numpy array.
    # def __len__(self):
    #     return self.nboot
    # Numerical operatoions

    def Log(self):
        result = BootStrap(self.nboot,name='log('+str(self)+')')
        result.cfgvals = logNA(self.cfgvals)
        result.bootvals = logNA(self.bootvals)
        if hasattr(self,'n_block') and self.n_block > 1:
            result.block_cfgvals = {}
            result.block_bootvals = {}
            for ikey,idata in self.block_cfgvals.items():
                result.block_cfgvals = logNA(self.cfgvals)
            for ikey,idata in self.block_cfgvals.items():
                result.block_bootvals = logNA(self.bootvals)

        return result

    def Exp(self):
        result = BootStrap(self.nboot,name='e^('+str(self)+')')
        result.cfgvals = np.exp(self.cfgvals)
        result.bootvals = np.exp(self.bootvals)
        if hasattr(self,'n_block') and self.n_block > 1:
            result.block_cfgvals = {}
            result.block_bootvals = {}
            for ikey,idata in self.block_cfgvals.items():
                result.block_cfgvals = np.exp(self.cfgvals)
            for ikey,idata in self.block_cfgvals.items():
                result.block_bootvals = np.exp(self.bootvals)
        return result

    def sin(self):
        result = BootStrap(self.nboot,name='sin('+str(self)+')')
        result.cfgvals = np.sin(self.cfgvals)
        result.bootvals = np.sin(self.bootvals)
        if hasattr(self,'n_block') and self.n_block > 1:
            result.block_cfgvals = {}
            result.block_bootvals = {}
            for ikey,idata in self.block_cfgvals.items():
                result.block_cfgvals = np.sin(self.cfgvals)
            for ikey,idata in self.block_cfgvals.items():
                result.block_bootvals = np.sin(self.bootvals)
        return result

    def sinh(self):
        result = BootStrap(self.nboot,name='sinh('+str(self)+')')
        result.cfgvals = np.sinh(self.cfgvals)
        result.bootvals = np.sinh(self.bootvals)
        if hasattr(self,'n_block') and self.n_block > 1:
            result.block_cfgvals = {}
            result.block_bootvals = {}
            for ikey,idata in self.block_cfgvals.items():
                result.block_cfgvals = np.sinh(self.cfgvals)
            for ikey,idata in self.block_cfgvals.items():
                result.block_bootvals = np.sinh(self.bootvals)
        return result


    def arcsinh(self):
        result = BootStrap(self.nboot,name='arcsinh('+str(self)+')')
        result.cfgvals = np.arcsinh(self.cfgvals)
        result.bootvals = np.arcsinh(self.bootvals)
        if hasattr(self,'n_block') and self.n_block > 1:
            result.block_cfgvals = {}
            result.block_bootvals = {}
            for ikey,idata in self.block_cfgvals.items():
                result.block_cfgvals = np.arcsinh(self.cfgvals)
            for ikey,idata in self.block_cfgvals.items():
                result.block_bootvals = np.arcsinh(self.bootvals)
        return result

    def cos(self):
        result = BootStrap(self.nboot,name='cos('+str(self)+')')
        result.cfgvals = np.cos(self.cfgvals)
        result.bootvals = np.cos(self.bootvals)
        if hasattr(self,'n_block') and self.n_block > 1:
            result.block_cfgvals = {}
            result.block_bootvals = {}
            for ikey,idata in self.block_cfgvals.items():
                result.block_cfgvals = np.cos(self.cfgvals)
            for ikey,idata in self.block_cfgvals.items():
                result.block_bootvals = np.cos(self.bootvals)
        return result

    def Sqrt(self):
        result = BootStrap(self.nboot,name='sqrt('+str(self)+')')
        result.cfgvals = np.sqrt(self.cfgvals)
        result.bootvals = np.sqrt(self.bootvals)
        if hasattr(self,'n_block') and self.n_block > 1:
            result.block_cfgvals = {}
            result.block_bootvals = {}
            for ikey,idata in self.block_cfgvals.items():
                result.block_cfgvals = np.sqrt(self.cfgvals)
            for ikey,idata in self.block_cfgvals.items():
                result.block_bootvals = np.sqrt(self.bootvals)
        return result


    def ResampleBS(self,new_nboot=nboot,rand_list=None):
        if new_nboot == self.nboot:
            return self
        if rand_list is None:
            rand_list = np.random.choice(self.nboot,size=new_nboot,replace=True)
        new_bootvals = self.bootvals.iloc[rand_list].reset_index()['bootvals']
        # new_bootvals = self.bootvals.iloc[rand_list].reset_index()['bootvals'].values
        new_bs = BootStrap(new_nboot,name=self.name+'_RS'+human_format(new_nboot),bootvals=new_bootvals)
        if hasattr(self,'n_block') and self.n_block > 1:
            new_bs.block_bootvals = {}
            for ikey,idata in self.block_bootvals.items():
                new_bs.block_bootvals[ikey] = idata.iloc[rand_list]
        return new_bs,rand_list


    def corr(self,bs2):
        if not isinstance(bs2,BootStrap):
            raise IOError('pass bootstrap variable into correlation')
        else:
            return self.bootvals.corr(bs2.bootvals)

    def cfg_corr(self,bs2):
        if not isinstance(bs2,BootStrap):
            raise IOError('pass bootstrap variable into correlation')
        else:
            return self.cfgvals.corr(bs2.cfgvals)


    def __add__(self,bs2):
        if type(bs2) in BootStrap.parentlist:
            return bs2.__radd__(self)
        else:
            result = BootStrap(self.nboot,name='+'.join([str(self).replace('BS_',''),str(bs2).replace('BS_','')]))
            if isinstance(bs2,BootStrap):
                result.cfgvals = self.cfgvals + bs2.cfgvals
                result.bootvals = self.bootvals + bs2.bootvals
                if hasattr(self,'n_block') and self.n_block > 1:
                    result.blocked_cfgvals = {}
                    result.blocked_bootvals = {}
                    for (ikey,iself),ibs2 in zip(self.blocked_cfgvals.items(),bs2.blocked_cfgvals.values()):
                        result.blocked_cfgvals[ikey] = iself + ibs2
                    for (ikey,iself),ibs2 in zip(self.blocked_bootvals.items(),bs2.blocked_bootvals.values()):
                        result.blocked_bootvals[ikey] = iself + ibs2
            else:
                try:
                    result.cfgvals = self.cfgvals + bs2
                    result.bootvals = self.bootvals + bs2
                    if hasattr(self,'n_block') and self.n_block > 1:
                        result.blocked_cfgvals = {}
                        result.blocked_bootvals = {}
                        for ikey,iself in self.blocked_cfgvals.items():
                            result.blocked_cfgvals[ikey] = iself + bs2
                        for ikey,iself in self.blocked_bootvals.items():
                            result.blocked_bootvals[ikey] = iself + bs2
                except:
                    print(type(bs2))
                    raise EnvironmentError('Invalid value to combine with BootStrap class')
            return result


    def __sub__(self,bs2):
        if type(bs2) in BootStrap.parentlist:
            return bs2.__rsub__(self)
        else:
            result = BootStrap(self.nboot,name='-'.join([str(self).replace('BS_',''),str(bs2).replace('BS_','')]))
            if isinstance(bs2,BootStrap):
                result.cfgvals = self.cfgvals - bs2.cfgvals
                result.bootvals = self.bootvals - bs2.bootvals
                if hasattr(self,'n_block') and self.n_block > 1:
                    result.blocked_cfgvals = {}
                    result.blocked_bootvals = {}
                    for (ikey,iself),ibs2 in zip(self.blocked_cfgvals.items(),bs2.blocked_cfgvals.values()):
                        result.blocked_cfgvals[ikey] = iself - ibs2
                    for (ikey,iself),ibs2 in zip(self.blocked_bootvals.items(),bs2.blocked_bootvals.values()):
                        result.blocked_bootvals[ikey] = iself - ibs2
            else:
                try:
                    result.cfgvals = self.cfgvals - bs2
                    result.bootvals = self.bootvals - bs2
                    if hasattr(self,'n_block') and self.n_block > 1:
                        result.blocked_cfgvals = {}
                        result.blocked_bootvals = {}
                        for ikey,iself in self.blocked_cfgvals.items():
                            result.blocked_cfgvals[ikey] = iself - bs2
                        for ikey,iself in self.blocked_bootvals.items():
                            result.blocked_bootvals[ikey] = iself - bs2
                except:
                    print(type(bs2))
                    raise EnvironmentError('Invalid value to combine with BootStrap class')
            return result


    def __mul__(self,bs2):
        if type(bs2) in BootStrap.parentlist:
            return bs2.__rmul__(self)
        else:
            result = BootStrap(self.nboot,name='*'.join([str(self).replace('BS_',''),str(bs2).replace('BS_','')]))
            if isinstance(bs2,BootStrap):
                result.cfgvals = self.cfgvals * bs2.cfgvals
                result.bootvals = self.bootvals * bs2.bootvals
                if hasattr(self,'n_block') and self.n_block > 1:
                    result.blocked_cfgvals = {}
                    result.blocked_bootvals = {}
                    for (ikey,iself),ibs2 in zip(self.blocked_cfgvals.items(),bs2.blocked_cfgvals.values()):
                        result.blocked_cfgvals[ikey] = iself * ibs2
                    for (ikey,iself),ibs2 in zip(self.blocked_bootvals.items(),bs2.blocked_bootvals.values()):
                        result.blocked_bootvals[ikey] = iself * ibs2
            else:
                try:
                    result.cfgvals = self.cfgvals * bs2
                    result.bootvals = self.bootvals * bs2
                    if hasattr(self,'n_block') and self.n_block > 1:
                        result.blocked_cfgvals = {}
                        result.blocked_bootvals = {}
                        for ikey,iself in self.blocked_cfgvals.items():
                            result.blocked_cfgvals[ikey] = iself * bs2
                        for ikey,iself in self.blocked_bootvals.items():
                            result.blocked_bootvals[ikey] = iself * bs2
                except:
                    print(type(bs2))
                    raise EnvironmentError('Invalid value to combine with BootStrap class')
            return result

    def __truediv__(self,bs2):
        return self.__div__(bs2)

    def __floordiv__(self,bs2):
        return self.__div__(bs2)

    def __div__(self,bs2):
        if type(bs2) in BootStrap.parentlist:
            return bs2.__rdiv__(self)
        else:
            result = BootStrap(self.nboot,name='/'.join([str(self).replace('BS_',''),str(bs2).replace('BS_','')]))
            if isinstance(bs2,BootStrap):
                result.cfgvals = self.cfgvals / bs2.cfgvals
                result.bootvals = self.bootvals / bs2.bootvals
                if hasattr(self,'n_block') and self.n_block > 1:
                    result.blocked_cfgvals = {}
                    result.blocked_bootvals = {}
                    for (ikey,iself),ibs2 in zip(self.blocked_cfgvals.items(),bs2.blocked_cfgvals.values()):
                        result.blocked_cfgvals[ikey] = iself / ibs2
                    for (ikey,iself),ibs2 in zip(self.blocked_bootvals.items(),bs2.blocked_bootvals.values()):
                        result.blocked_bootvals[ikey] = iself / ibs2
            else:
                try:
                    result.cfgvals = self.cfgvals / bs2
                    result.bootvals = self.bootvals / bs2
                    if hasattr(self,'n_block') and self.n_block > 1:
                        result.blocked_cfgvals = {}
                        result.blocked_bootvals = {}
                        for ikey,iself in self.blocked_cfgvals.items():
                            result.blocked_cfgvals[ikey] = iself / bs2
                        for ikey,iself in self.blocked_bootvals.items():
                            result.blocked_bootvals[ikey] = iself / bs2
                except Exception as err:
                    print(type(self.bootvals))
                    print(type(self.cfgvals))
                    print('in bootstrap ',type(bs2))
                    raise EnvironmentError(str(err) +'\nInvalid value to combine with BootStrap class')
            return result

    def __pow__(self,bs2):
        if type(bs2) in BootStrap.parentlist:
            return bs2.__rpow__(self)
        else:
            result = BootStrap(self.nboot,name='^'.join([str(self).replace('BS_',''),str(bs2).replace('BS_','')]))
            if isinstance(bs2,BootStrap):
                result.cfgvals = self.cfgvals ** bs2.cfgvals
                result.bootvals = self.bootvals ** bs2.bootvals
                if hasattr(self,'n_block') and self.n_block > 1:
                    result.blocked_cfgvals = {}
                    result.blocked_bootvals = {}
                    for (ikey,iself),ibs2 in zip(self.blocked_cfgvals.items(),bs2.blocked_cfgvals.values()):
                        result.blocked_cfgvals[ikey] = iself ** ibs2
                    for (ikey,iself),ibs2 in zip(self.blocked_bootvals.items(),bs2.blocked_bootvals.values()):
                        result.blocked_bootvals[ikey] = iself ** ibs2
            else:
                try:
                    result.cfgvals = self.cfgvals ** bs2
                    result.bootvals = self.bootvals ** bs2
                    if hasattr(self,'n_block') and self.n_block > 1:
                        result.blocked_cfgvals = {}
                        result.blocked_bootvals = {}
                        for ikey,iself in self.blocked_cfgvals.items():
                            result.blocked_cfgvals[ikey] = iself ** bs2
                        for ikey,iself in self.blocked_bootvals.items():
                            result.blocked_bootvals[ikey] = iself ** bs2
                except Exception as err:
                    print(type(bs2))
                    raise EnvironmentError(str(err) + '\nInvalid value to combine with BootStrap class')
            return result


    # Right multiplication functions

    def __radd__(self,bs2):
        if isinstance(bs2,BootStrap):
            return bs2.add(self)
        else:
            result = BootStrap(self.nboot,name='+'.join([str(bs2),str(self)]))
            try:
                result.cfgvals = bs2 + self.cfgvals
                result.bootvals = bs2 + self.bootvals
                if hasattr(self,'n_block') and self.n_block > 1:
                    result.blocked_cfgvals = {}
                    result.blocked_bootvals = {}
                    for ikey,iself in self.blocked_cfgvals.items():
                        result.blocked_cfgvals[ikey] = bs2 + iself
                    for ikey,iself in self.blocked_bootvals.items():
                        result.blocked_bootvals[ikey] = bs2 + iself
            except:
                print(type(bs2))
                raise EnvironmentError('Invalid value to combine with BootStrap class')
        return result


    def __rsub__(self,bs2):
        if isinstance(bs2,BootStrap):
            bs2.sub(self)
        else:
            result = BootStrap(self.nboot,name='-'.join([str(bs2),str(self)]))
            try:
                result.cfgvals = bs2 - self.cfgvals
                result.bootvals = bs2 - self.bootvals
                if hasattr(self,'n_block') and self.n_block > 1:
                    result.blocked_cfgvals = {}
                    result.blocked_bootvals = {}
                    for ikey,iself in self.blocked_cfgvals.items():
                        result.blocked_cfgvals[ikey] = bs2 - iself
                    for ikey,iself in self.blocked_bootvals.items():
                        result.blocked_bootvals[ikey] = bs2 - iself
            except:
                print(type(bs2))
                raise EnvironmentError('Invalid value to combine with BootStrap class')
        return result


    def __rmul__(self,bs2):
        if isinstance(bs2,BootStrap):
            bs2.mult(self)
        else:
            result = BootStrap(self.nboot,name='*'.join([str(bs2),str(self)]))
            try:
                result.cfgvals = bs2 * self.cfgvals
                result.bootvals = bs2 * self.bootvals
                if hasattr(self,'n_block') and self.n_block > 1:
                    result.blocked_cfgvals = {}
                    result.blocked_bootvals = {}
                    for ikey,iself in self.blocked_cfgvals.items():
                        result.blocked_cfgvals[ikey] = bs2 * iself
                    for ikey,iself in self.blocked_bootvals.items():
                        result.blocked_bootvals[ikey] = bs2 * iself
            except:
                print(type(bs2))
                raise EnvironmentError('Invalid value to combine with BootStrap class')
        return result


    def __rdiv__(self,bs2):
        if isinstance(bs2,BootStrap):
            bs2.div(self)
        else:
            result = BootStrap(self.nboot,name='/'.join([str(bs2),str(self)]))
            try:
                result.cfgvals = bs2 / self.cfgvals
                result.bootvals = bs2 / self.bootvals
                if hasattr(self,'n_block') and self.n_block > 1:
                    result.blocked_cfgvals = {}
                    result.blocked_bootvals = {}
                    for ikey,iself in self.blocked_cfgvals.items():
                        result.blocked_cfgvals[ikey] = bs2 / iself
                    for ikey,iself in self.blocked_bootvals.items():
                        result.blocked_bootvals[ikey] = bs2 / iself
            except Exception as err:
                print(type(bs2))
                raise EnvironmentError(str(err) +'\nInvalid value to combine with BootStrap class')
        return result

    def __rtruediv__(self,bs2):
        return self.__rdiv__(bs2)

    def __rfloordiv__(self,bs2):
        return self.__rdiv__(bs2)


    def __rpow__(self,bs2):
        if isinstance(bs2,BootStrap):
            bs2.pow(self)
        else:
            result = BootStrap(self.nboot,name='^'.join([str(self).replace('BS_',''),str(bs2).replace('BS_','')]))
            try:
                result.cfgvals = bs2 ** self.cfgvals
                result.bootvals = bs2 ** self.bootvals
                if hasattr(self,'n_block') and self.n_block > 1:
                    result.blocked_cfgvals = {}
                    result.blocked_bootvals = {}
                    for ikey,iself in self.blocked_cfgvals.items():
                        result.blocked_cfgvals[ikey] = bs2 ** iself
                    for ikey,iself in self.blocked_bootvals.items():
                        result.blocked_bootvals[ikey] = bs2 ** iself
            except:
                print(type(bs2))
                raise EnvironmentError('Invalid value to combine with BootStrap class')
        return result



    # unitary arithmetic operations


    def __neg__(self):
        result = BootStrap(self.nboot,name = '-'+self.name,cfgvals = -self.cfgvals,bootvals=-self.bootvals)
        if hasattr(self,'n_block') and self.n_block > 1:
            result.blocked_cfgvals = {}
            result.blocked_bootvals = {}
            for ikey,iself in self.blocked_cfgvals.items():
                result.blocked_cfgvals[ikey] = -iself
            for ikey,iself in self.blocked_bootvals.items():
                result.blocked_bootvals[ikey] = -iself
        return result

    def __pos__(self):
        return self

    def __abs__(self):
        result = BootStrap(self.nboot,name = '|'+self.name+'|',cfgvals = np.abs(self.cfgvals),bootvals=np.abs(self.bootvals))
        if hasattr(self,'n_block') and self.n_block > 1:
            result.blocked_cfgvals = {}
            result.blocked_bootvals = {}
            for ikey,iself in self.blocked_cfgvals.items():
                result.blocked_cfgvals[ikey] = np.abs(iself)
            for ikey,iself in self.blocked_bootvals.items():
                result.blocked_bootvals[ikey] = np.abs(iself)
        return result

    # Casting, just casts the bootstrap values

    def __complex__(self):
        self.cfgvals = self.cfgvals.astype(complex)
        self.bootvals = self.bootvals.astype(complex)
        if hasattr(self,'n_block') and self.n_block > 1:
            for ikey,iself in self.blocked_cfgvals.items():
                self.blocked_cfgvals[ikey] = iself.astype(complex)
            for ikey,iself in self.blocked_bootvals.items():
                self.blocked_bootvals[ikey] = iself.astype(complex)
        self.Stats()
        return self.Avg

    def __int__(self):
        self.cfgvals = self.cfgvals.astype(int)
        self.bootvals = self.bootvals.astype(int)
        if hasattr(self,'n_block') and self.n_block > 1:
            for ikey,iself in self.blocked_cfgvals.items():
                self.blocked_cfgvals[ikey] = iself.astype(int)
            for ikey,iself in self.blocked_bootvals.items():
                self.blocked_bootvals[ikey] = iself.astype(int)
        self.Stats()
        return int(round(self.Avg))

    def __float__(self):
        self.cfgvals = self.cfgvals.astype(np.float64)
        self.bootvals = self.bootvals.astype(np.float64)
        if hasattr(self,'n_block') and self.n_block > 1:
            for ikey,iself in self.blocked_cfgvals.items():
                self.blocked_cfgvals[ikey] = iself.astype(np.float64)
            for ikey,iself in self.blocked_bootvals.items():
                self.blocked_bootvals[ikey] = iself.astype(np.float64)
        self.Stats()
        return self.Avg

    def __round__(self):
        self.cfgvals = self.cfgvals.round
        self.bootvals = self.bootvals.round
        if hasattr(self,'n_block') and self.n_block > 1:
            for ikey,iself in self.blocked_cfgvals.items():
                self.blocked_cfgvals[ikey] = iself.round
            for ikey,iself in self.blocked_bootvals.items():
                self.blocked_bootvals[ikey] = iself.round
        self.Stats()
        return self.Avg

    def __iter__(self):
        return self

    def __next__(self):
        if self.current >= self.nboot:
            self.current = 0
            raise StopIteration
        else:
            self.current += 1
            if self.current >= self.nboot+1:
                return 0.0
            else:
                return self.bootvals[self.current-1]

    def __setitem__(self, key, value ):
        self.bootvals[key] = value

    def __getitem__(self, key):
        return self.bootvals[key]

    def __reversed__(self):
        return self.bootvals[::-1]

    #This would stuff up your statistics

    # def __delitem__(self,key):
    #     self.bootvals = np.delete(self.bootvals,key)
    #     self.nboot -= 1


    def BootHistogram(  self,plot_plane,histtype='stepfilled',label='',color=None,alpha=None,
                        key_val=True,norm=normed_histograms,Median=False):
        self.CheckBoot()
        self.Stats()
        if len(label) == 0:
            label = self.name
        if Median:
            plot_plane.axvline(self.Median,color=color)
            plot_plane.axvline(self.MedStd_Low,color=color,alpha=alpha/2.)
            plot_plane.axvline(self.MedStd_High,color=color,alpha=alpha/2.)
            if key_val and label is not None:
                label = label+' $Res_Median='+MakeValAndErr(self.Avg,self.Std)+'$'
        else:
            plot_plane.axvline(self.Avg,color=color)
            plot_plane.axvline(self.Avg-self.Std,color=color,alpha=alpha/2.)
            plot_plane.axvline(self.Avg+self.Std,color=color,alpha=alpha/2.)
            if key_val and label is not None:
                label = label+' $Res='+MakeValAndErr(self.Median,self.MedStd)+'$'
        plot_plane.hist(self.bootvals.values,histtype=histtype,label=label,color=color,alpha=alpha,density=norm)
        return plot_plane

    def Histogram(self,plot_plane,histtype='stepfilled',label='',color=None,alpha=None):
        self.CheckCfgVals()
        plot_plane.hist(self.cfgvals.values,histtype=histtype,label=self.name+ ' Cfgs '+label,color=color,alpha=alpha)

    def PlotBlockedError(self,plot_plane,label='',color=None,rel_err=True):
        self.CheckBoot_Blocked()
        y_data,x_data = [],[]
        for ikey,iblock in self.blocked_bootvals.items():
            x_data.append(int(ikey.replace('block','')))
            if rel_err:
                y_data.append(iblock.std()/iblock.mean())
            else:
                y_data.append(iblock.std())
        plot_plane.plot(x_data,y_data,color=color,label=self.name+' Berrs '+label)

def TestBoot(n_block=1):
    size = 20000
    this_nboot = 2000
    values = np.random.uniform(size=size)
    testdata = BootStrap(this_nboot,name='test_bootstrap_uniform',cfgvals=values,thisDelVal=False,n_block=n_block)

    values2 = np.arange(size)/size
    testdata2 = BootStrap(this_nboot,name='test_bootstrap_arange',cfgvals=values2,thisDelVal=False,n_block=n_block)

    values3 = np.random.normal(loc=0.5,scale=0.25,size=size)
    testdata3 = BootStrap(this_nboot,name='test_bootstrap_normal',cfgvals=values3,thisDelVal=False,n_block=n_block)
    return testdata,testdata2,testdata3

def TestBlock(n_block=1):
    size = 200
    this_nboot = 200
    values = np.random.uniform(size=size)

    testdata = BootStrap(this_nboot,name='test_bootstrap_uniform',cfgvals=values,thisDelVal=False,n_block=n_block)

    values2 = np.arange(size)/size
    testdata2 = BootStrap(this_nboot,name='test_bootstrap_arange',cfgvals=values2,thisDelVal=False,n_block=n_block)

    values3 = np.random.normal(loc=0.5,scale=0.25,size=size)
    testdata3 = BootStrap(this_nboot,name='test_bootstrap_normal',cfgvals=values3,thisDelVal=False,n_block=n_block)
    return testdata,testdata2,testdata3

if __name__ == '__main__':
    data,data2,data3 = TestBoot()
    print('testing data is in data, data2 and data3')
