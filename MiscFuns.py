#!/usr/bin/env python

from collections import OrderedDict,Callable
import os, errno, sys, re
import glob
import numpy as np
import copy
from types import ModuleType
import pandas as pa
from itertools import product
from copy import deepcopy
from TimeStuff import Timer
import operator as op
import imp
import string
import warnings

def mindex_iloc_md(df,index_list):
    for iindex in index_list:
        df = multindex_iloc(df,iindex)
    return df

def multindex_iloc(df, index):
    label = df.index.levels[0][index]
    return df.iloc[df.index.get_loc(label)]

myeps = np.finfo(0.0).eps

numeric_const_pattern = '[-+]? (?: (?: \d* \. \d+ ) | (?: \d+ \.? ) )(?: [Ee] [+-]? \d+ ) ?'
rx = re.compile(numeric_const_pattern, re.VERBOSE)

def map_str(this_list):
    return [[str(jval) for jval in ival] for ival in this_list]

def DEBUGprint():
    print('DEBUG')

def Series_fix_key(this_series,index,str_rep_from,str_rep_too):
    is_series = isinstance(this_series,pa.Series)
    if is_series:
        this_series = this_series.to_frame('data')
    this_list = this_series.index.names
    this_header = this_list[index]
    this_series.reset_index(inplace=True)
    this_series[this_header] = this_series[this_header].apply(lambda x : x.replace(str_rep_from,str_rep_too)).values
    this_series.set_index(this_list,inplace=True)
    if is_series:
        return this_series['data']
    else:
        return this_series


def StreamToNumb(val):
    val = val.replace('a','1')
    val = val.replace('b','2')
    val = val.replace('c','3')
    val = val.replace('d','4')
    return val.replace('-','')

def logNA_list(value):
    return np.ma.filled(np.log(np.ma.array(value)),-myeps)

def logNA(value):
    if isinstance(value,pa.Series):
        return value.apply(logNA)
    if isinstance(value,pa.DataFrame):
        return value.applymap(logNA)
    if isinstance(value,(list,tuple,np.ndarray)):
        return logNA_list(value)
    else:
        if isinstance(value,complex) :
            if value.imag > 0:
                return np.log(value)
            else:
                return float('NaN')
        elif value < 0:
            return  float('NaN')
        else:
            return np.log(value)
def human_format(num):
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    # add more suffixes if you need them
    return '%.0f%s' % (num, ['', 'K', 'M', 'G', 'T', 'P'][magnitude])

def MultiSplit(this_str,delim_list):
    for idelim in delim_list:
        this_str = this_str.replace(idelim,delim_list[0])
    return this_str.split(delim_list[0])

op_list = [op.add,op.sub,op.mul,op.truediv,op.pow]
op_dict = OrderedDict()
for iop in op_list:
    op_dict[iop.__name__] = iop
op_str_list = list(op_dict.keys())

op_dict_std = OrderedDict()
def add_std(x,y,xbar,ybar):
    return np.sqrt(x**2 + y**2)
def sub_std(x,y,xbar,ybar):
    return np.sqrt(x**2 + y**2)
def mul_std(x,y,xbar,ybar):
    return xbar*y + ybar*x
def div_std(x,y,xbar,ybar):
    return x/ybar - xbar*y/(y**2)
def pow_std(x,y,xbar,ybar):
    return ybar * x**(ybar-1) * x + np.log(xbar)*(xbar**ybar) * y
op_dict_std[op.add.__name__] = add_std
op_dict_std[op.sub.__name__] = sub_std
op_dict_std[op.mul.__name__] = mul_std
op_dict_std[op.truediv.__name__] = div_std
op_dict_std[op.pow.__name__] = pow_std


def GetCfgMissmatch(cfg_list1,cfg_list2,max_len_write = 10):
    # cfg_diff = list(set(cfg_list1).symmetric_difference(set(cfg_list2)))
    left_diff = set(cfg_list1) - set(cfg_list2)
    right_diff = set(cfg_list2) - set(cfg_list1)
    if len(left_diff)+ len(right_diff) < max_len_write:
        left_str,right_str = '',''
        if len(left_diff) > 0:
            left_str = 'cfgs not in cfg file: \n'+'\n'.join(sorted(left_diff))
        if len(right_diff) > 0:
            right_str = 'cfgs not in setup: \n'+'\n'.join(sorted(right_diff))
        return '\n'.join([left_str,right_str])
    else:
        len_diff = len(left_diff)- len(right_diff)
        if len_diff > 0:
            return 'setup has ' + str(len_diff)+ ' more configs'
        elif len_diff < 0:
            return 'cfg file has ' + str(-len_diff)+ ' more configs'


def fmt_file_type(this_type):
    this_type = this_type.replace('feather','fe')
    this_type = this_type.replace('parquet','par')
    this_type = this_type.replace('msgpack','msg')
    this_type = this_type.replace('pickle','p')
    if '.' in this_type:
        return this_type
    else:
        return '.'+this_type


def RemoveAllBoots(boot_data):
    if hasattr(boot_data,'bootvals'):
        boot_data.RemoveBoots()
        boot_data.RemoveVals()
    return boot_data


def convert_to_hex(rgba_color) :
    red = int(rgba_color[0]*255)
    green = int(rgba_color[1]*255)
    blue = int(rgba_color[2]*255)
    alpha = int(rgba_color[3]*255)
    # return '0x{r:02x}{g:02x}{b:02x}'.format(r=red,g=green,b=blue)
    return (red,green,blue,alpha)

def fmt_Qt5_col(this_class):
    if isinstance(this_class,(str,list,tuple,np.ndarray)):
        return this_class
    else:
        return [x/255. for x in this_class.getRgb()]


def get_val_float_0(val):
    this_val = get_val(val,this_type=float,suppress_warn=True)
    if isinstance(this_val,float):
        return this_val
    elif isinstance(this_val[0],float):
        return this_val[0]
    else:
        warnings.warn('float not found in get_val_float_0 call, using first character ordering')
        return string.ascii_lowercase.index(this_val[0])

def get_val_float(val):
    return get_val(val,this_type=float)

def get_val(val,Phys=True,this_type=int,suppress_warn=False):
    if not isinstance(val,str):
        return this_type(val)
    out_str = rx.findall(val)
    if len(out_str) == 0:
        if not suppress_warn: warnings.warn(this_type+' not found in "get_val" call')
        return val
    if len(out_str) == 1:
        return this_type(out_str[0])
    else:
        return list(map(this_type,out_str))
    # if not isinstance(val,str):
    #     return val
    # new_str = ''
    # for ichar in val:
    #     if ichar.isdigit() or ichar == '.':
    #         new_str += ichar
    # if len(new_str) > 0:
    #     try:
    #         return int(new_str)
    #     except:
    #         return float(new_str)
    # else:
    #     return 0

def ReductList(this_list,max_len=5):
    if isinstance(this_list,(list,tuple,np.ndarray)):
        if len(this_list) > max_len:
            return list(this_list[:max_len/2+1]) + ['...'] + list(this_list[-max_len/2:])
    return this_list

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    return False

def Remove_Empty_Str(thislist):
    out_list = []
    for ival in thislist:
        if len(ival) != 0:
            out_list.append(ival)
    return out_list

def MIndex_to_numpy(this_df):
    return this_df.reindex(list(product(*list(map(tuple, this_df.index.levels))))).values.reshape(list(map(len, this_df.index.levels)))


##
def flip_fun_2arg(this_fun):
    def result_fun(var1,var2):
        return this_fun(var2,var1)
    return result_fun

# def Intersect_DataFrame_2col(this_df):
#     out_df = pa.Series('False',index=this_df.index)
#     for istream,stream_df in this_df.groupby(level='stream'):
#         stream_series = stream_df.iloc[:,1][stream_df.isin(stream_df.iloc[:,0].values).iloc[:,1].values]
#         for ikey,ival in stream_series.iteritems():
#             out_df[ikey] = ival
#     return out_df[out_df != 'False']

## if only one dataframe has a xsrc list, it MUST be df_one
def Intersect_DataFrame_2col_xsrc(df_one,df_two):
    cfglist, xsrc_list = [],[[]]
    ilist = []
    no_two_xsrc = not ('xsrc_list' in df_two.columns)
    no_xsrc = not ('xsrc_list' in df_one.columns)
    for (istream,iccfg),idf in df_one.iterrows():
        if istream in df_two['configs'].index.levels[0]:
            icfg = idf['configs']
            this_where = np.where(df_two['configs'][istream].values == icfg)[0]
            if len(this_where) == 0: continue
            second_key = this_where[0]
            if 'xsrc_list' in idf:
                for ixsrc in idf['xsrc_list']:
                    if no_two_xsrc or (ixsrc in df_two['xsrc_list'][istream].values[second_key]):
                        xsrc_list[-1].append(ixsrc)
            if no_xsrc or len(xsrc_list[-1]) > 0:
                ilist.append((istream,iccfg))
                cfglist.append(icfg)
                xsrc_list.append([])
    if len(ilist) == 0:
        print('DataFrame 1:')
        print(df_one)
        print()
        print('DataFrame 2:')
        print(df_two)
        print()
        raise EnvironmentError('No common configs')
    indicies = pa.MultiIndex.from_tuples(ilist,names=['stream','config_number'])
    cfg_df = pa.DataFrame(cfglist,columns=['configs'],index=indicies)
    if not no_xsrc: cfg_df.loc[:,'xsrc_list'] = pa.Series(xsrc_list[:-1],index=indicies)
    return cfg_df

def Series_TO_ODict(this_series):
    if isinstance(this_series.index,pa.MultiIndex):
        outdict = ODNested()
        keysize = len(list(this_series.keys())[0])
        for ikey,ival in this_series.items():
            ikey = list(map(str,ikey))
            if keysize == 1:
                outdict[ikey[0]] = ival
            elif keysize == 2:
                outdict[ikey[0]][ikey[1]] = ival
            elif keysize == 3:
                outdict[ikey[0]][ikey[1]][ikey[2]] = ival
            elif keysize == 4:
                outdict[ikey[0]][ikey[1]][ikey[2]][ikey[3]] = ival
        return outdict
    else:
        return OrderedDict(list(zip(*[list(map(str,list(this_series.keys()))),this_series.values])))

def rreload(module, paths=[''], mdict={}):
    """Recursively reload modules."""
    if module not in mdict:
        # modules reloaded from this module
        mdict[module] = []
    imp.reload(module)
    for attribute_name in dir(module):
        attribute = getattr(module, attribute_name)
        if type(attribute) is ModuleType:
            if attribute not in mdict[module]:
                if attribute.__name__ not in sys.builtin_module_names:
                    if os.path.dirname(attribute.__file__) in paths:
                        mdict[module].append(attribute)
                        rreload(attribute, paths, mdict)
    imp.reload(module)

def intersect(*d):
    if len(d) == 1: return d[0]
    sets = iter(map(set, d))
    result = next(sets)
    for s in sets:
        result = result.intersection(s)
    return list(result)

class DOO(OrderedDict):
    # Source: http://stackoverflow.com/a/6190500/562769
    def __init__(self, default_factory=None, *a, **kw):
        if (default_factory is not None and
            not isinstance(default_factory, Callable)):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
            return type(self), args, None, None, iter(list(self.items()))

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        return type(self)(self.default_factory,
        copy.deepcopy(list(self.items())))

    def __repr__(self):
        return 'DOO(%s, %s)' % (self.default_factory,
                                OrderedDict.__repr__(self))



# class OrderedDefaultDict(OrderedDict):
#     def __init__(self, default_factory=None, *args, **kwargs):
#         OrderedDict.__init__(self,*args,**kwargs)
#         self.default_factory = default_factory
#     def __missing__(self, key):
#         if self.default_factory is None:
#             raise KeyError(key)
#         val = self[key] = self.default_factory()
#         return val



def ODNested():
    # return OrderedDefaultDict(ODNested)
    return DOO(ODNested)
    # return defaultdict(ODNested)


def FixDictArray(thislist,ikey):
    outDict = OrderedDict()
    for ic,ilist in enumerate(thislist):
        outDict[ikey+str(ic+1)] = ilist
    return outDict


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

def NumbToFileName(val):
    if isinstance(val,int):
        return str(val)
    else:
        if val > 0.01:
            return '{:.2f}'.format(val)
        else:
            return '0'


def GetFileNames(thisdir,strcomp):
    return glob.glob(thisdir+"/"+strcomp[0]+'*'+strcomp[1])

#Pulls out a class "atribute" out of a numpy tensor "data" of classes ##
#e.g. return data[:,:,...,:,:].atribute
def Pullflag(data,atribute):
    npdata = np.array(data)
    dataout = np.array([])
    if len(npdata.flatten()) == 0: return data
    for idata in npdata.flatten():
        flagdim = np.array(getattr(idata,atribute)).shape
        dataout = np.append(dataout,getattr(idata,atribute))
    return np.reshape(dataout,npdata.shape + flagdim)

def VecDelta(tvec,val):
    return list(map(int,val==np.array(tvec)))


## assuming class2 is what is being read in
## and class1 is what the parameters are
def CheckClass(class1,class2,checklist):
    sub_vector_list = ['tflowlist','checkmomlist','boot_tsum_list','tflow_list']
    fun_list = ['Fun','FunDer']
    for icheck in checklist:
        if icheck in fun_list:
            icheck_name = icheck+'_name'
            if icheck in list(class1.keys()):
                c1check = class1[icheck].__name__
            elif icheck_name in list(class1.keys()):
                c1check = class1[icheck_name]
            else:
                print(icheck , 'not in file ')
                return False
            if icheck in list(class2.keys()):
                c2check = class2[icheck].__name__
            elif icheck_name in list(class2.keys()):
                c2check = class2[icheck_name]
            else:
                print(icheck , 'not in setup class ')
                return False

            if  c1check != c2check:
                print('Fit Functions differe in classes:')
                print('class1:',c1check)
                print('class2:',c2check)
                return False
        elif (icheck not in list(class2.keys())):
            print(icheck , 'not in setup class ')
            return False
        elif (icheck not in list(class1.keys())):
            print(icheck , 'not in file ')
            return False
        elif icheck in sub_vector_list:
            for ival in class1[icheck]:
                if ival not in class2[icheck]:
                    if 'HumanFile' in list(class2.keys()):
                        print('value ' , ival , ' not found in ', class2['HumanFile'])
                    else:
                        print('value ' , ival , ' not found in class2')
                    return False
        elif (('_cfgs' in icheck or 'cfglist' in icheck) and 'fo_for_cfgs' not in icheck) and not isinstance(class1[icheck],bool):
            if 'xsrc_list' not in class1[icheck] or 'xsrc_list' not in class2[icheck]:
                for (istream,iccfg),icfg in class1[icheck]['configs'].items():
                    if istream not in class2[icheck].index.levels[0] or icfg not in class2[icheck]['configs'][istream].values:
                        if 'HumanFile' in list(class2.keys()):
                            print('config number ' , istream ,icfg , ' not found in ', class2['HumanFile'])
                        else:
                            print('config number ' , istream ,icfg , ' not found in class2')
                        return False
                for (istream,iccfg),icfg in class2[icheck]['configs'].items():
                    if istream not in class1[icheck].index.levels[0] or icfg not in class1[icheck]['configs'][istream].values:
                        if 'HumanFile' in list(class1.keys()):
                            print('config number ' , istream ,icfg , ' not found in ', class1['HumanFile'])
                        else:
                            print('config number ' , istream ,icfg , ' not found in class1')
                        return False
            else:
                for (istream,iccfg),icfg in class1[icheck]['configs'].items():
                    xsrc_test,cfg_test = False,False
                    stream_test = istream in class2[icheck].index.levels[0]
                    if stream_test:
                        cfg_test = icfg in class2[icheck]['configs'][istream].values
                        if cfg_test:
                            thisindex = np.where(class2[icheck]['configs'][istream].values == icfg)[0][0]
                            xsrc_set1 = set(class1[icheck]['xsrc_list'][istream,iccfg])
                            xsrc_set2 = set(class2[icheck]['xsrc_list'][istream].iloc[thisindex])
                            xsrc_test = xsrc_set1 == xsrc_set2
                    if not cfg_test:
                        if 'HumanFile' in list(class2.keys()):
                            print('config number missmatch for' , istream ,icfg , ' not the same as ', class2['HumanFile'])
                        else:
                            print('config number missmatch for' , istream ,icfg , ' not the same as class2')
                        return False
                    elif not xsrc_test:
                        print('x_source missmatch')
                        for ixsrc in sorted(xsrc_set1-xsrc_set2):
                            if 'HumanFile' in list(class2.keys()):
                                print(class2['HumanFile'] , 'is missing ',istream ,icfg ,ixsrc)
                            else:
                                print('class2 has is missing ',istream ,icfg ,ixsrc)
                        for ixsrc in sorted(xsrc_set2-xsrc_set1):
                            if 'HumanFile' in list(class2.keys()):
                                print(class1['HumanFile'] , 'is missing ',istream ,icfg ,ixsrc)
                            else:
                                print('class1 has is missing ',istream ,icfg ,ixsrc)
                        return False
                for (istream,iccfg),icfg in class2[icheck]['configs'].items():
                    cfg_test = False
                    stream_test = istream in class1[icheck].index.levels[0]
                    if stream_test:
                        if icfg not in class1[icheck]['configs'][istream].values:
                            if 'HumanFile' in list(class1.keys()):
                                print('config number missmatch for' , istream ,icfg , ' not the same as ', class1['HumanFile'])
                            else:
                                print('config number missmatch for' , istream ,icfg , ' not the same as class1')
                            return False

        elif isinstance(class1[icheck],list) or isinstance(class1[icheck],tuple) or isinstance(class1[icheck],np.ndarray):
            if len(class1[icheck]) != len(class2[icheck]):
                checkbool = False
            else:
                checkbool = class1[icheck] == class2[icheck]
            if isinstance(checkbool,np.ndarray): checkbool = all(checkbool)
            if not checkbool:
                print('missmatch' , icheck, class1[icheck],class2[icheck])
                return False
        else:
            # print 'DEBUG'
            # print icheck
            # print class1[icheck]
            # print class2[icheck]
            lcheck,rcheck = class1[icheck],class2[icheck]
            if isinstance(lcheck,str):
                lcheck = lcheck.replace('.bak','')
                lcheck = lcheck.replace(r'//',r'/')
            if isinstance(rcheck,str):
                rcheck = rcheck.replace('.bak','')
                rcheck = rcheck.replace(r'//',r'/')
            checkbool = rcheck == lcheck
            if isinstance(checkbool,np.ndarray): checkbool = all(checkbool)
            if not checkbool:
                print('missmatch' , icheck, rcheck,lcheck)
                return False
    return True

def CombineListCfgs(list_cfgs,this_attr):
    if len(list_cfgs) == 0:
        return None
    elif len(list_cfgs) == 1:
        return getattr(list_cfgs[0],this_attr)
    else:
        parselist = iter(list_cfgs)
        out_cfgs = CombineCfgs(getattr(next(parselist),this_attr),getattr(next(parselist),this_attr))
        for icfg in parselist:
            out_cfgs = CombineCfgs(out_cfgs,getattr(icfg,this_attr))
        return out_cfgs


def MakeStreamCfgs(cfg_df):
    cfg_df.loc[:,'stream_cfgs'] = pa.Series(['00'.join((istream,icfg)) for (istream,iccfg),icfg in cfg_df['configs'].items()],index=cfg_df.index)
    return cfg_df

def CombineCfgs(cfga,cfgb):
    if 'xsrc_list' in cfga:
        output = Intersect_DataFrame_2col_xsrc(cfga,cfgb)
    else:
        output = Intersect_DataFrame_2col_xsrc(cfgb,cfga)
    # print 'DEBUG'
    # print 'cfga'
    # print cfga
    # print
    # print 'cfgb'
    # print cfga
    # print
    # print 'comb_cfg'
    # print MakeStreamCfgs(output)
    # print
    return MakeStreamCfgs(output)

def ReduceCfgs(twopt,threept):
    output = deepcopy(threept)
    if 'stream_cfgs' not in threept.columns:
        threept.loc[:,'stream_cfgs'] = pa.Series(['00'.join((istream,icfg)) for (istream,iccfg),icfg in threept['configs'].items()],index=threept.index)
    if 'stream_cfgs' not in twopt.columns:
        twopt.loc[:,'stream_cfgs'] = pa.Series(['00'.join((istream,icfg)) for (istream,iccfg),icfg in twopt['configs'].items()],index=twopt.index)
    for ikey,stream_cfg in threept['stream_cfgs'].items():
        if stream_cfg not in twopt['stream_cfgs'].values:
            output.drop(ikey,inplace=True)
    return output
    #
    # twoptout = OrderedDict()
    # for icfg,isrclist in twopt.iteritems():
    #     if icfg in threept.keys():
    #         twoptout[icfg] = isrclist
    # for icfg in threept.keys():
    #     if icfg not in twopt.keys():
    #         print icfg
    #         raise IOError('cfg in 3pt not in 2pt, please delete three point correlator configuration')
    # return twoptout


def GetInfoList(thisInfo):
    if 'sm2ptRat' in list(thisInfo.keys()):
        if 'jsmCoeffs' in list(thisInfo.keys()):
            jsmcoeffs = thisInfo['jsmCoeffs']
        else:
            raise IOError('jsm coefficients are needed when passing in jsm2ptRat')
        if 'sinklab' in list(thisInfo.keys()):
            sinklab = thisInfo['sinklab']
        else:
            sinklab = 'Comb'

        sm2ptList = thisInfo['sm2ptRat']
        thisInfoList = []

        for jsm in thisInfo['sm2ptRat']:
            thisInfoPass = copy.deepcopy(thisInfo)
            # thisInfoPass['Interp'] = 'CPEven'
            thisInfoPass['jsm'] = jsm
            thisInfoList.append(thisInfoPass)
    else:
        jsmcoeffs = [1.0]
        if 'jsm' in list(thisInfo.keys()):
            sinklab = 'jsm'+str(thisInfo['jsm'])
            sm2ptList = thisInfo['jsm']
        else:
            sm2ptList = 0
            sinklab = ''
        thisInfoList = [thisInfo]
    return sinklab,jsmcoeffs,sm2ptList,thisInfoList


def GenSinkFun(jsmcoeffs,show_timer=True):
    def SinkCombFun(*c2list):
        if len(c2list) != len(jsmcoeffs) :
            print(jsmcoeffs)
            print(c2list)
            raise EnvironmentError('SinkCombFun expects '+str(len(jsmcoeffs))+' c2funs, got '+str(len(jsmcoeffs)))
        if show_timer: thistimer = Timer(linklist=c2list,name='Sink Comb: ')
        c2list = iter(c2list)
        this_jsmc = iter(jsmcoeffs)
        # if isinstance(coeff,basestring):
        #     print 'Warning, combining correlators had incorrect coefficient', coeff
        #     continue
        icoeff = next(this_jsmc)
        ic2 = next(c2list)
        if icoeff == 1.0:
            output = ic2
        elif icoeff != 0.0:
            output = ic2*icoeff
        else:
            output = 0.0
        if show_timer: thistimer.Lap()
        for ic2,icoeff in zip(c2list,this_jsmc):
            if isinstance(icoeff,str):
                print('Warning, combining correlators had incorrect coefficient', icoeff)
                continue
            # print
            # print 'test'
            # print 'coeff',coeff
            # if hasattr(ic2,'NNQ_Stats'):
            #     print 'NNQ'
            #     print [ival.bootvals[0] for ival in ic2.NNQ_Stats['boot'].values[:5]]
            # else:
            #     print 'C2'
            #     print [ival.bootvals[0] for ival in ic2.C2_Stats['boot'].values[:5]]

            # print
            # print 'test2'
            # print [ival.bootvals[0] for ival in output.values()]
            if icoeff == 1.0:
                output = output + ic2
            elif icoeff != 0.0:
                output = output + ic2*icoeff
            if show_timer: thistimer.Lap()
            # print 'testComb'
            # if hasattr(ic2,'NNQ_Stats'):
            #     print [ival.bootvals[0] for ival in output.NNQ_Stats['boot'].values[:5]]
            # else:
            #     print [ival.bootvals[0] for ival in output.C2_Stats['boot'].values[:5]]
            # print
        return output
    return SinkCombFun

def CreateNameAndLL(thisset,thissink,from_index=0):
    firstkey = list(thisset.keys())[from_index]
    thisname = thisset[firstkey].name
    thisleglab = thisset[firstkey].LegLab
    firstjsm = thisset[firstkey].jsm
    thisname = thisname.replace(firstjsm,thissink)
    thisleglab = thisleglab.replace(firstjsm,thissink)
    jsmlab = '_'.join([thissm.jsm for thissm in list(thisset.values())])
    return thisname,thisleglab,jsmlab
