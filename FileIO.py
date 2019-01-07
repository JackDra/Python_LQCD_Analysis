#!/usr/bin/env python
## THIS FILE IS TEMPORARY, FOR DEBUGGING AND SIMPLE DEFAULT PARAMETERS
## TODO: either default parameters in file, or even GUI?...

from collections import OrderedDict
from MiscFuns import is_number,fmt_file_type
from xmltodict import unparse,parse
from XmlFormatting import FormatToDictAvgStdChi,FormatToDictAvgStd
import pickle as pickle
# import dill as pickle
# import pickle as dill
import dill
import os,shutil
import pandas as pa
from Params import write_excel,ShowXmlWrite,functiondir
from Params import All_Classes_Names,ShowRead,cfg_file_type,Debug_Mode,test_file
import warnings
import fcntl

py2_encoding = 'latin1'

def WriteWithMeta(df_data,meta_data,this_file,file_type=cfg_file_type):
    if 'None' in cfg_file_type:
        # print
        # print 'forced code to not write cfgs'
        return
    if 'to_' not in file_type:
        file_type = 'to_'+file_type
    file_ext = fmt_file_type(file_type.replace('to_',''))
    if '.' not in file_ext:
        file_ext = '.'+file_ext
    if file_ext not in this_file:
        this_file += file_ext
    if file_ext in ['.csv','.json','.html']:
        getattr(df_data,file_type)(this_file)
        getattr(pa.DataFrame(meta_data),file_type)(this_file+'.meta')
    else:
        try:
            with open(this_file,'wb') as f:
                fcntl.flock(f, fcntl.LOCK_EX)
                pickle.dump(meta_data,f)
                getattr(df_data,file_type)(f)
                fcntl.flock(f, fcntl.LOCK_UN)
        except Exception as err:
            try:
                with open(this_file,'wb') as f:
                    pickle.dump(meta_data,f)
                    getattr(df_data,file_type)(f)
            except Exception as err2:
                raise IOError('\n' + '\n'.join([str(err),str(err2),'Error writing file',this_file]))

def ReadWithMetaPy2(this_file,file_type=cfg_file_type):
    if 'None' in cfg_file_type:
        # print
        # print 'forced code to not read cfgs'
        return [],pa.DataFrame()
    if 'read_' not in file_type:
        file_type = 'read_'+file_type
    file_ext = fmt_file_type(file_type.replace('read_',''))
    if '.' not in file_ext:
        file_ext = '.'+file_ext
    if file_ext not in this_file:
        this_file += file_ext
    if file_ext in ['.csv','.json','.html']:
        raise EnvironmentError('reading human readable types: csv, json, html not supported')
        # meta_data = getattr(pa,file_type)(this_file+'.meta')
        # df_data = getattr(pa,file_type)(this_file)
    else:
        try:
            with open(this_file,'rb') as f:
                fcntl.flock(f, fcntl.LOCK_EX)
                meta_data = pickle.load(f,encoding=py2_encoding)
                df_data = getattr(pa,file_type)(f)
                fcntl.flock(f, fcntl.LOCK_UN)
        except Exception as err:
            try:
                with open(this_file,'rb') as f:
                    meta_data = pickle.load(f,encoding=py2_encoding)
                    df_data = getattr(pa,file_type)(f)
            except Exception as err2:
                raise IOError('\n' + '\n'.join([str(err),str(err2),'Error loading file',this_file]))
    return meta_data,df_data


def ReadWithMeta(this_file,file_type=cfg_file_type):
    if 'None' in cfg_file_type:
        # print
        # print 'forced code to not read cfgs'
        return [],pa.DataFrame()
    if 'read_' not in file_type:
        file_type = 'read_'+file_type
    file_ext = fmt_file_type(file_type.replace('read_',''))
    if '.' not in file_ext:
        file_ext = '.'+file_ext
    if file_ext not in this_file:
        this_file += file_ext
    if file_ext in ['.csv','.json','.html']:
        raise EnvironmentError('reading human readable types: csv, json, html not supported')
        # meta_data = getattr(pa,file_type)(this_file+'.meta')
        # df_data = getattr(pa,file_type)(this_file)
    else:
        try:
            with open(this_file,'rb') as f:
                fcntl.flock(f, fcntl.LOCK_EX)
                meta_data = pickle.load(f)
                df_data = getattr(pa,file_type)(f)
                fcntl.flock(f, fcntl.LOCK_UN)
        except Exception as err:
            try:
                with open(this_file,'rb') as f:
                    meta_data = pickle.load(f)
                    df_data = getattr(pa,file_type)(f)
            except Exception as err2:
                warnings.warn('Error attempting to read pickle file, attempting to use encoding'+py2_encoding)
                ReadWithMetaPy2(this_file,file_type=file_type)
    return meta_data,df_data


## CURRENLTY LOCKING DOES NOT WORK WITH MULTICORE....

def WriteExcel(thisfile,df_list,flowcfglist=None,cfglist=None,params=None,Bak=True):
    if Bak: BackupFile(thisfile)
    if write_excel:
        with pa.ExcelWriter(thisfile) as writer:
            for ikey,idf in df_list.items():
                if 'boot' in idf.columns: del idf['boot']
                if 'Auto' in idf.columns: del idf['Auto']
                if 'Alphaboot' in idf.columns: del idf['Alphaboot']
                if 'AlphaAuto' in idf.columns: del idf['AlphaAuto']
                if 'EffM' in idf.columns: del idf['EffM']
                if 'FF' in idf.columns: del idf['FF']
                idf.to_excel(writer, sheet_name=ikey)

            if flowcfglist is not None:
                if not isinstance(flowcfglist,pa.Series):
                    flowcfglist = pa.Series(flowcfglist)
                flowcfglist.to_frame(name='flow_configurations').to_excel(writer,sheet_name='flow_cfg_list')
            if cfglist is not None:
                if not isinstance(cfglist,pa.DataFrame):
                    cfglist = pa.DataFrame(cfglist)
                outcfg_list = cfglist['configs'].to_frame(name='configurations')
                if 'xsrc_list' in list(cfglist.keys()): outcfg_list['n_xsrc'] = cfglist['xsrc_list'].apply(lambda x : len(x))
                outcfg_list.to_excel(writer,sheet_name='cfg_list')
            if params is not None:
                if not isinstance(params,pa.Series):
                    params = pa.Series(params,name='Parameters')
                params.to_frame(name='Parameters').to_excel(writer,sheet_name='Parameters')



def WriteXml(thisfile,outputdict,Bak=True):
    if Bak: BackupFile(thisfile)
    if ShowXmlWrite and 'Parameters' not in thisfile: print('Writing to ', thisfile)
    try:
        f = open(thisfile,'w')
        fcntl.flock(f, fcntl.LOCK_EX)
        f.write( unparse(outputdict,pretty=True).replace('\t','    '))
        fcntl.flock(f, fcntl.LOCK_UN)
        f.close()
    except Exception as err:
        # for ikey,iout in outputdict['Results'].iteritems():
        #     print ikey
        #     print type(iout)
        #     print iout
        #     print
        try:
            f.close()
        except:
            pass
        try:
            f = open(thisfile,'w')
            fcntl.flock(f, fcntl.LOCK_EX)
            f.write( unparse(outputdict,pretty=True).replace('\t','    '))
            fcntl.flock(f, fcntl.LOCK_UN)
            f.close()
        except Exception as err2:
            CheckFTDAS(outputdict)
            with open(thisfile,'w') as f:
                f.write(outputdict.__repr__())
            raise IOError(err + '\n'+err2 + '\n Error with outputting xml file \n' + thisfile)

def WriteFun_NoLock(this_fun):
    if isinstance(this_fun,str):
        this_file = functiondir+str(this_fun)+'.fun'
    elif isinstance(this_fun,int):
        this_file = functiondir+'par'+str(this_fun)+'.fun'
    elif hasattr(this_fun,'file_name'):
        this_file = functiondir+str(this_fun.file_name)+'.fun'
    else:
        this_file = functiondir+str(this_fun.__name__)+'.fun'
    if os.path.isfile(this_file):
        return
    try:
        with open(this_file,'wb') as f:
            dill.dump(this_fun,f)
    except Exception as err:
        raise IOError(str(err)+'\nFailed writing pickled file: \n'+this_file)


def WriteFuns(*these_funs):
    for ivar in these_funs:
        if isinstance(ivar,str):
            this_file = functiondir+str(ivar)+'.fun'
        elif isinstance(ivar,int):
            this_file = functiondir+'par'+str(ivar)+'.fun'
        elif hasattr(ivar,'file_name'):
            this_file = functiondir+str(ivar.file_name)+'.fun'
        else:
            this_file = functiondir+str(ivar.__name__)+'.fun'
        if os.path.isfile(this_file):
            continue
        try:
            with open(this_file,'wb') as f:
                fcntl.flock(f, fcntl.LOCK_EX)
                dill.dump(ivar,f)
                fcntl.flock(f, fcntl.LOCK_UN)
        except Exception:
            WriteFun_NoLock(ivar)
            # raise IOError(str(err)+'\nFailed writing pickled file: \n'+this_file)



# def ReadFuns_Py2(*fun_names):
#     fun_list = []
#     for ifun in fun_names:
#         fun_file = functiondir+ifun+'.fun'
#         if os.path.isfile(fun_file):
#             # try:
#             with open(fun_file,'rb') as f:
#                 fun_list.append(dill.load(f))
#             # except Exception as err:
#             #     raise IOError(str(err)+'\nFailed loading pickled file: \n'+fun_file)
#         else:
#             print('FNF: '+fun_file)
#             fun_list.append('FNF: '+fun_file)
#     return fun_list

def ReadFuns_NoLock(*fun_names):
    fun_list = []
    for ifun in fun_names:
        fun_file = functiondir+ifun+'.fun'
        if os.path.isfile(fun_file):
            try:
                with open(fun_file,'rb') as f:
                    fun_list.append(dill.load(f))
            except Exception as err:
                raise IOError(str(err)+'\nFailed loading pickled file: \n'+fun_file)
            # except Exception:
            #     fun_list.append(ReadFuns_Py2(ifun)[0])
        else:
            print('FNF: '+fun_file)
            fun_list.append('FNF: '+fun_file)
    return fun_list


def ReadFuns(*fun_names):
    fun_list = []
    for ifun in fun_names:
        fun_file = functiondir+ifun+'.fun'
        if os.path.isfile(fun_file):
            try:
                with open(fun_file,'rb') as f:
                    fcntl.flock(f, fcntl.LOCK_EX)
                    fun_list.append(dill.load(f))
                    fcntl.flock(f, fcntl.LOCK_UN)
            except Exception:
                fun_list.append(ReadFuns_NoLock(ifun)[0])
                # raise IOError(str(err)+'\nFailed loading pickled file: \n'+fun_file)
        else:
            print('FNF: '+fun_file)
            fun_list.append('FNF: '+fun_file)
    return fun_list


def WriteDill(thisfile,outputdict,this_err=''):
    try:
        with open(thisfile,'wb') as f:
            dill.dump(outputdict,f)
    except Exception as err2:
        CheckFTDAS(outputdict,this_err)
        raise IOError(str(this_err) + '\n'+str(err2)+ '\n Error with outputting pickled file \n' + str(thisfile))


def TestForFunction(dictin,this_string=[],temp_file=test_file):
    if isinstance(dictin,dict) or isinstance(dictin,OrderedDict) :
        for ikey,idata in dictin.items():
            pass_string = this_string + [ikey]
            TestForFunction(idata,this_string=pass_string,temp_file=temp_file)
    elif isinstance(dictin,pa.Series):
        for ikey,idata in dictin.items():
            pass_string = this_string + ['Series:'+str(ikey)]
            TestForFunction(idata,this_string=pass_string,temp_file=temp_file)
    elif isinstance(dictin,pa.DataFrame) :
        for icol,idata in dictin.items():
            for ikey,jdata in idata.items():
                pass_string = this_string + ['DF:'+str(icol)+';'+str(ikey)]
                TestForFunction(jdata,this_string=pass_string,temp_file=temp_file)
    elif any([iclass in str(type(dictin)) for iclass in All_Classes_Names]):
        for ikey,idata in dictin.__dict__.items():
            if '__' not in ikey:
                pass_string = this_string + ['Class:'+ikey]
                TestForFunction(idata,this_string=pass_string,temp_file=temp_file)
    else:
        try:
            with open(temp_file,'wb') as f:
                pickle.dump(dictin,f)
        except Exception as err:
            print()
            print('-------------------------------------------')
            print(str(err))
            print('An error occured at:')
            for itab,istr in enumerate(this_string):
                print(' '*itab,istr)
            print()
            print('value to print')
            print(dictin)
            print()
            print('has type')
            print(type(dictin))
            print('-------------------------------------------')
            print()
        # for itab,istr in enumerate(this_string):
        #     print '   '*itab , str(istr)

def WritePickle_NoLock(thisfile,outputdict,this_err=''):
    try:
        with open(thisfile,'wb') as f:
            pickle.dump(outputdict,f)
    except Exception as err:
        this_err += '\n'+str(err)
        print('WARNING: cannot write with cPickle, using dill for file:')
        print(thisfile)
        if Debug_Mode:
            print()
            print('DEBUGGING the cPickle problem...')
            TestForFunction(outputdict)
        WriteDill(thisfile,outputdict,this_err=this_err)


def WritePickle(thisfile,outputdict,Bak=True):
    if Bak: BackupFile(thisfile)
    try:
        with open(thisfile,'wb') as f:
            # fcntl.flock(f, fcntl.LOCK_EX)
            pickle.dump( outputdict , f )
            # fcntl.flock(f, fcntl.LOCK_UN)
    except Exception as err:
        print(err)
        WritePickle_NoLock(thisfile,outputdict,this_err=str(err))



def BackupFile(*thisfilelist):
    for thisfile in thisfilelist:
        if 'Parameters' in thisfile: continue
        try:
            if os.path.isfile(thisfile):
                shutil.move(thisfile, thisfile+'.bak')
        except:
            raise IOError('ERROR with backing up file '+thisfile)


def ReadDillWrap(thisfile,this_err = ''):
    try:
        with open( thisfile, "rb" ) as pfile:
            data = dill.load(pfile)
        # print 'Warning, dill used for reading file: \n',pfile
        return data
    except Exception as err2:
        this_err += '\n'+str(err2)
        print('WARNING: cannot read with cPickle, using dill for file:')
        print(thisfile)
        return ReadPicklePy2(thisfile,this_err)

def ReadPicklePy2(thisfile,this_err = ''):
    try:
        with open( thisfile, "rb" ) as pfile:
            data = pickle.load(pfile,encoding=py2_encoding)
        # print 'Warning, dill used for reading file: \n',pfile
        return data
    except Exception as err2:
        raise IOError(str(this_err)+'\n'+str(err2)+'\nFailed loading pickled file with encoding '+py2_encoding+': \n'+thisfile)

def ReadPickleWrap_NoLock(thisfile,this_err=''):
    try:
        with open( thisfile, "rb" ) as pfile:
            data = pickle.load(pfile)
        return data
    except Exception as err:
        this_err += '\n'+str(err)
        print('WARNING: cannot read with cPickle, using dill for file:')
        print(thisfile)
        return ReadDillWrap(thisfile,this_err)


def ReadPickleWrap(thisfile,thisShowRead=ShowRead):
    if thisShowRead and 'Parameters' not in thisfile: print('Reading Pickled file ', thisfile)
    try:
        with open( thisfile, "rb" ) as pfile:
            fcntl.flock(pfile, fcntl.LOCK_EX)
            data = pickle.load(pfile)
            fcntl.flock(pfile, fcntl.LOCK_UN)
        return data
    except Exception as err:
        return ReadPickleWrap_NoLock(thisfile,this_err = str(err))


def CheckFTDAS(dictin,prev_key='first',this_err=''):
    if isinstance(dictin,dict) or isinstance(dictin,OrderedDict):
        for ikey,idata in dictin.items():
            CheckFTDAS(idata,prev_key=ikey,this_err=this_err)
    else:
        try:
            str(dictin)
        except Exception as err:
            out_msg = 'final value in dictionary is causing xml write error'
            thisout =  '\n'+'\n'.join([str(this_err),str(err),prev_key,str(type(dictin)),out_msg])
            raise TypeError(thisout)

def RecFTDAS(dictin):
    if isinstance(dictin,dict) or isinstance(dictin,OrderedDict):
        for ikey,idata in dictin.items():
            dictin[ikey] = RecFTDAS(idata)
        return dictin
    else:
        try:
            dictout = float(dictin)
            if int(dictout) == dictout:
                dictout = int(dictout)
        except:
            try:
                if len(str(dictin).split()) == 3 and is_number(str(dictin).split()[0]):
                    dictout = FormatToDictAvgStdChi(str(dictin))
                elif len(str(dictin).split()) == 2 and is_number(str(dictin).split()[0]):
                    dictout = FormatToDictAvgStd(str(dictin))
                else:
                    dictout = str(dictin)
            except:
                print(dictin)
                print(type(dictin))
                raise TypeError('final value in dictionary is not string')
        return dictout

def FormatXml(data):
    data = data.replace(r'$32^{3}\times_64\_\quad_m_{\pi}=0.410_GeV_\_-def-\_<Q>$_','mpi410')
    data = data.replace(r'$32^{3}\times_64\_\quad_m_{\pi}=0.568_GeV_\_-def-\_<Q>$_','mpi570')
    data = data.replace(r'$32^{3}\times_64\_\quad_m_{\pi}=0.699_GeV_\_-def-\_<Q>$_','mpi700')
    data = data.replace(r'$28^{3}\times_56\_\quad_m_{\pi}=0.659_GeV_\_-def-\_<Q>$_','L28')
    data = data.replace(r'$16^{3}\times_32\_\quad_m_{\pi}=0.640_GeV_\_-def-\_<Q>$_','L16')
    data = data.replace(r'$20^{3}\times_40\_\quad_m_{\pi}=0.646_GeV_\_-def-\_<Q>$_','L20')
    data = data.replace(r'$24^{3}\times_48\_\quad_m_{\pi}=nan_GeV_\_-def-\_<Q>$_','qu_L24')
    for ival in ['-a-','-b-','-1-','-2-','-3-','-4-']:
        data = data.replace(ival,'stream_'+ival[1])
    for i in range(10000):
        data = data.replace('<'+str(i)+'>','<cfg_'+str(i)+'>')
        data = data.replace('</'+str(i)+'>','</cfg_'+str(i)+'>')
    data = data.replace(r'\alpha','alpha')
    data = data.replace(r'0.568','568')
    data = data.replace(r'0.699','699')
    data = data.replace(r'0.640','640')
    data = data.replace(r'0.646','646')
    data = data.replace(r'0.659','659')
    data = data.replace(r'0.410','410')
    data = data.replace('^','_pow_')
    data = data.replace('{','')
    data = data.replace('}','')
    data = data.replace(r'&lt;','_lt_')
    data = data.replace(r'&gt;','_gt_')
    data = data.replace(r'<Q>','Q')
    # data = data.replace(r'=','_eq_')
    # data = data.replace(r'-','_dash_')
    # data = data.replace(r'<$',r'<pie')
    # data = data.replace(r'</$',r'</pie')
    # data = data.replace(r'$_>',r'_>')
    # data = data.replace(r'\times',r'_times_')
    # data = data.replace(r'\quad',r'__')
    # data = data.replace(r'\pi',r'pi')
    # data = data.replace(r'\_',r'_')
    return data

def ReadXml(thisfile,thisShowRead=ShowRead):
    if thisShowRead and 'Parameters' not in thisfile: print('Reading Pickled file ', thisfile)
    # try:
    with open( thisfile, 'r' ) as xml_file:
        # fcntl.flock(xml_file, fcntl.LOCK_EX)
        parse_data = FormatXml(xml_file.read())
        # fcntl.flock(xml_file, fcntl.LOCK_UN)
    # print parse_data
    # print(parse_data)
    return RecFTDAS(parse(parse_data))
        # return RecFTDAS(parse(FormatXml(xml_file.read())))
    # except Exception as err:
    #     raise IOError(str(err)+'\n Failed loading xml file: \n'+thisfile)
