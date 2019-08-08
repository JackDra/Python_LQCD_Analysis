#!/usr/bin/env python


#IMPORT THIS FIRST, put in base import file
from Params import nboot,cfg_file_type,Wipe_All_Fits
from Params import cfundir,outputdir,cfgfmtdir
from MiscFuns import ODNested,mkdir_p,CheckClass,CombineCfgs,fmt_file_type
from XmlFormatting import AvgStdToFormat,tstr,unxmlfitr,untstr,tsumstr
from FlowOpps import FlowOp,FlowOpFullSquared
from PredefFitFuns import C3OneStateFitFunNoExp,C3MomTwoStateFitFunNoExp

from ReadBinaryCfuns import ReadFSCfunPickCHROMA,RC3Full
from BootStrapping import BootStrap
from FileIO import WriteXml,WritePickle,ReadPickleWrap,WriteExcel,ReadWithMeta,WriteWithMeta
from FileIO import Construct_File_Object
# import cPickle as pickle
# import dill as pickle
from DSDefines import DSfundict
import FitFunctions as ff
import MomParams as mp

import numpy as np
import re,glob,os
from copy import deepcopy,copy
from TwoPtCorrelators import TwoPointCorr
from TimeStuff import Timer

from collections import OrderedDict
import pandas as pa

# from FlowOpps import FlowOp




"""
How to use:

ThreePointCorr(cfglist=[], Info={})
creates class using:
 configuration list cfglist
 Information Dictionary Info (See below for what can be in Info)


.ReadAndWrite()
Does everything basically.


WARNING: only 1 stream is implemented now. only pass in config lists or have directories
         with 1 steam (-a- for example)

TODO:
also contains NJNQCorr(cfglist=[], Info={})

"""



class ThreePointCorr(object):
    """

    C3(q,tau,pp,t) uses bootstrap class
    pp  = sink momentum
    t   = source sink separation time
    q   = current insersion momentum
    tau = current insersion time

    Info is a dictionary containing information for the correlator:
    see comments below to see what is required.
    missing elements just sets field to default


    """
    def Construct_3pt_File(this_file):
        return Construct_File_Object(this_file,ThreePointCorr)


    ## TODO: create variational flag and interface with this function.

    ## array of all classes that inherit this class. Append this when adding new classes please!
    parentlist = ['RatioCorrelators.RatioCorr']

    def __init__(self, thisnboot=nboot, cfglist={}, Info={},name='',filename='',man_load_cfgs=False):

        self.thischecklist = ['nboot','nt','ism','jsm','tsink','DS','ppmom','Proj','tsrc','kud','ks']

        self.current = 0
        ## Initialising plotting parameters



        #self.C3 = [ igamma , ip , it , iconf ]
        self.C3_cfgs    = pa.DataFrame()
        '''
        self.C3_cfgs DataFrame:

        columns =   config_number , C3[igamma, ip, it]

        rows =      stream , configuration
        '''

        self.C3_Col_Names = ['current_operator','momentum','source_sink_separation']
        self.C3_Stats = pa.DataFrame()

        '''
        self.C3_Stats DataFrame:

        columns =   boot,   Avg,        Std

        rows =      current_operator,   momenta,    time separation
        '''

        # # C3boot [ igamma , ip , it ] bs
        # self.C3boot = ODNested()
        # self.C3Avg  = ODNested()
        # self.C3Std  = ODNested()
        self.C3_Fit_Col_Names = ['states','current_operator','momentum', 'time_fit_range']
        self.C3_Fit_Stats = pa.DataFrame()

        '''
        self.C3_Fit_Stats DataFrame:

        columns =   boot,   Avg,        Std

        rows =       states fitted to,current _operator, momenta, time fit range

        Avg/Std elements are pa.Series instances with:
            index = fit parameters
            values = Avg/Std values
        '''
        # # FitDict [ States_Fitted_To , Fit_range , ip , Fit_Parameter ] bs
        # self.FitDict    = ODNested()
        # self.FitDictAvg = ODNested()
        # self.FitDictStd = ODNested()

        self.TwoPtPars = []

        self.nboot = thisnboot


        if 'show_cfgs' in list(Info.keys()): self.show_cfgs = Info['show_cfgs']
        else: self.show_cfgs = True

        ## used to pick out different configuraiton/source lists.
        ## for example, if there are less flowed operators than 2 point correlators, passing 'TopCharge' will read that result
        if 'fo_for_cfgs' in list(Info.keys()): self.fo_for_cfgs = Info['fo_for_cfgs']
        else: self.fo_for_cfgs = False


        ## Info can contain fit ranges to be calculated
        ## must bein form Info['Fits'] = [ ('state#' , 'fit#-#' ) , ...]
        if 'Fits3pt' in list(Info.keys()): self.PredefFits = Info['Fits3pt']
        else: self.PredefFits = []

        ## Info can contain initial guess for fits, make sure it array of correct length
        if 'iGuessC3' in list(Info.keys()): self.iGuess = np.array(Info['iGuessC3'])
        else: self.iGuess = 'None' ## Defaults to value set for function in FitFunctions.py

        ## Source Smearing ism value ( in integer form)
        if 'ism' in list(Info.keys()): self.ism = 'ism'+str(Info['ism'])
        else: self.ism = ''

        if 'cfg_file_type' in list(Info.keys()): self.cfg_file_type = Info['cfg_file_type']
        else: self.cfg_file_type = cfg_file_type


        ## Sink Smearing type value.
        ## Distinguish from jsm since combined sink smearings is a thing.
        if 'jsm3pt' in list(Info.keys()): self.jsm = Info['jsm3pt']
        else: self.jsm = 'PoF0D1'

        ## t_source location ( in integer form)
        if 'tsrc' in list(Info.keys()): self.tsrc = 'tsrc'+str(Info['tsrc'])
        else: self.tsrc = ''

        ## t_sink location ( in integer form)
        if 'tsink' in list(Info.keys()): self.tsink = 'tsink'+str(Info['tsink'])
        else: self.tsink = ''

        # ## x_source list has been pushed to cfglist, this is depriciated
        # if 'xsrc' in Info.keys(): self.xsrc = ['xsrc'+str(isrc) for isrc in Info['xsrc']]
        # else: self.xsrc = ['xsrc1']

        ## current momentum list ( in form of [px py pz] )
        ## this now corresponds to the current insersion momentum q
        self.latparams = mp.LatticeParameters(Info=Info)
        self.nt = self.latparams.nt
        self.nxyz = self.latparams.nxyz
        self.dim_label = 'RC'+str(self.nxyz)+'x'+str(self.nt)
        self.dim_ll = str(self.nxyz)+r'^{3}\times ' + str(self.nt)

        self.latparams.LoadPickle()
        if 'qmom' in list(Info.keys()):
            if Info['qmom'][0] == 'All':
                self.qmom = self.latparams.GetMomList(list_type='str')
            else:
                self.qmom = Info['qmom']
        else:
            self.qmom = ['0 0 0']
        self.qform = [self.latparams.TOpform(ival,pref='q') for ival in self.qmom]
        self.checkmomlist = self.qform

        ## sink momentum ( in form of 'ppx ppy ppz' )
        if 'ppmom' in list(Info.keys()): self.ppmom = Info['ppmom']
        else: self.ppmom = '0 0 0'
        self.ppform = 'pp'+self.ppmom.replace(' ','')
        ## file has p###, not pp###
        self.filepp = self.ppform.replace('pp','p')


        ## kappa up down (in integer form)
        if 'kud' in list(Info.keys()): self.kud = 'kud'+str(Info['kud'])
        else: self.kud = ''

        ## kappa strange (in integer form)
        if 'ks' in list(Info.keys()): self.ks = 'ks'+str(Info['ks'])
        else: self.ks = ''

        ## these CurrOps contain a list of '_' separated parameters refering to:
        ## projectors, operator and real/complex part.\
        ## e.g. P3_g1g2_cmplx is z polarized projection of the operator g1g2, looking at the complex part
        if 'gammas' in list(Info.keys()): self.gammalist = Info['gammas']
        else:
            print('gammas not set, defaulting to vector current g4')
            self.gammalist = ['g4'] ## default to the vector current

        ## Projector tells whether to read the unpolarized (P4) or polarized (P1/P2/P3) projections
        if 'Projector' in list(Info.keys()): self.Proj = Info['Projector']
        else:
            print('Projector not set, defaulting to unpolarized projector P4')
            self.Proj = 'P4' ## default to unpolarized projector

        self.ProjNumb = self.Proj.replace('P','')

        ## DoubSing tells whether to read doublet (doub) or singlet (sing) contribution
        if 'DoubSing' in list(Info.keys()):
            if Info['DoubSing'] not in list(DSfundict.keys()):
                print('DSList (from DSDefines.py):')
                print(list(DSfundict.keys()))
                raise IOError("Doublet,Singlet type (Info['DoubSing']) not recognised")
            self.DS = Info['DoubSing']
        else: raise IOError('Please set DoubSing in Info for ThreePtCorrelators ')

        ## will attempt to find all replica streams in directory.
        ## set self.nstream to number of streams you want to find.
        if 'stream_list' in list(Info.keys()):
            if isinstance(Info['stream_list'],str):
                if Info['stream_list'].lower() != 'all':
                    self.stream_list = [Info['stream_list']]
                else:
                    self.stream_list = 'all'
            else:
                self.stream_list = [iinfo.lower() for iinfo in Info['stream_list']]
        else:
            self.stream_list = 'all'


        ## kappa strange (in integer form)
        if 'Csw' in list(Info.keys()): self.Csw = 'C'+str(Info['Csw'])
        else: self.Csw = 'C1715'
        ## kappa strange (in integer form)
        if 'Beta' in list(Info.keys()): self.Beta = 'B'+str(Info['Beta'])
        else: self.Beta = 'B1900'


        fileism = self.ism.replace('ism','sm')
        filejsm = self.jsm.replace('jsm','sm')
        self.kappafolder = 'Kud0'+self.kud.replace('kud','')+'Ks0'+self.ks.replace('ks','')
        self.kappafolder2 = 'Kud0'+self.kud.replace('kud','')[:-2]+'Ks0'+self.ks.replace('ks','')[:-2]

        self.readfilename = [self.dim_label+'_'+self.Beta+self.kappafolder+self.Csw,'_',
                            '_k'+self.kud.replace('kud','')+'_'+self.tsrc+fileism+'_'+filejsm+'GMA'
                            +self.ProjNumb+str(self.tsink)+self.filepp+self.DS+'NDer0.3cf']
        self.readfilename2 = [self.dim_label+'_'+self.Beta+self.kappafolder2+self.Csw,'_',
                            '_k'+self.kud.replace('kud','')
                            +'_'+self.tsrc+fileism+'_'+filejsm+'GMA'
                            +self.ProjNumb+str(self.tsink)+self.filepp+self.DS+'NDer0.3cf']
        self.readfilename3 = [self.dim_label+'_'+self.Beta+self.kappafolder2+self.Csw,'_',
                            '_k1382500_'+self.tsrc+fileism+'_'+filejsm+'GMA'
                            +self.ProjNumb+str(self.tsink)+self.filepp+self.DS+'NDer0.3cf']
        if os.path.isdir(cfundir+self.dim_label+'_'+self.kappafolder+'_AY'):
            self.filedir = cfundir+self.dim_label+'_'+self.kappafolder+'_AY/threept/'+fileism+filejsm+'GMA'+self.ProjNumb+str(self.tsink)+self.filepp+self.DS+'/'
        elif os.path.isdir(cfundir+self.dim_label+'_'+self.kappafolder):
            self.filedir = cfundir+self.dim_label+'_'+self.kappafolder+'/threept/'+fileism+filejsm+'GMA'+self.ProjNumb+str(self.tsink)+self.filepp+self.DS+'/'
        else:
            print(cfundir+self.dim_label+'_'+self.kappafolder)
            print('directory not found, using instead:')
            print(cfundir+self.kappafolder)
            self.filedir = cfundir+self.kappafolder+'/threept/'+fileism+filejsm+'GMA'+self.ProjNumb+str(self.tsink)+self.filepp+self.DS+'/'

        ## cfglist { configuration : [ x-source numbers ] }


        ## is needed?
        self.Info = Info
        self.cfglist = cfglist
        self.name = name
        self.filename = filename
        self.SetCustomName(string=name,stringfn=filename)

        if not man_load_cfgs:
            self.LoadCfgs(cfglist,name,filename)

    def LoadCfgs(self,cfglist,name='',filename=''):

        if len(cfglist) == 0:
            self.GetCfgList()
        else:
            self.ImportCfgList(cfglist)
        self.SetCustomName(string=name,stringfn=filename)


    def GetDSCfgList(self):
        thisInfo = deepcopy(self.Info)
        thisInfo['DoubSing'] = 'doub'
        doubclass = ThreePointCorr(thisnboot=self.nboot, cfglist=self.cfglist,
                                        Info=thisInfo,name=self.name.replace(self.Info['DoubSing'],'doub')
                                        ,filename=self.filename.replace(self.Info['DoubSing'],'doub'))

        thisInfo = deepcopy(self.Info)
        thisInfo['DoubSing'] = 'sing'
        singclass = ThreePointCorr(thisnboot=self.nboot, cfglist=self.cfglist,
                                        Info=thisInfo,name=self.name.replace(self.Info['DoubSing'],'sing')
                                        ,filename=self.filename.replace(self.Info['DoubSing'],'sing'))

        self.ImportCfgList(CombineCfgs(doubclass.C3_cfgs,singclass.C3_cfgs))


    def GetCfgList(self):
        def SplitCfgSrc(ifile):
            isrc = re.findall('_xsrc.*_k',ifile)
            icfg = re.findall('-.-00.*_xsrc',ifile)
            istream = re.findall('-.-',ifile)
            if len(isrc) == 0 or len(icfg) == 0:
                print()
                print('WARNING: file does not have xsrc or -.-####### (source or cfg number)')
                print(ifile)
                print()
                return '-ERROR-','xsrc-1',0
            isrc = isrc[0].split('_')[1]
            icfg = icfg[0].replace('_xsrc','')[5:]
            istream = istream[0]
            # return istream,isrc,icfg,int(icfg),int(isrc.replace('xsrc',''))
            return istream,isrc,icfg
        if self.DS not in ['doub','sing']:
            self.GetDSCfgList()
            return

        cfg_dict = ODNested()
        this_check = isinstance(self.stream_list,str) and self.stream_list == 'all'
        for ifile in glob.glob(self.filedir+'/*'):
            istream,ixsrc,icfg = SplitCfgSrc(ifile)
            if (this_check or istream in self.stream_list) and 'ERROR' not in istream:
                if istream not in list(cfg_dict.keys()) or icfg not in list(cfg_dict[istream].keys()):
                    cfg_dict[istream][icfg] = [ixsrc]
                else:
                    cfg_dict[istream][icfg].append(ixsrc)
        if len(list(cfg_dict.keys())) == 0:
            print('no configs found in ')
            print(self.filedir)

        cfglist,xsrclist,ilist = [],[],[]
        this_lambda = lambda ix : int(ix.replace('xsrc',''))
        for istream,stream_dict in cfg_dict.items():
            for icfg,(the_cfg,ixsrc_list) in enumerate(stream_dict.items()):
                ilist.append((istream,icfg))
                cfglist.append(the_cfg)
                ## sorting xsrc list with respect to its INTEGER value (not string sort!)
                ixsrc_int_list = list(map(this_lambda,ixsrc_list))
                xsrclist.append([ix for _,ix in sorted(zip(ixsrc_int_list,ixsrc_list))])
        ## finds cfgs and xsrcs conbinations that are in all streams.
        ##MAKE SURE TO ORDER CFG AND SRC

        if len(ilist) == 0:
            raise EnvironmentError('No configs found for \n'+self.filename)
        indicies = pa.MultiIndex.from_tuples(ilist,names=['stream','config_number'])
        cfg_df = pa.DataFrame(cfglist,columns=['configs'],index=indicies)
        cfg_df.loc[:,'xsrc_list'] = pa.Series(xsrclist,index=indicies)
        cfg_df.loc[:,'int_configs'] = pa.Series(list(map(int,cfglist)),index=indicies)
        cfg_df = cfg_df.sort_values('int_configs').sort_index(level='stream',sort_remaining=False)
        del cfg_df['int_configs']
        self.ImportCfgList(cfg_df)

    # def MakeOutputFile(self):
    #     self.name = '_'.join([self.tsrc,self.ism,self.jsm])
    #     self.HumanFile = '/'+self.kappafolder+'/G2/'+self.name+'.xml'
    #     self.PickleFile = '/'+self.kappafolder+'/G2/Pickle/'+self.name+'.py3p'

    #self.C3 = [ igamma , ip , it , iconf ]
    def WriteCfgsToFile(self,thisdir,show_timer=True):
        thisoutdir = thisdir + '/'+self.dim_label+'_'+self.kappafolder+'/'+'_'.join([self.tsrc,self.tsink,self.ism,self.jsm,self.ppform])+'/'
        # for igamma,gammadata in zip(self.gammalist,self.C3_cfgs['C3']):
        if show_timer: thistimer = Timer(linklist=[self.C3_cfgs['configs'].values,self.gammalist],name='Writing 3-ptCorrs ')
        for ((istream,iccfg),cfgdata),icfg in zip(iter(self.C3_cfgs['C3'].items()),self.C3_cfgs['configs'].values):
            for icgamma,gammadata in enumerate(cfgdata):
                igamma = self.gammalist[icgamma]
                insidedir = '_'.join([self.Proj,self.DS,igamma])
                mkdir_p(thisoutdir+'/'+insidedir)
                cfgstr = istream+'cfg'+str(icfg).zfill(5)
                thisfile = thisoutdir + '/'+insidedir + '/C3'+cfgstr+'.xml'
                dictout = ODNested()
                for iq,qdata in zip(self.qmom,gammadata):
                    for it,tdata in enumerate(qdata):
                        itstr = tstr(it+1)
                        dictout[iq][itstr] = tdata
                WriteXml(thisfile,{'Results':{insidedir+'_'+cfgstr:dictout}})
                if show_timer: thistimer.Lap((istream,icfg,igamma))

    def SetCustomName(self,string='',stringfn='',stringLL=''):
        if not hasattr(self,'stream_name'):
            self.stream_name = '-def-'
        if string == '':
            # self.name = '_'.join([self.kappafolder,self.Proj,'_'.join(self.gammalist),self.DS,
            #                       self.ppform,self.tsrc,self.tsink,self.ism,self.jsm])
            self.name = '_'.join([self.dim_label,self.kappafolder,self.stream_name,self.Proj,'_'.join(self.gammalist),'_'.join(self.qform),self.DS,self.ppform,self.tsrc,self.tsink,self.ism,self.jsm])
        else:
            self.name = string
        if stringfn == '':
            # self.filename = '_'.join([self.tsrc,self.tsink,self.ism,self.jsm])
            self.filename = '_'.join([self.dim_label,self.kappafolder,self.stream_name,self.Proj,'_'.join(self.gammalist),self.DS,
                                  self.ppform,self.qform[0],self.tsrc,self.tsink,self.ism,self.jsm])
        else:
            self.filename = stringfn
        if stringLL == '':
            self.LegLab = '$'+'\ '.join([self.dim_ll,self.latparams.GetPionMassLab(),self.stream_name,self.DS,self.ism,self.jsm,self.qform[0]])+'$' ## customise this how you like
        else:
            self.LegLab = stringLL
        self.SubDir = '_'.join([self.Proj,self.DS,self.ppform])
        if isinstance(self.fo_for_cfgs,str):
            self.G3dir = outputdir + '/'+self.dim_label+self.kappafolder+'/G3/'+self.fo_for_cfgs+'/'+self.SubDir+'/'
            self.G3Cfgdir = cfgfmtdir + '/'+self.dim_label+self.kappafolder+'/G3/'+self.fo_for_cfgs+'/'+self.SubDir+'/'
        else:
            self.G3dir = outputdir + '/'+self.dim_label+self.kappafolder+'/G3/'+self.SubDir+'/'
            self.G3Cfgdir = cfgfmtdir + '/'+self.dim_label+self.kappafolder+'/G3/'+self.SubDir+'/'
        mkdir_p(self.G3dir+'/Pickle/')
        mkdir_p(self.G3dir+'/Excel/')
        mkdir_p(self.G3Cfgdir)
        self.HumanFile = self.G3dir+self.filename+'.xml'
        self.ExcelFile = self.G3dir+'/Excel/'+self.filename+'.xlsx'
        self.PickleFile = self.G3dir+'/Pickle/'+self.filename+'.py3p'
        self.CfgFile = self.G3Cfgdir+'/'+self.filename+'_'+self.qform[0]+'_'+'_'.join(self.gammalist)+'.cfgs'


    def ImportCfgList(self,cfgsrclist,stream_list = 'PreDefine'):
        if not isinstance(cfgsrclist,pa.DataFrame):
            if not isinstance(stream_list,str):
                self.stream_list = stream_list
            ## cfgsrclist is [cfglist , xsrclist]
            ## cfglist must be 2-d array, first dimension is stream, second dimension is configuraiton number
            ## xsrclist is 3-d array, first dimension is stream, second dimension is cfg number, third is xsrc number
            ## xsrclist must have values 'xsrc##' with no zero pre buffering
            ilist = []
            xsrc_list = []
            for istream,(is_cfg,is_xsrc) in enumerate(zip(cfgsrclist[0],cfgsrclist[1])):
                for icf in range(is_cfg):
                    xsrc_list.append(is_xsrc)
                    ilist.append((istream,is_cfg))
            indicies = pa.MultiIndex.from_tuples(ilist,names=['stream','config_number'])
            cfg_df = pa.DataFrame(np.array(cfgsrclist[0]).flatten(),columns=['configs'],index=indicies)
            cfg_df.loc[:,'xsrc_list'] = pa.Series(xsrc_list,index=indicies)
            self.C3_cfgs = cfg_df
        else:
            self.C3_cfgs = cfgsrclist
        if 'configs' in self.C3_cfgs:
            self.C3_cfgs.loc[:,'stream_cfgs'] = pa.Series(['00'.join((istream,icfg)) for (istream,iccfg),icfg in self.C3_cfgs['configs'].items()],index=self.C3_cfgs.index)
            self.Update_Stream_List()
            self.ncfg_list = [icfg.size for istream,icfg in self.C3_cfgs['configs'].groupby(level='stream')]
            self.xsrcAvg = np.mean(list(map(len,self.C3_cfgs['xsrc_list'].values)))
            self.nMeas = self.xsrcAvg*sum(self.ncfg_list)
            self.cfglist = self.C3_cfgs

    def Update_Stream_List(self):
        if 'configs' in self.C3_cfgs:
            self.stream_list = self.C3_cfgs.index.levels[0]
            self.nstream = len(self.stream_list)
            if hasattr(self,'stream_name'):
                old_sn = self.stream_name
                self.stream_name = '-'+''.join([istr.replace('-','') for istr in self.stream_list])+'-'
                if hasattr(self,'filename'):
                    self.filename = self.filename.replace(old_sn,self.stream_name)
                if hasattr(self,'name'):
                    self.name = self.name.replace(old_sn,self.stream_name)
                if hasattr(self,'LegLab'):
                    self.LegLab = self.LegLab.replace(old_sn,self.stream_name)
            else:
                self.stream_name = '-'+''.join([istr.replace('-','') for istr in self.stream_list])+'-'

    def Stats(self):
        if 'boot' not in self.C3_Stats: return
        lAvg,lStd = [],[]
        for (igamma,ip,it),iC3 in self.items():
            iC3.Stats()
            lAvg.append(iC3.Avg)
            lStd.append(iC3.Std)
        indicies = pa.MultiIndex.from_tuples(self.C3_Stats.index,names=self.C3_Col_Names)
        self.C3_Stats.loc[:,'Avg'] = pa.Series(lAvg,index=indicies)
        self.C3_Stats.loc[:,'Std'] = pa.Series(lStd,index=indicies)



    def RemoveVals(self):
        if 'C3' in self.C3_cfgs: del self.C3_cfgs['C3']
        if 'C2' in self.C3_cfgs: del self.C3_cfgs['C2']

    ## self.C3 = [ igamma , ip, it , iconf ]
    def Read(self,show_timer=False,Full=False,file_type='PreDef'):
        if self.DS not in ['doub','sing']: return
        if self.Read_Cfgs(file_type=file_type,show_timer=show_timer,Full=Full): return
        ## making interpolator number more human readable
        thisdata = []
        prop_src = []
        ## this part isnt really nessisary, maybe later when we implement operator overloaing more precisly.
        # if len(self.filename.split('_')) > 4:
        #     raise IOError('Do not read after combining correlators, please recreate the correlator')
        fix_file_size = 'RC32x64' in self.dim_label and '13754' in self.kud and '1364' in self.ks
        if show_timer: thistimer = Timer(linklist=self.C3_cfgs['stream_cfgs'].values,name='Reading 3-ptCorrs '+ self.name.replace('_'.join(self.qform),''))
        for i_stream_cfg,ixsrc_list in zip(self.C3_cfgs['stream_cfgs'].values,self.C3_cfgs['xsrc_list'].values):
            thisxsrclist = [self.filedir+ self.readfilename[0]+i_stream_cfg+self.readfilename[1]+ixsrc+self.readfilename[2] for ixsrc in ixsrc_list]
            if not os.path.isfile(thisxsrclist[0]):
                thisxsrclist = [self.filedir+ self.readfilename2[0]+i_stream_cfg+self.readfilename2[1]+ixsrc+self.readfilename2[2] for ixsrc in ixsrc_list]
            if not os.path.isfile(thisxsrclist[0]):
                thisxsrclist = [self.filedir+ self.readfilename3[0]+i_stream_cfg+self.readfilename3[1]+ixsrc+self.readfilename3[2] for ixsrc in ixsrc_list]
            if Full:
                idata = RC3Full(thisxsrclist,self.qmom,self.gammalist,forcent=self.nt,check_fs=fix_file_size)
                thisdata.append(idata.data)
                prop_src.append([ifinfo['tsource'] for ifinfo in idata.f_info])
            else:
                thisdata.append(ReadFSCfunPickCHROMA(thisxsrclist,self.qmom,self.gammalist,forcent=self.nt,check_fs=fix_file_size).data)
            if show_timer: thistimer.Lap()
        if Full: self.C3_cfgs.loc[:,'tsrc_list'] = pa.Series(prop_src,index=self.C3_cfgs.index)
        self.C3_cfgs.loc[:,'C3'] = pa.Series(thisdata,index=self.C3_cfgs.index)
        self.Write_Cfgs(show_timer=show_timer,file_type=file_type)

    def Bootstrap(self,WipeData=True,DefWipe=False,show_timer=False):
        if self.DS not in ['doub','sing']:
            try:
                self.LoadDS(WipeData=WipeData,DefWipe=DefWipe,show_timer=show_timer)
            except Exception as err:
                print('reading doublet and singlet failed, reading with wipeing now')
                self.LoadDS(WipeData=WipeData,DefWipe=True,show_timer=show_timer)
        else:
            if 'C3' not in self.C3_cfgs: self.Read()
            lC3,lC3Avg,lC3Std = [],[],[]
            ilist = []
            rlist = None
            if show_timer: thistimer = Timer(   linklist=[self.gammalist,self.qform,list(range(self.latparams.nt))],
                                                name='BS C3 ')
            for icgamma,igamma in enumerate(self.gammalist):
                for icp,ip in enumerate(self.qform):
                    for it in range(self.latparams.nt):
                        strit = tstr(it+1)
                        ilist.append((igamma,ip,strit))
                        this_lambda = lambda x : x[icgamma][icp][it]
                        tdata = self.C3_cfgs['C3'].apply(this_lambda)
                        thisboot = BootStrap(   self.nboot, name='G_{3}('+','.join([igamma,ip,strit])+')',
                                                cfgvals=tdata,rand_list=rlist)
                        rlist = thisboot.Get_Rand_List()
                        lC3.append(thisboot)
                        lC3Avg.append(thisboot.Avg)
                        lC3Std.append(thisboot.Std)
                        if show_timer: thistimer.Lap((igamma,ip,strit))
            if len(ilist) > 0:
                indicies = pa.MultiIndex.from_tuples(ilist,names=self.C3_Col_Names)
                self.C3_Stats.loc[:,'boot'] = pa.Series(lC3,index=indicies)
                self.C3_Stats.loc[:,'Avg'] = pa.Series(lC3Avg,index=indicies)
                self.C3_Stats.loc[:,'Std'] = pa.Series(lC3Std,index=indicies)
            if WipeData: self.RemoveVals()




    # def ReadAndBoot(self):
    #     self.Read()
    #     self.Bootstrap()

    def ImportFitRanges(self,fit_ranges):
        self.fit_ranges = fit_ranges
        if isinstance(self.fit_ranges,str):
            self.fit_ranges = [self.fit_ranges]


    ## you can either just pass in the two point correlator class instance
    ## or pass in a bootstrapped parameters list
    ## or just pass in average values
    def ImportTwoPtFits(self,state,sourcepars,sinkpars):
        self.TwoPtPars = []
        for iq in self.qform:
            ip = iq.replace('q','p')
            self.TwoPtPars.append(np.append(list(self.GetTwoPtPars(ip,state,sourcepars).values()) ,list(self.GetTwoPtPars(ip,state,sinkpars).values())))


    def GetTwoPtPars(self,ip,state,pars,fitr='first'):
        if isinstance(pars,TwoPointCorr):
            if isinstance(state,int) or 'state' not in state:
                thisstate = 'state'+str(state)


            if fitr == 'first' and (state,ip) in pars.C2_Fit_Stats.index:
                fitr_list = list(pars.C2_Fit_Stats['boot'][thisstate,ip].keys())
                if len(fitr_list) > 1:
                    print('More than 1 fit parameter found in two point correlator fit parameters, ')
                    print('choosing ',fitr_list[0])
                return pars.C2_Fit_Stats['boot'][thisstate,ip,fitr_list[0]].fit_data['Params']
            elif (state,ip,fitr) in pars.C2_Fit_Stats.index:
                return pars.C2_Fit_Stats['boot'][thisstate,ip,fitr].fit_data['Params']
        elif isinstance(pars,ff.Fitting):
            return pars.fit_data['Params']
        elif isinstance(pars,list) or isinstance(pars,np.ndarray):
            return pa.Series([BootStrap(self.nboot,bootvals=[np.float64(ipar)]*self.nboot) for ipar in pars])
        else:
            raise IOError('two point fit parameters not recognised')

    ## fit_ranges is list of fit ranges to fit the correlator to
    ## formatted as fitr#-# where # are the min and max fit values.
    ## states correspond to how many states to fit to (currenlty only working for states = 1 or 2)
    ## C2pars is list of size 2, containing source and sink fit parameters used in the three point correlator fits
    ## C2pars can be class instances of TwoPointCorr, Fitting, or just a list of values.
    ## Word of caution, the bootstraps need to have been constructed with the same configurations.
    ## Best to have this run from RatioFunction class only.
    def Fit(self,C2pars='PreDef',state='PreDef',fit_ranges='PreDef',iGuess='PreDef',EstDir=False,WipeFit=True):
        def RunPredefFits():
            if state == 'PreDef' and fit_ranges == 'PreDef':
                for izip in self.PredefFits:
                    if len(izip) == 2:
                        self.Fit(C2pars=C2pars,state=izip[0],fit_ranges=[izip[1]],EstDir=EstDir,WipeFit=WipeFit)
                    else:
                        self.Fit(C2pars=C2pars,state=izip[0],fit_ranges=[izip[1]],iGuess=izip[2],EstDir=EstDir,WipeFit=WipeFit)
            elif state == 'PreDef':
                for izip in self.PredefFits:
                    if izip[1] == fit_ranges:
                        if len(izip) == 2:
                            self.Fit(C2pars=C2pars,state=izip[0],fit_ranges=[izip[1]],EstDir=EstDir,WipeFit=WipeFit)
                        else:
                            self.Fit(C2pars=C2pars,state=izip[0],fit_ranges=[izip[1]],iGuess=izip[2],EstDir=EstDir,WipeFit=WipeFit)
            elif fit_ranges == 'PreDef':
                for izip in self.PredefFits:
                    if izip[0] == state:
                        if len(izip) == 2:
                            self.Fit(C2pars=C2pars,state=izip[0],fit_ranges=[izip[1]],EstDir=EstDir,WipeFit=WipeFit)
                        else:
                            self.Fit(C2pars=C2pars,state=izip[0],fit_ranges=[izip[1]],iGuess=izip[2],EstDir=EstDir,WipeFit=WipeFit)
        if 'PreDef' in [state,fit_ranges]:
            RunPredefFits()
            return
        if iGuess != 'PreDef': self.iGuess = iGuess
        if isinstance(state,str):
            if 'state' in state:
                state = state.replace('state','')
            state = int(state)

        if C2pars != 'PreDef': self.ImportTwoPtFits(state,C2pars[0],C2pars[1])
        self.ImportFitRanges(fit_ranges)
        if state == 1:
            # thisfitfun = C3OneStateFitFunNoExp
            thisfitfun = C3OneStateFitFunNoExp
            npar = 1
        elif state == 2:
            # thisfitfun = C3TwoStateFitFunNoExp
            thisfitfun = C3MomTwoStateFitFunNoExp
            npar = 4
        elif state > 2:
            raise IOError('fitting state is only implemented up to 2 state so far. Cannot do '+str(state))
        ststr = 'state'+str(state)
        legstate = str(state)+'SF'
        lFit,lFitA,lFitS = [],[],[]
        ilist = []
        for igamma,gammadata in self.C3_Stats['boot'].groupby(level='current_operator'):
            for (ip,pdata),ptwopt in zip(gammadata.groupby(level='momentum'),self.TwoPtPars):
                # print 'Fitting state',str(state) , ' ' , ifit
                for ifit in self.fit_ranges:
                    if (ststr,igamma,ip,ifit) in self.C3_Fit_Stats.index:
                        if WipeFit:
                            self.C3_Fit_Stats.drop([(ststr,igamma,ip,ifit)],inplace=True)
                        else:
                            continue
                    ilist.append((ststr,igamma,ip,ifit))
                    ifitmin,ifitmax = list(map(int,unxmlfitr(ifit)))
                    # print '     Fitting ', ip
                    ## is this the right ordering of dictionary dimensions?
                    tdatarange = [[]]
                    ydata = []
                    for it,tdata in pdata.items():
                        tint = untstr(it)
                        if ifitmin <= tint <= ifitmax:
                            tcurrboot = BootStrap(self.nboot,bootvals=[np.float64(tint)]*self.nboot)
                            tdatarange[0].append(tcurrboot)
                            ydata.append(tdata)
                    thistlen = len(tdatarange[0])
                    ## this + 3 is just calabrated from looking at the conserved current two- and three- point correlators. Not sure why exactly.
                    tsinkboot = BootStrap(self.nboot,bootvals=[np.float64(self.tsink.replace('tsink',''))+3]*self.nboot)
                    tdatarange.append([tsinkboot]*thistlen)
                    for iparam in ptwopt:
                        tdatarange.append([iparam]*thistlen)

                    # print 'DEBUG'
                    # for it,iy in zip(np.rollaxis(np.array(tdatarange),1),ydata):
                    #     print it[0].Avg,it[1].Avg,it[2].Avg,it[3].Avg,it[4].Avg,it[5].Avg,iy.Avg

                    if EstDir:
                        thisff = ff.Fitting(Funs=[thisfitfun,npar,'Estimate'],data=[tdatarange,ydata],name=' '.join([legstate,ifit]),iGuess=self.iGuess)
                    else:
                        thisff = ff.Fitting(Funs=[thisfitfun,npar],data=[tdatarange,ydata],name=' '.join([legstate,ifit]),iGuess=self.iGuess)
                    thistest = thisff.FitBoots()
                    if not np.isnan(thistest):
                        lFit.append(thisff)
                        lFitA.append(thisff.fit_data['ParamsAvg'])
                        lFitS.append(thisff.fit_data['ParamsStd'])
        indicies = pa.MultiIndex.from_tuples(ilist,names=self.C3_Fit_Col_Names)
        if 'boot' in self.C3_Fit_Stats.columns:
            this_df = pa.DataFrame()
            this_df.loc[:,'boot'] = pa.Series(lFit,index=indicies)
            this_df.loc[:,'Avg'] = pa.Series(lFitA,index=indicies)
            this_df.loc[:,'Std'] = pa.Series(lFitS,index=indicies)
            self.C3_Fit_Stats = self.C3_Fit_Stats.append(this_df)
        else:
            self.C3_Fit_Stats.loc[:,'boot'] = pa.Series(lFit,index=indicies)
            self.C3_Fit_Stats.loc[:,'Avg'] = pa.Series(lFitA,index=indicies)
            self.C3_Fit_Stats.loc[:,'Std'] = pa.Series(lFitS,index=indicies)

                    # self.FitDict[ststr][igamma][ip][ifit] = thisff
                    # self.FitDictAvg[ststr][igamma][ip][ifit] = thisff.ParamsAvg
                    # self.FitDictStd[ststr][igamma][ip][ifit] = thisff.ParamsStd
                    #
                    # if state == 1:
                    #     self.FitDictAvg[ststr][igamma][ip][ifit]['A_Ep'] = ptwopt[0].Avg
                    #     self.FitDictStd[ststr][igamma][ip][ifit]['A_Ep'] = ptwopt[0].Std
                    #     self.FitDictAvg[ststr][igamma][ip][ifit]['Ep'] = ptwopt[1].Avg
                    #     self.FitDictStd[ststr][igamma][ip][ifit]['Ep'] = ptwopt[1].Std
                    #     self.FitDictAvg[ststr][igamma][ip][ifit]['A_Epp'] = ptwopt[2].Avg
                    #     self.FitDictStd[ststr][igamma][ip][ifit]['A_Epp'] = ptwopt[2].Std
                    #     self.FitDictAvg[ststr][igamma][ip][ifit]['Epp'] = ptwopt[3].Avg
                    #     self.FitDictStd[ststr][igamma][ip][ifit]['Epp'] = ptwopt[3].Std
                    # elif state == 2:
                    #     self.FitDictAvg[ststr][igamma][ip][ifit]['A_Ep'] = ptwopt[0].Avg
                    #     self.FitDictStd[ststr][igamma][ip][ifit]['A_Ep'] = ptwopt[0].Std
                    #     self.FitDictAvg[ststr][igamma][ip][ifit]['Ep'] = ptwopt[1].Avg
                    #     self.FitDictStd[ststr][igamma][ip][ifit]['Ep'] = ptwopt[1].Std
                    #     self.FitDictAvg[ststr][igamma][ip][ifit]['Ap_Ep'] = ptwopt[2].Avg
                    #     self.FitDictStd[ststr][igamma][ip][ifit]['Ap_Ep'] = ptwopt[2].Std
                    #     self.FitDictAvg[ststr][igamma][ip][ifit]['DEp'] = ptwopt[3].Avg
                    #     self.FitDictStd[ststr][igamma][ip][ifit]['DEp'] = ptwopt[3].Std
                    #     self.FitDictAvg[ststr][igamma][ip][ifit]['A_Epp'] = ptwopt[4].Avg
                    #     self.FitDictStd[ststr][igamma][ip][ifit]['A_Epp'] = ptwopt[4].Std
                    #     self.FitDictAvg[ststr][igamma][ip][ifit]['Epp'] = ptwopt[5].Avg
                    #     self.FitDictStd[ststr][igamma][ip][ifit]['Epp'] = ptwopt[5].Std
                    #     self.FitDictAvg[ststr][igamma][ip][ifit]['Ap_Epp'] = ptwopt[6].Avg
                    #     self.FitDictStd[ststr][igamma][ip][ifit]['Ap_Epp'] = ptwopt[6].Std
                    #     self.FitDictAvg[ststr][igamma][ip][ifit]['DEpp'] = ptwopt[7].Avg
                    #     self.FitDictStd[ststr][igamma][ip][ifit]['DEpp'] = ptwopt[7].Std
        self.Write()



    def Write(self):
        def FixDictArray(thislist,ikey):
            outDict = OrderedDict()
            for ic,ilist in enumerate(thislist):
                outDict[ikey+str(ic+1)] = ilist
            return outDict

        if not hasattr(self,'show_cfgs'):
            self.show_cfgs = True

        outDict = ODNested()
        outDict['name'] = self.name
        outDict['kud'] = self.kud
        outDict['ks'] = self.ks
        outDict['nt'] = self.nt
        outDict['nxyz'] = self.nxyz

        for istream,incfg in zip(self.stream_list,self.ncfg_list):
            outDict['ncfg_'+istream] = incfg
        outDict['nxsrc_avg'] = self.xsrcAvg
        outDict['nMeas'] = self.nMeas
        outDict['nboot'] = self.nboot
        outDict['tsrc'] = self.tsrc
        outDict['tsink'] = self.tsink
        # outDict['xsrc'] = FixDictArray(self.xsrc,'xsrc')
        outDict['ism'] = self.ism
        outDict['sink_type'] = self.jsm
        outDict['Projector'] = self.Proj
        outDict['Doub_Sing'] = self.DS
        outDict['ppmom'] = self.ppmom
        for icg,igamma in enumerate(self.gammalist):
            outDict['gamma'+str(icg+1)] = igamma
        for icq,iqmom in enumerate(self.qmom):
            outDict['qmom'+str(icq+1)] = iqmom

        excel_params = pa.Series(deepcopy(outDict))
        for icg,igamma in enumerate(self.gammalist):
            del outDict['gamma'+str(icg+1)]
        for icq,iqmom in enumerate(self.qmom):
            del outDict['qmom'+str(icq+1)]

        outDict['gamma_list'] = FixDictArray(self.gammalist,'gamma')
        outDict['qmom_list'] = FixDictArray(self.qmom,'qmom')


        for col_key,col_data in self.C3_Stats.items():
            if 'Avg' in col_key or 'Std' in col_key: continue
            for (igamma,ip,it),idata in col_data.items():
                outDict[col_key][igamma][ip][it] = AvgStdToFormat(idata.Avg,idata.Std)


        if 'boot' in self.C3_Fit_Stats:
            for (istate,igamma,ip,itr),fitdict in self.C3_Fit_Stats['boot'].items():
                if not isinstance(fitdict,ff.Fitting): continue
                outDict['Fits'][istate][igamma][ip][itr] = fitdict.GetOutputDict()


        # for ststr,statedict in self.FitDict.iteritems():
        #     if len(statedict.keys()) == 0: continue
        #     for igamma,gammadict in statedict.iteritems():
        #         for ip,pdict in gammadict.iteritems():
        #             for ifit,fitdict in pdict.iteritems():
        #                 thisfitAvg,thisfitStd = self.FitDictAvg[ststr][igamma][ip][ifit],self.FitDictStd[ststr][igamma][ip][ifit]
        #                 for (ipar,paravg),parstd in zip(thisfitAvg.iteritems(),thisfitStd.itervalues):
        #                     thisfmtflag = 'f'
        #                     if 'A_Ep' == ipar or 'A_Epp' == ipar: thisfmtflag = 'e'
        #                     outDict['Fits'][ststr][igamma][ip][ifit][ipar] = AvgStdToFormat(paravg,parstd,frmtflag=thisfmtflag)
        #                 outDict['Fits'][ststr][igamma][ip][ifit]['Chi^2_pdf'] = AvgStdToFormat(fitdict.Chi2DoF.Avg,fitdict.Chi2DoF.Std)
        #                 outDict['Fits'][ststr][igamma][ip][ifit]['MaxIter'] = fitdict.MaxIter
        #                 outDict['Fits'][ststr][igamma][ip][ifit]['Prec'] = fitdict.Prec

        if any(np.array(self.ncfg_list) > 0) and self.show_cfgs:
            for ((istream,iccfg),icfg),ixsrc_list in zip(iter(self.C3_cfgs['configs'].items()),self.C3_cfgs['xsrc_list'].values):
                outDict['cfglist'][istream][icfg] = 'n_xsrc_'+str(len(ixsrc_list))
        #
        # if self.ncfg > 0 and self.show_cfgs:
        #     for icfg,xsrclist in self.cfglist.iteritems():
        #         outDict['cfglist'][icfg] = FixDictArray(xsrclist,'xsrc')


        ## Write human readable file with data
        WriteXml(self.HumanFile,{'Results':outDict})


        if any(np.array(self.ncfg_list) > 0):
            for ((istream,iccfg),icfg),ixsrc_list in zip(iter(self.C3_cfgs['configs'].items()),self.C3_cfgs['xsrc_list'].values):
                outDict['cfglist'][istream][icfg] = 'n_xsrc_'+str(len(ixsrc_list))

        ## TODO, include fitting in the excell file
        WriteExcel(self.ExcelFile,{'C3_results':deepcopy(self.C3_Stats)},cfglist=self.C3_cfgs,params=excel_params)


        # if self.ncfg > 0:
        #     for icfg,xsrclist in self.cfglist.iteritems():
        #         outDict['cfglist'][icfg] = FixDictArray(xsrclist,'xsrc')

        ## pickles rest of data for reading
        WritePickle(self.PickleFile,self.__dict__)

    ##################################

    def LoadDS(self,DefWipe=False,WipeData=True,CheckCfgs=True,CheckMom=True,show_timer=False):
        thisInfo = deepcopy(self.Info)
        thisInfo['DoubSing'] = 'doub'
        self.doubclass = ThreePointCorr(thisnboot=self.nboot, cfglist=self.cfglist, Info=thisInfo,
                                        name=self.name.replace(self.Info['DoubSing'],'doub'),
                                        filename=self.filename.replace(self.Info['DoubSing'],'doub'))

        thisInfo = deepcopy(self.Info)
        thisInfo['DoubSing'] = 'sing'
        self.singclass = ThreePointCorr(thisnboot=self.nboot, cfglist=self.cfglist, Info=thisInfo,
                                        name=self.name.replace(self.Info['DoubSing'],'sing'),
                                        filename=self.filename.replace(self.Info['DoubSing'],'sing'))


        # self.doubclass.GetCfgList()
        # self.singclass.GetCfgList()
        if CheckCfgs:
            comb_cfgs = CombineCfgs(self.doubclass.C3_cfgs,self.singclass.C3_cfgs)
            self.doubclass.ImportCfgList(copy(comb_cfgs))
            self.singclass.ImportCfgList(copy(comb_cfgs))
        self.doubclass.LoadPickle(  DefWipe=DefWipe,WipeData=WipeData,
                                    CheckCfgs=CheckCfgs,CheckMom=CheckMom,show_timer=show_timer)
        self.singclass.LoadPickle(  DefWipe=DefWipe,WipeData=WipeData,
                                    CheckCfgs=CheckCfgs,CheckMom=CheckMom,show_timer=show_timer)
        if self.DS not in list(DSfundict.keys()):
            raise EnvironmentError('Doublet Singlet combination' + str(self.DS)+' not reconised, add to DSDefines.py')

        outlist,ilist = [],[]
        for (igamma,ip,it),idoubdata in self.doubclass.items():
            if (igamma,ip,it) in self.singclass.C3_Stats.index:
                ilist.append((igamma,ip,it))
                outlist.append(DSfundict[self.DS](idoubdata,self.singclass.C3_Stats['boot'][igamma,ip,it]))
        iindex = pa.MultiIndex.from_tuples(ilist,names=self.C3_Col_Names)
        self.C3_Stats = pa.Series(outlist,index=iindex).to_frame(name='boot')
        # self.Stats()
        if DefWipe:
            self.doubclass = None
            self.singclass = None
        self.Stats()

    def ReadAndWrite(self,DefWipe=False,WipeData=True,show_timer=False):
        if self.DS not in ['doub','sing']:
            try:
                self.LoadDS(WipeData=WipeData,DefWipe=DefWipe,show_timer=show_timer)
            except Exception as err:
                print('reading doublet and singlet failed, reading with wipeing now')
                print(self.filename)
                self.LoadDS(WipeData=WipeData,DefWipe=True,show_timer=show_timer)
            self.Bootstrap(WipeData=WipeData,DefWipe=DefWipe,show_timer=show_timer)
        else:
            self.Read(show_timer=show_timer)
            self.Bootstrap(WipeData=WipeData,DefWipe=DefWipe,show_timer=show_timer)
        # self.Fit()
        self.Write()

    def FixRead(self,readdict):
        readdict['Info']['outdir'] = self.Info['outdir']
        readdict['latparams'].outputdir = self.latparams.outputdir
        readdict['latparams'].PickleFile = self.latparams.PickleFile
        readdict['latparams'].HumanFile = self.latparams.HumanFile
        readdict['PickleFile'] = self.PickleFile
        readdict['ExcelFile'] = self.ExcelFile
        readdict['HumanFile'] = self.HumanFile
        readdict['CfgFile'] = self.CfgFile
        if 'cfg_file_type' not in list(readdict.keys()):
            readdict['cfg_file_type'] = self.cfg_file_type
        return readdict

    def LoadPickle(self,DefWipe=False,WipeData = True,CheckCfgs=False,CheckMom=True,show_timer=False):
        if os.path.isfile(self.PickleFile) and not DefWipe:
            # print 'Loading Pickle for ' , self.name
            loadeddict = ReadPickleWrap(self.PickleFile)
            loadeddict = self.FixRead(loadeddict)
            checklist = self.thischecklist
            if CheckCfgs: checklist = checklist + ['C3_cfgs']
            if CheckMom: checklist = checklist + ['checkmomlist']
            if CheckClass(self.__dict__,loadeddict,checklist):
                self.__dict__.update(loadeddict)
                self.read_cfgs_from_file = True
                if Wipe_All_Fits:
                    self.C3_Fit_Stats = pa.DataFrame()
                # if len(self.C3_Fit_Stats.values)> 0 and not isinstance(self.C3_Fit_Stats.values,sff.SetOfFitFuns):
                #     self.C3_Fit_Stats = pa.DataFrame()
            else:
                if os.path.isfile(self.PickleFile+'.bak'):
                    print('using backupfile for '+self.name)
                    self.PickleFile = self.PickleFile+'.bak'
                    self.HumanFile = self.HumanFile+'.bak'
                    self.ExcelFile = self.ExcelFile.replace('.xlsx','.bak.xlsx')
                    self.LoadPickle(DefWipe=DefWipe,WipeData=WipeData,CheckCfgs=CheckCfgs,CheckMom=CheckMom)
                    return
                print('Warning, file ' , self.HumanFile , ' has different parameters than this instance, reading and writing over file now')
                print()
                self.ReadAndWrite(WipeData=WipeData,DefWipe=DefWipe,show_timer=show_timer)
            # self.Write()
        else:
            self.ReadAndWrite(WipeData=WipeData,DefWipe=DefWipe,show_timer=show_timer)


    def Write_Cfgs(self,file_type='PreDef',show_timer=True):
        if hasattr(self,'read_cfgs_from_file') and self.read_cfgs_from_file:
            return
        if file_type == 'PreDef':
            file_type = self.cfg_file_type
        if 'to_' not in file_type:
            file_type = 'to_'+file_type
        this_file = self.CfgFile+fmt_file_type(file_type.replace('to_',''))
        if show_timer: thistimer = Timer(name='Writing Cfgs '+this_file)
        if len(self.C3_cfgs.index) == 0:
            self.Read(show_timer=show_timer)
        # this_df,dump = self.Get_Flattened_Cfgs()
        this_meta = [self.gammalist,self.qmom,self.nt]
        WriteWithMeta(self.C3_cfgs,this_meta,this_file,file_type=file_type)
        if show_timer: thistimer.Stop()

    def ChopDepVarList(self,this_df,read_meta,show_err=True):
        isin = True
        name_list = ['gamma_list','momentum_list']
        index = [slice(None),slice(None)]
        if self.nt != read_meta[2]:
            if show_err:
                print('\n missmatch in euclidian time, read_nt='+str(read_meta[2])+', setup_nt='+str(self.nt))
            return this_df,False
        gamma_list = read_meta[0]
        mom_list = read_meta[1]
        for ic,(ithis_list,icheck_list) in enumerate([(self.gammalist,gamma_list),(self.qmom,mom_list)]):
            isin = isin and set(ithis_list).issubset(icheck_list)
            if isin:
                same_list = len(ithis_list) == len(icheck_list)
                if not same_list:
                    index_list = [np.where(np.array(icheck_list) == it)[0][0] for it in ithis_list]
                    index[ic] = index_list
            elif show_err:
                print('\n missmatch in ',name_list[ic])
                print(ithis_list)
                print(icheck_list)
                break
        if isin:
            def Chop_fun(val):
                return np.array(val)[tuple(index)].tolist()
            this_df.loc[:,'C3'] = this_df.loc[:,'C3'].apply(Chop_fun)
        return this_df,isin




    def Read_Cfgs(self,file_type='PreDef',show_timer=True,CheckCfgs=True,Full=False):
        def CheckCfgList(cfg_data,read_meta):
            if cfg_data is None:
                return None,False
            if CheckCfgs and 'stream_cfgs' in self.C3_cfgs:
                second_class = pa.Series([cfg_data],index=['C3_cfgs'])
                cfg_bool = CheckClass(self.__dict__,second_class,['C3_cfgs'])
                # cfg_bool = np.array(self.C3_cfgs['stream_cfgs'].values) == np.array(cfg_data['stream_cfgs'].values)
                # if isinstance(cfg_bool,np.ndarray):
                #     cfg_bool = all(cfg_bool)
                # if show_timer and not cfg_bool:
                #     print '\n'+GetCfgMissmatch( np.array(self.C3_cfgs['stream_cfgs'].values),
                #                                 np.array(cfg_data['stream_cfgs'].values))
            else:
                cfg_bool = True
            if cfg_bool:
                cfg_data,outbool = self.ChopDepVarList(cfg_data,read_meta)
                if Full and 'tsrc_list' not in cfg_data.columns:
                    outbool = False
                return cfg_data,outbool
            else:
                return None,False

        if file_type == 'PreDef':
            file_type = self.cfg_file_type
        if 'read_' not in file_type:
            file_type = 'read_'+file_type
        this_file = self.CfgFile+fmt_file_type(file_type.replace('read_',''))
        if not os.path.isfile(this_file):
            if '-def-' in this_file:
                this_file = this_file.replace('-def-','-*-')
                file_list = glob.glob(this_file)
                if len(file_list) == 0:
                    print('Cfg FNF: \n', this_file)
                    return False
                else:
                    this_file = file_list[0]
        if not os.path.isfile(this_file):
            if show_timer: print('cfg FNF: ' + this_file)
            return False
        if show_timer: thistimer = Timer(name='Reading Cfgs '+this_file)
        meta_data,cfg_data = ReadWithMeta(this_file,file_type=file_type)
        cfg_data,check_bool = CheckCfgList(cfg_data,meta_data)
        # cfg_data = self.Format_Flattened_Cfgs(cfg_data)
        self.read_cfgs_from_file = check_bool
        if check_bool:
            self.ImportCfgList(cfg_data)
            if show_timer: thistimer = thistimer.Stop()
            self.ncfg_list = [icfg.size for istream,icfg in self.C3_cfgs['configs'].groupby(level='stream')]
            return True
        else:
            if show_timer: print('\n file incorrect,', end=' ')
            self.CfgFile = self.CfgFile.replace('.cfgs','.cfgs.bak')
            if os.path.isfile(this_file.replace('.cfgs','.cfgs.bak')):
                print(' found backup file')
                return self.Read_Cfgs(file_type=file_type,show_timer=show_timer,CheckCfgs=CheckCfgs,Full=Full)
            else:
                print(' reading unformatted')
                return False


    ## Comparisons

    def __str__(self):
        return '_'.join(self.name)


    def __iter__(self):
        return self

    def iteritems(self):
        return iter(self.C3_Stats['boot'].items())

    def items(self):
        return list(self.C3_Stats['boot'].items())

    def itemsAvgStd(self):
        outlist = []
        for ivals in self.C3_Stats.itertuples():
            outlist.append((ivals[0],ivals[2],ivals[3]))
        return outlist

    def values(self):
        return self.C3_Stats['boot'].values

    def itervalues(self):
        return iter(self.C3_Stats['boot'].values)

    def valuesAvgStd(self):
        outlist = []
        for ivals in self.C3_Stats.itertuples():
            outlist.append((ivals[2],ivals[3]))
        return outlist

    def keys(self):
        return self.C3_Stats.index


    def __setitem__(self, key, value ):
        self.C3_Stats.loc[key,'boot'] = value

    def __getitem__(self, key):
        return self.C3_Stats.loc[key,'boot']

    def __reversed__(self):
        self.C3_Stats = self.C3_Stats.iiloc[::-1]
        return self

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



    ## Operator overloading Starting to get tedious, leave to last or think of better way of implementing (maybe use predefined funtion passed in?)

    #Too Much work, only works for one operator. Please use SetCustomName after doing all algebraic manipulation


    # def UpdateName(self,opp,C3,C32):
    #     def GetOppSelf(opp):
    #         if opp == '+':
    #             return 'x2'
    #         elif opp == '-':
    #             return 'x0'
    #         elif opp == 'x':
    #             return 'pow2'
    #         elif opp == 'div':
    #             return 'divself'
    #         elif opp == 'pow':
    #             return 'powself'
    #     def GetLegSelf(opp):
    #         if opp == '+':
    #             return '\times 2'
    #         elif opp == '-':
    #             return '\times 0'
    #         elif opp == 'x':
    #             return '^{2}'
    #         elif opp == 'div':
    #             return 'divself'
    #         elif opp == 'pow':
    #             return 'powself'

    #     def LLFormatOp(opp):
    #         oppout = opp.replace('x',r'\times ')
    #         oppout = oppout.replace('pow','^{}')
    #         return oppout

    #     outname = []
    #     outLL = []
    #     oppLL = LLFormatOp(opp) ## \times looks better :)
    #     if isinstance(C32,ThreePointCorr) and isinstance(Op,ThreePointCorr):
    #         op1LL = C3.LegLab[1:-1] ## removes $ symbols, re adds at end
    #         op2LL = C32.LegLab[1:-1] ## removes $ symbols, re adds at end
    #         if C3.name == C32.name:
    #             outname = C3.name,GetOppSelf(opp)
    #             outLL = op1LL,' Comb ',GetLegSelf(oppLL)
    #         else:
    #             for iname,iname2 in zip(C3.name.split('_'),C32.name.split('_')):
    #                 if iname in iname2:
    #                     outname.append(iname2)
    #                 elif iname2 in iname:
    #                     outname.append(iname)
    #                 else:
    #                     outname.append(iname+opp+iname2)
    #             outname = '_'.join(outname)
    #             for iLegLab,iLegLab2 in zip(op1LL.split('_'),op2LL.split('\ ')):
    #                 if iLegLab in iLegLab2:
    #                     outLL.append(iLegLab2)
    #                 elif iLegLab2 in iLegLab:
    #                     outLL.append(iLegLab)
    #                 else:
    #                     outLL.append(iLegLab+oppLL+iLegLab2)
    #             outLL = '_'.join(outLL)
    #     elif isinstance(C3,ThreePointCorr):
    #         op1LL = C3.LegLab[1:-1] ## removes $ symbols, re adds at end
    #         if isinstance(C32,int):
    #             outname  =  C3.name,opp+str(C32)
    #             outLL  =  op1LL,oppLL+str(C32)
    #         else:
    #             if C32 > 0.01:
    #                 outname = C3.name,opp+'{:.2f}'.format(C32)
    #                 outLL = op1LL,oppLL+'{:.2f}'.format(C32)
    #     else:
    #         op2LL = C32.LegLab[1:-1] ## removes $ symbols, re adds at end
    #         if isinstance(C3,int):
    #             outname  =  str(C3) + opp ,C32.name
    #             outLL  =  str(C3) + oppLL,op2LL
    #         else:
    #             if C3 > 0.01:
    #                 outname = '{:.2f}'.format(C3) + opp,C32.name
    #                 outLL = '{:.2f}'.format(C3) + oppLL,op2LL
    #     # print r'$'+r'\ '.join(outLL)+'$'
    #     outLL = list(outLL)
    #     for ic,iLL in enumerate(outLL):
    #         if '^{}' in iLL:
    #             outLL[ic] = iLL.replace('^{}','^{')+'}'
    #     self.SetCustomName('_'.join(outname),r'$'+r'\ '.join(outLL)+'$')




    ## Numerical operatoions

    # def __add__(self,C32):
    #     if type(C32) in ThreePointCorr.parentlist:
    #         return C32.__radd__(self)
    #     else:
    #         result = ThreePointCorr(thisnboot=self.nboot, thisnt=self.nt, cfglist=self.cfglist,
    #                               Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
    #         result.UpdateName('+',self,C32)
    #         if isinstance(C32,ThreePointCorr):
    #             for (ip,it,iC3) in self.items():
    #                 result[(ip,it)] = iC3 + C32[(ip,it)]
    #         else:
    #             try:
    #                 for ip,it,iC3 in self.items():
    #                     result[(ip,it)] = iC3 + C32
    #             except Exception as err:
    #                 raise EnvironmentError('Invalid value to combine with BootStrap class')
    #         return result

    # def __sub__(self,C32):
    #     if type(C32) in ThreePointCorr.parentlist:
    #         return C32.__rsub__(self)
    #     else:
    #         result = ThreePointCorr(thisnboot=self.nboot, thisnt=self.nt, cfglist=self.cfglist,
    #                               Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift427=self.thisshift)
    #         result.UpdateName('-',self,C32)
    #         if isinstance(C32,ThreePointCorr):
    #             for (ip,it,iC3) in self.items():
    #                 result[(ip,it)] = iC3 - C32[(ip,it)]
    #         else:
    #             try:
    #                 for ip,it,iC3 in self.items():
    #                     result[(ip,it)] = iC3 - C32
    #             except Exception as err:
    #                 raise EnvironmentError('Invalid value to combine with BootStrap class')
    #         return result

    # def __mul__(self,C32):
    #     if type(C32) in ThreePointCorr.parentlist:
    #         return C32.__rmul__(self)
    #     else:
    #         result = ThreePointCorr(thisnboot=self.nboot, thisnt=self.nt, cfglist=self.cfglist,
    #                               Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
    #         result.UpdateName('x',self,C32)
    #         if isinstance(C32,ThreePointCorr):
    #             for (ip,it,iC3) in self.items():
    #                 result[(ip,it)] = iC3 * C32[(ip,it)]
    #         else:
    #             try:
    #                 for ip,it,iC3 in self.items():
    #                     result[(ip,it)] = iC3 * C32
    #             except Exception as err:
    #                 raise EnvironmentError('Invalid value to combine with BootStrap class')
    #         return result

    # def __div__(self,C32):
    #     if type(C32) in ThreePointCorr.parentlist:
    #         return C32.__rdiv__(self)
    #     else:
    #         result = ThreePointCorr(thisnboot=self.nboot, thisnt=self.nt, cfglist=self.cfglist,
    #                               Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
    #         result.UpdateName('div',self,C32)
    #         if isinstance(C32,ThreePointCorr):
    #             for (ip,it,iC3) in self.items():
    #                 try:
    #                     result[(ip,it)] = iC3 / C32[(ip,it)]
    #                 except Exception as err:
    #                     if any([ibs == 0 for ibs in C32[(ip,it)].bootvals]):
    #                         raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+C32.name + ' '+C32[(ip,it)].name)

    #         else:
    #             if C32 == 0:
    #                 raise ZeroDivisionError('Dividing by zero found in bootstrap division')
    #             try:
    #                 for ip,it,iC3 in self.items():
    #                     result[(ip,it)] = iC3 / C32
    #             except Exception as err:
    #                 raise EnvironmentError('Invalid value to combine with BootStrap class')
    #         return result


    # def __pow__(self,C32):
    #     if type(C32) in ThreePointCorr.parentlist:
    #         return C32.__rpow__(self)
    #     else:
    #         result = ThreePointCorr(thisnboot=self.nboot, thisnt=self.nt, cfglist=self.cfglist,
    #                               Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
    #         result.UpdateName('pow',self,C32)
    #         if isinstance(C32,ThreePointCorr):
    #             for (ip,it,iC3) in self.items():
    #                 result[(ip,it)] = iC3 ** C32[(ip,it)]
    #         else:
    #             try:
    #                 for ip,it,iC3 in self.items():
    #                     result[(ip,it)] = iC3 ** C32
    #             except Exception as err:
    #                 raise EnvironmentError('Invalid value to combine with BootStrap class')
    #         return result


    # ## Right multiplication functions


    # def __radd__(self,C32):
    #     result = ThreePointCorr(thisnboot=self.nboot, thisnt=self.nt, cfglist=self.cfglist,
    #                           Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
    #     result.UpdateName('+',C32,self)
    #     if isinstance(C32,ThreePointCorr):
    #         for (ip,it,iC3) in self.items():
    #             result[(ip,it)] = C32[(ip,it)] + iC3
    #     else:
    #         try:
    #             for ip,it,iC3 in self.items():
    #                 result[(ip,it)] = C32 + iC3
    #         except Exception as err:
    #             raise EnvironmentError('Invalid value to combine with BootStrap class')
    #     return result

    # def __rsub__(self,C32):
    #     result = ThreePointCorr(thisnboot=self.nboot, thisnt=self.nt, cfglist=self.cfglist,
    #                           Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
    #     result.UpdateName('-',C32,self)
    #     if isinstance(C32,ThreePointCorr):
    #         for (ip,it,iC3) in self.items():
    #             result[(ip,it)] =  C32[(ip,it)] - iC3
    #     else:
    #         try:
    #             for ip,it,iC3 in self.items():
    #                 result[(ip,it)] =  C32 - iC3
    #         except Exception as err:
    #             raise EnvironmentError('Invalid value to combine with BootStrap class')
    #     return result

    # def __rmul__(self,C32):
    #     result = ThreePointCorr(thisnboot=self.nboot, thisnt=self.nt, cfglist=self.cfglist,
    #                           Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
    #     result.UpdateName('x',C32,self)
    #     if isinstance(C32,ThreePointCorr):
    #         for (ip,it,iC3) in self.items():
    #             result[(ip,it)] =  C32[(ip,it)] * iC3
    #     else:
    #         try:
    #             for ip,it,iC3 in self.items():
    #                 result[(ip,it)] = C32 * iC3
    #         except Exception as err:
    #             raise EnvironmentError('Invalid value to combine with BootStrap class')
    #     return result

    # def __rdiv__(self,C32):
    #     result = ThreePointCorr(thisnboot=self.nboot, thisnt=self.nt, cfglist=self.cfglist,
    #                           Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
    #     result.UpdateName('div',C32,self)
    #     if isinstance(C32,ThreePointCorr):
    #         for (ip,it,iC3) in self.items():
    #             try:
    #                 result[(ip,it)] = C32[(ip,it)] / iC3
    #             except Exception as err:
    #                 if any([ibs == 0 for ibs in iC3.bootvals]):
    #                     raise ZeroDivisionError('Dividing by zero found in bootstrap division for '+self.name + ' '+iC3.name)

    #     else:
    #         for ip,it,iC3 in self.items():
    #             try:
    #                 result[(ip,it)] = C32 / iC3
    #             except Exception as err:
    #                 if iC3 == 0:
    #                     raise ZeroDivisionError('Dividing by zero found in bootstrap division')
    #                 else:
    #                     raise EnvironmentError('Invalid value to combine with BootStrap class')
    #     return result


    # def __rpow__(self,C32):
    #     result = ThreePointCorr(thisnboot=self.nboot, thisnt=self.nt, cfglist=self.cfglist,
    #                           Info=self.Info,thissym=self.thissym,thiscol=self.thiscol,thisshift=self.thisshift)
    #     result.UpdateName('pow',C32,self)
    #     if isinstance(C32,ThreePointCorr):
    #         for (ip,it,iC3) in self.items():
    #             result[(ip,it)] =  C32[(ip,it)] ** iC3
    #     else:
    #         try:
    #             for ip,it,iC3 in self.items():
    #                 result[(ip,it)] = C32 ** iC3
    #         except Exception as err:
    #             raise EnvironmentError('Invalid value to combine with BootStrap class')
    #     return result


    ## unitary arithmetic operations

    def __neg__(self):
        for ival in list(self.keys()):
            self[ival] = -self[ival]
        return self


    def __pos__(self):
        return self

    def __abs__(self):
        for ival in list(self.keys()):
            self[ival] = abs(self[ival])
        return self

    ## Casting, just casts the bootstrap values, ignore return value (only needed since return has to be same type).


    def __complex__(self):
        for ival in list(self.keys()):
            complex(self[ival])
        self.Stats()
        return 1+0j

    def __int__(self):
        for ival in list(self.keys()):
            int(self[ival])
        self.Stats()
        return 1

    def __float__(self):
        for ival in list(self.keys()):
            np.float64(self[ival])
        self.Stats()
        return 1.0

####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################


class NJNQCorr(object):
    """

    NJNQ(q,tau,pp,t,tf) uses bootstrap class
    pp  = sink momentum
    t   = source sink separation time
    q   = current insersion momentum
    tau = current insersion time
    tflow = flow time

    Info is a dictionary containing information for the correlator:
    see comments below to see what is required.
    missing elements just sets field to default


    """
    def Construct_NJNQ_File(this_file):
        return Construct_File_Object(this_file,NJNQCorr)


    ## TODO: create variational flag and interface with this function.

    ## array of all classes that inherit this class. Append this when adding new classes please!
    parentlist = ['RatioCorrelators.RatioCorr']

    def __init__(self, thisnboot=nboot, cfglist={}, Info={},name='',filename='',man_load_cfgs=False):


        self.current = 0
        ## initialising instance of three point correlator as NJN
        ## DoubSing tells whether to read doublet (doub) or singlet (sing) contribution
        if 'DoubSing' in list(Info.keys()):
            if Info['DoubSing'] not in list(DSfundict.keys()):
                print('DSList (from DSDefines.py):')
                print(list(DSfundict.keys()))
                raise IOError("Doublet,Singlet type (Info['DoubSing']) not recognised")
            self.DS = Info['DoubSing']
        else: raise IOError('Please set DoubSing in Info for ThreePtCorrelators ')

        self.isDScomb = self.DS not in ['doub','sing']

        self.FO = FlowOp(thisnboot,cfglist,Info)

        thisinfo = copy(Info)
        if self.FO.Observ == '':
            thisinfo['fo_for_cfgs'] = False
        else:
            thisinfo['fo_for_cfgs'] = self.FO.Observ
        if self.isDScomb:
            thisinfo['DoubSing'] = 'doub'
            self.NJN = ThreePointCorr(thisnboot=thisnboot,cfglist=cfglist,Info=thisinfo,name=name,filename=filename)

            thisinfo = copy(Info)
            thisinfo['DoubSing'] = 'sing'
            self.NJNsing = ThreePointCorr(thisnboot=thisnboot,cfglist=cfglist,Info=thisinfo,name=name,filename=filename)
        else:
            self.NJN = ThreePointCorr(thisnboot=thisnboot,cfglist=cfglist,Info=Info,name=name,filename=filename)



        self.thischecklist = [  'nboot','nt','ism','jsm','tsink','DS','ppmom','Proj',
                                'tsrc','kud','ks','NJNfile','FOfile','tflowlist']
        self.tflowlist = self.FO.tflowlist


        ## Grab some stuff (from .NJN and .FO)
        self.NJNfile = self.NJN.PickleFile
        self.FOfile = self.FO.flowPickleFile
        self.nboot = thisnboot
        self.nt = self.NJN.nt
        self.nxyz = self.NJN.nxyz
        self.dim_label = 'RC'+str(self.nxyz)+'x'+str(self.nt)
        self.dim_ll = str(self.nxyz)+r'^{3}\times ' + str(self.nt)
        self.ism = self.NJN.ism
        self.jsm = self.NJN.jsm
        self.tsink = self.NJN.tsink
        self.ppmom = self.NJN.ppmom
        self.Proj = self.NJN.Proj
        self.tsrc = self.NJN.tsrc
        self.kud = self.NJN.kud
        self.ks = self.NJN.ks
        self.Info = Info
        self.name = name
        self.filename = filename
        self.iGuess = self.NJN.iGuess
        self.checkmomlist = self.NJN.checkmomlist

        if 'show_cfgs' in list(Info.keys()): self.show_cfgs = Info['show_cfgs']
        else: self.show_cfgs = True


        if 'tflowfit' in list(Info.keys()): self.SetTflowFit(Info['tflowfit'])
        else: self.SetTflowFit([])

        self.NJNQ_Col_Names = [ 'current_operator','momentum',
                                'source_sink_separation','flow_times']
        self.NJNQ_Stats = pa.DataFrame()

        '''
        self.NJNQ_Stats DataFrame:

        columns =   boot,   Avg,        Std

        rows =      current_operator,   momenta,    time separation,    flow_times
        '''

        self.NJNQ_Fit_Col_Names = [ 'states','current_operator','momentum',
                                    'time_fit_range','flow_times']
        self.NJNQ_Fit_Stats = pa.DataFrame()

        '''
        self.NJNQ_Fit_Stats DataFrame:

        columns =   boot,   Avg,        Std

        rows =      states fitted to,   current _operator,  momenta,
                    time fit range ,    flow_times

        Avg/Std elements are pa.Series instances with:
            index = fit parameters
            values = Avg/Std values
        '''



        self.SetCustomName(string=name,stringfn=filename)

        ## is needed?
        self.Info = Info
        if not man_load_cfgs:
            self.LoadCfgs(cfglist,name,filename)

    def LoadCfgs(self,cfglist,name='',filename=''):

        if len(cfglist) == 0:
            self.GetCfgList()
        else:
            self.ImportCfgList(cfglist)
        self.SetCustomName(string=name,stringfn=filename)


    def GetCfgList(self):
        self.CombineCfgLists()

    def CombineCfgLists(self):
        cfglist = self.NJN.C3_cfgs
        if self.isDScomb:
            cfglist = CombineCfgs(cfglist,self.NJNsing.C3_cfgs)
        cfglist = CombineCfgs(cfglist,self.FO.Op_cfgs)
        self.ImportCfgList(copy(cfglist))


    def Stats(self):
        self.NJN.Stats()
        if self.isDScomb: self.NJNsing.Stats()
        self.FO.FlowStats()
        if 'boot' not in self.NJNQ_Stats: return
        lAvg,lStd = [],[]
        for (igamma,ip,it,itflow),iNJNQ in self.items():
            iNJNQ.Stats()
            lAvg.append(iNJNQ.Avg)
            lStd.append(iNJNQ.Std)
        indicies = pa.MultiIndex.from_tuples(self.NJNQ_Stats.index,names=self.NJNQ_Col_Names)
        self.NJNQ_Stats.loc[:,'Avg'] = pa.Series(lAvg,index=indicies)
        self.NJNQ_Stats.loc[:,'Std'] = pa.Series(lStd,index=indicies)

    def Bootstrap(self,WipeData=True,DefWipe=False,show_timer=False):
        ## Do not wipe the data, we need it for combining on the configuration level!
        self.NJN.Bootstrap(WipeData=False,DefWipe=DefWipe,show_timer=show_timer)
        if self.isDScomb: self.NJNsing.Bootstrap(WipeData=False,DefWipe=DefWipe,show_timer=show_timer)
        self.FO.FlowBootstrap(WipeData=False,show_timer=show_timer)
        lNJNQ,lNJNQAvg,lNJNQStd = [],[],[]
        ilist = []
        rlist = None
        if show_timer: thistimer = Timer(   linklist=[  self.NJN.gammalist,self.NJN.qform,
                                                        list(range(self.nt)),self.FO.tflowlist],
                                            name='BS NJNQ ')
        if self.isDScomb:
            for icg,igamma in enumerate(self.NJN.gammalist):
                for icp,ip in enumerate(self.NJN.qform):
                    for it in range(self.nt):
                        strit = tstr(it+1)
                        NJN_lambda = lambda x : x[icg][icp][it]
                        doub_data = self.NJN.C3_cfgs['C3'].apply(NJN_lambda).values
                        sing_data = self.NJNsing.C3_cfgs['C3'].apply(NJN_lambda).values
                        for ictflow,itflow in enumerate(self.FO.tflowlist):
                            FO_lambda = lambda x : x[ictflow]
                            FO_data = self.FO.Op_cfgs['Op'].apply(FO_lambda)
                            thisdata = DSfundict[self.DS](doub_data,sing_data) * FO_data
                            thisboot = BootStrap(   self.nboot, name='G_{3}('+','.join([igamma,ip,strit,itflow])+')',
                                                                cfgvals=thisdata,rand_list=rlist)
                            rlist = thisboot.Get_Rand_List()
                            ilist.append((igamma,ip,strit,itflow))
                            lNJNQ.append(thisboot)
                            lNJNQAvg.append(thisboot.Avg)
                            lNJNQStd.append(thisboot.Std)
                            if show_timer: thistimer.Lap((igamma,ip,strit,itflow))
        else:
            for icg,igamma in enumerate(self.NJN.gammalist):
                for icp,ip in enumerate(self.NJN.qform):
                    for it in range(self.nt):
                        strit = tstr(it+1)
                        NJN_lambda = lambda x : np.array(x[icg][icp][it])
                        doub_data = self.NJN.C3_cfgs['C3'].apply(NJN_lambda).values
                        for ictflow,itflow in enumerate(self.FO.tflowlist):
                            FO_lambda = lambda x : np.array(x[ictflow])
                            FO_data = self.FO.Op_cfgs['Op'].apply(FO_lambda)
                            thisdata = doub_data * FO_data
                            thisboot = BootStrap(   self.nboot, name='G_{3}('+','.join([igamma,ip,strit,itflow])+')',
                                                                                    cfgvals=thisdata,rand_list=rlist)
                            rlist = thisboot.Get_Rand_List()
                            ilist.append((igamma,ip,strit,itflow))
                            lNJNQ.append(thisboot)
                            lNJNQAvg.append(thisboot.Avg)
                            lNJNQStd.append(thisboot.Std)
                            if show_timer: thistimer.Lap((igamma,ip,strit,itflow))
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.NJNQ_Col_Names)
            self.NJNQ_Stats.loc[:,'boot'] = pa.Series(lNJNQ,index=indicies)
            self.NJNQ_Stats.loc[:,'Avg'] = pa.Series(lNJNQAvg,index=indicies)
            self.NJNQ_Stats.loc[:,'Std'] = pa.Series(lNJNQStd,index=indicies)
        self.RemoveVals(WipeData= WipeData)

    def ObsToName(self,thisop):
        thisobs = self.FO.FlowObsToName(thisop)
        return thisobs.replace('<',r'\alpha_{').replace('>',r'}')

    def SetCustomName(self,string='',stringfn='',stringLL=''):
        if string == '':
            self.name = '_'.join([self.dim_label,self.NJN.kappafolder,self.NJN.stream_name,self.NJN.Proj,'_'.join(self.NJN.gammalist),self.DS,
                                  self.NJN.ppform,self.NJN.tsrc,self.NJN.tsink,self.NJN.ism,self.NJN.jsm,self.FO.Observ])
        else:
            self.name = string
        if stringfn == '':
            self.filename = '_'.join([self.dim_label,self.NJN.kappafolder,self.NJN.stream_name,self.NJN.Proj,'_'.join(self.NJN.gammalist),self.NJN.DS,
                                  self.NJN.ppform,self.NJN.qform[0],self.NJN.tsrc,self.NJN.tsink,self.NJN.ism,self.NJN.jsm,self.FO.Observ])
        else:
            self.filename = stringfn
        if stringLL == '':
            self.LegLab = '$'+'\ '.join([   self.dim_ll,self.NJN.latparams.GetPionMassLab(),self.NJN.stream_name,
                                            self.DS,self.NJN.ism,self.NJN.jsm,self.ObsToName(self.FO.Observ),self.NJN.qform[0]])+'$' ## customise this how you like
        else:
            self.LegLab = stringLL
        self.SubDir = '_'.join([self.NJN.Proj,self.DS,self.NJN.ppform])
        self.NJNQdir = outputdir + '/'+self.dim_label+self.NJN.kappafolder+'/G3/'+self.FO.Observ+'/'+self.SubDir+'/'
        mkdir_p(self.NJNQdir+'/Pickle/')
        mkdir_p(self.NJNQdir+'/Excel/')
        self.HumanFile = self.NJNQdir+self.filename+'.xml'
        self.ExcelFile = self.NJNQdir+'/Excel/'+self.filename+'.xlsx'
        self.PickleFile = self.NJNQdir+'/Pickle/'+self.filename+'.py3p'


    # def MakeOutputFile(self):
    #     self.name = '_'.join([self.tsrc,self.ism,self.jsm])
    #     self.HumanFile = '/'+self.kappafolder+'/G2/'+self.name+'.xml'
    #     self.PickleFile = '/'+self.kappafolder+'/G2/Pickle/'+self.name+'.py3p'


    def ImportCfgList(self,cfgsrclist,stream_list = 'PreDefine'):
        if not isinstance(cfgsrclist,pa.DataFrame):
            if not isinstance(stream_list,str):
                self.stream_list = stream_list
            ## cfgsrclist is [cfglist , xsrclist]
            ## cfglist must be 2-d array, first dimension is stream, second dimension is configuraiton number
            ## xsrclist is 3-d array, first dimension is stream, second dimension is cfg number, third is xsrc number
            ## xsrclist must have values 'xsrc##' with no zero pre buffering
            ilist = []
            xsrc_list = []
            for istream,(is_cfg,is_xsrc) in enumerate(zip(cfgsrclist[0],cfgsrclist[1])):
                for icf in range(is_cfg):
                    xsrc_list.append(is_xsrc)
                    ilist.append((istream,is_cfg))
            indicies = pa.MultiIndex.from_tuples(ilist,names=['stream','config_number'])
            cfg_df = pa.DataFrame(np.array(cfgsrclist[0]).flatten(),columns=['configs'],index=indicies)
            cfg_df.loc[:,'xsrc_list'] = pa.Series(xsrc_list,index=indicies)
            self.NJNQ_cfgs = cfg_df
        else:
            self.NJNQ_cfgs = cfgsrclist
        self.Update_Stream_List()
        self.ncfg_list = [icfg.size for istream,icfg in self.NJNQ_cfgs['configs'].groupby(level='stream')]
        self.xsrcAvg = np.mean(list(map(len,self.NJNQ_cfgs['xsrc_list'].values)))
        self.nMeas = self.xsrcAvg*sum(self.ncfg_list)
        self.cfglist = self.NJNQ_cfgs
        if self.isDScomb:
            self.NJNsing.ImportCfgList(copy(self.cfglist))
        self.NJN.ImportCfgList(copy(self.cfglist))
        self.FO.FlowImportCfgList(copy(self.cfglist))

    def Update_Stream_List(self):
        self.stream_list = self.NJNQ_cfgs.index.levels[0]
        self.nstream = len(self.stream_list)
        if hasattr(self,'stream_name'):
            old_sn = self.stream_name
            self.stream_name = '-'+''.join([istr.replace('-','') for istr in self.stream_list])+'-'
            if hasattr(self,'filename'):
                self.filename = self.filename.replace(old_sn,self.stream_name)
            if hasattr(self,'name'):
                self.name = self.name.replace(old_sn,self.stream_name)
            if hasattr(self,'LegLab'):
                self.LegLab = self.LegLab.replace(old_sn,self.stream_name)
        else:
            self.stream_name = '-'+''.join([istr.replace('-','') for istr in self.stream_list])+'-'


    def RemoveVals(self,WipeData = True):
        if 'C3' == WipeData:
            self.NJN.RemoveVals()
            if self.isDScomb:  self.NJNsing.RemoveVals()
        elif 'Flow' ==  WipeData:
            self.FO.FlowRemoveVals()
        elif WipeData:
            self.NJN.RemoveVals()
            if self.isDScomb:  self.NJNsing.RemoveVals()
            self.FO.FlowRemoveVals()


    def Read(self,show_timer=False):
        self.NJN.Read(show_timer=show_timer)
        if self.isDScomb: self.NJNsing.Read(show_timer=show_timer)
        self.FO.FlowRead(show_timer=show_timer)

    # def ReadAndBoot(self):
    #     self.Read()
    #     self.Bootstrap()

    def ImportFitRanges(self,fit_ranges):
        self.fit_ranges = fit_ranges
        if isinstance(self.fit_ranges,str):
            self.fit_ranges = [self.fit_ranges]


    def SetTflowFit(self,tflow):
        # if tflow not in self.FO.tflowlist and tflow != 'Not Set':
        #     print 'Warning: ', tflow, ' not in flow list, tflowfit not added'
        # else:
        self.tflowfit = tflow

    def SetTFit(self,tsink):
        # if tsink not in ['t'+str(it) for it in xrange(1,self.nt+1)] and tsink != 'Not Set':
        #     print 'Warning: ', tsink, ' not in  list, tsinkfit not added'
        # else:
        self.tsinkfit = tsink


    def ImportTwoPtFits(self,state,sourcepars,sinkpars):
        self.NJN.ImportTwoPtFits(state,sourcepars,sinkpars)

    def GetTwoPtPars(self,ip,state,pars):
        self.NJN.GetTwoPtPars(ip,state,pars)


    ## fit_ranges is list of fit ranges to fit the correlator to
    ## formatted as fitr#-# where # are the min and max fit values.
    ## states correspond to how many states to fit to (currenlty only working for states = 1 or 2)
    ## C2pars is list of size 2, containing source and sink fit parameters used in the three point correlator fits
    ## C2pars can be class instances of TwoPointCorr, Fitting, or just a list of values.
    ## Word of caution, theb ootstraps need to have been constructed with the same configurations.
    ## Best to have this run from RatioFunction class only.

    def Fit(self,C2pars='PreDef',state='PreDef',fit_ranges='PreDef',iGuess='PreDef',EstDir=False,thistflowlist= 'PreDef',WipeFit=True):
        def RunPredefFits():
            if state == 'PreDef' and fit_ranges == 'PreDef':
                for izip in self.PredefFits:
                    if len(izip) == 2:
                        self.Fit(C2pars=C2pars,state=izip[0],fit_ranges=[izip[1]],EstDir=EstDir,WipeFit=WipeFit)
                    else:
                        self.Fit(C2pars=C2pars,state=izip[0],fit_ranges=[izip[1]],iGuess=izip[2],EstDir=EstDir,WipeFit=WipeFit)
            elif state == 'PreDef':
                for izip in self.PredefFits:
                    if izip[1] == fit_ranges:
                        if len(izip) == 2:
                            self.Fit(C2pars=C2pars,state=izip[0],fit_ranges=[izip[1]],EstDir=EstDir,WipeFit=WipeFit)
                        else:
                            self.Fit(C2pars=C2pars,state=izip[0],fit_ranges=[izip[1]],iGuess=izip[2],EstDir=EstDir,WipeFit=WipeFit)
            elif fit_ranges == 'PreDef':
                for izip in self.PredefFits:
                    if izip[0] == state:
                        if len(izip) == 2:
                            self.Fit(C2pars=C2pars,state=izip[0],fit_ranges=[izip[1]],EstDir=EstDir,WipeFit=WipeFit)
                        else:
                            self.Fit(C2pars=C2pars,state=izip[0],fit_ranges=[izip[1]],iGuess=izip[2],EstDir=EstDir,WipeFit=WipeFit)
        if 'PreDef' in [state,fit_ranges]:
            RunPredefFits()
            return
        if iGuess != 'PreDef': self.iGuess = iGuess
        if isinstance(state,str):
            if 'state' in state:
                state = state.replace('state','')
            state = int(state)

        if C2pars != 'PreDef': self.ImportTwoPtFits(state,C2pars[0],C2pars[1])
        self.ImportFitRanges(fit_ranges)
        if state == 1:
            # thisfitfun = C3OneStateFitFunNoExp
            thisfitfun = C3OneStateFitFunNoExp
            npar = 1
        elif state == 2:
            # thisfitfun = C3TwoStateFitFunNoExp
            thisfitfun = C3MomTwoStateFitFunNoExp
            npar = 4
        elif state > 2:
            raise IOError('fitting state is only implemented up to 2 state so far. Cannot do '+str(state))

        if thistflowlist == 'PreDef':
            if self.tflowfit == 'Not Set': raise IOError('tflowfit not set, please set by using .SetTflowFit()')
        else:
            if thistflowlist == 'All':
                self.tflowfit = self.FO.tflowlist
            else:
                self.tflowfit = thistflowlist

        ststr = 'state'+str(state)
        legstate = str(state)+'SF'
        lFit,lFitA,lFitS = [],[],[]
        ilist = []
        for (igamma,itflow),gammadata in self.NJNQ_Stats['boot'].groupby(level=('current_operator','flow_times')):
            for (ip,pdata),ptwopt in zip(gammadata.groupby(level='momentum'),self.NJN.TwoPtPars):
                for ifit in self.fit_ranges:
                    if (ststr,igamma,ip,ifit) in self.NJNQ_Fit_Stats.index:
                        if WipeFit:
                            self.NJNQ_Fit_Stats.drop([(ststr,igamma,ip,ifit)],inplace=True)
                        else:
                            continue
                    # print 'Fitting state',str(state) ,igamma,itflow,ifit
                    ilist.append((ststr,igamma,ip,ifit))
                    ifitmin,ifitmax = list(map(int,unxmlfitr(ifit)))
                    ## is this the right ordering of dictionary dimensions?
                    tdatarange = [[]]
                    ydata = []
                    for it,tdata in pdata.items():
                        if itflow not in list(tdata.keys()): continue
                        tint = untstr(it)
                        if ifitmin <= tint <= ifitmax:
                            tcurrboot = BootStrap(self.nboot,bootvals=[np.float64(tint)]*self.nboot)
                            tdatarange[0].append(tcurrboot)
                            ydata.append(tdata[itflow])
                    thistlen = len(tdatarange[0])
                    ## this + 3 is just calabrated from looking at the conserved current two- and three- point correlators. Not sure why exactly.
                    tsinkboot = BootStrap(self.nboot,bootvals=[np.float64(self.tsink.replace('tsink',''))+3]*self.nboot)
                    tdatarange.append([tsinkboot]*thistlen)
                    for iparam in ptwopt:
                        tdatarange.append([iparam]*thistlen)

                    # print 'DEBUG'
                    # for it,iy in zip(np.rollaxis(np.array(tdatarange),1),ydata):
                    #     print it[0].Avg,it[1].Avg,it[2].Avg,it[3].Avg,it[4].Avg,it[5].Avg,iy.Avg

                    if EstDir:
                        thisff = ff.Fitting(Funs=[thisfitfun,npar,'Estimate'],data=[tdatarange,ydata],name=' '.join([legstate,ifit,itflow]),iGuess=self.iGuess)
                    else:
                        thisff = ff.Fitting(Funs=[thisfitfun,npar],data=[tdatarange,ydata],name=' '.join([legstate,ifit,itflow]),iGuess=self.iGuess)
                    thistest = thisff.FitBoots()
                    if not np.isnan(thistest):
                        lFit.append(thisff)
                        lFitA.append(thisff.fit_data['ParamsAvg'])
                        lFitS.append(thisff.fit_data['ParamsStd'])
        indicies = pa.MultiIndex.from_tuples(ilist,names=self.NJNQ_Fit_Col_Names)
        if 'boot' in self.NJNQ_Fit_Stats.columns:
            this_df = pa.DataFrame()
            this_df.loc[:,'boot'] = pa.Series(lFit,index=indicies)
            this_df.loc[:,'Avg'] = pa.Series(lFitA,index=indicies)
            this_df.loc[:,'Std'] = pa.Series(lFitS,index=indicies)
            self.NJNQ_Fit_Stats = self.NJNQ_Fit_Stats.append(this_df)
        else:
            self.NJNQ_Fit_Stats.loc[:,'boot'] = pa.Series(lFit,index=indicies)
            self.NJNQ_Fit_Stats.loc[:,'Avg'] = pa.Series(lFitA,index=indicies)
            self.NJNQ_Fit_Stats.loc[:,'Std'] = pa.Series(lFitS,index=indicies)
        self.Write()



    def Write(self):
        def FixDictArray(thislist,ikey):
            outDict = OrderedDict()
            for ic,ilist in enumerate(thislist):
                outDict[ikey+str(ic+1)] = ilist
            return outDict

        self.NJN.Write()
        if self.isDScomb: self.NJNsing.Write()
        self.FO.FlowWrite()


        if not hasattr(self,'show_cfgs'):
            self.show_cfgs = True

        outDict = ODNested()


        if self.isDScomb:
            outDict['doub_ThreePC_File'] = self.NJN.HumanFile
            outDict['sing_ThreePC_File'] = self.NJNsing.HumanFile
        else:
            outDict['ThreePC_File'] = self.NJN.HumanFile
        outDict['FlowOp_File'] = self.FO.flowHumanFile

        excel_params = pa.Series(deepcopy(outDict))

        for col_key,col_data in self.NJNQ_Stats.items():
            if 'Avg' in col_key or 'Std' in col_key: continue
            for (igamma,ip,it,itflow),idata in col_data.items():
                outDict[col_key][igamma][ip][it][itflow] = AvgStdToFormat(idata.Avg,idata.Std)


        if 'boot' in self.NJNQ_Fit_Stats:
            for (istate,igamma,ip,itr,itflow),fitdict in self.NJNQ_Fit_Stats['boot'].items():
                if not isinstance(fitdict,ff.Fitting): continue
                outDict['Fits'][istate][igamma][ip][itr][itflow] = fitdict.GetOutputDict()



        if any(np.array(self.ncfg_list) > 0) and self.show_cfgs:
            for ((istream,iccfg),icfg),ixsrc_list in zip(iter(self.NJNQ_cfgs['configs'].items()),self.NJNQ_cfgs['xsrc_list'].values):
                outDict['cfglist'][istream][icfg] = 'n_xsrc_'+str(len(ixsrc_list))



        ## Write human readable file with data
        WriteXml(self.HumanFile,{'Results':outDict})

        if any(np.array(self.ncfg_list) > 0):
            for ((istream,iccfg),icfg),ixsrc_list in zip(iter(self.NJNQ_cfgs['configs'].items()),self.NJNQ_cfgs['xsrc_list'].values):
                outDict['cfglist'][istream][icfg] = 'n_xsrc_'+str(len(ixsrc_list))

        WriteExcel(self.ExcelFile,{'NJNQ_results':deepcopy(self.NJNQ_Stats)},cfglist=self.NJNQ_cfgs,params=excel_params)

        ## pickles rest of data for reading
        self.RemoveFuns()
        WritePickle(self.PickleFile,self.__dict__)
        self.GetFuns()
    ##################################

    def RemoveFuns(self):
        self.FO.RemoveFuns()

    def GetFuns(self):
        self.FO.GetFuns()

    def ReadAndWrite(self,WipeData=True,show_timer=False):
        self.Read(show_timer=show_timer)
        self.Bootstrap(WipeData=WipeData,show_timer=show_timer)
        self.Stats()
        # self.Fit()
        self.Write()

    def LoadPickle(self,DefWipe=False,WipeData = True,CheckCfgs=False,CheckMom=True,show_timer=False):
        if os.path.isfile(self.PickleFile) and not DefWipe:
            # print 'Loading Pickle for ' , self.name
            loadeddict = ReadPickleWrap(self.PickleFile)
            checklist = self.thischecklist
            if CheckCfgs: checklist = checklist + ['cfglist']
            if CheckMom: checklist = checklist + ['checkmomlist']
            if CheckClass(self.__dict__,loadeddict,checklist):
                self.__dict__.update(loadeddict)
                if Wipe_All_Fits:
                    self.NJNQ_Fit_Stats = pa.DataFrame()
                # self.read_cfgs_from_file = True
            else:
                if os.path.isfile(self.PickleFile+'.bak'):
                    print('using backupfile for '+self.name)
                    self.PickleFile = self.PickleFile+'.bak'
                    self.HumanFile = self.HumanFile+'.bak'
                    self.ExcelFile = self.ExcelFile.replace('.xlsx','.bak.xlsx')
                    self.LoadPickle(DefWipe=DefWipe,WipeData=WipeData,CheckCfgs=CheckCfgs,CheckMom=CheckMom)
                    return
                print('Warning, file ' , self.HumanFile , ' has different parameters than this instance, reading and writing over file now')
                print()
                self.ReadAndWrite(WipeData=WipeData,show_timer=show_timer)
            # self.Write()
        else:
            self.ReadAndWrite(WipeData=WipeData,show_timer=show_timer)



    ## Comparisons

    def __str__(self):
        return '_'.join(self.name)


    def __iter__(self):
        return self

    def iteritems(self):
        return iter(self.NJNQ_Stats['boot'].items())

    def items(self):
        return list(self.NJNQ_Stats['boot'].items())

    def itemsAvgStd(self):
        outlist = []
        for ivals in self.NJNQ_Stats.itertuples():
            outlist.append((ivals[0],ivals[2],ivals[3]))
        return outlist

    def values(self):
        return self.NJNQ_Stats['boot'].values

    def itervalues(self):
        return iter(self.NJNQ_Stats['boot'].values)

    def valuesAvgStd(self):
        outlist = []
        for ivals in self.NJNQ_Stats.itertuples():
            outlist.append((ivals[1],ivals[2]))
        return outlist

    def keys(self):
        return self.NJNQ_Stats.index


    def __setitem__(self, key, value ):
        self.NJNQ_Stats.loc[key,'boot'] = value

    def __getitem__(self, key):
        return self.NJNQ_Stats.loc[key,'boot']

    def __reversed__(self):
        self.NJNQ_Stats = self.NJNQ_Stats.iiloc[::-1]
        return self



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
        for ival in self.keys():
            self[ival] = -self[ival]
        return self


    def __pos__(self):
        return self

    def __abs__(self):
        for ival in self.keys():
            self[ival] = abs(self[ival])
        return self

    ## Casting, just casts the bootstrap values, ignore return value (only needed since return has to be same type).


    def __complex__(self):
        for ival in self.keys():
            complex(self[ival])
        self.Stats()
        return 1+0j

    def __int__(self):
        for ival in self.keys():
            int(self[ival])
        self.Stats()
        return 1

    def __float__(self):
        for ival in self.keys():
            np.float64(self[ival])
        self.Stats()
        return 1.0

####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################
####################################################################################################################


class NJNQFullCorr(object):
    """

    NJNQ(q,tau,pp,t,tf,ts) uses bootstrap class
    pp  = sink momentum
    t   = source sink separation time
    q   = current insersion momentum
    tau = current insersion time
    tflow = flow time
    ts = operator insersion time (  refered to "sum time", to keep in
                                    line with same notation for other quantities)

    Info is a dictionary containing information for the correlator:
    see comments below to see what is required.
    missing elements just sets field to default


    """

    def Construct_NJNQFull_File(this_file):
        return Construct_File_Object(this_file,NJNQFullCorr)

    ## TODO: create variational flag and interface with this function.

    ## array of all classes that inherit this class. Append this when adding new classes please!
    parentlist = ['RatioCorrelators.RatioCorr']

    def __init__(self, thisnboot=nboot, cfglist={}, Info={},name='',filename='',man_load_cfgs=False):


        self.current = 0
        ## initialising instance of three point correlator as NJN
        ## DoubSing tells whether to read doublet (doub) or singlet (sing) contribution
        if 'DoubSing' in list(Info.keys()):
            if Info['DoubSing'] not in list(DSfundict.keys()):
                print('DSList (from DSDefines.py):')
                print(list(DSfundict.keys()))
                raise IOError("Doublet,Singlet type (Info['DoubSing']) not recognised")
            self.DS = Info['DoubSing']
        else: raise IOError('Please set DoubSing in Info for ThreePtCorrelators ')

        self.isDScomb = self.DS not in ['doub','sing']

        self.FO = FlowOpFullSquared(thisnboot,cfglist,Info)

        thisinfo = copy(Info)
        if self.FO.Observ == '':
            thisinfo['fo_for_cfgs'] = False
        else:
            thisinfo['fo_for_cfgs'] = self.FO.Observ
        if self.isDScomb:
            thisinfo['DoubSing'] = 'doub'
            self.NJN = ThreePointCorr(thisnboot=thisnboot,cfglist=cfglist,Info=thisinfo,name=name)

            singinfo = copy(Info)
            singinfo['DoubSing'] = 'sing'
            self.NJNsing = ThreePointCorr(thisnboot=thisnboot,cfglist=cfglist,Info=singinfo,name=name)
        else:
            self.NJN = ThreePointCorr(thisnboot=thisnboot,cfglist=cfglist,Info=Info,name=name)



        self.thischecklist = [  'nboot','nt','ism','jsm','tsink','DS','ppmom','Proj',
                                'tsrc','kud','ks','NJNfile','FOfile','tflowlist','boot_tsum_list']
        self.tflowlist = self.FO.tflowlist


        ## Grab some stuff (from .NJN and .FO)
        self.NJNfile = self.NJN.PickleFile
        self.FOfile = self.FO.flowPickleFile
        self.nboot = thisnboot
        self.nt = self.NJN.nt
        self.nxyz = self.NJN.nxyz
        self.dim_label = 'RC'+str(self.nxyz)+'x'+str(self.nt)
        self.dim_ll = str(self.nxyz)+r'^{3}\times ' + str(self.nt)
        self.ism = self.NJN.ism
        self.jsm = self.NJN.jsm
        self.tsink = self.NJN.tsink
        self.ppmom = self.NJN.ppmom
        self.Proj = self.NJN.Proj
        self.tsrc = self.NJN.tsrc
        self.kud = self.NJN.kud
        self.ks = self.NJN.ks
        self.Info = Info
        self.name = name
        self.filename = filename
        self.iGuess = self.NJN.iGuess
        self.checkmomlist = self.NJN.checkmomlist

        if 'tsum_fit_list' in list(Info.keys()):
            self.tsum_fit_list = list(map(tsumstr,Info['tsum_fit_list']))
        else:
            self.tsum_fit_list = list(map(tsumstr,list(range(self.nt))))


        self.boot_tsum_list = list(range(self.nt))

        # if 'boot_tsum_list' in Info.keys():
        #     self.boot_tsum_list = map(untstr,Info['boot_tsum_list'])
        # else:
        #     self.boot_tsum_list = range(self.nt)
        #
        # new_blist = []
        # for itboot in self.boot_tsum_list:
        #     if 0 <= itboot < self.nt:
        #         new_blist.append(itboot)
        # self.boot_tsum_list = new_blist


        if 'show_cfgs' in list(Info.keys()): self.show_cfgs = Info['show_cfgs']
        else: self.show_cfgs = True


        if 'tflowfit' in list(Info.keys()): self.SetTflowFit(Info['tflowfit'])
        else: self.SetTflowFit([])

        '''
        I'm thinking types could be:
        from_src_sink_sum       is summing from source and sink, outwards
        from_src_sum            is summing from source
        from_sink_sum           is summing from sink
        from_current_sum        is summing from current point

        from_src_sym_sum                is summing symmetrically from source
        from_sink_sym_sum               is summing symmetrically from sink
        from_current_sym_sum            is summing symmetrically from current poin

        from_src_sink           is values from source and sink, outwards
        from_src                is values from source
        from_sink               is values from sink
        from_current            is values from current point

        from_src_sym                is values symmetrically from source
        from_sink_sym               is values symmetrically from sink
        from_current_sym            is values symmetrically from current point

        etc...
        '''
        if 'op_trange_sum_type_3pt' in list(Info.keys()):
            self.tr_sum_type = Info['op_trange_sum_type_3pt'].replace(')None','')
        else:
            self.tr_sum_type = 'from_src_sink'

        self.tr_sum_name = self.tr_sum_type.replace('src','sr')
        self.tr_sum_name = self.tr_sum_name.replace('sink','si')
        self.tr_sum_name = self.tr_sum_name.replace('sum','su')
        self.tr_sum_name = self.tr_sum_name.replace('current','cu')
        self.tr_sum_name = self.tr_sum_name.replace('sym','sy')
        self.tr_sum_name = self.tr_sum_name.replace('from','')

        self.NJNQFull_Col_Names = [ 'current_operator','momentum',
                                'source_sink_separation','flow_times','flow_op_trange']
        self.NJNQFull_Stats = pa.DataFrame()

        '''
        self.NJNQFull_Stats DataFrame:

        columns =   boot,   Avg,        Std

        rows =      current_operator,   momenta,    time separation,
                    flow_times,         flow op trange
        '''

        self.NJNQInt_Col_Names = [ 'current_operator','momentum',
                                'source_sink_separation','flow_times']
        self.NJNQInt_Stats = pa.DataFrame()

        '''
        self.NJNQFull_Stats DataFrame:

        columns =   boot,   Avg,        Std

        rows =      current_operator,   momenta,    time separation,    flow_times
        '''


        self.NJNQFull_Fit_Col_Names = [ 'states','current_operator','momentum',
                                    'time_fit_range','flow_times']
        self.NJNQFull_Fit_Stats = pa.DataFrame()

        '''
        self.NJNQFull_Fit_Stats DataFrame:

        columns =   boot,   Avg,        Std

        rows =      states fitted to,   current _operator,  momenta,
                    time fit range ,    flow_times

        Avg/Std elements are pa.Series instances with:
            index = fit parameters
            values = Avg/Std values
        '''



        self.SetCustomName(string=name,stringfn=filename)
        ## is needed?
        self.Info = Info


        if not man_load_cfgs:
            self.LoadCfgs(cfglist,name,filename)

    def LoadCfgs(self,cfglist,name='',filename=''):

        if len(cfglist) == 0:
            self.GetCfgList()
        else:
            self.ImportCfgList(cfglist)
        self.SetCustomName(string=name,stringfn=filename)




    def GetCfgList(self):
        self.CombineCfgLists()

    def CombineCfgLists(self):
        cfglist = self.NJN.C3_cfgs
        if self.isDScomb:
            cfglist = CombineCfgs(cfglist,self.NJNsing.C3_cfgs)
        cfglist = CombineCfgs(cfglist,self.FO.Op_cfgs)
        self.ImportCfgList(copy(cfglist))


    def Stats(self):
        self.NJN.Stats()
        if self.isDScomb: self.NJNsing.Stats()
        self.FO.FlowStats()
        if 'boot' not in self.NJNQFull_Stats: return
        lAvg,lStd = [],[]
        for (igamma,ip,it,itflow,its),iNJNQFull in self.items():
            iNJNQFull.Stats()
            lAvg.append(iNJNQFull.Avg)
            lStd.append(iNJNQFull.Std)
        indicies = pa.MultiIndex.from_tuples(self.NJNQFull_Stats.index,names=self.NJNQFull_Col_Names)
        self.NJNQFull_Stats.loc[:,'Avg'] = pa.Series(lAvg,index=indicies)
        self.NJNQFull_Stats.loc[:,'Std'] = pa.Series(lStd,index=indicies)

    def Create_NJNQ_fun(self,it_sum):
        def NJNQ_fun(this_df):
            # print 'DEBUG'
            # print np.array(this_df['iC3']).shape
            # print np.array(this_df['iFO'])[:,:].shape, it_sum
            try:
                return np.mean(this_df['iC3']*np.array(this_df['iFO'])[:,it_sum])
            except Exception as err:
                return float('NaN')
        return NJNQ_fun

    def Create_Tsrc_Fun(self,ictflow,it_curr,it_sink):
        if isinstance(it_sink,str):
            it_sink = untstr(it_sink.replace('tsink','t'))
        if isinstance(it_curr,str):
            it_curr = untstr(it_curr)
        def Tsrc_Fun(this_series):
            ## tsrc_list = [ ixsrc ]
            out_list,tflow_rolled = [],[]
            flow_op = this_series['Op'][ictflow]
            try:
                for itsrc in this_series['tsrc_list']:
                    if 'from_current' in self.tr_sum_type:
                        n_roll = (itsrc+it_curr/2)%self.nt
                    elif 'src' in self.tr_sum_type:
                        n_roll = itsrc
                    elif 'from_sink' in self.tr_sum_type:
                        n_roll = (itsrc+it_sink)%self.nt
                    else:
                        raise EnvironmentError(self.tr_sum_type + ' sum time not supported.')
                    ## op_series = [ tflow, t_src ]
                    tflow_rolled = np.roll(flow_op,-n_roll)
                    if '_src_sink' in self.tr_sum_type:
                        tflow_rolled_tsink = np.roll(np.roll(tflow_rolled,-it_sink)[::-1],1)
                        if '_sym' in self.tr_sum_type:
                            tflow_rolled = tflow_rolled+np.roll(tflow_rolled[::-1],1)
                            tflow_rolled_tsink = tflow_rolled_tsink+tflow_rolled_tsink[::-1]
                        elif '_rev' in self.tr_sum_type:
                            tflow_rolled = tflow_rolled[::-1]
                            tflow_rolled_tsink = tflow_rolled_tsink[::-1]
                        tflow_rolled = (tflow_rolled + tflow_rolled_tsink)/2.
                    elif '_rev' in self.tr_sum_type:
                        tflow_rolled = np.roll(tflow_rolled[::-1],1)
                    elif '_sym' in self.tr_sum_type:
                        tflow_rolled = tflow_rolled+np.roll(tflow_rolled[::-1],1)
                    if '_sum' in self.tr_sum_type:
                        if '_sym' in self.tr_sum_type and '_src_sink' not in self.tr_sum_type:
                            tflow_rolled[0] = tflow_rolled[0]/2.
                            if float(self.nt//2) == self.nt/2.:
                                tflow_rolled[int(self.nt/2-1)] = tflow_rolled[int(self.nt/2-1)]/2
                        out_list.append(np.cumsum(tflow_rolled))
                    else:
                        out_list.append(tflow_rolled*self.nt)
            except Exception as err:
                print('File Error:')
                print(self.NJN.CfgFile)
                print(self.FO.flowCfgFile)
                print('error contructing function for ',ictflow,it_curr,it_sink)
                print('this locations list is:' , ', '.join(map(str,this_series.name)))
                print('tsrc_list for this location is', this_series['tsrc_list'])
                print()
                print('Full series data:')
                print(this_series)
                print()
                print(err)
                return pa.Series('Failed',index=['iFO'])
            return pa.Series([out_list],index=['iFO'])
        return Tsrc_Fun


    def Bootstrap(self,WipeData=True,DefWipe=False,show_timer=False):
        ## Do not wipe the data, we need it for combining on the configuration level!
        # self.NJN.Bootstrap(WipeData=False,DefWipe=DefWipe)
        # if self.isDScomb: self.NJNsing.Bootstrap(WipeData=False,DefWipe=DefWipe)
        # self.FO.FlowBootstrap(WipeData=False,show_timer=False)
        if 'C3' not in self.NJN.C3_cfgs or 'FO' not in self.FO.Op_cfgs:
            self.Read(show_timer=show_timer)
        lNJNQFull,lNJNQFullAvg,lNJNQFullStd = [],[],[]
        lNJNQInt,lNJNQIntAvg,lNJNQIntStd = [],[],[]
        ilist,iIntlist = [],[]
        rlist = None
        self.FO.Op_cfgs['tsrc_list'] = self.NJN.C3_cfgs['tsrc_list']
        if show_timer: thistimer = Timer(   linklist=[  self.NJN.gammalist,self.NJN.qform,
                                                        list(range(self.nt)),self.FO.tflowlist,
                                                        self.boot_tsum_list],
                                            name='BS NJNQFull ')
        if self.isDScomb:
            for icg,igamma in enumerate(self.NJN.gammalist):
                for icp,ip in enumerate(self.NJN.qform):
                    for it in range(self.nt):
                        strit = tstr(it+1)
                        def C3_pick_fun(x):
                            return np.array(x[icg][icp][it])
                        doub_data = self.NJN.C3_cfgs['C3'].apply(C3_pick_fun)
                        sing_data = self.NJNsing.C3_cfgs['C3'].apply(C3_pick_fun)
                        thisdata = DSfundict[self.DS](doub_data,sing_data).to_frame('iC3')
                        for ictflow,itflow in enumerate(self.FO.tflowlist):
                            FO_lambda = lambda x : np.sum(x[ictflow])
                            FO_data = self.FO.Op_cfgs['Op'].apply(FO_lambda)
                            passdata = thisdata['iC3']* FO_data
                            passdata = passdata.apply(lambda x : np.mean(x))
                            thisboot = BootStrap(   self.nboot, name='G_{3}('+','.join([igamma,ip,strit,itflow])+')',
                                                                        cfgvals=passdata,rand_list=rlist)
                            rlist = thisboot.Get_Rand_List()
                            iIntlist.append((igamma,ip,strit,itflow))
                            lNJNQInt.append(thisboot)
                            lNJNQIntAvg.append(thisboot.Avg)
                            lNJNQIntStd.append(thisboot.Std)
                            thisdata['iFO'] = self.FO.Op_cfgs.apply(self.Create_Tsrc_Fun(ictflow,it,self.tsink),axis=1)['iFO']
                            if isinstance(thisdata['iFO'],str) and thisdata['iFO'] == 'Failed': continue
                            for it_sum in self.boot_tsum_list:
                                tsum_str = tstr(it_sum).replace('t','ts')
                                combdata = thisdata.apply(self.Create_NJNQ_fun(it_sum),axis=1)
                                if isinstance(combdata,(float)) and np.isnan(combdata): continue
                                ilist.append((igamma,ip,strit,itflow,tsum_str))
                                thisboot = BootStrap(   self.nboot, name='G_{3}('+','.join([igamma,ip,strit,itflow,tsum_str])+')',
                                                                                        cfgvals=combdata,rand_list=rlist)
                                lNJNQFull.append(thisboot)
                                lNJNQFullAvg.append(thisboot.Avg)
                                lNJNQFullStd.append(thisboot.Std)
                                if show_timer: thistimer.Lap((igamma,ip,strit,itflow,tsum_str))
        else:
            for icg,igamma in enumerate(self.NJN.gammalist):
                for icp,ip in enumerate(self.NJN.qform):
                    for it in range(self.nt):
                        strit = tstr(it+1)
                        def C3_pick_fun(x):
                            return np.array(x[icg][icp][it])
                        thisdata = self.NJN.C3_cfgs['C3'].apply(C3_pick_fun).to_frame('iC3')
                        for ictflow,itflow in enumerate(self.FO.tflowlist):
                            FO_lambda = lambda x : np.sum(x[ictflow])
                            FO_data = self.FO.Op_cfgs['Op'].apply(FO_lambda)
                            passdata = thisdata['iC3']* FO_data
                            passdata = passdata.apply(lambda x : np.mean(x))
                            thisboot = BootStrap(   self.nboot,name='G_{3}('+','.join([igamma,ip,strit,itflow])+')',
                                                    cfgvals=passdata,rand_list=rlist)
                            rlist = thisboot.Get_Rand_List()
                            iIntlist.append((igamma,ip,strit,itflow))
                            lNJNQInt.append(thisboot)
                            lNJNQIntAvg.append(thisboot.Avg)
                            lNJNQIntStd.append(thisboot.Std)
                            thisdata['iFO'] = self.FO.Op_cfgs.apply(self.Create_Tsrc_Fun(ictflow,it,self.tsink),axis=1)['iFO']
                            if isinstance(thisdata['iFO'],str) and thisdata['iFO'] == 'Failed': continue
                            for it_sum in self.boot_tsum_list:
                                tsum_str = tstr(it_sum).replace('t','ts')
                                combdata = thisdata.apply(self.Create_NJNQ_fun(it_sum),axis=1)
                                if isinstance(combdata,(float)) and np.isnan(combdata): continue
                                ilist.append((igamma,ip,strit,itflow,tsum_str))
                                thisboot = BootStrap(   self.nboot, name='G_{3}('+','.join([igamma,ip,strit,itflow,tsum_str])+')',
                                                                                        cfgvals=combdata,rand_list=rlist)
                                lNJNQFull.append(thisboot)
                                lNJNQFullAvg.append(thisboot.Avg)
                                lNJNQFullStd.append(thisboot.Std)
                                if show_timer: thistimer.Lap((igamma,ip,strit,itflow,tsum_str))
        if len(ilist) > 0:
            indicies = pa.MultiIndex.from_tuples(ilist,names=self.NJNQFull_Col_Names)
            self.NJNQFull_Stats.loc[:,'boot'] = pa.Series(lNJNQFull,index=indicies)
            self.NJNQFull_Stats.loc[:,'Avg'] = pa.Series(lNJNQFullAvg,index=indicies)
            self.NJNQFull_Stats.loc[:,'Std'] = pa.Series(lNJNQFullStd,index=indicies)
        if len(iIntlist) > 0:
            indicies = pa.MultiIndex.from_tuples(iIntlist,names=self.NJNQInt_Col_Names)
            self.NJNQInt_Stats.loc[:,'boot'] = pa.Series(lNJNQInt,index=indicies)
            self.NJNQInt_Stats.loc[:,'Avg'] = pa.Series(lNJNQIntAvg,index=indicies)
            self.NJNQInt_Stats.loc[:,'Std'] = pa.Series(lNJNQIntStd,index=indicies)
        self.RemoveVals(WipeData= WipeData)

    def ObsToName(self,thisop):
        thisobs = self.FO.FlowObsToName(thisop)
        return thisobs.replace('<',r'\alpha_{').replace('>',r'}')

    def SetCustomName(self,string='',stringfn='',stringLL=''):
        if string == '':
            self.name = '_'.join([  self.dim_label,self.NJN.kappafolder,self.NJN.stream_name,
                                    self.NJN.Proj,'_'.join(self.NJN.gammalist),'_'.join(self.NJN.qform),self.DS,
                                    self.NJN.ppform,self.NJN.tsrc,self.NJN.tsink,
                                    self.NJN.ism,self.NJN.jsm,self.FO.Observ])+self.tr_sum_name
        else:
            self.name = string
        if stringfn == '':
            self.filename = '_'.join([self.dim_label,self.NJN.kappafolder,self.NJN.stream_name,self.NJN.Proj,'_'.join(self.NJN.gammalist),self.NJN.DS,
                                  self.NJN.ppform,self.NJN.qform[0],self.NJN.tsrc,self.NJN.tsink,
                                  self.NJN.ism,self.NJN.jsm,self.FO.Observ])+self.tr_sum_name
        else:
            self.filename = stringfn
        if stringLL == '':
            self.LegLab = '$'+'\ '.join([   self.dim_ll,self.NJN.latparams.GetPionMassLab(),self.NJN.stream_name,
                                            self.DS,self.NJN.ism,self.NJN.jsm,self.ObsToName(self.FO.Observ),self.NJN.qform[0]]) \
                                            +self.tr_sum_name.replace('_','\ ')+'$' ## customise this how you like
        else:
            self.LegLab = stringLL
        self.SubDir = '_'.join([self.NJN.Proj,self.DS,self.NJN.ppform])
        self.NJNQFulldir = outputdir + '/'+self.dim_label+self.NJN.kappafolder+'/G3/'+self.FO.Observ+'/'+self.SubDir+'/'
        mkdir_p(self.NJNQFulldir+'/Pickle/')
        mkdir_p(self.NJNQFulldir+'/Excel/')
        self.HumanFile = self.NJNQFulldir+self.filename+'.xml'
        self.ExcelFile = self.NJNQFulldir+'/Excel/'+self.filename+'.xlsx'
        self.PickleFile = self.NJNQFulldir+'/Pickle/'+self.filename+'.py3p'


    # def MakeOutputFile(self):
    #     self.name = '_'.join([self.tsrc,self.ism,self.jsm])
    #     self.HumanFile = '/'+self.kappafolder+'/G2/'+self.name+'.xml'
    #     self.PickleFile = '/'+self.kappafolder+'/G2/Pickle/'+self.name+'.py3p'


    def ImportCfgList(self,cfgsrclist,stream_list = 'PreDefine'):
        if not isinstance(cfgsrclist,pa.DataFrame):
            if not isinstance(stream_list,str):
                self.stream_list = stream_list
            ## cfgsrclist is [cfglist , xsrclist]
            ## cfglist must be 2-d array, first dimension is stream, second dimension is configuraiton number
            ## xsrclist is 3-d array, first dimension is stream, second dimension is cfg number, third is xsrc number
            ## xsrclist must have values 'xsrc##' with no zero pre buffering
            ilist = []
            xsrc_list = []
            for istream,(is_cfg,is_xsrc) in enumerate(zip(cfgsrclist[0],cfgsrclist[1])):
                for icf in range(is_cfg):
                    xsrc_list.append(is_xsrc)
                    ilist.append((istream,is_cfg))
            indicies = pa.MultiIndex.from_tuples(ilist,names=['stream','config_number'])
            cfg_df = pa.DataFrame(np.array(cfgsrclist[0]).flatten(),columns=['configs'],index=indicies)
            cfg_df.loc[:,'xsrc_list'] = pa.Series(xsrc_list,index=indicies)
            self.NJNQFull_cfgs = cfg_df
        else:
            self.NJNQFull_cfgs = cfgsrclist
        self.Update_Stream_List()
        self.ncfg_list = [icfg.size for istream,icfg in self.NJNQFull_cfgs['configs'].groupby(level='stream')]
        self.xsrcAvg = np.mean(list(map(len,self.NJNQFull_cfgs['xsrc_list'].values)))
        self.nMeas = self.xsrcAvg*sum(self.ncfg_list)
        self.cfglist = self.NJNQFull_cfgs
        if self.isDScomb:
            self.NJNsing.ImportCfgList(copy(self.cfglist))
        self.NJN.ImportCfgList(copy(self.cfglist))
        self.FO.FlowImportCfgList(copy(self.cfglist))


    def Update_Stream_List(self):
        self.stream_list = self.NJNQFull_cfgs.index.levels[0]
        self.nstream = len(self.stream_list)
        if hasattr(self,'stream_name'):
            old_sn = self.stream_name
            self.stream_name = '-'+''.join([istr.replace('-','') for istr in self.stream_list])+'-'
            if hasattr(self,'filename'):
                self.filename = self.filename.replace(old_sn,self.stream_name)
            if hasattr(self,'name'):
                self.name = self.name.replace(old_sn,self.stream_name)
            if hasattr(self,'LegLab'):
                self.LegLab = self.LegLab.replace(old_sn,self.stream_name)
        else:
            self.stream_name = '-'+''.join([istr.replace('-','') for istr in self.stream_list])+'-'


    def RemoveVals(self,WipeData = True):
        if 'C3' == WipeData:
            self.NJN.RemoveVals()
            if self.isDScomb:  self.NJNsing.RemoveVals()
        elif 'Flow' ==  WipeData:
            self.FO.FlowRemoveVals()
        elif WipeData:
            self.NJN.RemoveVals()
            if self.isDScomb:  self.NJNsing.RemoveVals()
            self.FO.FlowRemoveVals()


    def Read(self,show_timer=False):
        if 'C3' not in self.NJN.C3_cfgs.columns: self.NJN.Read(Full=True,show_timer=show_timer)
        if self.isDScomb and 'C3' not in self.NJNsing.C3_cfgs.columns:
            self.NJNsing.Read(Full=True,show_timer=show_timer)
        if 'Op' not in self.FO.Op_cfgs.columns:
            self.FO.FlowRead(show_timer=show_timer)
        # print str(self.FO.Op_cfgs['Op'])

    # def ReadAndBoot(self):
    #     self.Read()
    #     self.Bootstrap()

    def ImportFitRanges(self,fit_ranges):
        self.fit_ranges = fit_ranges
        if isinstance(self.fit_ranges,str):
            self.fit_ranges = [self.fit_ranges]


    def SetTflowFit(self,tflow):
        # if tflow not in self.FO.tflowlist and tflow != 'Not Set':
        #     print 'Warning: ', tflow, ' not in flow list, tflowfit not added'
        # else:
        self.tflowfit = tflow

    def SetTFit(self,tsink):
        # if tsink not in ['t'+str(it) for it in xrange(1,self.nt+1)] and tsink != 'Not Set':
        #     print 'Warning: ', tsink, ' not in  list, tsinkfit not added'
        # else:
        self.tsinkfit = tsink


    def ImportTwoPtFits(self,state,sourcepars,sinkpars):
        self.NJN.ImportTwoPtFits(state,sourcepars,sinkpars)

    def GetTwoPtPars(self,ip,state,pars):
        self.NJN.GetTwoPtPars(ip,state,pars)


    ## fit_ranges is list of fit ranges to fit the correlator to
    ## formatted as fitr#-# where # are the min and max fit values.
    ## states correspond to how many states to fit to (currenlty only working for states = 1 or 2)
    ## C2pars is list of size 2, containing source and sink fit parameters used in the three point correlator fits
    ## C2pars can be class instances of TwoPointCorr, Fitting, or just a list of values.
    ## Word of caution, theb ootstraps need to have been constructed with the same configurations.
    ## Best to have this run from RatioFunction class only.

    # def Fit(self,C2pars='PreDef',state='PreDef',fit_ranges='PreDef',iGuess='PreDef',EstDir=False,thistflowlist= 'PreDef'):
    #     def RunPredefFits():
    #         if state == 'PreDef' and fit_ranges == 'PreDef':
    #             for izip in self.PredefFits:
    #                 if len(izip) == 2:
    #                     self.Fit(C2pars=C2pars,state=izip[0],fit_ranges=[izip[1]],EstDir=EstDir)
    #                 else:
    #                     self.Fit(C2pars=C2pars,state=izip[0],fit_ranges=[izip[1]],iGuess=izip[2],EstDir=EstDir)
    #         elif state == 'PreDef':
    #             for izip in self.PredefFits:
    #                 if izip[1] == fit_ranges:
    #                     if len(izip) == 2:
    #                         self.Fit(C2pars=C2pars,state=izip[0],fit_ranges=[izip[1]],EstDir=EstDir)
    #                     else:
    #                         self.Fit(C2pars=C2pars,state=izip[0],fit_ranges=[izip[1]],iGuess=izip[2],EstDir=EstDir)
    #         elif fit_ranges == 'PreDef':
    #             for izip in self.PredefFits:
    #                 if izip[0] == state:
    #                     if len(izip) == 2:
    #                         self.Fit(C2pars=C2pars,state=izip[0],fit_ranges=[izip[1]],EstDir=EstDir)
    #                     else:
    #                         self.Fit(C2pars=C2pars,state=izip[0],fit_ranges=[izip[1]],iGuess=izip[2],EstDir=EstDir)
    #     if 'PreDef' in [state,fit_ranges]:
    #         RunPredefFits()
    #         return
    #     if iGuess != 'PreDef': self.iGuess = iGuess
    #     if isinstance(state,basestring):
    #         if 'state' in state:
    #             state = state.replace('state','')
    #         state = int(state)
    #
    #     if C2pars != 'PreDef': self.ImportTwoPtFits(state,C2pars[0],C2pars[1])
    #     self.ImportFitRanges(fit_ranges)
    #     if state == 1:
    #         # thisfitfun = C3OneStateFitFunNoExp
    #         thisfitfun = C3OneStateFitFunNoExp
    #         npar = 1
    #     elif state == 2:
    #         # thisfitfun = C3TwoStateFitFunNoExp
    #         thisfitfun = C3MomTwoStateFitFunNoExp
    #         npar = 4
    #     elif state > 2:
    #         raise IOError('fitting state is only implemented up to 2 state so far. Cannot do '+str(state))
    #
    #     if thistflowlist == 'PreDef':
    #         if self.tflowfit == 'Not Set': raise IOError('tflowfit not set, please set by using .SetTflowFit()')
    #     else:
    #         if thistflowlist == 'All':
    #             self.tflowfit = self.FO.tflowlist
    #         else:
    #             self.tflowfit = thistflowlist
    #
    #     ststr = 'state'+str(state)
    #     legstate = str(state)+'SF'
    #     lFit,lFitA,lFitS = [],[],[]
    #     ilist = []
    #     for (igamma,itflow),gammadata in self.NJNQFull_Stats['boot'].groupby(level=('current_operator','flow_times')):
    #         for (ip,pdata),ptwopt in zip(gammadata.groupby(level='momentum'),self.NJN.TwoPtPars):
    #             for ifit in self.fit_ranges:
    #                 # print 'Fitting state',str(state) ,igamma,itflow,ifit
    #                 ifitmin,ifitmax = map(int,unxmlfitr(ifit))
    #                 ## is this the right ordering of dictionary dimensions?
    #                 tdatarange = [[]]
    #                 ydata = []
    #                 for it,tdata in pdata.iteritems():
    #                     if itflow not in tdata.keys(): continue
    #                     tint = untstr(it)
    #                     if ifitmin <= tint <= ifitmax:
    #                         tcurrboot = BootStrap(self.nboot,bootvals=[np.float64(tint)]*self.nboot)
    #                         tdatarange[0].append(tcurrboot)
    #                         ydata.append(tdata[itflow])
    #                 thistlen = len(tdatarange[0])
    #                 ## this + 3 is just calabrated from looking at the conserved current two- and three- point correlators. Not sure why exactly.
    #                 tsinkboot = BootStrap(self.nboot,bootvals=[np.float64(self.tsink.replace('tsink',''))+3]*self.nboot)
    #                 tdatarange.append([tsinkboot]*thistlen)
    #                 for iparam in ptwopt:
    #                     tdatarange.append([iparam]*thistlen)
    #
    #                 # print 'DEBUG'
    #                 # for it,iy in zip(np.rollaxis(np.array(tdatarange),1),ydata):
    #                 #     print it[0].Avg,it[1].Avg,it[2].Avg,it[3].Avg,it[4].Avg,it[5].Avg,iy.Avg
    #
    #                 if EstDir:
    #                     thisff = ff.Fitting(Funs=[thisfitfun,npar,'Estimate'],data=[tdatarange,ydata],name=' '.join([legstate,ifit,itflow]),iGuess=self.iGuess)
    #                 else:
    #                     thisff = ff.Fitting(Funs=[thisfitfun,npar],data=[tdatarange,ydata],name=' '.join([legstate,ifit,itflow]),iGuess=self.iGuess)
    #                 thisff.FitBoots()
    #
    #                 ilist.append((ststr,igamma,ip,ifit))
    #                 lFit.append(thisff)
    #                 lFitA.append(thisff.fit_data['ParamsAvg'])
    #                 lFitS.append(thisff.fit_data['ParamsStd'])
    #     indicies = pa.MultiIndex.from_tuples(ilist,names=self.NJNQFull_Fit_Col_Names)
    #     if 'boot' in self.NJNQFull_Fit_Stats.columns:
    #         this_df = pa.DataFrame()
    #         this_df.loc[:,'boot'] = pa.Series(lFit,index=indicies)
    #         this_df.loc[:,'Avg'] = pa.Series(lFitA,index=indicies)
    #         this_df.loc[:,'Std'] = pa.Series(lFitS,index=indicies)
    #         self.NJNQFull_Fit_Stats = self.NJNQFull_Fit_Stats.append(this_df)
    #     else:
    #         self.NJNQFull_Fit_Stats.loc[:,'boot'] = pa.Series(lFit,index=indicies)
    #         self.NJNQFull_Fit_Stats.loc[:,'Avg'] = pa.Series(lFitA,index=indicies)
    #         self.NJNQFull_Fit_Stats.loc[:,'Std'] = pa.Series(lFitS,index=indicies)
    #     self.Write()



    def Write(self):
        def FixDictArray(thislist,ikey):
            outDict = OrderedDict()
            for ic,ilist in enumerate(thislist):
                outDict[ikey+str(ic+1)] = ilist
            return outDict

        self.NJN.Write()
        if self.isDScomb: self.NJNsing.Write()
        self.FO.FlowWrite()


        if not hasattr(self,'show_cfgs'):
            self.show_cfgs = True

        outDict = ODNested()

        outDict['ThreePtCorr_File'] = self.NJN.HumanFile
        outDict['FlowOp_File'] = self.FO.flowHumanFile

        excel_params = pa.Series(deepcopy(outDict))

        for col_key,col_data in self.NJNQInt_Stats.items():
            if 'Avg' in col_key or 'Std' in col_key: continue
            for (igamma,ip,it,itflow),idata in col_data.items():
                outDict['NJNQInt_'+col_key][igamma][ip][it][itflow] = AvgStdToFormat(idata.Avg,idata.Std)


        for col_key,col_data in self.NJNQFull_Stats.items():
            if 'Avg' in col_key or 'Std' in col_key: continue
            for (igamma,ip,it,itflow,its),idata in col_data.items():
                outDict['NJNQFull_'+col_key][igamma][ip][it][itflow][its] = AvgStdToFormat(idata.Avg,idata.Std)


        if 'boot' in self.NJNQFull_Fit_Stats:
            for (istate,igamma,ip,itr,itflow),fitdict in self.NJNQFull_Fit_Stats['boot'].items():
                if not isinstance(fitdict,ff.Fitting): continue
                outDict['Fits'][istate][igamma][ip][itr][itflow] = fitdict.GetOutputDict()



        if any(np.array(self.ncfg_list) > 0) and self.show_cfgs:
            for ((istream,iccfg),icfg),ixsrc_list in zip(iter(self.NJNQFull_cfgs['configs'].items()),self.NJNQFull_cfgs['xsrc_list'].values):
                outDict['cfglist'][istream][icfg] = 'n_xsrc_'+str(len(ixsrc_list))



        ## Write human readable file with data
        WriteXml(self.HumanFile,{'Results':outDict})

        if any(np.array(self.ncfg_list) > 0):
            for ((istream,iccfg),icfg),ixsrc_list in zip(iter(self.NJNQFull_cfgs['configs'].items()),self.NJNQFull_cfgs['xsrc_list'].values):
                outDict['cfglist'][istream][icfg] = 'n_xsrc_'+str(len(ixsrc_list))

        WriteExcel(self.ExcelFile,{ 'NJNQFull_results':deepcopy(self.NJNQFull_Stats),
                                    'NJNQInt_results':deepcopy(self.NJNQInt_Stats)},
                                    cfglist=self.NJNQFull_cfgs,params=excel_params)

        ## pickles rest of data for reading
        self.RemoveFuns()
        WritePickle(self.PickleFile,self.__dict__)
        self.GetFuns()
    ##################################

    def GetFuns(self):
        self.FO.GetFuns()

    def RemoveFuns(self):
        self.FO.RemoveFuns()

    def ReadAndWrite(self,WipeData=True,show_timer=False):
        self.Read(show_timer=show_timer)
        self.Bootstrap(WipeData=WipeData,show_timer=show_timer)
        self.Stats()
        # self.Fit()
        self.Write()

    def LoadPickle(self,DefWipe=False,WipeData = True,CheckCfgs=False,CheckMom=True,show_timer=False):
        if os.path.isfile(self.PickleFile) and not DefWipe:
            # print 'Loading Pickle for ' , self.name
            loadeddict = ReadPickleWrap(self.PickleFile)
            checklist = self.thischecklist
            if CheckCfgs: checklist = checklist + ['cfglist']
            if CheckMom: checklist = checklist + ['checkmomlist']
            if CheckClass(self.__dict__,loadeddict,checklist):
                self.__dict__.update(loadeddict)
                if Wipe_All_Fits:
                    self.NJNQFull_Fit_Stats = pa.DataFrame()
            else:
                if os.path.isfile(self.PickleFile+'.bak'):
                    print('using backupfile for '+self.name)
                    self.PickleFile = self.PickleFile+'.bak'
                    self.HumanFile = self.HumanFile+'.bak'
                    self.ExcelFile = self.ExcelFile.replace('.xlsx','.bak.xlsx')
                    self.LoadPickle(DefWipe=DefWipe,WipeData=WipeData,CheckCfgs=CheckCfgs,CheckMom=CheckMom)
                    return
                print('Warning, file ' , self.HumanFile , ' has different parameters than this instance, reading and writing over file now')
                print()
                self.ReadAndWrite(WipeData=WipeData,show_timer=show_timer)
            # self.Write()
        else:
            self.ReadAndWrite(WipeData=WipeData,show_timer=show_timer)



    ## Comparisons

    def __str__(self):
        return '_'.join(self.name)


    def __iter__(self):
        return self

    def iteritems(self):
        return iter(self.NJNQFull_Stats['boot'].items())

    def items(self):
        return list(self.NJNQFull_Stats['boot'].items())

    def itemsAvgStd(self):
        outlist = []
        for ivals in self.NJNQFull_Stats.itertuples():
            outlist.append((ivals[0],ivals[2],ivals[3]))
        return outlist

    def values(self):
        return self.NJNQFull_Stats['boot'].values

    def itervalues(self):
        return iter(self.NJNQFull_Stats['boot'].values)

    def valuesAvgStd(self):
        outlist = []
        for ivals in self.NJNQFull_Stats.itertuples():
            outlist.append((ivals[1],ivals[2]))
        return outlist

    def keys(self):
        return self.NJNQFull_Stats.index


    def __setitem__(self, key, value ):
        self.NJNQFull_Stats.loc[key,'boot'] = value

    def __getitem__(self, key):
        return self.NJNQFull_Stats.loc[key,'boot']

    def __reversed__(self):
        self.NJNQFull_Stats = self.NJNQFull_Stats.iiloc[::-1]
        self.NJNQInt_Stats = self.NJNQInt_Stats.iiloc[::-1]
        return self



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
        for ival in self.keys():
            self[ival] = -self[ival]
        return self


    def __pos__(self):
        return self

    def __abs__(self):
        for ival in self.keys():
            self[ival] = abs(self[ival])
        return self

    ## Casting, just casts the bootstrap values, ignore return value (only needed since return has to be same type).


    def __complex__(self):
        for ival in self.keys():
            complex(self[ival])
        self.Stats()
        return 1+0j

    def __int__(self):
        for ival in self.keys():
            int(self[ival])
        self.Stats()
        return 1

    def __float__(self):
        for ival in self.keys():
            np.float64(self[ival])
        self.Stats()
        return 1.0




def Test3pt(DefWipe=False):
    from Params import defInfoFlow
    defInfoFlow['gammas'] = ['g4']
    defInfoFlow['qmom'] = ['0 0 0']
    defInfoFlow['Projector'] = 'P4'
    data = ThreePointCorr(Info=defInfoFlow)
    data.LoadPickle(DefWipe=DefWipe)
    return data

def TestNJNQ(DefWipe=False):
    from Params import defInfoFlow
    tflow = 6.01
    defInfoFlow['tflowlist'] = np.array([tflow])
    defInfoFlow['tflowfit'] = np.array([tflow])
    defInfoFlow['gammas'] = ['g4']
    defInfoFlow['qmom'] = ['0 0 1']
    defInfoFlow['Projector'] = 'P3'
    data = NJNQCorr(Info=defInfoFlow)
    data.LoadPickle(DefWipe=DefWipe)
    return data

def TestNJNQFull(DefWipe=False):
    from Params import defInfoFlow
    tsum = 32
    tflow = 6.01
    defInfoFlow['tsum_fit_list'] =  list(map(tsumstr,list(range(tsum))))
    defInfoFlow['boot_tsum_list'] =  list(map(tsumstr,list(range(tsum))))
    defInfoFlow['tflowlist'] = np.array([tflow])
    defInfoFlow['tflowfit'] = np.array([tflow])
    defInfoFlow['op_trange_sum_type_3pt'] = 'from_src_sym'
    defInfoFlow['gammas'] = ['g4']
    defInfoFlow['qmom'] = ['0 0 1']
    defInfoFlow['Projector'] = 'P3'
    data = NJNQFullCorr(Info=defInfoFlow)
    data.LoadPickle(DefWipe=DefWipe,show_timer=True)
    return data


if __name__ == '__main__':
    data3 = Test3pt()
    dataNJNQ = TestNJNQ()
    # dataFull = TestNJNQFull()
    print('Testing 3ptcorr complete, see data3 and dataNJNQ variables')
