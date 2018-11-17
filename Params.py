#!/usr/bin/env python
## THIS FILE IS TEMPORARY, FOR DEBUGGING AND SIMPLE DEFAULT PARAMETERS
## TODO: either default parameters in file, or even GUI?...

from collections import OrderedDict
from MiscFuns import mkdir_p
import numpy as np
from PredefFitFuns import ConstantFitFun
from copy import deepcopy
import pickle as pickle
# import dill as pickle
import os
import socket

import sys
if sys.version_info[0] < 3:
    raise Exception("This branch (Py3 or offshoot) is built with python3, please use older branch if python2 is needed")


this_machine = socket.gethostname()
guess_dir = os.getcwd()

def Def_pdict():
    paramdict = OrderedDict()
    paramdict['Multicore'] = False
    paramdict['Debug_Mode'] = False
    paramdict['num_of_procs'] = 4
    paramdict['max_job_per_thread'] = 20
    paramdict['write_excel'] = True
    paramdict['Show_Write'] = False
    paramdict['Show_Read'] = False
    paramdict['Wipe_All_Fits'] = False
    paramdict['cfg_file_type'] = 'msgpack'
    paramdict['Max_Iters'] = 10000
    paramdict['Solver_Prec'] = 10**(-8)
    # paramdict['data_directory'] = '/mnt/research/lqcd/CfunAnalysis/'
    # paramdict['cfun_prefix'] = '/cfunPChroma/'
    paramdict['data_directory'] = guess_dir+'/Data/'
    paramdict['cfg_dir'] = guess_dir+'/Configs/'
    paramdict['formatted_cfg_dir'] = guess_dir+'/FmtConfigs/'
    paramdict['cfun_prefix'] = '/cfuns/'
    paramdict['n_space'] = 32
    paramdict['n_time'] = 64
    paramdict['lattice_spacing'] = 0.0907 ## In fermi
    paramdict['n_boot'] = 200
    paramdict['nchi_threshold'] = 10
    paramdict['nchi_fit_threshold'] = 10
    paramdict['plot_the_ff'] = True
    paramdict['normed_histograms'] = True
    paramdict['def_min_fitr'] = 3
    return paramdict

try:
    if os.path.isfile('./Param_p3.py3p'):
        this_file = os.getcwd()+'/Param_p3.py3p'
        with open( this_file, "rb" ) as pfile:
            print('Reading Params from, ', this_file)
            paramdict =  pickle.load(pfile)
    elif os.path.isfile('../Param_p3.py3p'):
        this_file = os.getcwd()+'/../Param_p3.py3p'
        with open( this_file, "rb" ) as pfile:
            print('Reading Params from, ', this_file)
            paramdict =  pickle.load(pfile)
    elif len(guess_dir) > 0 and os.path.isfile(guess_dir+'Param_p3.py3p'):
        this_file = guess_dir+'Param_p3.py3p'
        with open( this_file, "rb" ) as pfile:
            print('Reading Params from, ', this_file)
            paramdict =  pickle.load(pfile)
    else:
        paramdict = Def_pdict()
        # raise IOError('Param_p3.p is not set, need to deal with this?')
except Exception as err:
    print('Error reading ' , this_file, ' Setting to default values.')
    print('The error was:')
    print(str(err))
    paramdict = Def_pdict()

if os.path.isfile('./JobParam_p3.py3p'):
    with open( './JobParam_p3.py3p', "rb" ) as pfile:
        def_job_params =  pickle.load(pfile)
elif os.path.isfile('../JobParam_p3.py3p'):
    with open( '../JobParam_p3.py3p', "rb" ) as pfile:
        def_job_params =  pickle.load(pfile)
else:
    def_job_params = OrderedDict()
    def_job_params['machine'] = this_machine
    def_job_params['python_exe'] = 'python'
    def_job_params['time'] = '4:50:00'
    def_job_params['email'] = '@gmail.com'
    def_job_params['nproc'] = 1
    def_job_params['RunPTG'] = False
    if this_machine == 'juqueen':
        def_job_params['quetype'] = 'bluegene'
        def_job_params['Scom'] = 'llsubmit'
    elif this_machine == 'laconia':
        def_job_params['quetype'] = 'dev-intel16'
        def_job_params['RunPTG'] = True
        def_job_params['Scom'] = 'qsub'
    else:
        def_job_params['quetype'] = 'batch'
        def_job_params['Scom'] = 'sbatch'
    def_job_params['memory'] = '12GB'


if os.path.isfile('./JobParamSolve.py3p'):
    with open( './JobParamSolve.py3p', "rb" ) as pfile:
        def_job_params_solve =  pickle.load(pfile)
elif os.path.isfile('../JobParamSolve.py3p'):
    with open( '../JobParamSolve.py3p', "rb" ) as pfile:
        def_job_params_solve =  pickle.load(pfile)
else:
    def_job_params_solve = OrderedDict()
    def_job_params_solve['machine'] = this_machine
    def_job_params_solve['python_exe'] = 'python'
    def_job_params_solve['time'] = '4:50:00'
    def_job_params_solve['email'] = '@gmail.com'
    def_job_params_solve['nproc'] = 1
    def_job_params_solve['RunPTG'] = False
    if this_machine == 'juqueen':
        def_job_params_solve['quetype'] = 'bluegene'
        def_job_params_solve['Scom'] = 'llsubmit'
    elif this_machine == 'laconia':
        def_job_params_solve['quetype'] = 'dev-intel16'
        def_job_params_solve['RunPTG'] = True
        def_job_params_solve['Scom'] = 'qsub'
    else:
        def_job_params_solve['quetype'] = 'batch'
        def_job_params_solve['Scom'] = 'sbatch'
    def_job_params_solve['memory'] = '12GB'

All_Classes_Names = [ 'AResAlpha','AResChi2','CalculationParams','RunParams','PlotParams',
    'JobParams','JobParamsSolve','TPCorr','FlowedOp','FlowedOpFull','RatioFun',
    'Alpha','AlphaFull','FormFacts','DictFlowEditor','DictEditor','BeginFrame',
    'AutoCorrelate','BootStrap','AlphaFrame','FlowOpFrame','BeginFrame','FFEquation',
    'Fitting','FlowOp','FlowOpFullSquared','FormFact','GammaMat','CorrSpinTrace','DOO',
    'DataFrame','InfoFrame','BeginFrame','LatticeParameters','Plotting','RatioCorr','RatioFOCorr',
    'RatioFOFullCorr','ReadFSCfunPickCHROMA','RC3Full','R2CChromaXMLFileList','RC2Full','NaNCfunError',
    'ReadFSCfunPickCHROMA','RC3Full','R2CChromaXMLFileList','RC2Full','NaNCfunError','SetOfFF',
    'SetOfFitFuns','SetOfFO','SetOfFOFull','SetOfRat','SetOfTwoPt','SetOfNNQ','SetOfNNQFull','ThreePointCorr',
    'NJNQCorr','NJNQFullCorr','Timer','TwoPointCorr','NNQCorr','NNQFullCorr','FlowedTwoPtCorr','VariationalTwoPt',
    'LegendFmt','KeyForamtting','ParsingInterrupted','_DictSAXHandler']

# G2_abs = False
# if G2_abs:
#     print 'G2_abs is ON! WARNING!'
myeps = np.finfo(0.0).eps
percentile = 0.158
Qeps = 10**(-5)
defMcore = False ## Not Needed
# defMcore = paramdict['Multicore'] ## changes all the default multicore LoadPickles
# defnProc = 12 ## default max processors to use
defnProc = paramdict['num_of_procs'] ## default max processors to use
defMaxTasks = paramdict['max_job_per_thread']
ShowXmlWrite = paramdict['Show_Write']
if 'cfg_file_type' in list(paramdict.keys()):
    cfg_file_type = paramdict['cfg_file_type']
else:
    cfg_file_type = 'msgpack'

if 'write_excel' not in list(paramdict.keys()):
    write_excel = True
else:
    write_excel = paramdict['write_excel']
# if 'Show_Read' not in paramdict.keys():
#     ShowRead = False
# else:
ShowRead = paramdict['Show_Read']
Wipe_All_Fits = False
if 'Wipe_All_Fits' in list(paramdict.keys()):
    Wipe_All_Fits = paramdict['Wipe_All_Fits']
plot_style = 'classic'

colourset8 = [ '#000080','#B22222','#00B800','#8B008B', '#0000FF','#FF0000','#000000','#32CD32','#FF0066']
markerset = ['o','s','^','v','d']

lineres = 100

defMaxIter = paramdict['Max_Iters']
defPrec = paramdict['Solver_Prec']## max number of iterations for least squares fitting.
fillalpha = 0.5
defSparam = 1.5

shiftmax = 3
shiftper = 0.05 ##R function shift
shiftperff = 0.01 ##Form Factor shift
shiftset = [0]
for ish in np.arange(1,shiftmax+1): shiftset += [-ish*shiftper,ish*shiftper]


## CURRENTLY THIS NEEDS TO CHANGE FOR YOUR DATA DIRECTORIES
# datadir = '/mnt/research/lqcd/CfunAnalysis/'
# cfundir = datadir + '/cfunPChroma/'
datadir = paramdict['data_directory']
cfg_dir = paramdict['cfg_dir']
if 'formatted_cfg_dir' in paramdict:
    cfgfmtdir = paramdict['formatted_cfg_dir']
else:
    cfgfmtdir = guess_dir+'/FmtConfigs/'
if not os.path.isdir(cfgfmtdir):
    print('Warning, config directory not found:')
    print(cfgfmtdir)
    print(' please change formatted_cfg_dir in General Parameters')

if 'output_directory' in paramdict:
    outputdir = paramdict['output_directory']
else:
    outputdir = guess_dir+'/results/'
if 'scratch_directory' in paramdict:
    scratchdir = paramdict['scratch_directory']
else:
    scratchdir = outputdir+'/scratch/'

test_file = scratchdir + '/debug_file.test'
if 'flow_new_format' in paramdict:
    flow_new_format = paramdict['flow_new_format']
else:
    flow_new_format = True
functiondir = scratchdir + '/functions/'
cfundir = cfg_dir + paramdict['cfun_prefix']
########################################################

Qdir = cfg_dir+'/topcharge/'
Wdir = cfg_dir+'/weinopp/'
graphdir = outputdir + '/graphs/'

try:
    mkdir_p(functiondir)
    mkdir_p(graphdir)
    mkdir_p(scratchdir)
    mkdir_p(functiondir)
except:
    print(os.getcwd())
    if os.path.isfile('./Param_p3.py3p'):
        print('./Param_p3.p found')
    elif os.path.isfile('../Param_p3.py3p'):
        print('../Param_p3.p found')
    elif len(guess_dir) > 0 and os.path.isfile(guess_dir+'Param_p3.py3p'):
        print(guess_dir+'Param_p3.p found')
    else:
        print('No parameter file found')
    errorstr = (graphdir +
                '\n error making directory, machine identified as \n '+
                this_machine)
    print(errorstr)
    print('THIS MUST BE FIXED BEFORE RUNNING')
# MomSqrdSet = [ '0','1','2','3','4','5','6','7','8','9' ]
MomSqrdSet = [ '0','1','2','3','4']
defmom2list = list(range(5))
defqsqrdMax = np.max(defmom2list)
defpmomlist = ['0 0 0','1 0 0','1 1 0','1 1 1','2 0 0']
nxyz = paramdict['n_space']
nt = paramdict['n_time']
nxyzt = [nxyz,nxyz,nxyz,nt]
# qunit = (2.0*np.pi)/float(nxyz)
hbarc = 0.1973269718 ## In GeV * fermi
cm_to_fm = 10**13
alpha_guess = -0.5 ## guess for what the nucleon mixing angle alpha is......
                   ## doesn't really matter what this is, just change it to see if Form Factors change.



latspace = paramdict['lattice_spacing'] ## In fermi
# hbarcdivlat = hbarc/latspace

defxlim = [4,10]
defxlimAlpha = [1,15]
defxlimAlphaTflow = [0,1] ## in sqrt(8t)
defxlimOpMore = np.arange(0.01,9.11,0.1)
defxlimOpMost = np.arange(0.01,10.00,0.01)
# defxlimOp = np.arange(0.01,10,1.0)
defxlimOp = np.array([4.01,5.01,6.01])
# defxlimOp = np.array([6.01])

def chitcoeff(nt,nxyz,latspace,obs):
    if 'QW' in obs or 'WQ' in obs:
        return  (hbarc/(latspace**(4.)*nxyz**(3.)*nt)) ## (1/V) [GeV]^{6}
    elif 'Q' in obs:
        return (hbarc/(latspace*nxyz**(0.75)*nt**(0.25))) ## (1/V)^1/4 [GeV]
    elif 'W' in obs:
        return (hbarc/((latspace*nxyz**(0.75)*nt**(0.25)))**0.5) ## (1/V)^1/8 [GeV]

nboot = paramdict['n_boot']
if 'nchi_threshold' in list(paramdict.keys()):
    nchi_threshold = paramdict['nchi_threshold']
else:
    nchi_threshold = 10

if 'nchi_fit_threshold' in list(paramdict.keys()):
    nchi_fit_threshold = paramdict['nchi_fit_threshold']
else:
    nchi_fit_threshold = 10

if 'plot_the_ff' in list(paramdict.keys()):
    plot_the_ff = paramdict['plot_the_ff']
else:
    plot_the_ff = True

if 'normed_histograms' in list(paramdict.keys()):
    normed_histograms = paramdict['normed_histograms']
else:
    normed_histograms = True

if 'def_min_fitr' in list(paramdict.keys()):
    def_min_fitr = paramdict['def_min_fitr']
else:
    def_min_fitr = 3


if 'Debug_Mode' in list(paramdict.keys()):
    Debug_Mode = paramdict['Debug_Mode']
else:
    Debug_Mode = False



MesICS1 = [10**-6,np.log(0.2)]
MesICS2 = [10**-6,10,np.log(0.2),np.log(0.2)]

defInfo = OrderedDict()
defInfo['outdir'] = outputdir
defInfo['show_cfgs'] = True
defInfo['nxyzt'] = [32,32,32,64]
defInfo['a'] = 0.0907 ## In fermi
defInfo['kud'] = 1375400
defInfo['ks'] = 1364000
defInfo['Csw'] = 1715
defInfo['Beta'] = 1900
defInfo['jsm'] = 64
defInfo['sm2ptRat'] = [16,32,64] ## remove this, jsmcoeffs and sinklab to get rid of sink combinations
if defInfo['kud'] == 1375400:
    defInfo['jsmCoeffs'] = np.array([0.0000005522,-0.0001589143,0.9999999874],np.float64)
    defInfo['sinklab'] = 'PoF0D1'
else:
    defInfo['jsmCoeffs'] = [0.0,0.0,1.0]
    defInfo['sinklab'] = 'jsm64'


defInfo['autocorr'] = True
defInfo['Sparam'] = 1.5
defInfo['MesOrBar'] = 'Baryon'
defInfo['Adjoint_Flow'] = True
defInfo['combine_sinks'] = False
defInfo['ism'] = 64
defInfo['tsrc'] = 0
defInfo['Fits'] = [('state1','fitr9-20'),('state2','fitr5-20')]
defInfo['tsink'] = 13
# defInfo['pmom'] = ['0 0 0','1 1 1']
defInfo['pmom'] = 'All' ## corresponds to defmom2list (above)
defInfo['ppmom'] = '0 0 0'
# defInfo['qmom'] = ['0 0 0','2 0 0']
defInfo['qmom'] = ['0 0 0']
# defInfo['qmom'] = ['0 0 0','0 0 1','0 1 1','1 1 1','0 0 2']
defInfo['DoubSing'] = 'doub'
# defInfo['DoubSing'] = 'Neutron'
# defInfo['DoubSing'] = 'Proton'
defInfo['Projector'] = 'P4'
# defInfo['gammas'] = ['g4Cons']
defInfo['gammas'] = ['g4']
# defInfo['Projector'] = 'P3'
# defInfo['gammas'] = ['g3g5cmplx']
defInfo['Interp'] = 'CPEven'



defInfo['Fits3pt'] = [('state1','fitr4-8'),('state2','fitr3-9')]
defInfo['FitsRF'] = ['fitr4-8']
defInfo['FitsRFFull'] = [['taufitr3-8','fitr4-8']]
defInfo['FFFits'] = ['All']

defInfo2 = deepcopy(defInfo)
defInfo2['ism'] = 32

defInfo3 = deepcopy(defInfo)
defInfo3['ism'] = 16


defInfoFlow = deepcopy(defInfo)
defInfoFlow['Interp'] = 'CPOdd'
defInfoFlow['Observable'] = 'TopCharge'
defInfoFlow['tflowlist'] = defxlimOp
# defInfoFlow['nrandt'] = 5
defInfoFlow['nrandt'] = 64
defInfoFlow['tsum_list'] = ['ts63']
defInfoFlow['FlowFits'] = ['fitr8-10']
defInfoFlow['FlowFitFun'] = (ConstantFitFun,1)
defInfoFlow['FitsAlpha'] = ['fitr10-20']
defInfoFlow['FitsAlphaFull'] = [('fitr10-20','tsumfitr8-32')]
defInfoFlow['FitFunAlpha'] = (ConstantFitFun,1)
defInfoFlow['tflowfit'] = ['t_f6.01','t_f5.01','t_f4.01']
defInfoFlow['tsinkfit'] = ['t10']
defInfoFlow['Projector'] = 'P3'
defInfoFlow['gammas'] = ['g4']
defInfoFlow['qmom'] = ['0 0 1']
defInfoFlow['tsum_picked'] = 'ts63'
defInfoFlow['op_trange_sum_type'] = 'from_src_None'

# defInfoFlow['tsum_fit_list'] = ['ts0','ts16']
defInfoFlow['tsum_fit_list'] = ['ts'+str(its) for its in range(0,33)]
defInfoFlow['boot_tsum_list'] = defInfoFlow['tsum_fit_list']
defInfoFlow['boot_nt_list'] = list(range(nt//2))
defInfoFlow['op_trange_sum_type'] = 'from_src_sum_sym_None'
defInfoFlow['fit_par_picked'] = 'First'
defInfoFlow['tsum_spec_fit_list'] = ['tss32']
# defInfoFlow['fit_min_picked'] = float('NaN')

# defInfoFlow['FitFun'] = (ChitFitFun,3)


defInfoFlow2 = deepcopy(defInfoFlow)
defInfoFlow2['ism'] = 32

defInfoFlow3 = deepcopy(defInfoFlow)
defInfoFlow3['ism'] = 16



defInfoMes = deepcopy(defInfo)
defInfoMes['MesOrBar'] = 'Meson'
defInfoMes['Interp'] = 'Pion' ## Pion is 15
defInfoMes['Fits'] = [('state1','fitr9-20',MesICS1),('state2','fitr5-20',MesICS2)]


defInfoMes2 = deepcopy(defInfoMes)
defInfoMes2['ism'] = 32

defInfoMes3 = deepcopy(defInfoMes)
defInfoMes3['ism'] = 16


defFitDict = OrderedDict()
defFitDict['Kud01375400Ks01364000_-a-_Baryon_CPEven_tsrc0_ism64_jsm64'] = ('state1','fitr10-20')
defFitDict['Kud01375400Ks01364000_-a-_Baryon_CPEven_tsrc0_ism32_jsm64'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_-a-_Baryon_CPEven_tsrc0_ism16_jsm64'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_-a-_Baryon_CPEven_tsrc0_ism64_PoF0D1'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_-a-_Baryon_CPEven_tsrc0_ism32_PoF0D1'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_-a-_Baryon_CPEven_tsrc0_ism16_PoF0D1'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_-a-_Baryon_CPOdd_tsrc0_ism64_jsm64'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_-a-_Baryon_CPOdd_tsrc0_ism32_jsm64'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_-a-_Baryon_CPOdd_tsrc0_ism16_jsm64'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_-a-_Baryon_CPOdd_tsrc0_ism64_PoF0D1'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_-a-_Baryon_CPOdd_tsrc0_ism32_PoF0D1'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_-a-_Baryon_CPOdd_tsrc0_ism16_PoF0D1'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_Meson_Pion_tsrc0_ism64_jsm64'] = ('state1','fitr9-20',MesICS1)
defFitDict['Kud01375400Ks01364000_Meson_Pion_tsrc0_ism32_jsm64'] = ('state1','fitr9-20',MesICS1)
defFitDict['Kud01375400Ks01364000_Meson_Pion_tsrc0_ism16_jsm64'] = ('state1','fitr9-20',MesICS1)
defFitDict['Kud01375400Ks01364000_-a-_P4_doub_pp000_tsrc0_tsink13_ism64_PoF0D1'] = 'fitr4-8'
defFitDict['Kud01375400Ks01364000_-a-_P4_doub_pp000_tsrc0_tsink13_ism32_PoF0D1'] = 'fitr4-8'
defFitDict['Kud01375400Ks01364000_-a-_P4_doub_pp000_tsrc0_tsink13_ism16_PoF0D1'] = 'fitr4-8'
defFitDict['Kud01375400Ks01364000_-a-_P3_doub_pp000_tsrc0_tsink13_ism64_PoF0D1'] = 'fitr4-8'
defFitDict['Kud01375400Ks01364000_-a-_P3_doub_pp000_tsrc0_tsink13_ism32_PoF0D1'] = 'fitr4-8'
defFitDict['Kud01375400Ks01364000_-a-_P3_doub_pp000_tsrc0_tsink13_ism16_PoF0D1'] = 'fitr4-8'
defFitDict['MassForFF'] = defFitDict['Kud01375400Ks01364000_-a-_Baryon_CPEven_tsrc0_ism64_jsm64']


defFitAllGamma = 'fitr4-8'

defFitDict2 = OrderedDict()
defFitDict2['Baryon_CPEven_tsrc0_ism64_jsm64'] = ('state2','fitr5-20')
defFitDict2['Baryon_CPEven_tsrc0_ism32_jsm64'] = ('state2','fitr5-20')
defFitDict2['Baryon_CPEven_tsrc0_ism16_jsm64'] = ('state2','fitr5-20')
defFitDict2['Meson_Pion_tsrc0_ism64_jsm64'] = ('state2','fitr5-20',MesICS2)
defFitDict2['Meson_Pion_tsrc0_ism32_jsm64'] = ('state2','fitr5-20',MesICS2)
defFitDict2['Meson_Pion_tsrc0_ism16_jsm64'] = ('state2','fitr5-20',MesICS2)


##HERE TODO
defFitAlpha = OrderedDict()
defFitAlpha['tsrc0_ism64_jsm64_TopCharge'] = 'fitr8-10'
defFitAlpha['tsrc0_ism32_jsm32_TopCharge'] = 'fitr8-10'
defFitAlpha['tsrc0_ism16_jsm16_TopCharge'] = 'fitr8-10'
defFitAlpha['RC32x64_Kud01375400Ks01364000_-a-_tsrc0_ism64_PoF0D1_TopCharge'] = 'fitr8-10'
defFitAlpha['Kud01375400Ks01364000_-a-_tsrc0_ism64_PoF0D1_TopCharge'] = 'fitr8-10'
defFitAlpha['Kud01375400Ks01364000_-a-_tsrc0_ism64_jsm64_TopCharge'] = 'fitr8-10'
defFitAlpha['Kud01375400Ks01364000_-a-_tsrc0_ism32_jsm64_TopCharge'] = 'fitr8-10'
defFitAlpha['Kud01375400Ks01364000_-a-_tsrc0_ism16_jsm64_TopCharge'] = 'fitr8-10'






defcfglist = OrderedDict()
defcfglist['-a-004310'] = list(map(str,[101,102]))
defcfglist['-a-004010'] = list(map(str,[101,102]))
defcfglist['-a-003310'] = list(map(str,[101,102]))


defSetsOfInfos = [defInfo,defInfo2,defInfo3]


## used to fix the stupid right add/sub/mult etc.. problem when overloading operators
def ParentClasses(thisclass,parentclasses):
    out = False
    for iparent in parentclasses:
        out = out or isinstance(thisclass,iparent)
    return out

## turns list of tflows into sqrt(8*tflow)
def TflowToPhys(tflowlist,thislatspace):
    if isinstance(tflowlist, list) or isinstance(tflowlist,np.ndarray):
        output = []
        for iflow in tflowlist:
            if isinstance(iflow,str):
                output.append(np.sqrt(8*np.float(iflow.replace('t_f','')))*thislatspace)
            else:
                output.append(np.sqrt(8*iflow)*thislatspace)
        return output ## in fermi
    elif isinstance(tflowlist,str):
        return np.sqrt(8*np.float(tflowlist.replace('t_f','')))*thislatspace
    else:
        return np.sqrt(8*tflowlist)*thislatspace

## CURRENLTY LOCKING DOES NOT WORK WITH MULTICORE....
