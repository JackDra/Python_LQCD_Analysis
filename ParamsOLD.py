#!/usr/bin/env python
## THIS FILE IS TEMPORARY, FOR DEBUGGING AND SIMPLE DEFAULT PARAMETERS
## TODO: either default parameters in file, or even GUI?...

from collections import OrderedDict
from MiscFuns import mkdir_p
import numpy as np
from PredefFitFuns import ConstantFitFun
from copy import deepcopy
from xmltodict import unparse
import dill as pickle

with open( './Param_p3.py3p', "rb" ) as pfile:
    paramdict =  pickle.load(pfile)


myeps = np.finfo(0.0).eps


defMcore = False ## changes all the default multicore LoadPickles
# defnProc = 12 ## default max processors to use 
defnProc = 4 ## default max processors to use 
ShowXmlWrite = True


colourset8 = [ '#000080','#B22222','#00B800','#8B008B', '#0000FF','#FF0000','#000000','#32CD32','#FF0066']
markerset = ['o','s','^','v','d']

lineres = 100

defMaxIter = 10000 ## max number of iterations for least squares fitting.
defPrec = 10**(-8) ## max number of iterations for least squares fitting.
fillalpha = 0.5
defMaxTasks = defnProc*5


shiftmax = 3
shiftper = 0.05 ##R function shift
shiftperff = 0.01 ##Form Factor shift
shiftset = [0]
for ish in np.arange(1,shiftmax+1): shiftset += [-ish*shiftper,ish*shiftper]


## CURRENTLY THIS NEEDS TO CHANGE FOR YOUR DATA DIRECTORIES
# datadir = '/mnt/research/lqcd/CfunAnalysis/'
# cfundir = datadir + '/cfunPChroma/'
datadir = '/home/jackdra/PHD/CHROMA/TestVar/scratch/'
cfundir = datadir + '/cfunsg5/'
########################################################

Qdir = datadir+'/topcharge/'
Wdir = datadir+'/weinopp/'
outputdir = datadir+'/resultsREWORK/'
graphdir = outputdir + '/graphs/'


mkdir_p(graphdir)

# MomSqrdSet = [ '0','1','2','3','4','5','6','7','8','9' ]
MomSqrdSet = [ '0','1','2','3','4']
defmom2list = list(range(5))
defqsqrdMax = np.max(defmom2list)
nxyz = 32
nt = 64
nxyzt = [nxyz,nxyz,nxyz,nt]
# qunit = (2.0*np.pi)/float(nxyz)
hbarc = 0.1973269718 ## In GeV * fermi
cm_to_fm = 10**13
alpha_guess = -0.5 ## guess for what the nucleon mixing angle alpha is......
                   ## doesn't really matter what this is, just change it to see if Form Factors change.

                   

latspace = 0.0907 ## In fermi
# hbarcdivlat = hbarc/latspace

defxlim = [2,20]
defxlimAlpha = [1,10]
defxlimAlphaTflow = [0,1] ## in sqrt(8t)
# defxlimOp = np.arange(0,10,0.1)
defxlimOp = np.arange(0,10,1.0)

def chitcoeff(nt,nxyz,latspace,obs):
    if 'QW' in obs or 'WQ' in obs:
        return  (hbarc/(latspace*nxyz**(0.75)*nt**(0.25))) ## (1/V)^(1/6) [GeV]
    elif 'Q' in obs:
        return (hbarc/(latspace*nxyz**(0.75)*nt**(0.25))) ## (1/V)^1/4 [GeV]
    elif 'W' in obs:
        return (hbarc/((latspace*nxyz**(0.75)*nt**(0.25)))**0.5) ## (1/V)^1/8 [GeV]

nboot = 200
# nboot = 20
print('nboot set to ' ,nboot)





            
MesICS1 = [10**-6,np.log(0.2)]
MesICS2 = [10**-6,10,np.log(0.2),np.log(0.2)]

defInfo = OrderedDict()
defInfo['outdir'] = outputdir
defInfo['nxyzt'] = [32,32,32,64]
defInfo['a'] = 0.0907 ## In fermi
defInfo['kud'] = 1375400
defInfo['ks'] = 1364000
defInfo['jsm'] = 64
defInfo['sm2ptRat'] = [16,32,64] ## remove this, jsmcoeffs and sinklab to get rid of sink combinations
if defInfo['kud'] == 1375400:
    defInfo['jsmCoeffs'] = np.array([0.0000005522,-0.0001589143,0.9999999874],np.float64)
    defInfo['sinklab'] = 'PoF0D1'
else:
    defInfo['jsmCoeffs'] = [0.0,0.0,1.0]
    defInfo['sinklab'] = 'jsm64'

    
defInfo['ism'] = 64
defInfo['tsrc'] = 0
defInfo['Fits'] = [('state1','fitr9-20'),('state2','fitr5-20')]
defInfo['tsink'] = 13
# defInfo['pmom'] = ['0 0 0','1 1 1']
defInfo['pmom'] = 'All' ## corresponds to defmom2list (above)
defInfo['ppmom'] = '0 0 0'
# defInfo['qmom'] = ['0 0 0','2 0 0']
defInfo['qmom'] = ['0 0 0']
# defInfo['DoubSing'] = 'doub'
defInfo['DoubSing'] = 'Neutron'
# defInfo['DoubSing'] = 'Proton'
defInfo['Projector'] = 'P4'
# defInfo['gammas'] = ['g4Cons']
defInfo['gammas'] = ['g4']
# defInfo['Projector'] = 'P3'
# defInfo['gammas'] = ['g3g5cmplx']

defInfo['Fits3pt'] = [('state1','fitr4-8'),('state2','fitr3-9')]
defInfo['FitsRF'] = ['fitr4-8']
defInfo['FFFits'] = ['All']

defInfo2 = deepcopy(defInfo)
defInfo2['ism'] = 32

defInfo3 = deepcopy(defInfo)
defInfo3['ism'] = 16


defInfoFlow = deepcopy(defInfo)
defInfoFlow['Interp'] = 'CPOdd'
defInfoFlow['Observable'] = 'TopCharge'
defInfoFlow['tflowlist'] = defxlimOp
defInfoFlow['FlowFits'] = ['fitr8-10']
defInfoFlow['FlowFitFun'] = (ConstantFitFun,1)
defInfoFlow['FitsAlpha'] = ['fitr8-10']
defInfoFlow['FitFunAlpha'] = (ConstantFitFun,1)
defInfoFlow['tflowfit'] = ['t_f4.0','t_f5.0','t_f6.0']
defInfoFlow['tsinkfit'] = ['t5']
defInfoFlow['Projector'] = 'P3'
defInfoFlow['gammas'] = ['g4']
defInfoFlow['qmom'] = ['0 0 1']
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
defFitDict['Kud01375400Ks01364000_Baryon_CPEven_tsrc0_ism64_jsm64'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_Baryon_CPEven_tsrc0_ism32_jsm64'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_Baryon_CPEven_tsrc0_ism16_jsm64'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_Baryon_CPEven_tsrc0_ism64_PoF0D1'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_Baryon_CPEven_tsrc0_ism32_PoF0D1'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_Baryon_CPEven_tsrc0_ism16_PoF0D1'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_Baryon_CPOdd_tsrc0_ism64_jsm64'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_Baryon_CPOdd_tsrc0_ism32_jsm64'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_Baryon_CPOdd_tsrc0_ism16_jsm64'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_Baryon_CPOdd_tsrc0_ism64_PoF0D1'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_Baryon_CPOdd_tsrc0_ism32_PoF0D1'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_Baryon_CPOdd_tsrc0_ism16_PoF0D1'] = ('state1','fitr9-20')
defFitDict['Kud01375400Ks01364000_Meson_Pion_tsrc0_ism64_jsm64'] = ('state1','fitr9-20',MesICS1)
defFitDict['Kud01375400Ks01364000_Meson_Pion_tsrc0_ism32_jsm64'] = ('state1','fitr9-20',MesICS1)
defFitDict['Kud01375400Ks01364000_Meson_Pion_tsrc0_ism16_jsm64'] = ('state1','fitr9-20',MesICS1)
defFitDict['Kud01375400Ks01364000_P4_doub_pp000_tsrc0_tsink13_ism64_PoF0D1'] = 'fitr4-8'
defFitDict['Kud01375400Ks01364000_P4_doub_pp000_tsrc0_tsink13_ism32_PoF0D1'] = 'fitr4-8'
defFitDict['Kud01375400Ks01364000_P4_doub_pp000_tsrc0_tsink13_ism16_PoF0D1'] = 'fitr4-8'
defFitDict['Kud01375400Ks01364000_P3_doub_pp000_tsrc0_tsink13_ism64_PoF0D1'] = 'fitr4-8'
defFitDict['Kud01375400Ks01364000_P3_doub_pp000_tsrc0_tsink13_ism32_PoF0D1'] = 'fitr4-8'
defFitDict['Kud01375400Ks01364000_P3_doub_pp000_tsrc0_tsink13_ism16_PoF0D1'] = 'fitr4-8'
defFitDict['MassForFF'] = defFitDict['Kud01375400Ks01364000_Baryon_CPEven_tsrc0_ism64_jsm64']


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
defFitAlpha['Kud01375400Ks01364000_tsrc0_ism64_PoF0D1_TopCharge'] = 'fitr8-10'
defFitAlpha['Kud01375400Ks01364000_tsrc0_ism64_jsm64_TopCharge'] = 'fitr8-10'
defFitAlpha['Kud01375400Ks01364000_tsrc0_ism32_jsm64_TopCharge'] = 'fitr8-10'
defFitAlpha['Kud01375400Ks01364000_tsrc0_ism16_jsm64_TopCharge'] = 'fitr8-10'






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
def TflowToPhys(tflowlist):
    if isinstance(tflowlist, list) or isinstance(tflowlist,np.ndarray):
        output = []
        for iflow in tflowlist:
            if isinstance(iflow,str):
                output.append(np.sqrt(8*np.float(iflow.replace('t_f','')))*latspace)
            else:
                output.append(np.sqrt(8*iflow)*latspace)
        return output ## in fermi
    elif isinstance(tflowlist,str):
        return np.sqrt(8*np.float(tflowlist.replace('t_f','')))*latspace
    else:
        return np.sqrt(8*tflowlist)*latspace 

def WriteXml(thisfile,outputdict):
    if ShowXmlWrite and 'Parameters' not in thisfile: print('Writing to ', thisfile)
    with open(thisfile,'w') as f:
        f.write( unparse(outputdict,pretty=True).replace('\t','    '))
