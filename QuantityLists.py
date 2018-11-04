#!/usr/bin/env python

import numpy as np
from copy import copy
from collections import OrderedDict
import pandas as pa
from Params import hbarc
from MiscFuns import op_str_list
from PlotData import null_series

## TODO, remember enum is str type i think.
# NxList = [32]
# NtList = [64]

# latsplist = [0.0907]

# KudList = [0.13754,0.13700]
# KsList = [0.1364]
r0 = 0.49 ## in fm

def npsumwrap(*args):
    return np.sum(args)

def npprodrap(*args):
    return np.prod(args)

def Ident0(*args):
    if len(args) > 0: return args[0]
    else: return 0.0

def Getnpfun(thisopp):
    if thisopp == '+':
        return npsumwrap
    elif thisopp == 'x':
        return npprodrap
    else:
        return Ident0

file_formats = ['msgpack','pickle','hdf','feather','parquet','csv','json','html','None']


CombineList = ['No','+','x']
CombineList_FF = op_str_list




all_gammas = ['I','g1','g2','g3','g4','g1g2','g1g3','g1g4','g2g3','g2g4','g3g4','g1g5','g2g5','g3g5','g4g5','g5']
all_gammas += [ival+'_cmplx' for ival in all_gammas]


alpha_sum_list = ['from_src','from_sink','from_center']
alpha_sum_list += [ival+'_mult' for ival in alpha_sum_list] + [ival+'_sum' for ival in alpha_sum_list]
alpha_sum_list += [ival+'_FO2' for ival in alpha_sum_list]
# alpha_sum_list += [ival+'_sum' for ival in alpha_sum_list]
alpha_sum_list += [ival+'_sym' for ival in alpha_sum_list]
# alpha_sum_list += [ival+'_rev' for ival in alpha_sum_list]

alpha_improv_types = ['None','P5NNQsub','P5NNQsubOnlySub','Qsub','NNQsub','NNQ2sub','NNSingQ2sub','Q2sub']


rat_sum_list = ['from_src','from_sink','from_current']
rat_sum_list += [ival+'_mult' for ival in rat_sum_list] + [ival+'_sum' for ival in rat_sum_list]
rat_sum_list += [ival+'_sym' for ival in rat_sum_list]
# rat_sum_list += [ival+'_rev' for ival in rat_sum_list]



ObsList = ['TopCharge','Weinberg']
ObsListAbr = ['Top','Wein']
ObsListExt = ['Q','W','QQ','QW','WW']

FOcfglist = ['None']+ObsList


LegLocList = ['best','upper right','upper left','lower left','lower right','center right','center left']

DSList = ['Neutron','Proton','Vector','IsoVector','doub','sing','PdivideN','NdivideP','NmultiplyP','PdivideNV2']

ProjList = ['P4','P3']

BMList = ['Baryon','Meson']

# BMNumbList = ['CPEven','CPOdd','Pion','P4','g5P4','P4g5','g5']
BarNumbList = ['CPEven','CPOdd','P4','g5P4','P4g5','g5','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17']
BarNumbListCPOdd = ['CPOdd','g5P4','P4g5','g5','CPEven','P4','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17']
MesNumbList = ['Pion','I','g1','g2','g3','g4','g1g2','g1g3','g1g4','g2g3','g2g4','g3g4','g1g5','g2g5','g3g5','g4g5','g5']
# BMNumbListOdd = ['CPOdd','g5P4','P4g5','g5','CPEven']
# BMNumbListOdd = ['CPOdd','g5P4','P4g5','g5']


CurrTypes = ['Scalar','Vector','VectorTop','VectorTopBNL','VectorTopBNL_mal','VectorWein','VectorWeinBNL','PsScalar','PsVector','Tensor','GeGm']

FFTypes = ['F_{1}','F_{2}','F_{3}','F_{3}/2m_{N}','tilde{F}_{2}','tilde{F}_{3}','tilde{F}_{3}/2m_{N}','G_{E}','G_{M}','G_{A}','G_{P}','G_{S}','H_{T}','E_{T}','tilde{H}_{T}']


mpi411_ens = {}
mpi411_ens['nxyzt'] = [32,64]
mpi411_ens['a'] = 0.0907
mpi411_ens['kud'] = 0.13754
mpi411_ens['ks'] = 0.1364
mpi411_ens['Csw'] = 1715
mpi411_ens['Beta'] = 1900
mpi411_ens['tsink'] = 13
mpi411_ens['Fits'] = [['state1','fitr9-20']]
mpi411_ens['FitsAlpha'] = ['fitr10-20']
mpi411_ens['FitsAlphaFull'] = [['fitr10-20','tsumfitr7-15']]
mpi411_ens['Rat_tsum_range'] = 'tsumfitr6-32'
mpi411_ens['jsmCoeffs'] = [5.522e-07,-0.0001589143,0.9999999874]



mpi411_a_ens = copy(mpi411_ens)
mpi411_a_ens['stream_list'] = ['-a-']

mpi411_b_ens = copy(mpi411_ens)
mpi411_b_ens['stream_list'] = ['-b-']


mpi570_ens = copy(mpi411_ens)
mpi570_ens['kud'] = 0.13727
mpi570_ens['FitsAlphaFull'] = [['fitr10-20','tsumfitr10-25']]
mpi570_ens['Rat_tsum_range'] = 'tsumfitr7-32'
mpi570_ens['jsmCoeffs'] = [0.0,0.0,1.0]

mpi701_ens = copy(mpi411_ens)
mpi701_ens['kud'] = 0.13700
mpi701_ens['FitsAlphaFull'] = [['fitr10-20','tsumfitr6-15']]
mpi701_ens['Rat_tsum_range'] = 'tsumfitr4-32'
mpi701_ens['jsmCoeffs'] = [0.0,0.0,1.0]

L16_ens = {}
L16_ens['nxyzt'] = [16,32]
L16_ens['a'] = 0.1215
L16_ens['kud'] = 0.13825
L16_ens['ks'] = 0.1371
L16_ens['Csw'] = 1761
L16_ens['Beta'] = 1830
L16_ens['tsink'] = 8
L16_ens['aL'] = L16_ens['a']*L16_ens['nxyzt'][0]
L16_ens['Fits'] = [['state1','fitr9-20']]
L16_ens['FitsAlpha'] = ['fitr5-11']
L16_ens['FitsAlphaFull'] = [['fitr5-11','tsumfitr10-15']]
L16_ens['Rat_tsum_range'] = 'tsumfitr3-16'
L16_ens['jsmCoeffs'] = [0.0,0.0,1.0]

L20_ens = {}
L20_ens['nxyzt'] = [20,40]
L20_ens['a'] = 0.098
L20_ens['kud'] = 0.13700
L20_ens['ks'] = 0.1364
L20_ens['Csw'] = 1715
L20_ens['Beta'] = 1900
L20_ens['tsink'] = 10
L20_ens['aL'] = L20_ens['a']*L20_ens['nxyzt'][0]
L20_ens['Fits'] = [['state1','fitr12-19']]
L20_ens['FitsAlpha'] = ['fitr7-17']
L20_ens['FitsAlphaFull'] = [['fitr7-17','tsumfitr10-20']]
L20_ens['Rat_tsum_range'] = 'tsumfitr4-20'
L20_ens['jsmCoeffs'] = [0.0,0.0,1.0]


L20_no2_ens = copy(L20_ens)
L20_no2_ens['stream_list'] = ['-1-','-3-','-4-']

L28_ens = {}
L28_ens['nxyzt'] = [28,56]
L28_ens['a'] = 0.0685
L28_ens['kud'] = 0.1356
L28_ens['ks'] = 0.1351
L28_ens['Csw'] = 1628
L28_ens['Beta'] = 2050
L28_ens['tsink'] = 18
L28_ens['aL'] = L28_ens['a']*L28_ens['nxyzt'][0]
L28_ens['Fits'] = [['state1','fitr17-25']]
L28_ens['FitsAlpha'] = ['fitr14-21']
L28_ens['FitsAlphaFull'] = [['fitr14-21','tsumfitr10-20']]
L28_ens['Rat_tsum_range'] = 'tsumfitr10-28'
L28_ens['jsmCoeffs'] = [0.0,0.0,1.0]

qu_L24_ens = {}
qu_L24_ens['nxyzt'] = [24,48]
qu_L24_ens['a'] = 0.0685
# qu_L24_ens['kud'] = 0.1356
# qu_L24_ens['ks'] = 0.1351
# qu_L24_ens['Csw'] = 1628
qu_L24_ens['Beta'] = 6.00
# qu_L24_ens['tsink'] = 18
# qu_L24_ens['aL'] = qu_L24_ens['a']*qu_L24_ens['nxyzt'][0]
# qu_L24_ens['Fits'] = [['state1','fitr17-25']]
# qu_L24_ens['FitsAlpha'] = ['fitr14-21']
# qu_L24_ens['FitsAlphaFull'] = [['fitr14-21','tsumfitr10-20']]
# qu_L24_ens['Rat_tsum_range'] = 'tsumfitr10-28'
# qu_L24_ens['jsmCoeffs'] = [0.0,0.0,1.0]


ens_dict = OrderedDict()
ens_dict['mpi411'] = mpi411_ens
ens_dict['mpi570'] = mpi570_ens
ens_dict['mpi701'] = mpi701_ens
ens_dict['mpi411-a-'] = mpi411_a_ens
ens_dict['mpi411-b-'] = mpi411_b_ens
ens_dict['L16'] = L16_ens
ens_dict['L20'] = L20_ens
ens_dict['L28'] = L28_ens
ens_dict['L20_no2'] = L20_no2_ens
ens_dict['L20_no2'] = L20_no2_ens
ens_dict['qu_L24'] = qu_L24_ens

setup_dict = OrderedDict()
setup_dict['pion_mass'] = [mpi411_ens,mpi570_ens,mpi701_ens]
setup_dict['lattice_spacing'] = [L16_ens,L20_ens,L28_ens]
setup_dict['box_size'] = [mpi701_ens,L20_ens]
setup_dict['pion_mass_split_stream'] = [mpi411_a_ens,mpi411_b_ens,mpi570_ens,mpi701_ens]
setup_dict['pion_mass_mpi411a'] = [mpi411_a_ens,mpi570_ens,mpi701_ens]
setup_dict['pion_mass_mpi411b'] = [mpi411_b_ens,mpi570_ens,mpi701_ens]
setup_dict['lattice_spacing_no2'] = [L16_ens,L20_no2_ens,L28_ens]
setup_dict['All_split'] = setup_dict['pion_mass_split_stream'] + setup_dict['lattice_spacing_no2']
setup_dict['All'] = setup_dict['pion_mass'] + setup_dict['lattice_spacing']

############### section of data from papers ##############################

# topological susceptibility taken from https://arxiv.org/pdf/1705.10906.pdf
# mobius domain wall fermions
# Nf = 2 + 1
# a > 0.119 fm

top_susc_1705_10906 = pa.DataFrame()
top_susc_1705_10906_label = r'S. Aoki et al, 2017, DWF, $N_f =2 + 1$, $a\in [0.043,0.081]fm$'
top_susc_1705_10906.name = top_susc_1705_10906_label
top_susc_1705_10906.loc[:,'$a[fm]$'] = pa.Series([hbarc/x for x in [2.453]*8 + [3.610]*6 + [4.496]])

top_susc_1705_10906.loc[:,'$L/a$'] = pa.Series([32]*7+[48]*7 + [64])
top_susc_1705_10906.loc[:,'$T/a$'] = pa.Series([64]*7+[96]*7 + [128])


top_susc_1705_10906.loc[:,'$m_{\pi}[MeV]$'] = pa.Series([
230,310,310,400,400,500,500,230,300,300,410,410,500,500,280])

top_susc_1705_10906.loc[:,'$\chi_{t}a^{4}$'] = pa.Series([
0.217*10**-5,0.400*10**-5,1.010*10**-5,1.593*10**-5,0.558*10**-5,
0.560*10**-5,1.201*10**-5,0.282*10**-5,0.91*10**-6,2.18*10**-6,
1.21*10**-6,0.59*10**-6,0.55*10**-6,1.70*10**-6])


top_susc_1705_10906.loc[:,'$\Delta \chi_{t}a^{4}$'] = pa.Series([
(0.039 + 0.014) *10**-5,(0.048 + 0.021) *10**-5,(0.192 + 0.13)  *10**-5,
(0.228 + 0.91)  *10**-5,(0.090 + 0.053) *10**-5,(0.112 + 0.22)  *10**-5,
(0.146 + 0.23)  *10**-5,(0.021 + 0.042) *10**-5,(0.23 + 0.03)   *10**-6,
(0.54 + 0.09)   *10**-6,(0.19 + 0.07)   *10**-6,(0.22 + 0.12)   *10**-6,
(0.30 + 0.13)   *10**-6,(0.31 + 0.21)   *10**-6,(0.118 + 0.002) *10**-6])


top_susc_1705_10906.loc[:,'$\chi_{t}^{1/4}[GeV]$'] = \
    top_susc_1705_10906.apply(lambda x : (x['$\chi_{t}a^{4}$']**0.25)*hbarc/x['$a[fm]$'],axis=1)
top_susc_1705_10906.loc[:,'$\Delta \chi_{t}^{1/4}[GeV]$'] = \
    top_susc_1705_10906.apply(lambda x : x['$\Delta \chi_{t}a^{4}$']*hbarc/((x['$\chi_{t}a^{4}$']**0.75)*4*x['$a[fm]$']),axis=1)

ylab = '$\chi_{t}^{1/4}[GeV]$'
yerrlab = '$\Delta \chi_{t}^{1/4}[GeV]$'
xlab = '$m_{\pi}[MeV]$'

top_susc_1705_10906_plot = pa.Series()
top_susc_1705_10906_plot['x_data'] = top_susc_1705_10906[xlab]
top_susc_1705_10906_plot['y_data'] = top_susc_1705_10906[ylab]
top_susc_1705_10906_plot['yerr_data'] = top_susc_1705_10906[yerrlab]
top_susc_1705_10906_plot['fit_class'] = None
top_susc_1705_10906_plot['shift'] = None
top_susc_1705_10906_plot['label'] = top_susc_1705_10906_label
top_susc_1705_10906_plot['color'] = None
top_susc_1705_10906_plot['symbol'] = None
top_susc_1705_10906_plot['type'] = 'error_bar'

#
# import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
# import matplotlib.pyplot as pl
# %matplotlib inline
# pl.errorbar(top_susc_1705_10906[xlab],top_susc_1705_10906[ylab],top_susc_1705_10906[yerrlab],label=top_susc_1705_10906_label)
# pl.ylabel(ylab)
# pl.xlabel(xlab)
# pl.legend()
# pl.show()


# topological susceptibility taken from https://arxiv.org/pdf/0710.1130.pdf
# overlap dirac operato
# Nf=2
# a > 0.119 fm

top_susc_0710_1130 = pa.DataFrame()
top_susc_0710_1130_label = r'S. Aoki et al, 2008, Ov-D, $16^{3} \times 32$, $N_f =2$, $a\in [0.119,0.128]fm$'
top_susc_0710_1130.loc[:,'$a[fm]$'] = pa.Series([0.1194,0.1206,0.1215,0.1236,0.1251,0.1272])
top_susc_0710_1130.loc[:,'$m_{\pi}r_{0}$'] = pa.Series([0.7096,0.8935,1.053,1.238,1.456,1.734])
top_susc_0710_1130.loc[:,'$\Delta m_{\pi}r_{0}$'] = pa.Series([0.0102,0.0137,0.013,0.014,0.015,0.017])
top_susc_0710_1130.loc[:,'$\chi_{t}r_{0}^{4}$'] = pa.Series([0.00559,0.00926,0.0123,0.0187,0.0263,0.0199])
top_susc_0710_1130.loc[:,'$\Delta \chi_{t}r_{0}^{4}$'] = pa.Series([0.00047,0.00129,0.0022,0.0034,0.0033,0.0033])


top_susc_0710_1130.loc[:,'$m_{\pi}[MeV]$'] = top_susc_0710_1130.apply(lambda x : x['$m_{\pi}r_{0}$']*hbarc*1000/r0,axis=1)
top_susc_0710_1130.loc[:,'$\Delta m_{\pi}[MeV]$'] = top_susc_0710_1130.apply(lambda x : x['$\Delta m_{\pi}r_{0}$']*hbarc*1000/r0,axis=1)
top_susc_0710_1130.loc[:,'$\chi_{t}[GeV]^{4}$'] = top_susc_0710_1130.apply(lambda x : x['$\chi_{t}r_{0}^{4}$']*(hbarc/r0)**4,axis=1)
top_susc_0710_1130.loc[:,'$\Delta \chi_{t}[GeV]^{4}$'] = top_susc_0710_1130.apply(lambda x : x['$\Delta \chi_{t}r_{0}^{4}$']*(hbarc/r0)**4,axis=1)


top_susc_0710_1130.loc[:,'$\chi_{t}^{1/4}[GeV]$'] = top_susc_0710_1130['$\chi_{t}[GeV]^{4}$'].apply(lambda x : x**(0.25))
top_susc_0710_1130.loc[:,'$\Delta \chi_{t}^{1/4}[GeV]$'] = top_susc_0710_1130.apply(lambda x : x['$\Delta \chi_{t}[GeV]^{4}$']/(4*x['$\chi_{t}[GeV]^{4}$']**(0.75)),axis=1)
ylab = '$\chi_{t}^{1/4}[GeV]$'
yerrlab = '$\Delta \chi_{t}^{1/4}[GeV]$'
xlab = '$m_{\pi}[MeV]$'

top_susc_0710_1130_plot = pa.Series()
top_susc_0710_1130_plot['x_data'] = top_susc_0710_1130[xlab]
top_susc_0710_1130_plot['y_data'] = top_susc_0710_1130[ylab]
top_susc_0710_1130_plot['yerr_data'] = top_susc_0710_1130[yerrlab]
top_susc_0710_1130_plot['fit_class'] = None
top_susc_0710_1130_plot['shift'] = None
top_susc_0710_1130_plot['label'] = top_susc_0710_1130_label
top_susc_0710_1130_plot['color'] = None
top_susc_0710_1130_plot['symbol'] = None
top_susc_0710_1130_plot['type'] = 'error_bar'


# import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
# import matplotlib.pyplot as pl
# pl.errorbar(top_susc_0710_1130[xlab],top_susc_0710_1130[ylab],top_susc_0710_1130[yerrlab],label=top_susc_0710_1130_label)
# pl.ylabel(ylab)
# pl.xlabel(xlab)


# topological susceptibility taken from https://arxiv.org/pdf/1406.5363.pdf
# alpha colab Mattia Bruno, Stefan Schaefer and Rainer Sommer
# Improved Wilson
# Nf=2
# \beta = 5.2 for first 7
# \beta = 5.3 for next 4
# \beta = 5.5 for last 2
tf1 =  0.061 # fm^2
top_susc_1406_5363 = pa.DataFrame()
top_susc_1406_5363_label = 'M. Bruno et al, 2014, Imp Wil, $N_f =2$'
top_susc_1406_5363.loc[:,'$\beta$'] = pa.Series([5.2, 5.2, 5.2, 5.2, 5.3, 5.3, 5.3, 5.3, 5.5, 5.5])
top_susc_1406_5363.loc[:,'$T/a$'] = pa.Series([64 ,64, 64 ,96 ,64 ,96 ,96 ,128,96 ,128])
top_susc_1406_5363.loc[:,'$L/a$'] = pa.Series([32, 32, 32, 48, 32, 48, 48, 64,  48, 64])
top_susc_1406_5363.loc[:,'$\kappa$'] = pa.Series([  0.13565, 0.13590, 0.13594, 0.13597,
                                                    0.13625, 0.13635, 0.13638, 0.13642,
                                                    0.13667, 0.13671])
top_susc_1406_5363.loc[:,'$m_{\pi}[MeV]$'] = pa.Series([    630, 380, 330, 280,
                                                            430, 310, 260, 190, 340, 260 ])



top_susc_1406_5363.loc[:,'$t_{1}^{2}\chi_{t}$'] = pa.Series([
1.48/1000, 1.06/1000, 0.88/1000, 0.81/1000, 1.02/1000, 1.00/1000, 0.72/1000, 0.58/1000, 0.67/1000, 0.40/1000])
top_susc_1406_5363.loc[:,'$\Delta t_{1}^{2}\chi_{t}$'] = pa.Series([
0.14/1000, 0.10/1000, 0.07/1000, 0.14/1000, 0.09/1000, 0.16/1000, 0.09/1000, 0.07/1000, 0.14/1000, 0.08/1000])


top_susc_1406_5363.loc[:,'$\chi_{t}^{1/4}[GeV]$'] = \
    top_susc_1406_5363['$t_{1}^{2}\chi_{t}$'].apply(lambda x : hbarc*((x/(tf1**2))**0.25))
top_susc_1406_5363.loc[:,'$\Delta \chi_{t}^{1/4}[GeV]$'] = \
    top_susc_1406_5363.apply(lambda x : x['$\Delta t_{1}^{2}\chi_{t}$']/((4*(tf1**0.5)*x['$t_{1}^{2}\chi_{t}$']**(0.75))),axis=1)
ylab = '$\chi_{t}^{1/4}[GeV]$'
yerrlab = '$\Delta \chi_{t}^{1/4}[GeV]$'
xlab = '$m_{\pi}[MeV]$'

top_susc_1406_5363_plot = pa.Series()
top_susc_1406_5363_plot['x_data'] = top_susc_1406_5363[xlab]
top_susc_1406_5363_plot['y_data'] = top_susc_1406_5363[ylab]
top_susc_1406_5363_plot['yerr_data'] = top_susc_1406_5363[yerrlab]
top_susc_1406_5363_plot['fit_class'] = None
top_susc_1406_5363_plot['shift'] = None
top_susc_1406_5363_plot['label'] = top_susc_1406_5363_label
top_susc_1406_5363_plot['color'] = None
top_susc_1406_5363_plot['symbol'] = None
top_susc_1406_5363_plot['type'] = 'error_bar'


# import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
# import matplotlib.pyplot as pl
# %matplotlib inline
# pl.errorbar(top_susc_1406_5363[xlab],top_susc_1406_5363[ylab],top_susc_1406_5363[yerrlab],label=top_susc_1406_5363_label)
# pl.ylabel(ylab)
# pl.xlabel(xlab)
# pl.legend()
# pl.show()




# topological susceptibility taken from https://arxiv.org/pdf/1312.5161.pdf
# twisted mass colab Krzysztof Cichy Elena Garcia-Ramos Karl Jansena
## N_f = 2 run
top_susc_1312_5161 = pa.DataFrame()
top_susc_1312_5161_label = 'K. Cichy et al, 2014, Twist, $N_f =2$, $a\in [0.083,0.0875]fm$'


top_susc_1312_5161.loc[:,'$L/a$'] = pa.Series([16,20,24,32,24,24])
top_susc_1312_5161.loc[:,'$T/a$'] = pa.Series([32,40,48,64,48,48])
top_susc_1312_5161.loc[:,'$L[fm]$'] = pa.Series([1.4,1.7,2.0,2.7,2.0,2.0])
top_susc_1312_5161.loc[:,'$m_{\pi}L$'] = pa.Series([2.5,2.8,3.3,4.3,4.1,4.7])
top_susc_1312_5161.loc[:,'$m_{\pi}[MeV]$'] = \
    top_susc_1312_5161.apply(lambda x : hbarc*1000*x['$m_{\pi}L$']/x['$L[fm]$'],axis=1)




top_susc_1312_5161.loc[:,'$r_{0}^{4}\chi_{t}$'] = pa.Series([0.0097,0.0092,0.0096,0.0082,0.0106,0.0118])
top_susc_1312_5161.loc[:,'$\Delta r_{0}^{4}\chi_{t}$'] = pa.Series([0.0020,0.0015,0.0015,0.0017,0.0021,0.0029])


top_susc_1312_5161.loc[:,'$\chi_{t}^{1/4}[GeV]$'] = \
    top_susc_1312_5161['$r_{0}^{4}\chi_{t}$'].apply(lambda x : (x**0.25)*hbarc/r0)
top_susc_1312_5161.loc[:,'$\Delta \chi_{t}^{1/4}[GeV]$'] = \
    top_susc_1312_5161.apply(lambda x : x['$\Delta r_{0}^{4}\chi_{t}$']*(hbarc/r0)/((x['$r_{0}^{4}\chi_{t}$']**0.75)*4),axis=1)
ylab = '$\chi_{t}^{1/4}[GeV]$'
yerrlab = '$\Delta \chi_{t}^{1/4}[GeV]$'
xlab = '$m_{\pi}[MeV]$'

top_susc_1312_5161_plot = pa.Series()
top_susc_1312_5161_plot['x_data'] = top_susc_1312_5161[xlab]
top_susc_1312_5161_plot['y_data'] = top_susc_1312_5161[ylab]
top_susc_1312_5161_plot['yerr_data'] = top_susc_1312_5161[yerrlab]
top_susc_1312_5161_plot['fit_class'] = None
top_susc_1312_5161_plot['shift'] = None
top_susc_1312_5161_plot['label'] = top_susc_1312_5161_label
top_susc_1312_5161_plot['color'] = None
top_susc_1312_5161_plot['symbol'] = None
top_susc_1312_5161_plot['type'] = 'error_bar'


# import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
# import matplotlib.pyplot as pl
# %matplotlib inline
# pl.errorbar(top_susc_1312_5161[xlab],top_susc_1312_5161[ylab],top_susc_1312_5161[yerrlab],label=top_susc_1312_5161_label)
# pl.ylabel(ylab)
# pl.xlabel(xlab)
# pl.legend()
# pl.show()


# topological susceptibility taken from https://arxiv.org/pdf/1312.5161.pdf
# twisted mass colab Krzysztof Cichy Elena Garcia-Ramos Karl Jansena
## N_f = 2 + 1 + 1 run
top_susc_1312_5161_211 = pa.DataFrame()
top_susc_1312_5161_211_label = 'K. Cichy et al, 2014, Twist, $N_f =2 + 1 + 1$, $a\in [0.059,0.0875]fm$'


top_susc_1312_5161_211.loc[:,'$\beta$'] = pa.Series([
1.90,1.90,1.90,1.90,1.90,1.90,1.90,1.95,1.95,1.95,1.95,1.95,2.10,2.10,2.10])

top_susc_1312_5161_211.loc[:,'$L/a$'] = pa.Series([
32,20,24,32,32,24,24,32,32,32,32,24,48,48,32])

top_susc_1312_5161_211.loc[:,'$T/a$'] = pa.Series([
64,40,48,64,64,48,48,64,64,64,64,48,96,96,64])

top_susc_1312_5161_211.loc[:,'$L[fm]$'] = pa.Series([
2.8,1.7,2.1,2.8,2.8,2.1,2.1,2.5,2.5,2.5,2.5,1.9,2.9,2.9,1.9])

top_susc_1312_5161_211.loc[:,'$m_{\pi}L$'] = pa.Series([
4.0,3.0,3.5,4.5,5.1,4.2,4.8,3.4,4.0,5.0,5.8,4.7,3.9,4.7,3.9])

top_susc_1312_5161_211.loc[:,'$m_{\pi}[MeV]$'] = \
    top_susc_1312_5161_211.apply(lambda x : hbarc*1000*x['$m_{\pi}L$']/x['$L[fm]$'],axis=1)

top_susc_1312_5161_211.loc[:,'$r_{0}^{4}\chi_{t}$'] = pa.Series([
0.0072,0.0130,0.0086,0.0074,0.0081,0.0092,0.0114,0.0070,0.0067,0.0080,0.0090,0.0106,0.0041,0.0073,0.0125])

top_susc_1312_5161_211.loc[:,'$\Delta r_{0}^{4}\chi_{t}$'] = pa.Series([
0.0015,0.0031,0.0018,0.0015,0.0017,0.0020,0.0024,0.0014,0.0012,0.0011,0.0016,0.0019,0.0013,0.0027,0.0027])

top_susc_1312_5161_211.loc[:,'$\chi_{t}^{1/4}[GeV]$'] = \
    top_susc_1312_5161_211['$r_{0}^{4}\chi_{t}$'].apply(lambda x : (x**0.25)*hbarc/r0)
top_susc_1312_5161_211.loc[:,'$\Delta \chi_{t}^{1/4}[GeV]$'] = \
    top_susc_1312_5161_211.apply(lambda x : x['$\Delta r_{0}^{4}\chi_{t}$']*hbarc/((x['$r_{0}^{4}\chi_{t}$']**0.75)*4*r0),axis=1)
ylab = '$\chi_{t}^{1/4}[GeV]$'
yerrlab = '$\Delta \chi_{t}^{1/4}[GeV]$'
xlab = '$m_{\pi}[MeV]$'

top_susc_1312_5161_211_plot = pa.Series()
top_susc_1312_5161_211_plot['x_data'] = top_susc_1312_5161_211[xlab]
top_susc_1312_5161_211_plot['y_data'] = top_susc_1312_5161_211[ylab]
top_susc_1312_5161_211_plot['yerr_data'] = top_susc_1312_5161_211[yerrlab]
top_susc_1312_5161_211_plot['fit_class'] = None
top_susc_1312_5161_211_plot['shift'] = None
top_susc_1312_5161_211_plot['label'] = top_susc_1312_5161_211_label
top_susc_1312_5161_211_plot['color'] = None
top_susc_1312_5161_211_plot['symbol'] = None
top_susc_1312_5161_211_plot['type'] = 'error_bar'


# import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
# import matplotlib.pyplot as pl
# %matplotlib inline
# pl.errorbar(top_susc_1312_5161_211[xlab],top_susc_1312_5161_211[ylab],top_susc_1312_5161_211[yerrlab],label=top_susc_1312_5161_211_label)
# pl.ylabel(ylab)
# pl.xlabel(xlab)
# pl.legend()
# pl.show()


# topological susceptibility taken from https://arxiv.org/pdf/1512.06746.pdf
# Claudio Bonati, Massimo D'Elia, Marco Mariti, Guido Martinelli, Michele
# Mesiti, Francesco Negro, Francesco Sanfilippo, Giovanni Villadorof
## N_f = 2 + 1 + 1 run
## the dataframe is usless, unless you wanted to plot lattice space dependance.
top_susc_1512_06746 = pa.DataFrame()
top_susc_1512_06746_label = 'C. Bonati et al, 2016, Stag, $N_{f}=2+1$, $a\in [0.0572,0.1249]fm$'

top_susc_1512_06746.loc[:,'$a[fm]$'] = pa.Series([0.1249,0.0989,0.0824,0.0707,0.0572])

top_susc_1512_06746.loc[:,'$m_{\pi}[MeV]$'] = pa.Series([135]*len(top_susc_1512_06746.loc[:,'$a[fm]$']))

top_susc_1512_06746.loc[:,'$\chi_{t}a^{4}$'] = pa.Series([
8.55*10**-5,2.22*10**-5,7.50*10**-6,2.60*10**-6,6.60*10**-7])

top_susc_1512_06746.loc[:,'$\Delta \chi_{t}a^{4}$'] = pa.Series([
0.32*10**-5,0.10*10**-5,0.40*10**-6,0.36*10**-6,0.82*10**-7])


top_susc_1512_06746.loc[:,'$\chi_{t}^{1/4}[GeV]$'] = \
    top_susc_1512_06746.apply(lambda x : (x['$\chi_{t}a^{4}$']**0.25)*hbarc/x['$a[fm]$'],axis=1)
top_susc_1512_06746.loc[:,'$\Delta \chi_{t}^{1/4}[GeV]$'] = \
    top_susc_1512_06746.apply(lambda x : x['$\Delta \chi_{t}a^{4}$']*hbarc/((x['$\chi_{t}a^{4}$']**0.75)*4*x['$a[fm]$']),axis=1)
ylab = '$\chi_{t}^{1/4}[GeV]$'
yerrlab = '$\Delta \chi_{t}^{1/4}[GeV]$'
xlab = '$m_{\pi}[MeV]$'
# xlab = '$a[fm]$'

val_1512_06746 = 0.069 ## GeV using a extrapolation O(a^{4})
err_1512_06746 = 0.009

# val_1512_06746 = 0.073 ## GeV using a extrapolation O(a^{2})
# val_1512_06746 = 0.009
#
# val_1512_06746 = top_susc_1512_06746[ylab].iloc[-1] ## GeV using smallest lattice
# val_1512_06746 = top_susc_1512_06746[yerrlab].iloc[-1]

top_susc_1512_06746_plot = pa.Series()
top_susc_1512_06746_plot['x_data'] = [top_susc_1512_06746[xlab][0]]
top_susc_1512_06746_plot['y_data'] = [val_1512_06746]
top_susc_1512_06746_plot['yerr_data'] = [err_1512_06746]
top_susc_1512_06746_plot['fit_class'] = None
top_susc_1512_06746_plot['shift'] = None
top_susc_1512_06746_plot['label'] = top_susc_1512_06746_label
top_susc_1512_06746_plot['color'] = None
top_susc_1512_06746_plot['symbol'] = None
top_susc_1512_06746_plot['type'] = 'error_bar'
#
#
# import matplotlib
# matplotlib.use('Agg') # Must be before importing matplotlib.pyplot or pylab!
# import matplotlib.pyplot as pl
# %matplotlib inline
# pl.errorbar([top_susc_1512_06746[xlab][0]],[val_1512_06746],[err_1512_06746],label=top_susc_1512_06746_label)
# pl.ylabel(ylab)
# pl.xlabel(xlab)
# pl.legend()
# pl.show()



############### Form Factor Others Results ###############################

# FF_index_names = [r'$nxyz$',r'$nt$',r'$a[fm]$',r'$a^{-1}[GeV]$',r'$Volume[fm]^3$',
#                r'$m_{\pi}[GeV]$',r'$m_{N} [GeV]$',r'$alpha$']




FF_1512_00566 = pa.DataFrame()
FF_1512_00566.name = r'E. Shintani et al, 2015, $N_{f}=2+1$, DWF'
FF_1512_00566.loc[:,r'$nxyz$'] = pa.Series([24,24,32])
FF_1512_00566.loc[:,r'$nt$'] = pa.Series([24*2,24*2,32*2])
FF_1512_00566.loc[:,r'$a[fm]$'] = pa.Series([0.1125,0.1125,0.14375])
FF_1512_00566.loc[:,r'$a^{-1}[GeV]$'] = hbarc/FF_1512_00566.loc[:,r'$a[fm]$']
FF_1512_00566.loc[:,r'$Volume[fm]^3$'] = (FF_1512_00566.loc[:,r'$a[fm]$'] * FF_1512_00566.loc[:,r'$nxyz$'])**3


FF_1512_00566.loc[:,r'$m_{\pi}[GeV]$'] = pa.Series([0.42,0.33,0.17])
FF_1512_00566.loc[:,r'$m_{\pi}[MeV]$'] = FF_1512_00566.loc[:,r'$m_{\pi}[GeV]$']*1000

FF_1512_00566.loc[:,r'$m_{N} fitr$'] = pa.Series(['fitr6-12',r'$fitr7-13',r'$fitr5-10$'])
FF_1512_00566.loc[:,r'$m_{N} [GeV]$'] = pa.Series([1.2641,1.1738,0.9746])
FF_1512_00566.loc[:,r'$\Delta m_{N} [GeV]$'] = pa.Series([0.0028,0.0025,0.0066])

FF_1512_00566.loc[:,r'$alpha fitr$'] = pa.Series(['fitr5-9',r'$fitr5-9',r'$fitr5-9$'])
FF_1512_00566.loc[:,r'$alpha$'] = pa.Series([-0.370,-0.356,-0.333])
FF_1512_00566.loc[:,r'$\Delta alpha$'] = pa.Series([0.022,0.022,0.128])

FF_1512_00566.loc[:,r'$d_{p}[efm] t_{sink}=1.32fm$'] = pa.Series([0.064,0.030,0.101])
FF_1512_00566.loc[:,r'$\Delta d_{p}[efm] t_{sink}=1.32fm$'] = pa.Series([0.027,0.025,0.090])
FF_1512_00566.loc[:,r'$d_{p}[efm] t_{sink}=0.9fm$'] = pa.Series([0.035,0.015,float('NaN')])
FF_1512_00566.loc[:,r'$\Delta d_{p}[efm] t_{sink}=0.9fm$'] = pa.Series([0.019,0.012,float('NaN')])
FF_1512_00566.loc[:,r'$d_{n}[efm] t_{sink}=1.32fm$'] = pa.Series([-0.021,-0.053,-0.093])
FF_1512_00566.loc[:,r'$\Delta d_{n}[efm] t_{sink}=1.32fm$'] = pa.Series([0.015,0.018,0.043])
FF_1512_00566.loc[:,r'$d_{n}[efm] t_{sink}=0.9fm$'] = pa.Series([-0.016,-0.029,float('NaN')])
FF_1512_00566.loc[:,r'$\Delta d_{n}[efm] t_{sink}=0.9fm$'] = pa.Series([0.011,0.008,float('NaN')])


xlab = r'$m_{\pi}[MeV]$'
xplot = FF_1512_00566.loc[:,xlab]


FF_1512_00566_plot = pa.DataFrame()
this_lab = FF_1512_00566.name+r' Neutron $t_{sink}=1.32fm$'
FF_1512_00566_plot.loc[:,this_lab] =  null_series
FF_1512_00566_plot.at['x_data',this_lab] = xplot
FF_1512_00566_plot.at['y_data',this_lab] = FF_1512_00566.loc[:,r'$d_{n}[efm] t_{sink}=1.32fm$']
FF_1512_00566_plot.at['yerr_data',this_lab] = FF_1512_00566.loc[:,r'$\Delta d_{n}[efm] t_{sink}=1.32fm$']
FF_1512_00566_plot.at['label',this_lab] = this_lab
FF_1512_00566_plot.at['type',this_lab] = 'error_bar'

this_lab = FF_1512_00566.name+r' Neutron $t_{sink}=0.9fm$'
FF_1512_00566_plot.at[:,this_lab] =  null_series
FF_1512_00566_plot.at['x_data',this_lab] = xplot
FF_1512_00566_plot.at['y_data',this_lab] = FF_1512_00566.loc[:,r'$d_{n}[efm] t_{sink}=0.9fm$']
FF_1512_00566_plot.at['yerr_data',this_lab] = FF_1512_00566.loc[:,r'$\Delta d_{n}[efm] t_{sink}=0.9fm$']
FF_1512_00566_plot.at['label',this_lab] = this_lab
FF_1512_00566_plot.at['type',this_lab] = 'error_bar'

this_lab = FF_1512_00566.name+r' Proton $t_{sink}=1.32fm$'
FF_1512_00566_plot.at[:,this_lab] =  null_series
FF_1512_00566_plot.at['x_data',this_lab] = xplot
FF_1512_00566_plot.at['y_data',this_lab] = FF_1512_00566.loc[:,r'$d_{p}[efm] t_{sink}=1.32fm$']
FF_1512_00566_plot.at['yerr_data',this_lab] = FF_1512_00566.loc[:,r'$\Delta d_{p}[efm] t_{sink}=1.32fm$']
FF_1512_00566_plot.at['label',this_lab] = this_lab
FF_1512_00566_plot.at['type',this_lab] = 'error_bar'

this_lab = FF_1512_00566.name+r' Proton $t_{sink}=0.9fm$'
FF_1512_00566_plot.at[:,this_lab] =  null_series
FF_1512_00566_plot.at['x_data',this_lab] = xplot
FF_1512_00566_plot.at['y_data',this_lab] = FF_1512_00566.loc[:,r'$d_{p}[efm] t_{sink}=0.9fm$']
FF_1512_00566_plot.at['yerr_data',this_lab] = FF_1512_00566.loc[:,r'$\Delta d_{p}[efm] t_{sink}=0.9fm$']
FF_1512_00566_plot.at['label',this_lab] = this_lab
FF_1512_00566_plot.at['type',this_lab] = 'error_bar'


########################################################################
FF_1510_05823 = pa.DataFrame()
FF_1510_05823.name = r'C. Alexandrou et al, 2016, $N_{f}=2+1+1$ TMF'
FF_1510_05823.loc[:,r'Action'] = pa.Series(['Wilson','Symanzik three-level','Iwasaki','Wilson','Symanzik three-level','Iwasaki','1701.07792'])
FF_1510_05823.loc[:,r'Q_Smearing'] = pa.Series(['Cooling','Cooling','Cooling','Gradient Flow','Gradient Flow','Gradient Flow','1701.07792'])
FF_1510_05823.loc[:,r'$nxyz$'] = pa.Series([32,32,32,32,32,32,32])
FF_1510_05823.loc[:,r'$nt$'] = pa.Series([32*2,32*2,32*2,32*2,32*2,32*2,32*2])
FF_1510_05823.loc[:,r'$a[fm]$'] = pa.Series([0.0823,0.0823,0.0823,0.0823,0.0823,0.0823,0.0823])
FF_1510_05823.loc[:,r'$\Delta a[fm]$'] = pa.Series([0.0010,0.0010,0.0010,0.0010,0.0010,0.0010,0.0010])
FF_1510_05823.loc[:,r'$a^{-1}[GeV]$'] = hbarc/FF_1510_05823.loc[:,r'$a[fm]$']
FF_1510_05823.loc[:,r'$Volume[fm]^3$'] = (FF_1510_05823.loc[:,r'$a[fm]$'] * FF_1510_05823.loc[:,r'$nxyz$'])**3


FF_1510_05823.loc[:,r'$m_{\pi}[GeV]$'] = pa.Series([0.373,0.373,0.373,0.373,0.373,0.373,0.373])
FF_1510_05823.loc[:,r'$m_{\pi}[MeV]$'] = FF_1510_05823.loc[:,r'$m_{\pi}[GeV]$']*1000

FF_1510_05823.loc[:,r'$m_{N} [GeV]$'] = pa.Series([1.220,1.220,1.220,1.220,1.220,1.220,1.220])
FF_1510_05823.loc[:,r'$\Delta m_{N} [GeV]$'] = pa.Series([0.005,0.005,0.005,0.005,0.005,0.005,0.005])

## selected using Iwasaki action,
FF_1510_05823.loc[:,r'$alpha fitr$'] = pa.Series(['fitr9-20','fitr9-20','fitr9-20','fitr9-20','fitr9-20','fitr9-20','fitr9-20'])
FF_1510_05823.loc[:,r'$alpha$'] = pa.Series([-0.217,-0.217,-0.217,-0.217,-0.217,-0.217,-0.217])
FF_1510_05823.loc[:,r'$\Delta alpha$'] = pa.Series([0.018,0.018,0.018,0.018,0.018,0.018,0.018])

paper_value = -0.555 / (2*1.220)
paper_value_err = 0.005/1.220
paper_value_err += 0.074/0.555
paper_value_err = paper_value_err * paper_value

paper_value_rot = 0.094 / (2*1.220)
paper_value_rot_err = paper_value_err
# paper_value_rot_err = 0.005/1.220
# paper_value_rot_err += 0.074/0.094
# paper_value_rot_err = paper_value_rot_err * paper_value_rot


## Just hacked it so that smaller tsink is F_{3} and larger tsink is \tilde{F}_{3}
FF_1510_05823.loc[:,r'$d_{n}[efm] t_{sink}=0.9876fm$'] = pa.Series([-0.035,-0.036,-0.035,-0.039,-0.046,-0.041,paper_value])
FF_1510_05823.loc[:,r'$\Delta d_{n}[efm] t_{sink}=0.9876fm$'] = pa.Series([0.009,0.010,0.010,0.010,0.012,0.012,paper_value_err])
FF_1510_05823.loc[:,r'$d_{n}[efm] t_{sink}=1.1522fm$'] = pa.Series([-0.056,-0.072,-0.049,-0.065,-0.042,-0.076,paper_value_rot])
FF_1510_05823.loc[:,r'$\Delta d_{n}[efm] t_{sink}=1.1522fm$'] = pa.Series([0.045,0.050,0.022,0.042,0.040,0.023,paper_value_rot_err])

### Derivative to the ratio technique (DRM)
FF_1510_05823.loc[:,r'$d_{n}[efm] t_{sink}=0.9876fm$ DRM'] = pa.Series([-0.044,-0.043,-0.044,-0.046,-0.043,-0.042,float('NaN')])
FF_1510_05823.loc[:,r'$\Delta d_{n}[efm] t_{sink}=0.9876fm$ DRM'] = pa.Series([0.012,0.009,0.011,0.015,0.010,0.010,float('NaN')])
FF_1510_05823.loc[:,r'$d_{n}[efm] t_{sink}=1.1522fm$ DRM'] = pa.Series([-0.044,-0.043,-0.043,-0.046,-0.042,-0.040,float('NaN')])
FF_1510_05823.loc[:,r'$\Delta d_{n}[efm] t_{sink}=1.1522fm$ DRM'] = pa.Series([0.016,0.015,0.021,0.023,0.019,0.019,float('NaN')])


### elimination of the momentum in the plateau region technique (EMPR)
FF_1510_05823.loc[:,r'$d_{n}[efm] t_{sink}=0.9876fm$ EMPR'] = pa.Series([-0.052,-0.066,-0.082,-0.056,-0.068,-0.082,float('NaN')])
FF_1510_05823.loc[:,r'$\Delta d_{n}[efm] t_{sink}=0.9876fm$ EMPR'] = pa.Series([0.017,0.018,0.021,0.017,0.017,0.021,float('NaN')])
FF_1510_05823.loc[:,r'$d_{n}[efm] t_{sink}=1.1522fm$ EMPR'] = pa.Series([-0.043,-0.053,-0.076,-0.043,-0.055,-0.073,float('NaN')])
FF_1510_05823.loc[:,r'$\Delta d_{n}[efm] t_{sink}=1.1522fm$ EMPR'] = pa.Series([0.026,0.028,0.032,0.026,0.028,0.032,float('NaN')])


# xlab = r'$m_{\pi}[GeV]$'
# xplot = FF_1510_05823.loc[:,xlab]
first_key = ['Wilson','Gradient Flow',slice(None)]
FF_1510_05823_MuIn = FF_1510_05823.set_index([r'Action',r'Q_Smearing',r'$m_{\pi}[MeV]$'])


FF_1510_05823_plot = pa.DataFrame()
this_lab = FF_1510_05823.name+r' Neutron $t_{sink}=0.9876fm$'
FF_1510_05823_plot.loc[:,this_lab] =  null_series
FF_1510_05823_plot.at['x_data',this_lab] = 'from_keys'
FF_1510_05823_plot.at['y_data',this_lab] = FF_1510_05823_MuIn.loc[:,r'$d_{n}[efm] t_{sink}=0.9876fm$']
FF_1510_05823_plot.at['yerr_data',this_lab] = FF_1510_05823_MuIn.loc[:,r'$\Delta d_{n}[efm] t_{sink}=0.9876fm$']
FF_1510_05823_plot.at['key_select',this_lab] = first_key
FF_1510_05823_plot.at['label',this_lab] = this_lab
FF_1510_05823_plot.at['type',this_lab] = 'error_bar_vary'

this_lab = FF_1510_05823.name+r' Neutron $t_{sink}=1.1522fm$'
FF_1510_05823_plot.at[:,this_lab] =  null_series
FF_1510_05823_plot.at['x_data',this_lab] = 'from_keys'
FF_1510_05823_plot.at['y_data',this_lab] = FF_1510_05823_MuIn.loc[:,r'$d_{n}[efm] t_{sink}=1.1522fm$']
FF_1510_05823_plot.at['yerr_data',this_lab] = FF_1510_05823_MuIn.loc[:,r'$\Delta d_{n}[efm] t_{sink}=1.1522fm$']
FF_1510_05823_plot.at['key_select',this_lab] = first_key
FF_1510_05823_plot.at['label',this_lab] = this_lab
FF_1510_05823_plot.at['type',this_lab] = 'error_bar_vary'

this_lab = FF_1510_05823.name+r' Neutron $t_{sink}=0.9876fm$ DRM'
FF_1510_05823_plot.loc[:,this_lab] =  null_series
FF_1510_05823_plot.at['x_data',this_lab] = 'from_keys'
FF_1510_05823_plot.at['y_data',this_lab] = FF_1510_05823_MuIn.loc[:,r'$d_{n}[efm] t_{sink}=0.9876fm$ DRM']
FF_1510_05823_plot.at['yerr_data',this_lab] = FF_1510_05823_MuIn.loc[:,r'$\Delta d_{n}[efm] t_{sink}=0.9876fm$ DRM']
FF_1510_05823_plot.at['key_select',this_lab] = first_key
FF_1510_05823_plot.at['label',this_lab] = this_lab
FF_1510_05823_plot.at['type',this_lab] = 'error_bar_vary'

this_lab = FF_1510_05823.name+r' Neutron $t_{sink}=1.1522fm$ DRM'
FF_1510_05823_plot.at[:,this_lab] =  null_series
FF_1510_05823_plot.at['x_data',this_lab] = 'from_keys'
FF_1510_05823_plot.at['y_data',this_lab] = FF_1510_05823_MuIn.loc[:,r'$d_{n}[efm] t_{sink}=1.1522fm$ DRM']
FF_1510_05823_plot.at['yerr_data',this_lab] = FF_1510_05823_MuIn.loc[:,r'$\Delta d_{n}[efm] t_{sink}=1.1522fm$ DRM']
FF_1510_05823_plot.at['key_select',this_lab] = first_key
FF_1510_05823_plot.at['label',this_lab] = this_lab
FF_1510_05823_plot.at['type',this_lab] = 'error_bar_vary'

this_lab = FF_1510_05823.name+r' Neutron $t_{sink}=0.9876fm$ EMPR'
FF_1510_05823_plot.loc[:,this_lab] =  null_series
FF_1510_05823_plot.at['x_data',this_lab] = 'from_keys'
FF_1510_05823_plot.at['y_data',this_lab] = FF_1510_05823_MuIn.loc[:,r'$d_{n}[efm] t_{sink}=0.9876fm$ EMPR']
FF_1510_05823_plot.at['yerr_data',this_lab] = FF_1510_05823_MuIn.loc[:,r'$\Delta d_{n}[efm] t_{sink}=0.9876fm$ EMPR']
FF_1510_05823_plot.at['key_select',this_lab] = first_key
FF_1510_05823_plot.at['label',this_lab] = this_lab
FF_1510_05823_plot.at['type',this_lab] = 'error_bar_vary'

this_lab = FF_1510_05823.name+r' Neutron $t_{sink}=1.1522fm$ EMPR'
FF_1510_05823_plot.at[:,this_lab] =  null_series
FF_1510_05823_plot.at['x_data',this_lab] = 'from_keys'
FF_1510_05823_plot.at['y_data',this_lab] = FF_1510_05823_MuIn.loc[:,r'$d_{n}[efm] t_{sink}=1.1522fm$ EMPR']
FF_1510_05823_plot.at['yerr_data',this_lab] = FF_1510_05823_MuIn.loc[:,r'$\Delta d_{n}[efm] t_{sink}=1.1522fm$ EMPR']
FF_1510_05823_plot.at['key_select',this_lab] = first_key
FF_1510_05823_plot.at['label',this_lab] = this_lab
FF_1510_05823_plot.at['type',this_lab] = 'error_bar_vary'

########################################################################
FF_1502_02295 = pa.DataFrame()
FF_1502_02295.name = r'F.-K. Guo et al, 2015, $N_{f}=2+1$ Clover'
FF_1502_02295.loc[:,r'$nxyz$'] = pa.Series([24,24])
FF_1502_02295.loc[:,r'$nt$'] = pa.Series([24*2,24*2])
FF_1502_02295.loc[:,r'$a[fm]$'] = pa.Series([0.074,0.074])
FF_1502_02295.loc[:,r'$\Delta a[fm]$'] = pa.Series([0.002,0.002])
FF_1502_02295.loc[:,r'$a^{-1}[GeV]$'] = hbarc/FF_1502_02295.loc[:,r'$a[fm]$']
FF_1502_02295.loc[:,r'$Volume[fm]^3$'] = (FF_1502_02295.loc[:,r'$a[fm]$'] * FF_1502_02295.loc[:,r'$nxyz$'])**3


FF_1502_02295.loc[:,r'$m_{\pi}[GeV]$'] = pa.Series([0.465,0.360])
FF_1502_02295.loc[:,r'$\Delta m_{\pi}[GeV]$'] = pa.Series([0.013,0.010])

FF_1502_02295.loc[:,r'$m_{\pi}[MeV]$'] = FF_1502_02295.loc[:,r'$m_{\pi}[GeV]$']*1000
FF_1502_02295.loc[:,r'$\Delta m_{\pi}[MeV]$'] = FF_1502_02295.loc[:,r'$\Delta m_{\pi}[GeV]$']*1000

FF_1502_02295.loc[:,r'$m_{N} [GeV]$'] = pa.Series([1.246,1.138])
FF_1502_02295.loc[:,r'$\Delta m_{N} [GeV]$'] = pa.Series([0.007,0.013])

## So dodgy..... esitmated from graph aparently...
FF_1502_02295.loc[:,r'$alpha$'] = pa.Series([-0.079,-0.092])
FF_1502_02295.loc[:,r'$\Delta alpha$'] = pa.Series([0.027,0.014])

## Taken from P. E. Shanahan et al., Phys. Rev. D89, 074511 (2014).
FF_1502_02295.loc[:,r'$F_{2}$'] = pa.Series([-1.491,-1.473])
FF_1502_02295.loc[:,r'$\Delta F_{2}$'] = pa.Series([0.022,0.037])

FF_1502_02295.loc[:,r'$F_{3}$'] = pa.Series([-0.375,-0.248])
FF_1502_02295.loc[:,r'$\Delta F_{3}$'] = pa.Series([0.048,0.029])

FF_1502_02295.loc[:,r'$\tilde{F}_{3}$'] = pa.Series([-0.130,0.020])
# FF_1502_02295.loc[:,r'$\Delta \tilde{F}_{3}$'] = pa.Series([0.076,0.058])
FF_1502_02295.loc[:,r'$\Delta \tilde{F}_{3}$'] = pa.Series([0.048,0.029])
# FF_1502_02295.loc[:,r'$\Delta \tilde{F}_{3}$'] = FF_1502_02295.loc[:,r'$\Delta F_{3}$']*FF_1502_02295.loc[:,r'$\tilde{F}_{3}$']/FF_1502_02295.loc[:,r'$F_{3}$']

FF_1502_02295.loc[:,r'$d_{n}[efm]$'] = FF_1502_02295.loc[:,r'$F_{3}$']/2*FF_1502_02295.loc[:,r'$m_{N} [GeV]$']
FF_1502_02295.loc[:,r'$\Delta d_{n}[efm]$'] = abs(FF_1502_02295.loc[:,r'$\Delta F_{3}$']/FF_1502_02295.loc[:,r'$F_{3}$'])
FF_1502_02295.loc[:,r'$\Delta d_{n}[efm]$'] += abs(FF_1502_02295.loc[:,r'$\Delta m_{N} [GeV]$']/FF_1502_02295.loc[:,r'$m_{N} [GeV]$'])
FF_1502_02295.loc[:,r'$\Delta d_{n}[efm]$'] = FF_1502_02295.loc[:,r'$\Delta d_{n}[efm]$'] * FF_1502_02295.loc[:,r'$d_{n}[efm]$']

FF_1502_02295.loc[:,r'$\tilde{d}{n}[efm]$'] = FF_1502_02295.loc[:,r'$\tilde{F}_{3}$']/2*FF_1502_02295.loc[:,r'$m_{N} [GeV]$']
FF_1502_02295.loc[:,r'$\Delta \tilde{d}{n}[efm]$'] = FF_1502_02295.loc[:,r'$\Delta \tilde{F}_{3}$']/FF_1502_02295.loc[:,r'$\tilde{F}_{3}$']
FF_1502_02295.loc[:,r'$\Delta \tilde{d}{n}[efm]$'] += FF_1502_02295.loc[:,r'$\Delta m_{N} [GeV]$']/FF_1502_02295.loc[:,r'$m_{N} [GeV]$']
FF_1502_02295.loc[:,r'$\Delta \tilde{d}{n}[efm]$'] = FF_1502_02295.loc[:,r'$\Delta \tilde{d}{n}[efm]$'] * FF_1502_02295.loc[:,r'$\tilde{d}{n}[efm]$']

xlab = r'$m_{\pi}[MeV]$'
xplot = FF_1502_02295.loc[:,xlab]


FF_1502_02295_plot = pa.DataFrame()
this_lab = FF_1502_02295.name+r' Non_rotated'
FF_1502_02295_plot.loc[:,this_lab] =  null_series
FF_1502_02295_plot.at['x_data',this_lab] = xplot
FF_1502_02295_plot.at['y_data',this_lab] = FF_1502_02295.loc[:,r'$d_{n}[efm]$']
FF_1502_02295_plot.at['yerr_data',this_lab] = FF_1502_02295.loc[:,r'$\Delta d_{n}[efm]$']
FF_1502_02295_plot.at['label',this_lab] = this_lab
FF_1502_02295_plot.at['type',this_lab] = 'error_bar'


this_lab = FF_1502_02295.name+r' Rotated'
FF_1502_02295_plot.loc[:,this_lab] =  null_series
FF_1502_02295_plot.at['x_data',this_lab] = xplot
FF_1502_02295_plot.at['y_data',this_lab] = FF_1502_02295.loc[:,r'$\tilde{d}{n}[efm]$']
FF_1502_02295_plot.at['yerr_data',this_lab] = FF_1502_02295.loc[:,r'$\Delta \tilde{d}{n}[efm]$']
FF_1502_02295_plot.at['label',this_lab] = this_lab
FF_1502_02295_plot.at['type',this_lab] = 'error_bar'



########################################################################
FF_0808_1428 = pa.DataFrame()
FF_0808_1428.name = 'R.Horsley et al, 2008, $N_{f}=2$ Clover'
FF_0808_1428.loc[:,r'$nxyz$'] = pa.Series([16])
FF_0808_1428.loc[:,r'$nt$'] = pa.Series([32*2])
FF_0808_1428.loc[:,r'$a[fm]$'] = pa.Series([0.11])
FF_0808_1428.loc[:,r'$\Delta a[fm]$'] = pa.Series([0.005])
FF_0808_1428.loc[:,r'$a^{-1}[GeV]$'] = hbarc/FF_0808_1428.loc[:,r'$a[fm]$']
FF_0808_1428.loc[:,r'$Volume[fm]^3$'] = (FF_0808_1428.loc[:,r'$a[fm]$'] * FF_0808_1428.loc[:,r'$nxyz$'])**3

## quoted as m_pi / m_rho, partical data base quotes m_rho approximatly 755 MeV
FF_0808_1428.loc[:,r'$m_{\pi}[GeV]$'] = pa.Series([0.755*0.8])
FF_0808_1428.loc[:,r'$\Delta m_{\pi}[GeV]$'] = pa.Series([0.0005*0.8])

FF_0808_1428.loc[:,r'$m_{\pi}[MeV]$'] = FF_0808_1428.loc[:,r'$m_{\pi}[GeV]$']*1000
FF_0808_1428.loc[:,r'$\Delta m_{\pi}[MeV]$'] = FF_0808_1428.loc[:,r'$\Delta m_{\pi}[GeV]$']*1000

FF_0808_1428.loc[:,r'$d_{n}[efm] \bar{\theta}=0.4$'] = pa.Series([-0.088])
FF_0808_1428.loc[:,r'$\Delta d_{n}[efm] \bar{\theta}=0.4$'] = pa.Series([0.008])
FF_0808_1428.loc[:,r'$d_{n}[efm] \bar{\theta}=0.2$'] = pa.Series([-0.190])
FF_0808_1428.loc[:,r'$\Delta d_{n}[efm] \bar{\theta}=0.2$'] = pa.Series([0.009])
FF_0808_1428.loc[:,r'$d_{n}[efm] \bar{\theta}\rightarrow 0$'] = pa.Series([-0.049])
FF_0808_1428.loc[:,r'$\Delta d_{n}[efm] \bar{\theta}\right0$'] = pa.Series([0.005])

xlab = r'$m_{\pi}[MeV]$'
xplot = FF_0808_1428.loc[:,xlab]


FF_0808_1428_plot = pa.DataFrame()
this_lab = FF_0808_1428.name + r'$\bar{\theta}=0.4$'
FF_0808_1428_plot.loc[:,this_lab] =  null_series
FF_0808_1428_plot.at['x_data',this_lab] = xplot
FF_0808_1428_plot.at['y_data',this_lab] = FF_0808_1428.loc[:,r'$d_{n}[efm] \bar{\theta}=0.4$']
FF_0808_1428_plot.at['yerr_data',this_lab] = FF_0808_1428.loc[:,r'$d_{n}[efm] \bar{\theta}=0.4$']
FF_0808_1428_plot.at['label',this_lab] = this_lab
FF_0808_1428_plot.at['type',this_lab] = 'error_bar'


this_lab = FF_0808_1428.name + r'$\bar{\theta}=0.2$'
FF_0808_1428_plot.loc[:,this_lab] =  null_series
FF_0808_1428_plot.at['x_data',this_lab] = xplot
FF_0808_1428_plot.at['y_data',this_lab] = FF_0808_1428.loc[:,r'$d_{n}[efm] \bar{\theta}=0.2$']
FF_0808_1428_plot.at['yerr_data',this_lab] = FF_0808_1428.loc[:,r'$d_{n}[efm] \bar{\theta}=0.2$']
FF_0808_1428_plot.at['label',this_lab] = this_lab
FF_0808_1428_plot.at['type',this_lab] = 'error_bar'


this_lab = FF_0808_1428.name + r'$\bar{\theta}\rightarrow 0$'
FF_0808_1428_plot.loc[:,this_lab] =  null_series
FF_0808_1428_plot.at['x_data',this_lab] = xplot
FF_0808_1428_plot.at['y_data',this_lab] = FF_0808_1428.loc[:,r'$d_{n}[efm] \bar{\theta}\rightarrow 0$']
FF_0808_1428_plot.at['yerr_data',this_lab] = FF_0808_1428.loc[:,r'$d_{n}[efm] \bar{\theta}\rightarrow 0$']
FF_0808_1428_plot.at['label',this_lab] = this_lab
FF_0808_1428_plot.at['type',this_lab] = 'error_bar'


########################################################################################################################
FF_0512004 = pa.DataFrame()
FF_0512004.name = 'F. Berruto et al, 2005, $N_{f}=2$, DWF'
FF_0512004.loc[:,r'$nxyz$'] = pa.Series([16,16])
FF_0512004.loc[:,r'$nt$'] = pa.Series([16*2,16*2])
FF_0512004.loc[:,r'$a^{-1}[GeV]$'] = pa.Series([1.7,1.7])
FF_0512004.loc[:,r'$a[fm]$'] = hbarc/FF_0512004.loc[:,r'$a^{-1}[GeV]$']
FF_0512004.loc[:,r'$Volume[fm]^3$'] = (FF_0512004.loc[:,r'$a[fm]$'] * FF_0512004.loc[:,r'$nxyz$'])**3

## only m_s is quoted.... shifting to -100
FF_0512004.loc[:,r'$m_{\pi}[GeV]$'] = pa.Series([0.690,0.605])
FF_0512004.loc[:,r'$m_{\pi}[MeV]$'] = FF_0512004.loc[:,r'$m_{\pi}[GeV]$']*1000

FF_0512004.loc[:,r'$m_{N} fitr$'] = pa.Series(['fitr7-12','fitr7-12'])
FF_0512004.loc[:,r'$m_{N}$'] = pa.Series([0.8646,0.9264])
FF_0512004.loc[:,r'$\Delta m_{N}$'] = pa.Series([0.0053,0.0054])
FF_0512004.loc[:,r'$m_{N} [GeV]$'] = FF_0512004.loc[:,r'$m_{N}$']*FF_0512004.loc[:,r'$a^{-1}[GeV]$']
FF_0512004.loc[:,r'$\Delta m_{N} [GeV]$'] = FF_0512004.loc[:,r'$\Delta m_{N}$']*FF_0512004.loc[:,r'$a^{-1}[GeV]$']
#
# FF_0512004.loc[:,r'$alpha fitr$'] = pa.Series(['fitr5-9',r'$fitr5-9',r'$fitr5-9$'])
# FF_0512004.loc[:,r'$alpha$'] = pa.Series([-0.370,-0.356,-0.333])
# FF_0512004.loc[:,r'$\Delta alpha$'] = pa.Series([0.022,0.022,0.128])

## so bad approximation.... they just took the smallest momentum and said its smaller than that
FF_0512004.loc[:,r'$d_{p}[efm]$ Eq.41'] = pa.Series([0.0,0.0])
FF_0512004.loc[:,r'$\Delta d_{p}[efm]$ Eq.41'] = pa.Series([0.20 + 0.04,0.087+0.095])
FF_0512004.loc[:,r'$d_{p}[efm]$ Eq.42'] = pa.Series([0.0,0.0])
FF_0512004.loc[:,r'$\Delta d_{p}[efm]$ Eq.42'] = pa.Series([0.49+0.45,0.12+0.27])

xlab = r'$m_{\pi}[MeV]$'
xplot = FF_0512004.loc[:,xlab]


FF_0512004_plot = pa.DataFrame()
this_lab = FF_0512004.name+r' Neutron Eq.41'
FF_0512004_plot.loc[:,this_lab] =  null_series
FF_0512004_plot.at['x_data',this_lab] = xplot
FF_0512004_plot.at['y_data',this_lab] = FF_0512004.loc[:,r'$d_{p}[efm]$ Eq.41']
FF_0512004_plot.at['yerr_data',this_lab] = FF_0512004.loc[:,r'$\Delta d_{p}[efm]$ Eq.41']
FF_0512004_plot.at['label',this_lab] = this_lab
FF_0512004_plot.at['type',this_lab] = 'error_bar'


this_lab = FF_0512004.name+r' Neutron Eq.42'
FF_0512004_plot.loc[:,this_lab] =  null_series
FF_0512004_plot.at['x_data',this_lab] = xplot
FF_0512004_plot.at['y_data',this_lab] = FF_0512004.loc[:,r'$d_{p}[efm]$ Eq.42']
FF_0512004_plot.at['yerr_data',this_lab] = FF_0512004.loc[:,r'$\Delta d_{p}[efm]$ Eq.42']
FF_0512004_plot.at['label',this_lab] = this_lab
FF_0512004_plot.at['type',this_lab] = 'error_bar'

########################################################################################################
FF_0505022 = pa.DataFrame()
FF_0505022.name = r'E. Shintani et al, 2005, $N_{f}=0$, DWF'
FF_0505022.loc[:,r'$nxyz$'] = pa.Series([16])
FF_0505022.loc[:,r'$nt$'] = pa.Series([16*2])
FF_0505022.loc[:,r'$a^{-1}[GeV]$'] = pa.Series([1.875])
FF_0505022.loc[:,r'$a[fm]$'] = hbarc/FF_0505022.loc[:,r'$a^{-1}[GeV]$']
FF_0505022.loc[:,r'$Volume[fm]^3$'] = (FF_0505022.loc[:,r'$a[fm]$'] * FF_0505022.loc[:,r'$nxyz$'])**3


FF_0505022.loc[:,r'$m_{\pi}[GeV]$'] = pa.Series([0.755*0.63])
FF_0505022.loc[:,r'$m_{\pi}[MeV]$'] = FF_0505022.loc[:,r'$m_{\pi}[GeV]$']*1000

FF_0505022.loc[:,r'$m_{N} fitr$'] = pa.Series(['fitr9-13'])
FF_0505022.loc[:,r'$am_{N}$'] = pa.Series([0.7114])
FF_0505022.loc[:,r'$\Delta am_{N}$'] = pa.Series([0.0043])

FF_0505022.loc[:,r'$m_{N} [GeV]$'] = FF_0505022.loc[:,r'$am_{N}$'] * FF_0505022.loc[:,r'$a^{-1}[GeV]$']
FF_0505022.loc[:,r'$\Delta m_{N} [GeV]$'] = FF_0505022.loc[:,r'$\Delta am_{N}$'] * FF_0505022.loc[:,r'$a^{-1}[GeV]$']

## quoted as f^1_N in paper
FF_0505022.loc[:,r'$alpha$'] = pa.Series([-0.247])
FF_0505022.loc[:,r'$\Delta alpha$'] = pa.Series([0.017])

FF_0505022.loc[:,r'$F_{2}$'] = pa.Series([-0.560])
FF_0505022.loc[:,r'$\Delta F_{2}$'] = pa.Series([0.040])

## Not Extrapolatied
FF_0505022.loc[:,r'$d_{n}[efm] (q^{2}\approx 0.58GeV^{2})$'] = pa.Series([-0.024])
FF_0505022.loc[:,r'$\Delta d_{n}[efm]$ (q^{2}\approx 0.58GeV^{2})'] = pa.Series([0.005])

dn = FF_0505022.loc[:,r'$d_{n}[efm] (q^{2}\approx 0.58GeV^{2})$']
d_dn = FF_0505022.loc[:,r'$\Delta d_{n}[efm]$ (q^{2}\approx 0.58GeV^{2})']
mn = FF_0505022.loc[:,r'$m_{N} [GeV]$']
d_mn = FF_0505022.loc[:,r'$\Delta alpha$']
alpha = FF_0505022.loc[:,r'$alpha$']
d_alpha = FF_0505022.loc[:,r'$\Delta alpha$']
F2 = FF_0505022.loc[:,r'$F_{2}$']
d_F2 = FF_0505022.loc[:,r'$\Delta F_{2}$']

FF_0505022.loc[:,r'$\tilde{d}_{n}[efm] (q^{2}\approx 0.58GeV^{2})$'] = dn*2*mn/hbarc + 2*alpha*F2

# FF_0505022.loc[:,r'$\Delta \tilde{d}_{n}[efm] (q^{2}\approx 0.58GeV^{2})$'] = \
#     d_dn*2*mn/hbarc + dn*2*d_mn/hbarc + 2*d_alpha*F2 + 2*alpha*d_F2
FF_0505022.loc[:,r'$\Delta \tilde{d}_{n}[efm] (q^{2}\approx 0.58GeV^{2})$'] = \
FF_0505022.loc[:,r'$\Delta d_{n}[efm]$ (q^{2}\approx 0.58GeV^{2})'] * FF_0505022.loc[:,r'$\tilde{d}_{n}[efm] (q^{2}\approx 0.58GeV^{2})$']/FF_0505022.loc[:,r'$d_{n}[efm] (q^{2}\approx 0.58GeV^{2})$']



xlab = r'$m_{\pi}[MeV]$'
xplot = FF_0505022.loc[:,xlab]


FF_0505022_plot = pa.DataFrame()
this_lab = FF_0505022.name+r'Neutron $q^{2}\approx 0.58GeV^{2}$ Unrotated'
FF_0505022_plot.loc[:,this_lab] =  null_series
FF_0505022_plot.at['x_data',this_lab] = xplot
FF_0505022_plot.at['y_data',this_lab] = FF_0505022.loc[:,r'$d_{n}[efm] (q^{2}\approx 0.58GeV^{2})$']
FF_0505022_plot.at['yerr_data',this_lab] = FF_0505022.loc[:,r'$\Delta d_{n}[efm]$ (q^{2}\approx 0.58GeV^{2})']
FF_0505022_plot.at['label',this_lab] = this_lab
FF_0505022_plot.at['type',this_lab] = 'error_bar'

this_lab = FF_0505022.name+r'Neutron $q^{2}\approx 0.58GeV^{2}$ Rotated'
FF_0505022_plot.loc[:,this_lab] =  null_series
FF_0505022_plot.at['x_data',this_lab] = xplot
FF_0505022_plot.at['y_data',this_lab] = FF_0505022.loc[:,r'$\tilde{d}_{n}[efm] (q^{2}\approx 0.58GeV^{2})$']
FF_0505022_plot.at['yerr_data',this_lab] = FF_0505022.loc[:,r'$\Delta \tilde{d}_{n}[efm] (q^{2}\approx 0.58GeV^{2})$']
FF_0505022_plot.at['label',this_lab] = this_lab
FF_0505022_plot.at['type',this_lab] = 'error_bar'
