
from XmlFormatting import MakeValAndErr

MakeValAndErr(1.2118170378   ,     0.8221137640,Dec=2)

'''
            Proton          Neutron
EDM         3.6(1.9)        1.21(82)
Schiff      2.7(1.5)        1.11(80)
'''

import pickle as pik
import numpy as np
from MomParams import hbarc
hbarcM = hbarc*1000

#### NEUTRON ####
mN = 939.5653 # neutron mass in MeV


def ReadCs(this_file):
    this_data = pik.load(open(this_file,'rb'))

    this_fit = this_data['plot_data']['fit _2']['fit_class']['t_f6.01','frac{tilde{F}_{3}}{2m_{N}}','fitr0-4','$F(0)$']

    par_list = this_fit.fit_data['Params']

    C1,C2,C3 = par_list.values
    C1 = -C1
    print(' In fm MeV^-2 (first table)')
    print('     C1: '+C1.MakeValAndErr(Dec=2))
    C2 = -C2
    print('     C2: '+C2.MakeValAndErr(Dec=2))



    C2_term = C2*np.log(mN**2)
    print('     C2_term: '+C2_term.MakeValAndErr(Dec=2))

    C1p = C1 + C2_term
    print('     C1\': '+C1p.MakeValAndErr(Dec=2))

    C1_fm = C1 *hbarcM**2
    C2_fm = C2 *hbarcM**2
    C2_term_fm = C2_term *hbarcM**2
    C1p_fm = C1p* hbarcM**2

    print()
    print(' In fm^3 (second table)')
    print('     C1: '+C1_fm.MakeValAndErr(Dec=2))
    print('     C2: '+C2_fm.MakeValAndErr(Dec=2))
    print('     C2*log(mN**2): '+C2_term_fm.MakeValAndErr(Dec=2))
    print('     C1\': '+C1p_fm.MakeValAndErr(Dec=2))
    print()


this_file = '/home/jackdra/LQCD/Scripts/EDM_paper/graphs/FF/FullFFFit/Neutron_ContFit_mpi.py3p'
print('Neutron')
ReadCs(this_file)

new_fit_file = '/home/jackdra/LQCD/Results/HPCCRES_FixedG2/graphs/FF/Set1_F0VsMpi.py3p'
nf_data = pik.load(open(new_fit_file,'rb'))
nf_data = nf_data['plot_data']['fit ']['fit_class']
nf_data = nf_data['t_f5.01','frac{tilde{F}_{3}}{2m_{N}}','fitr0-4','$F(0)$']

par_list = nf_data.fit_data['Params']

C1,C2,C3 = par_list.values

print('***')
print('Neutron fit with mN in log')
print(' In fm MeV^-2 (first table)')
print('     C1\': '+C1.MakeValAndErr(Dec=2))
print('     C2: '+C2.MakeValAndErr(Dec=2))

C1_fm = C1 *hbarcM**2
C2_fm = C2 *hbarcM**2

print()
print(' In fm^3 (second table)')
print('     C1\': '+C1_fm.MakeValAndErr(Dec=2))
print('     C2: '+C2_fm.MakeValAndErr(Dec=2))
print('***')
print()


this_file = '/home/jackdra/LQCD/Scripts/EDM_paper/graphs/FF/FullFFFit/Proton_ContFit_mpi.py3p'
print('Proton')
ReadCs(this_file)


M_labels = ['mpi410','mpi570','mpi700']
A_labels = ['L16','L20','L28']
# M_labels = ['mpi410']
# A_labels = []
all_ens = M_labels + A_labels
tsink_ens_dict = {}
tsink_ens_dict['mpi410'] =  64
tsink_ens_dict['mpi570'] =  64
tsink_ens_dict['mpi700'] =  64
tsink_ens_dict['L16'] =  32
tsink_ens_dict['L20'] =  40
tsink_ens_dict['L28'] =  56

a_dict = {}
a_dict['mpi410'] =  0.0907
a_dict['mpi570'] =  0.0907
a_dict['mpi700'] =  0.0907
a_dict['L16'] =  0.1095
a_dict['L20'] =  0.0936
a_dict['L28'] =  0.0684
t0_dict = {}
t0_dict['mpi410'] =     2.5377162290
t0_dict['mpi570'] =     2.3999592122
t0_dict['mpi700'] =     2.2590866059
t0_dict['L16'] =        1.3628954894
t0_dict['L20'] =        2.2386627909
t0_dict['L28'] =        4.9885618908

t0_dict_err = {}
t0_dict_err['mpi410'] =  0.0015529386
t0_dict_err['mpi570'] =  0.0011452810
t0_dict_err['mpi700'] =  0.0011742296
t0_dict_err['L16'] =     0.0015817442
t0_dict_err['L20'] =     0.0021913620
t0_dict_err['L28'] =     0.0064980950

for iens in all_ens:
    print(iens,np.sqrt(8*t0_dict[iens])*a_dict[iens])


def latTOphys(tf,this_ens):
    return np.sqrt(8*tf)*a_dict[this_ens]

def physTOlat(tf,this_ens):
    return (tf/a_dict[this_ens])**2 / 8

phys_fix = 0.60

lat_fix = []
phys_fix_close = []
for iens in all_ens:
    lat_fix.append(round(physTOlat(phys_fix,iens),2))
    phys_fix_close.append(latTOphys(lat_fix[-1],iens))
    print(iens,lat_fix[-1],phys_fix_close[-1])
