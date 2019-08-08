import numpy as np
from XmlFormatting import MakeValAndErr
from FitFunctions import Fitting
from MomParams import hbarc
from copy import deepcopy
from BootStrapping import BootStrap
from FileIO import ReadPicklePy2
this_path = '/home/jackdra/LQCD/Results/HPCCRES_FixedG2/test_data/ig3g5/VectorTop_t_f5.47_Rtsumfitr7-32_Afitr10-20_RFmin3_RChi50-100.p'
this_ff = ReadPicklePy2(this_path)
this_ff = this_ff['thisFF'].__dict__
alpha = list(this_ff['alpha'].values())[0]
alpha = list(alpha.values())[0]
a_mN = list(this_ff['mass'].values())[0]
this_coeff = this_ff['paramRed_Stats']
neg_P3 = True
q_list = {
'P3_g1':['q0-10'],
'P3_g3_Top_cmplx':['q-100','q00-1'],
'P3_g4_Top':['q00-1'],
'P4_g1_cmplx':['q-100'],
'P4_g4':['q-100']
}



A = []
A_avg = []
for ikey,iqlist in q_list.items():
    for iq in iqlist:
        this_key = (ikey,'pp000',iq,
                    'C2_state1_fitr9-20_alpha_fitr10-20',
                    't_f5.47')
        if neg_P3 and 'P3_' in ikey:
            A.append(list(-this_coeff['FF'][this_key].values))
        else:
            A.append(list(this_coeff['FF'][this_key].values))
        A_avg.append([ia.GetAvg() for ia in A[-1]])

A_avg[1] = [-0.004091203310362446, -0.0037999597699054553, -0.0015159981114171836]
A_avg[2] = [-0.004091203310362446, 0.0004806567910960839, -0.023532475648113317]
A_avg[3] = [0.02771624505955395, -0.0032562550430615955, 0.15942298938532015]

A[1] = [-0.004091203310362446, -0.0037999597699054553, -0.0015159981114171836]
A[2] = [-0.004091203310362446, 0.0004806567910960839, -0.023532475648113317]
A[3] = [0.02771624505955395, -0.0032562550430615955, 0.15942298938532015]

27.549757627311678
27.5561888037325

74.7513745345211
74.7568917413095

0.0003063246806996577
0.00020490889341879124
# alpha = -0.1898009109
# a_mN =  0.6506026910
a = 0.0907
renorm = 0.7354
this_qunit = 2*np.pi/32
this_Ep = (a_mN**2 + this_qunit**2)**(1/2)
thisQ2 = this_qunit**2 - (this_Ep - a_mN)**2
# from GammaMatricies import CorrSpinTrace
# pmu,qmu,ppmu = CorrSpinTrace().Create4Mom([-this_qunit,0,0],[0,0,0],a_mN.Avg)
#
# np.array(qmu[0])**2
# (-(this_Ep - a_mN)**2).MakeValAndErr(Dec=2)
#
# this_qunit**2
# np.sum(np.array(qmu[1:-1])**2)
#
# np.sum(np.array(qmu[:4])**2)
# thisQ2.MakeValAndErr(Dec=2)

Q2_GeV =  thisQ2*(hbarc/a)**2
GE_F2_coeff = - thisQ2/(4*a_mN**2)
GE_F2_coeff.Stats()
GE_F2_coeff.MakeValAndErr(Dec=2)
RQ_full_list = [1,2,3]
RQ_solve_list = [1,2,3]
print('------New Results, Rimp, sqrt(8t_f)= 0.6 fm------')
print('Q2_GeV = ',Q2_GeV.MakeValAndErr(Dec=2))
print('m_pi = 570 MeV')
print(' Solving with RQ solve list: ',RQ_solve_list)
exclude_list = [i for i in RQ_full_list if i not in RQ_solve_list]
def lin_fun(x,y):
    return sum([ix*ip for ix,ip in zip(x,y)])

fit_data = Fitting(Funs=(lin_fun,3))
fit_data_vec = Fitting(Funs=(lin_fun,2))

def AppendGeGm(this_ff):
    F1F2_list = this_ff[:2]
    GE = F1F2_list[0] + GE_F2_coeff*F1F2_list[1]
    GM = F1F2_list[0] + F1F2_list[1]
    GE.Stats()
    # gegm_list[0] = gegm_list[0].makevalanderr(dec=2)
    return np.array(list(this_ff[:2])+[GE,GM]+list(this_ff[2:]))


def alpha_rot(F_list):
    c2a = np.cos(2*alpha)
    s2a = np.sin(2*alpha)
    out = np.array([F_list[0],
                    c2a*F_list[1] - s2a*F_list[2],
                    s2a*F_list[1] + c2a*F_list[2]])
    if isinstance(out,BootStrap):
        out.Stats()
        return out.Avg
    else:
        return out

def alpha_rot_F2tilde(F_list):
    c2a = np.cos(2*alpha)
    s2a = np.sin(2*alpha)
    t2a = s2a/c2a
    out = np.array([F_list[0],F_list[1],
                     (c2a+t2a*s2a)*F_list[2] + t2a*F_list[1]])
    if isinstance(out,BootStrap):
        out.Stats()
        return out.Avg
    else:
        return out
def fmt_boot(val,doubbrak=False,width=13):
    if isinstance(val,(list,tuple,np.ndarray)):
        out_val = []
        for ival in val:
            if isinstance(ival,BootStrap):
                out_val.append(ival.MakeValAndErr(Dec=2))
                if doubbrak:
                    out_val[-1] = out_val[-1].replace('(','((').replace(')','))')
            elif isinstance(ival,float):
                out_val.append(f'{ival:.11}')
            else:
                out_val.append(ival)
                if doubbrak:
                    out_val[-1] = out_val[-1].replace('(','((').replace(')','))')
            this_val = out_val[-1]
            out_val[-1] = f'{this_val: <{width}}'
    else:
        if isinstance(val,BootStrap):
            out_val = val.MakeValAndErr(Dec=2)
            if doubbrak:
                out_val = out_val.replace('(','((').replace(')','))')
        elif isinstance(val,float):
            out_val = f'{val:.11}'
        else:
            out_val = val
            if doubbrak:
                out_val = out_val.replace('(','((').replace(')','))')
        out_val = f'{out_val: <{width}}'
    return out_val

## from file /mnt/scratch/dragosja/data/resultsFixedG2/RC32x64Kud01372700Ks01364000/FormFactors/Proton_qmax0_ppmax4/VectorTop_t_f5.47_Rtsumfitr7-32_Afitr10-20_RFmin3_RChi50-100.xml
## Q2 = 1

opp_list = [
'P3_g1_q0-10',
'P3_g3_Top_cmplx_q-100',
'P3_g3_Top_cmplx_q00-1',
'P3_g4_Top_q00-1',
'P4_g1_cmplx_q-100',
'P4_g4_q-100',
]

# A = [
# [-0.14603  , -0.14603  , 0.00000e+00 ],
# [-0.00409  , -0.00378  ,  0.00154     ],
# [-0.00409  ,  0.00050  ,  0.02356     ],
# [ 0.02772  , -0.00341  , -0.15958     ],
# [-0.14603  ,  0.00325  , 0.00000e+00 ],
# [ 0.98928  , -0.02204  , 0.00000e+00 ]
# ]
#
# A_err = [
# [0.00026, 0.00026, 0.00000e+00],
# [0.00057, 0.00044, 0.00042],
# [0.00057, 0.00019, 0.00041],
# [0.00389, 0.00127, 0.00277],
# [0.00026, 0.00002, 0.00000e+00],
# [0.00004, 0.00008, 0.00000e+00]
# ]



b_p = [
 0.38063,
 0.03822,
 0.01563,
-0.03962,
-0.16227,
 1.05168]

b_p_err = [
0.00149,
0.01120,
0.00977,
0.01239,
0.00080,
0.00082]

File_FF_p_list = [
 0.8073030761,
 1.1096350169,
 0.7825866763,
 1.9169380930,
-0.0814332110,
-0.0056583279]

File_FF_p_err = [
0.0006382802 ,
0.0087835107,
0.0005807771,
0.0090787129,
0.0497750749 ,
0.0034602218 ]
FileFF_p_str = [MakeValAndErr(iff,ifferr,Dec=2) for iff,ifferr in zip(File_FF_p_list,File_FF_p_err)]

## from file /mnt/scratch/dragosja/data/resultsFixedG2/RC32x64Kud01372700Ks01364000/FormFactors/Neutron_qmax0_ppmax4/VectorTop_t_f5.47_Rtsumfitr7-32_Afitr10-20_RFmin3_RChi50-100.xml
## Q2 = 1

b_n = [
     -0.23300,
     -0.02895,
     -0.01049,
    0.02368,
    0.00561,
     0.01299
]

b_n_err = [
    0.00093,
    0.00736,
    0.00657,
    0.00738,
    0.00051,
    0.00048
]

## F1, F2tilde, F3tilde,F3tild/2mN
File_FF_N_list= [
-0.0168334012   ,
-1.1566311772    ,
 0.0089298075   ,
-1.1734645784   ,
 0.1396013106   ,
 0.0096986761   ]

File_FF_N_err = [
    0.0003455014,
    0.0053874695,
    0.0003409781,
    0.0054708420,
    0.0321334321,
    0.0022335542]
FileFF_N_str = [MakeValAndErr(iff,ifferr,Dec=2) for iff,ifferr in zip(File_FF_N_list,File_FF_N_err)]

## No Rotation F1, F2tilde, F3tilde,F3tild/2mN
File_FF_N_list_norot= [
0.0343033761   ,
1.1391602858    ,
0.3070076873   ,
0.0213260732   ]

File_FF_N_err_norot = [
    0.0003524208,
    0.0054370368,
    0.0727008079,
    0.0050411796
]
FileFF_N_str_norot = [MakeValAndErr(iff,ifferr,Dec=2) for iff,ifferr in zip(File_FF_N_list_norot,File_FF_N_err_norot)]


## F1, F2 no topological charge
## from file: /mnt/scratch/dragosja/data/resultsFixedG2/RC32x64Kud01372700Ks01364000/FormFactors/Neutron_qmax0_ppmax4/Pickle/hold_g3g5Vec/Vector_RFmin3_RChi50-100.p.FFstats

File_FF_N_list_noT= [
    -0.0168322965,
    -1.1565815591,
     0.0089298070,
    -1.1734138556]

File_FF_N_err_noT = [
    0.0003454982,
    0.0053862758,
    0.0003409774,
    0.0054696330]
FileFF_N_str_noT = [MakeValAndErr(iff,ifferr,Dec=2) for iff,ifferr in zip(File_FF_N_list_noT,File_FF_N_err_noT)]


## F1, F2 no topological charge
## from file: /mnt/scratch/dragosja/data/resultsFixedG2/RC32x64Kud01372700Ks01364000/FormFactors/Proton_qmax0_ppmax4/Pickle/hold_g3g5Vec/Vector_RFmin3_RChi50-100.p.FFstats

File_FF_p_list_noT= [
    0.8073016486,
    1.1095708905,
    0.7825866772,
    1.9168725391]

File_FF_p_err_noT = [
    0.0006381381,
    0.0087797357,
    0.0005807777,
    0.0090745190]
FileFF_p_str_noT = [MakeValAndErr(iff,ifferr,Dec=2) for iff,ifferr in zip(File_FF_p_list_noT,File_FF_p_err_noT)]



A_np = np.array(A_avg)
equ_rat = A_np[2]/A_np[3]


for ie in exclude_list: del A[ie]



R_rat_n = b_n[2]/b_n[3]
R_rat_p = b_p[2]/b_p[3]
R_rat_n_err = np.abs(R_rat_n)*(np.abs(b_n_err[2]/b_n[2]) + np.abs(b_n_err[3]/b_n[3]))
R_rat_p_err = np.abs(R_rat_p)*(np.abs(b_p_err[2]/b_p[2]) + np.abs(b_p_err[3]/b_p[3]))

for ie in exclude_list: del b_p[ie]
for ie in exclude_list: del b_n[ie]

# FFs_p,covar,FFs_p_chi2 = fit_data.LSFit(np.array(A_avg).swapaxes(0,1),np.array(b_p),np.array(b_p_err))
FFs_p,covar,FFs_p_chi2 = fit_data.LSFit(np.array(A_avg).swapaxes(0,1),np.array(b_p),[1 for i in b_p])

FFs_p = np.append(FFs_p,[FFs_p[-1]*a/(2*a_mN)])
FFrot_p = alpha_rot(FFs_p)
FFrot_p = np.append(FFrot_p,FFrot_p[-1]*a/(2*a_mN))

FFrot_p_F2t = alpha_rot_F2tilde(FFs_p)
FFrot_p_F2t = np.append(FFrot_p_F2t,FFrot_p_F2t[-1]*a/(2*a_mN))

# print(FFs_p*renorm)

b_p_old = [-ib for ib in b_p[:4]] + b_p[4:]


# FFs_p_old = np.linalg.lstsq(np.array(A_avg),np.array(b_p_old))[0]
FFs_p_old = fit_data.LSFit(np.array(A_avg).swapaxes(0,1),np.array(b_p_old),np.array(b_p_err))[0]
FFs_p_old = np.append(FFs_p_old,[FFs_p_old[-1]*a/(2*a_mN)])
FFrot_p_old = alpha_rot(FFs_p_old)
FFrot_p_old = np.append(FFrot_p_old,FFrot_p_old[-1]*a/(2*a_mN))

FFrot_p_old_F2t = alpha_rot_F2tilde(FFs_p_old)
FFrot_p_old_F2t = np.append(FFrot_p_old_F2t,FFrot_p_old_F2t[-1]*a/(2*a_mN))

# print(FFs_p_old*renorm)


# FFs_n,FFs_n_chi2 = np.linalg.lstsq(np.array(A_avg),np.array(b_n))[0:2]
FFs_n,covar,FFs_n_chi2 = fit_data.LSFit(np.array(A_avg).swapaxes(0,1),np.array(b_n),np.array(b_n_err))
FFs_n,covar,FFs_n_chi2 = fit_data.LSFit(np.array(A_avg).swapaxes(0,1),np.array(b_n),[1 for i in b_n])
FFs_n = np.append(FFs_n,[FFs_n[-1]*a/(2*a_mN)])
FFrot_n = alpha_rot(FFs_n)
FFrot_n = np.append(FFrot_n,FFrot_n[-1]*a/(2*a_mN))
FFrot_n_F2t = alpha_rot_F2tilde(FFs_n)
FFrot_n_F2t = np.append(FFrot_n_F2t,FFrot_n_F2t[-1]*a/(2*a_mN))


b_n_old = [-ib for ib in b_n[:4]] + b_n[4:]
FFs_n_old = np.linalg.lstsq(np.array(A_avg),np.array(b_n_old))[0]
FFs_n_old = np.append(FFs_n_old,[FFs_n_old[-1]*a/(2*a_mN)])
FFrot_n_old = alpha_rot(FFs_n_old)
FFrot_n_old = np.append(FFrot_n_old,FFrot_n_old[-1]*a/(2*a_mN))

A_vec = [ia[:-1] for ia in A_avg[:1]] + [ia[:-1] for ia in A_avg[-2:]]
b_p_vec = b_p[:1] + b_p[-2:]
b_n_vec = b_n[:1] + b_n[-2:]
b_p_vec_err = b_p_err[:1] + b_p_err[-2:]
b_n_vec_err = b_n_err[:1] + b_n_err[-2:]
# A_vec = [ia[:-1] for ia in A_avg[:1]] + [ia[:-1] for ia in A_avg[-1:]]
# b_p_vec = b_p[:1] + b_p[-1:]
# b_n_vec = b_n[:1] + b_n[-1:]
# b_p_vec_err = b_p_err[:1] + b_p_err[-1:]
# b_n_vec_err = b_n_err[:1] + b_n_err[-1:]
print('DEBUG for Neutron')
for ia in A_vec:
    print(ia)
print(b_n_vec)
print(b_n_vec_err)
# np.matrix(A_vec)**-1
# F1,F2 = np.array(np.dot(np.matrix(A_vec)**-1 ,b_n_vec))[0]
# print('Explicit solve',fmt_boot(AppendGeGm(np.array([F1,F2])*renorm)))

# # FFs_p_vec,FFs_p_vec_chi2 = np.linalg.lstsq(np.array(A_vec),np.array(b_p_vec))[0:2]
# FFs_p_vec,covar,FFs_p_vec_chi2 = fit_data_vec.LSFit(np.array(A_vec).swapaxes(0,1),np.array(b_p_vec),np.array(b_p_vec_err))
# # FFs_n_vec,FFs_n_vec_chi2 = np.linalg.lstsq(np.array(A_vec),np.array(b_n_vec))[0:2]
# FFs_n_vec,covar,FFs_n_vec_chi2 = fit_data_vec.LSFit(np.array(A_vec).swapaxes(0,1),np.array(b_n_vec),np.array(b_n_vec_err))

# FFs_p_vec,FFs_p_vec_chi2 = np.linalg.lstsq(np.array(A_vec),np.array(b_p_vec))[0:2]
FFs_p_vec,covar,FFs_p_vec_chi2 = fit_data_vec.LSFit(np.array(A_vec).swapaxes(0,1),
                                                    np.array(b_p_vec),[1 for i in b_p_vec])
# FFs_n_vec,FFs_n_vec_chi2 = np.linalg.lstsq(np.array(A_vec),np.array(b_n_vec))[0:2]
FFs_n_vec,covar,FFs_n_vec_chi2 = fit_data_vec.LSFit(np.array(A_vec).swapaxes(0,1),
                                                    np.array(b_n_vec),[1 for i in b_n_vec])
print(fmt_boot(AppendGeGm(FFs_n_vec*renorm)),FFs_n_vec_chi2)
print()

def SolveF3_proton(F1,F2,eq_list=RQ_solve_list,Arot=False):
    out_list = []
    for eq in range(len(eq_list)):
        if Arot:
            C1,C2,C3 = tuple(A_rot[eq+2])
        else:
            C1,C2,C3 = tuple(A[eq+2])
        R = b_p[eq+2]
        if abs(C3) < 1e-7: continue
        out_list.append((R - C1*F1 - C2*F2)/C3)
    return np.mean(out_list)

def SolveF3_neutron(F1,F2,eq_list=RQ_solve_list,Arot=False):
    out_list = []
    for eq in range(len(eq_list)):
        if Arot:
            C1,C2,C3 = tuple(A_rot[eq+2])
        else:
            C1,C2,C3 = tuple(A[eq+2])
        if abs(C3) < 1e-7: continue
        R = b_n[eq+2]
        out_list.append((R - C1*F1 - C2*F2)/C3)
    return np.mean(out_list)

FF3_p_vec = [SolveF3_proton(FFs_p_vec[0],FFs_p_vec[1]),
             SolveF3_proton(FFs_p_vec[0],FFs_p_vec[1])*a/(2*a_mN)]
FF3_n_vec = [SolveF3_neutron(FFs_n_vec[0],FFs_n_vec[1]),
             SolveF3_neutron(FFs_n_vec[0],FFs_n_vec[1])*a/(2*a_mN)]
FFs_p_vec = np.append(FFs_p_vec,FF3_p_vec)
FFs_n_vec = np.append(FFs_n_vec,FF3_n_vec)

FFrot_p_vec = alpha_rot(FFs_p_vec)
FFrot_p_vec = np.append(FFrot_p_vec,FFrot_p_vec[-1]*a/(2*a_mN))
FFrot_p_vec_F2t = alpha_rot_F2tilde(FFs_p_vec)
FFrot_p_vec_F2t = np.append(FFrot_p_vec_F2t,FFrot_p_vec_F2t[-1]*a/(2*a_mN))

FFrot_n_vec = alpha_rot(FFs_n_vec)
FFrot_n_vec = np.append(FFrot_n_vec,FFrot_n_vec[-1]*a/(2*a_mN))
FFrot_n_vec_F2t = alpha_rot_F2tilde(FFs_n_vec)
FFrot_n_vec_F2t = np.append(FFrot_n_vec_F2t,FFrot_n_vec_F2t[-1]*a/(2*a_mN))


A_rot = deepcopy(A_avg)
for ie in range(len(RQ_solve_list)):
    A_rot[ie+1][1],A_rot[ie+1][2] = (   np.cos(2*alpha)*A[ie+1][1] - np.sin(2*alpha)*A[ie+1][2],
                                        np.sin(2*alpha)*A[ie+1][1] + np.cos(2*alpha)*A[ie+1][2])
    A_rot[ie+1][1] = A_rot[ie+1][1].GetAvg()
    A_rot[ie+1][2] = A_rot[ie+1][2].GetAvg()

# FFs_p_Arot,FFs_p_Arot_chi2 = np.linalg.lstsq(np.array(A_rot),np.array(b_p))[0:2]
FFs_p_Arot,covar,FFs_p_Arot_chi2 = fit_data.LSFit(np.array(A_rot).swapaxes(0,1),np.array(b_p),np.array(b_p_err))
FFs_p_Arot = np.append(FFs_p_Arot,[FFs_p_Arot[-1]*a/(2*a_mN)])

# FFs_n_Arot,FFs_n_Arot_chi2 = np.linalg.lstsq(np.array(A_rot),np.array(b_n))[0:2]
FFs_n_Arot,covar,FFs_n_Arot_chi2 = fit_data.LSFit(np.array(A_rot).swapaxes(0,1),np.array(b_n),np.array(b_n_err))
FFs_n_Arot = np.append(FFs_n_Arot,[FFs_n_Arot[-1]*a/(2*a_mN)])

A_vec_rot = [ia[:-1] for ia in A_rot[:1]] + [ia[:-1] for ia in A_rot[-2:]]
# FFs_p_vec_Arot,FFs_p_vec_Arot_chi2 = np.linalg.lstsq(np.array(A_vec_rot),np.array(b_p_vec))[0:2]
FFs_p_vec_Arot,covar,FFs_p_vec_Arot_chi2 = fit_data_vec.LSFit(np.array(A_vec_rot).swapaxes(0,1),np.array(b_p_vec),np.array(b_p_vec_err))
# FFs_n_vec_Arot,FFs_n_vec_Arot_chi2 = np.linalg.lstsq(np.array(A_vec_rot),np.array(b_n_vec))[0:2]
FFs_n_vec_Arot,covar,FFs_n_vec_Arot_chi2 = fit_data_vec.LSFit(np.array(A_vec_rot).swapaxes(0,1),np.array(b_n_vec),np.array(b_n_vec_err))

FF3_p_vec_Arot = [SolveF3_proton(FFs_p_vec_Arot[0],FFs_p_vec_Arot[1],Arot=True),
             SolveF3_proton(FFs_p_vec_Arot[0],FFs_p_vec_Arot[1],Arot=True)*a/(2*a_mN)]
FF3_n_vec_Arot = [SolveF3_neutron(FFs_n_vec_Arot[0],FFs_n_vec_Arot[1],Arot=True),
             SolveF3_neutron(FFs_n_vec_Arot[0],FFs_n_vec_Arot[1],Arot=True)*a/(2*a_mN)]
FFs_p_vec_Arot = np.append(FFs_p_vec_Arot,FF3_p_vec_Arot)
FFs_n_vec_Arot = np.append(FFs_n_vec_Arot,FF3_n_vec_Arot)

def percent_diff(val1,val2,fmt=True):
    out = np.abs((val1-val2)/val1)
    if isinstance(out,BootStrap):
        out.Stats()
        out = out.Avg
    if fmt:
        return f'{out:.7f}'
    else:
        return out

GeGm_n_direct = [0,0]
GeGm_n_direct[0] =  MakeValAndErr(  0.0089295673,  0.0003406836,Dec=2)
GeGm_n_direct[1] = MakeValAndErr(  -1.1734138556,  0.0054696330,Dec=2)


GeGm_p_direct = [0,0]
GeGm_p_direct[0] =  MakeValAndErr(  0.7825871135, 0.0006027896,Dec=2)
GeGm_p_direct[1] = MakeValAndErr(   1.9168725391, 0.0090745190,Dec=2)


print('Ratio of two equations with F3 in them: ')
print('neutron: ' , equ_rat,MakeValAndErr(R_rat_n,R_rat_n_err,Dec=2))
print('proton:  ' , equ_rat,MakeValAndErr(R_rat_p,R_rat_p_err,Dec=2))



FFs_p = AppendGeGm(FFs_p*renorm)
FFs_p_vec = AppendGeGm(FFs_p_vec*renorm)
FFs_n = AppendGeGm(FFs_n*renorm)
FFs_n_vec = AppendGeGm(FFs_n_vec*renorm)

print()
print('                         F1             F2              GE              Gm              F3              F3/2mn [efm]')
print('     Proton:')
# print('F_p old        ',FFs_p_old*renorm)
print('F_p Arot file    ',fmt_boot(FFs_p))
# print('F_p old rot    ',FFrot_p_old*renorm)
# print('F_p rot        ',fmt_boot(FFrot_p*renorm))
# print('F_p Arot       ',fmt_boot(FFs_p_Arot*renorm))
print('F_p file         ',fmt_boot(FileFF_p_str,doubbrak=True))
print('--------')
print('F_p vec          ',fmt_boot(FFs_p_vec[:4]))
print('F_p vec file     ',fmt_boot(FileFF_p_str_noT,doubbrak=True))
print('GeGm_p_direct      ',fmt_boot([0.,0.]+GeGm_p_direct,doubbrak=True))
# print('F_p rot vec      ',fmt_boot(FFrot_p_vec*renorm))
# print('F_p Arot vec     ',fmt_boot(FFs_p_vec_Arot*renorm))
# print(' rot,Arot % diff ',fmt_boot(percent_diff(FFrot_p[3],FFs_p_Arot[3])))
# print('F_p rotF3t       ',FFrot_p_F2t*renorm)
# print('F_p rotF3t vec   ',FFrot_p_vec_F2t*renorm)
print('     Neutron:    ')
# print('F_n old          ',FFs_n_old*renorm)
# print('F_n norot file   ',fmt_boot(FileFF_N_str_norot))
print('F_n Arot file    ',fmt_boot(FFs_n))
# print('F_n old rot      ',FFrot_n_old*renorm)
# print('F_n rot          ',fmt_boot(FFrot_n*renorm))
# print('F_n rotF3t       ',FFrot_n_F2t*renorm)
# print('F_n rotF3t vec   ',FFrot_n_vec_F2t*renorm)
# print('F_N Arot         ',fmt_boot(FFs_n_Arot*renorm))
print('F_n file         ',fmt_boot(FileFF_N_str,doubbrak=True))
print('--------')
print('F_n vec          ',fmt_boot(FFs_n_vec[:4]))
print('F_n vec File     ',fmt_boot(FileFF_N_str_noT,doubbrak=True))
print('GeGm_n_direct      ',fmt_boot([0.,0.]+GeGm_n_direct,doubbrak=True))

print('---------------')
print('A:')
for ia in A:
    print(fmt_boot(ia))
print()


print('A_rot:')
for ia in A_rot:
    print(ia)
print('A Arot diff:')
for ia,iar in zip(A,A_rot):
    print([fmt_boot(np.abs((ja-jarot)/ja)) for ja,jarot in zip(ia,iar) if np.abs(ja) > 0])
print()
print()
print('     Chi2s Proton: ')
print('F_p              ',FFs_p_chi2)
# print('F_p old rot      ',FFrot_p_old*renorm)
print('F_p Arot         ',FFs_p_Arot_chi2)
print('--------')
print('F_p vec          ',FFs_p_vec_chi2)
print('F_p Arot vec     ',FFs_p_vec_Arot_chi2)

print('     Chi2s Neutron: ')
print('F_n              ',FFs_n_chi2)
# print('F_n old rot      ',FFrot_n_old*renorm)
print('F_n Arot         ',FFs_n_Arot_chi2)
print('--------')
print('F_n vec          ',FFs_n_vec_chi2)
print('F_n Arot vec     ',FFs_n_vec_Arot_chi2)




using_solution_in_file = '''

Using ratio function results displayed in the output xml file (up to 5 sig fig.)

mpi = 570 MeV, Q^2 = 1 unit.

     Proton:
F_p               [ 0.72374922      -2.64064861     -0.64721004     -0.04524584]        No rotation at all, 6 equations
----------------------------------------------------------------------------------------
F_p rot           [0.723758855, '   -2.689(22)', '   0.38(14)', '    0.026(10)']        Post solving rotation (WRONG)
F_p rot file      ['0.72377((63))', -2.694((34))', ' 0.380((60))', ' 0.0264((42))']     Above result from full analsis code
----------------------------------------------------------------------------------------
F_p Arot          [ 0.72375894,     -2.6406779,      0.379685932, '  0.026466(50)']     Rotating Coeffs taken from F_p
F_p Arot file     [ 0.72375892,     -2.6406788,      0.376476541, '  0.026242(49)']     Rotating Coeffs in code (using bs samples)
F_p file          ['0.72377((63))', -2.6407((92))', '0.388((59))', ' 0.0270((41))']     Final result from correct rotation
----------------------------------------------------------------------------------------
F_p vec           [ 0.723760335,     -2.6406153  ]   Only solving for form factors F1 and F2
F_p vec File      ['0.72377((63))', '-2.6406((92))']   Above result from full analysis code


     Neutron:
F_n              [  0.0343115601,    1.1391135,       0.32282021, '   0.022502(42)']     No rotation at all, 6 equations
F_n norot file   [' 0.03430((35))', '1.1392((54))', ' 0.307((73))', ' 0.0213((50))']     Above result from full analsis
----------------------------------------------------------------------------------------
F_n rot           [ 0.034311561, '   1.1760(74)', '  -0.122(63)', '  -0.0085(44)']       Post solving rotation (WRONG)
F_n rot file      ['0.03430((35))', '1.173((23))', ' -0.140((35))', '-0.0097((24))']     Above result from full analsis code
----------------------------------------------------------------------------------------
F_N Arot          [0.0343115160,     1.13911150,     -0.12315946, '  -0.00858(2)']       Rotating Coeffs taken from F_n
F_n Arot file     [ 0.034311819,     1.13912510,     -0.13431995, '  -0.00936(2)']       Rotating Coeffs in code (using bs samples)
F_n file          ['0.03430((35))', '1.1392((54))', '-0.135((32))', '-0.0094((22))']     Final result from correct rotation
----------------------------------------------------------------------------------------
F_n vec           [ 0.034311,        1.13907861 ]      Only solving for form factors F1 and F2
F_n vec File      ['0.03430((35))',  1.1391((54))']    Above result from full analysis code

'''



UsingSameSolution = '''

Using the same ratio function results in all the files (up to 5 sig fig.)

mpi = 570 MeV, Q^2 = 1 unit.

     Proton:
F_p               [ 0.72374922      -2.64064861     -0.64721004     -0.04524584]        No rotation at all, 6 equations
----------------------------------------------------------------------------------------
F_p rot           [ 0.72375884, '   -2.690(22)', '   0.37(14)', '    0.026(10)']        Post solving rotation (WRONG)
F_p rot file      ['0.72377((63))', -2.694((34))', ' 0.380((60))', ' 0.0264((42))']     Above result from full analsis code
----------------------------------------------------------------------------------------
F_p Arot          [ 0.72375892,     -2.6406788,      0.376476541, '  0.026242(49)']     Rotating Coeffs taken from F_p
F_p Arot file     [ 0.72375892,     -2.6406788,      0.376476541, '  0.026242(49)']     Rotating Coeffs in code (using bs samples)
F_p file          ['0.72377((63))', -2.6407((92))', '0.388((59))', ' 0.0270((41))']     Final result from correct rotation
----------------------------------------------------------------------------------------
F_p vec           [ 0.723760335,     -2.6406153  ]   Only solving for form factors F1 and F2
F_p vec File      ['0.72377((63))', '-2.6406((92))']   Above result from full analysis code


     Neutron:
F_n              [  0.034311884,     1.13912802,      0.3111573, '    0.021689(41)']     No rotation at all, 6 equations
F_n norot file   [' 0.03430((35))', '1.1392((54))', ' 0.307((73))', ' 0.0213((50))']     Above result from full analsis
----------------------------------------------------------------------------------------
F_n rot           [ 0.0343118836, '  1.1717(79)', '  -0.133(62)', '  -0.0093(44)']       Post solving rotation (WRONG)
F_n rot file      ['0.03430((35))', '1.173((23))', ' -0.140((35))', '-0.0097((24))']     Above result from full analsis code
----------------------------------------------------------------------------------------
F_N Arot          [ 0.034311819,     1.13912510,     -0.13431995, '  -0.00936(2)']       Rotating Coeffs taken from F_n
F_n Arot file     [ 0.034311819,     1.13912510,     -0.13431995, '  -0.00936(2)']       Rotating Coeffs in code (using bs samples)
F_n file          ['0.03430((35))', '1.1392((54))', '-0.135((32))', '-0.0094((22))']     Final result from correct rotation
----------------------------------------------------------------------------------------
F_n vec           [ 0.034311,        1.13907861 ]    Only solving for form factors F1 and F2
F_n vec File      ['0.03430((35))', '1.1391((54))']  Above result from full analysis code

'''
