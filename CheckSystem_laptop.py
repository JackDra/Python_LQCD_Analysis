import numpy as np
from XmlFormatting import MakeValAndErr
from FitFunctions import Fitting

from copy import deepcopy
alpha = 0.5042365099
a_mN =  0.4817554883
a = 0.0907
renorm = 0.7354
RQ_full_list = [1,2,3]
RQ_solve_list = [1,2,3]
print('------New Results, Rimp, sqrt(8t_f)= 0.6 fm------')
print(' Solving with RQ solve list: ',RQ_solve_list)
exclude_list = [i for i in RQ_full_list if i not in RQ_solve_list]
def lin_fun(x,y):
    return sum([ix*ip for ix,ip in zip(x,y)])
fit_data = Fitting(Funs=(lin_fun,3))

def alpha_rot(F_list):
    c2a = np.cos(2*alpha)
    s2a = np.sin(2*alpha)
    return np.array([F_list[0],
                    c2a*F_list[1] - s2a*F_list[2],
                    s2a*F_list[1] + c2a*F_list[2]])

def alpha_rot_F2tilde(F_list):
    c2a = np.cos(2*alpha)
    s2a = np.sin(2*alpha)
    t2a = np.tan(2*alpha)
    return np.array([F_list[0],F_list[1],
                     (c2a+t2a*s2a)*F_list[2] + t2a*F_list[1]])


## from file /home/jackdra/LQCD/Results/DebugSystemV2/RC32x64Kud01375400Ks01364000/FormFactors/Neutron_qmax0_ppmax4/VectorTop_t_f5.47_Afitr10-20_RFmin3_RChi50-100.xml
## Q2 = 1

opp_list = [
'P3_g1_q0-10',
'P3_g3_Top_cmplx_q-100',
'P3_g3_Top_cmplx_q00-1',
'P3_g4_Top_q00-1',
'P4_g1_cmplx_q-100',
'P4_g4_q-100',
]


A = [
[-0.19231 ,-0.19231 ,0.00000e+00],
[ 0.01894 , 0.00831 , 0.01617   ],
[ 0.01894 ,-0.01467 , 0.05306   ],
[-0.09682 , 0.07513 ,-0.27091   ],
[-0.19231 , 0.00768 ,0.00000e+00],
[ 0.98133 ,-0.03919 ,0.00000e+00]
]
A_err = [
[0.00144, 0.00144, 0.00000e+00],
[0.00764, 0.00275, 0.01037    ],
[0.00764, 0.01378, 0.00707    ],
[0.03976, 0.07092, 0.03820    ],
[0.00144, 0.00018, 0.00000e+00],
[0.00028, 0.00062, 0.00000e+00]
]

b_n = [
-0.28434,
-0.34440,
-0.10205,
-0.72106,
-0.00987,
 0.03459]

b_n_err = [
0.00157,
0.10502,
0.07300,
0.06897,
0.00275,
0.00730]


## from file /home/jackdra/LQCD/Results/DebugSystemV2/RC32x64Kud01375400Ks01364000/FormFactors/Proton_qmax0_ppmax4/Pickle/VectorTop_t_f5.47_Afitr10-20_RFmin3_RChi50-100.p.Ratios
## Q2 = 1

b_p = [
0.457624  ,
-0.325935 ,
-0.753022 ,
-0.010671 ,
-0.183223 ,
1.056421
]

b_p_err = [
0.017612,
0.077270,
0.174658,
0.197883,
0.024995,
0.002661
]

## F1, F2tilde, F3tilde,F3tild/2mN
File_FF_N_list= [
0.0691131793   ,
1.0181905442    ,
1.4830612318   ,
0.1167313888   ]

File_FF_N_err = [
    0.0014611356,
    0.0006774452,
    0.1563818539,
    0.0120707698
]
FileFF_N_str = [MakeValAndErr(iff,ifferr,Dec=2) for iff,ifferr in zip(File_FF_N_list,File_FF_N_err)]

## F1, F2 no topological charge
File_FF_N_list_noT= [
0.0691129602   ,
1.0182155282   ]

File_FF_N_err_noT = [
    0.0014601343,
    0.0006922578
]
FileFF_N_str_noT = [MakeValAndErr(iff,ifferr,Dec=2) for iff,ifferr in zip(File_FF_N_list_noT,File_FF_N_err_noT)]


## F1, F2 no topological charge
File_FF_p_list_noT= [
 0.6940335750   ,
-2.4443253059   ]

File_FF_p_err_noT = [
    0.0036619942,
    0.0768310124
]
FileFF_p_str_noT = [MakeValAndErr(iff,ifferr,Dec=2) for iff,ifferr in zip(File_FF_p_list_noT,File_FF_p_err_noT)]



A_np = np.array(A)
equ_rat = A_np[2]/A_np[3]


for ie in exclude_list: del A[ie]



R_rat_n = b_n[2]/b_n[3]
R_rat_p = b_p[2]/b_p[3]
R_rat_n_err = np.abs(R_rat_n)*(np.abs(b_n_err[2]/b_n[2]) + np.abs(b_n_err[3]/b_n[3]))
R_rat_p_err = np.abs(R_rat_p)*(np.abs(b_p_err[2]/b_p[2]) + np.abs(b_p_err[3]/b_p[3]))

for ie in exclude_list: del b_p[ie]
for ie in exclude_list: del b_n[ie]

FFs_p,covar,FFs_p_chi2 = fit_data.LSFit(np.array(A).swapaxes(0,1),np.array(b_p),np.array(b_p_err))

FFs_p = np.append(FFs_p,[FFs_p[-1]*a/(2*a_mN)])
FFrot_p = alpha_rot(FFs_p)
FFrot_p = np.append(FFrot_p,FFrot_p[-1]*a/(2*a_mN))

FFrot_p_F2t = alpha_rot_F2tilde(FFs_p)
FFrot_p_F2t = np.append(FFrot_p_F2t,FFrot_p_F2t[-1]*a/(2*a_mN))

# print(FFs_p*renorm)

b_p_old = [-ib for ib in b_p[:4]] + b_p[4:]


# FFs_p_old = np.linalg.lstsq(np.array(A),np.array(b_p_old))[0]
FFs_p_old = fit_data.LSFit(np.array(A).swapaxes(0,1),np.array(b_p_old),np.array(b_p_err))[0]
FFs_p_old = np.append(FFs_p_old,[FFs_p_old[-1]*a/(2*a_mN)])
FFrot_p_old = alpha_rot(FFs_p_old)
FFrot_p_old = np.append(FFrot_p_old,FFrot_p_old[-1]*a/(2*a_mN))

FFrot_p_old_F2t = alpha_rot_F2tilde(FFs_p_old)
FFrot_p_old_F2t = np.append(FFrot_p_old_F2t,FFrot_p_old_F2t[-1]*a/(2*a_mN))

# print(FFs_p_old*renorm)


# FFs_n,FFs_n_chi2 = np.linalg.lstsq(np.array(A),np.array(b_n))[0:2]
FFs_n,covar,FFs_n_chi2 = fit_data.LSFit(np.array(A).swapaxes(0,1),np.array(b_n),np.array(b_n_err))
FFs_n = np.append(FFs_n,[FFs_n[-1]*a/(2*a_mN)])
FFrot_n = alpha_rot(FFs_n)
FFrot_n = np.append(FFrot_n,FFrot_n[-1]*a/(2*a_mN))
FFrot_n_F2t = alpha_rot_F2tilde(FFs_n)
FFrot_n_F2t = np.append(FFrot_n_F2t,FFrot_n_F2t[-1]*a/(2*a_mN))


b_n_old = [-ib for ib in b_n[:4]] + b_n[4:]
FFs_n_old = np.linalg.lstsq(np.array(A),np.array(b_n_old))[0]
FFs_n_old = np.append(FFs_n_old,[FFs_n_old[-1]*a/(2*a_mN)])
FFrot_n_old = alpha_rot(FFs_n_old)
FFrot_n_old = np.append(FFrot_n_old,FFrot_n_old[-1]*a/(2*a_mN))

# A_vec = [ia[:-1] for ia in A[:2]] + [ia[:-1] for ia in A[4:]]
# b_p_vec = b_p[:2] + b_p[4:]
# b_n_vec = b_n[:2] + b_n[4:]
A_vec = [ia[:-1] for ia in A[:1]] + [ia[:-1] for ia in A[-2:]]
b_p_vec = b_p[:1] + b_p[-2:]
b_n_vec = b_n[:1] + b_n[-2:]
b_p_vec_err = b_p_err[:1] + b_p_err[-2:]
b_n_vec_err = b_n_err[:1] + b_n_err[-2:]
# FFs_p_vec,FFs_p_vec_chi2 = np.linalg.lstsq(np.array(A_vec),np.array(b_p_vec))[0:2]
FFs_p_vec,covar,FFs_p_vec_chi2 = fit_data.LSFit(np.array(A_vec).swapaxes(0,1),np.array(b_p_vec),np.array(b_p_vec_err))
# FFs_n_vec,FFs_n_vec_chi2 = np.linalg.lstsq(np.array(A_vec),np.array(b_n_vec))[0:2]
FFs_n_vec,covar,FFs_n_vec_chi2 = fit_data.LSFit(np.array(A_vec).swapaxes(0,1),np.array(b_n_vec),np.array(b_n_vec_err))

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


A_rot = deepcopy(A)
for ie in range(len(RQ_solve_list)):
    A_rot[ie+1][1],A_rot[ie+1][2] = (   np.cos(2*alpha)*A_rot[ie+1][1] - np.sin(2*alpha)*A_rot[ie+1][2],
                                        np.sin(2*alpha)*A_rot[ie+1][1] + np.cos(2*alpha)*A_rot[ie+1][2])

# FFs_p_Arot,FFs_p_Arot_chi2 = np.linalg.lstsq(np.array(A_rot),np.array(b_p))[0:2]
FFs_p_Arot,covar,FFs_p_Arot_chi2 = fit_data.LSFit(np.array(A_rot).swapaxes(0,1),np.array(b_p),np.array(b_p_err))
FFs_p_Arot = np.append(FFs_p_Arot,[FFs_p_Arot[-1]*a/(2*a_mN)])

# FFs_n_Arot,FFs_n_Arot_chi2 = np.linalg.lstsq(np.array(A_rot),np.array(b_n))[0:2]
FFs_n_Arot,covar,FFs_n_Arot_chi2 = fit_data.LSFit(np.array(A_rot).swapaxes(0,1),np.array(b_n),np.array(b_n_err))
FFs_n_Arot = np.append(FFs_n_Arot,[FFs_n_Arot[-1]*a/(2*a_mN)])

A_vec_rot = [ia[:-1] for ia in A_rot[:1]] + [ia[:-1] for ia in A_rot[-2:]]
# FFs_p_vec_Arot,FFs_p_vec_Arot_chi2 = np.linalg.lstsq(np.array(A_vec_rot),np.array(b_p_vec))[0:2]
FFs_p_vec_Arot,covar,FFs_p_vec_Arot_chi2 = fit_data.LSFit(np.array(A_vec_rot).swapaxes(0,1),np.array(b_p_vec),np.array(b_p_vec_err))
# FFs_n_vec_Arot,FFs_n_vec_Arot_chi2 = np.linalg.lstsq(np.array(A_vec_rot),np.array(b_n_vec))[0:2]
FFs_n_vec_Arot,covar,FFs_n_vec_Arot_chi2 = fit_data.LSFit(np.array(A_vec_rot).swapaxes(0,1),np.array(b_n_vec),np.array(b_n_vec_err))

FF3_p_vec_Arot = [SolveF3_proton(FFs_p_vec_Arot[0],FFs_p_vec_Arot[1],Arot=True),
             SolveF3_proton(FFs_p_vec_Arot[0],FFs_p_vec_Arot[1],Arot=True)*a/(2*a_mN)]
FF3_n_vec_Arot = [SolveF3_neutron(FFs_n_vec_Arot[0],FFs_n_vec_Arot[1],Arot=True),
             SolveF3_neutron(FFs_n_vec_Arot[0],FFs_n_vec_Arot[1],Arot=True)*a/(2*a_mN)]
FFs_p_vec_Arot = np.append(FFs_p_vec_Arot,FF3_p_vec_Arot)
FFs_n_vec_Arot = np.append(FFs_n_vec_Arot,FF3_n_vec_Arot)

def percent_diff(val1,val2,fmt=True):
    out = np.abs((val1-val2)/val1)
    if fmt:
        return f'{out:.7f}'
    else:
        return out


print('Ratio of two equations with F3 in them: ')
print('neutron: ' , equ_rat,MakeValAndErr(R_rat_n,R_rat_n_err,Dec=2))
print('proton:  ' , equ_rat,MakeValAndErr(R_rat_p,R_rat_p_err,Dec=2))

print()
print('     Proton:')
# print('F_p old          ',FFs_p_old*renorm)
print('F_p              ',FFs_p*renorm)
# print('F_p old rot      ',FFrot_p_old*renorm)
print('F_p rot          ',FFrot_p*renorm)
# print('F_p file         ',FileFF_p_str) DONT HAVE!
print('F_p Arot         ',FFs_p_Arot*renorm)
print('--------')
print('F_p vec          ',FFs_p_vec*renorm)
print('F_p vec File     ',FileFF_p_str_noT)
print('F_p rot vec      ',FFrot_p_vec*renorm)
print('F_p Arot vec     ',FFs_p_vec_Arot*renorm)
print(' rot,Arot % diff ',percent_diff(FFrot_p[3],FFs_p_Arot[3]))
# print('F_p rotF3t       ',FFrot_p_F2t*renorm)
# print('F_p rotF3t vec   ',FFrot_p_vec_F2t*renorm)
print('     Neutron:    ')
# print('F_n old          ',FFs_n_old*renorm)
print('F_n              ',FFs_n*renorm)
# print('F_n old rot      ',FFrot_n_old*renorm)
print('F_n rot          ',FFrot_n*renorm)
print('F_n file         ',FileFF_N_str)
# print('F_n rotF3t       ',FFrot_n_F2t*renorm)
# print('F_n rotF3t vec   ',FFrot_n_vec_F2t*renorm)
print('F_N Arot         ',FFs_n_Arot*renorm)
print('--------')
print('F_n vec          ',FFs_n_vec*renorm)
print('F_n vec File     ',FileFF_N_str_noT)
print('F_n rot vec      ',FFrot_n_vec*renorm)
print('F_n Arotvec      ',FFs_n_vec_Arot*renorm)
print(' rot,Arot % diff ',percent_diff(FFrot_n[3],FFs_n_Arot[3]))

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
