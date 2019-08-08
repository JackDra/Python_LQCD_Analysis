import numpy as np
from XmlFormatting import MakeValAndErr

from copy import deepcopy
alpha = -0.1898009109
a_mN = 0.6487
a = 0.0907
renorm = 0.7354
RQ_full_list = [2,3]
RQ_solve_list = [2,3]
print('------New Results, Rimp, sqrt(8t_f)= 0.6 fm------')
print(' Solving with RQ solve list: ',RQ_solve_list)
exclude_list = [i for i in RQ_full_list if i not in RQ_solve_list]

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


## from file /mnt/scratch/dragosja/data/resultsFixedG2/RC32x64Kud01372700Ks01364000/FormFactors/Proton_qmax0_ppmax4/VectorTopBNL_t_f5.47_Rtsumfitr7-32_Afitr10-20_RFmin3_RChi50-100.xml
## Q2 = 1

A = [
    [-0.14603    ,-0.14603  ,0.00000e+00 ],
    [-0.00409    ,-0.00409  ,6.59947e-18 ],
    [-0.00409    ,-0.00827  , 0.02204    ],
    [0.02772     , 0.05605  ,-0.14928    ],
    [-0.14603    , 0.00325  ,0.00000e+00 ],
    [0.98928     ,-0.02204  ,0.00000e+00 ]
]

A_np = np.array(A)
equ_rat = A_np[2]/A_np[3]


for ie in exclude_list: del A[ie]

b_p = [
    0.38063,
    0.03955,
    0.01546,
    -0.04085,
    -0.16227,
    1.05168
]

b_p_err = [
    0.00149,
    0.01154,
    0.01022,
    0.01162,
    0.00080,
    0.00082
]

## F1, F2tilde, F3tilde,F3tild/2mN
File_FF_p_list= [
0.7237674820,
-2.6941844312 ,
0.3800439564,
0.0264021537 ]

File_FF_p_err = [
    0.0006255002,
    0.0337083441,
    0.0599716766,
    0.0041641810
]
FileFF_p_str = [MakeValAndErr(iff,ifferr,Dec=2) for iff,ifferr in zip(File_FF_p_list,File_FF_p_err)]

## from file /mnt/scratch/dragosja/data/resultsFixedG2/RC32x64Kud01372700Ks01364000/FormFactors/Neutron_qmax0_ppmax4/VectorTopBNL_t_f6.01_Rtsumfitr7-32_Afitr10-20_RFmin3_RChi50-100.xml
## Q2 = 1

b_n = [
    -0.23300,
    -0.02710,
    -0.01141,
    0.02072,
    0.00561,
    0.01299
]

b_n_err = [
    0.00093,
    0.00839,
    0.00611,
    0.00763,
    0.00051,
    0.00048
]

## F1, F2tilde, F3tilde,F3tild/2mN
File_FF_n_list= [
 0.0343035638,
 1.1725807648,
-0.1397728093,
-0.0097106977
]

File_FF_n_err = [
0.0003524230,
0.0225931607,
0.0351592296,
0.0024441459
]

FileFF_n_str = [MakeValAndErr(iff,ifferr,Dec=2) for iff,ifferr in zip(File_FF_n_list,File_FF_n_err)]



R_rat_n = b_n[2]/b_n[3]
R_rat_p = b_p[2]/b_p[3]
R_rat_n_err = np.abs(R_rat_n)*(np.abs(b_n_err[2]/b_n[2]) + np.abs(b_n_err[3]/b_n[3]))
R_rat_p_err = np.abs(R_rat_p)*(np.abs(b_p_err[2]/b_p[2]) + np.abs(b_p_err[3]/b_p[3]))

for ie in exclude_list: del b_p[ie]
for ie in exclude_list: del b_n[ie]

FFs_p = np.linalg.lstsq(np.array(A),np.array(b_p))[0]
FFs_p = np.append(FFs_p,[FFs_p[-1]*a/(2*a_mN)])
FFrot_p = alpha_rot(FFs_p)
FFrot_p = np.append(FFrot_p,FFrot_p[-1]*a/(2*a_mN))

FFrot_p_F2t = alpha_rot_F2tilde(FFs_p)
FFrot_p_F2t = np.append(FFrot_p_F2t,FFrot_p_F2t[-1]*a/(2*a_mN))

# print(FFs_p*renorm)

b_p_old = [-ib for ib in b_p[:4]] + b_p[4:]


FFs_p_old = np.linalg.lstsq(np.array(A),np.array(b_p_old))[0]
FFs_p_old = np.append(FFs_p_old,[FFs_p_old[-1]*a/(2*a_mN)])
FFrot_p_old = alpha_rot(FFs_p_old)
FFrot_p_old = np.append(FFrot_p_old,FFrot_p_old[-1]*a/(2*a_mN))

FFrot_p_old_F2t = alpha_rot_F2tilde(FFs_p_old)
FFrot_p_old_F2t = np.append(FFrot_p_old_F2t,FFrot_p_old_F2t[-1]*a/(2*a_mN))

# print(FFs_p_old*renorm)


FFs_n = np.linalg.lstsq(np.array(A),np.array(b_n))[0]
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
A_vec = [ia[:-1] for ia in A[:2]] + [ia[:-1] for ia in A[-2:]]
b_p_vec = b_p[:2] + b_p[-2:]
b_n_vec = b_n[:2] + b_n[-2:]
FFs_p_vec = np.linalg.lstsq(np.array(A_vec),np.array(b_p_vec))[0]
FFs_n_vec = np.linalg.lstsq(np.array(A_vec),np.array(b_n_vec))[0]

def SolveF3_proton(F1,F2,eq_list=RQ_solve_list,Arot=False):
    out_list = []
    for eq in range(len(eq_list)):
        if Arot:
            C1,C2,C3 = tuple(A_rot[eq+2])
        else:
            C1,C2,C3 = tuple(A[eq+2])
        R = b_p[eq+2]
        out_list.append((R - C1*F1 - C2*F2)/C3)
    return np.mean(out_list)

def SolveF3_neutron(F1,F2,eq_list=RQ_solve_list,Arot=False):
    out_list = []
    for eq in range(len(eq_list)):
        if Arot:
            C1,C2,C3 = tuple(A_rot[eq+2])
        else:
            C1,C2,C3 = tuple(A[eq+2])
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
    A_rot[ie+2][1],A_rot[ie+2][2] = (   np.cos(2*alpha)*A_rot[ie+2][1] - np.sin(2*alpha)*A_rot[ie+2][2],
                                        np.sin(2*alpha)*A_rot[ie+2][1] + np.cos(2*alpha)*A_rot[ie+2][2])

FFs_p_Arot = np.linalg.lstsq(np.array(A_rot),np.array(b_p))[0]
FFs_p_Arot = np.append(FFs_p_Arot,[FFs_p_Arot[-1]*a/(2*a_mN)])

FFs_n_Arot = np.linalg.lstsq(np.array(A_rot),np.array(b_n))[0]
FFs_n_Arot = np.append(FFs_n_Arot,[FFs_n_Arot[-1]*a/(2*a_mN)])

A_vec_rot = [ia[:-1] for ia in A_rot[:2]] + [ia[:-1] for ia in A_rot[-2:]]
FFs_p_vec_Arot = np.linalg.lstsq(np.array(A_vec_rot),np.array(b_p_vec))[0]
FFs_n_vec_Arot = np.linalg.lstsq(np.array(A_vec_rot),np.array(b_n_vec))[0]

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
print('F_p file         ',FileFF_p_str)
print('F_p Arot         ',FFs_p_Arot*renorm)
print('--------')
print('F_p vec          ',FFs_p_vec*renorm)
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
print('F_n file         ',FileFF_n_str)
# print('F_n rotF3t       ',FFrot_n_F2t*renorm)
# print('F_n rotF3t vec   ',FFrot_n_vec_F2t*renorm)
print('F_N Arot         ',FFs_n_Arot*renorm)
print('--------')
print('F_n vec          ',FFs_n_vec*renorm)
print('F_n rot vec      ',FFrot_n_vec*renorm)
print('F_n Arotvec      ',FFs_n_vec_Arot*renorm)
print(' rot,Arot % diff ',percent_diff(FFrot_n[3],FFs_n_Arot[3]))
