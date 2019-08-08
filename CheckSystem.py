import numpy as np
from XmlFormatting import MakeValAndErr

from copy import deepcopy
alpha = -0.1898009109
a_mN = 0.6487
a = 0.0907
renorm = 0.7354
RQ_full_list = [2,3]
RQ_solve_list = [3]
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


A =[[-0.14603  ,-0.14603 , 0.00000e+00  ],
    [-0.00409  ,-0.00409 , 6.59947e-18  ],
    [-0.00409  ,-0.00827 ,  0.02204     ],
    [ 0.02772  , 0.05605 , -0.14928     ],
    [-0.14603  , 0.00325 , 0.00000e+00  ],
    [ 0.98928  ,-0.02204 , 0.00000e+00  ]]
A_np = np.array(A)
equ_rat = A_np[2]/A_np[3]


for ie in exclude_list: del A[ie]

## from file /mnt/scratch/dragosja/data/resultsFixedG2/RC32x64Kud01372700Ks01364000/FormFactors/Proton_qmax0_ppmax4/VectorTopBNL_t_f5.47_Afitr10-20_RFmin3_RChi50-100.xml
## Q2 = 1


b_p = [  0.38063,
         0.031760,
         0.01396,
        -0.03675,
        -0.16227,
         1.05168]

b_p_err = [0.00149,
           0.01963,
           0.01637,
           0.02589,
           0.00080,
           0.00082]

## from file /mnt/scratch/dragosja/data/resultsFixedG2/RC32x64Kud01372700Ks01364000/FormFactors/Neutron_qmax0_ppmax4/VectorTopBNL_t_f5.47_Afitr10-20_RFmin3_RChi50-100.xml
## Q2 = 1

b_n = [   -0.23300,
        -0.02112,
        -0.01898,
         0.01094,
         0.00561,
         0.01299]
b_n_err = [ 0.00093,
            0.01546,
            0.01125,
            0.01732,
            0.00051,
            0.00048]


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


print('Ratio of two equations with F3 in them: ')
print('neutron: ' , equ_rat,MakeValAndErr(R_rat_n,R_rat_n_err,Dec=2))
print('proton:  ' , equ_rat,MakeValAndErr(R_rat_p,R_rat_p_err,Dec=2))

print()
print('     Proton:')
# print('F_p old          ',FFs_p_old*renorm)
print('F_p              ',FFs_p*renorm)
# print('F_p old rot      ',FFrot_p_old*renorm)
print('F_p rot          ',FFrot_p*renorm)
print('F_p Arot         ',FFs_p_Arot*renorm)
print('--------')
print('F_p vec          ',FFs_p_vec*renorm)
print('F_p rot vec      ',FFrot_p_vec*renorm)
print('F_p Arot vec     ',FFs_p_vec_Arot*renorm)
print(' % diff          ',np.abs((FFrot_p[3]-FFs_p_Arot[3])/FFrot_p[3]))
# print('F_p rotF3t       ',FFrot_p_F2t*renorm)
# print('F_p rotF3t vec   ',FFrot_p_vec_F2t*renorm)
print('     Neutron:    ')
# print('F_n old          ',FFs_n_old*renorm)
print('F_n              ',FFs_n*renorm)
# print('F_n old rot      ',FFrot_n_old*renorm)
print('F_n rot          ',FFrot_n*renorm)
# print('F_n rotF3t       ',FFrot_n_F2t*renorm)
# print('F_n rotF3t vec   ',FFrot_n_vec_F2t*renorm)
print('F_N Arot         ',FFs_n_Arot*renorm)
print('--------')
print('F_n vec          ',FFs_n_vec*renorm)
print('F_n rot vec      ',FFrot_n_vec*renorm)
print('F_n Arotvec      ',FFs_n_vec_Arot*renorm)
print(' % diff          ',np.abs((FFrot_n[3]-FFs_n_Arot[3])/FFrot_n[3]))


Full_System = '''

Ratio of two equations with F3 in them:
neutron:  [-0.1475469  -0.14754683 -0.14764202] -1.7(3.8)
proton:   [-0.1475469  -0.14754683 -0.14764202] -0.38(71)

     Proton:
F_p               [ 0.72364926 -2.64345484 -0.6710312  -0.04691115]
F_p rot           [ 0.72364926 -2.7039243   0.35627195  0.02490663]
F_p Arot          [ 0.72364926 -2.64345493  0.35503624  0.02482025]
--------
F_p vec           [ 0.72364924 -2.64345634 -0.53446375 -0.03736385]
F_p rot vec       [ 0.72364924 -2.65332054  0.48311803  0.03377432]
F_p Arot vec      [ 0.72364924 -2.64345634  0.48292114  0.03376056]
 % diff           0.0034684331426729333
     Neutron:
F_n               [ 0.03427114  1.14115847  0.36857303  0.02576659]
F_n rot           [ 0.03427114  1.19649695 -0.08052188 -0.00562921]
F_N Arot          [ 0.03427115  1.14115867 -0.07939051 -0.00555011]
--------
F_n vec           [ 0.03427118  1.14116153  0.09109836  0.0063686 ]
F_n rot vec       [ 0.03427118  1.09368137 -0.33824489 -0.02364638]
F_n Arotvec       [ 0.03427118  1.14116153 -0.33923102 -0.02371532]
 % diff           0.014050517173081166

 '''


RQ_1 = '''

     Proton:
F_p               [ 0.72364924 -2.64345634 -0.39180919 -0.02739101]
F_p rot           [ 0.72364924 -2.6004598   0.61561733  0.04303722]
F_p Arot          [ 0.72364924 -2.64345634  0.61650457  0.04309925]
--------
F_p vec           [ 0.72364924 -2.64345634 -0.39180919 -0.02739101]
F_p rot vec       [ 0.72364924 -2.6004598   0.61561733  0.04303722]
F_p Arot vec      [ 0.72364924 -2.64345634  0.61650457  0.04309925]
 % diff           0.0014412238386118318
     Neutron:
F_n               [ 0.03427118  1.14116153 -0.19874397 -0.013894  ]
F_n rot           [ 0.03427118  0.9862801  -0.607454   -0.04246653]
F_N Arot          [ 0.03427118  1.14116153 -0.61065001 -0.04268996]
--------
F_n vec           [ 0.03427118  1.14116153 -0.19874397 -0.013894  ]
F_n rot vec       [ 0.03427118  0.9862801  -0.607454   -0.04246653]
F_n Arotvec       [ 0.03427118  1.14116153 -0.61065001 -0.04268996]
 % diff           0.005261319768251815

'''

RQ_2 = '''

     Proton:
F_p               [ 0.72364924 -2.64345634 -0.67711831 -0.0473367 ]
F_p rot           [ 0.72364924 -2.70618128  0.35061873  0.02451142]
F_p Arot          [ 0.72364924 -2.64345634  0.3493377   0.02442187]
--------
F_p vec           [ 0.72364924 -2.64345634 -0.67711831 -0.0473367 ]
F_p rot vec       [ 0.72364924 -2.70618128  0.35061873  0.02451142]
F_p Arot vec      [ 0.72364924 -2.64345634  0.3493377   0.02442187]
 % diff           0.0036536218255132026
     Neutron:
F_n               [ 0.03427118  1.14116153  0.38094068  0.0266312 ]
F_n rot           [ 0.03427118  1.20108263 -0.06903579 -0.00482623]
F_N Arot          [ 0.03427118  1.14116153 -0.06781203 -0.00474067]
--------
F_n vec           [ 0.03427118  1.14116153  0.38094068  0.0266312 ]
F_n rot vec       [ 0.03427118  1.20108263 -0.06903579 -0.00482623]
F_n Arotvec       [ 0.03427118  1.14116153 -0.06781203 -0.00474067]
 % diff           0.017726538527234582

'''
