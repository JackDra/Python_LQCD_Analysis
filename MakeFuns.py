#!/usr/bin/env python

from FileIO import WriteFuns
import numpy as np

# for ival in [32,40,56,64]:
def C2atT(t,p):
    A0,Ep = p[0],p[1]
    Tmint = 33-t[0]
    try:
        return A0*(np.exp(-Ep*t[0])+np.exp(-Ep*Tmint))
    except Exception as err:
        return A0*((-Ep*t[0]).Exp() + (-Ep*Tmint).Exp())
C2atT.file_name = 'C2atT'+str(33)
this_fun = [C2atT,2]
def C2DeratT(t,p):
    A0,Ep = p[0],p[1]
    Tmint = 33-t[0]
    try:
        return [(np.exp(-Ep*t[0])+np.exp(-Ep*Tmint)),-A0*(t[0]*np.exp(-Ep*t[0])+Tmint*np.exp(-Ep*Tmint))]
    except Exception as err:
        return [((-Ep*t[0]).Exp()+(-Ep*Tmint).Exp()),-A0*(t[0]*(-Ep*t[0]).Exp()+Tmint*(-Ep*Tmint).Exp())]
C2DeratT.file_name = 'C2DeratT'+str(33)
WriteFuns(C2atT,C2DeratT)


def C2atT(t,p):
    A0,Ep = p[0],p[1]
    Tmint = 41-t[0]
    try:
        return A0*(np.exp(-Ep*t[0])+np.exp(-Ep*Tmint))
    except Exception as err:
        return A0*((-Ep*t[0]).Exp() + (-Ep*Tmint).Exp())
C2atT.file_name = 'C2atT'+str(41)
this_fun = [C2atT,2]
def C2DeratT(t,p):
    A0,Ep = p[0],p[1]
    Tmint = 41-t[0]
    try:
        return [(np.exp(-Ep*t[0])+np.exp(-Ep*Tmint)),-A0*(t[0]*np.exp(-Ep*t[0])+Tmint*np.exp(-Ep*Tmint))]
    except Exception as err:
        return [((-Ep*t[0]).Exp()+(-Ep*Tmint).Exp()),-A0*(t[0]*(-Ep*t[0]).Exp()+Tmint*(-Ep*Tmint).Exp())]
C2DeratT.file_name = 'C2DeratT'+str(41)
WriteFuns(C2atT,C2DeratT)



def C2atT(t,p):
    A0,Ep = p[0],p[1]
    Tmint = 57-t[0]
    try:
        return A0*(np.exp(-Ep*t[0])+np.exp(-Ep*Tmint))
    except Exception as err:
        return A0*((-Ep*t[0]).Exp() + (-Ep*Tmint).Exp())
C2atT.file_name = 'C2atT'+str(57)
this_fun = [C2atT,2]
def C2DeratT(t,p):
    A0,Ep = p[0],p[1]
    Tmint = 57-t[0]
    try:
        return [(np.exp(-Ep*t[0])+np.exp(-Ep*Tmint)),-A0*(t[0]*np.exp(-Ep*t[0])+Tmint*np.exp(-Ep*Tmint))]
    except Exception as err:
        return [((-Ep*t[0]).Exp()+(-Ep*Tmint).Exp()),-A0*(t[0]*(-Ep*t[0]).Exp()+Tmint*(-Ep*Tmint).Exp())]
C2DeratT.file_name = 'C2DeratT'+str(57)
WriteFuns(C2atT,C2DeratT)



def C2atT(t,p):
    A0,Ep = p[0],p[1]
    Tmint = 65-t[0]
    try:
        return A0*(np.exp(-Ep*t[0])+np.exp(-Ep*Tmint))
    except Exception as err:
        return A0*((-Ep*t[0]).Exp() + (-Ep*Tmint).Exp())
C2atT.file_name = 'C2atT'+str(65)
this_fun = [C2atT,2]
def C2DeratT(t,p):
    A0,Ep = p[0],p[1]
    Tmint = 65-t[0]
    try:
        return [(np.exp(-Ep*t[0])+np.exp(-Ep*Tmint)),-A0*(t[0]*np.exp(-Ep*t[0])+Tmint*np.exp(-Ep*Tmint))]
    except Exception as err:
        return [((-Ep*t[0]).Exp()+(-Ep*Tmint).Exp()),-A0*(t[0]*(-Ep*t[0]).Exp()+Tmint*(-Ep*Tmint).Exp())]
C2DeratT.file_name = 'C2DeratT'+str(65)
WriteFuns(C2atT,C2DeratT)
