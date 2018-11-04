#!/usr/bin/env python


"""

Contains all default Doublet Singlet combining functions.

You can add your own combinations here.

Remember Iso-Vector (doublet - singlet) quantities have no disconnected quark loop contributions :D

"""

upCharge = 2.0/3.0
downCharge = -1.0/3.0

udivdCharge = upCharge/downCharge
ddivuCharge = downCharge/upCharge





def doubDSFun(doublet,singlet):
    return doublet

def singDSFun(doublet,singlet):
    return singlet

def VectorDSFun(doublet,singlet):
    return doublet + singlet

def IsoVectorDSFun(doublet,singlet):
    return doublet - singlet

def ProtonDSFun(doublet,singlet):
    return upCharge*doublet + downCharge*singlet

def NeutronDSFun(doublet,singlet):
    return downCharge*doublet + upCharge*singlet


def PdivNDSFunV2(doublet,singlet):
    DdivS = doublet/singlet
    SdivD = DdivS**-1
    Rat1 = (udivdCharge + SdivD)**-1
    Rat2 = (ddivuCharge + DdivS)**-1
    return Rat1 + Rat2


# def NdivDDSFunV2(doublet,singlet):
#     DdivS = doublet/singlet
#     SdivD = DdivS**-1
#     Rat1 = (ddivuCharge + DdivS)**-1
#     Rat2 = (udivdCharge + SdivD)**-1
#     return Rat1 + Rat2

def PdivNDSFun(doublet,singlet):
    return ProtonDSFun(doublet,singlet)/NeutronDSFun(doublet,singlet)

def NdivPDSFun(doublet,singlet):
    return NeutronDSFun(doublet,singlet)/ProtonDSFun(doublet,singlet)


def NmulPDSFun(doublet,singlet):
    return NeutronDSFun(doublet,singlet)*ProtonDSFun(doublet,singlet)


DSfundict = {'doub':       doubDSFun,
             'sing':       singDSFun,
             'Iso-Vector': IsoVectorDSFun,
             'IsoVector':  IsoVectorDSFun,
             'Vector':     VectorDSFun,
             'Proton':     ProtonDSFun,
             'Neutron':    NeutronDSFun,
             'PdivideN':      PdivNDSFun,
             'NdivideP':      NdivPDSFun,
             'NmultiplyP':      NmulPDSFun,
             'PdivideNV2':    PdivNDSFunV2}
