#!/usr/bin/env python
## ADD YOUR OWN FUNCTIONS HERE
## Just remember, function(variable, parameter) where:
## variable is list (for multiple dimentional fitting) of dependant variables to fit over
## parameter is list of parameters to extract out of fit.

## you can also define derivatives of the function w.r.t the paremters to feed into the fitting method.
## just remember to add 2 lines to DerOfFun (above)

## you can also add initial conditions for your function in GuessFromFun (above)

import numpy as np
from MiscFuns import logNA

nM_phys = 939.5653
def makexunit(xin):
    try:
        return np.array(len(xin)*[1.0])
    except Exception as err:
        return 1.0


# def ChiDistribution(dof,chi):
#     return (np.exp(-chi/2.)*(chi**((dof/2.) - 1)))/(gamma(dof/2.)*2**(dof/2.))


def Ratio(*a):
    return a[0]/a[1]
    # return a[1]

def RatioDer(*a):
    return [1/a[1],-a[0]/(a[1]**2)]

def eyeFun(*a):
    return a[0]

def eyeFunDer(*a):
    return [1.]


# def IntChiDist(dof,alpha,chiList):
#     # chiList = np.append(chiList,chiList[-1]*10)
#     # if alpha == 0.5:
#     #     for ichi,ifun in zip(chiList,ChiDistribution(dof,np.array(chiList))):
#     #         print ichi/dof, ifun
#     for ichi,chival in enumerate(chiList):
#         intval = integrate.simps(ChiDistribution(dof,np.array(chiList[ichi:])),chiList[ichi:])
#         if intval < alpha:
#             return chiList[ichi]
#     return 1000.0

# def AlphaVsChiDOF(dof):
#     alphalist = np.arange(0.01,1.01,0.01)
#     chiList = np.arange(0.01,100.1,0.01)
#     chinu = []
#     for ialpha in alphalist:
#         chinu.append(IntChiDist(dof,ialpha,chiList))
#     return np.array(chinu)/dof,alphalist

def Chi_Prop_Fit(x,p):
    try:
        return p[0] + p[1]*(-p[2]*x[0]).Exp()
    except Exception as err:
        return p[0] + p[1]*np.exp(-p[2]*x[0])

def Chi_Prop_Fit_2exp(x,p):
    try:
        return p[0] + p[1]*(-p[2]*x[0]).Exp() + p[3]*(-p[4]*x[0]).Exp()
    except Exception as err:
        return p[0] + p[1]*np.exp(-p[2]*x[0]) + p[3]*np.exp(-p[4]*x[0])

def c3FitFun(x,p):
    try:
        tf = np.sqrt(x[0])
    except Exception as err:
        tf = x[0].Sqrt()
    return 1+p[0] + p[1]/tf + p[2]*tf

def c3FFDer(x,p):
    try:
        tf = np.sqrt(x[0])
    except Exception as err:
        tf = x[0].Sqrt()
    return [makexunit(tf),1/tf,tf]

def c3FitFun_log(x,p):
    try:
        tf = np.sqrt(x[0])
    except Exception as err:
        tf = x[0].Sqrt()
    return 1+p[0] +  p[1]/tf + p[2]*tf + p[3]*np.log(tf)

def c3FFDer_log(x,p):
    try:
        tf = np.sqrt(x[0])
    except Exception as err:
        tf = x[0].Sqrt()
    return [makexunit(tf),1/tf,tf,np.log(tf)]

def c3FitFun_nosqrt(x,p):
    return 1+p[0] + p[1]*(x[0]**-2) + p[2]*x[0]**2

def c3FFDer_nosqrt(x,p):
    return [makexunit(x[0]),x[0]**-2,x[0]**2]

def c3FitFun_nosqrt_log(x,p):
    return 1+p[0] + p[1]*(x[0]**-2) + p[2]*x[0]**2 + p[3]*np.log(x[0]**2)

def c3FFDer_nosqrt_log(x,p):
    return [makexunit(x[0]),x[0]**-2,x[0]**2,np.log(x[0]**2)]

def c3FitFun_V2(x,p):
    try:
        tf = np.sqrt(x[0])
    except Exception as err:
        tf = x[0].Sqrt()
    return 1+p[0] + p[1]/tf + p[2]*tf

def c3FFDer_V2(x,p):
    try:
        tf = np.sqrt(x[0])
    except Exception as err:
        tf = x[0].Sqrt()
    return [makexunit(tf),1/tf,tf]


def c3PolyFun(x,p):
    # try:
    #     p0,p1,p2 = np.exp(p[0]),np.exp(p[1]),np.exp(p[2])
    # except Exception as err:
    #     p0,p1,p2 = p[0].Exp(),p[1].Exp(),p[2].Exp()
    p0 = p[0]
    # return 1+p0*x[0]**2 + p1*x[0]**4+ p2*x[0]**6
    return makexunit(x[0])+p0*(x[0])

def c3PFDer(x,p):
    # try:
    #     p0,p1,p2 = np.exp(p[0]),np.exp(p[1]),np.exp(p[2])
    # except Exception as err:
    #     p0,p1,p2 = p[0].Exp(),p[1].Exp(),p[2].Exp()
    # return [p0*x[0]**2,p1*x[0]**4,p2*x[0]**6]
    return [x[0]]

def c3PolyFun2(x,p):
    # try:
    #     p0,p1,p2 = np.exp(p[0]),np.exp(p[1]),np.exp(p[2])
    # except Exception as err:
    #     p0,p1,p2 = p[0].Exp(),p[1].Exp(),p[2].Exp()
    p0,p1 = p[0],p[1]
    # return 1+p0*x[0]**2 + p1*x[0]**4+ p2*x[0]**6
    return makexunit(x[0])+p0*(x[0]) + p1*(x[0]**2)

def c3PFDer2(x,p):
    # try:
    #     p0,p1,p2 = np.exp(p[0]),np.exp(p[1]),np.exp(p[2])
    # except Exception as err:
    #     p0,p1,p2 = p[0].Exp(),p[1].Exp(),p[2].Exp()
    # return [p0*x[0]**2,p1*x[0]**4,p2*x[0]**6]
    return [x[0],x[0]**2]

def c3PolyFun3(x,p):
    # try:
    #     p0,p1,p2 = np.exp(p[0]),np.exp(p[1]),np.exp(p[2])
    # except Exception as err:
    #     p0,p1,p2 = p[0].Exp(),p[1].Exp(),p[2].Exp()
    p0,p1,p2 = p[0],p[1],p[2]
    # return 1+p0*x[0]**2 + p1*x[0]**4+ p2*x[0]**6
    return makexunit(x[0])+p0*(x[0]) + p1*(x[0]**2)+ p2*(x[0]**3)

def c3PFDer3(x,p):
    # try:
    #     p0,p1,p2 = np.exp(p[0]),np.exp(p[1]),np.exp(p[2])
    # except Exception as err:
    #     p0,p1,p2 = p[0].Exp(),p[1].Exp(),p[2].Exp()
    # return [p0*x[0]**2,p1*x[0]**4,p2*x[0]**6]
    return [x[0],x[0]**2,x[0]**3]


def c3PolyFun_skip1(x,p):
    # try:
    #     p0,p1,p2 = np.exp(p[0]),np.exp(p[1]),np.exp(p[2])
    # except Exception as err:
    #     p0,p1,p2 = p[0].Exp(),p[1].Exp(),p[2].Exp()
    p0 = p[0]
    # return 1+p0*x[0]**2 + p1*x[0]**4+ p2*x[0]**6
    return makexunit(x[0])+p0*(x[0]**2)

def c3PFDer_skip1(x,p):
    # try:
    #     p0,p1,p2 = np.exp(p[0]),np.exp(p[1]),np.exp(p[2])
    # except Exception as err:
    #     p0,p1,p2 = p[0].Exp(),p[1].Exp(),p[2].Exp()
    # return [p0*x[0]**2,p1*x[0]**4,p2*x[0]**6]
    return [x[0]**2]

def c3PolyFun2_skip1(x,p):
    # try:
    #     p0,p1,p2 = np.exp(p[0]),np.exp(p[1]),np.exp(p[2])
    # except Exception as err:
    #     p0,p1,p2 = p[0].Exp(),p[1].Exp(),p[2].Exp()
    p0,p1 = p[0],p[1]
    # return 1+p0*x[0]**2 + p1*x[0]**4+ p2*x[0]**6
    return makexunit(x[0])+p0*(x[0]**2) + p1*(x[0]**4)

def c3PFDer2_skip1(x,p):
    # try:
    #     p0,p1,p2 = np.exp(p[0]),np.exp(p[1]),np.exp(p[2])
    # except Exception as err:
    #     p0,p1,p2 = p[0].Exp(),p[1].Exp(),p[2].Exp()
    # return [p0*x[0]**2,p1*x[0]**4,p2*x[0]**6]
    return [x[0]**2,x[0]**4]

def c3PolyFun3_skip1(x,p):
    # try:
    #     p0,p1,p2 = np.exp(p[0]),np.exp(p[1]),np.exp(p[2])
    # except Exception as err:
    #     p0,p1,p2 = p[0].Exp(),p[1].Exp(),p[2].Exp()
    p0,p1,p2 = p[0],p[1],p[2]
    # return 1+p0*x[0]**2 + p1*x[0]**4+ p2*x[0]**6
    return makexunit(x[0])+p0*(x[0]**2) + p1*(x[0]**4)+ p2*(x[0]**6)

def c3PFDer3_skip1(x,p):
    # try:
    #     p0,p1,p2 = np.exp(p[0]),np.exp(p[1]),np.exp(p[2])
    # except Exception as err:
    #     p0,p1,p2 = p[0].Exp(),p[1].Exp(),p[2].Exp()
    # return [p0*x[0]**2,p1*x[0]**4,p2*x[0]**6]
    return [x[0]**2,x[0]**4,x[0]**6]

l8 = np.log(8)
def c3FitFun_NoA(x,p):
    try:
        logx = 2*np.log(x[0])-l8
    except Exception as err:
        logx = 2*x[0].Log()-l8
    try:
        p0,p1 = np.exp(p[0]),np.exp(p[1])
    except Exception as err:
        p0,p1 = p[0].Exp(),p[1].Exp()
    return 1+ p0/(logx - p1)

def c3FFDer_NoA(x,p):
    try:
        logx = 2*np.log(x[0])-l8
    except Exception as err:
        logx = 2*x[0].Log()-l8
    return [1/(p[1]+logx),-p[0]/(p[1]+logx)**2]

def c3ExpFitFun2(x,p):
    try:
        return (np.exp(-p[0]*x[0]) + np.exp(-p[1]*x[0]))/2.
    except Exception as err:
        return ((-p[0]*x[0]).Exp() + (-p[1]*x[0]).Exp())/2.

def c3EFFDer2(x,p):
    try:
        return [-x[0]*np.exp(-p[0]*x[0])/2,-x[0]*np.exp(-p[1]*x[0])/2]
    except Exception as err:
        return [-x[0]*(-p[0]*x[0]).Exp()/2,-x[0]*(-p[1]*x[0]).Exp()/2]

def c3ExpFitFun2_nosqrt(x,p):

    try:
        return (np.exp(-p[0]*(x[0]**2)/8) + np.exp(-p[1]*(x[0]**2)/8))/2.
    except Exception as err:
        return ((-p[0]*(x[0]**2)/8).Exp() + (-p[1]*(x[0]**2)/8).Exp())/2.

def c3EFFDer2_nosqrt(x,p):
    try:
        return [np.exp(-p[0]*(x[0]**2)/8)*(-x[0]**2)/16,np.exp(-p[1]*(x[0]**2)/8)*(-x[0]**2)/16]
    except Exception as err:
        return [(-p[0]*(x[0]**2)/8).Exp()*(-x[0]**2)/16,(-p[1]*(x[0]**2)/8).Exp()*(-x[0]**2)/16]

def c3ExpFitFun3_nosqrt(x,p):

    try:
        return (np.exp(-p[0]*(x[0]**2)/8) + np.exp(-p[1]*(x[0]**2)/8)+ np.exp(-p[2]*(x[0]**2)/8))/3.
    except Exception as err:
        return ((-p[0]*(x[0]**2)/8).Exp() + (-p[1]*(x[0]**2)/8).Exp()+ (-p[2]*(x[0]**2)/8).Exp())/3.

def c3EFFDer3_nosqrt(x,p):
    try:
        return [np.exp(-p[0]*(x[0]**2)/8)*(-x[0]**2)/24,
                np.exp(-p[1]*(x[0]**2)/8)*(-x[0]**2)/24,
                np.exp(-p[2]*(x[0]**2)/8)*(-x[0]**2)/24]
    except Exception as err:
        return [(-p[0]*(x[0]**2)/8).Exp()*(-x[0]**2)/24,
                (-p[1]*(x[0]**2)/8).Exp()*(-x[0]**2)/24,
                (-p[2]*(x[0]**2)/8).Exp()*(-x[0]**2)/24]


def c3ExpFitFun(x,p):
    try:
        return np.exp(-p[0]*x[0])
    except Exception as err:
        return (-p[0]*x[0]).Exp()

def c3EFFDer(x,p):
    try:
        return [-x[0]*np.exp(-p[0]*x[0])]
    except Exception as err:
        return [-x[0]*(-p[0]*x[0]).Exp()]

def c3ExpFitFun_nosqrt(x,p):

    try:
        return np.exp(-p[0]*(x[0]**2)/8)
    except Exception as err:
        return (-p[0]*(x[0]**2)/8).Exp()

def c3EFFDer_nosqrt(x,p):
    try:
        return [np.exp(-p[0]*(x[0]**2)/8)*(-x[0]**2)/8]
    except Exception as err:
        return [(-p[0](x[0]**2)/8).Exp()*(-x[0]**2)/8]


def Chi_Prop_Der(x,p):
    try:
        return [x[0]/x[0],(-p[2]*x[0]).Exp(),-x[0]*p[1]*(-p[2]*x[0]).Exp()]
    except Exception as err:
        return [x[0]/x[0],np.exp(-p[2]*x[0]),-x[0]*p[1]*np.exp(-p[2]*x[0])]

def Chi_Prop_Der_2exp(x,p):
    try:
        return [x[0]/x[0]   ,(-p[2]*x[0]).Exp(),-x[0]*p[1]*(-p[2]*x[0]).Exp()
                            ,(-p[4]*x[0]).Exp(),-x[0]*p[3]*(-p[4]*x[0]).Exp()]
    except Exception as err:
        return [x[0]/x[0]   ,np.exp(-p[2]*x[0]),-x[0]*p[1]*np.exp(-p[2]*x[0])
                            ,np.exp(-p[4]*x[0]),-x[0]*p[3]*np.exp(-p[4]*x[0])]


def Alpha_Prop_Fit(x,p):
    try:
        return p[0] + p[1]*((-p[2]*(x[2]-x[1])).Exp()-(-p[2]*x[1]).Exp())
    except Exception as err:
        return p[0] + p[1]*(np.exp(-p[2]*(x[2]-x[1])) - np.exp(-p[2]*x[1]))


def Alpha_Prop_Fit_plin(x,p):
    Tf = (x[1]-x[0]/2.)
    try:
        return p[0] + p[1]*(-p[2]*Tf).Exp() + p[3]*x[1]
    except Exception as err:
        return p[0] + p[1]*np.exp(-p[2]*Tf) + p[3]*x[1]



def Alpha_Prop_Fit_2exp(x,p):
    Tfmint = x[1]-x[0]
    # TminTf = x[2]-x[1]
    # try:
    #     return p[0] +   p[1]*((-p[2]*TminTf).Exp()-(-p[2]*Tfmint).Exp()) + \
    #                     p[3]*((-p[4]*TminTf).Exp()-(-p[4]*Tfmint).Exp())
    # except Exception as err:
    #     return p[0] +   p[1]*(np.exp(-p[2]*TminTf) - np.exp(-p[2]*Tfmint)) + \
    #                     p[3]*(np.exp(-p[4]*TminTf) - np.exp(-p[4]*Tfmint))
    try:
        return p[0] +   p[1]*(-p[2]*Tfmint).Exp() + \
                        p[3]*(-p[4]*Tfmint).Exp()
    except Exception as err:
        return p[0] +   p[1]*np.exp(-p[2]*Tfmint) + \
                        p[3]*np.exp(-p[4]*Tfmint)

def Alpha_Prop_Der(x,p):
    try:
        return [x[0]/x[0],(-p[2]*(x[2]-x[1])).Exp()-(-p[2]*x[1]).Exp(),
                (x[2]-x[1])*p[1]*(-p[2]*(x[2]-x[1])).Exp()+x[1]*p[1]*(-p[2]*x[1]).Exp()]
    except Exception as err:
        return [x[0]/x[0],(np.exp(-p[2]*(x[2]-x[1])) - np.exp(-p[2]*x[1])),
                (x[2]-x[1])*p[1]*np.exp(-p[2]*(x[2]-x[1])) +x[1]*p[1]*np.exp(-p[2]*x[1])]


def Alpha_Prop_Der_plin(x,p):
    Tf = (x[1]-x[0]/2.)
    try:
        return [x[0]/x[0],
                (-p[2]*Tf).Exp(),
                -Tf*p[1]*(-p[2]*Tf).Exp(),
                x[1]]
    except Exception as err:
        return [x[0]/x[0],
                np.exp(-p[2]*Tf),
                -Tf*p[1]*np.exp(-p[2]*Tf),
                x[1]]


def Alpha_Prop_Der_2exp(x,p):
    Tfmint = x[1]-x[0]
    # TminTf = x[2]-x[1]
    # try:
    #     return [x[0]/x[0],
    #             (-p[2]*TminTf).Exp()-(-p[2]*Tfmint).Exp(),
    #             TminTf*p[1]*(-p[2]*TminTf).Exp()+Tfmint*p[1]*(-p[2]*Tfmint).Exp(),
    #             (-p[4]*TminTf).Exp()-(-p[4]*Tfmint).Exp(),
    #             TminTf*p[3]*(-p[4]*TminTf).Exp()+Tfmint*p[3]*(-p[4]*Tfmint).Exp()]
    # except Exception as err:
    #     return [x[0]/x[0],
    #             (np.exp(-p[2]*TminTf) - np.exp(-p[2]*Tfmint)),
    #             TminTf*p[1]*np.exp(-p[2]*TminTf) +Tfmint*p[1]*np.exp(-p[2]*Tfmint),
    #             (np.exp(-p[3]*TminTf) - np.exp(-p[3]*Tfmint)),
    #             TminTf*p[3]*np.exp(-p[4]*TminTf) +Tfmint*p[3]*np.exp(-p[4]*Tfmint)]
    try:
        return [x[0]/x[0],
                (-p[2]*Tfmint).Exp(),
                -Tfmint*p[1]*(-p[2]*Tfmint).Exp(),
                (-p[4]*Tfmint).Exp(),
                -Tfmint*p[3]*(-p[4]*Tfmint).Exp()]
    except Exception as err:
        return [x[0]/x[0],
                np.exp(-p[2]*Tfmint),
                -Tfmint*p[1]*np.exp(-p[2]*Tfmint),
                np.exp(-p[3]*Tfmint),
                -Tfmint*p[3]*np.exp(-p[4]*Tfmint)]


## only difference here is that the current insersion time is included as x[0]
def RF_Prop_Fit(x,p):
    try:
        return p[0] + p[1]*((-p[2]*(x[2]-x[1])).Exp()-(-p[2]*x[1]).Exp())
    except Exception as err:
        return p[0] + p[1]*(np.exp(-p[2]*(x[2]-x[1])) - np.exp(-p[2]*x[1]))

def RF_Prop_Der(x,p):
    try:
        return [x[1]/x[1],(-p[2]*(x[2]-x[1])).Exp()-(-p[2]*x[1]).Exp(),
                (x[2]-x[1])*p[1]*(-p[2]*(x[2]-x[1])).Exp()+x[1]*p[1]*(-p[2]*x[1]).Exp()]
    except Exception as err:
        return [x[1]/x[1],(np.exp(-p[2]*(x[2]-x[1])) - np.exp(-p[2]*x[1])),
                (x[2]-x[1])*p[1]*np.exp(-p[2]*(x[2]-x[1])) +x[1]*p[1]*np.exp(-p[2]*x[1])]


def CreateFFO(thislen):
    def FormFactorO(x,p):
        return sum([ip*ix for ip,ix in zip(p[:thislen],x[:thislen])])
    def FormFactorODer(x,p):
        return np.array(x[:thislen]).tolist()
    return FormFactorO, FormFactorODer

def ParmDivX(x,p):
    return p[0] * (x[0]**(-1))

def ParmDivXDer(x,p):
    return [x[0]**(-1)]

## Sigma/(1/m_u + 1/m_d + 1/m_s) = Sigma/(2/m_ud + 1/m_s)
def Chipt_Chit(x,p):
    sum_part = 2*x[0]**-1 + x[1]**-1
    return p[0]*sum_part**-1

## A e^(-mt) + B(1-e^(-mt))
def RandTFitFun(x,p):
    try:
        thisexp = np.exp(-p[2]*x[0])
    except Exception as err:
        thisexp = (-p[2]*x[0]).Exp()
    return p[0]*thisexp + p[1]*(1-thisexp)


def ChitFitFun(x,p):
    Denom = (x[0]*(logNA(x[0])+p[2]))**(-1)
    return (p[0] *x[0] + p[1])*Denom

def ChitFitFunDer(x,p):
    Denom = (x[0]*(logNA(x[0])+p[2]))**(-1)
    return np.array([Denom*x[0],Denom,-x[0]*(Denom**2)])

## A*x*exp(-m*x)
def A_x_exp_minus_m_x(x,p):
    try:
        return p[0]*x[0]*np.exp(-p[1]*x[0])
    except Exception as err:
        return p[0]*x[0]*(-p[1]*x[0]).Exp()

## broken for some reason
def ParmDivXP(x,p):
    return p[0] * (x[0]**(-1)) + p[1]

def ParmDivXPDer(x,p):
    return np.array([x[0]**(-1),makexunit(x[0])])

def DPfitfun(x,p):
    return p[0]/((1+(x[0]/p[1]))**2)

def DPfitfunDer(x,p):
    # return [1/((1+(x[0]/p[1]))**2), (-2*p[0]*p[1]**2)/((x[0]+p[1])**3)]
    return [1/((1+(x[0]/p[1]))**2),(2*p[0]*x[0]/(p[1]**2))/((1+(x[0]/p[1]))**3)]

def DPfitfunOnePar(x,p):
    return 1/((1+(x[0]/p[0]))**2)

def DPfitfunOneParDer(x,p):
    # return [1/((1+(x[0]/p[1]))**2), (-2*p[0]*p[1]**2)/((x[0]+p[1])**3)]
    return [(-2/p[0])*1/((1+(x[0]/p[0]))**3)]

def DPfitfun2(x,p):
    return p[0]/(1+(x[0]/p[1])**2)**2

def DPfitfun2Der(x,p):
    return [1/(1+(x[0]/p[1])**2)**2,-4*p[0]*x[0]/(p[1]**2*(1+(x[0]/p[1])**2)**3)]

def FormFactorO1(x,p):return p[0]*x[0]
def FormFactorO1Der(x,p):return [x[0]]

def FormFactorO2(x,p):return p[0]*x[0]+p[1]*x[1]
def FormFactorO2Der(x,p): return [x[0],x[1]]

def FormFactorO3(x,p):return p[0]*x[0]+p[1]*x[1]+p[2]*x[2]
def FormFactorO3Der(x,p):return [x[0],x[1],x[2]]

def FormFactorO4(x,p):return p[0]*x[0]+p[1]*x[1]+p[2]*x[2]+p[3]*x[3]
def FormFactorO4Der(x,p):return [x[0],x[1],x[2],x[3]]

def FormFactorO5(x,p):return p[0]*x[0]+p[1]*x[1]+p[2]*x[2]+p[3]*x[3]+p[4]*x[4]
def FormFactorO5Der(x,p):return [x[0],x[1],x[2],x[3],x[4]]

def FormFactorO6(x,p):return p[0]*x[0]+p[1]*x[1]+p[2]*x[2]+p[3]*x[3]+p[4]*x[4]+p[5]*x[5]
def FormFactorO6Der(x,p):return [x[0],x[1],x[2],x[3],x[4],x[5]]

def FormFactorO7(x,p):return p[0]*x[0]+p[1]*x[1]+p[2]*x[2]+p[3]*x[3]+p[4]*x[4]+p[5]*x[5]+p[6]*x[6]
def FormFactorO7Der(x,p):return [x[0],x[1],x[2],x[3],x[4],x[5],x[6]]

def FormFactorO8(x,p):return p[0]*x[0]+p[1]*x[1]+p[2]*x[2]+p[3]*x[3]+p[4]*x[4]+p[5]*x[5]+p[6]*x[6]+p[7]*x[7]
def FormFactorO8Der(x,p):return [x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7]]

FormList = [FormFactorO8,
            FormFactorO7,
            FormFactorO6,
            FormFactorO5,
            FormFactorO4,
            FormFactorO3,
            FormFactorO2,
            FormFactorO1]


def ConstantFitFun(x,p):
    return p[0]*makexunit(x[0])

def ConstFFDer(x,p):
    return makexunit(x[0])



def CreateOORFF(Const):
    def OneOnRootNFitFun(x,p):
        return Const + p[0]*(np.array(x[0])**(-1/2.)-1)
    return OneOnRootNFitFun

def OORNFFDer(x,p):
    return [(-p[0]/2.)*(np.array(x[0])**(-3/2.))]

# def SchiffFitFun(x,p):
#     return p[0]*np.array(x[0])**2+p[1] + p[2]*np.array(x[1])
#
# def SchiffFFDer(x,p):
#     return [np.array(x[0])**2,makexunit(x[0]),np.array(x[1])]
#
#
# def SchiffFitFun(x,p):
#     return p[0]*np.array(x[0])**2 + p[1]*np.array(x[1]) + p[2]*np.array(x[1])**2
#
# def SchiffFFDer(x,p):
#     return [np.array(x[0])**2,np.array(x[1]),np.array(x[1])**2]
# def SchiffFitFun(x,p):
#     x = np.array(x)
#     return p[0]*x[0]**2 + p[1]*x[1]**2 + p[2]*np.log(x[1]**2)*x[1]**2
#
# def SchiffFFDer(x,p):
#     x = np.array(x)
#     return [x[0]**2,x[1]**2,np.log(x[1]**2)*x[1]**2]

def SchiffFitFun(x,p):
    return p[0]*np.array(x[0])**2 + p[1]

def SchiffFFDer(x,p):
    return [np.array(x[0])**2,makexunit(x[1])]

def SchiffFitFun_noa(x,p):
    return p[1]

def SchiffFFDer_noa(x,p):
    return [makexunit(x[1])]

def c3ContFitFun(x,p):
    return p[0] + p[1]*np.array(x[0])**2 + p[2]*np.array(x[1])

def c3ContFFDer(x,p):
    return [makexunit(x[0]),np.array(x[0])**2,np.array(x[1])]

def hx_fun(x):
    out = np.sqrt(1+1/x)
    out = out * np.log((out + 1)/(out-1))
    out = out - (2*(1+x/3))
    return (-15/4) * out

def fullq2FitFun(x,p):
    ## x[0] is q^2, x[1] is (m_\pi)^2
    ## both must be in same units, I suggest GeV!
    ## equation 7 in notes
    term_2 = hx_fun(x[0]/(4*x[1]))*8*x[1]/5
    return p[1] + p[0] * (x[0] -term_2)

def fullq2FFDer(x,p):
    term_2 = hx_fun(x[0]/(4*x[1]))*8*x[1]/5
    return [x[0] -term_2,makexunit(x[0])]

def LinearFitFun(x,p):
    return p[0]*np.array(x[0])+p[1]

def LinFFDer(x,p):
    return [x[0],makexunit(x[0])]


def LinearFitFun_x0(x,p):
    return p[0]*np.array(x[0])

def LinFFDer_x0(x,p):
    return [x[0]]

def SqPlLogFitFun(x,p):
    this_x = np.array(x[0])**2
    return p[0]*this_x+p[1] + p[2]*this_x * np.log(this_x)

def SqPlLogFFDer(x,p):
    this_x = np.array(x[0])**2
    return [this_x,makexunit(x[0]),this_x*np.log(this_x)]


def SqPlLogFitFun_x0(x,p):
    this_x = np.array(x[0])**2
    return p[0]*this_x + p[1]*this_x * np.log(this_x)

def SqPlLogFFDer_x0(x,p):
    this_x = np.array(x[0])**2
    return [this_x,this_x*np.log(this_x)]

def SquaredFitFun(x,p):
    return p[0]*np.array(x[0])**2+p[1]

def SquaredFFDer(x,p):
    return [np.array(x[0])**2,makexunit(x[0])]


def SquaredFitFun_x0(x,p):
    return p[0]*np.array(x[0])**2

def SquaredFFDer_x0(x,p):
    return [np.array(x[0])**2]


def ContPlLogFitFun(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    return p[0]*this_x+p[1] + p[2]*this_x * np.log(this_x) + p[3]*this_a

def ContPlLogFFDer(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    return [this_x,makexunit(x[0]),this_x*np.log(this_x),this_a]


def ContPlLogFitFun_x0(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    return p[0]*this_x + p[1]*this_x * np.log(this_x) + p[2]*this_a

def ContPlLogFFDer_x0(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    return [this_x,this_x*np.log(this_x),this_a]

def ContPlLogFitFun_x0_dN(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    ## physical neutron mass = 939.5653 MeV
    return p[0]*this_x + p[1]*this_x * np.log(this_x/(nM_phys**2)) + p[2]*this_a

def ContPlLogFFDer_x0_dN(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    return [this_x,this_x*np.log(this_x/(nM_phys**2)),this_a]

## third dim > 0 means neutron, < 0 means proton
def ContPlLogFitFun_NP_nompi3(x,p):
    if isinstance(x[2],(list,tuple,np.ndarray)):
        output = []
        for ix0,ix1,ix2 in zip(x[0],x[1],x[2]):
            if ix2 > 0:
                ##neutron case
                output.append(ContPlLogFitFun_x0_dN([ix0,ix1],[p[0],-p[2],p[3]]))
            elif ix2 < 1:
                output.append(ContPlLogFitFun_x0_dN([ix0,ix1],[p[1],p[2],p[4]]))
        return np.array(output)
    else:
        if x[2] > 0:
            return ContPlLogFitFun_x0_dN(x,[p[0],-p[2],p[3]])
        elif x[2] < 1:
            return ContPlLogFitFun_x0_dN(x,[p[1],p[2],p[4]])

def ContPlLogFFDer_NP_nompi3(x,p):
    if isinstance(x[2],(list,tuple,np.ndarray)):
        output = []
        for ix0,ix1,ix2 in zip(x[0],x[1],x[2]):
            out_list = np.array([makexunit(ix0)-makexunit(ix0) for i in p])
            if ix2 > 0:
                ##neutron case
                out_list[0],out_list[2],out_list[3] = ContPlLogFFDer_x0_dN([ix0,ix1],[p[0],p[2],p[3]])
                out_list[2] = -out_list[2]
            elif ix2 < 1:
                out_list[1],out_list[2],out_list[4] = ContPlLogFFDer_x0_dN([ix0,ix1],[p[1],p[2],p[4]])
            output.append(out_list)
        output = list(np.array(output).T)
    else:
        output = [0 for i in p]
        if x[2] > 0:
            output[0],output[2],output[3] =  ContPlLogFFDer_x0_dN(x,[p[0],p[2],p[3]])
            output[2] = -output[2]
        elif x[2] < 1:
            output[1],output[2],output[4] = ContPlLogFFDer_x0_dN(x,[p[1],p[2],p[4]])
    return output

def ContPlLogFitFun_x0_N(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    ## physical neutron mass = nM_phys MeV
    log_term = (np.log(this_x/(nM_phys**2))+(np.pi*np.array(x[1])/(2*nM_phys)))
    return p[0]*this_x - p[1]*this_x * log_term + p[2]*this_a

def ContPlLogFFDer_x0_N(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    log_term = (np.log(this_x/(nM_phys**2))+(np.pi*np.array(x[1])/(2*nM_phys)))
    return [this_x,-this_x*log_term,this_a]

def ContPlLogFitFun_x0_P(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    ## physical neutron mass = nM_phys MeV
    log_term = (np.log(this_x/(nM_phys**2))+(2*np.pi*np.array(x[1])/(nM_phys)))
    return p[0]*this_x + p[1]*this_x * log_term + p[2]*this_a

def ContPlLogFFDer_x0_P(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    log_term = (np.log(this_x/(nM_phys**2))+(2*np.pi*np.array(x[1])/(nM_phys)))
    return [this_x,this_x*log_term,this_a]



## third dim > 0 means neutron, < 0 means proton
def ContPlLogFitFun_NP(x,p):
    if isinstance(x[2],(list,tuple,np.ndarray)):
        output = []
        for ix0,ix1,ix2 in zip(x[0],x[1],x[2]):
            if ix2 > 0:
                output.append(ContPlLogFitFun_x0_N([ix0,ix1],[p[0],p[2],p[3]]))
            elif ix2 < 1:
                output.append(ContPlLogFitFun_x0_P([ix0,ix1],[p[1],p[2],p[4]]))
        return np.array(output)
    else:
        if x[2] > 0:
            return ContPlLogFitFun_x0_N(x,[p[0],p[2],p[3]])
        elif x[2] < 1:
            return ContPlLogFitFun_x0_P(x,[p[1],p[2],p[4]])

def ContPlLogFFDer_NP(x,p):
    if isinstance(x[2],(list,tuple,np.ndarray)):
        output = []
        for ix0,ix1,ix2 in zip(x[0],x[1],x[2]):
            out_list = np.array([makexunit(ix0)-makexunit(ix0) for i in p])
            if ix2 > 0:
                out_list[0],out_list[2],out_list[3] =  ContPlLogFFDer_x0_N([ix0,ix1],[p[0],p[2],p[3]])
            elif ix2 < 1:
                out_list[1],out_list[2],out_list[4] = ContPlLogFFDer_x0_P([ix0,ix1],[p[1],p[2],p[4]])
            output.append(out_list)
        output = list(np.array(output).T)
    else:
        output = [0 for i in p]
        if x[2] > 0:
            output[0],output[2],output[3] =  ContPlLogFFDer_x0_N(x,[p[0],p[2],p[3]])
        elif x[2] < 1:
            output[1],output[2],output[4] = ContPlLogFFDer_x0_P(x,[p[1],p[2],p[4]])
    return output

## third dim > 0 means neutron, < 0 means proton
def ContPlLogFitFun_NP_mC2(x,p):
    if isinstance(x[2],(list,tuple,np.ndarray)):
        output = []
        for ix0,ix1,ix2 in zip(x[0],x[1],x[2]):
            if ix2 > 0:
                output.append(ContPlLogFitFun_x0_N([ix0,ix1],[p[0],-p[2],p[3]]))
            elif ix2 < 1:
                output.append(ContPlLogFitFun_x0_P([ix0,ix1],[p[1],-p[2],p[4]]))
        return np.array(output)
    else:
        if x[2] > 0:
            return ContPlLogFitFun_x0_N(x,[p[0],-p[2],p[3]])
        elif x[2] < 1:
            return ContPlLogFitFun_x0_P(x,[p[1],-p[2],p[4]])

def ContPlLogFFDer_NP_mC2(x,p):
    if isinstance(x[2],(list,tuple,np.ndarray)):
        output = []
        for ix0,ix1,ix2 in zip(x[0],x[1],x[2]):
            out_list = np.array([makexunit(ix0)-makexunit(ix0) for i in p])
            if ix2 > 0:
                out_list[0],out_list[2],out_list[3] =  ContPlLogFFDer_x0_N([ix0,ix1],[p[0],-p[2],p[3]])
            elif ix2 < 1:
                out_list[1],out_list[2],out_list[4] = ContPlLogFFDer_x0_P([ix0,ix1],[p[1],p[2],p[4]])
            output.append(out_list)
        output = list(np.array(output).T)
    else:
        output = [0 for i in p]
        if x[2] > 0:
            output[0],output[2],output[3] =  ContPlLogFFDer_x0_N(x,[p[0],-p[2],p[3]])
        elif x[2] < 1:
            output[1],output[2],output[4] = ContPlLogFFDer_x0_P(x,[p[1],p[2],p[4]])
    return output

def ContPlLogFitFun_x0_N_latmN(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    ## mpi / mN
    this_mrat = np.array(x[2])**2
    log_term = np.log(this_mrat)-np.pi*np.array(x[2])/2
    return p[0]*this_x - p[1]*this_x * log_term + p[2]*this_a

def ContPlLogFFDer_x0_N_latmN(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    this_mrat = np.array(x[2])**2
    log_term = np.log(this_mrat)-np.pi*np.array(x[2])/2
    return [this_x,-this_x*log_term,this_a]

def ContPlLogFitFun_x0_P_latmN(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    ## physical neutron mass = this_mN MeV
    this_mrat = np.array(x[2])**2
    log_term = np.log(this_mrat)-2*np.pi*np.array(x[2])
    return p[0]*this_x + p[1]*this_x * log_term + p[2]*this_a

def ContPlLogFFDer_x0_P_latmN(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    this_mrat = np.array(x[2])**2
    log_term = np.log(this_mrat)-2*np.pi*np.array(x[2])
    return [this_x,this_x*log_term,this_a]

## third dim > 0 means neutron, < 0 means proton
## forth dim is nucleon mass in same units as pion mass
def ContPlLogFitFun_NP_latmN(x,p):
    if isinstance(x[2],(list,tuple,np.ndarray)):
        output = []
        for ix0,ix1,ix2,ix3 in zip(x[0],x[1],x[2],x[3]):
            if ix2 > 0:
                output.append(ContPlLogFitFun_x0_N_latmN([ix0,ix1,ix3],[p[0],p[2],p[3]]))
            elif ix2 < 1:
                output.append(ContPlLogFitFun_x0_P_latmN([ix0,ix1,ix3],[p[1],p[2],p[4]]))
        return np.array(output)
    else:
        if x[2] > 0:
            return ContPlLogFitFun_x0_N_latmN([x[0],x[1],x[3]],[p[0],p[2],p[3]])
        elif x[2] < 1:
            return ContPlLogFitFun_x0_P_latmN([x[0],x[1],x[3]],[p[1],p[2],p[4]])

def ContPlLogFFDer_NP_latmN(x,p):
    if isinstance(x[2],(list,tuple,np.ndarray)):
        output = []
        for ix0,ix1,ix2,ix3 in zip(x[0],x[1],x[2],x[3]):
            out_list = np.array([makexunit(ix0)-makexunit(ix0) for i in p])
            if ix2 > 0:
                out_list[0],out_list[2],out_list[3] =  ContPlLogFFDer_x0_N_latmN([ix0,ix1,ix3],[p[0],p[2],p[3]])
            elif ix2 < 1:
                out_list[1],out_list[2],out_list[4] = ContPlLogFFDer_x0_P_latmN([ix0,ix1,ix3],[p[1],p[2],p[4]])
            output.append(out_list)
        output = list(np.array(output).T)
    else:
        output = [0 for i in p]
        if x[2] > 0:
            output[0],output[2],output[3] =  ContPlLogFFDer_x0_N_latmN([x[0],x[1],x[3]],[p[0],p[2],p[3]])
        elif x[2] < 1:
            output[1],output[2],output[4] = ContPlLogFFDer_x0_P_latmN([x[0],x[1],x[3]],[p[1],p[2],p[4]])
    return output

def ContFitFun(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    return p[0]*this_x+p[1] + p[2]*this_a

def ContFFDer(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    return [this_x,makexunit(x[0]),this_a]


def ContFitFun_x0(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    return p[0]*this_x + p[1]*this_a

def ContFFDer_x0(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    return [this_x,this_a]


def ContPlLogFitFun_mix(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    this_mix = this_a * this_x
    return p[0]*this_x+p[1] + p[2]*this_x * np.log(this_x) + p[3]*this_a + p[4]*this_mix

def ContPlLogFFDer_mix(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    this_mix = this_a * this_x
    return [this_x,makexunit(x[0]),this_x*np.log(this_x),this_a,this_mix]


def ContPlLogFitFun_x0_mix(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    this_mix = this_a * this_x
    return p[0]*this_x + p[1]*this_x * np.log(this_x) + p[2]*this_a + p[3]*this_mix

def ContPlLogFFDer_x0_mix(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    this_mix = this_a * this_x
    return [this_x,this_x*np.log(this_x),this_a,this_mix]

def ContFitFun_mix(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    this_mix = this_a * this_x
    return p[0]*this_x+p[1] + p[2]*this_a + p[3]*this_mix

def ContFFDer_mix(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    this_mix = this_a * this_x
    return [this_x,makexunit(x[0]),this_a,this_mix]


def ContFitFun_x0_mix(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    this_mix = this_a * this_x
    return p[0]*this_x + p[1]*this_a + p[2]*this_mix

def ContFFDer_x0_mix(x,p):
    this_x = np.array(x[1])**2
    this_a = np.array(x[0])**2
    this_mix = this_a * this_x
    return [this_x,this_a,this_mix]




def LogFitFun(x,p):
    return p[0]*logNA(x[0])+p[1]

def LogFFDer(x,p):
    return [logNA(x[0]),makexunit(x[0])]

def LogFitFun_x0(x,p):
    return p[0]*logNA(x[0])

def LogFFDer_x0(x,p):
    return [logNA(x[0])]


def LinearFF_2D(x,p):
    return p[0]*np.array(x[0])+p[1]*np.array(x[1])+p[2]

def LinFFDer_2D(x,p):
    return [x[0],x[1],makexunit(x[0])]

def TestTwoVarFitFun(x,p):
    return p[0]*x[0]+ p[1]*x[1]**2

def TestTwoVarFFDer(x,p):
    return [x[0],x[1]**2]

##ONE STATE FIT FUNCTIONS###

def C2OSFAntiper(t,p,tsink):
    A0,Ep = p[0],p[1]
    Tmint = tsink-t[0]
    try:
        return A0*(np.exp(-Ep*t[0])+np.exp(-Ep*Tmint))
    except Exception as err:
        return A0*((-Ep*t[0]).Exp() + (-Ep*Tmint).Exp())

def C2OSFAntiperDer(t,p,tsink):
    A0,Ep = p[0],p[1]
    Tmint = tsink-t[0]
    try:
        return [(np.exp(-Ep*t[0])+np.exp(-Ep*Tmint)),-A0*(t[0]*np.exp(-Ep*t[0])+Tmint*np.exp(-Ep*Tmint))]
    except Exception as err:
        return [((-Ep*t[0]).Exp()+(-Ep*Tmint).Exp()),-A0*(t[0]*(-Ep*t[0]).Exp()+Tmint*(-Ep*Tmint).Exp())]


def C2OneStateFitFun(t,p):
    try:
        A0,Ep = p[0],np.exp(p[1])
        return A0*np.exp(-Ep*t[0])
    except Exception as err:
        A0,Ep = p[0],p[1].Exp()
        output = A0*(-Ep*t[0]).Exp()
        output.Stats()
        return output

def C2OneStateFitFunNoExp(t,p):
    A0,Ep = p[0],p[1]
    try:
        return A0*np.exp(-Ep*t[0])
    except Exception as err:
        output = A0*(-Ep*t[0]).Exp()
        output.Stats()
        return output


def C2OSFFDer(t,p):
    try:
        A0,Ep = p[0],np.exp(p[1])
        return [np.exp(-Ep*t[0]),-t[0]*Ep*A0*np.exp(-Ep*t[0])]
    except Exception as err:
        A0,Ep = p[0],p[1].Exp()
        output = [(-Ep*t[0]).Exp(),-t[0]*Ep*A0*(-Ep*t[0]).Exp()]
        output[0].Stats()
        output[1].Stats()
        return output

def C2OSFFNoExpDer(t,p):
    A0,Ep = p[0],p[1]
    try:
        return [np.exp(-Ep*t[0]),-t[0]*Ep*A0*np.exp(-Ep*t[0])]
    except Exception as err:
        output = [(-Ep*t[0]).Exp(),-t[0]*Ep*A0*(-Ep*t[0]).Exp()]
        output[0].Stats()
        output[1].Stats()
        return output


def C3OneStateFitFun(tf,p):
    FitA0 = tf[2]
    FitAp0 = tf[4]
    Fitm0 = np.exp(tf[3])
    FitEp0 = np.exp(tf[5])
    return np.sqrt(np.abs(FitA0*FitAp0))*np.exp(-Fitm0*tf[1])*np.exp(-(FitEp0-Fitm0)*tf[0])*p[0]

def C3OneStateFitFunNoExp(tf,p):
    FitA0 = tf[2]
    FitAp0 = tf[4]
    Fitm0 = tf[3]
    FitEp0 = tf[5]
    return np.sqrt(np.abs(FitA0*FitAp0))*np.exp(-Fitm0*tf[1])*np.exp(-(FitEp0-Fitm0)*tf[0])*p[0]

def C3OSFFDer(tf,p):
    FitA0 = tf[2]
    FitAp0 = tf[4]
    Fitm0 = np.exp(tf[3])
    FitEp0 = np.exp(tf[5])
    return [np.sqrt(np.abs(FitA0*FitAp0))*np.exp(-Fitm0*tf[1])*np.exp(-(FitEp0-Fitm0)*tf[0])]


##Momenta TwoStateFit##

def C2TwoStateFitFun(t,p):
    try:
        Am,Amp,m,dm = p[0],p[1],np.exp(p[2]),np.exp(p[3])
        return Am*(np.exp(m*-t[0]) + Amp*np.exp((m+dm)*-t[0]))
    except Exception as err:
        Am,Amp,m,dm = p[0],p[1],p[2].Exp(),p[3].Exp()
        output = Am*((m*-t[0]).Exp() + Amp*((m+dm)*-t[0]).Exp())
        output.Stats()
        return output

def C2TwoStateFitFunNoExp(t,p):
    Am,Amp,m,dm = p[0],p[1],p[2],p[3]
    try:
        return Am*(np.exp(m*-t[0]) + Amp*np.exp((m+dm)*-t[0]))
    except Exception as err:
        return Am*((m*-t[0]).Exp() + Amp*((m+dm)*-t[0]).Exp())

def C2TwoStateFitFunCM(t,p):
    Alen = (len(p)-2)/2
    Am,Amp,m,dm = np.array(p[:Alen]),np.array(p[Alen:-2]),np.exp(p[-2]),np.exp(p[-1])
    smpick = list(map(int,t[1]))
    return Am[smpick]*(np.exp(-m*t[0]) + Amp[smpick]*np.exp(-(m+dm)*t[0]))

def C2TwoStateFitFunCMNoExp(t,p):
    Alen = (len(p)-2)/2
    Am,Amp,m,dm = np.array(p[:Alen]),np.array(p[Alen:-2]),p[-2],p[-1]
    smpick = list(map(int,t[1]))
    return Am[smpick]*(np.exp(-m*t[0]) + Amp[smpick]*np.exp(-(m+dm)*t[0]))


def C2TSFLineFun(t,p):
    Am,Amp,m,dm = p[0],p[1],p[2],p[3]
    try:
        return Am*(np.exp(m*-t[0]) + Amp*np.exp((m+dm)*-t[0]))
    except Exception as err:
        return Am*((m*(-t[0])).exp(1) + Amp*((m+dm)*(-t[0])).exp(1))


# def C2TSFFCMDer(t,p):
#     Alen = (len(p)-2)/2
#     Am,Amp,m,dm = np.array(p[:Alen]),np.array(p[Alen:-2]),np.exp(p[-2]),np.exp(p[-1])
#     smpick = map(int,t[1])
#     pder = []
#     for ism,Ami in enumerate(Am): # Am derivatives
#         delta = VecDelta(smpick,ism)
#         pder.append(delta*(np.exp(-m*t[0]) + Amp[smpick]*np.exp(-(m+dm)*t[0])))
#     for ism,Ampi in enumerate(Amp): # Amp derivatives
#         delta = VecDelta(smpick,ism)
#         pder.append(delta*Am[smpick]*np.exp(-(m+dm)*t[0]))
#     pder.append(-t[0]*m*Am[smpick]*(np.exp(-m*t[0]) + Amp[smpick]*np.exp(-(m+dm)*t[0]))) #mass
#     pder.append(-t[0]*dm*Am[smpick]*Amp[smpick]*np.exp(-(m+dm)*t[0])) #Dm
#     return pder

def C2TSFFDer(t,p):
    try:
        Am,Amp,m,dm = p[0],p[1],np.exp(p[2]),np.exp(p[3])
        pder = []
        pder.append(np.exp(m*-t[0] )+ Amp*np.exp((m+dm)*-t[0])) ## dC/dAm
        pder.append(Am*np.exp((m+dm)*-t[0])) ## dC/dAmp
        pder.append(-t[0]*m*Am*(np.exp(-m*t[0]) + Amp*np.exp((m+dm)*-t[0]))) # mass
        pder.append(-t[0]*dm*Am*Amp*np.exp((m+dm)*-t[0])) # Dm
        return pder
    except Exception as err:
        Am,Amp,m,dm = p[0],p[1],p[2].Exp(),p[3].Exp()
        pder = []
        pder.append((m*-t[0] ).Exp()+ Amp*((m+dm)*-t[0]).Exp()) ## dC/dAm
        pder.append(Am*((m+dm)*-t[0]).Exp()) ## dC/dAmp
        pder.append(-t[0]*m*Am*((m*-t[0]).Exp() + Amp*((m+dm)*-t[0]).Exp())) # mass
        pder.append(-t[0]*dm*Am*Amp*((m+dm)*-t[0]).Exp()) # Dm
        for ip in pder: ip.Stats()
        return pder


def C2TSFFNoExpDer(t,p):
    Am,Amp,m,dm = p[0],p[1],p[2],p[3]
    pder = []
    try:
        pder.append(np.exp(m*-t[0] )) ## dC/dAm
        pder.append(np.exp((m+dm)*-t[0])) ## dC/dAmp
        pder.append(-t[0]*Am*(np.exp(m*-t[0]) + Amp*np.exp((m+dm)*-t[0]))) # mass
        pder.append(-t[0]*Am*Amp*np.exp((m+dm)*-t[0])) # Dm
        return pder
    except Exception as err:
        pder.append((m*-t[0] ).Exp()) ## dC/dAm
        pder.append(((m+dm)*-t[0]).Exp()) ## dC/dAmp
        pder.append(-t[0]*Am*((m*-t[0]).Exp() + Amp*((m+dm)*-t[0]).Exp())) # mass
        pder.append(-t[0]*Am*Amp*((m+dm)*-t[0]).Exp()) # Dm
        for ip in pder: ip.Stats()
        return pder


def C3TSFLineFun(Vals,thistcurr,thistsink):
    FitA0 = Vals[0]
    # FitA1 = Vals[1]
    Fitm0 = Vals[2]
    FitDm = Vals[3]
    FitB00 = Vals[4]
    FitB10 = Vals[5]
    # FitB01 = Vals[6]
    FitB11 = Vals[7]
    output =  (FitA0*np.exp(-Fitm0*thistsink) *
               (FitB00 + FitB10*(np.exp(-FitDm*thistcurr)+ np.exp(-FitDm*(thistsink-thistcurr))) +
                FitB11 *np.exp(-FitDm*thistsink)))
    return output



def C3MomTwoStateFitFun(tf,p):
    FitA0 = tf[2]
    # FitA1 = tf[3]
    FitAp0 = tf[6]
    # FitAp1 = tf[7]
    Fitm0 = np.exp(tf[4])
    FitDm = np.exp(tf[5])
    FitEp = np.exp(tf[8])
    FitDEp = np.exp(tf[9])
    if all(FitEp == Fitm0) and all(FitA0 == FitAp0):
        if all([tf[1][0]==itf for itf in tf[1]]):

            output = (FitA0*np.exp(-Fitm0*tf[1])*
                      (p[0] + p[1]*(np.exp(-FitDm*tf[0])+np.exp(-FitDm*(tf[1]-tf[0])))))
        else:
            output =  (FitA0*np.exp(-Fitm0*tf[1]) *
                       (p[0] + p[1]*(np.exp(-FitDm*tf[0])+ np.exp(-FitDm*(tf[1]-tf[0]))) +
                        p[3] *np.exp(-FitDm*tf[1])))
    else:
        if all([tf[1][0]==itf for itf in tf[1]]):
            output =  (np.sqrt(np.abs(FitA0*FitAp0))*np.exp(-Fitm0*tf[1]) * np.exp(-(FitEp-Fitm0)*tf[0])*
                       (p[0] + p[1]*np.exp(-FitDEp*tf[0])+ p[2]*np.exp(-FitDm*(tf[1]-tf[0]))))
        else:
            output =  (np.sqrt(np.abs(FitA0*FitAp0))*np.exp(-Fitm0*tf[1]) * np.exp(-(FitEp-Fitm0)*tf[0])*
                       (p[0] + p[1]*np.exp(-FitDEp*tf[0])+ p[2]*np.exp(-FitDm*(tf[1]-tf[0])) +
                        p[3] *np.exp(-FitDEp*tf[0])*np.exp(-FitDm*(tf[1]-tf[0]))))
    return output


def C3MomTwoStateFitFunNoExp(tf,p):
    FitA0 = tf[2]
    # FitA1 = tf[3]
    FitAp0 = tf[6]
    # FitAp1 = tf[7]
    Fitm0 = tf[4]
    FitDm = tf[5]
    FitEp = tf[8]
    FitDEp = tf[9]
    if all(FitEp == Fitm0) and all(FitA0 == FitAp0):
        if all([tf[1][0]==itf for itf in tf[1]]):

            output = (FitA0*np.exp(-Fitm0*tf[1])*
                      (p[0] + p[1]*(np.exp(-FitDm*tf[0])+np.exp(-FitDm*(tf[1]-tf[0])))))
        else:
            output =  (FitA0*np.exp(-Fitm0*tf[1]) *
                       (p[0] + p[1]*(np.exp(-FitDm*tf[0])+ np.exp(-FitDm*(tf[1]-tf[0]))) +
                        p[3] *np.exp(-FitDm*tf[1])))
    else:
        if all([tf[1][0]==itf for itf in tf[1]]):
            output =  (np.sqrt(np.abs(FitA0*FitAp0))*np.exp(-Fitm0*tf[1]) * np.exp(-(FitEp-Fitm0)*tf[0])*
                       (p[0] + p[1]*np.exp(-FitDEp*tf[0])+ p[2]*np.exp(-FitDm*(tf[1]-tf[0]))))
        else:
            output =  (np.sqrt(np.abs(FitA0*FitAp0))*np.exp(-Fitm0*tf[1]) * np.exp(-(FitEp-Fitm0)*tf[0])*
                       (p[0] + p[1]*np.exp(-FitDEp*tf[0])+ p[2]*np.exp(-FitDm*(tf[1]-tf[0])) +
                        p[3] *np.exp(-FitDEp*tf[0])*np.exp(-FitDm*(tf[1]-tf[0]))))
    return output


def C3MomTSFFDer(tf,p):
    FitA0 = tf[2]
    # FitA1 = tf[3]
    FitAp0 = tf[6]
    # FitAp1 = tf[7]
    Fitm0 = np.exp(tf[4])
    FitDm = np.exp(tf[5])
    FitEp = np.exp(tf[8])
    FitDEp = np.exp(tf[9])
    unitt = [0]*len(tf[1])
    if all(FitEp == Fitm0) and all(FitA0 == FitAp0):
        if all([tf[1][0]==itf for itf in tf[1]]):
            pder = [FitA0*np.exp(-Fitm0*tf[1])]
            pder.append(pder[0]*(np.exp(-FitDm*tf[0])+np.exp(-FitDm*(tf[1]-tf[0]))))
            pder.append(unitt)
            pder.append(unitt)
        else:
            pder = [FitA0*np.exp(-Fitm0*tf[1])]
            pder.append(pder[0]*(np.exp(-FitDm*tf[0])+np.exp(-FitDm*(tf[1]-tf[0]))))
            pder.append(unitt)
            pder.append(pder[0]*np.exp(-FitDm*tf[1]))
    else:
        if all([tf[1][0]==itf for itf in tf[1]]):
            pder = [np.sqrt(np.abs(FitA0*FitAp0))*np.exp(-Fitm0*tf[1]) * np.exp(-(FitEp-Fitm0)*tf[0])]
            pder.append(pder[0]*np.exp(-FitDEp*tf[0]))
            pder.append(pder[0]*np.exp(-FitDm*(tf[1]-tf[0])))
            pder.append(unitt)
        else:
            pder = [np.sqrt(np.abs(FitA0*FitAp0))*np.exp(-Fitm0*tf[1]) * np.exp(-(FitEp-Fitm0)*tf[0])]
            pder.append(pder[0]*np.exp(-FitDEp*tf[0]))
            pder.append(pder[0]*np.exp(-FitDm*(tf[1]-tf[0])))
            pder.append(pder[0]*np.exp(-FitDEp*tf[0])*np.exp(-FitDm*(tf[1]-tf[0])))
    return pder



## bool is whether it vanishes in chiral limit
## this paper  https://arxiv.org/pdf/0911.3981.pdf sais constant + log + linear
F3ChiralLimFun = {}
F3ChiralLimFun[(True,False,False,False)] = (SquaredFitFun_x0,1)
F3ChiralLimFun[(False,False,False,False)] = (SquaredFitFun,2)
F3ChiralLimFun[(True,True,False,False)] = (SqPlLogFitFun_x0,2)
F3ChiralLimFun[(False,True,False,False)] = (SqPlLogFitFun,3)

F3ChiralLimFun[(True,False,True,False)] = (ContFitFun_x0,2)
F3ChiralLimFun[(False,False,True,False)] = (ContFitFun,3)
F3ChiralLimFun[(True,True,True,False)] = (ContPlLogFitFun_x0_dN,3)
F3ChiralLimFun[(False,True,True,False)] = (ContPlLogFitFun,4)

F3ChiralLimFun[(True,False,True,True)] = (ContFitFun_x0_mix,3)
F3ChiralLimFun[(False,False,True,True)] = (ContFitFun_mix,4)
F3ChiralLimFun[(True,True,True,True)] = (ContPlLogFitFun_x0_mix,4)
F3ChiralLimFun[(False,True,True,True)] = (ContPlLogFitFun_mix,5)


F3LatSpaceFun = {}
F3LatSpaceFun[(False,False,False,False)] = (SquaredFitFun,2)

F3LatSpaceFun[(True,False,False,False)] = (ContFitFun,3)
F3LatSpaceFun[(True,True,False,False)] = (ContFitFun_mix,4)
F3LatSpaceFun[(True,False,True,False)] = (ContFitFun_x0,2)
F3LatSpaceFun[(True,True,True,False)] = (ContFitFun_x0_mix,3)

F3LatSpaceFun[(True,False,False,True)] = (ContPlLogFitFun,4)
F3LatSpaceFun[(True,True,False,True)] = (ContPlLogFitFun_mix,5)
F3LatSpaceFun[(True,False,True,True)] = (ContPlLogFitFun_x0_dN,3)
F3LatSpaceFun[(True,True,True,True)] = (ContPlLogFitFun_x0_mix,4)

# F3ChiralLimFun[True] = (LinearFitFun_x0,1)
# F3ChiralLimFun[False] = (LinearFitFun,2)
# F3ChiralLimFun[True] = (LogFitFun_x0,1)
# F3ChiralLimFun[False] = (LogFitFun,2)
