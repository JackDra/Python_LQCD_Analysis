#!/usr/bin/env python
import numpy as np
from copy import deepcopy,copy
from Params import myeps


def GammaTOChroma(Opp):
    if isinstance(Opp,int):
        Opp = 'g'+str(Opp)
    elif isinstance(Opp,(list,tuple,np.ndarray)):
        if len(Opp) == 2:
            Opp = 'g'+str(Opp[0])+'g'+str(Opp[1])
        elif len(Opp) == 3:
            Opp = np.log2(15-np.sum(2**(np.array(Opp)-1)))+1
            Opp = 'g5g'+str(int(Opp))
        else:
            raise IOError('gamma list too long '+str(len(Opp)))

    if Opp == 'Pion':
        return 15
    if len(Opp.split('g')[1:]) > 2:
        return GammaTOChroma(list(map(int,Opp.split('g')[1:])))
    for i in [1,2,3,4]:
        for j in [1,2,3,4]:
            if 'g'+str(i)+'g'+str(j) in Opp:
                return 2**(i-1) + 2**(j-1)
    for i in [1,2,3,4]:
        if 'g'+str(i)+'g5' in Opp or 'g5g'+str(i) in Opp:
            return 15-2**(i-1)
    for i in [1,2,3,4]:
        if 'g'+str(i) in Opp:
            return 2**(i-1)
    if 'g5' in Opp:
        return 15
    elif 'I' in Opp:
        return 0
    else:
        return -1



def ChromaTOGamma(Opp):
    Opp = int(Opp)
    if Opp == 0:
        return 'I'
    elif Opp in [1,2,4,8]:
        return 'g'+str(int(np.log2(Opp)+1))
    elif Opp in [3,5,9]:
        return 'g1g'+str(int(np.log2(Opp-1)+1))
    elif Opp in [6,10]:
        return 'g2g'+str(int(np.log2(Opp-2)+1))
    elif Opp == 12:
        return 'g3g4'
    elif Opp == 15:
        return 'g5'
    elif Opp in [7,11,13,14]:
        return 'g5g'+str(int(np.log2(15-Opp)+1))


mattype = np.complex

class GammaMat(object):

    """
    Class with basic gamma matricies.

    Initialise with representation.

    resulting gamma matricies are in np.matrix format

    """

    def __init__(self,Rep='Sakurai'):

        self.Rep = Rep
        ## sets up gamma matricies
        self.SetRep(Rep)

    def SetRep(self,Rep):
        if 'sakura' in Rep or 'Sakurai' in Rep:
            self.g0 = self.g4 = np.matrix([[  1,  0,  0,  0],
                                           [  0,  1,  0,  0],
                                           [  0,  0, -1,  0],
                                           [  0,  0,  0, -1]],mattype)

            self.g1 = np.matrix([[  0,  0,  0,-1j],
                                 [  0,  0,-1j,  0],
                                 [  0, 1j,  0,  0],
                                 [ 1j,  0,  0,  0]],mattype)

            self.g2 = np.matrix([[  0,  0,  0, -1],
                                 [  0,  0,  1,  0],
                                 [  0,  1,  0,  0],
                                 [ -1,  0,  0,  0]],mattype)

            self.g3 = np.matrix([[  0,  0,-1j,  0],
                                 [  0,  0,  0, 1j],
                                 [ 1j,  0,  0,  0],
                                 [  0,-1j,  0,  0]],mattype)

            self.g5 = np.matrix([[  0,  0, -1,  0],
                                 [  0,  0,  0, -1],
                                 [ -1,  0,  0,  0],
                                 [  0, -1,  0,  0]],mattype)
            self.g_mu_Tsign = [+1,-1,+1,-1,+1,+1]
        ## TODO, needs change for whole code tho
        # elif 'dirac' in Rep or 'Dirac' in Rep:
        else:
            raise IOError(Rep+' not a recognised representation for gamma matricies')

        ## include gamma 5 here?
        self.g_mu = [self.g0,self.g1,self.g2,self.g3,self.g4]
        ## is this right?
        self.sigma_mu_nu = [[np.zeros((4,4),mattype),-1j*self.g0*self.g1,-1j*self.g0*self.g2,-1j*self.g0*self.g3,np.zeros((4,4),mattype)],
                       [-1j*self.g1*self.g0,np.zeros((4,4),mattype),-1j*self.g1*self.g2,-1j*self.g1*self.g3,-1j*self.g1*self.g4],
                       [-1j*self.g2*self.g0,-1j*self.g2*self.g1,np.zeros((4,4),mattype),-1j*self.g2*self.g3,-1j*self.g2*self.g4],
                       [-1j*self.g3*self.g0,-1j*self.g3*self.g1,-1j*self.g3*self.g2,np.zeros((4,4),mattype),-1j*self.g3*self.g4],
                       [np.zeros((4,4),mattype),-1j*self.g4*self.g1,-1j*self.g4*self.g2,-1j*self.g4*self.g3,np.zeros((4,4),mattype)]]
        self.sigma_5_nu = [self.g5*self.g0,self.g5*self.g1,self.g5*self.g2,self.g5*self.g3,self.g5*self.g4]
        self.sigma_nu_5 = [self.g0*self.g5,self.g1*self.g5,self.g2*self.g5,self.g3*self.g5,self.g4*self.g5]

        self.G_unpol = (np.eye(4,dtype=mattype) + self.g4) / 2.
        self.G_pol = -1j*self.G_unpol* self.g5 * self.g3

        ## if you like dictionaries
        self.GammaDict = {}
        for ig in range(0,5):
            for jg in range(0,5):
                self.GammaDict['sig'+str(ig)+str(jg)] = self.sigma_mu_nu[ig][jg]
            self.GammaDict['g'+str(ig)] = self.g_mu[ig]
            self.GammaDict['sig5'+str(ig)] = self.sigma_5_nu[ig]
            self.GammaDict['sig'+str(ig)+'5'] = self.sigma_nu_5[ig]
        self.GammaDict['g5'] = self.g5
        self.GammaDict['P4'] = self.G_unpol
        self.GammaDict['P3'] = self.G_pol
        self.GammaDict['I'] = np.matrix(np.eye(4,dtype=mattype), copy=False)
        self.GammaDict['neg'] = -self.GammaDict['I']
        self.GammaDict['cmplx'] = 1.0j*self.GammaDict['I']

    def __getitem__(self,ikey):
        if isinstance(ikey,list) or isinstance(ikey,np.ndarray) or isinstance(ikey,tuple):
            output = self.GammaDict[ikey[0]]
            for jkey in ikey[1:]:
                if jkey in self.GammaDict:
                    output = output *self.GammaDict[jkey]
            return output
        else:
            return self.GammaDict[ikey]

    def SigMult(self,isig,jsig,thing):
        return -1j*self.GammaMult(isig,self.GammaMult(jsig,thing))

    ##TODO, go nuts with the predefined functions

    def ProjMult(self,pindex,thing):
        if pindex == 4 or pindex == 0:
            return (thing + self.GammaMult(4,thing))/2.
        else:
            return -1j * self.ProjMult(4,self.GammaMult(5,self.GammaMult(pindex,thing) ))

    ## either pass in string corresponding to above, or anything else
    ## this is for 1 gamma matrix, multiple gamma matricies to come!
    ## this is an efficiency thing btw, probably not needed but I did it before
    ## Infact, I havn't bothered to implement these efficiency routines in the rest of the code
    ## I guess TODO if you want...
    def MultGamma(self,thing,rgindex,T=False):
        if isinstance(thing,np.matrix):
            output = deepcopy(thing)
            if T: output = self.g_mu_Tsign[rgindex]*output
            if rgindex == 0 or rgindex == 4:
                output[:,(2,3)] = -output[:,(2,3)]
            elif rgindex == 1:
                output = np.rot90(output).T * 1j
                output[:,(2,3)] = -output[:,(2,3)]
            elif rgindex == 2:
                output = np.rot90(output).T
                output[:,(0,3)] = -output[:,(0,3)]
            elif rgindex == 3:
                output = np.roll(output,2,1)*1j
                output[:,(1,2)] = -output[:,(1,2)]
            elif rgindex == 5:
                output = -np.roll(output,2,1)
            return output
        else:
            if T: return output*self.g_mu[rgindex].T
            else: return output*self.g_mu[rgindex]


    def GammaMult(self,lgindex,thing):
        if isinstance(thing,np.matrix):
            return self.MultGamma(thing.T,lgindex,T=True).T
        else:
            return self.g_mu[lgindex]*thing

    ## default is (energy, px, py ,pz)
    ## can use (px,py,pz,energy) as well if Eindex = 4
    def Slashed(self,FVecObj,Eindex=0):
        if not isinstance(FVecObj,np.ndarray): raise IOError('Slashed requires numpy array to do element-wise multiplication')
        if Eindex == 0:
            return np.sum([ig * FVecObj for ig in self.g_mu[:-1]],0)
        elif Eindex == 4:
            return np.sum([ig * FVecObj for ig in self.g_mu[1:]],0)
        else:
            raise IOError('Eindex must be 0 or 4')


class CorrSpinTrace(object):

    """

    has routines which does correlator spin traces, for use in the form factor decomposition

    something proportional to: trace( Projector * (pp_slashed +m )* Opp *(pp_slashed + m) )
    Rodgers paper has good definitions for this :D

    """

    def __init__(self,Rep='Sakurai'):

        self.Rep = Rep
        ## sets up gamma matricies class for use
        self.Gammas = GammaMat(Rep=Rep)
        self.CurrFFs = {'Scalar'   : self.ScalarFF,
                        'Vector'   : self.VectorFF,
                        'GeGm'     : self.GeGmFF,
                        'VectorTop': self.VectorFFTop,
                        'VectorTopBNL': self.VectorFFTop,
                        'VectorTopBNL_mal': self.VectorFFTop,
                        'VectorWein': self.VectorFFTop,
                        'VectorWeinBNL': self.VectorFFTop,
                        'VectorTopNoV': self.VectorFFTop,
                        'PsScalar' : self.ScalarFF,
                        'PsVector' : self.PsVectorFF,
                        'Tensor'   : self.TensorFF}


    ## key of this class is just the form factor function
    def __getitem__(self,ikey):
        return self.CurrFFs[ikey]

    ## pmu is source 4-momentum
    ## ppmu is sink 4-momentum
    ## the energy component is just the energy, the -j is done in here
    ## Opp can be anyting you can use to __getitem__ from GammaMat (string, or list of string to multiply together etc...)
    def FFunOpp(self,Opp,pmu,ppmu,mass,Rfac=True):
        thisOpp = Opp
        if isinstance(thisOpp,str):thisOpp = self.Gammas[Opp]
        if isinstance(thisOpp,tuple) or isinstance(thisOpp,list):thisOpp = self.Gammas[Opp]
        Ep, Epp = -1.0j*pmu[0],-1.0j*ppmu[0]
        pplusm = (self.Gammas.g4 - (1j/Ep) * (pmu[1]*self.Gammas.g1 + pmu[2]*self.Gammas.g2 + pmu[3]*self.Gammas.g3) + (mass / Ep) * np.eye(4))
        pprimeplusm = (self.Gammas.g4 - (1j/Epp) * (ppmu[1]*self.Gammas.g1 + ppmu[2]*self.Gammas.g2 + ppmu[3]*self.Gammas.g3) + (mass / Epp )* np.eye(4))
        if Rfac:
            return (pplusm * thisOpp * pprimeplusm)*np.sqrt(Epp*Ep/((Epp+mass)*(Ep+mass))) * 1/4.
        else:
            return (pplusm * thisOpp * pprimeplusm) * 1/4.



    ## Form Factor decomposition needs alpha too when using CP voilating operator!!
    def FFunTopOpp(self,Opp,pmu,ppmu,mass,alpha,Rfac=True):
        thisOpp = Opp
        if isinstance(thisOpp,str):thisOpp = self.Gammas[Opp]
        if isinstance(thisOpp,tuple) or isinstance(thisOpp,list):thisOpp = self.Gammas[Opp]
        p,pp = pmu,ppmu
        m = mass
        Ep, Epp = -1.0j*p[0],-1.0j*pp[0]
        pplusm = (self.Gammas.g4 - (1j/Ep) * (p[1]*self.Gammas.g1 + p[2]*self.Gammas.g2 + p[3]*self.Gammas.g3) + (m / Ep) * np.eye(4))
        pprimeplusm = (self.Gammas.g4 - (1j/Epp) * (pp[1]*self.Gammas.g1 + pp[2]*self.Gammas.g2 + pp[3]*self.Gammas.g3) + (m / Epp )* np.eye(4))
        g5facp = (2.0*alpha*m/Ep)*self.Gammas.g5
        g5facpprime = (2.0*alpha*m/Epp)*self.Gammas.g5
        if Rfac:
            return ((pplusm * thisOpp * g5facpprime)*np.sqrt(Epp*Ep/((Epp+m)*(Ep+m))) * 1/4.,
                    (g5facp * thisOpp * pprimeplusm)*np.sqrt(Epp*Ep/((Epp+m)*(Ep+m))) * 1/4.)
        else:
            return ((pplusm * thisOpp * g5facpprime) * 1/4., (g5facp * thisOpp * pprimeplusm) * 1/4.)



    ## traced with projector of FFunOpp
    def TracedFFun(self,Opp,pmu,ppmu,mass,Rfac=True):
        thisProj,thisOpp = self.GetProjGamma(*Opp)
        return np.trace(thisProj * self.FFunOpp(thisOpp,pmu,ppmu,mass,Rfac=Rfac))

    ## see https://arxiv.org/pdf/1507.02343v2.pdf for definitin of CP odd form factor decomposition.
    def TracedFFunTop(self,Opp,pmu,ppmu,mass,alpha,Rfac=True):
        thisProj,thisOpp = self.GetProjGamma(*Opp)
        leftval,rightval = self.FFunTopOpp(thisOpp,pmu,ppmu,mass,alpha,Rfac=Rfac)
        return np.trace(thisProj * leftval)+np.trace(thisProj * rightval)



    ## these functions are expecting a list of gamma matricies
    def GetProjGamma(self,*Opp):
        return self.GetProj(*Opp),self.GetGamma(*Opp)

    ## only returns first projector found, if none found, returns identity
    def GetProj(self,*thisProj):
        for iProj in thisProj:
            if 'P' in iProj:
                return self.Gammas[iProj]
        return self.Gammas['I']

    def GetGamma(self,*Opp):
        outlist = []
        for iOpp in Opp:
            if 'P' not in iOpp:
                outlist.append(iOpp)
        if len(outlist) > 0:
            return self.Gammas[outlist]
        else:
            return self.Gammas['I']

    def CreateEs(self,momvec,thismass,curr=False):
        momvecsqrd = sum([iq**2 for iq in momvec])
        return np.sqrt(thismass**2+momvecsqrd)

    def Create4Mom(self,thisqvec,thisppvec,thismass):
        thispvec = [ipp-iq for ipp,iq in zip(thisppvec,thisqvec)]
        if isinstance(thisppvec,np.ndarray): thisppvec = thisppvec.tolist()
        if isinstance(thispvec,np.ndarray): thispvec = thispvec.tolist()
        thisp = [1.0j*self.CreateEs(thispvec,thismass)] + thispvec + [1.0j*self.CreateEs(thispvec,thismass)]
        thispp = [1.0j*self.CreateEs(thisppvec,thismass)] + thisppvec + [1.0j*self.CreateEs(thisppvec,thismass)]
        thisq = [ipp-ip for ipp,ip in zip(thispp,thisp)]
        return thisp,thisq,thispp


    def ScalarFF(self,opp,thisqvec,thisppvec,thismass,Rfac=True,alpha=1.0):
        thisp,thisq,thispp = self.Create4Mom(thisqvec,thisppvec,thismass)
        term1 = self.TracedFFun(opp,thisp,thispp,thismass,Rfac=Rfac)
        rcheck,ccheck = abs(term1.real)<myeps,abs(term1.imag)<myeps
        return [term1],not rcheck, not ccheck

    def VectorFF(self,opp,thisqvec,thisppvec,thismass,Rfac=True,PadF3=False,alpha=1.0):
        thisp,thisq,thispp = self.Create4Mom(thisqvec,thisppvec,thismass)
        thisopp,iscmplx,isneg = self.FormatOpp(opp)
        # print 'calculating form factor coeffs for ' ,thisopp
        term1 = self.TracedFFun(opp,thisp,thispp,thismass,Rfac=Rfac)
        term2 = 0.0j
        for i in [1,2,3,4]:
            if str(i) not in thisopp[-1]:
                if iscmplx:
                    # print 'calculating form factor coeffs for ' ,thisopp[:-1]+('g'+str(i),),'xq'
                    term2 += self.TracedFFun(thisopp+('g'+str(i)),thisp,thispp,thismass,Rfac=Rfac)*thisq[i]
                else:
                    ## sigma_mu_nu * q_mu = -i * g_mu * g_nu * q_nu
                    # print 'calculating form factor coeffs for ' ,thisopp+('g'+str(i),'cmplx','neg'),'xq'
                    term2 += self.TracedFFun(thisopp+('g'+str(i),'cmplx','neg'),thisp,thispp,thismass,Rfac=Rfac)*thisq[i]
        term2 = term2/(2.0*thismass)
        rcheck,ccheck = abs(term1.real)<myeps and abs(term2.real)<myeps,abs(term1.imag)<myeps and abs(term2.imag)<myeps
        if PadF3:
            return [term1,term2,0.0],not rcheck, not ccheck
        else:
            return [term1,term2],not rcheck, not ccheck


    def GeGmFF(self,opp,thisqvec,thisppvec,thismass,Rfac=True,PadF3=False,alpha=1.0):
        termlist,rcheck,ccheck = self.VectorFF(opp,thisqvec,thisppvec,thismass,Rfac=Rfac,PadF3=PadF3,alpha=alpha)
        thisp,thisq,thispp = self.Create4Mom(thisqvec,thisppvec,thismass)
        term1,term2 = termlist[:2]
        ## QM = Q^2/4m^2
        ## Ge = F1 - QM F2
        ## Gm = F1 + F2
        ## reranged
        ## F1 = (QM Gm + Ge) / (QM + 1)
        ## F2 = (Gm - Ge) / (QM + 1)

        ## Ratio = term1 * F1 + term2 * F2
        ## Ratio = [( QM*term1 +term2) Gm + (term1-term2) Ge]/(QM+1)
        QM = sum(np.array(thisq)**2)/(4*thismass**2)
        QMP1 = QM + 1.
        GeGmterm1 = (term1-term2)/QMP1
        if abs(term2.real) < myeps:
            GeGmterm2 = term2
        else:
            # if opp == ('P4','g4'):
            #     print 'debug'
            #     print opp,QM , term1, term2
            #     print
            GeGmterm2 = (QM*term1+term2)/QMP1
        rcheck,ccheck = abs(GeGmterm1.real)<myeps and abs(GeGmterm2.real)<myeps,abs(GeGmterm1.imag)<myeps and abs(GeGmterm2.imag)<myeps
        GeGmout = [GeGmterm1,GeGmterm2]
        if PadF3: GeGmout += [0.0]
        return GeGmout,not rcheck,not ccheck

    ## helper function used to pull out usefull stuff from opp
    def FormatOpp(self,opp,combProjGamma=True):
        thisProj,thisopp,iscmplx,isneg = (),(),False,False
        for iopp in opp:
            if iopp not in ['Top','Wein','cmplx','neg']:
                if 'P' in iopp:
                    thisProj = (iopp,)
                else:
                    thisopp += (iopp,)
            elif iopp == 'cmplx':
                iscmplx = True
            elif iopp == 'neg':
                isneg = True
        if combProjGamma:
            return thisProj+thisopp,iscmplx,isneg
        else:
            return thisProj,thisopp,iscmplx,isneg

    def VectorFFTop(self,opp,thisqvec,thisppvec,thismass,Rfac=True,alpha=1.0):
        thisp,thisq,thispp = self.Create4Mom(thisqvec,thisppvec,thismass)
        if 'Top' in opp or 'Wein' in opp:
            thisopp,iscmplx,isneg = self.FormatOpp(opp)
            term1 = self.TracedFFunTop(opp,thisp,thispp,thismass,alpha,Rfac=Rfac)
            term2 = 0.0j
            for i in [1,2,3,4]:
                if str(i) not in thisopp[-1]:
                    if iscmplx:
                        # print 'calculating form factor coeffs for ' ,thisopp[:-1]+('g'+str(i),),'xq'
                        term2 += self.TracedFFunTop(thisopp+('g'+str(i),),thisp,thispp,thismass,alpha,Rfac=Rfac)*thisq[i]
                    else:
                        ## sigma_mu_nu * q_mu = -i * g_mu * g_nu * q_nu
                        # print 'calculating form factor coeffs for ' ,thisopp+('g'+str(i),'cmplx','neg'),'xq'
                        term2 += self.TracedFFunTop(thisopp+('g'+str(i),'cmplx','neg'),thisp,thispp,thismass,alpha,Rfac=Rfac)*thisq[i]
            term2 = term2/(2.0*thismass)
            term3 = 0.0j
            for i in [1,2,3,4]:
                if str(i) not in thisopp[-1]:
                    if iscmplx:
                        ## sigma_mu_nu * q_mu = -i * g_mu * g_nu * q_nu
                        term3 += self.TracedFFun(thisopp+('g'+str(i),'g5'),thisp,thispp,thismass,Rfac=Rfac)*thisq[i]
                    else:
                        term3 += self.TracedFFun(thisopp+('g'+str(i),'g5','cmplx','neg'),thisp,thispp,thismass,Rfac=Rfac)*thisq[i]
            term3 = term3/(2.0*thismass)
            rcheck,ccheck = (abs(term1.real)<myeps and abs(term2.real)<myeps and abs(term3.real)<myeps,
                             abs(term1.imag)<myeps and abs(term2.imag)<myeps  and abs(term3.imag)<myeps)
            return [term1,term2,term3],not rcheck, not ccheck
        else:
            return self.VectorFF(opp,thisqvec,thisppvec,thismass,Rfac=Rfac,PadF3=True)

    def PsVectorFF(self,opp,thisqvec,thisppvec,thismass,Rfac=True,alpha=1.0):
        thisp,thisq,thispp = self.Create4Mom(thisqvec,thisppvec,thismass)
        term1 = self.TracedFFun(opp,thisp,thispp,thismass,Rfac=Rfac)
        oppNoG = list(copy(opp))
        oppNoG.remove('g5')
        for iopp in oppNoG:
            if 'g' in iopp:
                index1 = int(iopp[-1])
                break
        term2 = (1.0j*self.TracedFFun(oppNoG,thisp,thispp,thismass,Rfac=Rfac)*thisq[index1])/(2.0*thismass)
        rcheck,ccheck = abs(term1.real)<myeps and abs(term2.real)<myeps,abs(term1.imag)<myeps and abs(term2.imag)<myeps
        return [term1,term2],not rcheck, not ccheck


    def TensorFF(self,opp,thisqvec,thisppvec,thismass,Rfac=True,alpha=1.0):
        thisp,thisq,thispp = self.Create4Mom(thisqvec,thisppvec,thismass)
        ## P = pp + p = 2*pp -q
        thisP = [ip+ipp for ip,ipp in zip(thisp,thispp)]

        Proj,gammalist,iscmplx,isneg = self.FormatOpp(opp,combProjGamma=False)
        coeff = 1.0
        if iscmplx: coeff = 1.0j
        if isneg: coeff = -coeff
        index1,index2 = int(gammalist[0][-1]),int(gammalist[1][-1])
        ## sigma_mu_nu * q_mu = -i * g_mu * g_nu * q_nu
        ## term1 = i sigma_mu_nu H_T(Q^2)
        ## i and -i cancel when sub in sigma_mu_nu
        ## TODO, not sure of the sign of this!!!!
        ## most papers are in Minkowski space me thinks
        term1 = coeff*self.TracedFFun(Proj+gammalist,thisp,thispp,thismass,Rfac=Rfac)

        # print Proj,gammalist
        # print ''
        # print 'term1'
        # print self.TracedFFun(Proj+gammalist,thisp,thispp,thismass,Rfac=Rfac)
        # print ''


        ## term2 = (gamma_mu q_nu - gamma_nu * q_mu ) /2m * E_T(Q^2)
        ## TODO: not sure of the sign of this!!
        ## most papers are in Minkowski space me thinks
        term2 = 1.0j*coeff*(self.TracedFFun(Proj+(gammalist[0],),thisp,thispp,thismass,Rfac=Rfac)*thisq[index2] -
                 self.TracedFFun(Proj+(gammalist[1],),thisp,thispp,thismass,Rfac=Rfac)*thisq[index1] )/(2.0*thismass)
        # print ''
        # print 'term2'
        # print self.TracedFFun(Proj+(gammalist[0],),thisp,thispp,thismass,Rfac=Rfac), thisq[index2]
        # print self.TracedFFun(Proj+(gammalist[0],),thisp,thispp,thismass,Rfac=Rfac)*thisq[index2]
        # print self.TracedFFun(Proj+(gammalist[1],),thisp,thispp,thismass,Rfac=Rfac) , thisq[index1]
        # print self.TracedFFun(Proj+(gammalist[1],),thisp,thispp,thismass,Rfac=Rfac)*thisq[index1]
        # print ''

        ## term3 = (P_mu q_nu - P_nu q_nu)/2m^2 tilde{H}_1
        term3 = coeff*self.TracedFFun(Proj+('I',),thisp,thispp,thismass,Rfac=Rfac)*(thisP[index1]*thisq[index2] -
                                   thisP[index2]*thisq[index1] )/(2.0*thismass**2)
        # print 'term3'
        # print self.TracedFFun(Proj+('I',),thisp,thispp,thismass,Rfac=Rfac)
        # print thisP[index1]*thisq[index2]
        # print thisP[index2]*thisq[index1]
        # print (thisP[index1]*thisq[index2] - thisP[index2]*thisq[index1] )
        # print ''
        ##Discrepancy? i think
        rcheck,ccheck = (abs(term1.real)<myeps and abs(term2.real)<myeps and abs(term3.real)<myeps,
                         abs(term1.imag)<myeps and abs(term2.imag)<myeps and abs(term3.imag)<myeps)
        return [term1,term2,term3], not rcheck,not ccheck

if __name__ == '__main__':
    gdata = GammaMat()
    print(gdata)
