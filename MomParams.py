#!/usr/bin/env python

import numpy as np
import pandas as pa
from MiscFuns import mkdir_p,ODNested
from Params import defmom2list,nxyzt,hbarc,latspace
import os
# import dill as pickle
from Params import defInfo
from FileIO import WriteXml,WritePickle,ReadPickleWrap,Construct_File_Object
from collections import OrderedDict
from XmlFormatting import AvgStdToFormat


"""

Updated version of MomParams.py, now contains a class for specific momenta lists.

"""

physcial_pion_mass = 139.57061 ##MeV  from http://pdg.lbl.gov/2017/tables/rpp2017-tab-mesons-light.pdf
physcial_pion_mass_err = 0.00024 ##MeV  from http://pdg.lbl.gov/2017/tables/rpp2017-tab-mesons-light.pdf



class LatticeParameters(object):
    """

    Contains routines and parameters for momentum

    """
    def Construct_LParam_File(this_file):
        return Construct_File_Object(this_file,LatticeParameters)

    def __init__(self,mom2list=defmom2list,Info=defInfo,name = ''):
        self.mom2list = np.array(mom2list)
        self.Info = Info
        self.mom2max = np.max(self.mom2list)
        self.mom2min = np.min(self.mom2list)


        ## Box verticies in each dimension (nx,ny,nz,nt)
        if 'nxyzt' in list(Info.keys()): self.nxyzt = Info['nxyzt']
        else: self.nxyzt = nxyzt
        if len(self.nxyzt) == 2:
            self.nxyzt = [self.nxyzt[0],self.nxyzt[0],self.nxyzt[0],self.nxyzt[-1]]

        ## this is needed
        if 'outdir' in list(Info.keys()):
            self.outputdir = Info['outdir']
        else:
            print(Info)
            raise IOError("output directory is needed in information dictionary, (Info['outdir'] = ...) ")



        self.hbarc = hbarc
        self.nxyz = self.nxyzt[:-1]
        self.nt = self.nxyzt[-1]
        # self.nx,self.ny,self.nz = self.nxyz
        if self.nxyz != [self.nxyz[0],self.nxyz[0],self.nxyz[0]]:
            raise IOError('code only works for square spatial lattices')
        else:
            self.nxyz = self.nxyz[0]

        self.dim_label = 'RC'+str(self.nxyz)+'x'+str(self.nt)

        self.qunit = (2.0*np.pi)/float(self.nxyz)

        ## lattice spacing in fm
        if 'a' in list(Info.keys()): self.latspace = Info['a']
        else: self.latspace = latspace
        self.hbarcdivlat = hbarc/self.latspace ## in GeV
        self.qunitPhys = self.qunit*self.hbarcdivlat ## in GeV
        self.Lxyz = self.latspace*self.nxyz ## in fm
        self.Lt = self.latspace*self.nt ## FIX UNITS TODO SLD:HGSGH

        ## kappa up down (in integer form)
        if 'kud' in list(Info.keys()): self.kud = 'kud'+str(Info['kud'])
        else: self.kud = ''

        ## kappa strange (in integer form)
        if 'ks' in list(Info.keys()): self.ks = 'ks'+str(Info['ks'])
        else: self.ks = ''

        self.kappafolder = self.dim_label+'Kud0'+self.kud.replace('kud','')+'Ks0'+self.ks.replace('ks','')

        if self.kud in 'kud1370000' and self.ks in 'ks1364000' and 'RC32x64' in self.dim_label:
            self.NucleonMass = 0.7277 ## Lat Units
            self.NucleonMassErr = 0.0022 ## Lat Units
            self.PionMass = 0.3213 ## Lat Units
            self.PionMassErr = 0.00015 ## Lat Units
            self.UDMass = 0.030753 # am_ud in MSbar
            self.UDMassErr = 0.000110 # am_ud in MSbar
            self.SMass = 0.047142 # am_s in MSbar
            self.SMassErr = 0.0000110 # am_s in MSbar
        elif self.kud == 'kud1375400' and self.ks == 'ks1364000'and 'RC32x64' in self.dim_label:
            self.NucleonMass = 0.5584 ## Lat Units
            self.NucleonMassErr = 0.0053 ## Lat Units
            self.PionMass = 0.1883 ## Lat Units
            self.PionMassErr = 0.0003 ## Lat Units
            self.UDMass = 0.011028 # am_ud in MSbar
            self.UDMassErr = 0.000080 # am_ud in MSbar
            self.SMass = 0.042355 # am_s in MSbar
            self.SMassErr = 0.000079 # am_s in MSbar
        elif self.kud == 'kud1372700' and self.ks == 'ks1364000'and 'RC32x64' in self.dim_label:
            self.NucleonMass = 0.6487 ## Lat Units
            self.NucleonMassErr = 0.0056 ## Lat Units
            self.PionMass = 0.2609 ## Lat Units
            self.PionMassErr = 0.00015 ## Lat Units
            self.UDMass = 0.020834 # am_ud in MSbar
            self.UDMassErr = 0.000066 # am_ud in MSbar
            self.SMass = 0.044674 # am_s in MSbar
            self.SMassErr = 0.000072 # am_s in MSbar
        elif 'RC16x32' in self.dim_label:
            self.NucleonMass = float('NaN') ## Lat Units
            self.NucleonMassErr = float('NaN') ## Lat Units
            self.PionMass = 0.3940 ## Lat Units Taken from kud 13700
            self.PionMassErr = 0.00055 ## Lat Units
            self.UDMass = float('NaN') # am_ud in MSbar
            self.UDMassErr = float('NaN') # am_ud in MSbar
            self.SMass = float('NaN') # am_s in MSbar
            self.SMassErr = float('NaN') # am_s in MSbar
        elif 'RC20x40' in self.dim_label:
            self.NucleonMass = float('NaN') ## Lat Units
            self.NucleonMassErr = float('NaN') ## Lat Units
            self.PionMass = 0.3208 ## Lat Units Taken from kud 13700
            self.PionMassErr = 0.00033 ## Lat Units
            self.UDMass = float('NaN') # am_ud in MSbar
            self.UDMassErr = float('NaN') # am_ud in MSbar
            self.SMass = float('NaN') # am_s in MSbar
            self.SMassErr = float('NaN') # am_s in MSbar
        elif 'RC28x56' in self.dim_label:
            self.NucleonMass = float('NaN') ## Lat Units
            self.NucleonMassErr = float('NaN') ## Lat Units
            self.PionMass = 0.2289 ## Lat Units Taken from kud 13700
            self.PionMassErr = 0.00024 ## Lat Units
            self.UDMass = float('NaN') # am_ud in MSbar
            self.UDMassErr = float('NaN') # am_ud in MSbar
            self.SMass = float('NaN') # am_s in MSbar
            self.SMassErr = float('NaN') # am_s in MSbar
        else:
            print('ensemble not found for parameters, setting to NaN')
            self.NucleonMass = float('NaN') ## Lat Units
            self.NucleonMassErr = float('NaN') ## Lat Units
            self.PionMass = float('NaN') ## Lat Units Taken from kud 13700
            self.PionMassErr = float('NaN') ## Lat Units
            self.UDMass = float('NaN') # am_ud in MSbar
            self.UDMassErr = float('NaN') # am_ud in MSbar
            self.SMass = float('NaN') # am_s in MSbar
            self.SMassErr = float('NaN') # am_s in MSbar
        self.SetCustomName(name)


    def RemoveDup(self,thislist):
        output = []
        for ival in thislist:
            if ival not in output:
                output.append(ival)
        if isinstance(thislist,np.ndarray):
            return np.array(output)
        else:
            return output

    def Set_t0(self,this_to,cfglist=None):
        self.t0_values = this_to
        self.t0_cfglist = cfglist
        self.Write()

    def Get_t0(self):
        if not hasattr(self,'t0_values'):
            print('t0 not present in')
            print(self.PickleFile)
            self.Reload()
        if not hasattr(self,'t0_values'):
            raise EnvironmentError('t0 not previously computed')
        else:
            return self.t0_values


    def SetCustomName(self,name):
        if len(name) == 0:
            self.name = 'psqrd'+''.join(map(str,self.mom2list))
        else:
            self.name = name
        self.kappafolder = self.dim_label+'Kud0'+self.kud.replace('kud','')+'Ks0'+self.ks.replace('ks','')
        self.thisdir = self.outputdir + '/'+self.kappafolder+'/Parameters/'
        try:
            mkdir_p(self.thisdir)
        except Exception as err:
            pass

        self.HumanFile = self.thisdir+self.name+'.xml'
        self.PickleFile = self.thisdir+self.name+'.py3p'

    def FixRead(self,readdict):
        readdict['PickleFile'] = self.PickleFile
        readdict['HumanFile'] = self.HumanFile
        readdict['outputdir'] = self.outputdir
        return readdict

    def Reload(self):
        if os.path.isfile(self.PickleFile):
            loadeddict = ReadPickleWrap(self.PickleFile)
            # loadeddict = self.FixDict(loadeddict)
            self.__dict__.update(loadeddict)
        else:
            print('momparams Reload failed, cannot find file:')
            print(self.PickleFile)


    def LoadPickle(self,DefWipe=True,no_write=False):
        if os.path.isfile(self.PickleFile) and not DefWipe:
            # print 'Loading Pickle for ' , self.name
            loadeddict = ReadPickleWrap(self.PickleFile)
            # loadeddict = self.FixDict(loadeddict)
            self.__dict__.update(loadeddict)
        elif no_write:
            self.Avgmomlist()
            self.Chromamomlist()
        else:
            self.ReadAndWrite()

    def ReadAndWrite(self):
        self.Avgmomlist()
        self.Chromamomlist()
        self.Write()

    def Write(self):
        ## Write human readable file with data
        def FixDictArray(thislist,ikey):
            outDict = OrderedDict()
            for ic,ilist in enumerate(thislist):
                outDict[ikey+str(ic+1)] = ilist
            return outDict

        outDict = ODNested()
        outDict['kud'] = self.kud
        outDict['ks'] = self.ks
        outDict['q_unit'] = self.qunit
        outDict['q_unit_Phys_GeV'] = self.qunitPhys
        outDict['lattice_space_fm'] = self.latspace
        outDict['nt'] = self.nt
        outDict['nxyz'] = self.nxyz
        outDict['Lt_fm'] = self.Lt
        outDict['Lxyz_fm'] = self.Lxyz
        if hasattr(self,'t0_values'):
            outDict['t0_flow'] = AvgStdToFormat(self.t0_values['boot'].iloc[0].Avg,self.t0_values['boot'].iloc[0].Std)
        else:
            if os.path.isfile(self.PickleFile):
                file_data = ReadPickleWrap(self.PickleFile)
                if 't0_values' in file_data.keys():
                    self.t0_values = file_data['t0_values']
                    self.t0_cfglist = file_data['t0_cfglist']
                    outDict['t0_flow'] = AvgStdToFormat(self.t0_values['boot'].iloc[0].Avg,self.t0_values['boot'].iloc[0].Std)



        outDict['kappafolder'] = self.kappafolder
        outDict['PionMass'] = self.PionMass
        outDict['NucleonMass'] = self.NucleonMass
        outDict['PickleFile'] = self.PickleFile
        outDict['outputdir'] = self.outputdir
        outDict['momentum_list_average'] = FixDictArray(self.momformlistAvg,'pAvg_')
        outDict['momentum_list'] = FixDictArray(self.momformlist,'p_')



        WriteXml(self.HumanFile,{'Results':outDict})
        ## pickles rest of data for reading
        WritePickle(self.PickleFile,self.__dict__,Bak=False)

    def GetNucleonMass(self,Phys=True,MeV=False):
        if Phys:
            if MeV:
                return self.NucleonMass*self.hbarcdivlat*1000
            else:
                return self.NucleonMass*self.hbarcdivlat
        else:
            return self.NucleonMass

    def GetNucleonMassErr(self,Phys=True,MeV=False):
        if Phys:
            if MeV:
                return self.NucleonMassErr*self.hbarcdivlat*1000
            else:
                return self.NucleonMassErr*self.hbarcdivlat
        else:
            return self.NucleonMassErr


    def GetPionMass(self,Phys=True,MeV=False):
        if Phys:
            if MeV:
                return self.PionMass*self.hbarcdivlat*1000
            else:
                return self.PionMass*self.hbarcdivlat
        else:
            return self.PionMass

    def GetPionMassErr(self,Phys=True,MeV=False):
        if Phys:
            if MeV:
                return self.PionMassErr*self.hbarcdivlat*1000
            else:
                return self.PionMassErr*self.hbarcdivlat
        else:
            return self.PionMassErr

    def GetUDMass(self,Phys=True,MeV=False):
        if Phys:
            if MeV:
                return self.UDMass*self.hbarcdivlat*1000
            else:
                return self.UDMass*self.hbarcdivlat
        else:
            return self.UDMass

    def GetUDMassErr(self,Phys=True,MeV=False):
        if Phys:
            if MeV:
                return self.UDMassErr*self.hbarcdivlat*1000
            else:
                return self.UDMassErr*self.hbarcdivlat
        else:
            return self.UDMassErr

    def GetSMass(self,Phys=True,MeV=False):
        if Phys:
            if MeV:
                return self.SMass*self.hbarcdivlat*1000
            else:
                return self.SMass*self.hbarcdivlat
        else:
            return self.SMass

    def GetSMassErr(self,Phys=True,MeV=False):
        if Phys:
            if MeV:
                return self.SMassErr*self.hbarcdivlat*1000
            else:
                return self.SMassErr*self.hbarcdivlat
        else:
            return self.SMassErr




    def GetNucleonMassLab(self,Phys=True,MeV=False):
        if Phys:
            if MeV:
                return r'\quad m_{{N}}={:.3f} MeV '.format(self.GetNucleonMass(Phys=Phys,MeV=MeV))
            else:
                return r'\quad m_{{N}}={:.3f} GeV '.format(self.GetNucleonMass(Phys=Phys,MeV=MeV))
        else:
            return r'\quad m_{{N}}={:.3f} '.format(self.GetNucleonMass(Phys=Phys,MeV=MeV))

    def GetPionMassLab(self,Phys=True,MeV=False):
        if Phys:
            if MeV:
                return r'\quad m_{{\pi}}={:.3f} MeV '.format(self.GetPionMass(Phys=Phys,MeV=MeV))
            else:
                return r'\quad m_{{\pi}}={:.3f} GeV '.format(self.GetPionMass(Phys=Phys,MeV=MeV))
        else:
            return r'\quad m_{{\pi}}={:.3f} '.format(self.GetPionMass(Phys=Phys,MeV=MeV))


    def GetLatSpaceLab(self):
        return r'\quad a={:.3f} fm '.format(self.latspace)


    def GetBoxSizeLab(self):
        return r'\quad La={:.3f} fm '.format(self.Lxyz)

    def pvecTOpstr(self,ip):
        return ' '.join(map(str,ip))

    def pstrTOpvec(self,ip,actual=False):
        if actual:
            return [np.float(iq) * self.qunit for iq in ip.split(' ')[-3:]]
        else:
            return list(map(int,ip.split(' ')[-3:]))

    def pstrTOpform(self,ip,pref=''):
        return pref+ip.replace(' ','')

    def pformTOpstr(self,ip,nopref=False):
        if nopref:
            return ' '.join(ip.replace('p','').replace('q','')).replace('- ','-')
        else:
            return ' '.join(ip).replace('- ','-')

    def pformTOpvec(self,ip,actual=False):
        pout = ' '.join(ip).replace('- ','-').split(' ')
        if actual:
            return [np.float(iq)*self.qunit for iq in pout[-3:]]
        else:
            return list(map(int,pout[-3:]))

    def pvecTOpform(self,ip,pref=''):
        return pref+''.join(map(str,ip))

    ## use these instead of the above ones
    ## these work for any type
    def TOpstr(self,ip):
        if isinstance(ip,str):
            if ' ' in ip:
                return ip
            else:
                return self.pformTOpstr(ip)
        else:
            return self.pvecTOpstr(ip)

    def TOpvec(self,ip,actual=False):
        if isinstance(ip,str):
            if ' ' in ip:
                return self.pstrTOpvec(ip,actual=actual)
            else:
                return self.pformTOpvec(ip,actual=actual)
        else:
            if actual:
                return ip*self.qunit
            else:
                return ip


    def TOpform(self,ip,pref=''):
        if isinstance(ip,str):
            if ' ' in ip:
                return self.pstrTOpform(ip,pref)
            else:
                return pref+ip
        else:
            return self.pvecTOpform(ip,pref)

    def TOip(self,ip):
        if isinstance(ip,str):
            if ' ' in ip:
                return self.momstrlist.tolist().index(ip)
            else:
                return self.momformlist.tolist().index(ip)
        else:
            return self.momlist.tolist().index(ip)


    ## acutal=True corresponds to the actual p^2 ( in lattice units )
    ## actual=False corresponds to the integer amount of p^{2} [ in units of self.punit**2 = (2pi/L)**2]
    def Getpsqrd(self,ip,actual=False):
        thisip = list(self.TOpvec(ip))
        if actual:
            return self.Buffpsqrd[self.momlist.tolist().index(thisip)]*(self.qunit**2)
        else:
            return self.Buffpsqrd[self.momlist.tolist().index(thisip)]

    def Getpsqrdform(self,ip,pre='p'):
        return pre+'sqrd'+str(self.Getpsqrd(ip))

    ## converts the momenta to the version which has been averaged
    ## usefull for reading in corresponding 2 point correlator from the 3 point correlator
    def GetAvgmom(self,ip):
        thispsqrd = self.Getpsqrd(ip)
        return self.momlistAvg[self.mom2list.tolist().index(thispsqrd)]

    ## converts the momenta to the version which has been averaged
    ## usefull for reading in corresponding 2 point correlator from the 3 point correlator
    def GetAvgmomstr(self,ip):
        return self.TOpstr(self.GetAvgmom(ip))

    ## converts the momenta to the version which has been averaged
    ## usefull for reading in corresponding 2 point correlator from the 3 point correlator
    def GetAvgmomform(self,ip,pref=''):
        return self.TOpform(self.GetAvgmom(self.stripmom(ip)),pref=pref)




    ## for making p###, pp###, q### (see other routines
    def momlistPref(self,pref):
        return [pref+imom for imom in self.momformlist ]

    def momlistAvgPref(self,pref):
        return [pref+imom for imom in self.momformlistAvg ]

    ## for making p###, pp###, q### (see other routines
    def stripmom(self,thismom):
        if isinstance(thismom,str):
            return ' '.join(thismom.replace('p','').replace('q','')).replace('- ','-')
        else:
            return thismom[-3:]



    ## note, ordering of momenta is important
    def Avgmomlist(self):
        self.momlistAvg = []
        self.momstrlistAvg = []
        self.momformlistAvg = []
        mom2listhold = []
        for iq1 in range(self.mom2max+1):
            for iq2 in range(self.mom2max+1):
                for iq3 in range(self.mom2max+1):
                    if iq1**2 + iq2**2 + iq3**2 not in self.mom2list: continue
                    if iq1**2 + iq2**2 + iq3**2 in mom2listhold: continue
                    mom2listhold.append(iq1**2 + iq2**2 + iq3**2)
                    # qlist = np.append(qlist,'q = ' + str(iq3) + ' ' + str(iq2) + ' ' + str(iq1))
                    self.momlistAvg.append((iq3 , iq2 , iq1))
                    self.momstrlistAvg.append(' '.join(map(str,(iq3 , iq2 , iq1))))
                    self.momformlistAvg.append(''.join(map(str,(iq3 , iq2 , iq1))))
        self.momlistAvg = np.array(self.momlistAvg)
        self.momstrlistAvg = np.array(self.momstrlistAvg)
        self.momformlistAvg = np.array(self.momformlistAvg)
        self.mom2list = np.array(mom2listhold)
        self.mom2max = np.max(self.mom2list)
        self.mom2min = np.min(self.mom2list)


    def GetMomList(self,list_type='form',list_order='Def'):
        if list_order == 'Def':
            if 'Avg' in list_type:
                list_order = 'squared'
            else:
                list_order = 'flat'
        if 'form' in list_type:
            if 'Avg' in list_type:
                this_list = self.momformlistAvg
            else:
                this_list = self.momformlist
        elif 'number' in list_type:
            if 'Avg' in list_type:
                this_list = self.momlistAvg
            else:
                this_list = self.momlist
        elif 'str' in list_type:
            if 'Avg' in list_type:
                this_list = self.momstrlistAvg
            else:
                this_list = self.momstrlist
        else:
            raise IOError(str(list_type)+'not recognised as momentum list type, select form, number, or str (w/wo Avg)')
        if 'flat' in list_order:
            if 'rev' in list_order:
                return this_list[::-1]
            else:
                return this_list
        elif 'squared' in list_order:
            if 'rev' in list_order:
                return self.MomSqrdSort(this_list)[::-1]
            else:
                return self.MomSqrdSort(this_list)
        else:
            raise IOError(str(list_order)+'not recognised ordering for mom list, select flat, qsqrd, (w/wo rev)')

    def MomSqrdSort(self,this_mom_list):
        return np.array(sorted(this_mom_list, key=self.Getpsqrd))



    ## note, ordering of momenta is important
    def Chromamomlist(self):
        self.momlist = []
        self.momstrlist = []
        self.momformlist = []
        self.Buffpsqrd = []
        for iq1 in range(-self.mom2max,self.mom2max+1):
            for iq2 in range(-self.mom2max,self.mom2max+1):
                for iq3 in range(-self.mom2max,self.mom2max+1):
                    if iq1**2 + iq2**2 + iq3**2 not in self.mom2list: continue
                    self.Buffpsqrd.append(iq1**2 + iq2**2 + iq3**2)
                    self.momlist.append((iq3 , iq2 , iq1))
                    self.momstrlist.append(' '.join(map(str,(iq3 , iq2 , iq1))))
                    self.momformlist.append(''.join(map(str,(iq3 , iq2 , iq1))))
        self.momlist = np.array(self.momlist)
        self.momstrlist = np.array(self.momstrlist)
        self.momformlist = np.array(self.momformlist)
        self.Buffpsqrd = np.array(self.Buffpsqrd)


    def TOLatEBoot(self,ip,mass,Disp=0.):
        value = (mass/Disp)
        value.sinh()
        ipvec = self.TOpvec(ip)
        psqrdval = np.sum((np.sin(np.array(ipvec)*self.qunit)/Disp)**2)
        value = (psqrdval + value**2)
        value = value**0.5
        value.arcsinh()
        value = value* Disp
        if hasattr(value,'Stats'): value.Stats()
        return value



    def TOEBoot(self,ip,mass):
        value = (self.Getpsqrd(ip,actual=True) + mass**2)
        value = value**0.5
        if hasattr(value,'Stats'): value.Stats()
        return value

    ## expecting any momenta
    ## MassBoot BS or number
    ## output is in lattice units
    def EnergyFromRecRel(self,ip,MassBoot,DispIn=0.):
        if DispIn == 0.:
            return self.TOEBoot(ip,MassBoot)
        else:
            return self.TOLatEBoot(ip,MassBoot,Disp=DispIn)


    def Fmt_pmom_index(self,this_series,pmom_index,pre='p'):
        is_series = isinstance(this_series,pa.Series)
        if is_series:
            this_series = this_series.to_frame('data')
        pmom_list = this_series.index.names
        pmom_header = pmom_list[pmom_index]
        this_series.reset_index(inplace=True)
        this_series[pmom_header] = this_series[pmom_header].apply(lambda x : self.Getpsqrdform(x,pre=pre)).values
        this_series.set_index(pmom_list,inplace=True)
        if is_series:
            return this_series['data']
        else:
            return this_series


def TestMom(this_info = defInfo):
    data = LatticeParameters(mom2list=[1,2,3,4],Info=this_info)
    data.LoadPickle()
    return data

def ConvertQsqrd():
    from QuantityLists import ens_dict
    from GammaMatricies import CorrSpinTrace
    from Params import defInfoFlow
    mass_dict = {}
    mass_dict['mpi411'] = 0.5604128303
    mass_dict['mpi570'] = 0.6506026910
    mass_dict['mpi701'] = 0.7286377530
    mass_dict['L16'] = 0.9134165229
    mass_dict['L20'] = 0.7345331944
    mass_dict['L28'] = 0.5170724679

    a_dict = {}
    a_dict['mpi411'] =  0.0907
    a_dict['mpi570'] =  0.0907
    a_dict['mpi701'] =  0.0907
    a_dict['L16'] =  0.1095
    a_dict['L20'] =  0.0936
    a_dict['L28'] =  0.0684


    q_list = [[0,0,0],[1,0,0],[1,1,0],[1,1,1],[2,0,0]]
    q_keys = [0,1,2,3,4]
    old_df = pa.DataFrame()
    new_df = pa.DataFrame()
    calc_df = pa.DataFrame()
    phys_df = pa.DataFrame()
    for ikey,imass in mass_dict.items():
        defInfoFlow.update(ens_dict[ikey])
        defInfoFlow['a'] = a_dict[ikey]
        data = LatticeParameters(mom2list=[1,2,3,4],Info=defInfoFlow)
        # data.LoadPickle()
        this_oval,this_nval = [],[]
        this_calc,this_nphys = [],[]
        for iq in q_list:
            iqsqrd = np.sum(np.array(iq)**2)*data.qunit**2
            pmu,qmu,pmup = CorrSpinTrace().Create4Mom(np.array(iq)*data.qunit,
                                                         np.array([0,0,0]),imass)
            this_oval.append(np.abs(-np.sum(np.array(qmu)**2)))
            this_nval.append(np.abs(-np.sum(np.array(qmu[:4])**2)))
            this_calc.append(2*imass*(imass - np.sqrt(imass**2+iqsqrd)))
            this_nphys.append(this_nval[-1]*(data.hbarcdivlat**2))
        old_df[ikey] = pa.Series(this_oval,index=q_keys)
        new_df[ikey] = pa.Series(this_nval,index=q_keys)
        calc_df[ikey] = pa.Series(this_calc,index=q_keys)
        phys_df[ikey] = pa.Series(this_nphys,index=q_keys)
    return new_df,old_df,calc_df,phys_df


if __name__ == '__main__':
    data = TestMom()
    print('test parameters class is in variable data')
