#!/usr/bin/env python
import numpy as np
from Params import nt
import MomParams as mp
import os
import sys
from GammaMatricies import GammaTOChroma
# from MomParams import *

SeekIncSize = 8
ChromaSIS = 16
bar3ptChromaByteSwap = True
# complextype = np.complex64
complextype = np.complex128
CheckMagic = False


def GetBar3ptLoc(ig,ip,tsinklen,momlistlen,thiscomplextype,PolProj=False):
    intlen = np.dtype(np.int32).itemsize
    if PolProj:
        strlen = np.dtype('S11').itemsize
    else:
        strlen = np.dtype('S13').itemsize
    cmplxlen = np.dtype(thiscomplextype).itemsize
    cmplxblocklen = cmplxlen*tsinklen+intlen
    # cmplxblocklen = intlen*16+intlen
    headerlen = intlen*27+strlen
    magicNumshift = intlen*6
    magicNumshiftCons = intlen*5
    momblocklen = cmplxblocklen + magicNumshift
    momblocklenCons = cmplxblocklen*2 + magicNumshiftCons
    momblockList = [momblocklen,momblocklenCons,momblocklenCons,momblocklen,momblocklenCons,
                    momblocklen,momblocklen,momblocklenCons,momblocklenCons,momblocklen,
                    momblocklen,momblocklenCons,momblocklen,momblocklenCons,momblocklenCons,
                    momblocklen,momblocklen,momblocklen,momblocklen,momblocklen]

    gammablocklen = momblocklen*momlistlen + intlen*2
    gammablocklenCons = momblocklenCons*momlistlen + intlen*2
    gammablockList = np.array([gammablocklen,gammablocklenCons,gammablocklenCons,gammablocklen,gammablocklenCons,
                      gammablocklen,gammablocklen,gammablocklenCons,gammablocklenCons,gammablocklen,
                      gammablocklen,gammablocklenCons,gammablocklen,gammablocklenCons,gammablocklenCons,
                      gammablocklen,gammablocklen,gammablocklen,gammablocklen,gammablocklen])
    outloc = headerlen + np.sum(gammablockList[:ig]) + momblockList[ig]*ip
    return outloc , outloc + cmplxblocklen,magicNumshift


## noDer[igamma , ip , it]
## forcent is total length, not array index!
class ReadFSCfunPickCHROMA:
    def __init__(self,  thisfileList,thisMomList,thisGammaList,forcent=nt,check_fs=False,
                        force_cmplx_type=complextype):
        params = mp.LatticeParameters()
        params.LoadPickle(no_write=True)
        self.data = []
        for thisfile in thisfileList:
            # if thistsink == 'FileName':
            #     thistsink = int(re.findall('tsink.*?p',thisfile)[0].replace('tsink','').replace('p',''))
            datahold = []
            if check_fs:
                if os.path.getsize(thisfile) == 826009: thiscomplextype = np.complex128
                elif os.path.getsize(thisfile) < 600000:
                    # raise IOError(thisfile+ ' Is old 64 bit, please remove all smaller sized files')
                    thiscomplextype = np.complex64
                else: thiscomplextype = complextype
            else: thiscomplextype = force_cmplx_type
            if not os.path.isfile(thisfile):
                print('ERROR FILE:',thisfile)
                raise IOError('Could not read 3 point correlation function')
            with open(thisfile,'rb') as f:
                try:
                    for igamma,thisgamma in enumerate(thisGammaList):
                        igammaloc = GammaTOChroma(thisgamma)
                        datahold.append([])
                        for ip,ipstr in enumerate(thisMomList):
                            iploc = params.TOip(ipstr)
                            loc,loccons,magmin = GetBar3ptLoc(igammaloc,iploc,forcent,len(params.momlist),thiscomplextype,PolProj=('GMA3' in thisfile))
                            magicloc = loc-magmin
                            if 'Cons' in thisgamma:
                                loc = loccons
                            if CheckMagic:
                                f.seek(magicloc)
                                MagicList =  np.fromfile(f,dtype=np.int32,count=1).byteswap()
                                if MagicList[0] != 20301:
                                    f.seek(magicloc)
                                    MagicList =  np.fromfile(f,dtype=np.int32,count=100).byteswap()
                                    print(MagicList)
                                    raise IOError('Magic Number not found for ' +thisgamma+' '+params.TOpstr(iploc))
                            f.seek(loc)
                            tmpdata = np.fromfile(f,dtype=thiscomplextype,count=forcent)
                            if bar3ptChromaByteSwap:tmpdata = tmpdata.byteswap()
                            if 'cmplx' in thisgamma:
                                datahold[igamma].append(np.float64(tmpdata.imag))
                            else:
                                datahold[igamma].append(np.float64(tmpdata.real))
                            # if isshift> nt-min(AllTSinkList)-1:
                            # datahold[igamma][-1][-corr_shift-1] = -datahold[igamma][-1][-corr_shift-1]
                            # print
                            # print tmpdata.real
                            # print
                            # totsign = np.sign(np.sum(np.sign(datahold[igamma][-1][:CMTSinkList[0]])))
                            # datahold[igamma][-1] = totsign*np.abs(datahold[igamma][-1])
                            datahold[igamma][-1] = np.sign(datahold[igamma][-1][0])*np.abs(datahold[igamma][-1])
                            # print np.abs(datahold[igamma][-1][13]-datahold[igamma][-1][12])/datahold[igamma][-1][12]
                            # if np.abs(datahold[igamma][-1][13]-datahold[igamma][-1][12])/datahold[igamma][-1][12] > 1.0:
                            #     print
                            #     print thisfile.replace(xsrcList[0],xsrc)
                            #     print corr_shift
                            #     for it,tdata in enumerate(datahold[igamma][-1]):
                            #         print it, tdata
                            # print np.roll(datahold[igamma][-1],isshift)
                            # print np.roll(datahold[igamma][-1],-isshift)
                            if any(np.isnan(datahold[igamma][ip])):
                                print()
                                print(thisfile)
                                print('thisgamma ',thisgamma,' iploc ', iploc)
                                print('data:')
                                for it,idata in enumerate(datahold[igamma][ip]):
                                    print(it,idata)
                                # if DeleteNanCfgs:
                                #     raise NaNCfunError('NaN Values: '+thisgamma+' ' +qvecSet[iploc]  )
                except Exception as err:
                    print(sys.exc_info())
                    print(err)
                    print('Error with: \n',thisfile)
            if len(self.data) == 0:
                self.data = np.array(datahold)
                datalen = 1
            else:
                self.data += np.array(datahold)
                datalen += 1
        self.data = self.data/float(datalen)


## noDer[igamma , ip , it]
## forcent is total length, not array index!
class RC3Full:
    def __init__(self,thisfileList,thisMomList,thisGammaList,forcent=nt,check_fs=False):
        params = mp.LatticeParameters()
        params.LoadPickle(no_write=True)
        self.data = []
        self.f_info = []
        for thisfile in thisfileList:
            # if thistsink == 'FileName':
            #     thistsink = int(re.findall('tsink.*?p',thisfile)[0].replace('tsink','').replace('p',''))
            datahold = []
            if check_fs:
                if os.path.getsize(thisfile) == 826009: thiscomplextype = np.complex128
                elif os.path.getsize(thisfile) < 600000:
                    # raise IOError(thisfile+ ' Is old 64 bit, please remove all smaller sized files')
                    thiscomplextype = np.complex64
                else: thiscomplextype = complextype
            else: thiscomplextype = complextype
            if not os.path.isfile(thisfile):
                print('ERROR FILE:',thisfile)
                raise IOError('Could not read 3 point correlation function')
            with open(thisfile,'rb') as f:
                f,this_f_info = ReadC3Header(f)
                self.f_info.append(this_f_info)
                for igamma,thisgamma in enumerate(thisGammaList):
                    igammaloc = GammaTOChroma(thisgamma)
                    datahold.append([])
                    for ip,ipstr in enumerate(thisMomList):
                        iploc = params.TOip(ipstr)
                        loc,loccons,magmin = GetBar3ptLoc(  igammaloc,iploc,forcent,len(params.momlist),
                                                            thiscomplextype,PolProj=('GMA3' in thisfile))
                        magicloc = loc-magmin
                        if 'Cons' in thisgamma:
                            loc = loccons
                        if CheckMagic:
                            f.seek(magicloc)
                            MagicList =  np.fromfile(f,dtype=np.int32,count=1).byteswap()
                            if MagicList[0] != 20301:
                                f.seek(magicloc)
                                MagicList =  np.fromfile(f,dtype=np.int32,count=100).byteswap()
                                print(MagicList)
                                raise IOError('Magic Number not found for ' +thisgamma+' '+params.TOpstr(iploc))
                        f.seek(loc)
                        tmpdata = np.fromfile(f,dtype=thiscomplextype,count=forcent)
                        if bar3ptChromaByteSwap:tmpdata = tmpdata.byteswap()
                        if 'cmplx' in thisgamma:
                            datahold[igamma].append(np.float64(tmpdata.imag))
                        else:
                            datahold[igamma].append(np.float64(tmpdata.real))
                        # if isshift> nt-min(AllTSinkList)-1:
                        # datahold[igamma][-1][-corr_shift-1] = -datahold[igamma][-1][-corr_shift-1]
                        # print
                        # print tmpdata.real
                        # print
                        # totsign = np.sign(np.sum(np.sign(datahold[igamma][-1][:CMTSinkList[0]])))
                        # datahold[igamma][-1] = totsign*np.abs(datahold[igamma][-1])
                        datahold[igamma][-1] = np.sign(datahold[igamma][-1][0])*np.abs(datahold[igamma][-1])
                        # print np.abs(datahold[igamma][-1][13]-datahold[igamma][-1][12])/datahold[igamma][-1][12]
                        # if np.abs(datahold[igamma][-1][13]-datahold[igamma][-1][12])/datahold[igamma][-1][12] > 1.0:
                        #     print
                        #     print thisfile.replace(xsrcList[0],xsrc)
                        #     print corr_shift
                        #     for it,tdata in enumerate(datahold[igamma][-1]):
                        #         print it, tdata
                        # print np.roll(datahold[igamma][-1],isshift)
                        # print np.roll(datahold[igamma][-1],-isshift)
                        if any(np.isnan(datahold[igamma][ip])):
                            print()
                            print(thisfile)
                            print('thisgamma ',thisgamma,' iploc ', iploc)
                            print('data:')
                            for it,idata in enumerate(datahold[igamma][ip]):
                                print(it,idata)
                            # if DeleteNanCfgs:
                            #     raise NaNCfunError('NaN Values: '+thisgamma+' ' +qvecSet[iploc]  )

            self.data.append(datahold)
        self.data = np.array(self.data)
        self.data = np.rollaxis(self.data,0,self.data.ndim)


def ReadC3Header(f_object):
    d1 = np.fromfile(f_object,dtype=np.int32,count=10).byteswap()
    outdict = {}
    outdict['output_version'] = d1[0]
    outdict['mom2max'] = d1[1]
    outdict['j_decay'] = d1[2]
    outdict['nrow'] = d1[4:8]
    outdict['wilson_version'] = d1[8]
    outdict['Num_seqsrc'] = d1[9]

    outdict['nuc_type'] = np.fromfile(f_object,dtype='S11',count=1)[0]
    d2 = np.fromfile(f_object,dtype=np.int32,count=17).byteswap()
    outdict['tsource'] = d2[0]
    outdict['tsink'] = d2[1]
    outdict['sinkmom'] = d2[3:6]
    outdict['gammanumber'] = d2[6]
    outdict['FFoutput_version'] = d2[7]
    outdict['numFF'] = d2[8]
    outdict['gammaopp'] = d2[9]
    outdict['numMom'] = d2[10]
    outdict['MagicNumber'] = d2[11]
    outdict['qmom'] = d2[13:16]
    outdict['tcurrlen'] = d2[16]
    return f_object,outdict

## Just a reference function for the extra parameters in the 3 point function file.
def holder(thisfile,cmplx_dtype=np.complex64):
    with open(thisfile,'rb') as f:
        d1 = np.fromfile(f,dtype=np.int32,count=10).byteswap()
        output_version = d1[0]
        mom2max = d1[1]
        j_decay = d1[2]
        nrow = d1[4:8]
        wilson_version = d1[8]
        Num_seqsrc = d1[9]

        nuc_type = np.fromfile(f,dtype='S13',count=1)[0]
        d2 = np.fromfile(f,dtype=np.int32,count=17).byteswap()
        tsource = d2[0]
        tsink = d2[1]
        sinkmom = d2[3:6]
        gammanumber = d2[6]
        FFoutput_version = d2[7]
        numFF = d2[8]
        gammaopp = d2[9]
        numMom = d2[10]
        MagicNumber = d2[11]
        qmom = d2[13:16]
        tcurrlen = d2[16]
        values = np.fromfile(f,dtype=cmplx_dtype,count=32).byteswap()


    print('output_version',output_version)
    print('mom2max',mom2max)
    print('j_decay',j_decay)
    print('nrow',nrow)
    print('wilson_version',wilson_version)
    print('Num_seqsrc',Num_seqsrc)

    print('nuc_type',nuc_type)
    print('tsource',tsource)
    print('tsink ',tsink)
    print('sinkmom ',sinkmom)
    print('gammanumber ',gammanumber)
    print('FFoutput_version ',FFoutput_version)
    print('numFF ',numFF)
    print('gammaopp ',gammaopp)
    print('numMom ',numMom)
    print('MagicNumber',MagicNumber)
    print('qmom ',qmom)
    print('tcurrlen ',tcurrlen)
    print('values')
    for ic,ival in enumerate(values):
        print(ic,ival.real,ival.imag)






##TODO, make sure working with multibal moms in thisMomList
class R2CChromaXMLFileList:
    def __init__(self,thisfileList,thisMomList,InterpNumb='9',MesOrBar='Baryon',tfix=False):
        if MesOrBar == 'Baryon':
            InterpFlag = 'baryon_num'
            tfix_interp_ref = 9
        elif MesOrBar == 'Meson':
            InterpFlag = 'gamma_value'
            tfix_interp_ref = 15
        else: raise IOError('Pass Baryon or Meson into correlator read routine')
        self.data = []
        self.tshiftlist = []
        # params = mp.LatticeParameters()
        # params.LoadPickle(no_write=True)
        datalen = 0
        for thisfile in thisfileList:
            self.OutMomList = []
            TSRC_read = False
            datahold = []
            datatfixhold = []
            # print 'Reading ' ,thisfile
            if not os.path.isfile(thisfile):
                print('ERROR FILE:',thisfile)
                raise IOError('Could not read 2 point correlation function')
            with open(thisfile,'r') as f:
                BarPart,InterpPart,InterpParttfix,ReadMom = False,False,False,False
                for line in f:
                    strline = line.strip()
                    if 't_srce' in strline and not TSRC_read:
                        TSRC_read = True
                        thissrclist = strline.replace('<t_srce>','').replace('</t_srce>','').split()
                        # print thisfile.replace(xsrcList[0],xsrc)
                        # print int(thissrclist[-1])
                        self.tshiftlist.append(int(thissrclist[-1]))
                    elif strline == '<Shell_Shell_Wilson_'+MesOrBar+'s>':
                    # if strline == '<Shell_Shell_Wilson_Baryons>':
                        BarPart = True
                    elif InterpFlag in strline:
                        if int(strline.replace('<'+InterpFlag+'>','').replace('</'+InterpFlag+'>','')) == int(InterpNumb):
                            InterpPart = True
                        else:
                            InterpPart = False
                        if int(strline.replace('<'+InterpFlag+'>','').replace('</'+InterpFlag+'>','')) == int(tfix_interp_ref):
                            InterpParttfix = True
                        else:
                            InterpParttfix = False
                    elif BarPart and InterpPart:
                        if '<sink_mom>' in strline:
                            thismom = strline.replace('<sink_mom>','').replace('</sink_mom>','')
                            if thismom in thisMomList:
                                # if len(datahold) > 0 and tfix:
                                #     maxindex = datahold[-1].index(np.max(datahold[-1]))
                                #     datahold[-1] = np.roll(np.array(datahold[-1]),-maxindex).tolist()
                                datahold.append([])
                                self.OutMomList.append(thismom)
                                ReadMom = True
                            else:
                                ReadMom = False
                        elif '<re>' in strline and ReadMom:
                            datahold[-1].append(np.float64(strline.replace('<re>','').replace('</re>','')))
                            # if np.isnan(datahold[-1][-1]) and DeleteNanCfgs:
                            #     raise NaNCfunError('NaN Values: '+thisfile+'  ' +qvecSet[int(self.OutMomList[-1])]  )
                    if BarPart and InterpParttfix and tfix:
                        if '<sink_mom>' in strline:
                            thismom = strline.replace('<sink_mom>','').replace('</sink_mom>','')
                            if thismom in thisMomList:
                                # if len(datatfixhold) > 0:
                                #     maxindex = datatfixhold[-1].index(np.max(datatfixhold[-1]))
                                #     datatfixhold[-1] = maxindex
                                datatfixhold.append([])
                                ReadMom = True
                            else:
                                ReadMom = False
                        elif '<re>' in strline and ReadMom:
                            datatfixhold[-1].append(np.float64(strline.replace('<re>','').replace('</re>','')))
                            # if np.isnan(datatfixhold[-1][-1]) and DeleteNanCfgs:
                            #     raise NaNCfunError('NaN Values: '+thisfile+'  ' +qvecSet[int(self.OutMomList[-1])]  )
                    if strline == '</momenta>' and str(max(int(InterpPart),9)):
                        if len(datahold) > 0 and ((len(datatfixhold) > 0) or (not tfix )): break
            # print xsrc, ' data '
            # print np.array(datahold)
            # print self.data
            # print
            if tfix:
                for ic,ishift in enumerate(datatfixhold):
                    maxindex = ishift.index(np.max(ishift))
                    datahold[ic] = np.roll(datahold[ic],-maxindex).tolist()
            if len(self.data) == 0:
                self.data = np.rollaxis(np.array(datahold),0,1)
                datalen = 1
            else:
                rolldata = np.rollaxis(np.array(datahold),0,1)
                if self.data.shape != rolldata.shape: raise IOError('error with '+thisfile)
                self.data += rolldata
                datalen += 1
        self.data = self.data/np.float64(datalen)

        # indicies =  np.searchsorted(map(params.Getpsqrd,self.OutMomList),map(params.Getpsqrd,thisMomList))
        # indicies =  [thisMomList.index(imom) for imom in  self.OutMomList]
        indicies =  [self.OutMomList.index(imom) for imom in  thisMomList]


        # # if Debug:
        # #     print
        # #     print thisfile
        # #     print thisMomList
        # #     print self.OutMomList
        # #     print indicies
        # #     print self.data
        # #     print self.datatfix
        # print map(params.Getpsqrd,self.OutMomList),map(params.Getpsqrd,thisMomList)
        # print indicies
        # print thisMomList, self.OutMomList, np.array(self.OutMomList)[indicies]
        self.data = np.array(self.data)[indicies].tolist()
        # if tfix: self.datatfix = np.array(self.datatfix)[indicies].tolist()




##TODO, make sure working with multibal moms in thisMomList
class RC2Full:
    def __init__(self,thisfileList,thisMomList,InterpNumb='9',MesOrBar='Baryon',Dog5=False):
        if MesOrBar == 'Baryon': InterpFlag = 'baryon_num'
        elif MesOrBar == 'Meson': InterpFlag = 'gamma_value'
        else: raise IOError('Pass Baryon or Meson into correlator read routine')
        self.data = []
        self.datag5 = []
        self.tshiftlist = []
        # params = mp.LatticeParameters()
        # params.LoadPickle()
        for thisfile in thisfileList:
            self.OutMomList = []
            TSRC_read = False
            datahold = []
            datag5hold = []
            # print 'Reading ' ,thisfile
            if not os.path.isfile(thisfile):
                print('ERROR FILE:',thisfile)
                raise IOError('Could not read 2 point correlation function')
            with open(thisfile,'r') as f:
                BarPart,InterpPart,InterpPartg5,ReadMom = False,False,False,False
                for line in f:
                    strline = line.strip()
                    if 't_srce' in strline and not TSRC_read:
                        TSRC_read = True
                        thissrclist = strline.replace('<t_srce>','').replace('</t_srce>','').split()
                        # print thisfile.replace(xsrcList[0],xsrc)
                        # print int(thissrclist[-1])
                        self.tshiftlist.append(int(thissrclist[-1]))
                    elif strline == '<Shell_Shell_Wilson_'+MesOrBar+'s>':
                    # if strline == '<Shell_Shell_Wilson_Baryons>':
                        BarPart = True
                    elif InterpFlag in strline:
                        if int(strline.replace('<'+InterpFlag+'>','').replace('</'+InterpFlag+'>','')) == int(InterpNumb):
                            InterpPartg5 = False
                            InterpPart = True
                        # elif int(strline.replace('<'+InterpFlag+'>','').replace('</'+InterpFlag+'>','')) == int(INg5):
                        #     InterpPart = False
                        #     InterpPartg5 = True
                        else:
                            InterpPart = False
                            InterpPartg5 = False
                    elif BarPart and InterpPart:
                        if '<sink_mom>' in strline:
                            thismom = strline.replace('<sink_mom>','').replace('</sink_mom>','')
                            if thismom in thisMomList:
                                datahold.append([])
                                self.OutMomList.append(thismom)
                                ReadMom = True
                            else:
                                ReadMom = False
                        elif '<re>' in strline and ReadMom:
                            datahold[-1].append(np.float64(strline.replace('<re>','').replace('</re>','')))
                            # if np.isnan(datahold[-1][-1]) and DeleteNanCfgs:
                            #     raise NaNCfunError('NaN Values: '+thisfile+'  ' +qvecSet[int(self.OutMomList[-1])]  )
                    elif BarPart and InterpPartg5:
                        if '<sink_mom>' in strline:
                            thismom = strline.replace('<sink_mom>','').replace('</sink_mom>','')
                            if thismom in thisMomList:
                                datag5hold.append([])
                                ReadMom = True
                            else:
                                ReadMom = False
                        elif '<re>' in strline and ReadMom:
                            datag5hold[-1].append(np.float64(strline.replace('<re>','').replace('</re>','')))
                            # if np.isnan(datag5hold[-1][-1]) and DeleteNanCfgs:
                            #     raise NaNCfunError('NaN Values: '+thisfile+'  ' +qvecSet[int(self.OutMomList[-1])]  )
                    if strline == '</momenta>' and InterpPart:
                        if len(datahold) > 0 and ((len(datag5hold) > 0) or (not Dog5 )): break
            # print xsrc, ' data '
            # print np.array(datahold)
            # print self.data
            # print
            self.data.append(np.rollaxis(np.array(datahold),0,1))
            if Dog5:
                self.datag5.append(np.rollaxis(np.array(datag5hold),0,1))

        # indicies =  np.searchsorted(map(params.Getpsqrd,self.OutMomList),map(params.Getpsqrd,thisMomList))
        # indicies =  [thisMomList.index(imom) for imom in  self.OutMomList]
        indicies =  [self.OutMomList.index(imom) for imom in  thisMomList]


        # # if Debug:
        # #     print
        # #     print thisfile
        # #     print thisMomList
        # #     print self.OutMomList
        # #     print indicies
        # #     print self.data
        # #     print self.datag5
        # print map(params.Getpsqrd,self.OutMomList),map(params.Getpsqrd,thisMomList)
        # print indicies
        # print thisMomList, self.OutMomList, np.array(self.OutMomList)[indicies]
        self.data = np.rollaxis(np.array(self.data),0,2)[indicies]
        # self.data = np.rollaxis(self.data,1,2).tolist()
        self.data = np.rollaxis(self.data,1,3).tolist()
        if Dog5:
            self.datag5 = np.rollaxis(np.array(self.datag5),0,2)[indicies]
            self.datag5 = np.rollaxis(self.datag5,1,3).tolist()






class NaNCfunError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


if __name__ == "__main__":
    file_name = sys.argv[1]
    if not os.path.isfile(file_name):
        raise IOError('file not found \n'+file_name)
    with open(file_name,'rb') as f:
        f,data = ReadC3Header(f)
    for ikey,ival in data.items():
        print(ikey, ival)






# ## Der[igamma , ip , it]
# class ReadFSDerCfunPick:
#     def __init__(self,file,thisMomList,thisDGList):
#         self.data = []
#         f = open(file,'rb')
#         tmpdata = array('d')
#         tmpdata.read(f,DGBSize*nmom*nt*2)
#         tmpdata.byteswap()
#         for igd,thisgamma in enumerate(thisDGList):
#             cmplxflag = False
#             if 'cmplx' in thisgamma:
#                 cmplxflag = True
#                 thisgamma = thisgamma.replace('cmplx','')
#             igdloc = DGSet.index(thisgamma)
#             self.data.append([])
#             for ip,iploc in enumerate(thisMomList):
#                 self.data[igd].append([])
#                 for it in xrange(nt):
#                     loc = igdloc*2 + iploc*DGBSize*2 + it*DGBSize*nmom*2
#                     if cmplxflag: loc += 1
#                     self.data[igd][ip].append(tmpdata[loc])
#                     if np.isnan(self.data[igd][ip][it]):
#                         raise NaNCfunError('NaN Values: '+thisgamma+' ' +qvecSet[iploc] + ' it='+str(it+1) )

#         f.close()

# ## noDer[igamma , ip , it]
# class ReadFSCfunPick:
#     def __init__(self,file,thisMomList,thisGammaList):
#         self.data = []
#         f = open(file,'rb')
#         tmpdata = array('d')
#         try:
#             tmpdata.read(f,GBSize*nmom*nt*2)
#         except:
#             print file
#         tmpdata.byteswap()
#         for igamma,thisgamma in enumerate(thisGammaList):
#             cmplxflag = False
#             if 'cmplx' in thisgamma:
#                 cmplxflag = True
#                 thisgamma = thisgamma.replace('cmplx','')
#             igammaloc = GammaSet.index(thisgamma)
#             self.data.append([])
#             for ip,iploc in enumerate(thisMomList):
#                 self.data[igamma].append([])
#                 for it in xrange(nt):
#                     loc = igammaloc*2 + iploc*GBSize*2 + it*GBSize*nmom*2
#                     if cmplxflag: loc += 1
#                     self.data[igamma][ip].append(tmpdata[loc])
#                     if np.isnan(self.data[igamma][ip][it]):
#                         raise NaNCfunError('NaN Values: '+thisgamma+ ' ' +qvecSet[iploc] + ' it='+str(it+1) )

#         f.close()

# ##[ip,it]
# class Read2ptCfunPickOld:
#     def __init__(self,file,thisMomList):
#         self.data = []
#         f = open(file,'rb')
#         tmpdata = array('d')
#         tmpdata.read(f,nmom*ns*ns*nt*2)
#         tmpdata.byteswap()
#         for ip,iploc in enumerate(thisMomList):
#             self.data.append([])
#             for it in xrange(nt):
#         ##11 + 22 projected spin component ns = 4 works
#                 loc = it*2*ns**2 + iploc*nt*2*ns**2
#                 self.data[ip].append(tmpdata[loc]+tmpdata[loc+10])
#                 if np.isnan(self.data[ip][it]):
#                     raise NaNCfunError('NaN Values: ' +qvecSet[iploc] + ' it='+str(it+1) )
#                 # if 'cmplx' in RI:
#                 #     self.datai[ipread].append(tmpdata[loc+1]+tmpdata[loc+11])
#             if self.data[ip][tsource-1] < 0.0:
#                 self.data[ip] = np.negative(self.data[ip])
#         f.close()
