#!/usr/bin/env python

import SetsOfFlowOps as sfo
import SetsOfTwoPtCorrs as so
import SetsOfRatios as sr
import SetsOfFFs as ffs
import os
from FileIO import ReadPickleWrap
from TimeStuff import Timer

# from RunSetCommands import FOpRun,FOpFullRun
# from RunSetCommands import TPRun,ARun,AFullRun
# from RunSetCommands import RRun,FFRun

def GetJobFun(file_name):
    if '2pt' in file_name or 'TwoPt' in file_name:
        return TPRun
    elif 'Alpha' in file_name:
        if 'Full' in file_name:
            return AFullRun
        else:
            return ARun
    elif 'FOp' in file_name:
        if 'Full' in file_name:
            return FOpFullRun
        else:
            return FOpRun
    elif 'FFRun' in file_name:
        return FFRun
    elif 'FFSolve' in file_name:
        return FFSolve
    elif 'FFSetup' in file_name:
        return FFSetup
    elif 'FFFinish' in file_name:
        return FFFinish
    elif 'Rat' in file_name:
        return RRun
    elif '3pt' in file_name or 'ThreePt':
        return Out3ptRun
    else:
        raise IOError('file name is not associated with a run type \n'+file_name)


def RunJob(PickleFile):
    thisFun = GetJobFun(PickleFile)
    thisInput = ReadPickleWrap(PickleFile)
    thistimer = Timer(name='Job file ' + PickleFile)
    data = thisFun(*thisInput)
    thistimer.Stop()
    return data


def FOpRun(InfoFlowList,Ahmed_comp,corr_sets,check_cfgs,def_wipe):
    data = sfo.SetOfFO(InfoDict=InfoFlowList,CompAhmed=Ahmed_comp,corr_cfglists=corr_sets)
    data.LoadPickle(DefWipe=def_wipe,CheckCfgs=check_cfgs)
    data.FlowMakeChit(DefWipe=def_wipe)
    return data

def FOpFullRun(InfoFlowList,corr_sets,check_cfgs,def_wipe,Only_Q):
    data = sfo.SetOfFOFull(InfoDict=InfoFlowList,corr_cfglists=corr_sets)
    data.LoadPickle(DefWipe=def_wipe,CheckCfgs=check_cfgs,OnlyQ=Only_Q)
    return data

def TPRun(InfoList,corr_sets,check_cfgs,def_wipe):
    data = so.SetOfTwoPt(InfoDict=InfoList,corr_cfglists=corr_sets)
    data.LoadPickle(DefWipe=def_wipe,CheckCfgs=check_cfgs)
    print('do powers and do coeff turned off due to fitting issue')
    # data.DoPowers()
    # data.DoCoeff()
    return data

def ARun(InfoFlowList,a_noise_sub,corr_sets,check_cfgs,def_wipe):
    data = so.SetOfNNQ(InfoDict=InfoFlowList,corr_cfglists=corr_sets)
    data.LoadPickle(DefWipe=def_wipe,CheckCfgs=check_cfgs,CheckMom=True)
    if a_noise_sub:
        data.AppendNoiseSub([0,1])
    return data

def AFullRun(InfoFlowList,corr_sets,check_cfgs,def_wipe,save_mem):
    data = so.SetOfNNQFull(InfoDict=InfoFlowList,corr_cfglists=corr_sets)
    data.LoadPickle(DefWipe=def_wipe,CheckCfgs=check_cfgs,CheckMom=True,SaveMem=save_mem)
    return data

def RRun(InfoFlowList,corr_sets,check_cfgs,def_wipe,Mcore):
    data = sr.SetOfRat(InfoDict=InfoFlowList,corr_cfglists=corr_sets)
    data.LoadPickle(DefWipe=def_wipe,CheckCfgs=check_cfgs,Mcore=Mcore)
    return data

def Out3ptRun(InfoFlowList,corr_sets,this_dir):
    data = sr.SetOfRat(InfoDict=InfoFlowList,corr_cfglists=corr_sets)
    data.Reformat3pt(this_dir)
    return data


def Out2ptRun(InfoFlowList,corr_sets,this_dir):
    data = so.SetOfTwoPt(InfoDict=InfoFlowList,corr_cfglists=corr_sets)
    data.Reformat2pt(this_dir)
    return data

def OutBoot2ptRun(InfoFlowList,corr_sets,this_dir):
    data = so.SetOfTwoPt(InfoDict=InfoFlowList,corr_cfglists=corr_sets)
    data.WriteBoot2pt(this_dir)
    return data

def FFSetup(InfoFlowList,def_wipe,hard_wipe,cfg_list_ffile,Mcore,CommonRat):
    if cfg_list_ffile:
        data = ffs.SetOfFF(cfglist='fromfile',cfglist2pt='fromfile',InfoDict=InfoFlowList)
    else:
        data = ffs.SetOfFF(InfoDict=InfoFlowList)
    data.SetupPickle(DefWipe=def_wipe,HardWipe=hard_wipe,Mcore=Mcore,CommonRatios=CommonRat)
    return data

def FFSolve(InfoFlowList,def_wipe,hard_wipe,cfg_list_ffile,Mcore,CommonRat):
    if cfg_list_ffile:
        data = ffs.SetOfFF(cfglist='fromfile',cfglist2pt='fromfile',InfoDict=InfoFlowList)
    else:
        data = ffs.SetOfFF(InfoDict=InfoFlowList)
    data.SolvePickle(DefWipe=def_wipe,HardWipe=hard_wipe,Mcore=Mcore,CommonRatios=CommonRat)
    return data

def FFFinish(   InfoFlowList,def_wipe,hard_wipe,cfg_list_ffile,
                Mcore,WipeFit,OnlyFits,DoFits,CommonRat):
    if cfg_list_ffile:
        data = ffs.SetOfFF(cfglist='fromfile',cfglist2pt='fromfile',InfoDict=InfoFlowList)
    else:
        data = ffs.SetOfFF(InfoDict=InfoFlowList)
    data.LoadPickle(DefWipe=def_wipe,HardWipe=hard_wipe,Mcore=Mcore,
                    WipeFit=WipeFit,OnlyFits=OnlyFits,DoFits=DoFits,CommonRatios=CommonRat)
    return data

def FFRun(InfoFlowList,def_wipe,hard_wipe,cfg_list_ffile,Mcore,CommonRat):
    if cfg_list_ffile:
        data = ffs.SetOfFF(cfglist='fromfile',cfglist2pt='fromfile',InfoDict=InfoFlowList)
    else:
        data = ffs.SetOfFF(InfoDict=InfoFlowList)
    print('Setting up Form Factors')
    data.SetupPickle(DefWipe=def_wipe,HardWipe=hard_wipe,Mcore=Mcore,CommonRatios=CommonRat)
    print('Solving Form Factors')
    data.SolvePickle(DefWipe=def_wipe,HardWipe=hard_wipe,Mcore=Mcore,CommonRatios=CommonRat)
    print('Finishing up Form Factors')
    data.LoadPickle(DefWipe=def_wipe,HardWipe=hard_wipe,Mcore=Mcore,CommonRatios=CommonRat)
    return data


def WriteCSH(job_type,pickle_file,job_params,n_job,Next_Job=None):
    Jstring = job_type + '_'+str(n_job)
    outlist = GetCSHHeader(Jstring,job_params)
    outlist.append('')
    outlist.append(r'echo "Begining run of '+pickle_file+'"')
    outlist.append(r'echo "started "`date`')
    # outlist.append(r'mpirun -np 1 '+job_params['python_exe']+' '+os.getcwd()+'/AnalysisGUI.py '+pickle_file.replace('./',os.getcwd()+'/'))
    outlist.append(job_params['python_exe']+' '+os.getcwd()+'/AnalysisGUI.py '+pickle_file.replace('./',os.getcwd()+'/'))
    outlist.append(r'if ($? != 0) then')
    outlist.append(r'   echo "Job Crashed, error code $? , exiting"')
    outlist.append(r'   exit 1')
    outlist.append(r'endif')
    outlist.append(r'echo ""')
    outlist.append(r'echo "finished "`date`')
    if Next_Job is not None:

        outlist.append(r'echo "Submitting next job '+Next_Job+' "')
        outlist.append(r'if(`where qsub` == "") then')
        outlist.append(r'   '+Next_Job)
        outlist.append(r'else')
        outlist.append(r'   qsub '+Next_Job)
        outlist.append(r'endif')
    CreateCSHFile(pickle_file.replace('.py3p','.csh'),outlist)


def CreateCSHFile(thisfile,outputlist):
    with open(thisfile,'w') as f:
        for iout in outputlist:
            f.write("%s\n" % iout)
        f.write('\n')
    os.system("chmod u+x "+thisfile)

def GetCSHHeader(Jstring,job_params):
    outlist = []
    outlist.append(r'#! /bin/tcsh')
    outlist.append('')
    if job_params['machine'] == 'juqueen':
        outlist.append(r'# @ job_name = '+Jstring)
        outlist.append(r'# @ error = $(job_name).$(jobid).out')
        outlist.append(r'# @ output = $(job_name).$(jobid).out')
        outlist.append(r'# @ environment = COPY_ALL')
        outlist.append(r'# @ wall_clock_limit = '+job_params['time'])
        outlist.append(r'# @ notification = error')
        outlist.append(r'# @ notify_user = '+job_params['email'])
        outlist.append(r'# @ job_type = bluegene')
        outlist.append(r'# @ bg_size = '+str(job_params['nproc']))
        outlist.append(r'# @ queue')
    else:
        if job_params['Scom'] == 'sbatch':
            outlist.append(r'#SBATCH --nodes 1')
            outlist.append(r'#SBATCH --cpus-per-task=1')
            outlist.append(r'#SBATCH --time='+job_params['time'])
            outlist.append(r'#SBATCH --mem='+job_params['memory'])
            outlist.append(r'#SBATCH -J '+Jstring)
            if 'RunPTG' in job_params and job_params['RunPTG']:
                outlist.append(r'#SBATCH -A ptg')
            # outlist.append(r'#SBATCH -p '+job_params['quetype'])
            # outlist.append(r'#SBATCH -n '+str(job_params['nproc']))
            # outlist.append(r'#SBATCH --time='+job_params['time'])
            # outlist.append(r'#SBATCH --mem='+job_params['memory'])
        elif job_params['Scom'] == 'qsub':
            if 'email' in job_params:
                outlist.append(r'#PBS -m bea')
                outlist.append(r'#PBS -M '+job_params['email'])
            outlist.append(r'#PBS -l walltime='+job_params['time']+',nodes=1:ppn='+str(job_params['nproc'])+',mem='+job_params['memory'])
            outlist.append(r'#PBS -N '+Jstring)
            if job_params['machine'] == 'laconia':
                if 'RunPTG' in job_params and job_params['RunPTG']:
                    outlist.append(r'#PBS -A ptg')
                else:
                    outlist.append(r'#PBS -A hpccdefault')
    outlist.append('')
    if 'modules' in job_params:
        for imod in job_params['modules'].values():
            outlist.append(r'module load '+imod)
    outlist.append('')
    return outlist
