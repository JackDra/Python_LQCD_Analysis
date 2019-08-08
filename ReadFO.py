#!/usr/bin/env python
import numpy as np
from XmlFormatting import tflowstr
from Params import nt
import os
import warnings

def FlowReadTopCharge(thisfile,tf_list='All'):
    if os.path.isfile(thisfile):
        readfile = thisfile
    elif os.path.isfile(thisfile.replace('ng00','ng')):
        readfile = thisfile.replace('ng00','ng')
    else:
        print('Warning, file not found',thisfile)
        return [],[]

    with open(readfile,'r') as ifile:
        tflowlist,topchargelist = [],[]
        for iline,line in enumerate(ifile):
            val = list(map(np.float64,line.split()))
            if len(val) == 3:
                tflow,topcharge,cmplxdump = val
            elif len(val) == 2:
                tflow,topcharge = val
            else:
                raise IOError('error with file: \n'+readfile+'\n on line: '+str(iline+1))
            if tf_list == 'All' or  tflowstr(tflow) in tf_list:
                tflowlist.append(np.float(tflow))
                topchargelist.append(topcharge)
    return tflowlist,topchargelist


def FlowReadTopCharge_New(thisfile,tf_list='All'):
    if os.path.isfile(thisfile):
        readfile = thisfile
    elif os.path.isfile(thisfile.replace('ng00','ng')):
        readfile = thisfile.replace('ng00','ng')
    else:
        print('Warning, file not found',thisfile)
        return [],[]
    with open(readfile,'r') as ifile:
        tflowlist,topchargelist = [],[]
        flow_line = True
        for iline,line in enumerate(ifile):
            if flow_line and len(line.strip()) != 0:

                val = list(map(np.float64,line.split()))
                if len(val) == 3:
                    tflow,topcharge,cmplxdump = val
                elif len(val) == 2:
                    tflow,topcharge = val
                else:
                    raise IOError('error with file: \n'+readfile+'\n on line: '+str(iline+1))
                ## this is fix for 10.0 problem
                if tf_list == 'All' or  tflowstr(tflow) in tf_list:
                    tflowlist.append(np.float(tflow))
                    topchargelist.append(topcharge)
                flow_line = False
            elif len(line.strip()) == 0:
                flow_line = True
    tflowlist,topchargelist = zip(*[[y,x] for y,x in sorted(zip(tflowlist,topchargelist))])
    return list(tflowlist),list(topchargelist)

def ReadTsInFlowOp(thisfile):
    if os.path.isfile(thisfile):
        with open(thisfile,'r') as ifile:
            tlist,topchargelist = [],[]
            tpos_list = []
            for iline,line in enumerate(ifile):
                hold = list(map(np.float64,line.split()))
                if len(hold) < 3:
                    if len(hold) == 2 or len(hold)== 1 or len(tpos_list) < 65:
                        warnings.warn('Error in file:'+str(thisfile))
                        print('Error line:',iline)
                    tpos_list = []
                    continue
                tval,topcharge,cmplxdump = hold
                tlist.append(tval)
                topchargelist.append(topcharge)
                tpos_list.append(topcharge)
        return tlist,topchargelist,False
    else:
        print('Warning, file not found',thisfile)
        return [],[],True

def FmtReadTs(thisfile,nt=nt):
    tlist,topchargelist,err = ReadTsInFlowOp(thisfile)
    tlist = np.array(list(tlist))
    topchargelist = np.array(list(topchargelist))
    out_tlist = tlist[1:nt+1]
    out_toplist,out_sum_toplist = [],[]
    out_tflowlist = []
    for i,(it,itop) in enumerate(zip(tlist,topchargelist)):
        if i%(nt+1) == 0:
            out_sum_toplist.append(itop)
            out_tflowlist.append(it)
        else:
            out_toplist.append(itop)
    out_toplist = np.array(out_toplist).reshape(-1,64)
    out_sum_toplist = np.array(out_sum_toplist)
    return out_tlist,out_tflowlist,out_sum_toplist,out_toplist,err

if __name__ == "__main__":
    import sys,os
    if len(sys.argv) == 1:
        print('please pass in a file to read')
    file_name = sys.argv[1]
    if not os.path.isfile(file_name):
        raise IOError('file not found: \n'+file_name)
    if len(sys.argv) > 2 and sys.argv[2].lower() == 'full':
        tlist,tflow_list,file_data_summed,file_data,err = FmtReadTs(file_name)
        print('read full flowed operator in tlist tflow_list summed_file_data file_data dimensions are [ t flow , time ] for file_data_summed and [ t flow ] for summed_file_data')
    else:
        try:
            tflow_list,file_data = FlowReadTopCharge_New(file_name)
        except Exception as err:
            tflow_list,file_data = FlowReadTopCharge(file_name)
        print('read flowed operator in file_data and tflow_list, dimensions are [ t flow ]')
    file_data = np.array(file_data)
