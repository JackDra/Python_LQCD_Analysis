#!/usr/bin/env python

from FlowOpps import FlowOp
from Params import defInfoFlow,defxlimOpMost
from Params import outputdir,graphdir,scratchdir
from XmlFormatting import tflowTOLeg,untflowstr
from TimeStuff import Timer
import pandas as pa
from NullPlotData import null_series
import numpy as np
import os
from Autocorr import AutoCorrelate
from PredefFitFuns import LinearFitFun,LinFFDer,c3FitFun_nosqrt,c3FitFun_nosqrt_log
import SetsOfFits as sff
from copy import deepcopy
from PlotData import Plotting
from QuantityLists import ens_dict
import pickle as pik


# nx,nt,latspace = 16,32,0.1215
# kud,ks = 13825,13710
# beta = 1.83
# csw = 1761

# t0='t_f0.15'
# t0='t_f1.37'
# t0='t_f3.07'
# t0='t_f5.07'
force_t0 = False ## try running all lattices at to
mpi_latspace = 0.0907 # fm
ForceWipe = False
WipeBase = False
# t0_phys=float(untflowstr(t0))*mpi_latspace**2 ## in fm^2
# t0_phys=0.1 ## in fm^2
# t0_phys=0.3 ## in fm^2
# sqrt8t0_phys = 0.1 ## in fm
# sqrt8t0_phys = 0.3 ## in fm
my_eps_list = [f't_f{ival:.3}' for ival in np.arange(0,10.00,0.1)]
DoWW = False
DoAuto = False


def DoWWAnalysis(ens_name,data_plot,data_plot_NoE,data_Q,data_E,data_W,data_QW,InfoFlow,t0,this_a):
    ## this functions will have to do the bootstrapping, since it need to interpolate to t0_phys
    t0_phys = t0*this_a**2 ## in fm^2
    t0_phys_str = f't_fp{t0_phys:.3}'
    out_name_NoE = 'Test_WWrat_NoE_'+ens_name+'_'+t0_phys_str
    out_name = 'Test_WWrat_'+ens_name+'_'+t0_phys_str
    # t0_mpi = 'mpi' in ens_name or force_t0
    ens_latspace_2 = data_Q.latparams.latspace**2
    # data_QW_tt = FlowOp(Info=InfoFlow,man_load_cfgs=True)
    data_WW = FlowOp(Info=InfoFlow,man_load_cfgs=True)
    # data_QW_tt.FlowSetCustomName(string=ens_name+'_Test_QW_tt_'+t0_phys_str)
    data_WW.FlowSetCustomName(string=ens_name+'_Test_WW_'+t0_phys_str)
    WW_ratio = FlowOp(Info=InfoFlow,man_load_cfgs=True)
    WW_ratio.FlowSetCustomName(string=out_name)
    WW_NoE_ratio = FlowOp(Info=InfoFlow,man_load_cfgs=True)
    WW_NoE_ratio.FlowSetCustomName(string=out_name_NoE)
    if os.path.isfile(WW_ratio.flowPickleFile) and os.path.isfile(WW_NoE_ratio.flowPickleFile) and not ForceWipe:
        print('Found pickled file for Test_WWrat:')
        print(WW_ratio.flowPickleFile)
        print('Found pickled file for Test_WWrat_NoE:')
        print(WW_NoE_ratio.flowPickleFile)
        WW_ratio.FlowLoadPickle(CheckCfgs=False,DefWipe=False)
        WW_NoE_ratio.FlowLoadPickle(CheckCfgs=False,DefWipe=False)
        # data_QW.FlowSetCustomName(string=ens_name+'_Test_QW_'+t0_phys_str)
        data_WW.FlowSetCustomName(string=ens_name+'_Test_WW_'+t0_phys_str)
        # if os.path.isfile(data_QW.flowPickleFile) and not ForceWipe:
        #     data_QW.FlowLoadPickle(CheckCfgs=False,WipeData=False,DefWipe=False)
        # if not t0_mpi:
        #     if os.path.isfile(data_QW_tt.flowPickleFile) and not ForceWipe:
        #         data_QW_tt.FlowLoadPickle(CheckCfgs=False,WipeData=False,DefWipe=False)
        if os.path.isfile(data_WW.flowPickleFile) and not ForceWipe:
            data_WW.FlowLoadPickle(CheckCfgs=False,WipeData=False,DefWipe=False)
    else:
        # if os.path.isfile(data_QW_tt.flowPickleFile) and not ForceWipe:
        #     data_QW_tt.FlowLoadPickle(CheckCfgs=False,WipeData=False,DefWipe=False)
        # else:
        #     data_QW_tt = data_Q.FlowCombinePreBS(data_W,Operation='x')
        #     data_QW_tt.FlowSetCustomName(string=ens_name+'_QW_tt_'+t0_phys_str)
        #     data_QW_tt.FlowBootstrap(WipeData=False)
        #     data_QW_tt.FlowWrite()
        # if t0_mpi:
        #     data_QW.FlowSetCustomName(string=ens_name+'_Test_QW_'+t0_phys_str)
        #     data_QW.FlowWrite()
        #     if os.path.isfile(data_WW.flowPickleFile) and not ForceWipe:
        #         data_WW.FlowLoadPickle(CheckCfgs=False,WipeData=False,DefWipe=False)
        #     else:
        #         data_WW_low,data_WW_high,tlow,thigh = data_W.FlowCombinePreBS_FixFT_Interp(data_W,Operation='x',flow_time=t0_phys,Phys=True)
        #         data_WW_low.FlowSetCustomName(string=ens_name+'_Test_WW_low_'+str(tlow))
        #         data_WW_high.FlowSetCustomName(string=ens_name+'_Test_WW_high_'+str(thigh))
        #         data_WW_low.FlowBootstrap(WipeData=False)
        #         data_WW_high.FlowBootstrap(WipeData=False)
        #         data_WW_low.FlowWrite()
        #         data_WW_high.FlowWrite()
        #         fit_info = {}
        #         fit_info['Funs'] = (LinearFitFun,2)
        #         fit_info['name'] = data_WW_low.flowname + '_'+ data_WW_high.flowname + '_Interp'
        #         # tlow_index = list(data_WW_low.Op_Stats.index).index(tlow)
        #         # thigh_index = list(data_WW_high.Op_Stats.index).index(thigh)
        #         tm = t0_phys / (data_Q.latparams.latspace**2)
        #         t0_ratio_eval = (tm - untflowstr(tlow))/(untflowstr(thigh)-untflowstr(tlow))
        #         if 0>t0_ratio_eval >1:
        #             out_str = 'tm must be between the tlow, tight\n'
        #             out_str += 'tlow : '+untflowstr(tlow)+'\n'
        #             out_str += 'tmid : '+str(tm)+'\n'
        #             out_str += 'thigh : '+untflowstr(thigh)+'\n'
        #             raise EnvironmentError(out_str)
        #         def Fitt(value):
        #             data = pa.Series([value['low'],value['high']])
        #             this_fit = sff.SetOfFitFuns(data=data)
        #             this_test = this_fit.ScanRange(0,1,fit_info=fit_info,min_fit_len=0)
        #             if this_test:
        #                 this_fit.DoFits(show_timer=False)
        #                 this_fit.SortChi()
        #                 return this_fit.Fit_Stats.loc[:,'Fit'].iloc[0].Eval_Function(t0_ratio_eval)
        #             else:
        #                 return float('NaN')
        #         data_WW_df = data_WW_low.Op_Stats['boot'].to_frame('low')
        #         data_WW_df['high'] = data_WW_high.Op_Stats['boot']
        #         data_WW.FlowSetCustomName(string=ens_name+'_Test_WW_'+t0_phys_str)
        #         data_WW.Op_Stats['boot'] = data_WW_df.apply(Fitt,axis=1)
        #         data_WW.FlowWrite()
        #     QW_denom = data_QW.Op_Stats['boot'][t0]
        # else:
        flow_list = list(map(untflowstr,list(data_QW.Op_Stats.index)))
        def SortFun_pos(val):
            val = t0_phys-val*ens_latspace_2
            if val > 0:
                return abs(val)
            else:
                return 1000
        def SortFun_neg(val):
            val = t0_phys-val*ens_latspace_2
            if val < 0:
                return abs(val)
            else:
                return 1000
        left_flow = flow_list.index(min(flow_list,key=lambda x :SortFun_pos(x)))
        right_flow = flow_list.index(min(flow_list,key=lambda x :SortFun_neg(x)))
        tm = t0_phys / (data_Q.latparams.latspace**2)
        if left_flow > right_flow:
            out_str = 'problem deturmining closest flow times\n'
            out_str += 'left: '+str(left_flow) + '\n'
            out_str += 't0_phys: '+str(tm) + '\n'
            out_str += 'right: '+str(right_flow) + '\n'
            print(flow_list)
            raise EnvironmentError(out_str)
        flowfitr = 'fitr'+str(int(left_flow))+'-'+str(int(right_flow))
        t0_ratio_eval = (tm - untflowstr(flow_list[left_flow]))/(untflowstr(flow_list[right_flow])-untflowstr(flow_list[left_flow]))
        if 0>t0_ratio_eval >1:
            out_str = 'tm must be between the tlow, tight\n'
            out_str += 'tlow : '+str(flow_list[left_flow])+'\n'
            out_str += 'tmid : '+str(tm)+'\n'
            out_str += 'thigh : '+str(flow_list[right_flow])+'\n'
            raise EnvironmentError(out_str)
        data_QW.FlowSetFunction(LinearFitFun,2)
        data_QW.FlowFit(flowfit_range=flowfitr)
        QW_denom = data_QW.flowFit_Stats.Fit_Stats.loc[:,'Fit'].iloc[0].Eval_Function(t0_ratio_eval)

        if os.path.isfile(data_WW.flowPickleFile) and not ForceWipe:
            data_WW.FlowLoadPickle(CheckCfgs=False,WipeData=False,DefWipe=False)
        else:
            data_WW_low,data_WW_high,tlow,thigh = data_W.FlowCombinePreBS_FixFT_Interp(data_W,Operation='x',flow_time=t0_phys,Phys=True)
            data_WW_low.FlowSetCustomName(string=ens_name+'_Test_WW_low_'+str(tlow))
            data_WW_high.FlowSetCustomName(string=ens_name+'_Test_WW_high_'+str(thigh))
            data_WW_low.FlowBootstrap(WipeData=False)
            data_WW_high.FlowBootstrap(WipeData=False)
            data_WW_low.FlowWrite()
            data_WW_high.FlowWrite()
            fit_info = {}
            fit_info['Funs'] = (LinearFitFun,2)
            fit_info['name'] = data_WW_low.flowname + '_'+ data_WW_high.flowname + '_Interp'
            # tlow_index = list(data_WW_low.Op_Stats.index).index(tlow)
            # thigh_index = list(data_WW_high.Op_Stats.index).index(thigh)
            tm = t0_phys / (data_Q.latparams.latspace**2)
            t0_ratio_eval = (tm - untflowstr(tlow))/(untflowstr(thigh)-untflowstr(tlow))
            if 0>t0_ratio_eval >1:
                out_str = 'tm must be between the tlow, tight\n'
                out_str += 'tlow : '+untflowstr(tlow)+'\n'
                out_str += 'tmid : '+str(tm)+'\n'
                out_str += 'thigh : '+untflowstr(thigh)+'\n'
                raise EnvironmentError(out_str)
            def Fitt(value):
                data = pa.Series([value['low'],value['high']])
                this_fit = sff.SetOfFitFuns(data=data)
                this_test = this_fit.ScanRange(0,1,fit_info=fit_info,min_fit_len=0)
                if this_test:
                    this_fit.DoFits(show_timer=False)
                    this_fit.SortChi()
                    return this_fit.Fit_Stats.loc[:,'Fit'].iloc[0].Eval_Function(t0_ratio_eval)
                else:
                    return float('NaN')
            data_WW_df = data_WW_low.Op_Stats['boot'].to_frame('low')
            data_WW_df['high'] = data_WW_high.Op_Stats['boot']
            data_WW.FlowSetCustomName(string=ens_name+'_Test_WW_'+t0_phys_str)
            data_WW.Op_Stats['boot'] = data_WW_df.apply(Fitt,axis=1)
            data_WW.FlowWrite()
        print()
        print('QW_denom is',QW_denom.MakeValAndErr())
        print('evaulated at',t0_phys_str,' fm')
        tm = t0_phys / (data_Q.latparams.latspace**2)
        print('in lattice units','t_f'+f'{tm:.2}')
        print()
        WW_NoE_ratio = data_WW.__div__(QW_denom)
        WW_ratio = WW_NoE_ratio.__div__(data_E)
        WW_ratio.FlowSetCustomName(string=out_name)
        WW_ratio.FlowStats()
        WW_ratio.FlowWrite()
        WW_NoE_ratio.FlowSetCustomName(string=out_name_NoE)
        WW_NoE_ratio.FlowStats()
        WW_NoE_ratio.FlowWrite()
    data_plot = WW_ratio.FlowPlot_mul_tf2(data_plot)
    data_plot_NoE = WW_NoE_ratio.FlowPlot_mul_tf2(data_plot_NoE,mul=True)
    return WW_ratio,WW_NoE_ratio,data_plot,data_plot_NoE,data_WW


def GetCfgQuenched(ens_name,InfoFlow):
    nx,nt = int(InfoFlow['nxyzt'][0]),int(InfoFlow['nxyzt'][-1])
    beta_float = InfoFlow['Beta']
    beta = f'b{beta_float:.2f}'
    dim_str = f'{nx}x{nt}'
    cfglist_Q = pa.DataFrame()
    cfglist_W = pa.DataFrame()
    cfglist_E = pa.DataFrame()
    cval_Q,cval_W,cval_E,ival = [],[],[],[]
    master_folder = '/mnt/research/lqcd/weinberg_data/'
    master_folder += dim_str+beta + '/'
    if not os.path.isdir(master_folder):
        print('path not found:')
        print(master_folder)
        print()
        return None,None,None

    printed = False
    printedQNF,printedWNF,printedENF = False,False,False
    Q_folder = master_folder + 'Q/'
    W_folder = master_folder+ 'W/'
    E_folder = master_folder+ 'E/'
    for icfg in range(0,1500):
        file_name = f'conf{icfg:04}.dat'
        this_file_Q = Q_folder + file_name
        this_file_W = W_folder + file_name
        this_file_E = E_folder + file_name
        if not os.path.isfile(this_file_Q):
            if not printedQNF:
                print('FNF:')
                print(this_file_Q)
            printedQNF=True
            continue
        elif not os.path.isfile(this_file_W):
            if not printedWNF:
                print('FNF:')
                print(this_file_W)
            printedWNF=True
            continue
        elif not os.path.isfile(this_file_E):
            if not printedENF:
                print('FNF:')
                print(this_file_E)
            printedENF=True
            continue
        else:
            if not printed:
                print('adding files')
                print(this_file_Q)
                print(this_file_W)
                print(this_file_E)
                print()
            printed=True
            cval_Q.append(this_file_Q)
            cval_W.append(this_file_W)
            cval_E.append(this_file_E)
            ival.append(('-q-',str(icfg)))
    if len(cval_Q) == 0:
        print('No CFGS found')
        return None,None,None
    mi_ival = pa.MultiIndex.from_tuples(ival,names=['stream','number'])
    cfglist_Q['file_names'] = pa.Series(cval_Q,index=mi_ival)
    cfglist_W['file_names'] = pa.Series(cval_W,index=mi_ival)
    cfglist_E['file_names'] = pa.Series(cval_E,index=mi_ival)
    return cfglist_Q,cfglist_W,cfglist_E


def GetCfgEns(ens_name,InfoFlow):
    nx,nt = InfoFlow['nxyzt'][0],InfoFlow['nxyzt'][-1]
    kud,ks = InfoFlow['kud']//100,InfoFlow['ks']//100
    beta = str(InfoFlow['Beta'])
    beta = beta[0]+'.'+beta[1:]
    if beta[:3] == '1.9':
        beta = beta[:3]
    else:
        beta = beta[:4]

    csw = InfoFlow['Csw']

    beta_str = 'B'+str(InfoFlow['Beta'])
    RC_str = 'RC'+str(nx)+'x'+str(nt)
    RC_str2 = str(nx)+'X'+str(nt)
    kud_str = 'Kud0'+str(kud)+'00Ks0'+str(ks)+'00'
    if nx == 20:
        kud_str2 = 'b'+str(beta)+'L20kl0.'+str(kud)+'ks0.'+str(ks)[:-1]
        kud_str3 = 'Kud0'+str(kud)+'Ks0'+str(ks)
    elif nx == 32:
        kud_str2 = 'b'+str(beta)+'kl0.'+str(kud)+'ks0.'+str(ks)[:-1]
        kud_str3 = 'Kud0'+str(kud)+'00Ks0'+str(ks)+'00'
    else:
        kud_str3 = 'Kud0'+str(kud)+'Ks0'+str(ks)
        kud_str2 = 'b'+str(beta)+'kl0.'+str(kud)+'ks0.'+str(ks)


    Csw = 'C'+str(csw)
    cfglist_Q = pa.DataFrame()
    cfglist_W = pa.DataFrame()
    cfglist_E = pa.DataFrame()
    cval_Q,cval_W,cval_E,ival = [],[],[],[]
    for istream in ['-1-','-2-','-3-','-4-','-a-','-b-']:
    # for istream in ['-1-','-2-','-a-','-b-']:
        if not os.path.isdir('/mnt/research/lqcd/CfunAnalysis/FlowNewForm/'+RC_str+'_'+kud_str+istream):
            print('path not found:')
            print('/mnt/research/lqcd/CfunAnalysis/FlowNewForm/'+RC_str+'_'+kud_str+istream)
            print()
            continue
        print('Making file list for stream',istream)
        printed = False
        printedQNF,printedWNF,printedENF = False,False,False
        for icfg in range(0,700):
            icfg_pad = f'{icfg:03}'
            this_file_Q = '/mnt/research/lqcd/CfunAnalysis/FlowNewForm/'+RC_str+'_'+kud_str+istream+'/PerGF/q_flow_b1.90_ng@0.out'.replace('@',icfg_pad)

            this_file_W = '/mnt/research/lqcd/qcEDM/'+kud_str2+'_'+RC_str2+'/data/weinberg/'+RC_str+'_'+beta_str+kud_str3+Csw+istream+'00@0/weinberg_T_0.00.dat'.replace('@',icfg_pad)

            this_file_E = '/mnt/research/lqcd/qcEDM/'+kud_str2+'_'+RC_str2+'/data/energy/'+RC_str+'_'+beta_str+kud_str3+Csw+istream+'00@0/energy_T_0.00.dat'.replace('@',icfg_pad)
            if not os.path.isfile(this_file_Q):
                if not printedQNF:
                    print('FNF:')
                    print(this_file_Q)
                printedQNF=True
                continue
            elif not os.path.isfile(this_file_W):
                if not printedWNF:
                    print('FNF:')
                    print(this_file_W)
                printedWNF=True
                continue
            elif not os.path.isfile(this_file_E):
                if not printedENF:
                    print('FNF:')
                    print(this_file_E)
                printedENF=True
                continue
            else:
                if not printed:
                    print('adding files')
                    print(this_file_Q)
                    print(this_file_W)
                    print(this_file_E)
                    print()
                printed=True
                cval_Q.append(this_file_Q)
                cval_W.append(this_file_W)
                cval_E.append(this_file_E)
                ival.append((istream,str(icfg)))
    if len(cval_Q) == 0:
        print('No CFGS found')
        return None,None,None
    mi_ival = pa.MultiIndex.from_tuples(ival,names=['stream','number'])
    cfglist_Q['file_names'] = pa.Series(cval_Q,index=mi_ival)
    cfglist_W['file_names'] = pa.Series(cval_W,index=mi_ival)
    cfglist_E['file_names'] = pa.Series(cval_E,index=mi_ival)
    return cfglist_Q,cfglist_W,cfglist_E




def AlphaRatAnalysis(ens_name,InfoFlow,data_plot,data_plot_NoE,
                     data_plot_tauint,data_plot_covar,t0_scale,nblock=1):
    t0_scale_str = f't_s{t0_scale:.3}'
    leg_lab_ext = r' \frac{t_{s}}{t_{0}}='+t0_scale_str.replace('t_s','')
    InfoFlow['tflowlist'] = defxlimOpMost
    InfoFlow['n_block'] = nblock

    data_Q = FlowOp(Info=InfoFlow,man_load_cfgs=True)
    data_W = FlowOp(Info=InfoFlow,man_load_cfgs=True)
    data_E = FlowOp(Info=InfoFlow,man_load_cfgs=True)
    data_QW = FlowOp(Info=InfoFlow,man_load_cfgs=True)
    data_QQ = FlowOp(Info=InfoFlow,man_load_cfgs=True)
    data_QQ0 = FlowOp(Info=InfoFlow,man_load_cfgs=True)
    if nblock > 1:
        nbstr = '_nb'+str(nblock)
    else:
        nbstr = ''
    out_name_NoE = 'Test_alpha_NoE_'+ens_name+'_'+t0_scale_str+nbstr
    out_name_Q0= 'Test_alpha_Q0_'+ens_name+'_'+t0_scale_str+nbstr
    out_name = 'Test_alpha_'+ens_name+'_'+t0_scale_str+nbstr



    alpha_ratio = FlowOp(Info=InfoFlow,man_load_cfgs=True)
    alpha_ratio.FlowSetCustomName(string=out_name)
    alpha_ratio_NoE = FlowOp(Info=InfoFlow,man_load_cfgs=True)
    alpha_ratio_NoE.FlowSetCustomName(string=out_name_NoE)
    alpha_ratio_Q0 = FlowOp(Info=InfoFlow,man_load_cfgs=True)
    alpha_ratio_Q0.FlowSetCustomName(string=out_name_Q0)
    DoRead = True
    if  os.path.isfile(alpha_ratio.flowPickleFile) and \
        os.path.isfile(alpha_ratio_NoE.flowPickleFile) and \
        os.path.isfile(alpha_ratio_Q0.flowPickleFile) and not ForceWipe:
        print('Found pickled file for Test_alpha:')
        print(alpha_ratio.flowPickleFile)
        alpha_ratio.FlowLoadPickle(CheckCfgs=False,DefWipe=False)
        if 'Auto' in alpha_ratio.Op_Stats.columns:
            alpha_ratio.FlowStats()
            alpha_ratio.FlowWrite()
        print('Found pickled file for Test_alpha_NoE:')
        print(alpha_ratio_NoE.flowPickleFile)
        alpha_ratio_NoE.FlowLoadPickle(CheckCfgs=False,DefWipe=False)
        if 'Auto' in alpha_ratio_NoE.Op_Stats.columns:
            alpha_ratio_NoE.FlowStats()
            alpha_ratio_NoE.FlowWrite()

        print('Found pickled file for Test_alpha_Q0:')
        print(alpha_ratio_Q0.flowPickleFile)
        alpha_ratio_Q0.FlowLoadPickle(CheckCfgs=False,DefWipe=False)
        if 'Auto' in alpha_ratio_Q0.Op_Stats.columns:
            alpha_ratio_Q0.FlowStats()
            alpha_ratio_Q0.FlowWrite()

        data_Q.FlowSetCustomName(string=ens_name+'_Test_Q'+nbstr)
        data_E.FlowSetCustomName(string=ens_name+'_Test_E'+nbstr)
        data_W.FlowSetCustomName(string=ens_name+'_Test_W'+nbstr)
        data_QQ.FlowSetCustomName(string=ens_name+'_Test_QQ'+nbstr)
        data_QQ0.FlowSetCustomName(string=ens_name+'_Test_QQ0'+nbstr)
        data_QW.FlowSetCustomName(string=ens_name+'_Test_QW_'+t0_scale_str+nbstr)
        data_QW.FlowSetCustomName(string=ens_name+'_Test_QW_'+t0_scale_str+nbstr,stringLL=data_QW.flowLegLab+leg_lab_ext)
        if os.path.isfile(data_Q.flowPickleFile):
            data_Q.FlowLoadPickle(CheckCfgs=False,WipeData=False,DefWipe=WipeBase)
        if os.path.isfile(data_E.flowPickleFile):
            data_E.FlowLoadPickle(CheckCfgs=False,WipeData=False,DefWipe=WipeBase)
            data_E.Compute_t0()
            this_t0 = data_E.t0_values['boot'].iloc[0]
            ts = t0_scale*this_t0
            ts.Stats()
            ts = ts.Avg
        if os.path.isfile(data_W.flowPickleFile):
            data_W.FlowLoadPickle(CheckCfgs=False,WipeData=False,DefWipe=WipeBase)
        if os.path.isfile(data_QQ.flowPickleFile) and not ForceWipe:
            data_QQ.FlowLoadPickle(CheckCfgs=False,WipeData=False,DefWipe=False)
        if os.path.isfile(data_QQ0.flowPickleFile) and not ForceWipe:
            data_QQ0.FlowLoadPickle(CheckCfgs=False,WipeData=False,DefWipe=False)
        if os.path.isfile(data_QW.flowPickleFile) and not ForceWipe:
            data_QW.FlowLoadPickle(CheckCfgs=False,WipeData=False,DefWipe=False)
        DoRead = False
    if DoRead:
        if 'qu_' in ens_name:
            cfglist_Q,cfglist_W,cfglist_E = GetCfgQuenched(ens_name,InfoFlow)
        else:
            cfglist_Q,cfglist_W,cfglist_E = GetCfgEns(ens_name,InfoFlow)
        if cfglist_Q is None:
            raise IOError('failed to find configurations')
        cfg_list_fileQ = outputdir + ens_name+'_cfglist_Q.csv'
        cfg_list_fileW = outputdir + ens_name+'_cfglist_W.csv'
        cfg_list_fileE = outputdir + ens_name+'_cfglist_E.csv'
        print('outputting cfglist to files:')
        print(cfg_list_fileQ)
        print(cfg_list_fileW)
        print(cfg_list_fileE)

        cfglist_Q.to_csv(cfg_list_fileQ)
        cfglist_W.to_csv(cfg_list_fileW)
        cfglist_E.to_csv(cfg_list_fileE)

        data_Q.FlowImportCfgList(cfglist_Q)
        data_Q.FlowSetCustomName(string=ens_name+'_Test_Q'+nbstr)
        data_Q.FlowLoadPickle(CheckCfgs=False,WipeData=False,DefWipe=WipeBase)

        data_W.FlowImportCfgList(cfglist_W)
        data_W.FlowSetCustomName(string=ens_name+'_Test_W'+nbstr)
        data_W.FlowLoadPickle(CheckCfgs=False,WipeData=False,DefWipe=WipeBase)

        data_E.FlowImportCfgList(cfglist_E)
        data_E.FlowSetCustomName(string=ens_name+'_Test_E'+nbstr)
        data_E.FlowLoadPickle(CheckCfgs=False,WipeData=False,DefWipe=WipeBase)
        data_E.Compute_t0()
        this_t0 = data_E.t0_values['boot'].iloc[0]
        ts = t0_scale*this_t0
        ts.Stats()
        ts = ts.Avg
        data_QW = FlowOp(Info=InfoFlow,man_load_cfgs=True)
        data_QQ = FlowOp(Info=InfoFlow,man_load_cfgs=True)
        data_QQ0 = FlowOp(Info=InfoFlow,man_load_cfgs=True)
        data_QW.FlowSetCustomName(string=ens_name+'_Test_QW_'+t0_scale_str+nbstr)
        data_QW.FlowSetCustomName(string=ens_name+'_Test_QW_'+t0_scale_str+nbstr,stringLL=data_QW.flowLegLab+leg_lab_ext)
        data_QQ.FlowSetCustomName(string=ens_name+'_Test_QQ'+nbstr)
        data_QQ0.FlowSetCustomName(string=ens_name+'_Test_QQ0'+nbstr)

        if os.path.isfile(data_QW.flowPickleFile) and not ForceWipe:
            data_QW.FlowLoadPickle(CheckCfgs=False,WipeData=False,DefWipe=False)
        else:
            data_QW_low,data_QW_high,tlow,thigh = data_W.FlowCombinePreBS_FixFT_Interp(data_Q,Operation='x',flow_time=ts,Phys=False)
            data_QW_low.FlowSetCustomName(string=ens_name+'_Test_QW_low_'+str(tlow)+nbstr)
            data_QW_high.FlowSetCustomName(string=ens_name+'_Test_QW_high_'+str(thigh)+nbstr)
            data_QW_low.FlowBootstrap(WipeData=False)
            data_QW_high.FlowBootstrap(WipeData=False)
            data_QW_low.FlowWrite()
            data_QW_high.FlowWrite()
            fit_info = {}
            fit_info['Funs'] = (LinearFitFun,2)
            fit_info['name'] = data_QW_low.flowname + '_'+ data_QW_high.flowname + '_Interp'
            # tlow_index = list(data_QW_low.Op_Stats.index).index(tlow)
            # thigh_index = list(data_QW_high.Op_Stats.index).index(thigh)
            t0_ratio_eval = (ts - untflowstr(tlow))/(untflowstr(thigh)-untflowstr(tlow))
            if 0>t0_ratio_eval >1:
                out_str = 'tm must be between the tlow, tight\n'
                out_str += 'tlow : '+untflowstr(tlow)+'\n'
                out_str += 'tmid : '+str(ts)+'\n'
                out_str += 'thigh : '+untflowstr(thigh)+'\n'
                raise EnvironmentError(out_str)
            def Fitt(value):
                data = pa.Series([value['low'],value['high']])
                this_fit = sff.SetOfFitFuns(data=data)
                this_test = this_fit.ScanRange(0,1,fit_info=fit_info,min_fit_len=0)
                if this_test:
                    this_fit.DoFits(show_timer=False)
                    this_fit.SortChi()
                    return this_fit.Fit_Stats.loc[:,'Fit'].iloc[0].Eval_Function(t0_ratio_eval)
                else:
                    return float('NaN')
            data_QW_df = data_QW_low.Op_Stats['boot'].to_frame('low')
            data_QW_df['high'] = data_QW_high.Op_Stats['boot']
            data_QW.FlowSetCustomName(string=ens_name+'_Test_QW_'+t0_scale_str+nbstr)
            data_QW.FlowSetCustomName(string=ens_name+'_Test_QW_'+t0_scale_str+nbstr,
                                      stringLL=data_QW.flowLegLab+leg_lab_ext)
            data_QW.Op_Stats['boot'] = data_QW_df.apply(Fitt,axis=1)
            data_QW.FlowWrite()
        if os.path.isfile(data_QQ.flowPickleFile) and not ForceWipe:
            data_QQ.FlowLoadPickle(CheckCfgs=False,WipeData=False,DefWipe=False)
        else:
            data_QQ = data_Q.FlowCombinePreBS(data_Q,Operation='x')
            data_QQ.FlowSetCustomName(string=ens_name+'_Test_QQ'+nbstr)
            data_QQ.FlowBootstrap(WipeData=False)
            data_QQ.FlowWrite()

        if os.path.isfile(data_QQ0.flowPickleFile) and not ForceWipe:
            data_QQ0.FlowLoadPickle(CheckCfgs=False,WipeData=False,DefWipe=False)
        else:
            data_QQ0 = data_Q.FlowCombinePreBS_FixFT(data_Q,Operation='x',flow_time=0)
            data_QQ0.FlowSetCustomName(string=ens_name+'_Test_QQ0'+nbstr)
            data_QQ0.FlowBootstrap(WipeData=False)
            data_QQ0.FlowWrite()

        # alpha_ratio = data_QW.__div__(data_E)
        # alpha_ratio.Op_Stats['boot'] = alpha_ratio.Op_Stats['boot'] / data_QQ.Op_Stats['boot'][t0]
        alpha_ratio_Q0 = data_QW.__div__(data_E*data_QQ0.GetInterp(ts))
        alpha_ratio_NoE = data_QW.__div__(data_QQ.GetInterp(ts))
        alpha_ratio = alpha_ratio_NoE.__div__(data_E)
        # this_df = pa.DataFrame()
        # this_df['Op_QQ'] = data_QQ.Op_cfgs['Op']
        # this_df['Op_QW'] = data_QW.Op_cfgs['Op']
        # this_df['Op_E'] = data_E.Op_cfgs['Op']
        # box_len_3 = alpha_ratio.latparams.nxyz**3
        if DoAuto:
            raise NotImplementedError('Auto not implemented since we interpolating')
            # Alist = []
            # thistimer = Timer(linklist=alpha_ratio.tflowfloat,name='Autocorrelating '+ens_name)
            # for ictflow,strit in enumerate(alpha_ratio.Op_Stats.index):
            #     itf = float(untflowstr(strit))
            #     # bl3_div_itf = box_len_3/itf
            #     bl3_div_itf = itf
            #     def RatAutoFun(QW,E,QQ):
            #         return bl3_div_itf * QW/(E*QQ)
            #     def RatAutFunDer(QW,E,QQ):
            #         return [bl3_div_itf * np.reciprocal(E*QQ),-bl3_div_itf * QW/(QQ*np.power(E,2)),-bl3_div_itf * QW/(E*np.power(QQ,2))]
            #     tdata = pa.DataFrame()
            #     tdata['Op_QW'] = data_QW.Op_cfgs['Op'].apply(lambda x: x[ictflow])
            #     tdata['Op_E'] = data_E.Op_cfgs['Op'].apply(lambda x: x[ictflow])
            #     tdata['Op_QQ'] = data_QQ.Op_cfgs['Op'].apply(lambda x: x[list(alpha_ratio.tflowlist).index(ts)])
            #     Alist.append(AutoCorrelate( Fun=[RatAutoFun,RatAutFunDer],Sparam=alpha_ratio.Sparam,
            #                                 name=alpha_ratio.flowname + ' $' +tflowTOLeg(strit)+'$',data=tdata,save_covar=True))
            #     thistimer.Lap(strit)
            # alpha_ratio.Op_Stats['Auto'] = pa.Series(Alist,index=alpha_ratio.Op_Stats.index)
        # alpha_ratio = data_QW.__div__(data_QQ*data_E)
        alpha_ratio.FlowSetCustomName(string=out_name)
        alpha_ratio.FlowSetCustomName(string=out_name,stringLL=alpha_ratio.flowLegLab+leg_lab_ext)
        alpha_ratio.FlowStats()
        alpha_ratio.FlowWrite()
        alpha_ratio_NoE.FlowSetCustomName(string=out_name_NoE)
        alpha_ratio_NoE.FlowSetCustomName(string=out_name_NoE,stringLL=alpha_ratio_NoE.flowLegLab+leg_lab_ext)
        alpha_ratio_NoE.FlowStats()
        alpha_ratio_NoE.FlowWrite()
        alpha_ratio_Q0.FlowSetCustomName(string=out_name_Q0)
        alpha_ratio_Q0.FlowSetCustomName(string=out_name_Q0,stringLL=alpha_ratio_Q0.flowLegLab+leg_lab_ext)
        alpha_ratio_Q0.FlowStats()
        alpha_ratio_Q0.FlowWrite()
    # data_plot_NoE = alpha_ratio_NoE.FlowPlot_mul_tf2(data_plot_NoE,mul=True)
    data_plot_NoE = alpha_ratio_NoE.FlowPlot_mul_tf2(data_plot_NoE,mul=True,forcet0=t0_dict[ens_name])
    # data_plot_NoE = alpha_ratio_NoE.FlowPlot(data_plot_NoE)
    data_plot = alpha_ratio.FlowPlot_mul_tf2(data_plot)
    # if DoAuto:
    #     data_plot = alpha_ratio.FlowPlotAuto(data_plot)
    #     data_plot_tauint =  alpha_ratio.PlotTauIntWOpt(data_plot_tauint)
    #     data_plot_covar = PlotCovar(alpha_ratio.Op_Stats['Auto'],data_plot_covar,ens_name)
    return alpha_ratio,alpha_ratio_Q0,alpha_ratio_NoE,data_plot,data_plot_NoE,data_plot_tauint, \
            data_Q,data_E,data_W,data_QW,data_QQ,data_QQ0,data_plot_covar,InfoFlow,ts


def PlotCovar(Auto_Series,data_plot,label):
    this_key = [slice(None),'C11']
    covarl,indexl = [],[]
    for ikey,idata in Auto_Series.items():
        for ic,icovar in enumerate(idata.covar):
            for jc,jcovar in enumerate(icovar):
                if jc >= ic:
                    indexl.append((ikey,'C'+str(ic+1)+str(jc+1)))
                    covarl.append(jcovar)
        for ic,icovar in enumerate(idata.covar_NoNorm):
            for jc,jcovar in enumerate(icovar):
                if jc >= ic:
                    indexl.append((ikey,'C'+str(ic+1)+str(jc+1)+'_NoNorm'))
                    covarl.append(jcovar/np.sqrt(idata.Avg_paras[ic]*idata.Avg_paras[jc]))
    if len(indexl) > 0:
        indicies = pa.MultiIndex.from_tuples(indexl,names=['FlowTime','Covar_matrix'])
        ploty = pa.Series(covarl,index=indicies)
    else:
        ploty = pa.Series()
    hold_series = null_series
    hold_series['x_data'] = 'from_keys'
    hold_series['key_select'] = this_key
    hold_series['y_data'] = ploty
    hold_series['yerr_data'] = None
    hold_series['type'] = 'scatter_vary'
    hold_series['label'] = label
    data_plot.AppendData(hold_series)
    return data_plot


def PlotErrRat(this_data,this_data2):
    print('plotting ratio of relative error WW/WQ')
    this_info = pa.Series()
    this_info['save_file'] = graphdir+'/flowops_WWdivWQ_errrat.pdf'
    this_info['title'] = r'$\frac{\Delta<W(t)W(t_0)>}{\Delta<W(t)Q(t_0)>}$ '
    this_info['xlabel'] = '$t/a$'
    this_info['ylabel'] = 'ratio'
    data_plot = Plotting(plot_info=this_info)
    for (ikey,idata),idata2 in zip(this_data.items(),this_data2.values()):
        if len(idata.Op_Stats['boot'].index) == 0:
            print('Warning, ',ikey,' did not pass properly, skipping')
            continue
        this_data = idata.Op_Stats['boot'].apply(lambda x : np.abs(x.Std/x.Avg))
        this_data = this_data/idata2.Op_Stats['boot'].apply(lambda x : np.abs(x.Std/x.Avg))
        hold_series = null_series
        hold_series['x_data'] = list(map(untflowstr,this_data.index))
        hold_series['y_data'] = this_data.values
        hold_series['label'] = ikey
        hold_series['shift'] = 0.0
        hold_series['type'] = 'plot'
        data_plot.AppendData(hold_series)
    data_plot.PlotAll()
    return data_plot


def PlotEt2(this_data):
    print('plotting ratio of relative error WW/WQ')
    this_info = pa.Series()
    this_info['save_file'] = graphdir+'/flowops_Et2.pdf'
    this_info['title'] = r'$(\frac{6}{N_{xyzt}})t^{2}<E(t)>$'
    this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
    this_info['ylabel'] = r'$(\frac{6}{N_{xyzt}}t^{2}<E(t)>$'
    data_plot = Plotting(plot_info=this_info)

    for ikey,idata in this_data.items():
        if len(idata.Op_Stats['boot'].index) == 0:
            print('Warning, ',ikey,' did not pass properly, skipping')
            continue
        def this_fun(val):
            return (untflowstr(val.name)**2)*6*val['boot']/(idata.latparams.nt*idata.latparams.nxyz**3)
        this_data = idata.Op_Stats.apply(this_fun,axis=1)
        this_data = this_data.apply(lambda x : x.Stats())
        this_s8t0 = (8*idata.latparams.Get_t0()['boot'].iloc[0]).Sqrt()
        tflowphys,tflowphys_err = [],[]
        for itflow in map(untflowstr,list(idata.Op_Stats.index)):
            it = np.sqrt(8*itflow)/this_s8t0
            it.Stats()
            tflowphys.append(it.Avg)
            tflowphys_err.append(it.Std)
        tflowphys = np.array(tflowphys)
        tflowphys_err = np.array(tflowphys_err)

        hold_series = null_series
        hold_series['x_data'] = tflowphys
        hold_series['xerr_data'] = tflowphys_err
        hold_series['y_data'] = this_data.apply(lambda x : x.Avg)
        hold_series['yerr_data'] = this_data.apply(lambda x : x.Std)
        hold_series['label'] = ikey
        hold_series['shift'] = 0.0
        hold_series['type'] = 'error_bar'
        data_plot.AppendData(hold_series)
    data_plot.PlotAll()
    return data_plot


def PlotEpsWrap(this_data,file_name):
    print('plotting eps ratio for '+file_name)
    op1,op2 = file_name[0],file_name[1]
    this_info = pa.Series()
    this_info['save_file'] = graphdir+'/flowops_'+file_name+'_errrat_.pdf'
    this_info['title'] = r'$\frac{\Delta<'+op1+'(t+eps)'+op2+'(t_0)>}{\Delta<'+op1+'(t)'+op2+'(t_0)>}$'
    this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
    this_info['ylabel'] = 'ratio'
    data_plot = Plotting(plot_info=this_info)
    for ikey,idata in ErrEpsWrap(this_data).items():
        if len(idata.index) == 0:
            print('Warning, ',file_name,ikey,' did not pass properly, skipping')
            continue
        hold_series = null_series
        hold_series['x_data'] = 'from_keys'
        hold_series['y_data'] = idata
        hold_series['key_select'] = (slice(None),idata.index[0][1])
        hold_series['label'] = ikey
        hold_series['type'] = 'plot_vary'
        data_plot.AppendData(hold_series)
    data_plot.PlotAll()
    return data_plot


def ErrEpsWrap(data_list):
    data_out = {}
    for ikey,idata in data_list.items():
        data_out[ikey] = ErrEpsShift(idata.Op_Stats['boot'].apply(lambda x : np.abs(x.Std/x.Avg)))
    return data_out


def ErrEpsShift(this_data):
    vall,indexl = [],[]
    for ieps,ikey in enumerate(my_eps_list):
        indexl += [('eps'+ikey.replace('t_f',''),jkey) for jkey in this_data.keys()]
        vall += list((np.roll(this_data,-ieps)/this_data).values)
    if len(indexl)>0:
        index_list = pa.MultiIndex.from_tuples(indexl,names=['shift','flow_time'])
        out_series = pa.Series(vall,index=index_list)
    else:
        out_series = pa.Series()
    return out_series



import sys
master_ens_list = ['mpi411','mpi570','mpi701','L16','L20','L28','qu_L24','qu_L28','qu_L32']
t0_dict = {}
t0_dict['mpi411'] =     2.5377162290
t0_dict['mpi570'] =     2.3999592122
t0_dict['mpi701'] =     2.2590866059
t0_dict['L16'] =        1.3628954894
t0_dict['L20'] =        2.2386627909
t0_dict['L28'] =        4.9885618908
t0_dict['qu_L24'] =     3.2060027740
t0_dict['qu_L28'] =     4.4433376418
t0_dict['qu_L32'] =     6.0324704599


a_dict = {}
a_dict['mpi411'] =  0.0907
a_dict['mpi570'] =  0.0907
a_dict['mpi701'] =  0.0907
a_dict['L16'] =  0.1095
a_dict['L20'] =  0.0936
a_dict['L28'] =  0.0684
# these lattice spacing are fixed such that sqrt(8t_0) is 0.4 fm
a_dict['qu_L24'] =  0.07898289551435124
a_dict['qu_L28'] =  0.06709039362546686
a_dict['qu_L32'] =  0.057579434581241117
fit_info_list = []
fit_info = {}
# fit_info['Funs'] = [c3FitFun_nosqrt,3]
fit_info['Funs'] = [c3FitFun_nosqrt_log,4]
fit_info_list.append(fit_info)
tflow_fit_min = 0.1
tflow_fit_max = 0.4
min_par_tflow = 1
min_hold = tflow_fit_min
max_hold = tflow_fit_max
tflow_fit_min = {}
tflow_fit_max = {}
use_t0rat = True
if use_t0rat:
    sqrt8t_key = 'sqrt8t_t0rat'
else:
    sqrt8t_key = 'sqrt8t_fm'
if use_t0rat:
    print(min_hold,' to ',max_hold, ' in fm corresponds to')
    for ikey in a_dict.keys():
        print(ikey)
        tflow_fit_min[ikey] = min_hold/(a_dict[ikey]*np.sqrt(8*t0_dict[ikey]))
        tflow_fit_max[ikey] = max_hold/(a_dict[ikey]*np.sqrt(8*t0_dict[ikey]))
        print(tflow_fit_min[ikey],' to ',tflow_fit_max[ikey], ' in units of sqrt(8t_{0})')
        print()
else:
    for ikey in a_dict.keys():
        tflow_fit_min[ikey] = min_hold
        tflow_fit_max[ikey] = max_hold

run_nb = 1
if len(sys.argv)> 1:
    ens_list = []
    for iarg in sys.argv[1:]:
        if 'mpi' in iarg and 'ens' in iarg:
            ens_list.append(['mpi411','mpi570','mpi701'])
        elif 'latspace' in iarg and 'ens' in iarg:
            ens_list.append(['L16','L20','L28'])
        elif 'quenched' in iarg and 'ens' in iarg:
            ens_list.append(['qu_L24','qu_L28','qu_L32'])
        elif 'nblock' in iarg:
            run_nb = int(iarg.replace('nblock',''))
        elif iarg not in master_ens_list:
            print('Warning, ensemble name not found: ',iarg)
        else:
            ens_list.append(iarg)
    if len(ens_list) == 0:
        ens_list = master_ens_list
else:
    ens_list = master_ens_list
print()
print('Running over ensembles')
print(', '.join(ens_list))
print()

DefWipe=True
def FitRat(this_data,this_key,nbstr=''):
    fit_data = this_data.Op_Stats[['boot']].reset_index()
    this_file = scratchdir+'/fit_WQ_'+this_key+nbstr+'.py3p'
    this_fit = False
    if os.path.isfile(this_file) and not DefWipe:
        with open(this_file,'rb') as f:
            this_fit,dummy = pik.load(f)
    if this_fit is False or len(this_fit.Fit_Stats_fmt['Fit']) == 0:
        if use_t0rat:
            fit_data[sqrt8t_key] = fit_data['Flow Times'].apply(lambda x : np.sqrt(untflowstr(x)/t0_dict[this_key]))
        else:
            fit_data[sqrt8t_key] = fit_data['Flow Times'].apply(lambda x : np.sqrt(untflowstr(x))*a_dict[this_key])
        fit_data['boot'] = fit_data.apply(lambda x : x['boot']*untflowstr(x['Flow Times']),axis=1)
        fit_data['boot'].apply(lambda x : x.Stats())
        fit_data = fit_data.set_index(sqrt8t_key)['boot']
        fit_less_data = fit_data[::5]
        this_fit = sff.SetOfFitFuns(data = fit_less_data,name=this_key)
        this_fit.ScanRange(tflow_fit_min[ikey],tflow_fit_max[ikey],
                           fit_info_list,min_fit_len=min_par_tflow,from_xdata=True)
        this_fit.DoFits()
        this_fit.SortChi()
        this_fit.Fitr_Len_Sort()
        # this_fit.Cut_min_len(30)
        this_fit.RemoveFuns()
        with open(this_file,'wb') as f:
            pik.dump((this_fit,fit_less_data),f)
        this_fit.GetFuns()
    return this_fit


def PlotFlowTime(t0_scale,this_nb=1):
    if this_nb > 1:
        nbstr = '_nb'+str(this_nb)
    else:
        nbstr = ''

    t0_scale_str = f't_s{t0_scale:.3}'
    this_info = pa.Series()
    this_info['save_file'] = graphdir+'/flowops_test_NoE_'+t0_scale_str+nbstr+'.pdf'
    this_info['title'] = r'$\frac{t<Q(t_{s})W(t)>}{<Q(t_{s})^{2}>}$ $t_{s}='+t0_scale_str+ r' t_{0}$'+nbstr
    this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
    this_info['ylabel'] = 'ratio'
    data_plot_NoE = Plotting(plot_info=this_info)

    this_info = pa.Series()
    this_info['save_file'] = graphdir+'/flowops_test_'+t0_scale_str+nbstr+'.pdf'
    this_info['title'] = r'$\frac{<Q(t_{s})W(t)>}{t<E(t)><Q(t_{s})^{2}>}$ $t_{s}='+t0_scale_str+ r' t_{0}$'+nbstr
    this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
    this_info['ylabel'] = 'ratio'
    data_plot = Plotting(plot_info=this_info)


    # this_info = pa.Series()
    # this_info['save_file'] = graphdir+'/flowops_R28_'+t0_scale_str+'.pdf'
    # this_info['title'] = r'$\frac{<Q(t_{s})W(t)>}{t<E(t)><Q(t_{s})^{2}>_{28}}$ $t_{s}='+t0_scale_str+ r' t_{0}$'
    # this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
    # this_info['ylabel'] = 'ratio'
    # data_plot_28 = Plotting(plot_info=this_info)
    #
    # this_info = pa.Series()
    # this_info['save_file'] = graphdir+'/flowops_test_Q0_'+t0_scale_str+'.pdf'
    # this_info['title'] = r'$\frac{<Q(t_{s})W(t)>}{t<E(t)><Q(t_{s})Q(0)>}$ $t_{s}='+t0_scale_str+ r' t_{0}$'
    # this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
    # this_info['ylabel'] = 'ratio'
    # data_plot_Q0 = Plotting(plot_info=this_info)
    #
    this_info = pa.Series()
    this_info['save_file'] = graphdir+'/flowops_test_WWrat_'+t0_scale_str+nbstr+'.pdf'
    this_info['title'] = r'$\frac{<W(t_{s})W(t)>}{t<E(t)><Q(t_{s})W(t_{s})>}$ $t_{s}='+t0_scale_str+ r' t_{0}$'+nbstr
    this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
    this_info['ylabel'] = 'ratio'
    WW_data_plot = Plotting(plot_info=this_info)

    this_info = pa.Series()
    this_info['save_file'] = graphdir+'/flowops_test_WWrat_NoE_'+t0_scale_str+nbstr+'.pdf'
    this_info['title'] = r'$\frac{t<W(t_{s})W(t)>}{<Q(t_{s})W(t_{s})>}$ $t_{s}='+t0_scale_str+ r' t_{0}$'+nbstr
    this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
    this_info['ylabel'] = 'ratio'
    WW_data_plot_NoE = Plotting(plot_info=this_info)

    this_info = pa.Series()
    this_info['save_file'] = graphdir+'/flowops_test_'+t0_scale_str+nbstr+'_tauint.pdf'
    this_info['title'] = r'$\tau_{int}$ of alpha ratio $t_{s}='+t0_scale_str+ r' t_{0}$'+nbstr
    this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
    this_info['ylabel'] = '$\tau_{int}$'
    data_plot_tauint = Plotting(plot_info=this_info)

    this_info = pa.Series()
    this_info['save_file'] = graphdir+'/flowops_test_covar'+nbstr+'.pdf'
    this_info['title'] = 'Covariance matrix of alpha ratio '+nbstr
    this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
    this_info['ylabel'] = '$Covar$'
    data_plot_covar = Plotting(plot_info=this_info)
    #
    # this_info = pa.Series()
    # this_info['save_file'] = graphdir+'/flowops_W.pdf'
    # this_info['title'] = r'$<W(t)>$'
    # this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
    # this_info['ylabel'] = '$W$'
    # data_plot_W = Plotting(plot_info=this_info)
    #
    # this_info = pa.Series()
    # this_info['save_file'] = graphdir+'/flowops_Q.pdf'
    # this_info['title'] = r'$<Q(t)>$'
    # this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
    # this_info['ylabel'] = '$Q$'
    # data_plot_Q = Plotting(plot_info=this_info)
    #
    # this_info = pa.Series()
    # this_info['save_file'] = graphdir+'/flowops_E.pdf'
    # this_info['title'] = r'$<\frac{6*t^{2}}{V}*E(t)>$'
    # this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
    # this_info['ylabel'] = '$E$'
    # data_plot_E = Plotting(plot_info=this_info)
    #
    # this_info = pa.Series()
    # this_info['save_file'] = graphdir+'/flowops_QQ_'+t0_scale_str+'.pdf'
    # this_info['title'] = r'$<Q(t)^{2}>$'
    # this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
    # this_info['ylabel'] = r'$<QQ>$'
    # data_plot_QQ = Plotting(plot_info=this_info)
    #
    # this_info = pa.Series()
    # this_info['save_file'] = graphdir+'/flowops_QQ0_'+t0_scale_str+'.pdf'
    # this_info['title'] = r'$<Q(t)Q(0)>$'
    # this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
    # this_info['ylabel'] = r'$<QQ>$'
    # data_plot_QQ0 = Plotting(plot_info=this_info)
    #
    #
    #
    # this_info = pa.Series()
    # this_info['save_file'] = graphdir+'/flowops_QW_'+t0_scale_str+'.pdf'
    # this_info['title'] = r'$<Q(t_{s})W(t)>$ $t_{s}='+t0_scale_str+ r' t_{0}$'
    # this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
    # this_info['ylabel'] = r'$<QW>$'
    # data_plot_QW = Plotting(plot_info=this_info)

    # this_info = pa.Series()
    # this_info['save_file'] = graphdir+'/flowops_WW_'+t0_scale_str+'.pdf'
    # this_info['title'] = r'$<W(t_{s})W(t)>$ $t_{s}='+t0_scale_str+ r' t_{0}$'
    # this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
    # this_info['ylabel'] = r'$<WW>$'
    # data_plot_WW = Plotting(plot_info=this_info)

    Rat_list = {}
    Rat_28_list = {}
    Rat_Q0_list = {}
    Rat_NoE_list = {}
    data_Q_list = {}
    data_W_list = {}
    data_E_list = {}
    data_QW_list = {}
    data_QQ_list = {}
    data_QQ0_list = {}
    data_WW_list = {}
    # data_QW2_list = {}
    WWrat_ens_list = {}
    ts_list = {}
    Rat_Fit_list = {}
    for ikey,iens in ens_dict.items():
        if ikey not in ens_list: continue
        if '-' in ikey or '_no2' in ikey: continue
        parse_dict = defInfoFlow.copy()
        parse_dict.update(iens)
        parse_dict['nxyzt'] = [parse_dict['nxyzt'][0],parse_dict['nxyzt'][0],parse_dict['nxyzt'][0],parse_dict['nxyzt'][1]]
        parse_dict['kud'] = int(parse_dict['kud']*10**7)
        parse_dict['ks'] = int(parse_dict['ks']*10**7)
        print('         ANALYSING CONFIGURATION ',ikey)
        this_alpha,this_alpha_Q0,this_alpha_NoE,data_plot, \
        data_plot_NoE,data_plot_tauint,data_Q,data_E,data_W,data_QW, \
        data_QQ,data_QQ0,data_plot_covar,InfoFlow,its = AlphaRatAnalysis(ikey,parse_dict,data_plot,data_plot_NoE,
                                                            data_plot_tauint,data_plot_covar,t0_scale,nblock=this_nb)

        Rat_Fit_list[ikey] = FitRat(this_alpha_NoE,ikey,nbstr=nbstr)
        ts_list[ikey] = its
        Rat_list[ikey] = this_alpha
        Rat_Q0_list[ikey] = this_alpha_Q0
        Rat_NoE_list[ikey] = this_alpha_NoE
        data_Q_list[ikey] = data_Q
        data_E_list[ikey] = data_E
        data_W_list[ikey] = data_W
        data_QQ_list[ikey] = data_QQ
        data_QQ0_list[ikey] = data_QQ0
        data_QW_list[ikey] = deepcopy(data_QW)

        hold_series = null_series
        hold_series['type'] = 'fit_vary'
        # plot_fit = A0Fit_data[ikey].ReduceViaChi(nchi_cutoff=15)[0]
        plot_fit = Rat_Fit_list[ikey]
        hold_series['key_select'] = 'First'
        hold_series['fit_class'] = plot_fit
        hold_series['label'] = ikey+' Fit'
        # hold_series['xaxis'] = 0
        # hold_series['xdatarange'] = 'Data'
        # hold_series['otherXvals'] = [1]
        data_plot_NoE.AppendData(hold_series)

        # data_plot_Q0 = this_alpha_Q0.FlowPlot_mul_tf2(data_plot_Q0)
        # data_plot_Q = data_Q.FlowPlot(data_plot_Q)
        # data_plot_W = data_W.FlowPlot(data_plot_W)
        # data_plot_E = data_E.FlowPlot_mul_tf2(data_plot_E,just_E=True)
        # data_plot_QQ = data_QQ.FlowPlot(data_plot_QQ)
        # data_plot_QQ0 = data_QQ0.FlowPlot(data_plot_QQ0)
        # data_plot_QW = data_QW.FlowPlot_mul_tf2(data_plot_QW,mul=True)
        if DoWW:
            this_WWrat,this_WWrat_NoE,WW_data_plot,WW_data_plot_NoE,data_WW = DoWWAnalysis(ikey,WW_data_plot,WW_data_plot_NoE,data_Q,data_E,data_W,data_QW,InfoFlow,t0_scale)
            data_WW_list[ikey] = data_WW
            # data_QW2_list[ikey] = data_QW2
            WWrat_ens_list[ikey] = this_WWrat
            # data_plot_WW = data_WW.FlowPlot(data_plot_WW)
            # data_plot_QWcomp = data_QW.FlowPlot(data_plot_QWcomp)
            # # data_plot_QWcomp = data_QW2.FlowPlot(data_plot_QWcomp)
            # data_plot_WWvsQW = data_QW.FlowPlot_mul_tf2(data_plot_WWvsQW,mul=True)
            # data_plot_WWvsQW = data_WW.FlowPlot_mul_tf2(data_plot_WWvsQW,mul=True)
    if 'L28' in data_QQ_list.keys() and 'qu_' not in data_QQ_list.keys():
        QQ_28 = data_QQ_list['L28'].GetInterp(ts_list['L28'])
    else: QQ_28 = None
    leg_lab_ext = r' \frac{t_{s}}{t_{0}}='+t0_scale_str.replace('t_s','')
    for ikey,iens in ens_dict.items():
        if ikey not in ens_list: continue
        if '-' in ikey or '_no2' in ikey: continue
        if 'mpi' in ikey: continue
        if QQ_28 is not None:
            this_data = data_QW_list[ikey].__div__(data_E_list[ikey]*QQ_28)
            out_name = 'Test_alpha_28_'+ikey+'_'+t0_scale_str+nbstr
            this_data.FlowSetCustomName(string=out_name)
            this_data.FlowSetCustomName(string=out_name,stringLL=this_data.flowLegLab+leg_lab_ext)
            this_data.FlowStats()
            this_data.FlowWrite()
            # data_plot_28 = this_data.FlowPlot_mul_tf2(data_plot_28)
            Rat_28_list[ikey] = this_data
    # data_plot_28.PlotAll()
    data_plot.PlotAll()
    # data_plot_Q0.PlotAll()
    data_plot_NoE.PlotAll()
    # data_plot_tauint.PlotAll()
    # data_plot_covar.PlotAll()
    # data_plot_W.PlotAll()
    # data_plot_E.PlotAll()
    # data_plot_Q.PlotAll()
    # data_plot_QQ.PlotAll()
    # data_plot_QQ0.PlotAll()
    # data_plot_QW.PlotAll()
    # data_plot_WW.PlotAll()
    # data_plot_QWcomp.PlotAll()

    if DoWW:
        # data_plot_WWvsQW.PlotAll()
        WW_data_plot.PlotAll()
        WW_data_plot_NoE.PlotAll()
    # data_plot_QQ_errrat = PlotEpsWrap(data_QQ_list,'QQ')
    # data_plot_QW_errrat = PlotEpsWrap(data_QW_list,'QW')
    # data_plot_WW_errrat = PlotEpsWrap(data_WW_list,'WW')
    # data_plot_QW2_errrat = PlotEpsWrap(data_QW2_list,'QW2')
    # data_plot_Et2 = PlotEt2(data_E_list)
    # data_plot_WWdivQW = PlotErrRat(data_WW_list,data_QW_list)
    out_dict = {}
    out_dict['Rat_Fit_data'] = Rat_Fit_list
    out_dict['Rat_data'] = Rat_list
    out_dict['Rat_28_data'] = Rat_28_list
    out_dict['Rat_Q0_data'] = Rat_Q0_list
    out_dict['Rat_NoE_data'] = Rat_NoE_list
    out_dict['E_data'] = data_E_list
    out_dict['Q_data'] = data_Q_list
    out_dict['W_data'] = data_W_list
    out_dict['QQ_data'] = data_QQ_list
    out_dict['QQ0_data'] = data_QQ0_list
    out_dict['QW_data'] = data_QW_list
    out_dict['WW_data'] = data_WW_list
    # out_dict['QW2_data'] = data_QW2_list
    out_dict['data_plot'] = data_plot
    out_dict['data_plot_NoE'] = data_plot_NoE
    out_dict['data_plot_tauint'] = data_plot_tauint
    out_dict['data_plot_covar'] = data_plot_covar
    # out_dict['data_plot_W'] = data_plot_W
    # out_dict['data_plot_E'] = data_plot_E
    # out_dict['data_plot_Q'] = data_plot_Q
    # out_dict['data_plot_QQ'] = data_plot_QQ
    # out_dict['data_plot_QQ0'] = data_plot_QQ0
    # out_dict['data_plot_QW'] = data_plot_QW
    # out_dict['data_plot_WW'] = data_plot_WW
    # out_dict['data_plot_QWcomp'] = data_plot_QWcomp
    # out_dict['data_plot_WWvsQW'] = data_plot_WWvsQW
    # out_dict['WW_data_plot'] = WW_data_plot
    # out_dict['WW_data_plot_NoE'] = WW_data_plot_NoE
    # out_dict['data_plot_Et2'] = data_plot_Et2
    # out_dict['data_plot_WWdivQW'] = data_plot_WWdivQW
    return out_dict

## sqrt(8t_s) = "sqrt_t0_scale_list" X sqrt(8t_0)
# sqrt_t0_scale_list = [0.2, 0.4, 0.6, 0.8, 1.0]
sqrt_t0_scale_list = [0.6]
t0_scale_list = [it0**2 for it0 in sqrt_t0_scale_list]
pow_list = 3,2,3/2
all_data = {}
for it0s in t0_scale_list:
    t0_scale_str = f't_s{it0s:.3f}'
    data_list = PlotFlowTime(it0s,this_nb=run_nb)
    all_data[t0_scale_str] = data_list

# for ikey,iens in ens_dict.items():
#     if ikey not in ens_list: continue
#     if '-' in ikey or '_no2' in ikey: continue
#
#     this_info = pa.Series()
#     this_info['save_file'] = graphdir+'/flowops_test_NoE_'+ikey+'.pdf'
#     this_info['title'] = r'$\frac{t<Q(t_{s})W(t)>}{<Q(t_{s})^{2}>}$ '+ikey
#     this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
#     this_info['ylabel'] = 'ratio'
#     data_plot_NoE = Plotting(plot_info=this_info)
#
#     this_info = pa.Series()
#     this_info['save_file'] = graphdir+'/flowops_test_'+ikey+'.pdf'
#     this_info['title'] = r'$\frac{<Q(t_{s})W(t)>}{t<E(t)><Q(t_{s})^{2}>}$ '+ikey
#     this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
#     this_info['ylabel'] = 'ratio'
#     data_plot = Plotting(plot_info=this_info)
#
#     # this_info = pa.Series()
#     # this_info['save_file'] = graphdir+'/flowops_R28_'+ikey+'.pdf'
#     # this_info['title'] = r'$\frac{<Q(t_{s})W(t)>}{t<E(t)><Q(t_{s})^{2}>_{28}}$ '+ikey
#     # this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
#     # this_info['ylabel'] = 'ratio'
#     # data_plot_28 = Plotting(plot_info=this_info)
#
#     this_info = pa.Series()
#     this_info['save_file'] = graphdir+'/flowops_test_Q0_'+ikey+'.pdf'
#     this_info['title'] = r'$\frac{<Q(t_{s})W(t)>}{t<E(t)><Q(t_{s})Q(0)>}$ '+ikey
#     this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
#     this_info['ylabel'] = 'ratio'
#     data_plot_Q0 = Plotting(plot_info=this_info)
#
#     # this_info = pa.Series()
#     # this_info['save_file'] = graphdir+'/flowops_test_WWrat_'+ikey+'.pdf'
#     # this_info['title'] = r'$\frac{<W(t_{s})W(t)>}{t<E(t)><Q(t_{s})W(t_{s})>}$ '+ikey
#     # this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
#     # this_info['ylabel'] = 'ratio'
#     # WW_data_plot = Plotting(plot_info=this_info)
#     #
#     # this_info = pa.Series()
#     # this_info['save_file'] = graphdir+'/flowops_test_WWrat_NoE_'+ikey+'.pdf'
#     # this_info['title'] = r'$\frac{t<W(t_{s})W(t)>}{<Q(t_{s})W(t_{s})>}$ '+ikey
#     # this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
#     # this_info['ylabel'] = 'ratio'
#     # WW_data_plot_NoE = Plotting(plot_info=this_info)
#
#     this_info = pa.Series()
#     this_info['save_file'] = graphdir+'/flowops_QWdiv28_'+ikey+'.pdf'
#     this_info['title'] = r'$<Q(t_{s})W(t)>_{'+ikey+'}/<Q(t_{s})W(t)>_{L28}$ '
#     this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
#     this_info['ylabel'] = r'$<QW>/<QW>_{28}$'
#     data_plot_QW_rat28 = Plotting(plot_info=this_info)
#
#     this_info = pa.Series()
#     this_info['save_file'] = graphdir+'/flowops_QWdiv20_'+ikey+'.pdf'
#     this_info['title'] = r'$<Q(t_{s})W(t)>_{'+ikey+'}/<Q(t_{s})W(t)>_{L20}$ '
#     this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
#     this_info['ylabel'] = r'$<QW>/<QW>_{20}$'
#     data_plot_QW_rat20 = Plotting(plot_info=this_info)
#
#
#     this_info = pa.Series()
#     this_info['save_file'] = graphdir+'/flowops_QWdiv16_'+ikey+'.pdf'
#     this_info['title'] = r'$<Q(t_{s})W(t)>_{'+ikey+'}/<Q(t_{s})W(t)>_{L16}$ '
#     this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
#     this_info['ylabel'] = r'$<QW>/<QW>_{16}$'
#     data_plot_QW_rat16 = Plotting(plot_info=this_info)
#
#     this_info = pa.Series()
#     this_info['save_file'] = graphdir+'/flowops_QW_'+ikey+'.pdf'
#     this_info['title'] = r'$<Q(t_{s})W(t)>$ '+ikey
#     this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
#     this_info['ylabel'] = r'$<QW>$'
#     data_plot_QW = Plotting(plot_info=this_info)
#
#
#     # this_info = pa.Series()
#     # this_info['save_file'] = graphdir+'/flowops_WW_'+ikey+'.pdf'
#     # this_info['title'] = r'$<W(t_{s})W(t)>$ '+ikey
#     # this_info['xlabel'] = r'$\sqrt{8t}/\sqrt{8t_{0}}$'
#     # this_info['ylabel'] = r'$<WW>$'
#     # data_plot_WW = Plotting(plot_info=this_info)
#
#
#     for itc,(its,ts_data) in enumerate(all_data.items()):
#         # if ikey in ts_data['Rat_28_data'].keys():
#         #     this_rat_28 = ts_data['Rat_28_data'][ikey]
#         #     data_plot_28 = this_rat_28.FlowPlot_mul_tf2(data_plot_28)
#
#         this_QW = ts_data['QW_data'][ikey]
#         # this_WW = ts_data['WW_data'][ikey]
#         this_rat = ts_data['Rat_data'][ikey]
#         this_rat_NoE = ts_data['Rat_NoE_data'][ikey]
#         this_rat_Q0 = ts_data['Rat_Q0_data'][ikey]
#         if 'L28' in ts_data['QW_data'].keys():
#             this_QW_rat = this_QW.__div__(ts_data['QW_data']['L28'])
#             this_QW_rat.FlowStats()
#             data_plot_QW_rat28 = this_QW_rat.FlowPlot(data_plot_QW_rat28,lab_append=its)
#         if 'L20' in ts_data['QW_data'].keys():
#             this_QW_rat = this_QW.__div__(ts_data['QW_data']['L20'])
#             this_QW_rat.FlowStats()
#             data_plot_QW_rat20 = this_QW_rat.FlowPlot(data_plot_QW_rat20,lab_append=its)
#         if 'L16' in ts_data['QW_data'].keys():
#             this_QW_rat = this_QW.__div__(ts_data['QW_data']['L16'])
#             this_QW_rat.FlowStats()
#             data_plot_QW_rat16 = this_QW_rat.FlowPlot(data_plot_QW_rat16,lab_append=its)
#         data_plot_QW = this_QW.FlowPlot(data_plot_QW)
#         if 'L28' in ikey:
#             for ipow in pow_list:
#                 data_plot_QW = this_QW.FlowPlot_mul_tf2(data_plot_QW,mul=True,tpow=ipow,lab_append='$<QW>^{'+str(ipow)+'}$')
#                 data_plot = this_rat.FlowPlot_mul_tf2(data_plot,tpow=ipow,lab_append='$<QW>^{'+str(ipow)+'}$')
#         else:
#             data_plot_QW = this_QW.FlowPlot_mul_tf2(data_plot_QW,mul=True)
#             data_plot = this_rat.FlowPlot_mul_tf2(data_plot)
#         data_plot_NoE = this_rat_NoE.FlowPlot_mul_tf2(data_plot_NoE,mul=True)
#         data_plot_Q0 = this_rat_Q0.FlowPlot_mul_tf2(data_plot_Q0)
#     data_plot_QW_rat28.PlotAll()
#     data_plot_QW_rat20.PlotAll()
#     data_plot_QW_rat16.PlotAll()
#     data_plot_QW.PlotAll()
#     data_plot.PlotAll()
#     # data_plot_28.PlotAll()
#     data_plot_NoE.PlotAll()
#     data_plot_Q0.PlotAll()
