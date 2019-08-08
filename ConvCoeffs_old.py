

import pickle as pik
import numpy as np
from MomParams import hbarc
hbarcM = hbarc*1000
from Params import this_dir
import pandas as pa
import matplotlib.pyplot as pl

#### NEUTRON ####
mN = 939.5653 # neutron mass in MeV



chi_lab = r'$\chi^{2}_{pdf}$'


def ReadCs(this_file):
    this_data = pik.load(open(this_file,'rb'))

    this_fit = this_data['plot_data']
    par_list,name_list = [],[]
    chi_list = []
    chi_list_f3 = []
    fun_list = []
    for ifit,fit_data in this_fit.iteritems():
        if 'fit_class' in fit_data.keys() and hasattr(fit_data['fit_class'],'iteritems'):
            for ikey,iffs in fit_data['fit_class'].iteritems():
                if '$F(0)$' not in ikey: continue
                if 'frac{tilde{F}_{3}}{2m_{N}}' not in ikey and \
                    'frac{F_{3}}{2m_{N}}' not in ikey: continue
                if 'fitr0-4' not in ikey: continue
                par_list.append(iffs.fit_data['Params'].values)
                chi_list_f3.append(CreateChi2f3(iffs))
                chi_list.append(iffs.fit_data['Chi2DoF'].values[0])
                this_lab = fit_data['label']
                while this_lab in name_list:
                    this_lab += '_1'
                name_list.append(this_lab)
                iffs.GetFuns()
                fun_list.append(iffs.Fun)
    out_df = pa.DataFrame()
    for C_list,iname,ifun,ichi,ichi3 in zip(par_list,name_list,fun_list,chi_list,chi_list_f3):
        if len(C_list) == 5:
            C1,C2,C3,C4,C5 = C_list
        else:
            C1,C2,C3 = C_list
        this_series = pa.Series()
        print('Parameters for ',iname)
        print('Uses fit function:',ifun.__name__)
        print('Chi2 pdf = ',ichi.MakeValAndErr(Dec=2))
        print('Chi2 pdf first 3 points = ',ichi3.MakeValAndErr(Dec=2))
        # C1,C2,C3 = par_list.values
        this_series[chi_lab] = ichi
        this_series[chi_lab+ ' first 3'] = ichi3
        this_series['C1'] = C1
        this_series['C2'] = C2
        this_series['C3'] = C3
        if len(C_list) == 5:
            this_series['C4'] = C4
            this_series['C5'] = C5
        print(' In fm MeV^-2 (first table)')
        if len(C_list) == 5:
            print('     C1_n: '+C1.MakeValAndErr(Dec=2))
            print('     C1_p: '+C2.MakeValAndErr(Dec=2))
            print('     C2: '+C3.MakeValAndErr(Dec=2))
            print('     C3_n: '+C4.MakeValAndErr(Dec=2))
            print('     C3_p: '+C5.MakeValAndErr(Dec=2))
        else:
            print('     C1: '+C1.MakeValAndErr(Dec=2))
            print('     C2: '+C2.MakeValAndErr(Dec=2))
            print('     C3: '+C3.MakeValAndErr(Dec=2))




        this_series['C1_fm'] = C1 *hbarcM**2
        this_series['C2_fm'] = C2 *hbarcM**2
        if len(C_list) == 5:
            this_series['C3_fm'] = C3 *hbarcM**2

        print()
        print(' In fm^3 (second table)')
        if len(C_list) == 5:
            print('     C1_n: '+this_series['C1_fm'].MakeValAndErr(Dec=2))
            print('     C1_p: '+this_series['C2_fm'].MakeValAndErr(Dec=2))
            print('     C2: '+this_series['C3_fm'].MakeValAndErr(Dec=2))
            print('     C3_n: '+C4.MakeValAndErr(Dec=2))
            print('     C3_p: '+C5.MakeValAndErr(Dec=2))
        else:
            print('     C1: '+this_series['C1_fm'].MakeValAndErr(Dec=2))
            print('     C2: '+this_series['C2_fm'].MakeValAndErr(Dec=2))
            print('     C3: '+this_series['C3'].MakeValAndErr(Dec=2))
        print()
        out_df[iname] = this_series
    return out_df,fun_list

def CreateChi2f3(this_fit,points=3,dof=3):
    this_fit.GetFuns()
    this_chi_fun = this_fit.ChiFun
    this_params = this_fit.fit_data['Params'].values
    this_chi_list = []
    for ix,iy,iyerr in zip((this_fit.xdata.T)[:points],this_fit.ydata[:points],this_fit.yerr[:points]):
        # print()
        # print('DEBUG')
        # print(list(ix)+[iy,iyerr])
        # print(this_params)
        # print(this_chi_fun(this_params,list(ix)+[iy,iyerr]))
        # print()
        this_chi_list.append(this_chi_fun(this_params,list(ix)+[iy,iyerr]))
    return sum(np.array(this_chi_list)**2)/dof

def ReadCs_Schiff(this_file):
    this_data = pik.load(open(this_file,'rb'))

    this_fit = this_data['plot_data']
    par_list,name_list = [],[]
    chi_list_f3 = []
    chi_list = []
    fun_list = []
    for ifit,fit_data in this_fit.iteritems():
        if 'fit_class' in fit_data.keys() and hasattr(fit_data['fit_class'],'iteritems'):
            for ikey,iffs in fit_data['fit_class'].iteritems():
                if "$Q^{2}F^{'}(0)$" not in ikey: continue
                if 'frac{tilde{F}_{3}}{2m_{N}}' not in ikey and \
                    'frac{F_{3}}{2m_{N}}' not in ikey: continue
                if 'fitr0-4' not in ikey: continue
                par_list.append(iffs.fit_data['Params'].values)
                chi_list.append(iffs.fit_data['Chi2DoF'].values[0])
                chi_list_f3.append(CreateChi2f3(iffs))
                this_lab = fit_data['label']
                while this_lab in name_list:
                    this_lab += '_1'
                name_list.append(this_lab)
                iffs.GetFuns()
                fun_list.append(iffs.Fun)
    out_df = pa.DataFrame()
    for C_list,iname,ifun,ichi,ichi3 in zip(par_list,name_list,fun_list,chi_list,chi_list_f3):
        C5,C4 = C_list
        this_series = pa.Series()
        print('Parameters for ',iname)
        print('Uses fit function:',ifun.__name__)
        print('Chi2 pdf = ',ichi.MakeValAndErr(Dec=2))
        print('Chi2 pdf first 3 points = ',ichi3.MakeValAndErr(Dec=2))
        # C1,C2,C3 = par_list.values
        this_series[chi_lab] = ichi
        this_series[chi_lab+ ' first 3'] = ichi3
        this_series['C4'] = C4
        this_series['C5'] = C5
        print(' In fm GeV^-2 (first table)')
        print('     C4: '+C4.MakeValAndErr(Dec=2))
        print('     C5: '+C5.MakeValAndErr(Dec=2))

        this_series['C4_fm'] = C4 *hbarc**2
        this_series['C5_fm'] = C5 *hbarc**2
        print()
        print(' In fm^3 (second table)')
        print('     C4: '+this_series['C4_fm'].MakeValAndErr(Dec=2))
        print('     C5: '+this_series['C5_fm'].MakeValAndErr(Dec=2))
        print()
        out_df[iname] = this_series
    return out_df,fun_list

fix7_file = '/home/jackdra/LQCD/Results/HPCCRES_FixedG2/graphs/FF/CorrectSystem//d_mpi_fix700.py3p'
print('------d Fixed 700 MeV for L ens file------')
fix7_res,fix7_funs= ReadCs(fix7_file)


n_file = '/home/jackdra/LQCD/Results/HPCCRES_FixedG2/graphs/FF/CorrectSystem/mpi3fit/d_n_mpi3.py3p'
print('------d_n_mpi3 file------')
n_res,n_funs= ReadCs(n_file)

p_file = '/home/jackdra/LQCD/Results/HPCCRES_FixedG2/graphs/FF/CorrectSystem/mpi3fit/d_p_mpi3.py3p'
print('------d_p_mpi3 file-------')
p_res,p_funs = ReadCs(p_file)

np_file = '/home/jackdra/LQCD/Results/HPCCRES_FixedG2/graphs/FF/CorrectSystem/mpi3fit/d_mpi3_NPcomb.py3p'
print('------NP fit combine, Common C2------')
np_res,np_funs = ReadCs(np_file)
np_res.loc[chi_lab,:] = np_res.loc[chi_lab,:]*9/7


npnm_file = '/home/jackdra/LQCD/Results/HPCCRES_FixedG2/graphs/FF/CorrectSystem/mpi3fit/d_mpi3_NPcomb_nompi3.py3p'
print('------NP fit combine, Common C2  No mpi3------')
npnm_res,npnm_funs = ReadCs(npnm_file)

Schiff_file = '/home/jackdra/LQCD/Results/HPCCRES_FixedG2/graphs/FF/CorrectSystem/S_mpi.py3p'
print('------Shiff Moment Results------')
ReadCs_Schiff(Schiff_file)

# npMC2_file = '/home/jackdra/LQCD/Results/HPCCRES_FixedG2/graphs/FF/CorrectSystem/mpi3fit/d_mpi3_NPcomb_mc2.py3p'
# print('------NP fit combine, Common C2, same sign------')
# npMC2_res,npMC2_funs = ReadCs(npMC2_file)

# prl_file = '/home/jackdra/LQCD/Scripts/EDM_paper/PRL_graphs/FF/EDM_mpi.py3p'
# print('------old EDM file (PRL)------')
# prl_res,prl_funs = ReadCs(prl_file)


funtolab = {
    'ContPlLogFitFun_x0_N':     ' $m_{\pi}^{3}$ term',
    'ContPlLogFitFun_x0_P':     ' $m_{\pi}^{3}$ term',
    'ContPlLogFitFun_x0_dN':    ' no $m_{\pi}^{3}$ term',
    'ContPlLogFitFun_NP':       ' Combine Proton and Neutron',
    'ContPlLogFitFun_NP_nompi3':       ' Combine Proton and Neutron',
    'ContPlLogFitFun_NP_latmN': ' Using lattice $m_{N}$',
    'ContPlLogFitFun_x0':       ' Not dividing by any $m_{N}$ (Old Fit)'
}

def GetLabels(this_funs,this_df,pref='',post = ''):
    new_cols = []
    neg_C2_list = []
    for ifun in this_funs:
        if isinstance(ifun,str): continue
        new_cols.append(pref + funtolab[ifun.__name__]+post)
        if 'ContPlLogFitFun_x0_N' == ifun.__name__:
            neg_C2_list.append(-1)
        else:
            neg_C2_list.append(1)
    this_df.columns = new_cols
    for icol in new_cols:
        if 'ContPlLogFitFun_x0_N' in icol:
            this_df[icol]['C2'] = -this_df[icol]['C2']



GetLabels(fix7_funs,fix7_res,pref='Fixed mpi700')
GetLabels(n_funs,n_res,pref='Neutron ')
GetLabels(p_funs,p_res,pref='Proton ')
GetLabels(np_funs,np_res)
GetLabels(npnm_funs,npnm_res)

# GetLabels(prl_funs,prl_res,post=' Same C2 sign')
output_columns = [ 'C1 $[e fm^{3}]$', 'C2 $[e fm^{3}]$','C3 $[e fm^{-1}]$' ,chi_lab, chi_lab + ' first 3' ]

def FormatForLatex(this_res):
    if 'C4' in this_res.index:
        Neutron_res = this_res.loc[['C1_fm','C3_fm','C4',chi_lab, chi_lab + ' first 3'],:].T
        Proton_res = this_res.loc[['C2_fm','C3_fm','C5',chi_lab, chi_lab + ' first 3'],:].T
        Neutron_res.columns = output_columns
        Proton_res.columns = output_columns
        Neutron_res.loc[:,'C2 $[e fm^{3}]$'] = -Neutron_res.loc[:,'C2 $[e fm^{3}]$']
        return Neutron_res,Proton_res
    else:
        output_res = this_res.loc[['C1_fm','C2_fm','C3',chi_lab, chi_lab + ' first 3'],:].T
        output_res.columns = output_columns
        return output_res

fmt_fix7res =   FormatForLatex(fix7_res)
fmt_nres =   FormatForLatex(n_res)
fmt_pres =   FormatForLatex(p_res)
fmt_nres_comb,fmt_pres_comb =  FormatForLatex(np_res)
fmt_nres_comb_nm,fmt_pres_comb_nm =  FormatForLatex(npnm_res)
# fmt_prlres = FormatForLatex(prl_res)
fmt_pres = fmt_pres.iloc[[0,2],:]
fmt_nres_comb,fmt_pres_comb = fmt_nres_comb.iloc[[0],:],fmt_pres_comb.iloc[[0],:]
fmt_nres_comb_nm,fmt_pres_comb_nm = fmt_nres_comb_nm.iloc[[0],:],fmt_pres_comb_nm.iloc[[0],:]
fmt_no_mpi3 = fmt_pres.iloc[[1],:].append(fmt_nres.iloc[[1],:])
fmt_mpi3 = fmt_pres.iloc[[0],:].append(fmt_nres.iloc[[0],:])
fmt_res_comb = fmt_pres_comb.append(fmt_nres_comb)
fmt_res_comb_nm = fmt_pres_comb_nm.append(fmt_nres_comb_nm)
fmt_res_comb.index = ['Combined Fit, Proton Cs','Combined Fit, Neutron Cs']
fmt_res_comb_nm.index = ['Comb Fit (no mpi3), Proton Cs','Comb Fit (no mpi3), Neutron Cs']

fmt_fix7res.applymap(lambda x : x.MakeValAndErr(Dec=2))
fmt_no_mpi3.applymap(lambda x : x.MakeValAndErr(Dec=2))
fmt_mpi3.applymap(lambda x : x.MakeValAndErr(Dec=2))
fmt_res_comb.applymap(lambda x : x.MakeValAndErr(Dec=2))
fmt_res_comb_nm.applymap(lambda x : x.MakeValAndErr(Dec=2))
print(fmt_fix7res.applymap(lambda x : '$'+x.MakeValAndErr(Dec=2,latex=True)+'$').to_latex(escape=False))
print(fmt_no_mpi3.applymap(lambda x : '$'+x.MakeValAndErr(Dec=2,latex=True)+'$').to_latex(escape=False))
print(fmt_mpi3.applymap(lambda x : '$'+x.MakeValAndErr(Dec=2,latex=True)+'$').to_latex(escape=False))
print(fmt_res_comb.applymap(lambda x : '$'+x.MakeValAndErr(Dec=2,latex=True)+'$').to_latex(escape=False))
print(fmt_res_comb_nm.applymap(lambda x : '$'+x.MakeValAndErr(Dec=2,latex=True)+'$').to_latex(escape=False))


def plot_FF(this_series,square_mpi=False,fn=this_dir+'/TestGraphs/testCs.pdf',title=''):
    C1,C2,C3,chi2,chi2f = this_series.values
    xvals = pa.Series(np.arange(700))
    y_C1 = xvals.apply(lambda x : C1*x**2/hbarcM**2)
    y_C2 = xvals.apply(lambda x : C2*(x**2)*np.log(x**2/mN**2)/hbarcM**2)
    y_C12 = y_C1 + y_C2
    xplot = xvals
    if square_mpi:
        xplot = xplot**2
    y_C1_Avg = y_C1.apply(lambda x : x.GetAvg())
    y_C2_Avg = y_C2.apply(lambda x : x.GetAvg())
    y_C12_Avg = y_C12.apply(lambda x : x.GetAvg())
    y_C1_Up = y_C1.apply(lambda x : x.GetAvg()+x.GetStd())
    y_C2_Up = y_C2.apply(lambda x : x.GetAvg()+x.GetStd())
    y_C12_Up = y_C12.apply(lambda x : x.GetAvg()+x.GetStd())
    y_C1_Down = y_C1.apply(lambda x : x.GetAvg()-x.GetStd())
    y_C2_Down = y_C2.apply(lambda x : x.GetAvg()-x.GetStd())
    y_C12_Down = y_C12.apply(lambda x : x.GetAvg()-x.GetStd())
    pl.plot(xvals,y_C1_Avg.values,color='blue',label='$C_{1}m_{\pi}^{2}$')
    pl.plot(xvals,y_C2_Avg.values,color='red',label='$C_{2}m_{\pi}^{2}log(m_{\pi}^{2}/m_{N,phys}^{2})$')
    pl.plot(xvals,y_C12_Avg.values,color='green',label='blue+red')
    pl.fill_between(xvals,y_C1_Up.values,y_C1_Down.values,color='blue',alpha=0.3)
    pl.fill_between(xvals,y_C2_Up.values,y_C2_Down.values,color='red',alpha=0.3)
    pl.fill_between(xvals,y_C12_Up.values,y_C12_Down.values,color='green',alpha=0.3)
    pl.legend()
    pl.grid()
    pl.savefig(fn)
    pl.clf()
plot_FF(fmt_no_mpi3.iloc[0,:],fn=this_dir+'/TestGraphs/testCs_proton.pdf',title='Proton')
plot_FF(fmt_no_mpi3.iloc[1,:],fn=this_dir+'/TestGraphs/testCs_neutron.pdf',title='Neutron')
