

import pickle as pik
import numpy as np
from MomParams import hbarc
hbarcM = hbarc*1000
import pandas as pa
import matplotlib.pyplot as pl
from Params import this_dir,plot_style

#### NEUTRON ####
mN = 939.5653 # neutron mass in MeV

# tpad,tsize,xsize,ysize = 20,40,60,60
# legcol,legsize,legalpha = 1,25,0.7
# linewidth,dotsize,errorbarlen = 2,10,6
# def_tick_size = 30
# def_grid = True
# params = {'legend.fontsize': 25,
#           'xtick.labelsize':def_tick_size,
#           'ytick.labelsize':def_tick_size,
#           'legend.numpoints': 1,
#           'axes.labelsize' : 35,
#           'axes.titlesize' : 35,
#           'figure.autolayout': True,
#           'axes.grid': def_grid,
#           'errorbar.capsize': errorbarlen,
#           'lines.markersize': dotsize,
#           'lines.linewidth': linewidth,
#           'font.size': 30,
#           # 'lines.markeredgewidth':2,
#           'axes.xmargin':0.01,
#           'axes.ymargin':0.01}
# params = {'errorbar.capsize': 5}
pl.style.use(plot_style)
# pl.rcParams.update(params)

chi_lab = r'$\chi^{2}_{pdf}$'


def ReadCs(this_file,F3neg=True):
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
        if F3neg:
            C1,C2,C3 = -C1,-C2,-C3
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


funtolab = {
    'ContPlLogFitFun_x0_N':     ' $m_{\pi}^{3}$ term',
    'ContPlLogFitFun_x0_P':     ' $m_{\pi}^{3}$ term',
    # 'ContPlLogFitFun_x0_dN':    ' no $m_{\pi}^{3}$ term',
    'ContPlLogFitFun_x0_dN':    ' EDM fit',
    'ContPlLogFitFun_NP':       ' Combine Proton and Neutron',
    'ContPlLogFitFun_NP_nompi3':       ' Combine Proton and Neutron',
    'ContPlLogFitFun_NP_latmN': ' Using lattice $m_{N}$',
    'ContPlLogFitFun_x0':       ' Not dividing by any $m_{N}$ (Old Fit)',
    'SchiffFitFun':             ' Schiff Fit'
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


# npMC2_file = '/home/jackdra/LQCD/Results/HPCCRES_FixedG2/graphs/FF/CorrectSystem/mpi3fit/d_mpi3_NPcomb_mc2.py3p'
# print('------NP fit combine, Common C2, same sign------')
# npMC2_res,npMC2_funs = ReadCs(npMC2_file)

# prl_file = '/home/jackdra/LQCD/Scripts/EDM_paper/PRL_graphs/FF/EDM_mpi.py3p'
# print('------old EDM file (PRL)------')
# prl_res,prl_funs = ReadCs(prl_file)


EDM_cols = ['$C_{1}$ $\\left[\\bar{\\theta}\\, e\\text{~fm}^{3}\\right]$',
            '$C_{2}$ $\\left[\\bar{\\theta}\\, e\\text{~fm}^{3}\\right]$',
            '$C_{3}$ $\\left[\\bar{\\theta}\\, e\\text{~fm}^{-1} \\right]$' ,'EDM '+chi_lab ]

Schiff_cols = ['$C_{4}$ $\\left[\\bar{\\theta}\\, e\\text{~fm}^{3}\\right]$',
               '$C_{5}$ $\\left[\\bar{\\theta}\\, e\\text{~fm}\\right]$' ,'Schiff '+chi_lab ]

def Get_Schiff_res(Schiff_res,Schiff_funs,swap_np=False):
    GetLabels(Schiff_funs,Schiff_res)
    Schiff_res = Schiff_res.T
    Schiff_res = Schiff_res[['C4_fm','C5_fm','$\chi^{2}_{pdf}$']]
    Schiff_res.columns = Schiff_cols
    if swap_np:
        Schiff_res = Schiff_res.iloc[::-1,:]
    Schiff_res.index = ['neutron','proton']
    return Schiff_res


def Get_EDM_res(EDM_res,EDM_funs,swap_np=False):
    GetLabels(EDM_funs,EDM_res)
    EDM_res = EDM_res.T
    EDM_res = EDM_res[['C1_fm','C2_fm','C3','$\chi^{2}_{pdf}$']]
    EDM_res.columns = EDM_cols
    if swap_np:
        EDM_res = EDM_res.iloc[::-1,:]
    EDM_res.index = ['proton','neutron']
    return EDM_res


def ReadAndFormat_EDM(EDM_file,swap_np=False):
    print('------EDM fit read Cs------')
    EDM_res,EDM_funs = ReadCs(EDM_file)
    return Get_EDM_res(EDM_res,EDM_funs,swap_np=swap_np)


def ReadAndFormat_Schiff(Schiff_file,swap_np=False):
    print('------Schiff fit read Cs------')
    Schiff_res,Schiff_funs = ReadCs_Schiff(Schiff_file)
    return Get_Schiff_res(Schiff_res,Schiff_funs,swap_np=swap_np)


## file should have neutron then proton, swap_np=True will assume the file is flipped tho.
EDM_file = '/home/jackdra/LQCD/Results/HPCCRES_FixedG2/graphs/FF/AfterNotes/PRL_cont/d_mpi.py3p'
EDM_res = ReadAndFormat_EDM(EDM_file,swap_np=True)

Schiff_file = '/home/jackdra/LQCD/Results/HPCCRES_FixedG2/graphs/FF/AfterNotes/PRL_cont/S_mpi.py3p'
Schiff_res = ReadAndFormat_Schiff(Schiff_file,swap_np=False)

# GetLabels(prl_funs,prl_res,post=' Same C2 sign')


# fmt_prlres = FormatForLatex(prl_res)
both_res = pa.concat([EDM_res,Schiff_res],axis=1)

EDM_res.applymap(lambda x : x.MakeValAndErr(Dec=2))
Schiff_res.applymap(lambda x : x.MakeValAndErr(Dec=2))
both_res.applymap(lambda x : x.MakeValAndErr(Dec=2))
print(EDM_res.applymap(lambda x : '$'+x.MakeValAndErr(Dec=2,latex=True)+'$').to_latex(escape=False))
print(Schiff_res.applymap(lambda x : '$'+x.MakeValAndErr(Dec=2,latex=True)+'$').to_latex(escape=False))
print(both_res.applymap(lambda x : '$'+x.MakeValAndErr(Dec=2,latex=True)+'$').to_latex(escape=False))



def plot_EDM(this_series,square_mpi=False,fn=this_dir+'/TestGraphs/testCs.pdf',title=''):
    C1,C2,C3,chi2 = this_series.values
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
    pl.fill_between(xvals,y_C1_Up.values,y_C1_Down.values,color='blue',alpha=0.3,edgecolor='none')
    pl.fill_between(xvals,y_C2_Up.values,y_C2_Down.values,color='red',alpha=0.3,edgecolor='none')
    pl.fill_between(xvals,y_C12_Up.values,y_C12_Down.values,color='green',alpha=0.3,edgecolor='none')
    pl.title(title)
    pl.xlabel(r'$m_{\pi}^{2} [MeV]^{2}$')
    pl.ylabel(r'$d [e fm]$')
    pl.legend(loc='upper left')
    pl.grid()
    pl.savefig(fn)
    pl.clf()


def plot_Schiff(this_series,square_mpi=False,fn=this_dir+'/TestGraphs/testCs.pdf',title=''):
    C4,C5,chi2,chi2f = this_series.values
    xvals = pa.Series(np.arange(700))
    xplot = xvals
    if square_mpi:
        xplot = xplot**2
    C4_Avg = C4.apply(lambda x : x.GetAvg())
    C4_Up = C4.apply(lambda x : x.GetAvg()+x.GetStd())
    C4_Down = C4.apply(lambda x : x.GetAvg()-x.GetStd())
    pl.plot(xvals,C4_Avg.values,color='blue',label='$C_{4}$')
    pl.fill_between(xvals,C4_Up.values,C4_Down.values,color='blue',alpha=0.3)
    pl.legend()
    pl.savefig(fn)
    pl.clf()


plot_EDM(EDM_res.loc['proton',:],fn=this_dir+'/TestGraphs/testCs_EDM_proton.pdf',
         title=r'$d_{p}$ C parameters')
plot_EDM(EDM_res.loc['neutron',:],fn=this_dir+'/TestGraphs/testCs_EDM_neutron.pdf',
         title=r'$d_{n}$ C parameters')
# plot_FF(Schiff_pres,fn=this_dir+'/TestGraphs/testCs_Schiff_proton.pdf',title='Schiff Proton')
# plot_FF(Schiff_nres,fn=this_dir+'/TestGraphs/testCs_Schiff_neutron.pdf',title='Schiff Neutron')
