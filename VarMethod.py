#!/usr/bin/env python

import numpy as np
import xarray as xr
import scipy.linalg as lalg
from warnings import warn
import PlotData as jpl
import pandas as pa
from MiscFuns import map_str
from BootStrapping import BootStrap


"""
How to use:

takes a matrix of objects and creates performs the variational method on it


"""



class VariationalMethod(object):

    '''

    '''

    def __init__(self, matrix_corr = None , t0=2,dt=2,name='',symetrize=True):
        self.symetrize = symetrize
        self.name = name
        self.t0,self.dt = t0,dt
        self.left_evec = self.right_evec = None,None
        self.state_corrs = None
        if matrix_corr is not None:
            self.ImportMatrixCorr(matrix_corr)

    def ImportMatrixCorr(self,this_mc):
        ## matrix_corr is [..., ism, jsm, t]
        if isinstance(this_mc,(list,np.ndarray)):
            this_shape = np.array(this_mc).shape
            key_list = []
            if len(this_shape) > 3:
                for i in range(len(this_shape)-3):
                    key_list.append('extra_dim_'+str(i+1))
            key_list += ['ism','jsm','t_sep']
            self.matrix_corr = Make_DataArray(this_mc,key_list,name=self.name)
        elif isinstance(this_mc,xr.DataArray):
            self.matrix_corr = this_mc
        elif isinstance(this_mc,xr.Dataset):
            self.matrix_corr = this_mc.to_array()
        else:
            raise ImportError('datatype for this_mc not recognised '+str(type(this_mc)))
        if self.matrix_corr.shape[-3] != self.matrix_corr.shape[-2]:
            raise NotImplementedError('non-square matricies not implemented for variational method')
        self.ndim_other = self.matrix_corr.dims
        self.other_coords = self.matrix_corr.coords
        self.has_extra_dims = len(self.ndim_other) > 3
        this_dict = {}
        for idim in self.ndim_other[:-3]:
            this_dict[idim] = self.other_coords[idim].values
        self.other_coords = this_dict
        if self.has_extra_dims:
            self.matrix_corr = self.matrix_corr.stack(higher_dims=self.matrix_corr.dims[:-3])
        elif self.matrix_corr.ndim < 3:
            raise IOError('matrix of correlators must be atleast 3 dims for [ism,jsm,t]')
        if self.matrix_corr.dtype == object:
            self.booted_mcorr = self.matrix_corr
            self.nboot = this_mc.flatten()[0].nboot
            self.matrix_corr = xr.apply_ufunc(lambda x : x.Avg,self.matrix_corr,vectorize=True)
        else:
            self.booted_mcorr = None

        self.nism = self.matrix_corr.shape[-3]
        self.njsm = self.matrix_corr.shape[-2]
        self.nt = self.matrix_corr.shape[-1]


    def VectorizedGEVP(self,Amat,Bmat):
        if self.symetrize:
            Amat = (Amat.T + Amat)/2
            Bmat = (Bmat.T + Bmat)/2
            eval,right_evec  = lalg.eigh(Amat,b=Bmat)
            # Cmat = lalg.inv(Bmat)*Amat
            # eval,right_evec  = lalg.eigh(Cmat)
            left_evec = right_evec
        else:
            eval,left_evec,right_evec  = lalg.eig(Amat,b=Bmat,left=True,right=True)
        eval,left_evec,right_evec = self.SortEs(eval,left_evec.swapaxes(0,1),right_evec.swapaxes(0,1))
        self.eval.append(np.array(eval))
        self.eEnergy.append(-np.log(np.array(eval))/self.dt)
        self.left_evec.append(np.array(left_evec).swapaxes(0,1))
        self.right_evec.append(np.array(right_evec).swapaxes(0,1))

    def SortEs(self,el,lev,rev):
        return list(zip(*[[iel,ilev,irev] for _,iel,ilev,irev in sorted(zip(-np.log(el),el,lev,rev))]))


    def PerformVarMeth(self,t0=None,dt=None):
        if t0 is not None:
            self.t0 = t0
        if dt is not None:
            self.dt = dt
        these_coords = self.matrix_corr.coords['t_sep'].values
        if isinstance(these_coords[0],tuple):
            these_coords = np.array([ival[0] for ival in these_coords])
        if self.t0 not in these_coords:
            raise IOError('t0='+str(self.t0) + ' not in correlation matrix')
        if self.dt not in these_coords:
            raise IOError('dt='+str(self.dt) + ' not in correlation matrix')
        if self.t0+self.dt not in these_coords:
            raise IOError('t0+dt='+str(self.t0+self.dt) + ' not in correlation matrix')
        Amat = self.matrix_corr.sel(t_sep=self.t0+self.dt)
        Bmat = self.matrix_corr.sel(t_sep=self.t0)
        self.eval = []
        self.eEnergy = []
        self.left_evec = []
        self.right_evec = []
        eval_keys = {}
        left_evec_keys = {}
        right_evec_keys = {}
        if self.has_extra_dims:
            xr.apply_ufunc(self.VectorizedGEVP,Amat,Bmat,
                           input_core_dims=[['ism','jsm'],['ism','jsm']],vectorize=True)
            eval_keys.update(self.other_coords)
            left_evec_keys.update(self.other_coords)
            right_evec_keys.update(self.other_coords)
        else:
            self.VectorizedGEVP(Amat,Bmat)
            self.eval = self.eval[0]
            self.eEnergy = self.eEnergy[0]
            self.left_evec = self.left_evec[0]
            self.right_evec = self.right_evec[0]
        eval_keys['state'] = range(np.array(self.eval).shape[-1])
        left_evec_keys['ism'] = range(np.array(self.left_evec).shape[-1])
        left_evec_keys['state'] = range(np.array(self.left_evec).shape[-2])
        right_evec_keys['jsm'] = range(np.array(self.right_evec).shape[-1])
        right_evec_keys['state'] = range(np.array(self.right_evec).shape[-2])
        self.eval = Make_DataArray_Keys(self.eval,eval_keys,name='evals')
        self.eEnergy = Make_DataArray_Keys(self.eEnergy,eval_keys,name='eEnergy')
        self.left_evec = Make_DataArray_Keys(self.left_evec,left_evec_keys,name='left_evecs')
        self.right_evec = Make_DataArray_Keys(self.right_evec,right_evec_keys,name='right_evacs')
        if self.has_extra_dims:
            self.eval= self.eval.stack(higher_dims=self.eval.dims[:-1])
            self.eEnergy= self.eEnergy.stack(higher_dims=self.eEnergy.dims[:-1])
            self.left_evec= self.left_evec.stack(higher_dims=self.left_evec.dims[:-2])
            self.right_evec= self.right_evec.stack(higher_dims=self.right_evec.dims[:-2])
        self.ProjectCorrs()


    def ProjectCorrs(self):
        if not hasattr(self,'left_evec') or not hasattr(self,'right_evec'):
            warn('PerformVarMeth need to be called to project, attempting now with default t0='+str(self.t0)+' and dt='+str(self.dt))
            self.PerformVarMeth()
        self.proj_corr = xr.dot(self.left_evec,self.matrix_corr,dims='ism')
        self.proj_corr = xr.dot(self.right_evec,self.proj_corr,dims='jsm')
        if hasattr(self,'booted_mcorr') and self.booted_mcorr is not None:
            self.booted_proj_corr = self.proj_corr.copy()
            self.booted_proj_corr = xr.apply_ufunc(lambda x : BootStrap(thisnboot=self.nboot),self.booted_proj_corr,vectorize=True)
            for iboot in range(self.nboot):
                this_boot_mcorr = xr.apply_ufunc(lambda x : x.bootvals[iboot],self.booted_mcorr,vectorize=True)
                this_booted_proj_corr = xr.dot(self.left_evec,this_boot_mcorr,dims='ism')
                this_booted_proj_corr = xr.dot(self.right_evec,self.proj_corr,dims='jsm')
                def this_fun(itemp,boot_val):
                    itemp.bootvals = itemp.bootvals.append(pa.Series([boot_val]))
                    return itemp
                self.booted_proj_corr = xr.apply_ufunc(this_fun,self.booted_proj_corr,this_booted_proj_corr,vectorize=True)



    def MakeEffMass(self):
        self.eff_mass = np.log(self.matrix_corr / self.matrix_corr.roll(t_sep=-1,roll_coords=False))
        if not hasattr(self,'proj_corr'):
            self.ProjectCorrs()
        self.proj_eff_mass = np.log(self.proj_corr / self.proj_corr.roll(t_sep=-1,roll_coords=False))
        if hasattr(self,'booted_mcorr')  and self.booted_mcorr is not None:
            def boot_log_fun(num,denom):
                this_out = (num/denom).Log()
                return this_out
            self.booted_eff_mass = xr.apply_ufunc(boot_log_fun,self.booted_mcorr,self.booted_mcorr.roll(t_sep=-1,roll_coords=False),vectorize=True)
            self.booted_proj_eff_mass = xr.apply_ufunc(boot_log_fun,self.booted_proj_corr,self.booted_proj_corr.roll(t_sep=-1,roll_coords=False),vectorize=True)

    def PlotEffMass_Boot(self,plot_plane,**kwargs):
        if not hasattr(self,'booted_eff_mass') or self.booted_eff_mass is None:
            print('effective mass not computed, computing now')
            self.MakeEffMass()
        hold_series = jpl.null_series
        if 'color' in kwargs.keys():
            hold_series['color'] = kwargs['color']
            del kwargs['color']
        if 'shift' in kwargs.keys():
            hold_series['shift'] = kwargs['shift']
            del kwargs['shift']
        try:
            data_plot = self.booted_eff_mass.unstack().sel(**kwargs)
        except:
            data_plot = self.booted_eff_mass.sel(**kwargs).unstack()
        data_plot = data_plot.stack(all=data_plot.dims).to_pandas()
        this_index = pa.MultiIndex.from_tuples(list(map_str(data_plot.index)),names=data_plot.index.names)
        data_plot = pa.Series(data_plot.values,index =this_index)
        data_plot.apply(lambda x : x.Stats())
        data_avg = data_plot.apply(lambda x : x.Avg)
        data_std = data_plot.apply(lambda x : x.Std)
        hold_series['x_data'] = 'from_keys'
        hold_series['type'] = 'error_bar'
        if isinstance(data_plot.index[0],(tuple,list,np.ndarray)):
            if len(data_plot.index[0])>1:
                hold_series['key_select'] = data_plot.index[0]
                hold_series['key_select'] = list(hold_series['key_select'])[:-1]+[slice(None)]
                hold_series['type'] += '_vary'
        hold_series['y_data'] = data_avg
        hold_series['yerr_data'] = data_std
        hold_series['label'] = self.name
        plot_plane.AppendData(hold_series)
        return plot_plane

    def PlotEffMass_VarMeth_Boot(self,plot_plane,**kwargs):
        if not hasattr(self,'booted_proj_eff_mass')  or self.booted_proj_eff_mass is None:
            print('effective mass not computed, computing now')
            self.MakeEffMass()
        hold_series = jpl.null_series
        if 'color' in kwargs.keys():
            hold_series['color'] = kwargs['color']
            del kwargs['color']
        if 'shift' in kwargs.keys():
            hold_series['shift'] = kwargs['shift']
            del kwargs['shift']
        try:
            data_plot = self.booted_proj_eff_mass.unstack().sel(**kwargs)
        except:
            data_plot = self.booted_proj_eff_mass.sel(**kwargs).unstack()
        data_plot = data_plot.stack(all=self.booted_proj_eff_mass.dims).to_pandas()
        this_index = pa.MultiIndex.from_tuples(list(map_str(data_plot.index)),names=data_plot.index.names)
        data_plot = pa.Series(data_plot.values,index =this_index)
        data_plot.apply(lambda x : x.Stats())
        data_avg = data_plot.apply(lambda x : x.Avg)
        data_std = data_plot.apply(lambda x : x.Std)
        hold_series['x_data'] = 'from_keys'
        hold_series['key_select'] = data_plot.index[0]
        hold_series['type'] = 'error_bar'
        if isinstance(data_plot.index[0],(tuple,list,np.ndarray)):
            if len(data_plot.index[0])>1:
                hold_series['key_select'] = data_plot.index[0]
                hold_series['key_select'] = list(hold_series['key_select'])[:-1]+[slice(None)]
                hold_series['type'] += '_vary'
        hold_series['y_data'] = data_avg
        hold_series['yerr_data'] = data_std
        hold_series['label'] = self.name+' VarMeth'
        plot_plane.AppendData(hold_series)
        return plot_plane

    def PlotEffMass(self,plot_plane,**kwargs):
        if hasattr(self,'booted_mcorr')  and self.booted_mcorr is not None:
            return self.PlotEffMass_Boot(plot_plane,**kwargs)
        if not hasattr(self,'eff_mass'):
            print('effective mass not computed, computing now')
            self.MakeEffMass()
        hold_series = jpl.null_series
        if 'color' in kwargs.keys():
            hold_series['color'] = kwargs['color']
            del kwargs['color']
        if 'shift' in kwargs.keys():
            hold_series['shift'] = kwargs['shift']
            del kwargs['shift']
        try:
            data_plot = self.eff_mass.unstack().sel(**kwargs)
        except:
            data_plot = self.eff_mass.sel(**kwargs).unstack()
        data_plot = data_plot.stack(all=data_plot.dims).to_pandas()
        this_index = pa.MultiIndex.from_tuples(list(map_str(data_plot.index)),names=data_plot.index.names)
        data_plot = pa.Series(data_plot.values,index =this_index)
        hold_series['x_data'] = 'from_keys'
        hold_series['type'] = 'plot'
        if isinstance(data_plot.index[0],(tuple,list,np.ndarray)):
            if len(data_plot.index[0])>1:
                hold_series['key_select'] = data_plot.index[0]
                hold_series['key_select'] = list(hold_series['key_select'])[:-1]+[slice(None)]
                hold_series['type'] += '_vary'
        hold_series['y_data'] = data_plot
        hold_series['label'] = self.name
        plot_plane.AppendData(hold_series)
        return plot_plane

    def PlotEffMass_VarMeth(self,plot_plane,**kwargs):
        if hasattr(self,'booted_proj_corr')  and self.booted_proj_corr is not None:
            return self.PlotEffMass_VarMeth_Boot(plot_plane,**kwargs)
        if not hasattr(self,'eff_mass'):
            print('effective mass not computed, computing now')
            self.PlotEffMass_VarMeth()
        hold_series = jpl.null_series
        if 'color' in kwargs.keys():
            hold_series['color'] = kwargs['color']
            del kwargs['color']
        if 'shift' in kwargs.keys():
            hold_series['shift'] = kwargs['shift']
            del kwargs['shift']
        try:
            data_plot = self.proj_eff_mass.unstack().sel(**kwargs)
        except:
            data_plot = self.proj_eff_mass.sel(**kwargs).unstack()
        data_plot = data_plot.stack(all=self.proj_eff_mass.dims).to_pandas()
        this_index = pa.MultiIndex.from_tuples(list(map_str(data_plot.index)),names=data_plot.index.names)
        data_plot = pa.Series(data_plot.values,index =this_index)
        hold_series['x_data'] = 'from_keys'
        hold_series['key_select'] = data_plot.index[0]
        hold_series['type'] = 'plot'
        if isinstance(data_plot.index[0],(tuple,list,np.ndarray)):
            if len(data_plot.index[0])>1:
                hold_series['key_select'] = data_plot.index[0]
                hold_series['key_select'] = list(hold_series['key_select'])[:-1]+[slice(None)]
                hold_series['type'] += '_vary'
        hold_series['y_data'] = data_plot
        hold_series['label'] = self.name+' VarMeth'
        plot_plane.AppendData(hold_series)
        return plot_plane


def Make_DataArray(this_data,key_list,name=''):
    this_data = np.array(this_data)
    this_shape = this_data.shape
    this_coords = {}
    for ic,ikey in enumerate(key_list):
        this_coords[ikey] = range(this_shape[ic])
    return xr.DataArray(this_data,coords=this_coords,dims=list(this_coords.keys()),name=name)

def Make_DataArray_Keys(this_data,key_list,name=''):
    this_data = np.array(this_data)
    this_coords = {}
    for ikey,ival in key_list.items():
        this_coords[ikey] = []
        for jval in ival:
            this_coords[ikey].append(jval)
    if this_data.shape[0] != len(list(this_coords.values())[0]):
        this_data = this_data.reshape(*map(len,this_coords.values()))
    return xr.DataArray(this_data,coords=this_coords,dims=list(this_coords.keys()),name=name)

def CreateTestData():
    energy_states = np.array([0.2,0.4,0.8,1.6])
    coeff_list = np.array([[1,0.5,0.25,0.125],
                           [0.5,0.5,0.25,0.25],
                           [0.1,0.2,0.4,0.8],
                           [0.1,1,10,100]]).T

    coeff_err_per = 0
    energy_err_per = 2
    this_nboot = 20
    this_ncfg = 20
    nt = 50
    booted = True
    coeff_list = np.array([np.outer(ival,ival) for ival in coeff_list])
    coeff_list[1]
    coeff_err = coeff_list*coeff_err_per
    coeff_err = np.random.normal(loc=0,scale=coeff_err_per,size=list(coeff_list.shape)+[this_ncfg])
    # energy_states_err = np.random.normal(loc=1,scale=energy_err_per,size=list(energy_states.shape)+[this_ncfg])
    energy_states_err = np.random.normal(loc=0,scale=energy_err_per,size=list(energy_states.shape)+[this_ncfg])
    def exp_fun(val,ism,jsm,size):
        rand_coeff = np.array([iblist+icoeff for iblist,icoeff in zip(coeff_list[:,ism,jsm],coeff_err[:,ism,jsm,:])])
        rand_err = np.array([np.exp(-iblist*val)*np.exp(-ierr) for iblist,ierr in zip(energy_states,energy_states_err)])
        # rand_err_p1 = np.array([ierr + iblist*(val+1) for iblist,ierr in zip(energy_states,energy_states_err)])
        #     numer = rand_coeff[0]* np.exp(-rand_err[0])
        #     denom = rand_coeff[0]* np.exp(-rand_err_p1[0])
        #     print(np.log(np.mean(numer)/
        #                  np.mean(denom)),np.mean(rand_err[0]))
        return np.sum( rand_coeff* rand_err,axis=0)
    # data_shape = [coeff_list.shape[1],coeff_list.shape[2],6]
    data_shape = [coeff_list.shape[1],coeff_list.shape[2],nt]
    values = []
    if len(data_shape)>3:
        loop_shape = [np.prod(data_shape[:-3])]+data_shape[-3:]
        for iouter in range(loop_shape[0]):
            for ism in range(loop_shape[1]):
                for jsm in range(loop_shape[2]):
                    for it in range(loop_shape[3]):
                        if booted:
                            this_vals = exp_fun(it,ism,jsm,this_ncfg)
                            values.append(BootStrap(thisnboot = this_nboot,cfgvals=this_vals))
                        else:
                            values.append(exp_fun(it,ism,jsm,1))
                        ## not testing for bootstrapping yet
    else:
        for ism in range(data_shape[0]):
            for jsm in range(data_shape[1]):
                for it in range(data_shape[2]):
                    if booted:
                        this_vals = exp_fun(it,ism,jsm,this_nboot)
                        values.append(BootStrap(thisnboot = this_nboot,cfgvals=this_vals))
                    else:
                        values.append(exp_fun(it,ism,jsm,1))
                    # values.append(exp_fun(it,ism,jsm))
    return np.array(values).reshape(data_shape)

if __name__ == '__main__':
    raw_data = CreateTestData()
    test_class = VariationalMethod(matrix_corr=raw_data,symetrize=False)
    test_class.PerformVarMeth(t0=1,dt=5)
    # test_class.booted_mcorr
    # test_class.eEnergy
    # test_class.left_evec
    this_info = pa.Series()
    this_info['save_file'] = '/home/jackdra/LQCD/Scripts/Python_Analysis/TestGraphs/TestVar.pdf'
    this_info['title'] = 'Var Meth'
    this_info['x_label'] = 't'
    this_info['y_label'] = 'EffMass'
    data_plot = jpl.Plotting(plot_info=this_info)
    data_plot = test_class.PlotEffMass(data_plot)
    data_plot = test_class.PlotEffMass_VarMeth(data_plot)
    data_plot.PlotAll()
