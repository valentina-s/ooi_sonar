#!/usr/local/bin/python

'''
Echogram decomposition
'''
import os, sys, re
import h5py
import datetime
import math
import numpy as np
from matplotlib.dates import date2num, num2date
from calendar import monthrange
import matplotlib.pyplot as plt
from modest_image import imshow
from sklearn import decomposition
from mpl_toolkits.axes_grid1 import make_axes_locatable


def find_nearest_time_idx(all_timestamp_num,time_wanted):
    ''' Function to find nearest element '''
    time_wanted_num = date2num(time_wanted)
    idx = np.searchsorted(all_timestamp_num, time_wanted_num, side="left")
    if idx > 0 and (idx == len(all_timestamp_num) or \
        math.fabs(time_wanted_num - all_timestamp_num[idx-1]) < math.fabs(time_wanted_num - all_timestamp_num[idx])):
        return idx-1
    else:
        return idx


def reshape_into_3freq(mtx,vec_len):
    ''' Manipulate 3 freq Sv data into 1 long vector for each day '''
    long_vec_3freq = [np.reshape(np.transpose(mtx[band,:,:]),(-1,vec_len)) for band in range(3)]
    # convert from list to numpy array, dim: (3, 30 or 31, vec_len_each_day))
    long_vec_3freq = np.asarray(long_vec_3freq)
    # Swal dim 0 and 1 and reshape triple long vector for each day
    return np.reshape(np.swapaxes(long_vec_3freq,0,1),(-1,3*vec_len))


def reshape_into_1freq(mtx,vec_len,fidx):
    ''' Manipulate 1 freq Sv data into 1 long vector for each day '''
    return np.reshape(np.transpose(mtx[fidx,:,:]),(-1,vec_len))


def sep_into_freq(v_comps,pings_per_day,depth_bin_num):
    ''' Manipulate decomposed output (each a long vector of each day) into mtx
        organized into dim: (components,freq,depth,time) --> 3 freq data
        or             dim: (components,depth,time) --> 1 freq data
    '''
    # 3 freq data
    if v_comps.shape[1] == pings_per_day*3*depth_bin_num:
        comp_sep = np.reshape(v_comps,(-1,pings_per_day*3,depth_bin_num))      # dim: (6, 720, 1046)
        comp_sep_freq = np.reshape(np.swapaxes(comp_sep,1,2),(comp_sep.shape[0],depth_bin_num,-1,pings_per_day))  # dim: (6, 1046, 3, 240)
        return np.swapaxes(comp_sep_freq,1,2)  #  dim: (6, 3, 1046, 240)
    # 1 freq data
    elif v_comps.shape[1] == pings_per_day*depth_bin_num:
        comp_sep = np.reshape(v_comps,(-1,pings_per_day,depth_bin_num))  # dim: (6, 240, 1046)
        return np.swapaxes(comp_sep,1,2)  # dim: (6, 1046, 240)
    else:
        print 'Input invalid!'
        return []


def plot_decomp_v(v_comps,save_path,plot_params):
    ''' Plot each components of the decomposition '''
    # 1 freq results
    if v_comps.ndim==3:
        fig_v,ax_v = plt.subplots(v_comps.shape[0],sharey=True,sharex=True)
        for comp in range(v_comps.shape[0]):
            im = imshow(ax_v[comp],v_comps[comp,:,:],aspect='auto')
            divider = make_axes_locatable(ax_v[comp])
            cax = divider.append_axes("right", size="2%", pad=0.05)
            cbar = plt.colorbar(im,cax=cax)
            ax_v[comp].set_xticks(plot_params['time_tick'])
            ax_v[comp].set_xticklabels(plot_params['time_label'])
            ax_v[comp].set_yticks(plot_params['depth_tick'])
            ax_v[comp].set_yticklabels(plot_params['depth_label'])
            ax_v[comp].set_ylabel('Depth',fontsize=14)
            if comp==v_comps.shape[0]-1:
                ax_v[comp].set_xlabel('Time of day',fontsize=14)
        fig_v.set_figwidth(5)
        fig_v.set_figheight(10)
    # 3 freq results
    if v_comps.ndim==4:
        fig_v,ax_v = plt.subplots(v_comps.shape[0],v_comps.shape[1],sharey=True,sharex=True)
        for comp in range(v_comps.shape[0]):  # component loop
            for ff in range(v_comps.shape[1]):  # frequency loop
                im = imshow(ax_v[comp,ff],v_comps[comp,ff,:,:],aspect='auto')
                divider = make_axes_locatable(ax_v[comp,ff])
                cax = divider.append_axes("right", size="2%", pad=0.05)
                cbar = plt.colorbar(im,cax=cax)
                ax_v[comp,ff].set_xticks(plot_params['time_tick'])
                ax_v[comp,ff].set_xticklabels(plot_params['time_label'])
                ax_v[comp,ff].set_yticks(plot_params['depth_tick'])
                ax_v[comp,ff].set_yticklabels(plot_params['depth_label'])
                ax_v[comp,ff].set_ylabel('Depth',fontsize=14)
                if comp==v_comps.shape[0]-1:
                    ax_v[comp,ff].set_xlabel('Time of day',fontsize=14)
                elif comp==0:
                    if ff==0:
                        ax_v[comp,ff].set_title('38 kHz',fontsize=16)
                    elif ff==1:
                        ax_v[comp,ff].set_title('120 kHz',fontsize=16)
                    elif ff==2:
                        ax_v[comp,ff].set_title('200 kHz',fontsize=16)
                if ff==0:
                    ax_v[comp,ff].set_ylabel('Depth',fontsize=14)
        fig_v.set_figwidth(16)
        fig_v.set_figheight(10)
    save_fname = '%s%02d_%dcomps_vec.png' % (plot_params['year'],plot_params['month'],v_comps.shape[0])
    plt.savefig(os.path.join(save_path,save_fname))


def plot_decomp_transform(r_mtx,save_path,plot_params,freq='all'):
    ''' Plot transformed data '''
    plt.imshow(r_mtx,aspect='auto')
    plt.xlabel('Components')
    plt.ylabel('Day of month')
    plt.colorbar()
    plt.title(str(freq)+' kHz, %d components' % r_mtx.shape[1],fontsize=16)
    save_fname = '%d_%02d_%dcomps_transformed.png' % (plot_params['year'],plot_params['month'],r_mtx.shape[1])
    plt.savefig(os.path.join(save_path,save_fname))
