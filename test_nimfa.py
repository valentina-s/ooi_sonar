#!/usr/local/bin/python
'''
Testing NMF initialization and various methods using NIMFA
'''

import os, sys, glob, re
import datetime
from matplotlib.dates import date2num, num2date
import numpy as np
from calendar import monthrange
import h5py
from sklearn import decomposition
import nimfa

sys.path.insert(0,'/home/wu-jung/code_git/mi-instrument/')
from concat_raw import get_num_days_pings, get_data_from_h5
from echogram_decomp import find_nearest_time_idx,reshape_into_3freq,reshape_into_1freq,\
    sep_into_freq,plot_decomp_v,plot_decomp_transform

import matplotlib.pyplot as plt
from modest_image import imshow
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Set default colormap                               
plt.rcParams['image.cmap'] = 'jet'

# Get info of all files and set path                 
data_path = '/home/wu-jung/internal_2tb/ooi_sonar/figs/20170310_monthly_smpl'
save_path = '/home/wu-jung/internal_2tb/ooi_sonar/figs/20170310_nmf_options'
fname_form = '*.h5'
fname_all = glob.glob(os.path.join(data_path,fname_form))
n_comp = 5                # number of components for NMF

# Constant params
fmt = '%Y%m%d'
all_hr = range(24)         # list of all hour
all_min = range(1,11)      # list of all minutes
pings_per_day = len(all_hr)*len(all_min)  # number of pings per day
max_missing_ping = 5       # maximum number of missing pings in a day

# Load data
fname = 'CE04OSPS_20150901-20150930_smpl.h5'
f = h5py.File(os.path.join(data_path,fname),'r')
Sv_mtx0 = np.asarray(f['Sv_mtx'])
depth_tick = np.asarray(f['depth_tick'])
depth_label = np.asarray(f['depth_label'])
time_tick = np.asarray(f['time_tick'])
time_label = np.asarray(f['time_label'])
depth_bin_num = np.asarray(f['depth_bin_num'])
pings_per_day = np.asarray(f['pings_per_day'])
f.close()

fname = 'CE04OSPS_20151001-20151031_smpl.h5'
#fname = 'CE04OSPS_20160501-20160601_smpl.h5'
f = h5py.File(os.path.join(data_path,fname),'r')
Sv_mtx1 = np.asarray(f['Sv_mtx'])
f.close()

# Set plotting params
vec_len_each_day = pings_per_day*depth_bin_num  # length of vector for 1 day

# Prep data
if Sv_mtx0.shape[2]>Sv_mtx1.shape[2]:
    Sv_mtx0 = np.delete(Sv_mtx0,range(Sv_mtx1.shape[2],Sv_mtx0.shape[2]),axis=2)
else:
    Sv_mtx1 = np.delete(Sv_mtx1,range(Sv_mtx0.shape[2],Sv_mtx1.shape[2]),axis=2)

Sv_vec0 = reshape_into_3freq(Sv_mtx0,vec_len_each_day)
Sv_vec1 = reshape_into_3freq(Sv_mtx1,vec_len_each_day)

# Try out NIMFA
# Default nmf
nmf = nimfa.Nmf(Sv_vec0-Sv_vec0.min(),rank=5)
nmf_fit = nmf()
W = nmf_fit.basis()
H = nmf_fit.coef()
V = np.empty((5,1440,1046))
for Hcomp in range(H.shape[0]):
    V[Hcomp,:,:] = H[Hcomp,:].reshape((1440,1046))
    
# SNMF
snmf = nimfa.Snmf(Sv_vec0-Sv_vec0.min(),rank=5)
snmf_fit = snmf()
W_snmf = snmf_fit.basis()
H_snmf = snmf_fit.coef()
V_snmf = np.empty((5,1440,1046))
for c in range(H_snmf.shape[0]):
    V_snmf[c,:,:] = H_snmf[c,:].reshape((1440,1046))

snmf1 = nimfa.Snmf(Sv_vec1-Sv_vec1.min(),rank=5)
snmf1_fit = snmf1()
W_snmf1 = snmf1_fit.basis()
H_snmf1 = snmf1_fit.coef()
V_snmf1 = np.empty((5,1440,1046))
for c in range(H_snmf1.shape[0]):
    V_snmf1[c,:,:] = H_snmf1[c,:].reshape((1440,1046))

snmf2 = nimfa.Snmf(Sv_vec1-Sv_vec1.min(),rank=5,\
                   seed='fixed',W=np.copy(W_snmf),H=np.copy(H_snmf))
snmf2_fit = snmf2()
W_snmf2 = snmf2_fit.basis()
H_snmf2 = snmf2_fit.coef()
V_snmf2 = np.empty((5,1440,1046))
for c in range(H_snmf2.shape[0]):
    V_snmf2[c,:,:] = H_snmf2[c,:].reshape((1440,1046))

    
erank = snmf.estimate_rank(rank_range=[2,3,4,5,6,7,8],\
                           n_run=3, idx=0, what='all')
erank1 = snmf1.estimate_rank(rank_range=[2,3,4,5,6,7,8],\
                             n_run=3, idx=0, what='all')
erank2 = snmf2.estimate_rank(rank_range=[2,3,4,5,6,7,8],\
                             n_run=3, idx=0, what='all')
    
fig,ax=plt.subplots(5,2)
for c in range(H_snmf.shape[0]):
    ax[c,0].imshow(V_snmf[c,0:480,:].T,aspect='auto')
    ax[c,1].imshow(V_snmf1[c,0:480,:].T,aspect='auto')
fig.set_figwidth(8)
fig.set_figheight(6)

    
    
    
    
    
    
    
    
    
    
    
    
    
# Assemblge plot_param
plot_params = dict(zip(['depth_tick','depth_label','time_tick','time_label'],\
                       [depth_tick,depth_label,time_tick,time_label]))

fig_comp = plot_decomp_v(v_comps0,plot_params)
fig_comp.set_figheight(3)
fig_comp.set_figwidth(8)
fig_comp.suptitle('Set 0')

fig_comp = plot_decomp_v(v_comps1,plot_params)
fig_comp.set_figheight(3)
fig_comp.set_figwidth(8)
fig_comp.suptitle('Set 1 no init')

fig_comp = plot_decomp_v(v_comps2,plot_params)
fig_comp.set_figheight(3)
fig_comp.set_figwidth(8)
fig_comp.suptitle('Set 1 with init')

save_fname = 'test1.png'
comp_num = 1;
fig,ax=plt.subplots(1,3)
ax[0].imshow(v_comps0[comp_num,0,:,:],aspect='auto')
ax[0].set_title('Set 1')
ax[1].imshow(v_comps1[comp_num,0,:,:],aspect='auto')
ax[1].set_title('Set 2 no init')
ax[2].imshow(v_comps2[comp_num,0,:,:],aspect='auto')
ax[2].set_title('Set 2 with init')
fig.set_figwidth(12)
fig.set_figheight(2)
fig.savefig(os.path.join(save_path,save_fname),dpi=200)

    
save_fname = 'test1_new2.png'
comp_num = 1;
fig,ax=plt.subplots(1,2)
ax[0].imshow(v_comps0[comp_num,0,:,:],aspect='auto')
ax[0].set_title('Set 1')
ax[1].imshow(v_comps1[comp_num,0,:,:],aspect='auto')
ax[1].set_title('Set 2 no init')
fig.set_figwidth(8)
fig.set_figheight(2)
fig.savefig(os.path.join(save_path,save_fname),dpi=200)
