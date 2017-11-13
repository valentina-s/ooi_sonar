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
save_path = '/home/wu-jung/internal_2tb/ooi_sonar/figs/20170313_figs_incubator_final'
fname_form = '*.h5'
fname_all = glob.glob(os.path.join(data_path,fname_form))
n_comp = 5                # number of components for NMF

# Constant params
all_hr = range(24)         # list of all hour
all_min = range(1,11)      # list of all minutes
pings_per_day = len(all_hr)*len(all_min)  # number of pings per day
max_missing_ping = 5       # maximum number of missing pings in a day

# Load data
num_file = 11
for iF in range(num_file):
    fname = fname_all[iF+12]
    print 'Processing '+ fname
    f = h5py.File(os.path.join(data_path,fname),'r')
    if iF==0:
        Sv_mtx = np.asarray(f['Sv_mtx'])
        depth_tick = np.asarray(f['depth_tick'])
        depth_label = np.asarray(f['depth_label'])
        time_tick = np.asarray(f['time_tick'])
        time_label = np.asarray(f['time_label'])
        depth_bin_num = np.asarray(f['depth_bin_num'])
        pings_per_day = np.asarray(f['pings_per_day'])
        print 'number of days: '+str(Sv_mtx.shape[2]/pings_per_day)
    else:
        Sv_mtx = np.concatenate((Sv_mtx,np.asarray(f['Sv_mtx'])),axis=2)
        print 'number of days: '+str(Sv_mtx.shape[2]/pings_per_day)
    f.close()

bin_size = np.int_(depth_label[-1])/depth_tick[-1]
    
# Take out surface reflection
num_del_idx = 40
Sv_mtx = np.delete(Sv_mtx,range(num_del_idx),axis=1)
depth_bin_num = depth_bin_num-num_del_idx
depth_max = np.int_(depth_label[-1])
depth_min = bin_size*(num_del_idx+1)
    
# Set plotting params
vec_len_each_day = pings_per_day/2*depth_bin_num  # length of vector for 1 day

# Sliding window-- scikit-learn NMF
days_win = 30     # number of days per window
days_overlap = 5  # number of days overlapping when sliding
num_win = np.int_(np.floor((Sv_mtx.shape[2]/np.float(pings_per_day)-days_win)/days_overlap))
num_comp = 5


# Plot a daily echogram
for n_day in range(0,61):
    #n_day = 35  # which day in the whole sequence
    idx_n_day = range(pings_per_day)+n_day*pings_per_day
    fig,ax = plt.subplots(1)
    cax = ax.imshow(Sv_mtx[1,:,idx_n_day].T,aspect='auto',vmax=-30,vmin=-80)
    ax.set_yticks(depth_tick)
    ax.set_yticklabels(depth_label)
    ax.set_xticks(time_tick)
    ax.set_xticklabels(['00:00','06:00','12:00','18:00','24:00'])
    ax.set_xlabel('Time of day',fontsize=13)
    ax.set_ylabel('Depth (m)',fontsize=13)
    ax.set_title('Day %02d' % n_day,fontsize=13)
    fig.set_figwidth(4)
    fig.set_figheight(2.5)
    fig.colorbar(cax)
    plt.tight_layout()
    fig.savefig(os.path.join(save_path,'echogram_day%02d.png' % n_day),dpi=150)
    plt.close()

    
# Plot echogram for a few days
n_day = np.arange(35,40)  # which day in the whole sequence
idx_n_day = range(pings_per_day*n_day.shape)+n_day[0]*pings_per_day
fig,ax = plt.subplots(1)
cax = ax.imshow(Sv_mtx[1,:,idx_n_day].T,aspect='auto',vmax=-30,vmin=-80)
ax.set_yticks(depth_tick)
ax.set_yticklabels(depth_label)
ax.set_xticks(np.arange(0,pings_per_day*n_day.shape,pings_per_day))
ax.set_xticklabels([str(x) for x in n_day])
ax.set_xlabel('Days',fontsize=13)
ax.set_ylabel('Depth (m)',fontsize=13)
fig.set_figwidth(5)
fig.set_figheight(2.5)
fig.colorbar(cax)
plt.tight_layout()
fig.savefig(os.path.join(save_path,'day%02d-%02d_echogram.png' % (n_day[0],n_day[-1])),dpi=150)



# PCA, ICA, NMF Comparison =========================
# Get data
iS = 0
num_comp = 5
idx_slide = iS*pings_per_day*days_overlap+range(pings_per_day*days_win)
idx_slide = idx_slide[::2]
Sv_vec_tmp = reshape_into_3freq(Sv_mtx[:,:,idx_slide],vec_len_each_day)

# PCA ===========================
pca = decomposition.PCA(n_components=num_comp)
#W_pca = pca.fit_transform(Sv_vec_tmp-Sv_vec_tmp.min())
W_pca = pca.fit_transform(Sv_vec_tmp)
H_pca = pca.components_

# Plot each components
nan_insert = np.empty((20,depth_bin_num))
nan_insert[:] = np.nan
V_pca = np.empty((num_comp,pings_per_day/2*3,depth_bin_num))
for c in range(num_comp):
    V_pca[c,:,:] = H_pca[c,:].reshape((pings_per_day/2*3,depth_bin_num))

fig,ax=plt.subplots(1,num_comp,sharey=True)
for c in range(num_comp):
    cmean = np.mean(V_pca[c,:,:])
    cstd = np.std(V_pca[c,:,:])
    cmin = max((0,cmean-2*cstd))
    cmax = min((np.max(V_pca[c,:,:]),cmean+4*cstd))
    ax[c].imshow(np.concatenate((V_pca[c,0:240,:],nan_insert,\
                                 V_pca[c,240:480,:],nan_insert,\
                                 V_pca[c,480:-1,:])).T,\
                 aspect='auto',vmax=cmax,vmin=cmin)
    ax[c].set_xticks((120,360,600))
    ax[c].set_xticklabels(['38k','120k','200k'])
    ax[c].tick_params('both', length=0)
    ax[c].set_yticks((0,depth_bin_num))
    ax[c].set_yticklabels([str(np.int_(depth_min)),str(depth_max)])
    if c==0:
        ax[c].set_ylabel('Depth (m)')
fig.set_figwidth(16)
fig.set_figheight(1.7)
fig.suptitle('Window %d' % iS)
fig_fname = 'PCA_comp_win%02d.png' % iS
fig.savefig(os.path.join(save_path,fig_fname),dpi=150)
plt.close()


# ICA ======================
ica = decomposition.FastICA(n_components=num_comp)
#W_pca = pca.fit_transform(Sv_vec_tmp-Sv_vec_tmp.min())
W_ica = ica.fit_transform(Sv_vec_tmp)
H_ica = ica.components_

# Plot each components
nan_insert = np.empty((20,depth_bin_num))
nan_insert[:] = np.nan
V_ica = np.empty((num_comp,pings_per_day/2*3,depth_bin_num))
for c in range(num_comp):
    V_ica[c,:,:] = H_ica[c,:].reshape((pings_per_day/2*3,depth_bin_num))

fig,ax=plt.subplots(1,num_comp,sharey=True)
for c in range(num_comp):
    cmean = np.mean(V_ica[c,:,:])
    cstd = np.std(V_ica[c,:,:])
    cmin = max((0,cmean-2*cstd))
    cmax = min((np.max(V_ica[c,:,:]),cmean+4*cstd))
    ax[c].imshow(np.concatenate((V_ica[c,0:240,:],nan_insert,\
                                 V_ica[c,240:480,:],nan_insert,\
                                 V_ica[c,480:-1,:])).T,\
                 aspect='auto',vmax=cmax,vmin=cmin)
    ax[c].set_xticks((120,360,600))
    ax[c].set_xticklabels(['38k','120k','200k'])
    ax[c].tick_params('both', length=0)
    ax[c].set_yticks((0,depth_bin_num))
    ax[c].set_yticklabels([str(np.int_(depth_min)),str(depth_max)])
    if c==0:
        ax[c].set_ylabel('Depth (m)')
fig.set_figwidth(16)
fig.set_figheight(1.7)
fig.suptitle('Window %d' % iS)
fig_fname = 'ICA_comp_win%02d.png' % iS
fig.savefig(os.path.join(save_path,fig_fname),dpi=150)
plt.close()


# NMF ===============
nmf = decomposition.NMF(n_components=num_comp)
W_nmf = nmf.fit_transform(Sv_vec_tmp-Sv_vec_tmp.min())
H_nmf = nmf.components_

# Plot each components
nan_insert = np.empty((20,depth_bin_num))
nan_insert[:] = np.nan
V_nmf = np.empty((num_comp,pings_per_day/2*3,depth_bin_num))
for c in range(num_comp):
    V_nmf[c,:,:] = H_nmf[c,:].reshape((pings_per_day/2*3,depth_bin_num))

fig,ax=plt.subplots(1,num_comp,sharey=True)
for c in range(num_comp):
    cmean = np.mean(V_nmf[c,:,:])
    cstd = np.std(V_nmf[c,:,:])
    cmin = max((0,cmean-2*cstd))
    cmax = min((np.max(V_nmf[c,:,:]),cmean+4*cstd))
    ax[c].imshow(np.concatenate((V_nmf[c,0:240,:],nan_insert,\
                                 V_nmf[c,240:480,:],nan_insert,\
                                 V_nmf[c,480:-1,:])).T,\
                 aspect='auto',vmax=cmax,vmin=cmin)
    ax[c].set_xticks((120,360,600))
    ax[c].set_xticklabels(['38k','120k','200k'])
    ax[c].tick_params('both', length=0)
    ax[c].set_yticks((0,depth_bin_num))
    ax[c].set_yticklabels([str(np.int_(depth_min)),str(depth_max)])
    if c==0:
        ax[c].set_ylabel('Depth (m)')
fig.set_figwidth(16)
fig.set_figheight(1.7)
fig.suptitle('Window %d' % iS)
fig_fname = 'NMF_comp_win%02d.png' % iS
fig.savefig(os.path.join(save_path,fig_fname),dpi=150)
plt.close()




# ========== Sliding window: sci-kit learn NMF ==========
H_nmf_all = np.empty((num_win,num_comp,pings_per_day/2*3*depth_bin_num))
W_nmf_all = np.empty((num_win,days_win,num_comp))
for iS in range(num_win):  # slide through all windows
    print 'window '+str(iS)
    
    idx_slide = iS*pings_per_day*days_overlap+range(pings_per_day*days_win)
    idx_slide = idx_slide[::2]
    Sv_vec_tmp = reshape_into_3freq(Sv_mtx[:,:,idx_slide],vec_len_each_day)
    
    # Plot echogram of current window
    print 'Plotting window %d' % iS
    #cmean = np.mean(Sv_mtx[:,:,idx_slide])
    #cstd = np.std(Sv_mtx[:,:,idx_slide])
    #cmin = max((np.min(Sv_mtx[:,:,idx_slide]),cmean-2*cstd))
    #cmax = min((np.max(Sv_mtx[:,:,idx_slide]),cmean+4*cstd))
    fig,ax = plt.subplots(3,1,sharex=True)
    for ff in range(Sv_mtx.shape[0]):
        #ax[ff].imshow(Sv_mtx[ff,:,idx_slide].T,aspect='auto',vmax=cmax,vmin=cmin)
        im = imshow(ax[ff],Sv_mtx[ff,:,idx_slide].T,aspect='auto',vmax=-30,vmin=-80)
        #ax[ff].imshow(Sv_mtx[ff,:,idx_slide].T,aspect='auto',vmax=-30,vmin=-80)
        ax[ff].set_xticks(range(0,pings_per_day/2*days_win,pings_per_day/2))
        ax[ff].set_xticklabels([str(x) for x in range(days_win)])
        ax[ff].set_yticks((0,depth_bin_num))
        ax[ff].set_yticklabels([str(np.int_(depth_min)),str(depth_max)])
        ax[ff].set_ylabel('Depth (m)')
        if ff==2:
            ax[ff].set_xlabel('Days within sliding window')
        for dd in range(days_win):
            ax[ff].plot((dd*pings_per_day/2,dd*pings_per_day/2),(0,depth_bin_num-0.5),'k--',linewidth=0.5)
        divider = make_axes_locatable(ax[ff])
        cax = divider.append_axes("right", size="1%", pad=0.05)
        cbar = plt.colorbar(im,cax=cax)
        cbar.ax.tick_params(labelsize=14)
    fig.set_figheight(3)
    fig.set_figwidth(16)
    fig_fname = 'echogram_win%02d.png' % iS
    fig.savefig(os.path.join(save_path,fig_fname),dpi=150)
    plt.close()
    
    
    # Scikit learn NMF
    nmf_min = decomposition.NMF(n_components=n_comp)
    W_nmf_all[iS,:,:] = nmf_min.fit_transform(Sv_vec_tmp-Sv_vec_tmp.min())
    H_nmf_all[iS,:,:] = nmf_min.components_

# Plot each components
nan_insert = np.empty((20,depth_bin_num))
nan_insert[:] = np.nan
for iS in range(num_win):
    V_nmf = np.empty((num_comp,pings_per_day/2*3,depth_bin_num))
    for c in range(num_comp):
        V_nmf[c,:,:] = H_nmf_all[iS,c,:].reshape((pings_per_day/2*3,depth_bin_num))

    fig,ax=plt.subplots(1,5,sharey=True)
    for c in range(num_comp):
        cmean = np.mean(V_nmf[c,:,:])
        cstd = np.std(V_nmf[c,:,:])
        cmin = max((0,cmean-2*cstd))
        cmax = min((np.max(V_nmf[c,:,:]),cmean+4*cstd))
        ax[c].imshow(np.concatenate((V_nmf[c,0:240,:],nan_insert,\
                                     V_nmf[c,240:480,:],nan_insert,\
                                     V_nmf[c,480:-1,:])).T,\
                     aspect='auto',vmax=cmax,vmin=cmin)
        ax[c].set_xticks((120,360,600))
        ax[c].set_xticklabels(['38k','120k','200k'])
        ax[c].tick_params('both', length=0)
        ax[c].set_yticks((0,depth_bin_num))
        ax[c].set_yticklabels([str(np.int_(depth_min)),str(depth_max)])
        if c==0:
            ax[c].set_ylabel('Depth (m)')
    fig.set_figwidth(16)
    fig.set_figheight(1.7)
    fig.suptitle('Window %d' % iS)
    fig_fname = 'comp_win%02d.png' % iS
    fig.savefig(os.path.join(save_path,fig_fname),dpi=150)
    plt.close()

# Correlation matrix
R = np.corrcoef(H_nmf_all.reshape((-1,H_nmf_all.shape[-1])))
np.fill_diagonal(R, np.nan)

# R for part of the windows
num_win_incl = 20
fig,ax = plt.subplots(1)
im=ax.imshow(R[0:num_win_incl*num_comp,0:num_win_incl*num_comp],cmap=plt.get_cmap('Spectral_r'))
for iS in range(num_comp,num_comp*num_win_incl,num_comp):
    plt.plot((-0.5,num_win_incl*num_comp-0.5),(iS,iS),'k--',linewidth=0.5)
for iS in range(num_comp,num_comp*num_win_incl,num_comp):
    plt.plot((iS,iS),(-0.5,num_win_incl*num_comp-0.5),'k--',linewidth=0.5)
ax.set_xticks(np.arange(0,num_win_incl*num_comp,num_comp))
ax.set_xticklabels(np.arange(num_win_incl))
ax.set_yticks(np.arange(0,num_win_incl*num_comp,num_comp))
ax.set_yticklabels(np.arange(num_win_incl))
ax.set_xlabel('Window number',fontsize=16)
ax.set_ylabel('Window number',fontsize=16)
fig.set_figwidth(12)
fig.set_figheight(12)
fig.colorbar(im)
fig_fname = 'R_win0-19.png'
fig.savefig(os.path.join(save_path,fig_fname),dpi=150)
plt.close()

# R for all windows
fig,ax = plt.subplots(1)
im=ax.imshow(R,cmap=plt.get_cmap('Spectral_r'))
for iS in range(num_comp,num_comp*num_win,num_comp):
    plt.plot((-0.5,num_win*num_comp-0.5),(iS,iS),'k--',linewidth=0.5)
for iS in range(num_comp,num_comp*num_win,num_comp):
    plt.plot((iS,iS),(-0.5,num_win*num_comp-0.5),'k--',linewidth=0.5)
ax.set_xticks(np.arange(0,num_win*num_comp,num_comp))
ax.set_xticklabels(np.arange(num_win))
ax.set_yticks(np.arange(0,num_win*num_comp,num_comp))
ax.set_yticklabels(np.arange(num_win))
fig.set_figwidth(12)
fig.set_figheight(12)
fig.colorbar(im)
fig_fname = 'R_win_all.png'
fig.savefig(os.path.join(save_path,fig_fname),dpi=150)
plt.close()

# Save results
f = h5py.File(os.path.join(save_path,'11files_5comps.h5'))
f.create_dataset("Sv_mtx",data=Sv_mtx)
f.create_dataset("pings_per_day",data=pings_per_day)
f.create_dataset("depth_bin_num",data=depth_bin_num)
f.create_dataset("depth_tick",data=depth_tick)
f.create_dataset("depth_label",data=depth_label)
f.create_dataset("time_tick",data=time_tick)
f.create_dataset("time_label",data=time_label)
f.create_dataset("H_nmf_all",data=H_nmf_all)
f.create_dataset("W_nmf_all",data=W_nmf_all)
f.create_dataset("num_del_idx",data=num_del_idx)
f.create_dataset("vec_len_each_day",data=vec_len_each_day)
f.create_dataset("days_win",data=days_win)
f.create_dataset("days_overlap",data=days_overlap)
f.create_dataset("num_win",data=num_win)
f.create_dataset("num_comp",data=num_comp)
f.create_dataset("R",data=R)
f.close()



    
    
    
if 0:
    # ========== Sliding window: NIMFA SNMF ==========
    days_win = 30     # number of days per window
    days_overlap = 5  # number of days overlapping when sliding
    num_win = np.int_(np.floor((Sv_mtx.shape[2]/np.float(pings_per_day)-days_win)/days_overlap))
    num_comp = 5

    H_snmf_all = np.empty((num_win,num_comp,pings_per_day/2*3*depth_bin_num))
    for iS in range(num_win):  # slide through all windows
        print 'window '+str(iS)

        idx_slide = iS*pings_per_day*days_overlap+range(pings_per_day*days_win)
        idx_slide = idx_slide[::2]
        Sv_vec_tmp = reshape_into_3freq(Sv_mtx[:,:,idx_slide],vec_len_each_day)

        # NIMFA SNMF
        snmf = nimfa.Snmf(Sv_vec_tmp-Sv_vec_tmp.min(),rank=num_comp)
        snmf_fit = snmf()
        W_snmf = snmf_fit.basis()
        H_snmf = snmf_fit.coef()
        V_snmf = np.empty((num_comp,pings_per_day/2*3,depth_bin_num))
        for c in range(H_snmf.shape[0]):
            V_snmf[c,:,:] = H_snmf[c,:].reshape((pings_per_day/2*3,depth_bin_num))
        H_snmf_all[iS,:,:] = H_snmf

        ## Scikit learn NMF
        #nmf_min = decomposition.NMF(n_components=n_comp)
        #W_nmf_min = nmf_min.fit_transform(Sv_vec_tmp-Sv_vec_tmp.min())
        #H_nmf_min = nmf_min.components_
        #V_nmf_min = np.empty((num_comp,pings_per_day/2*3,depth_bin_num))
        #for c in range(H_nmf_min.shape[0]):
        #    V_nmf_min[c,:,:] = H_nmf_min[c,:].reshape((pings_per_day/2*3,depth_bin_num))
        #H_nmf_all[iS,:,:] = H_nmf_min

        # Plot NMF results
        if 0:
            nan_insert = np.empty((20,depth_bin_num))
            nan_insert[:] = np.nan
            fig,ax=plt.subplots(1,5,sharey=True)
            for c in range(H_nmf_min.shape[0]):
                cmean = np.mean(V_nmf_min[c,:,:])
                cstd = np.std(V_nmf_min[c,:,:])
                cmin = max((0,cmean-2*cstd))
                cmax = min((np.max(V_nmf_min[c,:,:]),cmean+4*cstd))
                ax[c].imshow(np.concatenate((V_nmf_min[c,0:240,:],nan_insert,\
                                             V_nmf_min[c,240:480,:],nan_insert,\
                                             V_nmf_min[c,480:-1,:])).T,\
                             aspect='auto',vmax=cmax,vmin=cmin)
                ax[c].set_xticks((120,360,600))
                ax[c].set_xticklabels(['38k','120k','200k'])
            fig.set_figwidth(10)
            fig.set_figheight(1.8)
            fig.suptitle('Window %d' % iS)

        # Plot SNMF results
        nan_insert = np.empty((20,depth_bin_num))
        nan_insert[:] = np.nan
        fig,ax=plt.subplots(1,5,sharey=True)
        for c in range(H_snmf.shape[0]):
            cmean = np.mean(V_snmf[c,:,:])
            cstd = np.std(V_snmf[c,:,:])
            cmin = max((0,cmean-2*cstd))
            cmax = min((np.max(V_snmf[c,:,:]),cmean+4*cstd))
            ax[c].imshow(np.concatenate((V_snmf[c,0:240,:],nan_insert,\
                                         V_snmf[c,240:480,:],nan_insert,\
                                         V_snmf[c,480:-1,:])).T,\
                         aspect='auto',vmax=cmax,vmin=cmin)
            ax[c].set_xticks((120,360,600))
            ax[c].set_xticklabels(['38k','120k','200k'])
        fig.set_figwidth(11)
        fig.set_figheight(1.8)
        fig.suptitle('Window %d' % iS)    

    # Correlation matrix
    R = np.corrcoef(H_snmf_all.reshape((-1,H_snmf_all.shape[-1])))
    np.fill_diagonal(R, np.nan)
    #R = np.tril(R)
    #R[R==0] = np.nan

    fig,ax = plt.subplots(1)
    im=ax.imshow(R,cmap=plt.get_cmap('Spectral_r'))#,vmax=1,vmin=-1)
    ax.set_xticks(np.arange(0,50,5))
    ax.set_xticklabels(np.arange(10)+1)
    ax.set_yticks(np.arange(0,50,5))
    ax.set_yticklabels(np.arange(10)+1)
    fig.set_figwidth(10)
    fig.set_figheight(10)
    fig.colorbar(im)




