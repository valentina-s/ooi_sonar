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
    '''
    Function to find nearest element
    time_wanted is a datetime object
    '''
    time_wanted_num = date2num(time_wanted)
    idx = np.searchsorted(all_timestamp_num, time_wanted_num, side="left")
    if idx > 0 and (idx == len(all_timestamp_num) or \
        math.fabs(time_wanted_num - all_timestamp_num[idx-1]) < math.fabs(time_wanted_num - all_timestamp_num[idx])):
        idx -= 1

    # If interval between the selected index and time wanted > 2 seconds
    sec_diff = datetime.timedelta(all_timestamp_num[idx]-time_wanted_num).total_seconds()
    if math.fabs(sec_diff)>2:
        return np.nan
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


def plot_decomp_v(v_comps,plot_params,freq='all'):
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
            elif comp==0:
                ax_v[comp].set_title('%s kHz' % freq, fontsize=16)
        fig_v.set_figwidth(5)
        fig_v.set_figheight(10)
    # 3 freq results
    if v_comps.ndim==4:
        fig_v,ax_v = plt.subplots(v_comps.shape[0],v_comps.shape[1],sharey=True,sharex=True)
        for comp in range(v_comps.shape[0]):  # component loop
            for ff in range(v_comps.shape[1]):  # frequency loop
                # Get color axis limtis
                cmean = np.mean(v_comps[comp,ff,:,:].reshape((-1,1)))
                cstd = np.std(v_comps[comp,ff,:,:].reshape((-1,1)))
                cmax = np.max(v_comps[comp,ff,:,:].reshape((-1,1)))
                # Plot
                im = imshow(ax_v[comp,ff],v_comps[comp,ff,:,:],aspect='auto',vmax=cmax,vmin=cmean-cstd*3)
                divider = make_axes_locatable(ax_v[comp,ff])
                cax = divider.append_axes("right", size="2%", pad=0.05)
                cbar = plt.colorbar(im,cax=cax)
                # Set labels
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
    return fig_v


def plot_decomp_transform(r_mtx,plot_params,freq='all'):
    ''' Plot transformed data '''
    fig,ax = plt.subplots(1)
    cax = ax.imshow(r_mtx.T,aspect='auto')
    ax.set_xlabel('Day of month')
    ax.set_ylabel('Components')
    fig.colorbar(cax)
    fig.set_figwidth(8)
    fig.set_figheight(2)
    ax.set_title(str(freq)+' kHz, %d components' % r_mtx.shape[1],fontsize=16)
    return fig
    
    
def get_data_based_on_day(every_ping,SITE_CODE,data_path):
    ''' Get data for pings specified using the datetime structure every_ping
        SITE_CODE is the location code for the instrument, right now it should be 'CE04OSPS' or 'CE02SHBP'
        data_path is the path where the data are stored
        The function will assemble the correct filename based on dates to get to the ping data
    '''
    cnt = 0  # counter for location in unpacked array
    y_want = np.unique(np.asarray([every_ping[x].year for x in range(len(every_ping))]))
    for y in y_want:  # loop through all years included
        ping_sel_y = [x for x in every_ping if x.year==y]  # all pings from year y               
        m_want = np.unique(np.asarray([ping_sel_y[x].month for x in range(len(ping_sel_y))]))

        for m in m_want:
            # y,m are the filename   
            ym_fname = datetime.datetime.strftime(datetime.date(y,m,1),'%Y%m');
            ping_sel_m = [x for x in ping_sel_y if x.month==m]  # all pings from month m

            # Open h5 file
            f = h5py.File(os.path.join(data_path,SITE_CODE+'_'+ym_fname+'.h5'),'r')

            # Get f['data_times'] idx for every hour in all days in the month
            all_idx = [find_nearest_time_idx(f['data_times'],t) for t in ping_sel_m]
            all_idx = np.array(all_idx)  # to allow numpy operation                                                                               
            # Inititalize mtx if this is the first file read
            if m==m_want[0] and y==y_want[0]:
                ping_time = np.empty(len(every_ping))  # pinging time                                                                            
                ping_time[:] = np.nan
                Sv_tmp = f['Sv'][:,:,0]                 # Sv array                                                                                
                Sv_mtx = np.empty((Sv_tmp.shape[0],Sv_tmp.shape[1],len(every_ping)))
                Sv_mtx[:] = np.nan
                bin_size = f['bin_size'][0]             # size of each depth bin
                
            # Fill in array
            notnanidx = np.ndarray.flatten(np.argwhere(~np.isnan(all_idx)))
            if notnanidx.any():  # if notnanidx not empty, then fill in values
                ping_time[notnanidx+cnt] = f['data_times'][all_idx[notnanidx].tolist()]
                Sv_mtx[:,:,notnanidx+cnt] = f['Sv'][:,:,all_idx[notnanidx].tolist()]

            # Increment counter
            cnt = cnt+len(all_idx)
            
            # Close h5 file
            f.close()
        
    return Sv_mtx,ping_time,bin_size


def clean_days(max_missing_ping,pings_per_day,Sv_mtx,ping_time):
    ''' Clean up Sv_mtx: delete days with too many missing pings '''
#     nan_ping_idx = np.isnan(Sv_mtx[0,0,:]).reshape((-1,pings_per_day))  # get ping index where Sv_mtx is NaN
    nan_ping_idx = np.isnan(ping_time).reshape((-1,pings_per_day))  # get ping index where Sv_mtx is NaN
    num_nan_of_day = np.sum(nan_ping_idx,1)  # number of missing pings of each day
    
    # Determine days and idx with lots of missing pings                                                             
    day_to_delete = []
    idx_to_delete = []
    for day in range(len(num_nan_of_day)):
        if num_nan_of_day[day]>max_missing_ping:
            day_to_delete.append(day)
            idx_to_delete.append(day*pings_per_day+np.arange(0,pings_per_day))
    day_to_delete = np.array(day_to_delete)
    idx_to_delete = np.ndarray.flatten(np.array(idx_to_delete))
    Sv_mtx_clean = np.delete(Sv_mtx,idx_to_delete,axis=2)  # delete bad idx
    ping_time_clean = np.delete(ping_time,idx_to_delete)  # delete bad idx
    return Sv_mtx_clean,ping_time_clean,idx_to_delete,day_to_delete


def clean_pings(Sv_mtx,ping_time):
    ''' Clean up Sv_mtx: fill in missing pings '''
#     nan_ping_idx = np.ndarray.flatten(np.argwhere(np.isnan(Sv_mtx[0,0,:])))  # ping index where Sv_mtx is NaN
    nan_ping_idx = np.ndarray.flatten(np.argwhere(np.isnan(ping_time)))  # ping index where Sv_mtx is NaN
    fill_idx = np.copy(nan_ping_idx)
    for idx in nan_ping_idx:
        if idx!=0:
            Sv_mtx[:,:,idx] = Sv_mtx[:,:,idx-1]
            ping_time[idx] = ping_time[idx-1]
        else:
            good_idx = np.ndarray.flatten(np.argwhere(~np.isnan(ping_time)))  # ping index where Sv_mtx is not NaN
            Sv_mtx[:,:,idx] = Sv_mtx[:,:,good_idx[0]]
            ping_time[idx] = ping_time[good_idx[0]]
    return Sv_mtx, ping_time
