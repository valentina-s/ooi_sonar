#!/usr/local/bin/python
''' NMF decomposition for each month '''

import os, sys, glob, re
import datetime
from matplotlib.dates import date2num, num2date
import numpy as np
from calendar import monthrange
import h5py
from sklearn import decomposition

sys.path.insert(0,'/home/wu-jung/code_git/mi-instrument/')
from concat_raw import get_num_days_pings, get_data_from_h5
from echogram_decomp import find_nearest_time_idx,reshape_into_3freq,reshape_into_1freq,\
    sep_into_freq,plot_decomp_v,plot_decomp_transform,\
    get_data_based_on_day,clean_days,clean_pings

import matplotlib.pyplot as plt
from modest_image import imshow
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Set default colormap                     
plt.rcParams['image.cmap'] = 'jet'

# Set path
data_path = '/home/wu-jung/whoi_bk/OOI/ooi_all_data_h5'
fname_form = '*.h5'
fname_all = glob.glob(os.path.join(data_path,fname_form))
save_path = '/home/wu-jung/internal_2tb/ooi_sonar/figs/20170308_nmf_monthly'


# Set params
SITE_CODE = 'CE04OSPS'
#SITE_CODE = 'CE02SHBP'
# Constant params
all_hr = range(24)         # list of all hour: 0-23
all_min = range(1,11)      # list of all minutes
pings_per_day = len(all_hr)*len(all_min)  # number of pings per day
max_missing_ping = 5       # maximum number of missing pings in a day
n_comp = 5                 # number of NMF components

# Run through all months
for yy in (2015,2016):
    if yy==2015:
        mm_range = range(8,13,1)
    elif yy==2016:
        mm_range = range(1,8,1)
        
    for mm in mm_range:
        print 'Processing data from site %s %04d%02d' % (SITE_CODE,yy,mm)
        
        _,daynum = monthrange(yy,mm)  # number of days in the month

        # Get datetime object for on the hour every hour in all days in the month
        all_day = range(1,daynum+1)  # list of all days
        every_ping = [datetime.datetime(yy,mm,x_day,x_hr,x_min,0) \
                      for x_day in all_day for x_hr in all_hr for x_min in all_min]
        
        # correct for time shift after Nov 2016
        curr_ym = '%04d%02d' % (yy,mm)
        if curr_ym>='201511' and SITE_CODE=='CE04OSPS':
            every_ping = [x+datetime.timedelta(hours=6) for x in every_ping]

        # Load data
        Sv_mtx,ping_time,bin_size = get_data_based_on_day(every_ping,SITE_CODE,data_path)

        # Clean up data (remove NaNs and fill missing pings)
        Sv_mtx,ping_time,idx_to_delete,day_to_delete = \
            clean_days(max_missing_ping,pings_per_day,Sv_mtx,ping_time)  # Discard days with too many missing pings
        Sv_mtx,ping_time = clean_pings(Sv_mtx,ping_time)  # Fill missing pings
        Sv_mtx = np.delete(Sv_mtx,np.arange(50),axis=1)   # Take out surface layer

        # Set plotting params
        depth_bin_num = Sv_mtx.shape[1]  # number of depth bins
        max_depth = np.round(Sv_mtx.shape[1]*bin_size)
        vec_len_each_day = pings_per_day*depth_bin_num  # length of vector for 1 day

        # Get time and depth stamps
        depth_tick = np.linspace(0,depth_bin_num,5)
        depth_label = ['%d' % d for d in np.linspace(0,max_depth,5)]
        time_tick = np.linspace(0,pings_per_day,5)
        time_label = ['%02d:00' % x for x in np.arange(0,24,6)]
        plot_params = dict(zip(['depth_tick','depth_label','time_tick','time_label'],\
                               [depth_tick,depth_label,time_tick,time_label]))

        # Run NMF
        Sv_vec_3freq = reshape_into_3freq(Sv_mtx,vec_len_each_day)
        nmf_3freq = decomposition.NMF(n_components=n_comp)
        nmf_3freq.fit(Sv_vec_3freq-Sv_vec_3freq.min())
        r_mtx = nmf_3freq.transform(Sv_vec_3freq-Sv_vec_3freq.min())  # tranformation
        v_comps = sep_into_freq(nmf_3freq.components_,pings_per_day,depth_bin_num)  # get nmf components

        # Filename base
        date_range = datetime.datetime.strftime(every_ping[0],'%Y%m%d')+\
                    '-'+datetime.datetime.strftime(every_ping[-1],'%Y%m%d')

        # Plot results and save figs
        fig_comp = plot_decomp_v(v_comps,plot_params)  # Plot all components
        save_fname = '%s_%s_%dcomps_vec.png' % (SITE_CODE,date_range,v_comps.shape[0])
        fig_comp.savefig(os.path.join(save_path,save_fname),dpi=150)
        plt.close(fig_comp)

        fig_tfm = plot_decomp_transform(r_mtx,plot_params)
        save_fname = '%s_%s_%dcomps_tfm.png' % (SITE_CODE,date_range,v_comps.shape[0])
        fig_tfm.savefig(os.path.join(save_path,save_fname),dpi=150)
        plt.close(fig_tfm)










