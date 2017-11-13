
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
from echogram_decomp import find_nearest_time_idx,reshape_into_3freq,reshape_into_1freq,sep_into_freq,plot_decomp_v,plot_decomp_transform


data_path = '/media/wu-jung/My Passport/OOI/ooi_zplsc_h5_figs'  # path to folder with all data
save_path = '/home/wu-jung/internal_2tb/ooi_sonar/figs'

SITE_CODE = 'CE04OSPS'
h5_fname_all = [SITE_CODE+'_'+datetime.datetime(2015,x,1).strftime('%Y%m')+'.h5' for x in range(8,13)]
H5_FILENAME_MATCHER = re.compile('(?P<SITE_CODE>\S*)_(?P<YearMonth>\S*)\.\S*')

# Initialize big matrix that will hold everything
month_pings_days = np.array([get_num_days_pings(x) for x in h5_fname_all])
month_num_pings = [month_pings_days[ii,0]*month_pings_days[ii,1] for ii in range(month_pings_days.shape[0])]
total_num_pings = np.sum([month_pings_days[ii,0]*month_pings_days[ii,1] for ii in range(month_pings_days.shape[0])])
f = h5py.File(os.path.join(data_path,h5_fname_all[0]),'r')
Sv_tmp = f['Sv'][:,:,0]
Sv_mtx = np.empty((Sv_tmp.shape[0],Sv_tmp.shape[1],total_num_pings))
data_times = np.empty(total_num_pings)
all_idx = np.empty(total_num_pings)


# Fill in data from files
file_ping_idx = np.insert(np.cumsum(month_num_pings[:-1]),0,0)
for h5,h5_fname in enumerate(h5_fname_all[0:]):
    print 'Loading file: ' + h5_fname + ' ...'
    Sv_mtx_tmp,data_times_tmp,bin_size_tmp,pings_per_day_tmp,all_idx_tmp = get_data_from_h5(data_path,h5_fname)
    fill_idx = file_ping_idx[h5] + np.arange(Sv_mtx_tmp.shape[2])
    Sv_mtx[:,:,fill_idx] = Sv_mtx_tmp
    data_times[fill_idx] = data_times_tmp
    all_idx[fill_idx] = all_idx_tmp

# Figure out what days of data to discard and fill in short missing pings
all_idx_rshp = np.reshape(all_idx,(-1,pings_per_day))
num_nan_of_day = np.sum(np.isnan(all_idx_rshp),1)
max_missing_ping = 5;

# Determine days and idx with lots of missing pings
day_to_delete = []
idx_to_delete = []
for day in range(len(num_nan_of_day)):
    if num_nan_of_day[day]>max_missing_ping:
        day_to_delete.append(day)
        idx_to_delete.append(day*pings_per_day_tmp+np.arange(0,pings_per_day_tmp))
day_to_delete = np.array(day_to_delete)
idx_to_delete = np.ndarray.flatten(np.array(idx_to_delete))

# Generate idx array with short missing pings filled and bad days removed
univ_idx = np.arange(len(all_idx))

nanidx = np.where(np.isnan(all_idx))
for ii in nanidx[0]:
    if ii==0:
        univ_idx[ii] = univ_idx[ii+1]
    else:
        univ_idx[ii] = univ_idx[ii-1]

mask = np.ones(len(univ_idx), dtype=bool)
if idx_to_delete:
    mask[idx_to_delete] = False
univ_idx = univ_idx[mask]

# Params of Sv data
bin_size = bin_size_tmp      # size of each depth bin
depth_bin_num = Sv_mtx.shape[1]  # number of depth bins
max_depth = np.round(Sv_mtx.shape[1]*bin_size)
vec_len_each_day = pings_per_day_tmp*depth_bin_num  # length of vector for 1 day

# Get time and depth stamps
depth_tick = np.linspace(0,depth_bin_num,5)
depth_label = ['%d' % d for d in np.linspace(0,max_depth,5)]
time_tick = np.linspace(0,pings_per_day,5)
time_label = ['%02d:00' % x for x in np.arange(0,24,6)]
plot_params = dict(zip(['year','month','depth_label','time_label','depth_tick','time_tick'],\
                        [2015,8,depth_label,time_label,depth_tick,time_tick]))

# Reshape array and perform NMF
Sv_vec_3freq = reshape_into_3freq(Sv_mtx[:,:,univ_idx],vec_len_each_day)
nmf_3freq = decomposition.NMF(n_components=6)
nmf_3freq.fit(Sv_vec_3freq-Sv_vec_3freq.min())
Sv_vec_3freq_r_nmf = nmf_3freq.transform(Sv_vec_3freq-Sv_vec_3freq.min())  # tranformation 
nmf_3freq_comps = sep_into_freq(nmf_3freq.components_,pings_per_day,depth_bin_num)  # get nmf components

# Plot NMF results
plot_decomp_v(nmf_3freq_comps,save_path,plot_params)
plot_decomp_transform(Sv_vec_3freq_r_nmf,save_path,plot_params)



