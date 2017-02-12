
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
mask[idx_to_delete] = False
univ_idx = univ_idx[mask]

# Params of Sv data
bin_size = bin_size_tmp      # size of each depth bin
depth_bin_num = Sv_mtx.shape[1]  # number of depth bins
max_depth = np.round(Sv_mtx.shape[1]*bin_size)
vec_len_each_day = pings_per_day_tmp*depth_bin_num  # length of vector for 1 day

# Reshape array and perform NMF
Sv_vec_3freq = reshape_into_3freq(Sv_mtx[:,:,univ_idx],vec_len_each_day)

nmf_3freq = decomposition.NMF(n_components=6)
nmf_3freq.fit(Sv_vec_3freq-Sv_vec_3freq.min())
Sv_vec_3freq_r_nmf = nmf_3freq.transform(Sv_vec_3freq-Sv_vec_3freq.min())  # tranformation 
nmf_3freq_comps = sep_into_freq(nmf_3freq.components_,pings_per_day,depth_bin_num)  # get nmf components

plot_decomp_v(nmf_3freq_comps,save_path,plot_params)

plot_decomp_transform(Sv_vec_3freq_r_nmf,save_path,plot_params)


# Get time and depth stamps
depth_tick = np.linspace(0,depth_bin_num,5)
depth_label = ['%d' % d for d in np.linspace(0,max_depth,5)]
time_tick = np.linspace(0,pings_per_day,5)
time_label = ['%02d:00' % x for x in np.arange(0,24,6)]
plot_params = dict(zip(['year','month','depth_label','time_label','depth_tick','time_tick'],\
                        [2015,8,depth_label,time_label,depth_tick,time_tick]))


def get_num_days_pings(h5_fname):
    ''' Get the total number of days and number of pings per day for the given year and month '''
    # Get month and day range
    ym = datetime.datetime.strptime(H5_FILENAME_MATCHER.match(h5_fname).group('YearMonth'),'%Y%m')
    year = ym.year
    month = ym.month
    _,daynum = monthrange(year,month)

    # Get datetime object for on the hour every hour in all days in the month
    all_day = range(1,daynum+1)  # list of all days
    all_hr = range(24)  # list of all hour: 0-23
    all_minutes = range(1,11)  # list of all minutes: 0-9
    every_ping = [datetime.datetime(year,month,day,hr,minutes,0) \
                for day in all_day for hr in all_hr for minutes in all_minutes]
    pings_per_day = len(all_hr)*len(all_minutes)
    return pings_per_day,daynum


def get_data_from_h5(data_path,h5_fname):
    ''' Retrieve data from h5 files '''
    f = h5py.File(os.path.join(data_path,h5_fname),'r')

    # Get month and day range
    ym = datetime.datetime.strptime(H5_FILENAME_MATCHER.match(h5_fname).group('YearMonth'),'%Y%m')
    year = ym.year
    month = ym.month
    _,daynum = monthrange(year,month)

    # Get datetime object for on the hour every hour in all days in the month
    all_day = range(1,daynum+1)  # list of all days
    all_hr = range(24)  # list of all hour: 0-23
    all_minutes = range(1,11)  # list of all minutes: 0-9
    every_ping = [datetime.datetime(year,month,day,hr,minutes,0) \
                for day in all_day for hr in all_hr for minutes in all_minutes]
    pings_per_day = len(all_hr)*len(all_minutes)

    # Get f['data_times'] idx for every hour in all days in the month
    all_idx = [find_nearest_time_idx(f['data_times'],hr) for hr in every_ping]
    all_idx = np.array(all_idx)  # to allow numpy operation

    # Clean up all_idx
    # --> throw away days with more than 5 pings missing
    # --> fill in occasional missing pings with neighboring values
    all_idx_rshp = np.reshape(all_idx,(-1,pings_per_day))
    num_nan_of_day = np.sum(np.isnan(all_idx_rshp),1)
    for day in range(len(num_nan_of_day)):
        if num_nan_of_day[day]>5:
            all_idx_rshp[day,:] = np.nan
        elif num_nan_of_day[day]<5 and num_nan_of_day[day]!=0:
            nanidx = np.where(np.isnan(all_idx_rshp[day,:]))
    #        for ii in nanidx[0]:
    #            if ii==0:
    #                all_idx_rshp[day,ii] = all_idx_rshp[day,ii+1]
    #            else:
    #                all_idx_rshp[day,ii] = all_idx_rshp[day,ii-1]
    all_idx = np.reshape(all_idx_rshp,-1)

    # Extract timing and Sv data
    notnanidx = np.int_(all_idx[~np.isnan(all_idx)])
    data_times = np.empty(all_idx.shape)  # initialize empty array
    data_times[~np.isnan(all_idx)] = f['data_times'][notnanidx.tolist()]
  
    Sv_tmp = f['Sv'][:,:,0]
    Sv_mtx = np.empty((Sv_tmp.shape[0],Sv_tmp.shape[1],all_idx.shape[0]))
    Sv_mtx[:,:,~np.isnan(all_idx)] = f['Sv'][:,:,notnanidx.tolist()]

    bin_size = f['bin_size'][0]      # size of each depth bin

    return Sv_mtx,data_times,bin_size,pings_per_day,all_idx


# Params of Sv data
depth_bin_num = Sv_mtx.shape[1]  # number of depth bins
max_depth = np.round(Sv_mtx.shape[1]*bin_size)
vec_len_each_day = len(all_hr)*len(all_minutes)*depth_bin_num  # length of vector for 1 day

# Get time and depth stamps
depth_tick = np.linspace(0,depth_bin_num,5)
depth_label = ['%d' % d for d in np.linspace(0,max_depth,5)]
time_tick = np.linspace(0,pings_per_day,5)
time_label = ['%02d:00' % x for x in np.arange(0,24,6)]
plot_params = dict(zip(['year','month','depth_label','time_label','depth_tick','time_tick'],\
                        [year,month,depth_label,time_label,depth_tick,time_tick]))

# Decomposition using 3 freq
Sv_vec_3freq=reshape_into_3freq(Sv_mtx,vec_len_each_day)
nmf_3freq = decomposition.NMF(n_components=6)
nmf_3freq.fit(Sv_vec_3freq-Sv_vec_3freq.min())
Sv_vec_3freq_r_nmf = nmf_3freq.transform(Sv_vec_3freq-Sv_vec_3freq.min())  # tranformation
nmf_3freq_comps = sep_into_freq(nmf_3freq.components_,pings_per_day,depth_bin_num)  # get nmf components
plot_decomp_v(nmf_3freq_comps,save_path,plot_params)
plot_decomp_transform(Sv_vec_3freq_r_nmf,save_path,plot_params)

Sv_vec_38k=reshape_into_1freq(Sv_mtx,vec_len_each_day,0)
nmf_38k = decomposition.NMF(n_components=6)
nmf_38k.fit(Sv_vec_38k-Sv_vec_38k.min())
Sv_vec_38k_r_nmf = nmf_38k.transform(Sv_vec_38k-Sv_vec_38k.min())  # tranformation
nmf_38k_comps = sep_into_freq(nmf_38k.components_,pings_per_day,depth_bin_num)  # get nmf components
plot_decomp_v(nmf_38k_comps,save_path,plot_params,38)
plot_decomp_transform(Sv_vec_38k_r_nmf,save_path,plot_params,38)

Sv_vec_120k=reshape_into_1freq(Sv_mtx,vec_len_each_day,1)
nmf_120k = decomposition.NMF(n_components=6)
nmf_120k.fit(Sv_vec_120k-Sv_vec_120k.min())
Sv_vec_120k_r_nmf = nmf_120k.transform(Sv_vec_120k-Sv_vec_120k.min())  # tranformation
nmf_120k_comps = sep_into_freq(nmf_120k.components_,pings_per_day,depth_bin_num)  # get nmf components
plot_decomp_v(nmf_120k_comps,save_path,plot_params,120)
plot_decomp_transform(Sv_vec_120k_r_nmf,save_path,plot_params,120)

Sv_vec_200k=reshape_into_1freq(Sv_mtx,vec_len_each_day,2)
nmf_200k = decomposition.NMF(n_components=6)
nmf_200k.fit(Sv_vec_200k-Sv_vec_200k.min())
Sv_vec_200k_r_nmf = nmf_200k.transform(Sv_vec_200k-Sv_vec_200k.min())  # tranformation
nmf_200k_comps = sep_into_freq(nmf_200k.components_,pings_per_day,depth_bin_num)  # get nmf components
plot_decomp_v(nmf_200k_comps,save_path,plot_params,120)
plot_decomp_transform(Sv_vec_200k_r_nmf,save_path,plot_params,120)
