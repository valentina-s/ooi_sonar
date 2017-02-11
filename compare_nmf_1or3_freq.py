
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


data_path = '~/internal_2tb/ooi_sonar/ooi_zplsc_h5_figs'  # path to folder with all data
save_path = '~/internal_2tb/ooi_sonar/analysis_figs'
h5_fname = 'CE04OSPS_201509.h5'  # file to analyze
f = h5py.File(os.path.join(data_path,h5_fname),'r')

# Get month and day range
H5_FILENAME_MATCHER = re.compile('(?P<SITE_CODE>\S*)_(?P<YearMonth>\S*)\.\S*')
ym = datetime.datetime.strptime(H5_FILENAME_MATCHER.match(h5_fname).group('YearMonth'),'%Y%m')
year = ym.year
month = ym.month
_,daynum = monthrange(year,month)

# Get datetime object for on the hour every hour in all days in the month
all_day = range(1,daynum+1)  # list of all days
all_hr = range(24)  # list of all hour: 0-23
all_minutes = range(1,11)  # list of all minutes: 0-9
every_hr = [datetime.datetime(year,month,day,hr,minutes,0) \
            for day in all_day for hr in all_hr for minutes in all_minutes]
pings_per_day = len(all_hr)*len(all_minutes)

# Get f['data_times'] idx for every hour in all days in the month
all_idx = [find_nearest_time_idx(f['data_times'],hr) for hr in every_hr]

# Extract Sv and timing data
Sv_mtx = f['Sv'][:,:,all_idx]
data_times = f['data_times'][all_idx]

# Params of Sv data
bin_size = f['bin_size'][0]      # size of each depth bin
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
