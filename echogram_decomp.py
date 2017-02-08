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
import matplotlib.pylab as plt
from modest_image import imshow
from calendar import monthrange

data_path = '/Volumes/wjlee_apl_2/ooi_zplsc_h5_figs'
h5_fname = 'CE04OSPS_201509.h5'

f = h5py.File(os.path.join(data_path,h5_fname))

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

# Function to find nearest element
def find_nearest_time_idx(all_timestamp_num,time_wanted):
    time_wanted_num = date2num(time_wanted)
    idx = np.searchsorted(all_timestamp_num, time_wanted_num, side="left")
    if idx > 0 and (idx == len(all_timestamp_num) or \
        math.fabs(time_wanted_num - all_timestamp_num[idx-1]) < math.fabs(time_wanted_num - all_timestamp_num[idx])):
        return idx-1
    else:
        return idx

# Get f['data_times'] idx for every hour in all days in the month
all_idx = [find_nearest_time_idx(f['data_times'],hr) for hr in every_hr]

# Extract Sv and timing data
Sv_mtx = f['Sv'][:,:,all_idx]
data_times = f['data_times'][all_idx]

# Manipulate Sv data into 1 column for each day
vec_len_each_day = len(all_hr)*len(all_minutes)*Sv_mtx.shape[1]  # length of vector for 1 day
Sv_mtx_rshp = [np.transpose(Sv_mtx[band,:,:]),(-1,vec_len_each_day)) for band in range(3)]
Sv_mtx_rshp = np.asarray(Sv_mtx_rshp)  # convert from list to numpy array
