#!/usr/local/bin/python

'''
Concatenate multiple files and save into hdf5 files
'''

import os, sys, glob
import datetime
sys.path.insert(0,'/Users/wujung/code/mi-instrument/')

from mi.instrument.kut.ek60.ooicore.zplsc_b import *
from concat_raw_output import *
from calendar import monthrange
import numpy as np

month_first_last = []
year_want = 2016
for month in np.arange(12)+1:
    _,daynum = monthrange(year_want,month)
    month_first_last.append([datetime.date(year_want,month,1).strftime('%Y%m%d'),\
                             datetime.date(year_want,month,daynum).strftime('%Y%m%d')])

# 80m site
data_path = '/Volumes/wjlee_apl_2/ooi_zplsc_80m/'
SITE_CODE = 'CE02SHBP'
# SITE_CODE = 'CE04OSPS';
fname_form = 'OOI-D*.raw'  # index all files in 2015
fname_all = glob.glob(os.path.join(data_path,fname_form))

for month in np.arange(7,12):
    idx_date = get_date_idx(month_first_last[month],fname_all)
    h5_fname = os.path.join('/Volumes/wjlee_apl_2/ooi_zplsc_h5_figs',
                SITE_CODE+'_'+datetime.date(year_want,month+1,1).strftime('%Y%m')+'.h5')
    for idx in idx_date:
        unpack_raw_to_h5(fname_all[idx],h5_fname)


# date_wanted = ['20150911','20150920']
# idx_date = get_date_idx(date_wanted,fname_all)
# h5_fname = '/Volumes/wjlee_apl_2/2015_0911-20_80m_sub5.h5'
# for idx in idx_date:
#     unpack_raw_to_h5(fname_all[idx],h5_fname,deci_len)

# date_wanted = ['20150921','20150930']
# idx_date = get_date_idx(date_wanted,fname_all)
# h5_fname = '/Volumes/wjlee_apl_2/2015_0921-30_80m_sub5.h5'
# for idx in idx_date:
#     unpack_raw_to_h5(fname_all[idx],h5_fname,deci_len)
#
#
# ''' 600m site '''
# data_path = '/Volumes/wjlee_apl_2/ooi_zplsc_600m/'
# fname_form = 'OOI-D*.raw'  # index all files in 2015
# fname_all = glob.glob(os.path.join(data_path,fname_form))
#
# date_wanted = ['20150901','20150910']
# idx_date = get_date_idx(date_wanted,fname_all)
# h5_fname = '/Volumes/wjlee_apl_2/2015_0901-10_600m_sub5.h5'
# for idx in idx_date:
#     unpack_raw_to_h5(fname_all[idx],h5_fname,deci_len)
#
# date_wanted = ['20150911','20150920']
# idx_date = get_date_idx(date_wanted,fname_all)
# h5_fname = '/Volumes/wjlee_apl_2/2015_0911-20_600m_sub5.h5'
# for idx in idx_date:
#     unpack_raw_to_h5(fname_all[idx],h5_fname,deci_len)
#
# date_wanted = ['20150921','20150930']
# idx_date = get_date_idx(date_wanted,fname_all)
# h5_fname = '/Volumes/wjlee_apl_2/2015_0921-30_600m_sub5.h5'
# for idx in idx_date:
#     unpack_raw_to_h5(fname_all[idx],h5_fname,deci_len)
