#!/Users/wujung/anaconda/envs/py27/bin/python

'''
Plot echogram for 1 day
'''

import os, sys, glob, re
import datetime
import numpy as np
from calendar import monthrange
import h5py

sys.path.insert(0,'/home/wu-jung/code_git/mi-instrument/')
from concat_raw import get_num_days_pings, get_data_from_h5
from echogram_decomp import find_nearest_time_idx

import matplotlib.pyplot as plt
from modest_image import imshow
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Set default colormap
plt.rcParams['image.cmap'] = 'jet'


# Get info of all files and set path
#data_path = '/media/wu-jung/My Passport/OOI/ooi_all_data_h5'
data_path = '/Volumes/wjlee_apl_2/ooi_zplsc_h5_figs/'
fname_form = 'CE04OSPS_2017*.h5'                                                   
fname_all = glob.glob(os.path.join(data_path,fname_form))
#save_path = '/home/wu-jung/internal_2tb/ooi_sonar/figs/plot_echogram_month'
save_path = '/Volumes/wjlee_apl_2/ooi_zplsc_h5_figs/'

H5_FILENAME_MATCHER = re.compile('(?P<SITE_CODE>\S*)_(?P<YearMonth>\S*)\.\S*')


for fname in fname_all:
    print 'Processing '+fname+' ...'

    # Open h5 file
    f = h5py.File(fname,'r')
    
    # Get month and day range
    file_datetime = H5_FILENAME_MATCHER.match(os.path.basename(fname)).group('YearMonth')
    ym = datetime.datetime.strptime(file_datetime,'%Y%m')
    year = ym.year
    month = ym.month
    _,daynum = monthrange(year,month)

    # Save fig name
    save_fname = os.path.basename(fname).split('.')[0]
    
    # Get datetime object for on the hour every hour in all days in the month
    all_day = range(1,daynum+1)  # list of all days         
    all_hr = range(24)  # list of all hour: 0-23
    all_minutes = range(1,21)  # list of all minutes: 1-10
    every_ping = [datetime.datetime(year,month,day,hr,minutes,0) \
                  for day in all_day for hr in all_hr for minutes in all_minutes]
    pings_per_day = len(all_hr)*len(all_minutes)

    # Get f['data_times'] idx for every hour in all days in the month
    all_idx = [find_nearest_time_idx(f['data_times'],hr) for hr in every_ping]
    all_idx = np.array(all_idx)  # to allow numpy operation                                                                 
    # Extract timing and Sv data
    notnanidx = np.int_(all_idx[~np.isnan(all_idx)])
    data_times = np.empty(all_idx.shape)  # initialize empty array
    data_times[~np.isnan(all_idx)] = f['data_times'][notnanidx.tolist()]
    
    Sv_tmp = f['Sv'][:,:,0]
    Sv_mtx = np.empty((Sv_tmp.shape[0],Sv_tmp.shape[1],all_idx.shape[0]))
    Sv_mtx[:] = np.nan
    Sv_mtx[:,:,~np.isnan(all_idx)] = f['Sv'][:,:,notnanidx.tolist()]

    bin_size = f['bin_size'][0]      # size of each depth bin                      

    f.close()

    # Get plotting params
    depth_bin_num = Sv_mtx.shape[1]  # number of depth bins
    max_depth = np.round(Sv_mtx.shape[1]*bin_size)
    day_tick = np.arange(daynum)*pings_per_day
    day_label = [str(x) for x in np.arange(daynum)]
    depth_tick = np.linspace(0,depth_bin_num,5)
    depth_label = ['%d' % d for d in np.linspace(0,max_depth,5)]
    plot_params = dict(zip(['year','month','depth_label','day_label','depth_tick','day_tick'],\
                           [year,month,depth_label,day_label,depth_tick,day_tick]))

    # Plot and save figure
    fig,ax = plt.subplots(3,sharex=True,sharey=True)
    for ff in range(Sv_mtx.shape[0]):
        im = imshow(ax[ff],Sv_mtx[ff,:,:],aspect='auto',vmax=-30,vmin=-80)
        ax[ff].set_yticks(plot_params['depth_tick'])
        ax[ff].set_yticklabels(plot_params['depth_label'],fontsize=16)
        ax[ff].set_xticks(plot_params['day_tick'])
        ax[ff].set_xticklabels(plot_params['day_label'],fontsize=16)
        ax[ff].set_ylabel('Depth (m)',fontsize=20)
        if ff==0:
            ax[ff].set_title(save_fname+', 38 kHz',fontsize=20)
        elif ff==1:
            ax[ff].set_title('120 kHz',fontsize=20)
        elif ff==2:
            ax[ff].set_title('200 kHz',fontsize=20)
        divider = make_axes_locatable(ax[ff])
        cax = divider.append_axes("right", size="1%", pad=0.05)
        cbar = plt.colorbar(im,cax=cax)
        cbar.ax.tick_params(labelsize=14)
        ax[2].set_xlabel('Day',fontsize=20)
    fig.set_figwidth(24)
    fig.set_figheight(8)
    plt.savefig(os.path.join(save_path,save_fname),dpi=300)
    plt.close(fig)


