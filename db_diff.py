

#import os, sys, glob
#from datetime import datetime
#sys.path.insert(0,'/home/wu-jung/code_git/mi-instrument')

#from mi.instrument.kut.ek60.ooicore.zplsc_b import *
#from concat_raw import *

import numpy as np
import datetime as dt
from matplotlib.dates import date2num
import matplotlib.colors as colors


# Colormap from Jech & Michaels 2006
multifreq_th_colors = np.array([[0,0,0],\
                                [86,25,148],\
                                [28,33,179],\
                                [0,207,239],\
                                [41,171,71],\
                                [51,204,51],\
                                [255,239,0],\
                                [255,51,0]])/255.
mf_cmap = colors.ListedColormap(multifreq_th_colors)
mf_bounds = range(9)
mf_norm = colors.BoundaryNorm(mf_bounds,mf_cmap.N)



# Colormap: standard EK60
e_cmap_colors = np.array([[255, 255, 255],\
                          [159, 159, 159],\
                          [ 95,  95,  95],\
                          [  0,   0, 255],\
                          [  0,   0, 127],\
                          [  0, 191,   0],\
                          [  0, 127,   0],\
                          [255, 255,   0],\
                          [255, 127,   0],\
                          [255,   0, 191],\
                          [255,   0,   0],\
                          [166,  83,  60],\
                          [120,  60,  40]])/255.
e_cmap_th = [-80,-30]
e_cmap = colors.ListedColormap(e_cmap_colors)
e_bounds = np.linspace(e_cmap_th[0],e_cmap_th[1],e_cmap.N+1)
e_norm = colors.BoundaryNorm(e_bounds,e_cmap.N)



def get_cal_params(power_data_dict,particle_data,config_header,config_transducer):
    """
    Get calibration params from the unpacked file
    Parameters come from config_header and config_transducer (both from header),
                                      as well as particle_data (from .RAW file)
    """
    cal_params = []
    for ii in range(len(power_data_dict)):
        cal_params_tmp = {}
        cal_params_tmp['soundername'] = config_header['sounder_name'];
        cal_params_tmp['frequency'] = config_transducer[ii]['frequency'];
        cal_params_tmp['soundvelocity'] = particle_data[0]['zplsc_sound_velocity'][ii];
        cal_params_tmp['sampleinterval'] = particle_data[0]['zplsc_sample_interval'][ii];
        cal_params_tmp['absorptioncoefficient'] = particle_data[0]['zplsc_absorption_coeff'][ii];
        cal_params_tmp['gain'] = config_transducer[ii]['gain']    # data.config(n).gain;
        cal_params_tmp['equivalentbeamangle'] = config_transducer[ii]['equiv_beam_angle']   # data.config(n).equivalentbeamangle;
        cal_params_tmp['pulselengthtable'] = config_transducer[ii]['pulse_length_table']   # data.config(n).pulselengthtable;
        cal_params_tmp['gaintable']  = config_transducer[ii]['gain_table']   # data.config(n).gaintable;
        cal_params_tmp['sacorrectiontable'] = config_transducer[ii]['sa_correction_table']   # data.config(n).sacorrectiontable;
        cal_params_tmp['transmitpower'] = particle_data[0]['zplsc_transmit_power'][ii]   # data.pings(n).transmitpower(pingNo);
        cal_params_tmp['pulselength'] = particle_data[0]['zplsc_pulse_length'][ii]   # data.pings(n).pulselength(pingNo);
        cal_params_tmp['anglesensitivityalongship'] = config_transducer[ii]['angle_sensitivity_alongship']  # data.config(n).anglesensitivityalongship;
        cal_params_tmp['anglesensitivityathwartship'] = config_transducer[ii]['angle_sensitivity_athwartship']   #data.config(n).anglesensitivityathwartship;
        cal_params_tmp['anglesoffsetalongship'] = config_transducer[ii]['angle_offset_alongship']   # data.config(n).anglesoffsetalongship;
        cal_params_tmp['angleoffsetathwartship'] = config_transducer[ii]['angle_offset_athwart']   # data.config(n).angleoffsetathwartship;
        cal_params_tmp['transducerdepth'] = particle_data[0]['zplsc_transducer_depth'][ii]   # data.pings(n).transducerdepth(pingNo);
        cal_params.append(cal_params_tmp)
    return cal_params



def get_noise(power_data,depth_bin_size,ping_bin_range,depth_bin_range,tvgCorrectionFactor=2):
    '''
    INPUT:
        power_data            2D mtx of power data [depth x ping num]
        ping_bin_range        average over M pings
        depth_bin_range       average over depth_bin_range [m]
        tvgCorrectionFactor   default (=2) is to apply TVG correction with offset of 2 samples
                              note this factor is important in TVG compensation
                              and therefore in how power_bin is obtained as well
    OUTPUT:
        minimum value for bins of averaged ping
    '''
    N = int(np.floor(depth_bin_range/depth_bin_size))
    
    # Average uncompensated power over M pings and N depth bins
    depth_bin_num = int(np.floor((power_data.shape[0]-tvgCorrectionFactor)/N))
    ping_bin_num = int(np.floor(power_data.shape[1]/ping_bin_range))
    power_bin = np.empty([depth_bin_num,ping_bin_num])
    for iD in range(depth_bin_num):
        for iP in range(ping_bin_num):
            depth_idx = np.arange(N)+N*iD+tvgCorrectionFactor  # match the 2-sample offset
            ping_idx = np.arange(ping_bin_range)+ping_bin_range*iP
            power_bin[iD,iP] = np.mean(10**(power_data[np.ix_(depth_idx,ping_idx)]/10))

    # Noise = minimum value for each averaged ping
    return np.min(power_bin,0),ping_bin_num



def remove_noise(power_data,cal,noise_est,ping_bin_range=40,tvg_correction_factor=2):
    '''
    Function for noise removal and TVG + absorption compensation
    Ref: De Robertis et al. 2010
    
    INPUT:
        power_data      2D mtx of power data [depth x ping num]
        noise_est       results from `get_noise`
                        if pass in only a scalar, it will be treated as the noise estimate for all data
                        if pass in a vecotr, the noise will be removed adaptively using values in the vector
        ping_bin_range          average over M pings
        depth_bin_range         average over depth_bin_range [m]
        tvg_correction_factor   default(=2) for converting power_data to Sv
    OUTPUT:
        Sv_raw       TVG and absorption compensated Sv data, no noise removal
        Sv_corr      TVG and absorption compensated Sv data, no noise removal        
        Sv_noise     TVG and absorption compensated noise estimation
    '''

    # Get cal params
    f = cal['frequency']
    c = cal['soundvelocity']
    t = cal['sampleinterval']
    alpha = cal['absorptioncoefficient']
    G = cal['gain']
    phi = cal['equivalentbeamangle']
    pt = cal['transmitpower']
    tau = cal['pulselength']

    # key derived params
    dR = c*t/2   # sample thickness
    wvlen = c/f  # wavelength

    # Calc gains
    CSv = 10 * np.log10((pt * (10**(G/10))**2 * wvlen**2 * c * tau * 10**(phi/10)) / (32 * np.pi**2))

    # calculate Sa Correction
    idx = [i for i,dd in enumerate(cal['pulselengthtable']) if dd==tau]
    Sac = 2 * cal['sacorrectiontable'][idx]

    # Get TVG
    range_vec = np.arange(power_data.shape[0]) * dR
    rangeCorrected = range_vec - (tvg_correction_factor * dR)
    rangeCorrected[rangeCorrected<0] = 0

    TVG = np.empty(rangeCorrected.shape)
    TVG[rangeCorrected!=0] = \
        np.real( 20*np.log10(rangeCorrected[rangeCorrected!=0]) )  # TVG = real(20 * log10(rangeCorrected));
    TVG[rangeCorrected==0] = 0

    # Get absorption
    ABS = 2*alpha*rangeCorrected

    # Compensate measurement for noise and corrected for transmission loss
    # also estimate Sv_noise component for subsequent SNR check

    # Noise removal
    if noise_est.shape==():  # if noise_est is single element
        subtract = 10**(power_data/10)-noise_est
        tmp = 10*np.log10(np.ma.masked_less_equal(subtract,0))
        tmp.set_fill_value(-999)
        Sv_corr = (tmp.T+TVG+ABS-CSv-Sac).T
        Sv_noise = 10*np.log10(noise_est)+TVG+ABS-CSv-Sac

    else:
        ping_bin_num = int(np.floor(power_data.shape[1]/ping_bin_range))
        Sv_corr = np.ma.empty(power_data.shape)   # log domain corrected Sv
        Sv_noise = np.empty(power_data.shape)  # Sv_noise
        for iP in range(ping_bin_num):
            ping_idx = np.arange(ping_bin_range) +iP*ping_bin_range
            subtract = 10**(power_data[:,ping_idx]/10) -noise_est[iP]
            tmp = 10*np.log10(np.ma.masked_less_equal(subtract,0))
            tmp.set_fill_value(-999)
            Sv_corr[:,ping_idx] = (tmp.T +TVG+ABS-CSv-Sac).T
            Sv_noise[:,ping_idx] = np.array([10*np.log10(noise_est[iP])+TVG+ABS-CSv-Sac]*ping_bin_range).T
    
    # Raw Sv withour noise removal but with TVG/absorption compensation
    Sv_raw = (power_data.T+TVG+ABS-CSv-Sac).T
    
    return Sv_raw,Sv_corr,Sv_noise



def get_MVBS(Sv,depth_bin_size,ping_bin_range,depth_bin_range):
    '''
    Obtain mean MVBS
    
    INPUT:
        th                Sv threshold: discard Sv values below th during averaging
        depth_bin_size    depth bin size from unpacked data
        ping_bin_range    average over M pings
        depth_bin_range   average over depth_bin_range [m]
    OUTPUT:
        smoothed Sv data
    '''

    N = int(np.floor(depth_bin_range/depth_bin_size))  # total number of depth bins
    
    # Average Sv over M pings and N depth bins
    depth_bin_num = int(np.floor(Sv.shape[1]/N))
    ping_bin_num = int(np.floor(Sv.shape[2]/ping_bin_range))
    MVBS = np.ma.empty([Sv.shape[0],depth_bin_num,ping_bin_num])
    for iF in range(Sv.shape[0]):
        for iD in range(depth_bin_num):
            for iP in range(ping_bin_num):
                depth_idx = np.arange(N) + N*iD
                ping_idx = np.arange(ping_bin_range) + ping_bin_range*iP
                MVBS[iF,iD,iP] = 10*np.log10( np.nanmean(10**(Sv[np.ix_((iF,),depth_idx,ping_idx)]/10)) )            
    return MVBS



def multifreq_color_code(Sv38,Sv120,Sv200):
    '''
    Multi-frequency color-coding regiem
    Ref: Jech and Michaels 2006
    
    INPUT:
        Sv at 38, 120, and 200 kHz
        
    OUTPUT:
        Sv_mf   numerical indicator mtx of combination of
                presence/absence of multiple freq
    '''
    yes_1 = ~np.isnan(Sv38)
    yes_2 = ~np.isnan(Sv120)
    yes_3 = ~np.isnan(Sv200)

    Sv_mf = np.empty(Sv38)
    Sv_mf[~yes_1 & ~yes_2 & ~yes_3] = 0
    Sv_mf[ yes_1 & ~yes_2 & ~yes_3] = 1
    Sv_mf[ yes_1 &  yes_2 & ~yes_3] = 2
    Sv_mf[ yes_1 &  yes_2 &  yes_3] = 3
    Sv_mf[ yes_1 & ~yes_2 &  yes_3] = 4
    Sv_mf[~yes_1 &  yes_2 & ~yes_3] = 5
    Sv_mf[~yes_1 &  yes_2 &  yes_3] = 6
    Sv_mf[~yes_1 & ~yes_2 &  yes_3] = 7

    return Sv_mf



def get_power_data_mtx(data_dict,frequencies):
    '''
    Convert data_dict to numpy array
    Input:
        data_dict     power_data_dict or Sv from power2Sv()
        frequencies   unpacked dict from parse_echogram_file()
    '''
    fval = frequencies.values()
    fidx = sorted(range(len(fval)), key=lambda k: fval[k])   # get key sequence for low to high freq
    fidx = [x+1 for x in fidx]
    return np.array((data_dict[fidx[0]],\
                     data_dict[fidx[1]],\
                     data_dict[fidx[2]]))  # organize all values into matrix



def find_nearest_time_idx(all_timestamp_num,time_wanted,tolerance):
    '''
    Function to find nearest element
    time_wanted is a datetime object
    tolerance is the max tolerance in second allowed between `time_wanted` and `all_timestamp`
    all_timestamp_num is a numerical date object (i.e., output from `date2num`)
    '''
    time_wanted_num = date2num(time_wanted)
    idx = np.searchsorted(all_timestamp_num, time_wanted_num, side="left")
    if idx > 0 and (idx == len(all_timestamp_num) or \
        np.abs(time_wanted_num - all_timestamp_num[idx-1]) < np.abs(time_wanted_num - all_timestamp_num[idx])):
        idx -= 1

    # If interval between the selected index and time wanted > `tolerance` seconds
    sec_diff = dt.timedelta(all_timestamp_num[idx]-time_wanted_num).total_seconds()
    if np.abs(sec_diff)>tolerance:
        return np.nan
    else:
        return idx



