from mi.instrument.kut.ek60.ooicore.zplsc_b import *
import glob
# import numpy as np   # already imported in zplsc_b
from datetime import datetime
from matplotlib.dates import date2num, num2date
import h5py


def get_data_mtx(data_dict,frequencies):
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


def get_date_idx(date_wanted,fname_all):
    '''
    Index the files in the wanted date range
    '''
    raw_file_times = [FILE_NAME_MATCHER.match(x) for x in fname_all];
    # idx = [x for x in range(len(X)) if X[x].group('Date')=='20150912']  # solution 1
    date_list = map(lambda x: x.group('Date'),raw_file_times)  # if in python 3 need to do list(map()) to convert map object to list
    time_list = map(lambda x: x.group('Time'),raw_file_times)

    if len(date_wanted)==1:
        idx_date = [i for i,dd in enumerate(date_list) if dd==date_wanted[0]]  # solution 2
    elif len(date_wanted)==2:
        idx_date = [i for i,dd in enumerate(date_list) if dd>=date_wanted[0] and dd<=date_wanted[-1]]  # solution 2
    else:
        print 'Invalid range: date_wanted!'
        idx_date = []

    if len(idx_date)!=0:
        if time_list[idx_date[0]]>'000000':  # if midnight was recorded in the previous file
            idx_date.insert(0,idx_date[0]-1)
    else:
        print 'date wanted does not exist!'
    return idx_date


def unpack_raw_to_h5(fname,h5_fname):
    '''
    Unpack *.raw files and save directly into h5 files
    '''
    # unpack data
    particle_data, data_times, power_data_dict, freq, bin_size, config_header, config_transducer = \
        parse_echogram_file(fname)

    # convert from power to Sv
    cal_params = get_cal_params(power_data_dict,particle_data,config_header,config_transducer)
    Sv = power2Sv(power_data_dict,cal_params)  # convert from power to Sv
    Sv_mtx = get_data_mtx(Sv,freq)
    sz = Sv_mtx.shape

    # write to hdf5 file
    f = h5py.File(h5_fname,"a")
    if "Sv" in f:  # if file alread contains Sv mtx
        print 'append new data mtx'
        # append new data
        sz_exist = f['Sv'].shape  # shape of existing Sv mtx
        f['Sv'].resize((sz_exist[0],sz_exist[1],sz_exist[2]+sz[2]))
        f['Sv'][:,:,sz_exist[2]:] = Sv_mtx
        f['data_times'].resize((sz_exist[2]+sz[2],))
        # f['data_times'].resize((1,sz_exist[2]+sz[2]))
        f['data_times'][sz_exist[2]:] = data_times
    else:
        print 'create new dataset'
        # create dataset and save Sv
        f.create_dataset("Sv", sz, maxshape=(sz[0],sz[1],None), data=Sv_mtx, chunks=True)
        f.create_dataset("data_times", (sz[2],), maxshape=(None,), data=data_times, chunks=True)
        # label dimensions
        f['Sv'].dims[0].label = 'frequency'
        f['Sv'].dims[1].label = 'depth'
        f['Sv'].dims[2].label = 'time'
        # create frequency dimension scale, use f['Sv'].dims[0][0][0] to access
        f['frequency'] = freq.values()
        f['Sv'].dims.create_scale(f['frequency'])
        f['Sv'].dims[0].attach_scale(f['frequency'])
    f.close()


def concat_raw(data_path,fname_all,date_wanted):
    ''' Unpack *.raw files and concatenate in memory '''
    # Index the files in the wanted date range
    raw_file_times = [FILE_NAME_MATCHER.match(os.path.join(data_path,x)) for x in fname_all];
    # idx = [x for x in range(len(X)) if X[x].group('Date')=='20150912']  # solution 1
    date_list = map(lambda x: x.group('Date'),raw_file_times)  # if in python 3 need to do list(map()) to convert map object to list
    time_list = map(lambda x: x.group('Time'),raw_file_times)

    if len(date_wanted)==1:
        idx_date = [i for i,dd in enumerate(date_list) if dd==date_wanted[0]]  # solution 2
    elif len(date_wanted)==2:
        idx_date = [i for i,dd in enumerate(date_list) if dd>=date_wanted[0] and dd<=date_wanted[-1]]  # solution 2
    else:
        print 'Invalid range: date_wanted!'
        idx_date = []

    if len(idx_date)!=0:
        if time_list[idx_date[0]]>'000000':  # if midnight was recorded in the previous file
            idx_date.insert(0,idx_date[0]-1)

        # Read first file for initialization
        particle_data, data_times, power_data_dict, freq, bin_size, config_header, config_transducer = \
            parse_echogram_file(os.path.join(data_path,fname_all[idx_date[0]]))
        data_mtx = get_data_mtx(power_data_dict,freq)

        # print 'original size is ' + str({k: v.shape for (k,v) in power_data_dict.items()})
        # freq = frequencies.value()  # get freuqency values
        # freq = {k: str(int(v/1E3))+'k' for (k,v) in frequencies.items()}  # get frequencies

        for fnum in idx_date[1:]:
            print '====== Processing file: ' + fname_all[fnum]
            # Read next file
            _, data_times_tmp, power_data_dict_tmp, freq_tmp, _, _, _ = parse_echogram_file(os.path.join(data_path,fname_all[fnum]))  # load new file
            data_mtx_tmp = get_data_mtx(power_data_dict_tmp,freq_tmp)
            data_mtx = np.concatenate((data_mtx,data_mtx_tmp),axis=1)

    else:
        print 'date wanted does not exist!'

    return data_mtx

        # if power_data_dict[fidx[0]].shape==power_data_dict[fidx[1]].shape==power_data_dict[fidx[2]].shape:
        # power_mtx_temp = np.array((power_data_dict[fidx[0]],\
        #                           power_data_dict[fidx[1]],\
        #                           power_data_dict[fidx[2]]))  # organize all values into matrix
        # else:
        #    sz = np.array(map(lambda x: x.shape,power_mtx_temp))

        # concatenate power data
        # power_data_dict = {k: numpy.concatenate((v,power_data_dict_temp[k]),axis=1) \
        #                            for (k,v) in power_data_dict.items() }   # keep the original dict structure
