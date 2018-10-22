import numpy as np
import os

def read_srad_file(filefullpath):
    '''
    Read txt file downloaded from **reatlime solar radiation data**
    from http://www.ndbc.noaa.gov/station_page.php?station=46098
    '''
    f = open(filefullpath, 'rU')
    header1 = f.readline().split()
    header2 = f.readline().split()
    data_block = f.readlines()
    
    data = {}
    for col_name in header1:
        data[col_name] = np.ma.zeros(len(data_block), 'f',\
                                     fill_value = -999.999)
        
    for (line_count, line) in enumerate(data_block):
        items = line.split()    
        for (col_count, col_name) in enumerate(header1):
            value = items[col_count]
            if value == "MM":
                value = np.ma.masked
            else:
                value = float(value)
            data[col_name][line_count] = value
            
    f.close()
    return data