#!/usr/local/bin/python

# Change filename predix for OOI 80m site

import glob
import os
import re

data_path = '/Volumes/wjlee_apl_2/ooi_zplsc_80m'
data_fname_all = glob.glob(os.path.join(data_path,'*'))

test_str = '(\S*_|^)(?P<FileName>\S*)\.(?P<ExName>\S*)'
test_str_matcher = re.compile(test_str)

fname_all = map(lambda x: os.path.basename(x), data_fname_all)

for fname in fname_all:
    fmatch = test_str_matcher.match(fname)
    print fmatch.group('FileName')+'.'+fmatch.group('ExName')
    new_fname = fmatch.group('FileName')+'.'+fmatch.group('ExName')
    os.rename(os.path.join(data_path,fname), os.path.join(data_path, new_fname))
