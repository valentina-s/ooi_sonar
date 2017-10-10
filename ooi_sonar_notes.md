# OOI Sonar Data Processing Project



## 2017/01/16
* Found out [Scrapy](https://scrapy.org) is a nice tool for scrapying websites ([tutorial](https://doc.scrapy.org/en/1.3/intro/overview.html)). Seems straightforward to use and make crawling the website tree very easy. Later on when we need to make the developed methods more widely applicable we can consider using this.
* For now since we just need to download a large chunk of data for exploration, it seems sufficient to just use `wget`.
	```
	$ wget -r --no-parent -nd -e robots=off -A raw,bot,idx,png [URL]
	```
	* `-r`: recursive operation
	* `--no-parent`: not indexing the parent directories. Note this will only work if the URL ends with '/' and therefore is recognized as a directory.
	* `-nd`: no directory structure in the downloaded file. This is convenient so that we can just sort and select through the whole pool of files for data from specific dates.
	* `-e robots=off`: force to crawl the website
	* `-A acclist`: same as `--accept`; `acclist` can be a list of file suffixes or patterns, can use `*`, `?`, etc.
* Resource: an [example code](https://github.com/billhowe/ooifetch/blob/master/fetchmovies.py) to look and fetch movies from OOI raw data server.



## 2017/01/24
* `conda info --envs` to see which env is on
* `conda activate py27` switch to python 2.7 environment
* Do `pip install` under the respective environment since the packages are separate
* modify `generate_relative_file_path` in `zplsc_b.py` to work with filename convention
* Get date strings: `from matplotlib.dates import date2num, num2date`
* Time recorded in `data_times` needs to be converted using the following:
   ```
   # Reference time "seconds since 1900-01-01 00:00:00"
   REF_TIME = date2num(datetime(1900, 1, 1, 0, 0, 0))
   (data_times / (60 * 60 * 24)) + REF_TIME
   ```
* About finding peaks in python: [discussion](https://blog.ytotech.com/2015/11/01/findpeaks-in-python/) and associated [github page](https://github.com/MonsieurV/py-findpeaks)
* Getting data loaded and echogram plotted using OOI repository code! :)



## 2017/01/26
* `os.getcwd()` gets current direcotry
* Indicator function: `[int(x == '20150911') for x in date_list]`
* map and lambda: `map(lambda x: x.group('Date'),raw_file_times)`
* ways to iterate over dictionary: `freq = {k: str(int(v/1E3))+'k' for (k,v) in frequencies.items()}  # get frequencies`



## 2017/01/30-02/01
* Use memory profiler:
	```
	%load_ext memory_profiler
	%load_ext line_profiler
	%mprun -f concat_raw concat_raw(data_path,fname_all,date_wanted)
	```
* The xarray format is based on netCDF, which is based on HDF5 in the recent version.
* Matlab has functions to directly import and export HDF5, so we'll use HDF5 for now and devise another structure to store metadata.
* Found out that function `ZPLSPlot` actually transpose and flip `power_data_dict` during `__init__` by calling `_transpose_and_flip(power_data_dict)`, and that's why with and without the plotting routine the data matrix dimensions are changed
* `matplotlib` [subplot demo](http://matplotlib.org/examples/pylab_examples/subplots_demo.html)
* Use `f['Sv2'].dims[0][0][0]` to access the attached dimension scale



## 2017/02/03
* Discussion with Valentina on analysis:
	* PCA, ICA, NMF
	* scikit learn decomposition [tutorial](http://scikit-learn.org/stable/modules/decomposition.html)
	* GIL **"lock"**



## 2017/02/04-05
* Code to index `data_times` and `Sv` of specific times
* Found out the first ping of each day seems to have smaller amplitude than the rest of the day. **[CONFIRM THIS WITH SKIP]**



## 2017/02/07-12
* Get PCA, ICA, and NMF working for sonar data
* Spent lot of time figuring out how to use `numpy.reshape()`. It's not very intuitive in that the elements are taken to fill the *last* instead of the *first* dimension of the array.
* `h5py` doesn't allow duplicated indexing as in `numpy`, see [here](https://github.com/h5py/h5py/issues/8)



## 2017/02/13
* Set up jupyter noteboook server on the APL machine (see resources below for how-to)
* Comparison of NMF analysis on either single frequencies or all 3 frequencies together for 30 days of data: it seems like the decomposition is more stable using all 3 frequencies together.
* Also tried NMF analysis on 80 consecutive days of data: it seems like the patterns are a little “washed away” probably due to larger variability across longer time
* Noticed that the echosounders sometimes skipped pinging at particular instances



## 2017/02/20
Summarize discussions this and last week:

* Discussed with Sasha about rPCA/rNMF: he used different cost functions and different formulation than the PCA pursuit which doesn't allow any parameter tuning. He's currently debugging the Julia code, which will allow us to tune the parameter to get different "sparse" components.
* The ["multi-resolution dynamic mode decomposition"](https://arxiv.org/pdf/1506.00564.pdf) can potentially be used to detect changes in our data. In the example they use the SST data to discover the El Niño event in 1997. Can we use it to find some usual/unusual patterns in the echogram?
* Thinking about using [subspace tracking](http://web.eecs.umich.edu/~girasole/?page_id=190) to detect changes in the echogram? maybe in the long run...



## 2017/02/28
* rPCA to separate background and sparse components, and decompose background component
* remove surface wave from the decomposition
* use sliding window for NMF and use previous components to "seed" the following decomposition --> use results of this to see if can derive/observe longer temporal scale patterns



## 2017/03/02
* correlation with external variable using NMF?



## 2017/03/07
* Modules to eliminate nan pings and spit out nan locations and fixed matrix:
	* `get_data_based_on_day`
	* `clean_days`
	* `clean_pings`
* Do sliding window NMF



## 2017/03/09
* Discussion with Valentina and Sasha on how to structure the optimization problem
* Make list of specific goals for follow-up after incubator
* Need to: plot components across all 3 frequencies with same caxis



## 2017/09/11
* Pick back up the EK60 .RAW unpacking code. Check what was changed/added from the forked mi-instrument package at that time.
* Changed `parse_echogram_file` and `read_header`: added output `config_header` and `config_transducer`, which are used in new function `get_cal_params`.
* `config_header` outputs are from function `read_config_header` from `zplsc_echogrampy` that was supplied in the forked mi-instrument package.
* `config_transducer` were added by WJL and read transducer parameters using `read_config_transducer` supplied by the forked mi-instrument package
* The various parameters in `cal_params` are unpacked from `config_header`, `config_transducer`, and `particle_data`. `particle_data` were unpacked from .raw file using `parse_echogram_file`
* Also added `data_times = (data_times / (60 * 60 * 24)) + REF_TIME` toward the end of `parse_echogram_file` and transpose the variable `power_data_dict`. This conversion was produced referencing `ZPLSPlot` in `zplsc_echogram.py`
* Params to add (? or not??):
	* to `cal_param` (corresponds to the `ping` structure array in Matlab Echolab):
		- bandwidth
		- heave/roll/pitch
		- temperature
		- trawl*
		- offset
		- count
		- power
		- alongship
		- athwartship
		- samplerange
		- seg
	* **The above data are stored in `sample_data` unpacked by function `process_sample`, can be seen from `sample_dtype`, just not extracted by the current `parse_echogram_file`**
* Variable correspondence:
	* `config_transducer` -- `data.config`
	* `cal_param` -- part of `data.ping`
* The mi-instrument implementation `parse_echogram_file` only includes the `RAW0` datagram and ignore all others found in `readEKRaw.m`
* The various parameters in `cal_params` are unpacked from `config_header`, `config_transducer`, and `particle_data`. `particle_data` were unpacked from .raw file using `parse_echogram_file`
* In `parse_echogram_file`, `sample_data` refers to the params and `power_data` refers to the actual echo return data
* **CHECK** **(V)**: in `get_cal_params` if using `particle_data[0]` is adequate or if need data from other index --> perhaps the values are identical? --> No, this is because somehow `particle_data` was saved as a tuple with only 1 element, params for different frequencies were accessed using `particle_data[0]['param_name'][freq_num]`
* **CHECK** **(V)**: power2Sv routine between Matlab and Python



## 2017/09/12
* Download all data from CE04OSPS and CE02SHBP sites (use `wget -c -N` options to do incremental mirroring and fix any previous incomplete files)
* Finished checking power2Sv routine in Matlab and Python --> good match!
* **CHECK** **(V)** what are the actual ranges returned by the following section in Matlab `readEKRaw_ConvertPower.m`
```matlab
%  create range vector (in m)
data.pings(n).range = double((0:pSize(1) - 1) + ...
	double(data.pings(n).samplerange(1)) - 1)' * dR;
```

The above corresponds to `range_vec = np.arange(pSize[0]) * dR` in Python `power2Sv` function.

* **Simplify** **(V)** the following section in Python `power2Sv` using broadcasting
```
Sv[n+1] = power_data_dict[n+1] + np.transpose(np.matlib.repmat(TVG,pSize[1],1)) +\
	2*cal_params[n]['absorptioncoefficient']*np.transpose(np.matlib.repmat(rangeCorrected,pSize[1],1)) - CSv - Sac
```
to
```
Sv[n+1] = (power_data_dict[n+1].T \
             +TVG +2*cal_params[n]['absorptioncoefficient']*rangeCorrected\
              -CSv -Sac).T
```
* Temperature readings in the .RAW files are all identical --> need to find other sources to estimate seawater absorption
* **NEXT TO-DOs**:
	* Find data from profiler and estimate seawater absorption at OOI locations --> go to [this page](http://oceanobservatories.org/site/ce04osps/) to use CTD data from both the profiler and the platform
	* Implement background noise estimation/remove from De Robertis & Higginbottom 2007
	* Think about the bin size to average data for MVBS



## 2017/09/14-17
* Finished reading De Robertis & Higginbottom 2007 yesterday, will test implementing the noise removal algorithm today (20170914)
* It took quite a few steps to make it right, including:
    1. Determine the number of pings used for getting averaged noise power
    2. making sure the 2 methods (substract then compensate or the other way around) are equivalent to each other
* The results and functions are in notebooks:
    * `Noise removal (during development, results not correct).ipynb` contains correct part 1. results, but part 2. algorithms were somehow broken
    * `Noise removal.ipynb` has everything in the correct form
* Need to determine a Sv thresholding method before calculating $$\Delta S_V$$: tested 2 methods:
    * De Robertis & Higginbottom 2007: using the noise estimated as in `Noise removal.ipynb`
    * Logerwill & Wilson 2004: simple thresholding approach, also used in McKelvey 2004, Jech & Michaels 2006, Sato et al. 2015
* Lower SNR threshold needs to be used to retain more echogram points for data from 20150910 (SNR threshould=6) compared to data from 20170910 (SNR threshold=10). However, based on the method by Logermill & Wilson, the noise threshold should also be higher (-67 dB for 20150910 and -79 dB for 20170910).



## 2017/09/18
* First kind of successful try of dB-differencing: in `Freq-differencing.ipynb`
	* Spent lots of time figuring out how to handle pixles with NaN values in one or more frequencies
	* Use 2 methods:
		1. Use thresholding between 200 and 38 kHz (following in Sato et al. 2015)
		2. Use color-coding by presence of particular combination of frequencies (follwoing Jech and Michaels 2006)
* **CHECK** Need to check the simple Sv thresholding code again: see if can have a threshold selection criteria based on achieving stable NASC values (method used in Sato et al. 2015)
* **UPDATE** **(V)** functions related to noise estimation and MVBS calculation using masked array functionality


## 2017/10/09
* Update noise estimation functions and related operations using masked array functionality. The new frequency-differencing results are in `Freq-differencing 20171009.ipynb`.
* Add `db_diff.py`, which contains the functions and colormap definitions in `Freq-differencing 20171009.ipynb`.





## TO-DO
* Check if decimated power (envelope?) series give the same statistics as the original time series
* Need to take care of the divided by zero warning for TVG
* Histogram statistics for the echogram
* Do decomposition without the surface layer
* Topic modeling--> finding the best number of components



## Goals
* Methods to determine the best number of components
* Continuity of components (within and across windows)
* Constraints for slow variation in both components and weights
* Constraints for force using spatial relationship in the image pixels (i.e., location is important)
* Optimize extracting data from H5 files
* Run PCP before NMF?


## RESOURCES
* Raw data:
	* [80m site benthic experiment pacakge](https://rawdata.oceanobservatories.org/files/CE02SHBP/MJ01C/ZPLSCB101_10.33.13.7/)
	* [600m site shallow water profiler](https://rawdata.oceanobservatories.org/files/CE04OSPS/PC01B/ZPLSCB102_10.33.10.143/)
* Python:
	* [timing and memory profiling](http://pynash.org/2013/03/06/timing-and-profiling/)
	* [image tutorial](http://matplotlib.org/users/image_tutorial.html)
	* [matplotlib subplots](http://matplotlib.org/examples/pylab_examples/subplots_demo.html)
	* [colormaps](http://matplotlib.org/users/colormaps.html)
* [Extracting and Enriching Ocean Biogeographic Information System (OBIS) Data with R](https://ropensci.org/blog/blog/2017/01/25/obis)
* Scikit learn decomposition:
	* [2.5. Decomposing signals in components (matrix factorization problems)](http://scikit-learn.org/stable/modules/decomposition.html)
	* [Practical Guide to Principal Component Analysis (PCA) in R & Python](https://www.analyticsvidhya.com/blog/2016/03/practical-guide-principal-component-analysis-python/)
* Principal Component Pursuit algorithm [python implementation](https://github.com/dfm/pcp)
* Running ipython notebook on remote machine:
	* first set up a public server on remote machine using [these steps](http://jupyter-notebook.readthedocs.io/en/latest/public_server.html):
		* Prerequisite: A notebook configuration file
		* Preparing a hashed password
		* Using SSL for encrypted communication
		* Running a public notebook server
	* then connect using `https://YOUR_IP:PORT/`
	* also tried [this](https://coderwall.com/p/ohk6cg/remote-access-to-ipython-notebooks-via-ssh) but didn't work
	* Note the self-generated SSL **won't work with Safari** but is fine with other browsers
* Machine learning [next steps](https://www.quora.com/I-have-completed-Andrew-Ngs-Coursera-class-on-Machine-Learning-What-should-I-do-next-What-*can*-I-do-next)
* [Quick HDF5 with Pandas](http://glowingpython.blogspot.com/2014/08/quick-hdf5-with-pandas.html)
* Various NMFs:
	* [NIMFA](http://nimfa.biolab.si/)
	* [Decomposing Time Series Data by a Non-negative Matrix Factorization Algorithm with Temporally Constrained Coefficients](http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7319146)
	* [Exploting long-term temporal dependencies in NMF using recurrent neural networks with application to source separation](https://ccrma.stanford.edu/~gautham/Site/Publications_files/lewandowski-icassp2014.pdf)
	* [NMF with spectral and temporal continuity criteria for monaural sound source separation](http://www.eurasip.org/Proceedings/Eusipco/Eusipco2014/HTML/papers/1569921699.pdf)
	* [Dynamic NMF slides](https://www.uni-oldenburg.de/fileadmin/user_upload/mediphysik/ag/sigproc/audio/nmoh/downloads/D-NMF.pdf)
* Dynamic NMF and temporal continuity:
	* Derek Greene: [Dynamic Topic Modeling](https://github.com/derekgreene/dynamic-nmf), [Temporal stability](https://github.com/derekgreene/topic-stability),[How Many Topics? Stability Analysis for Topic Models](https://arxiv.org/abs/1404.4606)
* [A Unified View of Static and Dynamic Source
Separation Using Non-Negative Factorizations](https://ccrma.stanford.edu/~gautham/Site/Publications_files/smaragdis-spm2014.pdf)
* [23 types of regressions](http://www.datasciencecentral.com/profiles/blogs/23-types-of-regression)
* DMD in [pyrunner](http://www.pyrunner.com/weblog/2016/07/25/dmd-python/)
* [Sparse DMD](https://github.com/aaren/sparse_dmd/blob/master/sparse_dmd/dmd.py)

## Possible funding sources:
* NSF [ABI](https://www.nsf.gov/bio/dbi/about.jsp)
* NSF [BIGDATA](https://www.nsf.gov/funding/pgm_summ.jsp?pims_id=504767)


## REFERENCES
* Recent paper about krill Sv38 and Sv120 ratio: [Volume backscattering strength of ice krill (Euphausia crystallorophias) in the Amundsen Sea coastal polynya](http://www.sciencedirect.com/science/article/pii/S0967064515002106)
* Correlate backscatteing with oceanographic background changes [Acoustic backscatter observations with implications for seasonal and vertical migrations of zooplankton and nekton in the Amundsen shelf (Antarctica)](http://www.sciencedirect.com/science/article/pii/S0272771414003485)
* [Multi-resolution dynamic mode decomposition](https://arxiv.org/pdf/1506.00564.pdf) 
* [Subspace tracking](http://web.eecs.umich.edu/~girasole/?page_id=190) algorithms

## Misc Python notes
* the `take` and `put` functions in general have better performance than their fancy indexing equivalents by a significant margin
* `ravel` does not produce a copy of the underlying data if it does not have to (more on this below). The `flatten` method behaves like ravel except it always returns a copy of the data
* Broadcasting has better performance than using `repeat`.
* Broadcasting rule: Two arrays are compatible for broadcasting if for each trailing dimension (that is, starting from the end), the axis lengths match or if either of the lengths is 1. Broadcasting is then performed over the missing and / or length 1 dimensions.
* `np.newaxes` can be used to add a new axis with length 1 specifically for broadcasting purposes, especially in generic algorithms
* Use Pandas:
	* `from pandas import Series, DataFrame`
	* `import pandas as pd`
