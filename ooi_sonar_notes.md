# OOI Sonar Data Processing Project

############################################
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

############################################
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

############################################
## 2017/01/26
* `os.getcwd()` gets current direcotry
* Indicator function: `[int(x == '20150911') for x in date_list]`
* map and lambda: `map(lambda x: x.group('Date'),raw_file_times)`
* ways to iterate over dictionary: `freq = {k: str(int(v/1E3))+'k' for (k,v) in frequencies.items()}  # get frequencies`

############################################
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

############################################
## 2017/02/03
* Discussion with Valentina on analysis:
	* PCA, ICA, NMF
	* scikit learn decomposition [tutorial](http://scikit-learn.org/stable/modules/decomposition.html)
	* GIL **"lock"**

############################################
## 2017/02/04-05
* Code to index `data_times` and `Sv` of specific times
* Found out the first ping of each day seems to have smaller amplitude than the rest of the day. **[CONFIRM THIS WITH SKIP]**

############################################
## 2017/02/07-12
* Get PCA, ICA, and NMF working for sonar data
* Spent lot of time figuring out how to use `numpy.reshape()`. It's not very intuitive in that the elements are taken to fill the *last* instead of the *first* dimension of the array.
* `h5py` doesn't allow duplicated indexing as in `numpy`, see [here](https://github.com/h5py/h5py/issues/8)

############################################
## 2017/02/13
* Set up jupyter noteboook server on the APL machine (see resources below for how-to)
* Comparison of NMF analysis on either single frequencies or all 3 frequencies together for 30 days of data: it seems like the decomposition is more stable using all 3 frequencies together.
* Also tried NMF analysis on 80 consecutive days of data: it seems like the patterns are a little “washed away” probably due to larger variability across longer time
* Noticed that the echosounders sometimes skipped pinging at particular instances
* Questions:
	* How to deal with missing data? 1 ping vs multiple hours in a day
	* 




## TO-DO
* Check if decimated power (envelope?) series give the same statistics as the original time series
* Need to take care of the divided by zero warning for TVG



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


## REFERENCES
* Recent paper about krill Sv38 and Sv120 ratio: [Volume backscattering strength of ice krill (Euphausia crystallorophias) in the Amundsen Sea coastal polynya](http://www.sciencedirect.com/science/article/pii/S0967064515002106)
* Correlate backscatteing with oceanographic background changes [Acoustic backscatter observations with implications for seasonal and vertical migrations of zooplankton and nekton in the Amundsen shelf (Antarctica)](http://www.sciencedirect.com/science/article/pii/S0272771414003485)


## Misc Python notes
* the `take` and `put` functions in general have better performance than their fancy indexing equivalents by a significant margin
* `ravel` does not produce a copy of the underlying data if it does not have to (more on this below). The `flatten` method behaves like ravel except it always returns a copy of the data
* Broadcasting has better performance than using `repeat`.
* Broadcasting rule: Two arrays are compatible for broadcasting if for each trailing dimension (that is, starting from the end), the axis lengths match or if either of the lengths is 1. Broadcasting is then performed over the missing and / or length 1 dimensions.
* `np.newaxes` can be used to add a new axis with length 1 specifically for broadcasting purposes, especially in generic algorithms
* Use Pandas:
	* `from pandas import Series, DataFrame`
	* `import pandas as pd`