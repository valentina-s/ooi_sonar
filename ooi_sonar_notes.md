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


## TO-DO
* Fix the python script for getting Sv properly calibrated
* Check if decimated power (envelope?) series give the same statistics as the original time series
* Need to take care of the divided by zero warning for TVG

## RESOURCES
* Raw data:
	* [80m site benthic experiment pacakge](https://rawdata.oceanobservatories.org/files/CE02SHBP/MJ01C/ZPLSCB101_10.33.13.7/)
	* [600m site shallow water profiler](https://rawdata.oceanobservatories.org/files/CE04OSPS/PC01B/ZPLSCB102_10.33.10.143/)
* [timing and memory profiling](http://pynash.org/2013/03/06/timing-and-profiling/)
* python [image tutorial](http://matplotlib.org/users/image_tutorial.html)
