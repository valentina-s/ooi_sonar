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
* 
