# OOI Sonar for Ecosystem Monitoring

## Data sources
* Contact point: Friedrich Knuth (knuth@marine.rutgers.edu)
* Official OOI data help: help@oceanobservatories.org

The sonar data files are current only accessible from the [OOI Raw Data Archive](http://oceanobservatories.org/data/raw-data/) and _not_ on the [OOI Data Portal](http://oceanobservatories.org/data-portal/). Individual instruments are indexed using the Reference Designator instrument code shown below:

| Array   | Reference       | Designator             |
|:-------:|:---------------:|:----------------------:|
|Coastal  | Pioneer         | CP01CNSM-MFD37-07-ZPLSC|
|Coastal  | Pioneer	       | CP03ISSM-MFD37-07-ZPLSC|
|Coastal  | Pioneer	       | CP04OSSM-MFD37-07-ZPLSC|
|Coastal  | Endurance	       | CE01ISSM-MFD37-07-ZPLSC|
|Coastal  | Endurance	       | CE06ISSM-MFD37-07-ZPLSC|
|Coastal  | Endurance       | CE07SHSM-MFD37-07-ZPLSC|
|Coastal  | Endurance       | CE09OSSM-MFD37-07-ZPLSC|
|Global   | Argentine Basin | GA02HYPM-MPM01-02-ZPLSG|
|Global   | Irminger        | GI02HYPM-MPM01-02-ZPLSG|
|Global   | Station Pap     | GP02HYPM-MPM01-02-ZPLSG|
|Global   | Southern Ocean  | GS02HYPM-MPM01-02-ZPLSG|
|Coastal  | Endurance	       | CE02SHBP-MJ01C-07-ZPLSC|
|Coastal  | Endurance	       | CE04OSPS-PC01B-05-ZPLSC|

Additional information about the instruments can be found at these locations by searching for the "zpls" instrument code. The data management team at Rutgers are working on centralizing these resources.

Below are a few related websites:

* [OOI Main Website](http://oceanobservatories.org/instruments/)
* [Data Portal Asset Management Page](https://ooinet.oceanobservatories.org/assets/management/)
* [Data Team QA/QC Testing Database](https://ooi.visualocean.net/instruments/all)
* National Center for Environmental Information (NCEI) - [Water column sonar data](Link to the sonar data repository http://www.ngdc.noaa.gov/mgg/wcd/)

## Data type and features
In this project we will focus on the sonar echo data collected using modified Kongberg Simrad EK60 echosounders (those with instrument code ZPLSC-**B**). The files are typically 50 MB in size, each containing ~11 hours of echo data. The data are in _.raw_ format, which is a compressed data format from the manufacturer. We currently have a Matlab and [a Python package](https://github.com/oceanobservatories/mi-instrument/tree/master/mi/instrument/kut/ek60/ooicore) that can unpack the data and plot an echogram (see below).

The data are echo time series of multiple sonar _pings_ transmitted roughly every second. For each ping, the returning echoes form a time series, which can be converted to _depth_ by calculating the two-way travel time of sound from any particular point in the water colume to the sonar transducer. By plotting the returning echo time series across pings, we arrive at the _echograms_ as shown below.

Echogram at 38 kHz:
<figure>
  <img src=".\img\ooi_ex_38k.png width="300">
</figure>

Echogram at 120 kHz:
<figure>
  <img src=".\img\ooi_ex_120k.png width="300">
</figure>

Note the two freuqencies of data were collected nearly simultaneously.


## Other resources
* [Ocean Network Canada](www.oceannetworks.ca/)'s wiki page about their data product [ASL Acoustic Profiler Time Series](https://wiki.oceannetworks.ca/display/DP/24). The ASL echosounders are also used on un-cabled OOI nodes.





