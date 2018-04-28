# starmap-northern

## Overview and features
Star map of the northern hemisphere made with Python packages for astronomy: astropy, ephem, pyorbital. It uses matplotlib to render the graphics. Constellations (more accurately, asterisms) and their stars to plot are determined by user and defines as Python dictionary object. Then, dedicated script automatically download necessary data (stars coordinates). Main program puts all constellations forms on the plot and draws over the main plane the field of view: the circle showing region that observer can see at this time in this position. So there are 2 coordinate systems: bigger outer circle corresponds to RA-dec (ICRS, equatorial coordinate system) and smaller inner circle goes with alt-az (horizontal coordinate system). Additionally, Sun, Moon, Solar system planets and ISS are placed on the plot (if located in the northern hemisphere). The app can be used for educational purposes.

## Requirements
Requirements are listed in `requirements.txt` file so you can run
```bash
pip3 install -r requirements.txt
```
to install them.
  - astropy: SIMBAD queries, main coordinates transformations, database file storage
  - ephem (pyephem): TLE processing
  - pyorbital: TLE querying and processing
  - matplotlib: currently only 2.0.0 version is supported
  - geopy: retrieve geolocation by string from the Internet
  - BeautifulSoup4: HTML files handling

Desired (the following packages were used while developing):
  - PyQt5
  - IPython and Jupyter
  - Spyder IDE

## Usage
Edit `constellations` dictionary in `dl_constellations.py` file to include/exclude desired stars/constellations.

To download stars' coordinates and form database files run
```bash
rm stars.html .constellations.html  # start from scratch
python3 dl_constellations.py
```
and wait for completing (it may take tens of minutes). In the end you will see a report that will indicate if some stars hadn't been downloaded. Currently, you must have all stars to be downloaded to successfully make the map. After writing `stars.html` and `.constellations.html` files you can use them for plotting and there is no more needs to start `dl_constellations.py`.

To run main script execute:
```bash

```
