# starmap-northern

## Overview and features
Star map of the northern hemisphere made with Python packages for astronomy: astropy, ephem, pyorbital. It uses matplotlib to render the graphics. Constellations (more accurately, asterisms) and their stars to plot are determined by a user and defined as Python dictionary object. Then, the dedicated script automatically download necessary data (stars coordinates) (SIMBAD database). The main program puts all constellations forms on the plot and draws over the main plane the field of view: the circle showing the region that observer can see at this time in this geolocation. So, 2 coordinate systems are present: bigger outer circle corresponds to RA-dec (ICRS, equatorial coordinate system) and smaller inner circle goes with alt-az (horizontal coordinate system). Additionally, Sun, Moon, Solar system planets and ISS are placed on the plot (if currently located under the northern hemisphere sky). The app can be used for educational purposes.

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
  - geopy: retrieve geolocation by a string from the Internet
  - BeautifulSoup4: HTML files handling

Desired (the following packages were used while developing):
  - PyQt5
  - IPython and Jupyter
  - Spyder IDE

## Usage
Edit `constellations` dictionary in `dl_constellations.py` file to include/exclude desired stars/constellations. Form of the constellation (asterism) is set by the sequence of stars in which that constellation needs to be drawn to get its form. Think of it like about a task where you need to draw some figure, not taking away the pencil from the paper. So of course you can go through the same points and lines where you've already been. For example, let's see Ursa Major ('UMa' code) constellation. For telling the program to draw its Big Dipper (Plough) asterism we need to add following line at `constellations` `OrderedDict`:
```python
'UMa': ['η','ζ','ε','δ','γ','β','α','δ']
```
We see that `δ UMa` star is included twice as it is point where "dipper" is closing itself.

Also, before starting the main script you can specify some settings e.g. online/local ISS TLE, geolocation, light/dark style, fonts and markers sizes and so on.

To download stars' coordinates and form database files run
```bash
rm stars.html .constellations.html  # start from scratch
python3 dl_constellations.py
```
and wait for completing (it may take tens of minutes). In the end you will see a report that will indicate if some stars hadn't been downloaded. Currently, you must have all stars to be downloaded to successfully make the map. After writing `stars.html` and `.constellations.html` files you can use them for plotting and there is no more needs to start `dl_constellations.py`.

To run main script execute
```bash
python3 sky_map.py
```
but it's recommended to use IPython (in a browser or Spyder IDE) instead of pure Python.
