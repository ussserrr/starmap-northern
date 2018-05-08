#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# TODO: different stars sizes dependending on their flux
# TODO: map where only FOV and center is a current 'polar star' analogue (zenith)
# TODO: case when not all stars were downloaded
# TODO: implement constellations as graphs
# TODO: ask for lateset matplotlib fix


import sys, os

import numpy as np

import matplotlib.pyplot as plt

import astropy.time
from astropy.coordinates import EarthLocation
import astropy.units as u
import astropy.table

from datetime import datetime, timedelta
import time

from util import extract_forms, extract_constellations, prepare_skymap, plot_fov, plot_moon,\
                 plot_sun, plot_solarsystem, plot_iss, fonts, plot_size
from dl_constellations import constellations as constellations_dict


#
# Set desired location here. Specify name (required for caption), choose online/offline
# flag and in the last case also set latitude, longitude and altitude in the class Loc
#
location_str = "Syktyvkar"
internet_geolocation = 1
if not internet_geolocation:
    class Loc:
        """
        Hard-coded location. You can use it instead of Internet geolocationing
        """
        latitude = 61.40
        longitude = 50.49
        altitude = 100
        address = location_str


#
# Leave custom_time=None or custom_time=0 to use current OS time or specify
# the custom time as datetime(year, month, day, hours, minutes, seconds, ms)
#
custom_time = None  # datetime(2018, 4, 29, 0, 38, 42, 00000)


#
# Choose what TLE to use (only one key at a time can have '1' value).
#   'internet': use pyobject to retrieve TLE
#   'internet_by_url': use custom URL to download TLE
#   'local': specify list for TLE by yourself
# Set all keys to '0' values to skip ISS plotting
#
tle_type = { 'internet': 1,
             'internet_by_url': 0,
             'local': 0 }
tle_url = 'https://www.celestrak.com/NORAD/elements/stations.txt'
tle_local = [ '1 25544U 98067A   18116.84910257  .00002058  00000-0  38261-4 0  9992',
              '2 25544  51.6422 274.6376 0002768  23.1132 122.6984 15.54167281110534' ]



#
# get necessary data
#
# read stars table database from the file to Table object
if os.path.isfile('stars.html'):
    stars = astropy.table.Table.read('stars.html')
else:
    print("'stars.html' not found. Run 'dl_constellations.py'")
    sys.exit()

# generate constellations forms if there is no one
if not os.path.isfile('.constellations.html'):
    constellations_forms = extract_forms(constellations_dict, stars)
    constellations_forms.write('.constellations.html', format='ascii.html', overwrite=True)
else:
    constellations_forms = astropy.table.Table.read('.constellations.html')

# exctract each constellation
constellations = extract_constellations(stars)


#
# make plot field
#
fig,ax = prepare_skymap(fontsize=fonts['skymap'])


#
# plotting constellations
#
for constellation,constellation_form in zip(constellations,constellations_forms['idxs']):
    constellation_form = [int(s) for s in constellation_form.split()]

    # put stars on plot and connect them by pairs
    for i,iplus1 in zip(constellation_form[:-1], constellation_form[1:]):
        ax.plot( np.radians([ constellation['RA'][i], constellation['RA'][iplus1] ]),
                 [ -constellation['dec'][i]+45, -constellation['dec'][iplus1]+45 ],
                 color='skyblue', marker=None, linewidth=2.0 )

    # put name of constellation beside its first star
    ax.text( np.radians(constellation['RA'][0]),
             -constellation['dec'][0]+45,
             "{}".format(constellation['Constellation'][0]),
             fontsize=fonts['const_name'], weight='bold' )

    # print progress
    print("\r{} plotted".format(constellation['Constellation'][0]).ljust(80), end='')
print('\rconstellations are plotted'.ljust(80))


#
# handle observation geolocation
#
location = None

if internet_geolocation:
    # request geolocation by string
    print("request geolocation via geopy...")
    from geopy.geocoders import Nominatim as geocoder
    location = geocoder().geocode(location_str, timeout=5)
else:
    location = Loc()

if location is not None:
    print("observation location: " + location.address)
    obs_loc = EarthLocation( lat=location.latitude*u.deg,
                             lon=location.longitude*u.deg,
                             height=location.altitude*u.m  )
else:
    print("NO LOCATION")  # we can't proceed without location
    sys.exit()


#
# handle observation time
#
obs_time = None

# Moscow time: 3*u.hour (or datetime.now()-datetime.utcnow())
utc_offset = (time.timezone if (time.localtime().tm_isdst==0) else time.altzone)/60/60*-1

if custom_time:
    from math import modf
    obs_time = astropy.time.Time( custom_time - timedelta( hours=modf(utc_offset)[1],
                                                           minutes=modf(utc_offset)[0] ) )
else:
    obs_time = astropy.time.Time( datetime.utcnow() )

if not obs_time:
    print("NO TIME")  # we can't proceed without time
    sys.exit()


#
# plot current field-of-view, planets, Moon, Sun
#
plot_fov(ax, obs_time, obs_loc, fontsize=fonts['fov'])
plot_moon(ax, obs_time, obs_loc)
plot_sun(ax, obs_time)
plot_solarsystem(ax, obs_time, obs_loc)


#
# plot ISS (or just skip)
#
tle = None

if 1 in tle_type.values():
    if tle_type['internet']:
        pass
    elif tle_type['internet_by_url']:
        import urllib.request
        print("request ISS' two-line elements from URL...")
        with urllib.request.urlopen(tle_url) as tle_file:
            tle = [line.decode('utf-8')[:-2] for line in tle_file][1:3]
    elif tle_type['local']:
        tle = tle_local

    plot_iss(ax, obs_time, obs_loc, tle=tle)

else:
    pass



#
# title, capture with location & time and legend
#
fig.suptitle( "Star map of the\nnorthern semisphere",
              fontsize=fonts['title'], fontname="Andale Mono",
              x=plot_size+(plot_size/7),y=(plot_size/2)+(plot_size/6) )

fig.text( plot_size+(plot_size/7), (plot_size/2)-0.025,
          "{0}, ( {1:.2f}, {2:.2f} )\n\
           {3}\n".format( location_str, obs_loc.lat, obs_loc.lon,
                          str(obs_time + (utc_offset*u.hour)) ),
          fontsize=20, fontname="Andale Mono",
          horizontalalignment='center', verticalalignment='center' )

ax.legend( labelspacing=4, fontsize=fonts['legend'], handletextpad=2, borderpad=2,
           ncol=5, bbox_to_anchor=(1+0.75, 0.425) )



plt.show()
