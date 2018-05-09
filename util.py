#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np

import matplotlib.pyplot as plt

from astropy.coordinates import SkyCoord, get_moon, get_sun, get_body
from astropy.coordinates.name_resolve import NameResolveError
import astropy.table

from collections import OrderedDict



#
# settings section
#
fonts = { 'skymap': 10,  # fonts sizes
          'const_name': 10,
          'fov': 10,
          'legend': 14,
          'title': 48 }
markers = { 'moon': 50,  # markers sizes
            'sun': 50,
            'iss': 20,
            'planet': 20 }
style = 'light'  # 'light' or 'dark'
plot_size = 4
dpi = 120



if style == 'light':
    colors = { 'fov_outer': 'orange' }
    plt.style.use('fivethirtyeight')
elif style == 'dark':
    colors = { 'fov_outer': 'lightseagreen' }
    plt.style.use('dark_background')



def download_stars(constellations_dict):
    """
    Form astropy.Table with constellations, their stars and (RA,dec) coordinates

    input:
        dictionary in format {'constellation_3-letter_code': [list_of_stars'_names]}

    returns:
        Table instance
    """

    print('downloading {} constellations...'.format(len(constellations_dict.keys())))

    stars_table = astropy.table.Table( names=['Constellation', 'Star', 'RA', 'dec'],
                                       meta={'name': "constellations' stars"},
                                       dtype=['object', 'object', 'float', 'float'] )

    # fill this table with data
    skipped_stars = 0
    for name,stars in constellations_dict.items():
        unique_stars = list(set(stars))

        for letter in unique_stars:
            search_request = letter + ' ' + name
            print("\r{}".format(search_request.ljust(80)), end='')
            try:
                star = SkyCoord.from_name(search_request)
                stars_table.add_row([name, letter, star.ra, star.dec])
            except NameResolveError:
                print("\rWarning: {} not found!".format(search_request).ljust(80))
                skipped_stars += 1
                continue

    print( "\r{} stars were downloaded, {} were skipped"
           .format(len(stars_table['Star']), skipped_stars).ljust(80) )
    if skipped_stars > 0:
        print("You will not be able to plot constellations. Try to re-download")

    return stars_table



def extract_constellations(stars):
    """
    Extract individual constellations from a one table and put them into separate
    tables stored in a list (sorted by constellations' names)

    returns:
        list with constellations' astropy.Table tables
    """

    constellations = []

    # find unique constellations in db
    for i,constellation in enumerate(sorted(np.unique(stars['Constellation']))):
        constellations.append( astropy.table.Table( names=['Constellation', 'Star', 'RA', 'dec'],
                                                    dtype=['object', 'object', 'float', 'float'] ) )
        for star in stars:
            if star['Constellation'] == constellation:
                constellations[i].add_row(star)

    print("extracted {} stars from {} constellations".format(len(stars), len(constellations)))

    return constellations



def extract_forms(constellations_dict, stars_table):
    """
    Define indexes of stars in the database table that corresponds to stars in constellations' forms

    returns:
        astropy.Table with 3 columns: constellation name,
                                      path (with stars' letters) (for human-reading),
                                      path (with stars' indexes) (for machine-reading)
    """

    constellations_forms = astropy.table.Table( names=['constellation', 'path', 'idxs'],
                                                dtype=['object', 'object', 'object'] )

    for (name,path),constellation in zip( constellations_dict.items(),
                                          extract_constellations(stars_table) ):
        idxs = ''
        for letter in path:
            idxs = idxs + str(list(constellation['Star']).index(letter)) + ' '

        path_for_table = ''
        for elem in path:
            path_for_table = path_for_table + elem + ' '

        constellations_forms.add_row([name, path_for_table, idxs])

    return constellations_forms



def prepare_skymap(fontsize=10):
    """
    Form outer circle for ICRS coordinate system

    returns:
        matplotlib Figure and Axes instances
    """

    fig = plt.figure(dpi=dpi)
    fig.canvas.set_window_title('starmap-northern')

    # docs quote: add an axes at position [left, bottom, width, height] where
    # all quantities are in fractions of figure width and height
    ax = fig.add_axes([0, 0, plot_size, plot_size], polar=True)

    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)  # anti-clockwise

    ax.set_ylim(-45, 90)
    ax.set_yticks(np.arange(-45, 90+0.1, 15))
    ax.set_yticklabels(ax.get_yticks()[::-1])  # reverse axis

    ax.grid(True)
    gridlines = ax.get_xgridlines() + ax.get_ygridlines()
    for line in gridlines:
        line.set_linestyle(':')
        line.set_linewidth(0.5)

    for yticklabel in ax.get_yticklabels():
        yticklabel.set_fontsize(fontsize)
    for xticklabel in ax.get_xticklabels():
        xticklabel.set_fontsize(1.75*fontsize)

    ax.tick_params(axis='x', which='major', pad=20)

    ax.get_yticklabels()[0].set_visible(False)
    ax.get_yticklabels()[-1].set_visible(False)

    fig.text( plot_size+(plot_size/2.7), 0,
              '.', color=fig.get_facecolor() )  # for additional space at right

    print("sky map is prepared")

    return fig,ax



def plot_moon(ax, obs_time, obs_loc):
    """
    Put Moon on a given Axes instance
    """

    moon = get_moon(obs_time, location=obs_loc)
    moon = SkyCoord(moon.ra, moon.dec, frame='gcrs').transform_to('icrs')

    ax.plot( [moon.ra.radian], [-moon.dec.value+45], label='Moon', linestyle='',
             color='indianred', marker='$☽$', markersize=markers['moon'] )

    print("Moon is plotted")



def plot_sun(ax, obs_time):
    """
    Put Sun on a given Axes instance
    """

    sun = get_sun(obs_time)
    sun = SkyCoord(sun.ra, sun.dec, frame='gcrs').transform_to('icrs')

    ax.plot( [sun.ra.radian], [-sun.dec.value+45], label='Sun', linestyle='',
             color='yellow', marker='$☀︎$', markersize=markers['sun'] )

    print("Sun is plotted")



def plot_iss(ax, obs_time, obs_loc, tle=None):
    """
    Put International Space Station on the plot. The function uses TLE list if
    it was given or tries to retrieve it via PyOrbital package
    """

    if tle is None:
        from pyorbital.orbital import Orbital
        print("request ISS' two-line elements via PyOrbital...")
        orb = Orbital("ISS (ZARYA)")
        iss = orb.get_observer_look( obs_time.value,
                                     obs_loc.lon.value,
                                     obs_loc.lat.value,
                                     obs_loc.height.value )
        iss_coord = SkyCoord( iss[0], iss[1], unit='deg', frame='altaz',
                              obstime=obs_time, location=obs_loc )
    else:
        import ephem
        iss = ephem.readtle('ISS', tle[0], tle[1])
        location = ephem.Observer()
        location.lat = str(obs_loc.lat.value)
        location.lon = str(obs_loc.lon.value)
        location.date = obs_time.value
        iss.compute(location)
        iss_coord = SkyCoord( iss.az, iss.alt, unit='rad', frame='altaz',
                              obstime=obs_time, location=obs_loc )

    iss_coord = iss_coord.transform_to('icrs')

    ax.plot( [iss_coord.ra.radian], [-iss_coord.dec.value+45], label='ISS', linestyle='',
             color='green', marker='$⋈$', markersize=markers['iss'] )

    print("ISS is plotted")



def plot_solarsystem(ax, obs_time, obs_loc):
    """
    Put planets of Solar System on a given Axes instance. Designation of each planet
    is its roman symbol
    """

    # for Uranus must be another symbol, but at Unicode U+2645, which renders as ♅ (Wiki)
    planets = OrderedDict([ ('Mercury', '☿'),
                            ('Venus', '♀'),
                            ('Mars', '♂'),
                            ('Jupiter', '♃'),
                            ('Saturn', '♄'),
                            ('Uranus', '♅'),
                            ('Neptune', '♆') ])

    planets_coords = [get_body(planet, obs_time, location=obs_loc) for planet in planets.keys()]
    for coords,(name,symbol) in zip(planets_coords,planets.items()):
        planet = SkyCoord(coords.ra, coords.dec, frame='gcrs').transform_to('icrs')

        ax.plot( [planet.ra.radian], [-planet.dec.value+45], label=name, linestyle='',
                 color='violet', marker='$'+symbol+'$', markersize=markers['planet'] )

    print("planets are plotted")



def plot_fov(ax, obs_time, obs_loc, fontsize=10):
    """
    Plot current field-of-view (FOV)
    """

    # coordinates of field-of-view-circle' points in (alt,az) frame
    fov_az = np.arange(0, 360+0.1, 1)
    fov_alt = np.zeros(len(fov_az))
    fov = SkyCoord(fov_az, fov_alt, unit='deg', frame='altaz', obstime=obs_time, location=obs_loc)
    # converting this coordinates to (RA,dec) format for plotting them onto plot
    fov = fov.transform_to('icrs')
    # plotting field-of-view circle
    ax.plot(fov.ra.radian, -fov.dec.value+45, '-', linewidth=0.5, color=colors['fov_outer'])

    # fill the area that we cannot observe now
    shared_ax = fov.ra.radian
    fov_circle = -fov.dec.value+45
    outer_circle = len(fov_circle) * [-(-45)+45]
    ax.fill_between( shared_ax, fov_circle, outer_circle, where=outer_circle>=fov_circle,
                     facecolor=colors['fov_outer'], alpha=0.25 )

    #
    # putting on plot ticks of circle axis (axis of azimuth) in the same way as outer circle
    #
    fov_ticks_az = np.arange(0, 345+0.1, 15)
    fov_ticks_alt = np.zeros(len(fov_ticks_az))
    fov_ticks = SkyCoord( fov_ticks_az, fov_ticks_alt, unit='deg', frame='altaz',
                          obstime=obs_time, location=obs_loc )
    fov_ticks = fov_ticks.transform_to('icrs')
    cardinal_directions = ['N', 'E', 'S', 'W']  # anti-clockwise
    cnt = 0
    for (tick_coord, label) in zip(fov_ticks, fov_ticks_az):
        ax.plot( [tick_coord.ra.radian], [-tick_coord.dec.value+45],
                 '.', color=colors['fov_outer'] )
        if int(label) % 90 == 0:
            ax.text( tick_coord.ra.radian, -tick_coord.dec.value+45,
                     cardinal_directions[cnt], fontsize=fontsize*2,
                     fontname="Apple Chancery", fontweight='bold' )
            cnt += 1
        else:
            ax.text( tick_coord.ra.radian, -tick_coord.dec.value+45,
                     "{}°".format(int(label)), fontsize=fontsize )

    #
    # plot straight axis - from South to North - of field of view circle (similarly)
    #
    SN_ax_alt = [0, 0]
    SN_ax_az = [0, 180]
    SN_ax = SkyCoord( SN_ax_az, SN_ax_alt, unit='deg', frame='altaz',
                      obstime=obs_time, location=obs_loc )
    SN_ax = SN_ax.transform_to('icrs')
    ax.plot(SN_ax.ra.radian, -SN_ax.dec.value+45, '-', linewidth=0.5, color=colors['fov_outer'])

    SN_ticks_alt = [alt for alt in np.arange(75, 15-0.1, -15)]
    SN_ticks_az = len(SN_ticks_alt)*[0]
    SN_ticks_alt = SN_ticks_alt + SN_ticks_alt[::-1]
    SN_ticks_az = SN_ticks_az + len(SN_ticks_az)*[180]
    SN_ticks = SkyCoord( SN_ticks_az, SN_ticks_alt, unit='deg', frame='altaz',
                         obstime=obs_time, location=obs_loc )
    SN_ticks = SN_ticks.transform_to('icrs')
    for (tick_coord, label) in zip(SN_ticks, SN_ticks_alt):
        ax.plot( [tick_coord.ra.radian], [-tick_coord.dec.value+45],
                 '.', color=colors['fov_outer'] )
        ax.text( tick_coord.ra.radian, -tick_coord.dec.value+45,
                 "{}°".format(int(label)), fontsize=fontsize )

    #
    # plot curved axis - from West to East - of field of view circle (similarly)
    #
    WE_ax_alt = [alt for alt in np.arange(0, 90+0.1, 1)]
    WE_ax_az = len(WE_ax_alt)*[90] + (len(WE_ax_alt)-1)*[270]
    WE_ax_alt = WE_ax_alt + WE_ax_alt[:-1][::-1]
    WE_ax = SkyCoord (WE_ax_az, WE_ax_alt, unit='deg', frame='altaz',
                     obstime=obs_time, location=obs_loc )
    WE_ax = WE_ax.transform_to('icrs')
    ax.plot(WE_ax.ra.radian, -WE_ax.dec.value+45, '-', linewidth=0.5, color=colors['fov_outer'])

    WE_ticks_alt = [alt for alt in np.arange(15, 75+0.1, 15)]
    WE_ticks_az = len(WE_ticks_alt)*[90] + len(WE_ticks_alt)*[270]
    WE_ticks_alt = WE_ticks_alt + WE_ticks_alt[::-1]
    WE_ticks = SkyCoord( WE_ticks_az, WE_ticks_alt, unit='deg', frame='altaz',
                         obstime=obs_time, location=obs_loc )
    WE_ticks = WE_ticks.transform_to('icrs')
    for (tick_coord, label) in zip(WE_ticks, WE_ticks_alt):
        ax.plot( [tick_coord.ra.radian], [-tick_coord.dec.value+45],
                 '.', color=colors['fov_outer'] )
        ax.text( tick_coord.ra.radian, -tick_coord.dec.value+45,
                 "{}°".format(int(label)), fontsize=fontsize )

    print("field-of-view is plotted")



if __name__ == '__main__':
    prepare_skymap()
