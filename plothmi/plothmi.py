#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Read FITS data & Plot
[Import] numpy, matplotlib, astropy, sunpy, usr_sunpy
[Example Data] https://pan.baidu.com/s/1nwsIcDr (pswd: s5re)
'''
# 2017-12-11 written by Lydia
# 2018-04-01 modified by Lydia

from __future__ import division, print_function

# import matplotlib
# matplotlib.use('Qt5Agg')  # 'Qt5Agg', 'TkAgg', 'Agg', ...

import matplotlib.pyplot as plt
import numpy as np

from astropy.coordinates import SkyCoord
import astropy.units as u
import sunpy.map

from copy import deepcopy
import gc, os, time
gc.disable()

# [usr_sunpy]
# Funcions: read_sdo, plot_map, plot_vmap, image_to_helio, ...
path = os.path.split(os.path.abspath(__file__))[0]
import sys
sys.path.append(path + '/../modules')
from usr_sunpy import *

#======================================================================|
# Global Parameters

# Data
fname1 = path + '/' + 'data/hmi.B_720s.20150827_052400_TAI.field.fits'
fname2 = path + '/' + 'data/hmi.B_720s.20150827_052400_TAI.inclination.fits'
fname3 = path + '/' + 'data/hmi.B_720s.20150827_052400_TAI.azimuth.fits'
fname4 = path + '/' + 'data/hmi.B_720s.20150827_052400_TAI.disambig.fits'

# Range of submap (arcsec)
xmin, xmax = (300., 800.)
ymin, ymax = (-500., -100.)

#======================================================================|
# Read data

print('[Path] %s' % path)
print('Reading data...')
mapb = read_sdo(fname1)
mapi = read_sdo(fname2)
mapa = read_sdo(fname3)
mapd = read_sdo(fname4)
# Disambiguate
# mapa.data[np.isfinite(mapd.data) & (mapd.data > 3)] += 180.
mapa.data[mapd.data > 3] += 180.

t0 = time.time()
dtor = np.pi/180.
mapbx = deepcopy(mapb)
mapby = deepcopy(mapb)
mapbz = deepcopy(mapb)
mapbx.data[:] = mapb.data * np.sin(mapi.data * dtor) * np.cos((mapa.data + 270.) * dtor)
mapby.data[:] = mapb.data * np.sin(mapi.data * dtor) * np.sin((mapa.data + 270.) * dtor)
mapbz.data[:] = mapb.data * np.cos(mapi.data * dtor)
print('(Time of getting Bvec: %f sec)' % (time.time() - t0))

# Rotate(CCW)
# This function will remove old CROTA keywords from the header.
order = 3  # Test: 1 or 3 is ok
print('Correcting image axes...')
t0 = time.time()
with np.errstate(invalid='ignore'):  # Suppress warnings of NaNs
    mapbx = mapbx.rotate(order=order)
    mapby = mapby.rotate(order=order)
    mapbz = mapbz.rotate(order=order)
print('Rotation angle = %f deg (CCW)' % -mapb.meta['crota2'])
print('(Time of rotation: %f sec)' % (time.time() - t0))

# Check the center ('crpix1', 'crpix2') - First pixel is number 1.
# mask = np.isfinite(mapbz.data)  # Get the center by finite values
# jmin, jmax = (mask.nonzero()[0].min(), mask.nonzero()[0].max())
# imin, imax = (mask.nonzero()[1].min(), mask.nonzero()[1].max())
# pcenter = (0.5 * (imin + imax) * u.pix, 0.5 * (jmin + jmax) * u.pix)
pcenter = ((mapbz.meta['crpix1'] - 1) * u.pix, (mapbz.meta['crpix2'] - 1) * u.pix)
center = mapbz.pixel_to_world(*pcenter)
print('[Image_center]\n\t(%.3f, %.3f) pixel = (%7.4f, %7.4f) arcsec\n\t(lon, lat) = (%8.5f, %8.5f) deg' %
      ((mapbz.dimensions.x.value-1.)/2., (mapbz.dimensions.y.value-1.)/2.,
        mapbz.center.Tx.value, mapbz.center.Ty.value,
        mapbz.center.heliographic_stonyhurst.lon.value,
        mapbz.center.heliographic_stonyhurst.lat.value))
print('[Disk_center]\n\t(%.3f, %.3f) pixel = (%7.4f, %7.4f) arcsec\n\t(lon, lat) = (%8.5f, %8.5f) deg' %
      (pcenter[0].value, pcenter[1].value,
       center.Tx.value, center.Ty.value,
       center.heliographic_stonyhurst.lon.value,
       center.heliographic_stonyhurst.lat.value))
print('[Observation]\n\t(lon, lat, radius) = (%g deg, %g deg, %g m)' %
      (mapbz.heliographic_longitude.value,
       mapbz.heliographic_latitude.value,
       mapbz.observer_coordinate.radius.value))

# Submap
bl = SkyCoord(xmin*u.arcsec, ymin*u.arcsec, frame=mapbz.coordinate_frame)
tr = SkyCoord(xmax*u.arcsec, ymax*u.arcsec, frame=mapbz.coordinate_frame)
smapbx = mapbx.submap(bl, tr)
smapby = mapby.submap(bl, tr)
smapbz = mapbz.submap(bl, tr)
print('\nSubmap: %s = %s arcsec' %
      (list(map(int, u.Quantity(smapbz.dimensions).value)), [[xmin, xmax], [ymin, ymax]]))

#======================================================================|
# Plot

fig1 = plt.figure(figsize=(8, 6), dpi=100)
ax1 = fig1.add_subplot(111, projection=mapbz)
plot_map(ax1, mapbz)

# Properties
mapbz.draw_grid(axes=ax1, grid_spacing=20*u.deg, color='w', linestyle=':')
# mapbz.draw_limb(axes=ax1, color='b', linewidth=1.5)
mapbz.draw_rectangle(bl, (xmax-xmin)*u.arcsec, (ymax-ymin)*u.arcsec,
                     axes=ax1, color='yellow', linewidth=1.5)
# ax1.set_title(mapbz.latex_name, y=1.05)
plt.clim(-2000., 2000.)
fig1.savefig(path+'/'+'plothmi_disk.png', dpi=200)

#----------------------------------------------------------------------|
iskip, jskip = (12, 12)

fig2 = plt.figure(figsize=(9, 6), dpi=100)
ax2 = fig2.add_subplot(111, projection=smapbz)
im2 = plot_map(ax2, smapbz)
plot_vmap(ax2, smapbx, smapby, smapbz, iskip=iskip, jskip=jskip, cmin=100., vmax=500., cmap='binary',
          scale_units='xy', scale=1/0.05, minlength=0.02)

# Properties
smapbz.draw_grid(axes=ax2, grid_spacing=10*u.deg, color='yellow', linestyle=':')
ax2.set_title(mapbz.latex_name+' (submap)', y=1.1)
plt.subplots_adjust(right=0.8)  # Reduce the value to move the colorbar to the right
im2.set_clim(-2000., 2000.)

fig2.savefig(path+'/'+'plothmi_sub.png', dpi=200)
#----------------------------------------------------------------------|
plt.show()
