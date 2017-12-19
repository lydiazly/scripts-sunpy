#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Read FITS data & Plot
[Import] usr_sunpy
'''
# 2017-12-11 written by Lydia
# 2017-12-19 modified by Lydia

import matplotlib.pyplot as plt
import numpy as np

from astropy.coordinates import SkyCoord
import astropy.units as u
import sunpy.map

from copy import deepcopy
import gc, os

try:
    import usr_sunpy
except ImportError, e:
    print("ERROR: Module 'usr_sunpy' doesn't exist!")
    exit(1)
from usr_sunpy import *

path = os.path.split(os.path.abspath(__file__))[0]
fname1 = path + '/' + 'data/hmi.B_720s.20150827_052400_TAI.field.fits'
fname2 = path + '/' + 'data/hmi.B_720s.20150827_052400_TAI.inclination.fits'
fname3 = path + '/' + 'data/hmi.B_720s.20150827_052400_TAI.azimuth.fits'
fname4 = path + '/' + 'data/hmi.B_720s.20150827_052400_TAI.disambig.fits'

xmin, xmax = (300., 800.)  # deg
ymin, ymax = (-500., -100.)

#======================================================================|
# Read data

print('[Path] %s' % path)
print('Reading files...')
mapb = read_sdo(fname1)
mapi = read_sdo(fname2)
mapa = read_sdo(fname3)
mapd = read_sdo(fname4)
# Disambiguate
mapa.data[np.isfinite(mapd.data) * (mapd.data > 3)] += 180.

dtor = np.pi/180.
mapbx = deepcopy(mapb)
mapby = deepcopy(mapb)
mapbz = deepcopy(mapb)
mapbx.data[:] = mapb.data * np.sin(mapi.data * dtor) * np.cos((mapa.data + 270.) * dtor)
mapby.data[:] = mapb.data * np.sin(mapi.data * dtor) * np.sin((mapa.data + 270.) * dtor)
mapbz.data[:] = mapb.data * np.cos(mapi.data * dtor)  # ~3s

# Rotate(CCW)
ang = mapb.meta['crota2']  # ~ 180
order = 3  # Test: 3 is the best
# rotation_matrix:
# [[ cos(-ang), sin(-ang)],
#  [-sin(-ang), cos(-ang)]]
mapbx = mapbx.rotate(order=order)
mapby = mapby.rotate(order=order)
mapbz = mapbz.rotate(order=order)
print('Rotation angle = %f deg (CCW)' % -ang)

# Get the center ('crpix1', 'crpix2') - First pixel is number 1.
# mask = np.isfinite(mapbz.data)  # non-nan values
# jmin, jmax = (mask.nonzero()[0].min(), mask.nonzero()[0].max())
# imin, imax = (mask.nonzero()[1].min(), mask.nonzero()[1].max())
# pcen = (0.5 * (imin + imax) * u.pix, 0.5 * (jmin + jmax) * u.pix)
pcenter = ((mapbz.meta['crpix1'] - 1) * u.pix, (mapbz.meta['crpix2'] - 1) * u.pix)
center = mapbz.pixel_to_world(*pcenter)
print('[Image_center] (%.3f, %.3f)pixel = (%7.4f, %7.4f)arcsec  (lon, lat) = (%8.5f, %8.5f)deg' %
      ((mapbz.dimensions.x.value-1.)/2., (mapbz.dimensions.y.value-1.)/2.,
         mapbz.center.Tx.value, mapbz.center.Ty.value,
         mapbz.center.heliographic_stonyhurst.lon.value, mapbz.center.heliographic_stonyhurst.lat.value))
print('[ Disk_center] (%.3f, %.3f)pixel = (%7.4f, %7.4f)arcsec  (lon, lat) = (%8.5f, %8.5f)deg' %
      (pcenter[0].value, pcenter[1].value, center.Tx.value, center.Ty.value,
       center.heliographic_stonyhurst.lon.value, center.heliographic_stonyhurst.lat.value))
print('[ Observation] (lon, lat, radius) = (%g deg, %g deg, %g m)' %
      (mapbz.heliographic_longitude.value, mapbz.heliographic_latitude.value,
       mapbz.observer_coordinate.radius.value))

# Release
del mapb, mapi, mapa, mapd; gc.collect()

# Submap
bl = SkyCoord(xmin*u.arcsec, ymin*u.arcsec, frame=mapbz.coordinate_frame)
tr = SkyCoord(xmax*u.arcsec, ymax*u.arcsec, frame=mapbz.coordinate_frame)
smapbx = mapbx.submap(bl, tr)
smapby = mapby.submap(bl, tr)
smapbz = mapbz.submap(bl, tr)
print('Submap: %s = %s arcsec' %
      (map(int, u.Quantity(smapbz.dimensions).value), [[xmin, xmax], [ymin, ymax]]))

#======================================================================|
# Plot

fig1 = plt.figure(figsize=(8, 6), dpi=100)
ax1 = fig1.add_subplot(111, projection=mapbz)
plot_map(ax1, mapbz)

# Properties
mapbz.draw_grid(axes=ax1, grid_spacing=20*u.deg, color='w', linestyle=':')
mapbz.draw_rectangle(bl, (xmax-xmin)*u.arcsec, (ymax-ymin)*u.arcsec, axes=ax1, color='yellow', linewidth=1.5)
# ax1.set_title(mapbz.latex_name, y=1.05);
plt.clim(-2000., 2000.)
fig1.savefig(path+'/'+'plothmi.png', dpi=200)

#----------------------------------------------------------------------|
iskip, jskip = (12, 12)

# Threshold for vectors
bzsmall = 100.
smapbx.data[abs(smapbz.data) < bzsmall] = 0.
smapby.data[abs(smapbz.data) < bzsmall] = 0.

fig2 = plt.figure(figsize=(9, 6), dpi=100)
ax2 = fig2.add_subplot(111, projection=smapbz)
im2 = plot_map(ax2, smapbz)
plot_vmap(ax2, smapbx, smapby, smapbz, iskip=iskip, jskip=jskip, limit=500., cmap='binary',
          scale_units='xy', scale=1/0.04, minlength=0.02)

# Properties
smapbz.draw_grid(axes=ax2, grid_spacing=10*u.deg, color='yellow', linestyle=':')
ax2.set_title(mapbz.latex_name, y=1.09)
plt.subplots_adjust(right=0.8)  # Reduce the value to shift the colorbar right
im2.set_clim(-2000., 2000.)

fig2.savefig(path+'/'+'plothmi_sub.png', dpi=200)
#----------------------------------------------------------------------|
plt.show()
