#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Read FITS data & Projection & Plot
[Packages] numpy, matplotlib, astropy, sunpy, usr_sunpy
[Example Data] https://pan.baidu.com/s/1nwsIcDr (pswd: s5re)
'''
# 2017-12-19 written by Lydia
# 2018-01-24 modified by Lydia

from __future__ import division, print_function

import matplotlib.pyplot as plt
import numpy as np

from astropy.coordinates import SkyCoord
import astropy.units as u
import sunpy.map

from copy import deepcopy
import gc, os

try:
    import usr_sunpy
except ImportError as e:
    print("Import Error:", e)
    exit(1)
from usr_sunpy import *

path = os.path.split(os.path.abspath(__file__))[0]
fname1 = path + '/' + 'data/hmi.B_720s.20150827_052400_TAI.field.fits'
fname2 = path + '/' + 'data/hmi.B_720s.20150827_052400_TAI.inclination.fits'
fname3 = path + '/' + 'data/hmi.B_720s.20150827_052400_TAI.azimuth.fits'
fname4 = path + '/' + 'data/hmi.B_720s.20150827_052400_TAI.disambig.fits'

xmin, xmax = (500.,800.)  # arcsec
ymin, ymax = (-450.,-200.)

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
order = 3  # Test: 3 is the best
print('Correcting image axes...')
with np.errstate(invalid='ignore'):  # Suppress warnings of NaNs
    mapbx = mapbx.rotate(order=order)
    mapby = mapby.rotate(order=order)
    mapbz = mapbz.rotate(order=order)
print('Rotation angle = %f deg (CCW)' % -mapb.meta['crota2'])

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
      (list(map(int, u.Quantity(smapbz.dimensions).value)), [[xmin, xmax], [ymin, ymax]]))

#======================================================================|
# Projection

hx, hy = image_to_helio(smapbz)
smapbx_h, smapby_h, smapbz_h = image_to_helio(smapbx, smapby, smapbz)

print('(xmin, xmax) = (%9.3f, %9.3f) arcsec\n(ymin, ymax) = (%9.3f, %9.3f) arcsec' %
      (hx.min(), hx.max(), hy.min(), hy.min()))

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
fig1.savefig(path+'/'+'plothmi_projection_disk.png', dpi=200)

#----------------------------------------------------------------------|
iskip, jskip = (12, 12)

fig2 = plt.figure(figsize=(9, 6), dpi=100)
ax2 = fig2.add_subplot(111, projection=smapbz)
im2 = plot_map(ax2, smapbz)
plot_vmap(ax2, smapbx, smapby, smapbz, iskip=iskip, jskip=jskip, cmin=100., vmax=500., cmap='binary',
          scale_units='xy', scale=1/0.05, minlength=0.02)

# Properties
smapbz.draw_grid(axes=ax2, grid_spacing=10*u.deg, color='yellow', linestyle=':')
ax2.set_title(mapbz.latex_name, y=1.09)
plt.subplots_adjust(right=0.8)  # Reduce the value to shift the colorbar right
im2.set_clim(-2000., 2000.)

fig2.savefig(path+'/'+'plothmi_projection_sub.png', dpi=200)

#----------------------------------------------------------------------|
iskip, jskip = (10, 10)

fig3 = plt.figure(figsize=(9, 5), dpi=100)
ax3 = fig3.add_subplot(111)
im3 = plot_map_p(ax3, hx, hy, smapbz_h)
plot_vmap_p(ax3, hx, hy, smapbx_h, smapby_h, smapbz_h, 
            iskip=iskip, jskip=jskip, cmin=100., vmax=300., cmap='binary',
            scale_units='xy', scale=1/0.03, minlength=0.05)

# Properties
#ax3.grid(True, ls=':', color='yellow', alpha=0.8)
ax3.set_title(mapbz.latex_name+' (projected)', y=1.02);
plt.subplots_adjust(right=0.9)  # Reduce the value to shift the colorbar right
ax3.set_xlim((-170,170))
ax3.set_ylim((-110,100))
im3.set_clim((-2500,2500))

fig3.savefig(path+'/'+'plothmi_projection_sub.png', dpi=200)
#----------------------------------------------------------------------|
plt.show()
