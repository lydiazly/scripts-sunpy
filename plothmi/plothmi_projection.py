#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Read FITS data & Projection & Plot
[Import] numpy, matplotlib, astropy, sunpy, usr_sunpy
[Example Data] https://pan.baidu.com/s/1nwsIcDr (pswd: s5re)
'''
# 2017-12-19 written by Lydia
# 2018-04-25 modified by Lydia

from __future__ import division, print_function

# import matplotlib
# matplotlib.use('Qt5Agg')  # 'Qt5Agg', 'TkAgg', 'Agg', ...

import matplotlib.pyplot as plt
import numpy as np

from astropy.coordinates import SkyCoord
import astropy.units as u
import sunpy.map

from copy import deepcopy
import os, time

# [usr_sunpy] Funcions: read_sdo, plot_map, plot_vmap, image_to_helio, ...
import sys
sys.path.append('../modules')
from usr_sunpy import read_sdo, plot_map, plot_vmap, image_to_helio

#======================================================================|
# Global Parameters

# Data
fname1 = 'data/hmi.B_720s.20150827_052400_TAI.field.fits'
fname2 = 'data/hmi.B_720s.20150827_052400_TAI.inclination.fits'
fname3 = 'data/hmi.B_720s.20150827_052400_TAI.azimuth.fits'
fname4 = 'data/hmi.B_720s.20150827_052400_TAI.disambig.fits'

# Range of submap (arcsec)
xmin, xmax = (500.,800.)
ymin, ymax = (-450.,-200.)

#======================================================================|
# Read data

print('[Path] %s' % os.getcwd())
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
# `rotate` function will remove old CROTA keywords.
order = 1  # Test: 1 or 3 is ok
# Suppress metadata warnings if sunpy >= 0.9.0:
mapbx.meta['hgln_obs'] = 0.; mapby.meta['hgln_obs'] = 0.; mapbz.meta['hgln_obs'] = 0.
print('Correcting image axes...')
t0 = time.time()
# Suppress warnings of NaNs:
with np.errstate(invalid='ignore'):
    mapbx = mapbx.rotate(order=order)
    mapby = mapby.rotate(order=order)
    mapbz = mapbz.rotate(order=order)
print('Rotation angle = %f deg (CCW)' % -mapb.meta['crota2'])
print('(Time of rotation: %f sec)' % (time.time() - t0))

# Check the center ('crpix1', 'crpix2') - First pixel is number 1.
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
      (tuple(map(int, u.Quantity(smapbz.dimensions).value)), ((xmin, xmax), (ymin, ymax))))

#======================================================================|
# Projection

t0 = time.time()
hx, hy = image_to_helio(smapbz)
smapbx_h, smapby_h, smapbz_h = image_to_helio(smapbx, smapby, smapbz)
print('(Time of projection: %f sec)' % (time.time() - t0))
print('\nProjected:')
print('(xmin, xmax) = (%9.3f, %9.3f) arcsec\n(ymin, ymax) = (%9.3f, %9.3f) arcsec' %
      (hx.min(), hx.max(), hy.min(), hy.min()))

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
fig1.savefig('projection_disk.png', dpi=200)

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

fig2.savefig('projection_sub.png', dpi=200)

#----------------------------------------------------------------------|
iskip, jskip = (10, 10)

fig3 = plt.figure(figsize=(9, 5), dpi=100)
ax3 = fig3.add_subplot(111)
im3 = plot_map(ax3, smapbz_h, coords=(hx, hy))
plot_vmap(ax3, smapbx_h, smapby_h, smapbz_h, coords=(hx, hy),
          iskip=iskip, jskip=jskip, cmin=100., vmax=300., cmap='binary',
          scale_units='xy', scale=1/0.03, minlength=0.05)

# Properties
ax3.grid(True, ls=':', alpha=0.8)
ax3.set_title(mapbz.latex_name+' (projected)', y=1.02);
plt.subplots_adjust(right=0.9)  # Reduce the value to shift the colorbar right
ax3.set_xlim((-170,170))
ax3.set_ylim((-110,100))
im3.set_clim((-2000,2000))

fig3.savefig('projection_sub_projected.png', dpi=200)
#----------------------------------------------------------------------|
plt.show()
