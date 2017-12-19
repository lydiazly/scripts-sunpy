#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Funtions for sunpy.
[Import] matplotlib, numpy, sunpy
[Reference] http://docs.sunpy.org/en/stable/code_ref/map.html
'''
# 2017-12-11 written by Lydia
# 2017-12-19 modified by Lydia

import astropy.units as u
from copy import deepcopy
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
import sunpy.map
from sunpy.visualization import axis_labels_from_ctype
import os

#======================================================================|
def read_sdo(filename):
    '''
    -------------------------------------------------------------------
    Read from a FITS file.
    
    [Return] A Map object
    
    [Properties]
    * data: ndarray - A 2d list or ndarray containingthe map data.
                          data[i,j]: j from the bottom, i from te left.
    * meta: dict - A dictionary of the original image headr tags.
    * dimensions: PixelPair - x axis first, y axis second.
                              Use 'u.Quantity(mapbz.dimensions).value'
                              to get the numerical values.
    [Notes]
    A number of the properties of this class are returned as two-value named
    tuples that can either be indexed by position ([0] or [1]) or be
    accessed by the names (.x and .y) or (.axis1 and .axis2). Things that
    refer to pixel axes use the .x, .y convention, where x and y refer to
    the FITS axes (x for columns y for rows). Spatial axes use .axis1 and
    .axis2 which correspond to the first and second axes in the header.
    axis1 corresponds to the coordinate axis for x and axis2 corresponds to y.
    
    [Reference]
    http://docs.sunpy.org/en/stable/code_ref/map.html#using-map-objects
    --------------------------------------------------------------------
    '''
    smap = sunpy.map.Map(filename)
    print('%s\t%s' %
          (os.path.basename(filename),
           map(int, u.Quantity(smap.dimensions).value)))
    return smap

#======================================================================|
def plot_map(ax, smap, grid=False, **kwargs):
    '''
    --------------------------------------------------------------------
    Plot the map object using matplotlib, in a method equivalent to imshow()
    using nearest neighbour interpolation.
    
    [Return] matplotlib image
    
    [Parameters]
    * smap: Map
    kwargs:
    - plot_kwargs(sunpy): annotate, axes, title
    - imshow_kwargs(matplotlib): kwargs
    
    [Reference]
    http://docs.sunpy.org/en/stable/code_ref/map.html#sunpy.map.mapbase.GenericMap.plot
    --------------------------------------------------------------------
    '''
    im = smap.plot(annotate=True, **kwargs)
    ax.set_aspect(1)
    ax.grid(grid)
    ax.set_autoscale_on(False)  # Disable autoscaling from other plots
    plt.subplots_adjust(left=0.1, bottom=0.11, right=0.9, top=0.88)  # x0, y0, x1, y1 (default = [0.125, 0.11, 0.9, 0.88])
    cax = plt.axes([0.9-0.015, 0.11, 0.03*0.8, 0.77])  # [x=(x1+pad), y=y0, w=ratio*(x1-x0), h=(y1-y0)]
    plt.colorbar(im, cax=cax)
    return im

#======================================================================|
def plot_map_p(ax, X, Y, smap, grid=False, cmap='gray', **kwargs):
    '''
    --------------------------------------------------------------------
    Plot the map object using matplotlib, in a method equivalent to imshow()
    using nearest neighbour interpolation.
    
    [Return] matplotlib image
    
    [Parameters]
    * X, Y: 2D ndarrays
    * smap: Map
    kwargs:
    - plot_kwargs(sunpy): annotate, axes, title
    - pcolormesh_kwargs(matplotlib): kwargs
    
    [Reference]
    http://docs.sunpy.org/en/stable/code_ref/map.html#sunpy.map.mapbase.GenericMap.plot
    --------------------------------------------------------------------
    '''
    im = ax.pcolormesh(X, Y, smap.data.T, cmap=cmap, **kwargs)
    ax.set_aspect(1)
    ax.grid(grid, ls=':')
    ax.set_autoscale_on(False)  # Disable autoscaling from other plots
    ax.set_title(smap.latex_name, y=1.02)
    ax.set_xlabel('X (arcsec)')
    ax.set_ylabel('Y (arcsec)')
    plt.subplots_adjust(left=0.1, bottom=0.11, right=0.9, top=0.88)  # x0, y0, x1, y1 (default = [0.125, 0.11, 0.9, 0.88])
    cax = plt.axes([0.9-0.025, 0.11, 0.03*0.8, 0.77])  # [x=(x1+pad), y=y0, w=ratio*(x1-x0), h=(y1-y0)]
    plt.colorbar(im, cax=cax)
    return im

#======================================================================|
def plot_vmap(ax, mapu, mapv, mapc,
              iskip=10, jskip=10, limit=1000., cmap='binary',
              scale_units='xy', scale=1/0.05, minlength=0.05, width=0.003,
              headlength=6, headwidth=5, headaxislength=3,
              **kwargs):
    '''
    --------------------------------------------------------------------
    Vector plot using quiver in matplotlib.
    
    [Return] matplotlib image
    
    [Parameters]
    * mapu: Map of Vx
    * mapv: Map of Vy
    * mapc: Map for setting color
            (reset mapc.data first to change to color function)
    * limit: max length of (Vx, Vy), same unit as Vx
    * kwargs: quiver_kwargs(matplotlib)
    
    [Reference]
    https://matplotlib.org/devdocs/api/_as_gen/matplotlib.axes.Axes.quiver.html#matplotlib-axes-axes-quiver
    --------------------------------------------------------------------
    '''
    iskip = int(iskip); jskip = int(jskip)
    dimy, dimx = mapc.data.shape
    
    pixmin = (0., 0.)  # pixel
    pixmax = u.Quantity(mapc.dimensions).value
    x0 = (mapc.pixel_to_world(pixmin[0]*u.pix,pixmin[1]*u.pix).Tx).to(u.deg).value  # deg
    x1 = (mapc.pixel_to_world(pixmax[0]*u.pix,pixmax[1]*u.pix).Tx).to(u.deg).value
    y0 = (mapc.pixel_to_world(pixmin[0]*u.pix,pixmin[1]*u.pix).Ty).to(u.deg).value
    y1 = (mapc.pixel_to_world(pixmax[0]*u.pix,pixmax[1]*u.pix).Ty).to(u.deg).value
    # Resample
#    rmapu = mapu.resample(u.Quantity([nx, ny], u.pix)).data.T
#    rmapv = mapv.resample(u.Quantity([nx, ny], u.pix)).data.T
#    rmapc = np.sign(mapc.resample(u.Quantity([nx, ny], u.pix)).data.T)
    rmapu = mapu.data.T[::iskip, ::jskip]
    rmapv = mapv.data.T[::iskip, ::jskip]
    rmapc = np.sign(mapc.data.T)[::iskip, ::jskip]
    # Clip
    mag = np.sqrt(rmapu**2 + rmapv**2)
    mask = np.where(mag > limit)
    rmapu[mask] = rmapu[mask] * (limit/mag[mask])
    rmapv[mask] = rmapv[mask] * (limit/mag[mask])

    X, Y = np.mgrid[x0:x1:dimx*1j, y0:y1:dimy*1j][:, ::iskip, ::jskip]  # deg
    im = ax.quiver(X, Y, rmapu, rmapv,
                   rmapc, cmap=cmap,
                   angles='uv', pivot='tail',
                   scale_units=scale_units, scale=scale,
                   minlength=minlength, width=width,
                   headlength=headlength, headwidth=headwidth, headaxislength=headaxislength,
                   transform=ax.get_transform('world'),
                   **kwargs)
    im.set_clim(-1, 1)
    return im

#======================================================================|
def plot_vmap_p(ax, X, Y, mapu, mapv, mapc,
                iskip=10, jskip=10, limit=1000., cmap='binary',
                scale_units='xy', scale=1/0.05, minlength=0.05, width=0.003,
                headlength=6, headwidth=5, headaxislength=3,
                **kwargs):
    '''
    --------------------------------------------------------------------
    Vector plot using quiver in matplotlib(for projected maps).
    
    [Return] matplotlib image
    
    [Parameters]
    * X, Y: 2D ndarrays
    * mapu: map of Vx
    * mapv: map of Vy
    * mapc: map for setting color
            (reset mapc.data first to change to color function)
    * limit: max length of (Vx, Vy), same unit as Vx
    * kwargs: quiver_kwargs(matplotlib)
    
    [Reference]
    https://matplotlib.org/devdocs/api/_as_gen/matplotlib.axes.Axes.quiver.html#matplotlib-axes-axes-quiver
    --------------------------------------------------------------------
    '''
    iskip = int(iskip); jskip = int(jskip)
    
    P = 0.  # The angle of the northern extremity, CCW from the north point of the disk.
    L0 = np.deg2rad(mapc.heliographic_longitude.value)  # The longitude of the center of the disk.
    B0 = np.deg2rad(mapc.heliographic_latitude.value)  # The latitude of the center of the disk.
    dimy, dimx = mapc.data.shape
    dx, dy = (mapc.meta['cdelt1'], mapc.meta['cdelt2'])  # arcsec/pix
    x0 = -dx * (dimx-1)/2.  # Set x=0, y=0 at the center of the image.
    y0 = -dy * (dimy-1)/2.
    x1 = dx * (dimx-1)/2.
    y1 = dy * (dimy-1)/2.
    tmp = mapc.pixel_to_world(*((dimx-1)/2., (dimy-1)/2.)*u.pix).transform_to('heliographic_stonyhurst')
    Lc = np.deg2rad(tmp.lon.value)  # The longitude of the center of the image.
    Bc = np.deg2rad(tmp.lat.value)  # The latitude of the center of the image.
    
    # Resample
    fmapu = interpolate.RectBivariateSpline(np.linspace(x0, x1, dimx), np.linspace(y0, y1, dimy), mapu.data.T)
    fmapv = interpolate.RectBivariateSpline(np.linspace(x0, x1, dimx), np.linspace(y0, y1, dimy), mapv.data.T)
    fmapc = interpolate.RectBivariateSpline(np.linspace(x0, x1, dimx), np.linspace(y0, y1, dimy), mapc.data.T)
    ax1, ax2, ay1, ay2 = proj_matrix(P, L0, B0, Bc, Lc, 2)
    hx = np.linspace(X.min(), X.max(), dimx)
    hy = np.linspace(Y.min(), Y.max(), dimy)
    hx, hy = np.mgrid[X.min():X.max():dimx*1j, Y.min():Y.max():dimy*1j]
    ix = ax1 * hx + ay1 * hy
    iy = ax2 * hx + ay2 * hy
    
    rmapu = np.array(map(fmapu, ix.flatten(), iy.flatten())).reshape((dimx, dimy))
    rmapv = np.array(map(fmapv, ix.flatten(), iy.flatten())).reshape((dimx, dimy))
    rmapc = np.array(map(fmapc, ix.flatten(), iy.flatten())).reshape((dimx, dimy))
    rmapu = rmapu[::iskip, ::jskip]
    rmapv = rmapv[::iskip, ::jskip]
    rmapc = np.sign(rmapc)[::iskip, ::jskip]
    
    # Clip
    mag = np.sqrt(rmapu**2 + rmapv**2)
    mask = np.where(mag > limit)
    rmapu[mask] = rmapu[mask] * (limit/mag[mask])
    rmapv[mask] = rmapv[mask] * (limit/mag[mask])

    im = ax.quiver(hx[::iskip, ::jskip], hy[::iskip, ::jskip], rmapu, rmapv,
                   rmapc, cmap=cmap,
                   angles='uv', pivot='tail',
                   scale_units=scale_units, scale=scale,
                   minlength=minlength, width=width,
                   headlength=headlength, headwidth=headwidth, headaxislength=headaxislength,
                   **kwargs)
    im.set_clim(-1, 1)
    return im

#======================================================================|
def image_to_helio(*smap):
    '''
    -------------------------------------------------------------------
    Vector: (U, V, W)_helio = A33 dot (U, V, W)_image
    Coords: (x, y)_helio = A22.T.I dot (x, y)_image
    
    [Return] list
    smap -> (x_h, y_h)
    smapx, smapy, smapz -> (smapx_h, smapy_h, smapz_h)
    
    [Reference]
    http://link.springer.com/10.1007/BF00158295
    --------------------------------------------------------------------
    '''
    P = 0.  # The angle of the northern extremity, CCW from the north point of the disk.
    L0 = np.deg2rad(smap[-1].heliographic_longitude.value)  # The longitude of the center of the disk.
    B0 = np.deg2rad(smap[-1].heliographic_latitude.value)  # The latitude of the center of the disk.
    dimy, dimx = smap[-1].data.shape
    dx, dy = (smap[-1].meta['cdelt1'], smap[-1].meta['cdelt2'])  # arcsec/pix
    x0 = -dx * (dimx-1)/2.  # Set x=0, y=0 at the center of the image.
    y0 = -dy * (dimy-1)/2.
    tmp = smap[-1].pixel_to_world(*((dimx-1)/2., (dimy-1)/2.)*u.pix).transform_to('heliographic_stonyhurst')
    Lc = np.deg2rad(tmp.lon.value)  # The longitude of the center of the image.
    Bc = np.deg2rad(tmp.lat.value)  # The latitude of the center of the image.
    
    ix, iy = np.mgrid[0:dimx-1:dimx*1j, 0:dimy-1:dimy*1j]
    ix *= dx; ix += x0
    iy *= dy; iy += y0
    
    ax1, ax2, ax3, ay1, ay2, ay3, az1, az2, az3 = proj_matrix(P, L0, B0, Bc, Lc)
    if len(smap) == 1:  # Image transformation
        [[cx1, cx2], [cy1, cy2]] = np.array(np.mat([[ax1, ay1], [ax2, ay2]]).I)
        hx = cx1 * ix + cx2 * iy
        hy = cy1 * ix + cy2 * iy
        return hx, hy
    elif len(smap) == 3:  # Vector transformation
        hmapx = deepcopy(smap[0])
        hmapy = deepcopy(smap[1])
        hmapz = deepcopy(smap[2])
        hmapx.data[:] = ax1 * smap[0].data + ax2 * smap[1].data + ax3 * smap[2].data
        hmapy.data[:] = ay1 * smap[0].data + ay2 * smap[1].data + ay3 * smap[2].data
        hmapz.data[:] = az1 * smap[0].data + az2 * smap[1].data + az3 * smap[2].data
        return hmapx, hmapy, hmapz
    else:
        print('Warning: The number of variables must be 1 or 3! Doing nothing.')
        return

#======================================================================|
def proj_matrix(P, L0, B0, Bc, Lc, *dim):
    '''
    Vector: (U, V, W)_helio = A33 dot (U, V, W)_image
    Coords: (x, y)_image = A22.T dot (x, y)_helio
    
    [Return] dim=2: ax1, ax2, ay1, ay2
             dim=3: ax1, ax2, ax3, ay1, ay2, ay3, az1, az2, az3
    
    [Parameters]
    *  P: The angle of the northern extremity, CCW from the north point of the disk.
    * L0: The longitude of the center of the disk.
    * B0: The latitude of the center of the disk.
    * Lc: The longitude of the the referenced point.
    * Bc: The latitude of the referenced point.
    '''
    if len(dim) == 0: dim = (3,)
    ax1 =  -np.sin(B0) * np.sin(P) * np.sin(Lc - L0) + np.cos(P) * np.cos(Lc - L0)
    ax2 =   np.sin(B0) * np.cos(P) * np.sin(Lc - L0) + np.sin(P) * np.cos(Lc - L0)
    ay1 = (-np.sin(Bc) * (np.sin(B0) * np.sin(P) * np.cos(Lc - L0) + np.cos(P) * np.sin(Lc - L0))
          - np.cos(Bc) *  np.cos(B0) * np.sin(P))
    ay2 = ( np.sin(Bc) * (np.sin(B0) * np.cos(P) * np.cos(Lc - L0) - np.sin(P) * np.sin(Lc - L0))
          + np.cos(Bc) *  np.cos(B0) * np.cos(P))
    if dim[0] == 2:
        return ax1, ax2, ay1, ay2
    else:
        ax3 =  -np.cos(B0) * np.sin(Lc - L0)
        ay3 =  -np.cos(B0) *  np.sin(Bc) * np.cos(Lc - L0) + np.sin(B0) * np.cos(Bc)
        az1 = ( np.cos(Bc) * (np.sin(B0) * np.sin(P) * np.cos(Lc - L0) + np.cos(P) * np.sin(Lc - L0))
              - np.sin(Bc) *  np.cos(B0) * np.sin(P))
        az2 = (-np.cos(Bc) * (np.sin(B0) * np.cos(P) * np.cos(Lc - L0) - np.sin(P) * np.sin(Lc - L0))
              + np.sin(Bc) *  np.cos(B0) * np.cos(P))
        az3 =   np.cos(Bc) *  np.cos(B0) * np.cos(Lc - L0) + np.sin(Bc) * np.sin(B0)
        return ax1, ax2, ax3, ay1, ay2, ay3, az1, az2, az3

