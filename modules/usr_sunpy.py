'''
User functions.
- SunPy Version: 0.9.0
- Reference http://docs.sunpy.org/en/stable/code_ref/map.html
'''
# 2017-12-11 written by Lydia
# 2018-07-01 modified by Lydia
from __future__ import division, print_function
import astropy.units as u
import numpy as np
import sunpy.map
#======================================================================|
def read_sdo(filename):
    '''
    ----------------------------------------------------------------------------
    Example of reading from a FITS file. Print the filename and dimensions.
    
    [Return] a sunpy `GenericMap` or subclass(e.g. `HMIMap`) object
    - data - a 2D numpy `ndarray`
      data[i, j]: i from the bottom(y), j from te left(x).
    - meta - a `dict` of the original image headr tags.
    
    [See also]
    - help(sunpy.map.GenericMap)
    - http://docs.sunpy.org/en/stable/code_ref/map.html#using-map-objects

    ----------------------------------------------------------------------------
    '''
    import os
    smap = sunpy.map.Map(filename)  # Read & return a `GenericMap`
    print('%s\t%s' % (os.path.basename(filename),
          list(map(int, u.Quantity(smap.dimensions).value))))
    return smap

#======================================================================|
def plot_map(ax, smap, coords=None, grid=False, cmap='gray', **kwargs):
    '''
    ----------------------------------------------------------------------------
    Plot image.
    
    [Plot Function]
    - plot_map(ax, smap, **kwargs)       -> `pcolormesh` from matplotlib
    - plot_map(ax, smap, X, Y, **kwargs) -> `imshow` from matplotlib

    [Parameters]
    - ax: a matplotlib `Axes` object
    - smap: a sunpy `GenericMap`
    - coords: two 2D numpy `ndarray`s: (X, Y)
    - grid: draw grids or not
    - cmap: name of color map
    - kwargs:
      - sunpy_kwargs: annotate, axes, title
      - matplotlib_kwargs:
        If coords is None(default), use kwargs of `imshow`,
        else use kwargs of `pcolormesh`.

    [Return] a matplotlib image object

    [See also]
    http://docs.sunpy.org/en/stable/code_ref/map.html#sunpy.map.mapbase.GenericMap.plot

    ----------------------------------------------------------------------------
    '''
    if not isinstance(smap, sunpy.map.mapbase.GenericMap):
        raise TypeError("smap should be a sunpy GenericMap.")
    if coords and (len(coords) != 2 or not any(isinstance(i, np.ndarray) for i in coords)):
        raise ValueError("coords should be a list of two 2D ndarrays.")
    
    # Plot
    import matplotlib.pyplot as plt
    if not coords:
        im = smap.plot(annotate=True, cmap=cmap, **kwargs)
    else:
        X, Y = coords
        im = ax.pcolormesh(X, Y, smap.data.T, cmap=cmap, **kwargs)
        ax.set_title(smap.latex_name, y=1.02)
        ax.set_xlabel('X (arcsec)')
        ax.set_ylabel('Y (arcsec)')
    ax.set_aspect(1)
    ax.grid(grid)
    ax.set_autoscale_on(False)  # Disable autoscaling from other plots
    plt.subplots_adjust(left=0.1, bottom=0.11, right=0.9, top=0.88)  # xmin,ymin,xmax,ymax (default=[0.125, 0.11, 0.9, 0.88])
    cax = plt.axes([0.9-0.015, 0.11, 0.03*0.8, 0.77])  # [x=(xmax+pad),y=ymin,w=ratio*(xmax-xmin),h=(ymax-ymin)]
    plt.colorbar(im, cax=cax)
    return im

#======================================================================|
def plot_vmap(ax, mapu, mapv, mapc, coords=None,
              iskip=10, jskip=10, cmin=0., vmax=1000., cmap='binary',
              scale_units='xy', scale=1/0.05, minlength=0.05, width=0.003,
              headlength=6, headwidth=5, headaxislength=3,
              **kwargs):
    '''
    ----------------------------------------------------------------------------
    Vector plot.
    
    [Plot Function] `quiver` from matplotlib
    
    [Parameters]
    - ax: matplotlib axes object
    - mapu: a sunpy `GenericMap` of Vector_x
    - mapv: a sunpy `GenericMap` of Vector_y
    - mapc: a sunpy `GenericMap` to set color values
    - coords: two 2D numpy `ndarray`s: (X, Y)
    - iskip, jskip: number of skipped points in both dimensions
    - cmin: where mapc.data < cmin => set to zero
    - vmax: where norm(Vector) > vmax => set to vmax
    - cmap: name of a color map
    - scale_units, ..., kwargs: kwargs of `quiver`
    
    [Return] a matplotlib `artist`(image object)
    
    [See also]
    https://matplotlib.org/devdocs/api/_as_gen/matplotlib.axes.Axes.quiver.html#matplotlib-axes-axes-quiver
    
    ----------------------------------------------------------------------------
    '''
    if not any(isinstance(i, sunpy.map.mapbase.GenericMap) for i in (mapu, mapv, mapc)):
        raise TypeError("mapu, mapv, mapc should be sunpy GenericMaps.")
    if coords and (len(coords) != 2 or not any(isinstance(i, np.ndarray) for i in coords)):
        raise ValueError("coords should be a list of two 2D ndarrays.")
    
    iskip = int(iskip); jskip = int(jskip)
    # Resample
    if not coords:
        dimy, dimx = mapc.data.shape
        pixmin = (0., 0.)  # pixel
        pixmax = u.Quantity(mapc.dimensions).value
        xmin = (mapc.pixel_to_world(pixmin[0]*u.pix,pixmin[1]*u.pix).Tx).to(u.deg).value  # deg
        xmax = (mapc.pixel_to_world(pixmax[0]*u.pix,pixmax[1]*u.pix).Tx).to(u.deg).value
        ymin = (mapc.pixel_to_world(pixmin[0]*u.pix,pixmin[1]*u.pix).Ty).to(u.deg).value
        ymax = (mapc.pixel_to_world(pixmax[0]*u.pix,pixmax[1]*u.pix).Ty).to(u.deg).value
        rmapu = mapu.data.T[::iskip, ::jskip].copy()
        rmapv = mapv.data.T[::iskip, ::jskip].copy()
        rmapc = mapc.data.T[::iskip, ::jskip].copy()
        X, Y = np.mgrid[xmin:xmax:dimx*1j, ymin:ymax:dimy*1j][:, ::iskip, ::jskip]  # deg
    else:
        from scipy.interpolate import RectBivariateSpline
        X, Y = coords
        P, L0, B0, Bc, Lc, xmin, xmax, ymin, ymax, dimx, dimy, *_ = _get_image_params(mapc)
        ax1, ax2, ay1, ay2 = proj_matrix(P, L0, B0, Bc, Lc, 2)
        hx = np.linspace(X.min(), X.max(), dimx)
        hy = np.linspace(Y.min(), Y.max(), dimy)
        hx, hy = np.mgrid[X.min():X.max():dimx*1j, Y.min():Y.max():dimy*1j]
        ix = ax1 * hx + ay1 * hy
        iy = ax2 * hx + ay2 * hy
        # Interpolate to 2D spline functions (over a rectangular mesh)
        fu = RectBivariateSpline(np.linspace(xmin, xmax, dimx), np.linspace(ymin, ymax, dimy), mapu.data.T)
        fv = RectBivariateSpline(np.linspace(xmin, xmax, dimx), np.linspace(ymin, ymax, dimy), mapv.data.T)
        fc = RectBivariateSpline(np.linspace(xmin, xmax, dimx), np.linspace(ymin, ymax, dimy), mapc.data.T)
        # Resample
        rmapu = np.array(list(map(fu, ix.flatten(), iy.flatten()))).reshape((dimx, dimy))[::iskip, ::jskip].copy()
        rmapv = np.array(list(map(fv, ix.flatten(), iy.flatten()))).reshape((dimx, dimy))[::iskip, ::jskip].copy()
        rmapc = np.array(list(map(fc, ix.flatten(), iy.flatten()))).reshape((dimx, dimy))[::iskip, ::jskip].copy()
        X = hx[::iskip, ::jskip]; Y = hy[::iskip, ::jskip]
    # Clip
    mag = np.sqrt(rmapu**2 + rmapv**2)
    mask = np.where(abs(rmapc) < cmin)
    rmapu[mask] = 0.
    rmapv[mask] = 0.
    mask = np.where(mag > vmax)
    rmapu[mask] = rmapu[mask] * (vmax / mag[mask])
    rmapv[mask] = rmapv[mask] * (vmax / mag[mask])
    rmapc = np.sign(rmapc)
    # Plot
    if not coords: kwargs['transform'] = ax.get_transform('world')
    im = ax.quiver(X, Y, rmapu, rmapv, rmapc,
                   cmap=cmap,
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
    ----------------------------------------------------------------------------
    Transform maps from image-coordinate to helio-coordinate.
    - Helo-coordinate: Helioprojective(Cartesian) system
    - Matrix: A22 or A33 get from `usr_sunpy.proj_matrix`
    - Unit: arcsec

    [Parameters]
    - smap: list of sunpy `GenericMap`, one or three elements.
      - For scalar: `image_to_helio(smap)`
      - For vectors: `image_to_helio(smapx, smapy, smapz)`

    [Return]
    - For scalar: a list of numpy `ndarray`: x_h, y_h
    - For vectors: a list of sunpy `GenericMap`: smapx_h, smapy_h, smapz_h

    [See also]
    http://docs.sunpy.org/en/stable/code_ref/coordinates.html#sunpy-coordinates

    ----------------------------------------------------------------------------
    '''
    if not any(isinstance(i, sunpy.map.mapbase.GenericMap) for i in smap):
        raise TypeError("*smap should be 1 or 3 sunpy GenericMaps.")
    
    from copy import deepcopy
    
    P, L0, B0, Bc, Lc, xmin, _, ymin, _, dimx, dimy, dx, dy = _get_image_params(smap[-1])
    
    ix, iy = np.mgrid[0:dimx-1:dimx*1j, 0:dimy-1:dimy*1j]  # image-coordinate
    ix *= dx; ix += xmin
    iy *= dy; iy += ymin
    
    ax1, ax2, ax3, ay1, ay2, ay3, az1, az2, az3 = proj_matrix(P, L0, B0, Bc, Lc)
    
    if len(smap) == 1:  # Scalar
        [[cx1, cx2], [cy1, cy2]] = np.array(np.mat([[ax1, ay1], [ax2, ay2]]).I)
        hx = cx1 * ix + cx2 * iy
        hy = cy1 * ix + cy2 * iy
        return hx, hy
    elif len(smap) == 3:  # Vector
        hmapx = deepcopy(smap[0])
        hmapy = deepcopy(smap[1])
        hmapz = deepcopy(smap[2])
        hmapx.data[:] = ax1 * smap[0].data + ax2 * smap[1].data + ax3 * smap[2].data
        hmapy.data[:] = ay1 * smap[0].data + ay2 * smap[1].data + ay3 * smap[2].data
        hmapz.data[:] = az1 * smap[0].data + az2 * smap[1].data + az3 * smap[2].data
        return hmapx, hmapy, hmapz
    else:
        raise ValueError('The number of arguments must be 1 or 3.')

#======================================================================|
def proj_matrix(P, L0, B0, Bc, Lc, *dim):
    '''
    ----------------------------------------------------------------------------
    - For coords: (x, y)_helio = A22.T.I * (x, y)_image
    - For vectors: (U, V, W)_helio = A33 * (U, V, W)_image

     A33 = [[ax1, ay1, az1],
            [ax2, ay2, az2],
            [ax3, ay3, az3]]
    
    [Return] elements of a 2D array
            (use values instead of arrays just for easy reading & comparison)
    - default: ax1, ax2, ax3, ay1, ay2, ay3, az1, az2, az3
    -   dim=2: ax1, ax2, ay1, ay2
    
    [Parameters]
    -  P: the angle of the northern extremity, CCW from the north point of the disk.
    - L0: the longitude of the center of the disk.
    - B0: the latitude of the center of the disk.
    - Lc: the longitude of the the referenced point.
    - Bc: the latitude of the referenced point.
    
    [Reference]
    http://link.springer.com/10.1007/BF00158295

    ----------------------------------------------------------------------------
    '''
    ax1 =  -np.sin(B0) * np.sin(P) * np.sin(Lc - L0) + np.cos(P) * np.cos(Lc - L0)
    ax2 =   np.sin(B0) * np.cos(P) * np.sin(Lc - L0) + np.sin(P) * np.cos(Lc - L0)
    ay1 = (-np.sin(Bc) * (np.sin(B0) * np.sin(P) * np.cos(Lc - L0) + np.cos(P) * np.sin(Lc - L0))
          - np.cos(Bc) *  np.cos(B0) * np.sin(P))
    ay2 = ( np.sin(Bc) * (np.sin(B0) * np.cos(P) * np.cos(Lc - L0) - np.sin(P) * np.sin(Lc - L0))
          + np.cos(Bc) *  np.cos(B0) * np.cos(P))
    if dim and dim[0] == 2:
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

#======================================================================|
def _get_image_params(smap):
    '''
    [Return] P, L0, B0, Bc, Lc, xmin, xmax, ymin, ymax, dimy, dimx, dx, dy
    '''
    dimy, dimx = smap.data.shape  # dimy(vertical) goes first
    P = 0.  # The angle of the northern extremity, CCW from the north point of the disk.
    L0 = np.deg2rad(smap.heliographic_longitude.value)  # The longitude of the center of the disk.
    B0 = np.deg2rad(smap.heliographic_latitude.value)  # The latitude of the center of the disk.
    dx, dy = (smap.meta['cdelt1'], smap.meta['cdelt2'])  # arcsec/pix
    xmin = - dx * (dimx - 1) / 2.  # Set (0, 0) at the center of the image.
    ymin = - dy * (dimy - 1) / 2.
    xmax = dx * (dimx - 1) / 2.
    ymax = dy * (dimy - 1) / 2.
    _tmp = smap.pixel_to_world(*((dimx-1)/2.,(dimy-1)/2.)*u.pix).transform_to('heliographic_stonyhurst')
    Lc = np.deg2rad(_tmp.lon.value)  # The longitude of the center of the image.
    Bc = np.deg2rad(_tmp.lat.value)  # The latitude of the center of the image.
    return P, L0, B0, Bc, Lc, xmin, xmax, ymin, ymax, dimx, dimy, dx, dy
