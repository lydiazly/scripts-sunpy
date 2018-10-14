'''
- SunPy Version: 0.9.3
- Reference http://docs.sunpy.org/en/stable/code_ref/map.html
- Change log:
  - 2018-10-11
    - Change parameter list of `plot_map()`, `plot_vmap()`, `proj_matrix()`.
    - Change return value of `_get_image_params()` from a tuple to a dict.
    - Adjust positions of colorbar & title.
    - Fix docs.

------------------------------------------------------------------------
'''
# 2017-12-11 written by Lydia
# 2018-10-13 modified by Lydia
from __future__ import absolute_import, division, print_function
import astropy.units as u
import numpy as np
import sunpy.map
from sunpy.visualization import wcsaxes_compat

__all__ = ['plot_map', 'plot_vmap', 'image_to_helio', 'proj_matrix']

#======================================================================|
def plot_map(smap, ax=None, coords=None, annotate=True, title=True, colorbar=True,
             grid=True, grid_color='yellow', grid_ls=':', grid_lw=0.8, grid_alpha=0.5, **kwargs):
    '''
    Plot image.
    
    [Plot Function]
    - plot_map(smap, **kwargs) -> use smap.plot(), `imshow` from matplotlib
    - plot_map(smap, coords=(X, Y), **kwargs) -> `pcolormesh` from matplotlib

    [Parameters]
    - smap: a sunpy `GenericMap`
    - ax: a matplotlib `Axes` object
    - coords: tuple or list of `ndarray`s: (X, Y)
    - grid: bool or `~astropy.units.Quantity`(spacing), zorder = 90
    - grid_color: grid color
    - grid_ls: grid line stlye
    - grid_lw: grid line width
    - grid_alpha: grid alpha
    - cmap: name of color map
    - kwargs: dict, matplotlib kwargs
      If coords is None(default), use kwargs of `imshow`,
      else of `pcolormesh`.

    [Returns] a matplotlib image object

    [See also]
    http://docs.sunpy.org/en/stable/code_ref/map.html#sunpy.map.mapbase.GenericMap.plot

    --------------------------------------------------------------------
    '''
    # Compatible with old syntax `plot_map(ax, smap, ...)`
    if ax and isinstance(ax, sunpy.map.mapbase.GenericMap):
        smap, ax = ax, smap
    
    if not isinstance(smap, sunpy.map.mapbase.GenericMap):
        raise TypeError("smap should be a sunpy GenericMap.")
    if coords is not None and (len(coords) != 2 or not any(isinstance(i, np.ndarray) for i in coords)):
        raise ValueError("coords should be a list of two 2D ndarrays.")
    if ax is None:
        ax = wcsaxes_compat.gca_wcs(smap.wcs)
    
    if coords is not None and isinstance(grid, bool):
        grid = False
    if (isinstance(grid, bool) and not grid) or all(np.isnan(smap.data[:, -1])):
        # No right yaxis label
        _pad = 0.045  # colorbar
        _right = 0.97  # subplots_adjust
    else:
        _pad = 0.12  # colorbar
        _right = 0.99  # subplots_adjust
    if (isinstance(grid, bool) and not grid) or all(np.isnan(smap.data[-1, :])):
        # No top xaxis label
        _y = 1.02  # title
        _top = 0.9  # subplots_adjust
    else:
        _y = 1.12  # title
        _top = 0.84  # subplots_adjust
    # Plot
    import matplotlib.pyplot as plt
    if coords is None:
        im = smap.plot(axes=ax, annotate=annotate, title=title, **kwargs)
    else:
        X, Y = coords
        im = ax.pcolormesh(X, Y, smap.data.T, **kwargs)
        if annotate:
            if title is True:
                title = smap.latex_name
            if title:
                ax.set_title(title)
            ax.set_xlabel('X (arcsec)')
            ax.set_ylabel('Y (arcsec)')
    ax.set_aspect(1)
    ax.grid(False)  # rectangular grids
    grid_kw = dict(color=grid_color, ls=grid_ls, lw=grid_lw, alpha=grid_alpha)
    if isinstance(grid, bool):
        if grid:
            smap.draw_grid(axes=ax, zorder=90, **grid_kw)
    elif isinstance(grid, u.Quantity):
        smap.draw_grid(axes=ax, grid_spacing=grid, **grid_kw)
    else:
        raise TypeError("grid should be a bool or an astropy Quantity.")
    if annotate and title:
        ax.set_title(ax.get_title(), y=_y)
    ax.set_autoscale_on(False)  # Disable autoscaling from other plots
    plt.subplots_adjust(left=0.11, bottom=0.1, right=_right, top=_top)  # xmin,ymin,xmax,ymax (default=[0.125, 0.1, 0.9, 0.9])
    # cax = plt.axes([0.9+0.1, 0.1, 0.03*0.8, 0.78])  # [x=(xmax+pad),y=ymin,w=ratio*(xmax-xmin),h=(ymax-ymin)]
    # plt.colorbar(im, ax=ax, cax=cax)
    if colorbar:
        plt.colorbar(im, ax=ax, pad=_pad)
    # pip install mpldatacursor
    # https://github.com/joferkington/mpldatacursor
    # try:
    #     from mpldatacursor import datacursor
    #     datacursor(im)
    # except:
    #     pass
    return im

#======================================================================|
def plot_vmap(mapu, mapv, mapc, ax=None, coords=None, cmap='binary',
              iskip='auto', jskip='auto',
              cmin=0., cmax=None, vmin=None, vmax=1000.,
              scale_units='xy', scale=1/0.05, minlength=0.05, width=0.003,
              headlength=6.5, headwidth=5, headaxislength=3.5,
              **kwargs):
    '''
    Quiver plot of (U, V).
    
    [Plot Function] `quiver` from matplotlib, zorder(default) = 100
    
    [Parameters]
    - mapu: a sunpy `GenericMap` of U
    - mapv: a sunpy `GenericMap` of V
    - mapc: a sunpy `GenericMap` to set color values
    - ax: matplotlib axes object
    - coords: tuple or list of 2D `ndarray`s: (X, Y)
    - iskip, jskip: 'auto' or number of skipped points
    - cmin: where mapc.data < cmin => set U, V to zero
    - cmax: where mapc.data > cmax => set U, V to zero
    - vmin: clip norm(U, V) to vmin
    - vmax: clip norm(U, V) to vmax
    - cmap: name of a color map
    - scale_units, ..., **kwargs: kwargs of `quiver`
    
    [Returns] a matplotlib `artist`(image object)
    
    [See also]
    https://matplotlib.org/devdocs/api/_as_gen/matplotlib.axes.Axes.quiver.html#matplotlib-axes-axes-quiver

    --------------------------------------------------------------------
    '''
    # Compatible with old syntax `plot_vmap(ax, mapu, mapv, mapc, ...)`
    if ax and isinstance(ax, sunpy.map.mapbase.GenericMap):
        _ax, _mapu, _mapv, _mapc = mapu, mapv, mapc, ax
        mapu, mapv, mapc, ax = _mapu, _mapv, _mapc, _ax
        del _mapu, _mapv, _mapc, _ax
    
    if not any(isinstance(i, sunpy.map.mapbase.GenericMap) for i in (mapu, mapv, mapc)):
        raise TypeError("mapu, mapv, mapc should be sunpy GenericMaps.")
    if coords is not None and (len(coords) != 2 or not any(isinstance(i, np.ndarray) for i in coords)):
        raise ValueError("coords should be a list of two 2D ndarrays.")
    if ax is None:
        ax = wcsaxes_compat.gca_wcs(smap.wcs)
    
    dimy, dimx = mapc.data.shape
    _ngx = 52  # default points number of xaxis
    iskip = dimx // _ngx if iskip == 'auto' else int(iskip)
    jskip = dimx // _ngx if jskip == 'auto' else int(jskip)
    # Resample
    if coords is None:
        pixmin = (0., 0.)  # pix
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
        _par = _get_image_params(mapc)
        a11, a12, a21, a22 = proj_matrix(*[_par[_] for _ in ('P', 'L0', 'B0', 'Lc', 'Bc')], 2)
        hx, hy = np.mgrid[X.min():X.max():dimx*1j, Y.min():Y.max():dimy*1j]
        # (x, y)_image = A22.T * (x, y)_helio
        ix = a11 * hx + a21 * hy
        iy = a12 * hx + a22 * hy
        # Interpolate to 2D spline functions (over a rectangular mesh)
        fu = RectBivariateSpline(np.linspace(_par['xmin'], _par['xmax'], dimx), np.linspace(_par['ymin'], _par['ymax'], dimy), mapu.data.T)
        fv = RectBivariateSpline(np.linspace(_par['xmin'], _par['xmax'], dimx), np.linspace(_par['ymin'], _par['ymax'], dimy), mapv.data.T)
        fc = RectBivariateSpline(np.linspace(_par['xmin'], _par['xmax'], dimx), np.linspace(_par['ymin'], _par['ymax'], dimy), mapc.data.T)
        # Resample
        rmapu = np.array(list(map(fu, ix.flatten(), iy.flatten()))).reshape((dimx, dimy))[::iskip, ::jskip].copy()
        rmapv = np.array(list(map(fv, ix.flatten(), iy.flatten()))).reshape((dimx, dimy))[::iskip, ::jskip].copy()
        rmapc = np.array(list(map(fc, ix.flatten(), iy.flatten()))).reshape((dimx, dimy))[::iskip, ::jskip].copy()
        X = hx[::iskip, ::jskip]; Y = hy[::iskip, ::jskip]
    # Clip
    if cmin is not None or cmax is not None:
        mask = True
        if cmin is not None:
            mask = mask & (abs(rmapc) < cmin)
        if cmax is not None:
            mask = mask & (abs(rmapc) > cmax)
        rmapu[mask] = rmapv[mask] = 0.
    if vmin is not None or vmax is not None:
        mag = np.sqrt(rmapu**2 + rmapv**2)
        if vmin is not None:
            mask = mag < vmin
            rmapu[mask] = rmapu[mask] * (vmin / mag[mask])
            rmapv[mask] = rmapv[mask] * (vmin / mag[mask])
        if vmax is not None:
            mask = mag > vmax
            rmapu[mask] = rmapu[mask] * (vmax / mag[mask])
            rmapv[mask] = rmapv[mask] * (vmax / mag[mask])
    # Plot
    rmapc = np.sign(rmapc)
    if not coords: kwargs['transform'] = ax.get_transform('world')
    if 'zorder' in kwargs:
        _zorder = kwargs.pop('zorder')
    else:
        _zorder = 100
    im = ax.quiver(X, Y, rmapu, rmapv, rmapc,
                   cmap=cmap,
                   clim=(-1, 1),
                   angles='uv', pivot='tail',
                   scale_units=scale_units, scale=scale,
                   minlength=minlength, width=width,
                   headlength=headlength, headwidth=headwidth, headaxislength=headaxislength,
                   zorder=_zorder,
                   **kwargs)
    return im

#======================================================================|
def image_to_helio(*smap):
    '''
    Transform maps from image-coordinate to helio-coordinate.
    - Helo-coordinate: Helioprojective(Cartesian) system
    - Matrix: A22 or A33 get from `usr_sunpy.proj_matrix`
    - Unit: arcsec

    [Parameters]
    - smap: `GenericMap`, one or three elements.
      - For scalar: `image_to_helio(smap)`
      - For vectors: `image_to_helio(smapx, smapy, smapz)`

    [Returns]
    - For scalar: x_h, y_h (type: `ndarray`)
    - For vectors: smapx_h, smapy_h, smapz_h (type: `GenericMap`)

    [See also]
    http://docs.sunpy.org/en/stable/code_ref/coordinates.html#sunpy-coordinates

    --------------------------------------------------------------------
    '''
    if not any(isinstance(i, sunpy.map.mapbase.GenericMap) for i in smap):
        raise TypeError("*smap should be 1 or 3 sunpy GenericMaps.")
    from copy import deepcopy
    _par = _get_image_params(smap[-1])
    ix, iy = np.mgrid[0:_par['dimx']-1:_par['dimx']*1j, 0:_par['dimy']-1:_par['dimy']*1j]  # image-coordinate
    ix *= _par['dx']; ix += _par['xmin']
    iy *= _par['dy']; iy += _par['ymin']
    a11, a12, a13, a21, a22, a23, a31, a32, a33 = proj_matrix(*[_par[_] for _ in ('P', 'L0', 'B0', 'Lc', 'Bc')])
    if len(smap) == 1:  # Scalar
        [[cx1, cx2], [cy1, cy2]] = np.array(np.mat([[a11, a21], [a12, a22]]).I)
        hx = cx1 * ix + cx2 * iy
        hy = cy1 * ix + cy2 * iy
        return hx, hy
    elif len(smap) == 3:  # Vector
        hmapx = deepcopy(smap[0])
        hmapy = deepcopy(smap[1])
        hmapz = deepcopy(smap[2])
        hmapx.data[:] = a11 * smap[0].data + a12 * smap[1].data + a13 * smap[2].data
        hmapy.data[:] = a21 * smap[0].data + a22 * smap[1].data + a23 * smap[2].data
        hmapz.data[:] = a31 * smap[0].data + a32 * smap[1].data + a33 * smap[2].data
        return hmapx, hmapy, hmapz
    else:
        raise ValueError('The number of arguments must be 1 or 3.')

#======================================================================|
def proj_matrix(P, L0, B0, Lc, Bc, *dim):
    '''
    - For coords: (x, y)_helio = A22.T.I * (x, y)_image
    - For vectors: (U, V, W)_helio = A33 * (U, V, W)_image

     A33 = [[a11, a12, a13],
            [a21, a22, a23],
            [a31, a32, a33]]
    
    [Parameters]
    -  P: the angle of the northern extremity, CCW from the north point of the disk.
    - L0: the longitude of the center of the disk.
    - B0: the latitude of the center of the disk.
    - Lc: the longitude of the the referenced point.
    - Bc: the latitude of the referenced point.
    
    [Returns] elements of a 2D array
            (use values instead of arrays just for easy reading & comparison)
    - default: a11, a12, a13, a21, a22, a23, a31, a32, a33
    -   dim=2: a11, a12, a21, a22

    [Reference]
    http://link.springer.com/10.1007/BF00158295

    --------------------------------------------------------------------
    '''
    a11 =  -np.sin(B0) * np.sin(P) * np.sin(Lc - L0) + np.cos(P) * np.cos(Lc - L0)
    a12 =   np.sin(B0) * np.cos(P) * np.sin(Lc - L0) + np.sin(P) * np.cos(Lc - L0)
    a21 = (-np.sin(Bc) * (np.sin(B0) * np.sin(P) * np.cos(Lc - L0) + np.cos(P) * np.sin(Lc - L0))
          - np.cos(Bc) *  np.cos(B0) * np.sin(P))
    a22 = ( np.sin(Bc) * (np.sin(B0) * np.cos(P) * np.cos(Lc - L0) - np.sin(P) * np.sin(Lc - L0))
          + np.cos(Bc) *  np.cos(B0) * np.cos(P))
    if dim and dim[0] == 2:
        return a11, a12, a21, a22
    else:
        a13 =  -np.cos(B0) * np.sin(Lc - L0)
        a23 =  -np.cos(B0) *  np.sin(Bc) * np.cos(Lc - L0) + np.sin(B0) * np.cos(Bc)
        a31 = ( np.cos(Bc) * (np.sin(B0) * np.sin(P) * np.cos(Lc - L0) + np.cos(P) * np.sin(Lc - L0))
              - np.sin(Bc) *  np.cos(B0) * np.sin(P))
        a32 = (-np.cos(Bc) * (np.sin(B0) * np.cos(P) * np.cos(Lc - L0) - np.sin(P) * np.sin(Lc - L0))
              + np.sin(Bc) *  np.cos(B0) * np.cos(P))
        a33 =   np.cos(Bc) *  np.cos(B0) * np.cos(Lc - L0) + np.sin(Bc) * np.sin(B0)
        return a11, a12, a13, a21, a22, a23, a31, a32, a33

#======================================================================|
def _get_image_params(smap):
    '''
    [Returns] a dict with keys:
    'P', 'L0', 'B0', 'Bc', 'Lc', 'xmin', 'xmax', 'ymin', 'ymax', 'dimy', 'dimx', 'dx', 'dy'

    --------------------------------------------------------------------
    '''
    _par = {}
    _par['dimy'], _par['dimx'] = smap.data.shape  # dimy(vertical) goes first
    _par['P'] = 0.  # The angle of the northern extremity, CCW from the north point of the disk.
    _par['L0'] = np.deg2rad(smap.heliographic_longitude.value)  # The longitude of the center of the disk.
    _par['B0'] = np.deg2rad(smap.heliographic_latitude.value)  # The latitude of the center of the disk.
    _par['dx'], _par['dy'] = (smap.meta['cdelt1'], smap.meta['cdelt2'])  # arcsec/pix, smap.meta['cdelt1'] == smap.scale.axis1.value
    _par['xmin'] = - _par['dx'] * (_par['dimx'] - 1) / 2.  # Set (0, 0) at the center of the image.
    _par['ymin'] = - _par['dy'] * (_par['dimy'] - 1) / 2.
    _par['xmax'] = _par['dx'] * (_par['dimx'] - 1) / 2.
    _par['ymax'] = _par['dy'] * (_par['dimy'] - 1) / 2.
    _tmp = smap.pixel_to_world(*((_par['dimx']-1)/2.,(_par['dimy']-1)/2.)*u.pix).transform_to('heliographic_stonyhurst')
    _par['Lc'] = np.deg2rad(_tmp.lon.value)  # The longitude of the center of the image.
    _par['Bc'] = np.deg2rad(_tmp.lat.value)  # The latitude of the center of the image.
    return _par
