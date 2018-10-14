'''
- SunPy Version: 0.9.3
- Reference http://docs.sunpy.org/en/stable/code_ref/map.html
- Change log:
  - 2018-10-11
    - Add `aiaprep_usr()`.

------------------------------------------------------------------------
'''
# 2017-12-11 written by Lydia
# 2018-10-12 modified by Lydia
from __future__ import absolute_import, division, print_function

__all__ = ['read_sdo', 'tai', 'aiaprep_usr']

#======================================================================|
def read_sdo(filename):
    '''
    Example of reading from a FITS file. Print the filename and dimensions.
    
    [Returns] a sunpy `GenericMap` or subclass(e.g. `HMIMap`) object
    - data - a 2D numpy `ndarray`
      data[i, j]: i from the bottom(y), j from te left(x).
    - meta - a `dict` of the original image headr tags.
    
    [See also]
    - help(sunpy.map.GenericMap)
    - http://docs.sunpy.org/en/stable/code_ref/map.html#using-map-objects

    --------------------------------------------------------------------
    '''
    import astropy.units as u
    import sunpy.map
    import os
    smap = sunpy.map.Map(filename)  # Read & return a `GenericMap`
    print('%s\t%s' % (os.path.basename(filename),
          list(map(int, u.Quantity(smap.dimensions).value))))
    return smap

#======================================================================|
def tai(*timestr):
    '''
    Warp time strings as TAI time.
    
    [Parameters]
    - timestr: time strings
      e.g. '2010-01-01T00:00:00', '2010.01.01_00:00:00_TAI'

    [Returns]
    - len(timestr) == 1: a TAI `Time` object
    - len(timestr) > 1: a list of TAI `Time` objects

    [See also]
    - http://docs.astropy.org/en/stable/time/
    - http://docs.sunpy.org/en/stable/guide/time.html

    --------------------------------------------------------------------
    '''
    from sunpy.time import parse_time
    try:
        tmp = [parse_time(i) for i in timestr]
    except ValueError as e:
        print('ValueError:', e)
    import astropy.time
    taitime = [astropy.time.Time(i, scale='tai') for i in tmp]
    return taitime if len(taitime) > 1 else taitime[0]
#======================================================================|
def aiaprep_usr(aiamap, order=3):
    """
    **Modified [`sunpy.instr.aia.aiaprep()`](https://docs.sunpy.org/en/latest/api/sunpy.instr.aia.aiaprep.html)**
    
    Processes a level 1 `~sunpy.map.sources.sdo.AIAMap` into a level 1.5
    `~sunpy.map.sources.sdo.AIAMap`. Rotates, scales and
    translates the image so that solar North is aligned with the y axis, each
    pixel is 0.6 arcsec across, and the center of the sun is at the center of
    the image. The actual transformation is done by Map's
    :meth:`~sunpy.map.mapbase.GenericMap.rotate` method.

    This function is similar in functionality to aia_prep() in SSWIDL, but
    it does not use the same transformation to rotate the image and it handles
    the meta data differently. It should therefore not be expected to produce
    the same results.

    Parameters
    ----------
    aiamap : `~sunpy.map.sources.sdo.AIAMap` instance
        A `sunpy.map.Map` from AIA

    Returns
    -------
    newmap : A level 1.5 copy of `~sunpy.map.sources.sdo.AIAMap`

    Notes
    -----
    This routine makes use of Map's :meth:`~sunpy.map.mapbase.GenericMap.rotate`
    method, which modifies the header information to the standard PCi_j WCS
    formalism.
    The FITS header resulting in saving a file after this procedure will
    therefore differ from the original file.

    --------------------------------------------------------------------
    """
    import numpy as np
    import astropy.units as u
    from sunpy.map.sources.sdo import AIAMap, HMIMap
    if not isinstance(aiamap, (AIAMap, HMIMap)):
        raise ValueError("Input must be an AIAMap")

    # Target scale is 0.6 arcsec/pixel, but this needs to be adjusted if the map
    # has already been rescaled.
    if (aiamap.scale[0] / 0.6).round() != 1.0 * u.arcsec and aiamap.data.shape != (4096, 4096):
        scale = (aiamap.scale[0] / 0.6).round() * 0.6 * u.arcsec
    else:
        scale = 0.6 * u.arcsec  # pragma: no cover # can't test this because it needs a full res image
    scale_factor = aiamap.scale[0] / scale

    # tempmap = aiamap.rotate(recenter=True, scale=scale_factor.value, missing=aiamap.min())
    tempmap = aiamap.rotate(order=order, recenter=True, scale=scale_factor.value, missing=np.nan)

    # extract center from padded aiamap.rotate output
    # crpix1 and crpix2 will be equal (recenter=True), as aiaprep does not work with submaps
    center = np.floor(tempmap.meta['crpix1'])
    range_side = (center + np.array([-1, 1]) * aiamap.data.shape[0] / 2) * u.pix
    newmap = tempmap.submap(u.Quantity([range_side[0], range_side[0]]),
                            u.Quantity([range_side[1], range_side[1]]))

    newmap.meta['r_sun'] = newmap.meta['rsun_obs'] / newmap.meta['cdelt1']
    newmap.meta['lvl_num'] = 1.5

    return newmap
