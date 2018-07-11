'''
- SunPy Version: 0.9.0
- Reference http://docs.sunpy.org/en/stable/code_ref/map.html
'''
# 2017-12-11 written by Lydia
# 2018-07-11 modified by Lydia
from __future__ import absolute_import, division, print_function

__all__ = ['read_sdo', 'tai']

#======================================================================|
def read_sdo(filename):
    '''
    ----------------------------------------------------------------------------
    Example of reading from a FITS file. Print the filename and dimensions.
    
    [Returns] a sunpy `GenericMap` or subclass(e.g. `HMIMap`) object
    - data - a 2D numpy `ndarray`
      data[i, j]: i from the bottom(y), j from te left(x).
    - meta - a `dict` of the original image headr tags.
    
    [See also]
    - help(sunpy.map.GenericMap)
    - http://docs.sunpy.org/en/stable/code_ref/map.html#using-map-objects

    ----------------------------------------------------------------------------
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
    '''
    from sunpy.time import parse_time
    try:
        tmp = [parse_time(i) for i in timestr]
    except ValueError as e:
        print('ValueError:', e)
    import astropy.time
    taitime = [astropy.time.Time(i, scale='tai') for i in tmp]
    return taitime if len(taitime) > 1 else taitime[0]