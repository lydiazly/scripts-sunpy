<h1 id="usr_sunpy">usr_sunpy</h1>


User functions.
- SunPy Version: 0.9.0
- Reference http://docs.sunpy.org/en/stable/code_ref/map.html

<h2 id="usr_sunpy.read_sdo">read_sdo</h2>

```python
read_sdo(filename)
```

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

<h2 id="usr_sunpy.plot_map">plot_map</h2>

```python
plot_map(ax, smap, coords=None, grid=False, cmap='gray', **kwargs)
```

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

<h2 id="usr_sunpy.plot_vmap">plot_vmap</h2>

```python
plot_vmap(ax, mapu, mapv, mapc, coords=None, iskip=10, jskip=10, cmin=0.0, vmax=1000.0, cmap='binary', scale_units='xy', scale=20.0, minlength=0.05, width=0.003, headlength=6, headwidth=5, headaxislength=3, **kwargs)
```

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

<h2 id="usr_sunpy.image_to_helio">image_to_helio</h2>

```python
image_to_helio(*smap)
```

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

<h2 id="usr_sunpy.proj_matrix">proj_matrix</h2>

```python
proj_matrix(P, L0, B0, Bc, Lc, *dim)
```

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

