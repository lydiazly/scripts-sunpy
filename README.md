# SunPy scripts & examples

> Tested ( *Update: 2018-06-01* )<br>
> Ubuntu 16.04 LTS<br>
> Python [2.7.12] Sunpy [0.9.0] Astropy [2.0.6]<br>
> Python [3.6.5] Sunpy [0.9.0] Astropy [3.0.2]<br>

---

## Download

    git clone https://git.coding.net/lydiazly/scripts-sunpy.git

## Install [SunPy](http://sunpy.org)

( See also: [sunpy_troubleshooting.txt](https://coding.net/u/lydiazly/p/scripts-sunpy/git/blob/master/sunpy_troubleshooting.txt) )

## Examples

Scripts:

    python plothmi.py

IPython notebooks:

    jupyter notebook example_plothmi.ipynb

## Get help

Python:

``` python
>>> import usr_sunpy
>>> help(usr_sunpy)
```

Shell:

``` sh
pydoc usr_sunpy
```

## Use user modules

e.g. *modules/usr_sunpy*

[Funcions] *read_sdo, plot_map, plot_vmap, image_to_helio, ...*

Append *./modules* to *$PYTHONPATH* and import:

``` python
>>> from usr_sunpy import *
```

Or

``` python
>>> import sys
>>> sys.path.append('../modules')  # If current directory is plothmi/
>>> from usr_sunpy import *
```
