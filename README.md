## SunPy scripts & examples
> Tested:<br />
> Ubuntu 16.04 LTS<br />
> Python 2.7.12: SunPy 0.9.0, Astropy 2.0.6<br />
> Python 3.5.2: SunPy 0.9.0, Astropy 3.0.2<br />
> Python 3.6.5: SunPy 0.9.0, Astropy 3.0.2

* *Update: 2018-04-01*

<br />

###### Download

``` sh
git clone https://git.coding.net/lydiazly/scripts-sunpy.git
```

---

#### Install [<u>SunPy</u>](http://sunpy.org)

(See also: [sunpy_troubleshooting.txt](https://coding.net/u/lydiazly/p/scripts-sunpy/git/blob/master/sunpy_troubleshooting.txt))

#### Use user modules

e.g. *modules/usr_sunpy*

Funcions: *read_sdo, plot_map, plot_vmap, image_to_helio, ...*

Append&nbsp;&nbsp;<span style="color:#000000">*./modules*</span>&nbsp;&nbsp;
to&nbsp;&nbsp;<span style="color:#445eac">$PYTHONPATH</span>&nbsp;&nbsp;and import:

``` python
>>> from usr_sunpy import *
```

Or:

``` python
>>> import sys
>>> sys.path.append('../modules')  # If current directory is plothmi/
>>> from usr_sunpy import *
```

Get help:

Python:

``` python
>>> import usr_sunpy
>>> help(usr_sunpy)
```

Shell:

``` sh
pydoc usr_sunpy
```

(See example:
[plothmi.py](https://coding.net/u/lydiazly/p/scripts-sunpy/git/blob/master/plothmi/plothmi.py))

---

#### Examples

Scripts:

    python plothmi.py

IPython notebooks:

    jupyter notebook example_plothmi.ipynb
