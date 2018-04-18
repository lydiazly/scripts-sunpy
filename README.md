## SunPy scripts & examples
> Tested:<br />
> Ubuntu 14.04 LTS<br />
> Python 2.7 - SunPy(0.8.5), Astropy(2.0.4)<br />
> Python 3.4 - SunPy(0.8.5), Astropy(2.0.4)<br />
> Python 3.5 - SunPy(0.8.5), Astropy(3.0.1)<br />
> Python 3.6 - SunPy(0.8.5), Astropy(3.0.1)


* *Update: 2018-04-01*

<br />

###### Download

``` sh
git clone https://git.coding.net/lydiazly/scripts-sunpy.git
```

---

#### Install [<u>SunPy</u>](http://sunpy.org)

(See also: [sunpy_troubleshooting.txt](https://coding.net/u/lydiazly/p/scripts-sunpy/git/blob/master/sunpy_troubleshooting.txt))

#### Use modules

[modules/usr_sunpy]

Funcions: *read_sdo, plot_map, plot_vmap, image_to_helio, ...*

Append&nbsp;&nbsp;<span style="color:#000000">*./modules*</span>&nbsp;&nbsp;
to&nbsp;&nbsp;<span style="color:#445eac">$PYTHONPATH</span>&nbsp;&nbsp;and import:

``` python
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
