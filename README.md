## SunPy scripts & examples
> Tested:<br />
> Ubuntu 14.04 LTS<br />
> Python 2.7 & Python 3.4 - SunPy(0.8.5), Astropy(2.0.4)<br />
> Python 3.5 - SunPy(0.8.5), Astropy(3.0.1)

* *Update: 2018-03-14*

<br />

###### Download

``` sh
git clone https://git.coding.net/lydiazly/scripts-sunpy.git
```

---

#### Install [<u>SunPy</u>](http://sunpy.org)

    pip install sunpy[all]

(See [<u>sunpy_install_troubleshooting</u>](https://coding.net/u/lydiazly/p/scripts-sunpy/git/blob/master/sunpy_install_troubleshooting.txt))

#### Import user modules

Append&nbsp;&nbsp;<span style="color:#000000">*scripts-sunpy/modules*</span>&nbsp;&nbsp;
into&nbsp;&nbsp;<span style="color:#445eac">*$PYTHONPATH*</span>.

&nbsp;&nbsp;*e.g.* For Bash users, in ~/.bashrc:

``` sh
export PYTHONPATH=$PYTHONPATH:<YOUR_PATH>/scripts-sunpy/modules
```

Import:

``` python
>>> from usr_sunpy import *
```

To get help:

Python:

``` python
>>> import usr_sunpy
>>> help(usr_sunpy)
```

Shell:

``` sh
pydoc usr_sunpy
```

---

#### Examples

Scripts:

    python plothmi.py

IPython notebooks:

    ipython notebook example_plothmi.ipynb

&nbsp;&nbsp;&nbsp;&nbsp;or

    jupyter notebook example_plothmi.ipynb
