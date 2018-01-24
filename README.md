## SunPy scripts & examples
> Tested: Python 2.7 & Python 3.4, Ubuntu 14.04 LTS

* *Update: 2018-01-24*

<br />

###### Download

``` sh
git clone https://git.coding.net/lydiazly/scripts-sunpy.git
```

---

#### Install [<u>SunPy</u>](http://sunpy.org)

    pip install sunpy[all]

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

``` python
>>> import usr_sunpy
>>> help(usr_sunpy)
```

---

#### Examples

Scripts:

    python plothmi.py

IPython notebooks:

    ipython notebook example_plothmi.ipynb

&nbsp;&nbsp;&nbsp;&nbsp;or

    jupyter notebook example_plothmi.ipynb
