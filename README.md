# SunPy scripts & examples

> Tested ( *Update: 2018-06-01* )<br>
> Ubuntu 16.04 LTS<br>
> Python [2.7.12] Sunpy [0.9.0] Astropy [2.0.6]<br>
> Python [3.6.5] Sunpy [0.9.0] Astropy [3.0.3]<br>

---

## Download

    git clone https://git.coding.net/lydiazly/scripts-sunpy.git

## Install [SunPy](http://sunpy.org)

( See also: [sunpy_troubleshooting.txt](https://coding.net/u/lydiazly/p/scripts-sunpy/git/blob/master/sunpy_troubleshooting.txt) )

## Examples

Run scripts:

    python plothmi.py

Open IPython notebooks:

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

## Import user modules

e.g. Import functions in *modules/usr_sunpy*

First append *./modules* to PYTHONPATH in your ~/.bashrc:

``` sh
export PYTHONPATH=$PYTHONPATH:<your_path>/scripts-sunpy/modules
```

Then import in python:

``` python
from usr_sunpy import *
```

Or add the path temporarily in python:

``` python
import sys
sys.path.append('../modules')  # If current location is plothmi/
from usr_sunpy import *
```
