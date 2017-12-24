#SunPy scripts & examples
[ Update: 2017-12-19 ]


``` sh
git clone https://git.coding.net/lydiazly/scripts-sunpy.git
```

###Python version

Tested with python 2.7.

###Install [**SunPy**](http://sunpy.org/)

``` sh
pip install sunpy[all] pytest
```

###User modules

Append **modules/** to **PYTHONPATH**.

e.g. Bash and sh:

Add

``` sh
export PYTHONPATH=$PYTHONPATH:<SOME_PATH>/scripts-sunpy/modules
```
 
to ~/.bashrc
<br /><br />
######Usage

``` python
>>> from usr_sunpy import *
```

######Get help

``` python
>>> import usr_sunpy
>>> help(usr_sunpy)
```