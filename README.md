# SunPy scripts & examples

[ Tested ]

Ubuntu 16.04 LTS
Python|Sunpy|Astropy|Date
:---:|:---:|:---:|:---:
2.7.12| 0.9.0|2.0.6|2018-06-28
3.6.5| 0.9.0|3.0.3|2018-06-28

[ Update ]
> *2018-07-01*&emsp;上传了 pdf 版本
<br>

---

## Download

`$ git clone https://git.coding.net/lydiazly/scripts-sunpy.git`

---

## Install [SunPy](http://sunpy.org)

( See also: [sunpy_troubleshooting.txt](https://coding.net/u/lydiazly/p/scripts-sunpy/git/blob/master/sunpy_troubleshooting.txt) )

---

## Examples

* plothmi/example_plothmi.ipynb<br>
Download (
[Jupyter Notebook](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/plothmi/example_plothmi.ipynb)
|
[HTML](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/plothmi/example_plothmi.html)
|
[PDF](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/plothmi/example_plothmi.pdf)
)&ensp;
<a href="http://htmlpreview.github.io/?https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/plothmi/example_plothmi.html" target="_blank">
Preview
</a>

* plothmi/example_projection.ipynb<br>
Download (
[Jupyter Notebook](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/plothmi/example_projection.ipynb)
|
[HTML](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/plothmi/example_projection.html)
|
[PDF](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/plothmi/example_projection.pdf)
)&ensp;
<a href="http://htmlpreview.github.io/?https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/plothmi/example_projection.html" target="_blank">
Preview
</a>

Run scripts:

    $ python plothmi.py

Open Jupyter notebooks:

    $ jupyter notebook example_plothmi.ipynb

---

## User modules

### Get help

Python:

``` python
>>> import usr_sunpy
>>> help(usr_sunpy)
```

Shell:

    $ pydoc usr_sunpy

See this 
<a href="http://htmlpreview.github.io/?https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/usr_sunpy.html" target="_blank">
documentation preview
</a>
.

### Import

e.g. Import functions in [modules/usr_sunpy.py](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/modules/usr_sunpy.py)

First append&ensp;[modules/](https://coding.net/u/lydiazly/p/scripts-sunpy/git/tree/master/modules)&ensp;to&ensp;**PYTHONPATH**&ensp;in your&ensp;*~/.bashrc* :

``` sh
export PYTHONPATH=$PYTHONPATH:<your_path>/scripts-sunpy/modules
```

Then import in python:

``` python
>>> import usr_sunpy
```

Or add the path temporarily in python:

``` python
>>> import sys
>>> sys.path.append('../modules')  # If current location is plothmi/
>>> import usr_sunpy
```

See examples in [plothmi/](https://coding.net/u/lydiazly/p/scripts-sunpy/git/tree/master/plothmi) .
