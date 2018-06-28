# SunPy scripts & examples

> Tested ( *Update: 2018-06-28* )<br>
> Ubuntu 16.04 LTS<br>
> Python 2.7.12 | Sunpy 0.9.0 | Astropy 2.0.6<br>
> Python 3.6.5 | Sunpy 0.9.0 | Astropy 3.0.3<br>

---

## Download

    $ git clone https://git.coding.net/lydiazly/scripts-sunpy.git

## Install [SunPy](http://sunpy.org)

( See also: [sunpy_troubleshooting.txt](https://coding.net/u/lydiazly/p/scripts-sunpy/git/blob/master/sunpy_troubleshooting.txt) )

## Examples

* plothmi/example_plothmi.ipynb<br>
Download (
[HTML](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/plothmi/example_plothmi.html)
|
[Jupyter Notebook](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/plothmi/example_plothmi.ipynb)
)&ensp;
<a href="http://htmlpreview.github.io/?https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/plothmi/example_plothmi.html" target="_blank">
Preview
</a>

* plothmi/example_projection.ipynb<br>
Download (
[HTML](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/plothmi/example_projection.html)
|
[Jupyter Notebook](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/plothmi/example_projection.ipynb)
)&ensp;
<a href="http://htmlpreview.github.io/?https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/plothmi/example_projection.html" target="_blank">
Preview
</a>

Run scripts:

    $ python plothmi.py

Open Jupyter notebooks:

    $ jupyter notebook example_plothmi.ipynb

## Get help

Python:

``` python
>>> import usr_sunpy
>>> help(usr_sunpy)
```

Shell:

    $ pydoc usr_sunpy

## Import user modules

e.g. Import functions in <font color=navy>modules/usr_sunpy.py</font>

First append <font color=navy>./modules</font> to <font color=blue>PYTHONPATH</font> in your <font color=navy>~/.bashrc</font>:

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

See examples in <font color=navy>plothmi/</font> .
