# SunPy scripts & examples




[ Blog ] **NEW!** https://lydiazly.coding.me/python-notes

[ Update ]
> *2018-10-13*&emsp;Upgrade to SunPy 0.9.3<br>
> *2018-07-12*&emsp;Fix links in notebooks.<br>
> *2018-07-11*&emsp;Add functions in [usr_sunpy](https://coding.net/u/lydiazly/p/scripts-sunpy/git/tree/master/modules)<br>
> &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; Add examples of Fido.<br>
> &emsp;&emsp;&emsp;&emsp;&emsp;&emsp; Move files to blog.<br>
> *2018-07-01*&emsp;Uploaded pdf version; modified doc
<br>

Please let me know if you have any questions: lydiazly@nju.edu.cn

[ Tested ]

Ubuntu 16.04 LTS

Date|Python|Sunpy|NumPy|SciPy|Matplotlib|Astropy|Pandas
:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:
2018-10-13|3.6.6|**0.9.3**|1.15.2|1.1.0|3.0.0|3.0.4|0.23.4
2018-07-08|3.6.5|0.9.0|1.14.5|1.1.0|2.2.2|3.0.3|0.23.0
2018-06-28|2.7.12|0.9.0|1.14.4|1.1.0|2.2.2|2.0.6|0.22.0

---

## Download

    $ git clone https://git.coding.net/lydiazly/scripts-sunpy.git

---

## Install [SunPy](http://sunpy.org)

( See also: [sunpy_troubleshooting.txt](https://coding.net/u/lydiazly/p/scripts-sunpy/git/blob/master/sunpy_troubleshooting.txt) )

SunPy Installation Information:

``` python
>>> import sunpy
>>> sunpy.util.system_info()
```

---

## Examples

* plothmi/example_plothmi.ipynb<br>
[
[Jupyter Notebook](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/plothmi/example_plothmi.ipynb)
|
[HTML](https://lydiazly.coding.me/python-notes/_pages/example_plothmi.html)
|
[PDF](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/plothmi/example_plothmi.pdf)
]

* plothmi/example_projection.ipynb<br>
[
[Jupyter Notebook](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/plothmi/example_projection.ipynb)
|
[HTML](https://lydiazly.coding.me/python-notes/_pages/example_projection.html)
|
[PDF](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/plothmi/example_projection.pdf)
]

&emsp;&emsp;**+++ NEW! +++**

* fido/fido.ipynb<br>
[
[Jupyter Notebook](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/fido/fido.ipynb)
|
[HTML](https://lydiazly.coding.me/python-notes/_pages/fido.html)
|
[PDF](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/fido/fido.pdf)
]

* fido/jsoc.ipynb<br>
[
[Jupyter Notebook](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/fido/jsoc.ipynb)
|
[HTML](https://lydiazly.coding.me/python-notes/_pages/jsoc.html)
|
[PDF](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/fido/jsoc.pdf)
]

Run scripts:

    $ python plothmi.py

Open Jupyter notebooks:

    $ jupyter notebook example_plothmi.ipynb

---

## User modules

### Get help

**+++ NEW! +++**

Python:

``` python
>>> import usr_sunpy
>>> help(usr_sunpy.plot)
```

Shell:

    $ pydoc usr_sunpy.plot

See this 
<a href="https://lydiazly.coding.me/python-notes/usr_sunpy.html" target="_blank">
documentation preview
</a>
.

### Import

First append&ensp;[modules/](https://coding.net/u/lydiazly/p/scripts-sunpy/git/tree/master/modules)&ensp;to&ensp;**PYTHONPATH**&ensp;in your&ensp;*~/.bashrc* :

``` sh
export PYTHONPATH=$PYTHONPATH:<your_path>/scripts-sunpy/modules
```

Then import in python:

**+++ NEW! +++**

``` python
>>> import usr_sunpy
>>> import usr_sunpy.plot
>>> from usr_sunpy.plot import <function_name>
```

Or add the path temporarily in python:

``` python
>>> import sys
>>> sys.path.append('../modules')  # If current location is plothmi/
>>> import usr_sunpy
```

See examples in [plothmi/](https://coding.net/u/lydiazly/p/scripts-sunpy/git/tree/master/plothmi) ,&ensp;[fido/jsoc.ipynb](https://coding.net/u/lydiazly/p/scripts-sunpy/git/raw/master/fido/jsoc.ipynb).
