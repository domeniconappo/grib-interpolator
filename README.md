Intro
=====

The grib_interpolator package is a python API to read GRIB files and interpolate values for any target grid.

It's a light and generalized version of part of python code that was developed under the Copernicus EFAS operational center @ECMWF.

Installation
------------

With pip tool:

    pip install https://bitbucket.org/nappodo/grib-interpolator/get/master.zip
    
Or just download code/clone repository and run

    python setup.py
    
Setup will install dependencies as numpy, scipy and numexpr (if not already present in your python/virtualenv). 

You also need to install GRIB API python interface. 
Please refer to GRIB API docs for details: [GRIB API docs](https://software.ecmwf.int/wiki/display/GRIB)

Usage
-----

```python
s = "Python syntax highlighting"
print s
```




