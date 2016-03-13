"""
This software comes as Open Source and licensed via AGPL v3.
It was developed under the initiative Copernicus, EFAS operational center @ECMWF (Reading, UK).
"""

from distutils.core import setup

packages_deps = ['numpy>=1.10.1', 'scipy>=0.16.0', 'numexpr>=2.4.6', 'dask[bag]', 'dask[array]', 'toolz']

setup(
    name='interpolator',
    version='1.0',
    packages=['grib_interpolator', 'grib_interpolator.tests'],
    url='https://bitbucket.org/nappodo/grib-interpolator',
    license='AGPL v3. Copernicus @ECMWF',
    author='Domenico Nappo',
    author_email='domenico.nappo@gmail.com',
    description='A python package to read and interpolate GRIB data',
    install_requires=packages_deps,
    keywords="GRIB interpolation Copernicus EFAS ECMWF",
)
