from distutils.core import setup

packages_deps = ['numpy>=1.10.1', 'scipy>=0.16.0', 'numexpr>=2.4.6']

setup(
    name='interpolator',
    version='1.0',
    packages=['grib_interpolator', 'grib_interpolator.tests'],
    url='',
    license='',
    author='dominik',
    author_email='domenico.nappo@gmail.com',
    description='A python package to read and interpolate GRIB data',
    install_requires=packages_deps,
    keywords="GRIB interpolation",
)
