Intro
=====

The grib_interpolator package is a python API to read GRIB files and interpolate values for any target grid.

It's a light and generalized version of part of python code that was developed under the Copernicus EFAS operational center @ECMWF.

It handles GRIB files up to two resolutions messages, and interpolates values using multiple methods with Interpolator class.

Available methods are:

+ mode = grib method = nearest Nearest neighbour n=1 with GRIB API
+ mode = grib method = invdist Inverse Distance method n=4 with GRIB API
+ mode = scipy method = nearest Nearest neighbour n=1 with scipy.kdtree
+ mode = scipy method = invdist Inverse Distance method n=4 with scipy.kdtree

Known problems
--------------

* GRIB API doesn't interpolate rotated grids yet so you have to use scipy methods. In this case, target grid must be regular
to avoid problems at extreme regions or grid borders.

* Missing values for input GRIBs are not properly handled yet.



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

The API provides an easy way to read messages from a GRIB file and to interpolate values against a target grid.
Interpolation process can take really long so calculated indexes and weights are being stored 
in a pickled numpy array under a folder of choice and automatically selected for future interpolations 
(filenaming is based on source grid metadata). Use one folder per each target grid you have, 
otherwise intertables will be overwritten and you can have unexpected results!

If your target grid is rotated, include the flag _rotated_target_ when instatiate Interpolator.
There isn't a flag for rotated source grids since you should never have such input grids, as coordinates are always transformed to 
regular grid when read through GRIB API.

```python
import os
import numpy as np

from grib_interpolator import Interpolator
from grib_interpolator import GRIBReader


if __name__ == '__main__':
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # shortName of variable to extract, as stored in GRIB messages
    variable = '2t'

    input_file = '/dataset/test_2013330702/EpsN320-2013063000.grb'
    print 'Opening {}'.format(input_file)
    # Opening input file, using indexes if present in GRIB file
    reader = GRIBReader(input_file, indexes=('shortName', 'perturbationNumber'))
    # you can filter messages with any key,
    # as startStep, endStep and so on. Filters can also be lists of values to select

    # >>> messages = reader.select_messages(shortName=['2t', 'e'], perturbationNumber=[8, 10])

    # Filters can also be functions (very useful to select messages between a range)
    start_step = lambda s: s >= 6
    end_step = lambda s: s <= 120
    messages = reader.select_messages(shortName=variable, perturbationNumber=10, startStep=start_step, endStep=end_step)

    # get coordinates from GRIB
    lats, lons = messages.latlons
    # grid_details is a GridDetails object with some geodetic metadata describing grid, used in interpolation
    grid_details = messages.grid_details

    # we need these gids (GRIB message ids) from reader because they are required from GRIB API 'find nearest' methods.

    aux_g, aux_v, aux_g2, aux_v2 = reader.get_gids_for_intertable()

    # Intertables will be saved/loaded using this folder.
    # The intertable is a numpy binary file with indexes and weights that will be reused
    # for future interpolations
    # Important Note: Use one folder per target grid,
    # otherwise intertables will be overwritten as file naming is based only on source grid metadata.
    store = '/dataset/interpolator_intertables/europe_5km'

    # Loading target grid. In this example, files are pickled numpy array representing Europe grid 5Km
    target_lats = np.load(current_dir + '/grib_interpolator/tests/target_lats.npy')
    target_lons = np.load(current_dir + '/grib_interpolator/tests/target_lons.npy')

    # mode can be 'nearest', 'invdist'. method can be 'grib' or 'scipy'
    interpolator = Interpolator(source_lons=lons, source_lats=lats, source_grid_details=grid_details,
                                gid=aux_g, mode='nearest', method='grib', store=store)

    # Note: intertable creation can take from minutes to several hours or days,
    # depending on source and target sizes and from CPU speed. Once intertable is created, the whole process
    # will last a few seconds
    # In this example result will be saved as numpy binary files
    print 'The intertable {} will be created if not existing yet'.format(interpolator.intertable_path)
    for timestep, values in messages.first_resolution_values().iteritems():
        print 'Interpolating timestep {}'.format(timestep)
        out_file = '{}_{}_{}.npy'.format(timestep.start_step, timestep.end_step, variable)
        out_file = os.path.join('/dataset/interpolator_tests/EpsN320', out_file)
        interpolated_values = interpolator.interpolate(values, target_lons, target_lats)
        print 'Saving result in numpy file {}'.format(out_file)
        np.save(out_file, interpolated_values.data)  # result is a MaskedArray so we only save real array
    reader.close()

    ################################
    # Same test for COSMO input data
    input_file = '/dataset/cosmo/2012111912_pf10_t2.grb'
    print 'Opening {}'.format(input_file)
    reader = GRIBReader(input_file, indexes=('shortName', 'perturbationNumber'))
    messages = reader.select_messages(shortName=variable, perturbationNumber=10)
    lats, lons = messages.latlons
    grid_details = messages.grid_details

    interpolator = Interpolator(source_lons=lons, source_lats=lats, source_grid_details=grid_details,
                                mode='nearest', method='scipy', store='/dataset/interpolator_intertables/cosmo')

    for timestep, values in messages.first_resolution_values().iteritems():
        print 'Interpolating timestep {}'.format(timestep)
        out_file = '{}_{}_{}.npy'.format(timestep.start_step, timestep.end_step, variable)
        out_file = os.path.join('/dataset/interpolator_tests/cosmo', out_file)
        interpolated_values = interpolator.interpolate(values, target_lons, target_lats)
        np.save(out_file, interpolated_values.data)
    reader.close()

```




