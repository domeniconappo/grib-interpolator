import os
import numpy as np


if __name__ == '__main__':
    current_dir = os.path.dirname(os.path.abspath(__file__))
    variable = '2t'
    from grib_interpolator.base import Interpolator
    from grib_interpolator.gribreader import GRIBReader

    input_file = '/dataset/test_2013330702/EpsN320-2013063000.grb'
    print 'Opening {}'.format(input_file)
    reader = GRIBReader(input_file, True)
    target_lats = np.load(current_dir + '/grib_interpolator/tests/target_lats.npy')
    target_lons = np.load(current_dir + '/grib_interpolator/tests/target_lons.npy')
    args = {'shortName': variable, 'perturbationNumber': 10}
    messages = reader.select_messages(**args)
    lats, lons = messages.latlons
    grid_details = messages.grid_details
    aux_g, aux_v, aux_g2, aux_v2 = reader.get_gids_for_intertable()

    kwargs = {'source_lons': lons, 'source_lats': lats, 'store': '/dataset/interpolator_intertables/EpsN320',
              'source_grid_details': grid_details, 'gid': aux_g, 'mode': 'nearest', 'method': 'grib'}

    interpolator = Interpolator(**kwargs)
    print 'The intertable {} will be created if not existing yet'.format(interpolator.intertable_path)
    for timestep, values in messages.first_resolution_values().iteritems():
        print 'Interpolating timestep {}'.format(timestep)
        out_file = '{}_{}_{}.npy'.format(timestep.start_step, timestep.end_step, variable)
        out_file = os.path.join('/dataset/interpolator_tests/EpsN320', out_file)
        interpolated_values = interpolator.interpolate(values, target_lons, target_lats)
        print 'Saving result in numpy file {}'.format(out_file)
        np.save(out_file, interpolated_values.data)
    reader.close()

    # COSMO tests
    input_file = '/dataset/cosmo/2012111912_pf10_t2.grb'
    print 'Opening {}'.format(input_file)
    reader = GRIBReader(input_file, True)
    target_lats = np.load(current_dir + '/grib_interpolator/tests/target_lats.npy')
    target_lons = np.load(current_dir + '/grib_interpolator/tests/target_lons.npy')
    args = {'shortName': variable, 'perturbationNumber': 10}
    messages = reader.select_messages(**args)
    lats, lons = messages.latlons
    grid_details = messages.grid_details
    aux_g, aux_v, aux_g2, aux_v2 = reader.get_gids_for_intertable()

    kwargs = {'source_lons': lons, 'source_lats': lats, 'store': '/dataset/interpolator_intertables/cosmo',
              'source_grid_details': grid_details, 'gid': aux_g, 'mode': 'nearest', 'method': 'scipy'}

    interpolator = Interpolator(**kwargs)
    print 'The intertable {} will be created if not existing yet'.format(interpolator.intertable_path)
    for timestep, values in messages.first_resolution_values().iteritems():
        print 'Interpolating timestep {}'.format(timestep)
        out_file = '{}_{}_{}.npy'.format(timestep.start_step, timestep.end_step, variable)
        out_file = os.path.join('/dataset/interpolator_tests/cosmo', out_file)
        interpolated_values = interpolator.interpolate(values, target_lons, target_lats)
        print 'Saving result in numpy file {}'.format(out_file)
        np.save(out_file, interpolated_values.data)
    reader.close()
