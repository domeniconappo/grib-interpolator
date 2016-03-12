import os
import abc

import numpy as np

from grib_interpolator.griblib import grib_nearest, grib_invdist
from grib_interpolator.scipylib import InverseDistance
from grib_interpolator.utils import mask_it


class _Interpolator(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, source_lons, source_lats, source_grid_details, source_mv, target_mv,
                 rotated_target=False, parallel=True, gid=1):
        self.source_lons = source_lons
        self.source_lats = source_lats
        self.grid_details = source_grid_details
        self.source_mv = source_mv
        self.target_mv = target_mv
        self.parallel = parallel
        self.rotated_target = rotated_target
        self.gid = gid

    @abc.abstractmethod
    def interpolate(self, source_values, target_lons, target_lats):
        raise NotImplementedError()

    @abc.abstractmethod
    def interpolate_with_table(self, intertable, source_values, target_lons, target_lats):
        raise NotImplementedError()


class ScipyNearest(_Interpolator):

    def interpolate_with_table(self, intertable, source_values, target_lons, target_lats):
        indexes = intertable['indexes']
        z = source_values
        z = np.append(z, self.target_mv)
        result = z[indexes]
        return result.reshape(target_lons.shape)

    def __init__(self, *args, **kwargs):
        super(ScipyNearest, self).__init__(*args, **kwargs)
        self.scipy_interpolator = InverseDistance(self.source_lons, self.source_lats,
                                                  self.grid_details, nnear=1, target_mv=self.target_mv,
                                                  source_mv=self.source_mv, rotated_target=self.rotated_target,
                                                  parallel=self.parallel)

    def interpolate(self, source_values, target_lons, target_lats):
        result, indexes, weights = self.scipy_interpolator.interpolate(source_values, target_lons, target_lats)
        intertable = np.rec.fromarrays((indexes, weights), names=('indexes', 'coeffs'))
        result = result.reshape(target_lons.shape)
        return result, intertable


class ScipyInvdist(ScipyNearest):

    def interpolate_with_table(self, intertable, source_values, target_lons, target_lats):
        indexes = intertable['indexes']
        coeffs = intertable['coeffs']
        z = source_values
        z = np.append(z, self.target_mv)
        result = np.einsum('ij,ij->i', coeffs, z[indexes])
        return result.reshape(target_lons.shape)

    def __init__(self, *args, **kwargs):
        super(ScipyNearest, self).__init__(*args, **kwargs)
        self.scipy_interpolator = InverseDistance(self.source_lons, self.source_lats,
                                                  self.grid_details, nnear=4, target_mv=self.target_mv,
                                                  source_mv=self.source_mv, rotated_target=self.rotated_target,
                                                  parallel=self.parallel)


class GribNearest(_Interpolator):

    def interpolate_with_table(self, intertable, source_values, target_lons, target_lats):
        xs, ys, idxs = intertable[0], intertable[1], intertable[2]
        result = np.empty(target_lons.shape)
        result.fill(self.target_mv)
        result = mask_it(result, self.target_mv)
        result[xs, ys] = source_values[idxs]
        return result

    def __init__(self, *args, **kwargs):
        super(GribNearest, self).__init__(*args, **kwargs)
        self.gid = kwargs.get('gid', -1)

    def interpolate(self, source_values, target_lons, target_lats):
        result = np.empty(target_lons.shape)
        result.fill(self.target_mv)
        result = mask_it(result, self.target_mv)
        xs, ys, idxs = grib_nearest(self.gid, target_lats, target_lons, self.target_mv)
        intertable = np.asarray([xs, ys, idxs], dtype=np.int32)
        result[xs, ys] = source_values[idxs]
        return result, intertable


class GribInvdist(GribNearest):

    def interpolate_with_table(self, intertable, source_values, target_lons, target_lats):
        indexes = intertable['indexes']  # first two arrays of this group are target xs and ys indexes
        coeffs = intertable['coeffs']
        v = source_values
        result = np.empty(target_lons.shape)
        result.fill(self.target_mv)
        result = mask_it(result, self.target_mv)
        xs, ys, idxs1, idxs2, idxs3, idxs4, coeffs1, coeffs2, coeffs3, coeffs4 = indexes[0], indexes[1], indexes[2], indexes[3], indexes[4], indexes[5], coeffs[0], coeffs[1], coeffs[2], coeffs[3]
        result[xs, ys] = v[idxs1] * coeffs1 + v[idxs2] * coeffs2 + v[idxs3] * coeffs3 + v[idxs4] * coeffs4
        return result

    def interpolate(self, source_values, target_lons, target_lats):
        v = source_values
        result = np.empty(target_lons.shape)
        result.fill(self.target_mv)
        result = mask_it(result, self.target_mv)
        xs, ys, idxs1, idxs2, idxs3, idxs4, coeffs1, coeffs2, coeffs3, coeffs4 = grib_invdist(self.gid, target_lats,
                                                                                              target_lons,
                                                                                              self.target_mv)
        indexes = np.asarray([xs, ys, idxs1, idxs2, idxs3, idxs4], dtype=np.int32)
        coeffs = np.asarray([coeffs1, coeffs2, coeffs3, coeffs4, np.zeros(coeffs1.shape), np.zeros(coeffs1.shape)])
        intertable = np.rec.fromarrays((indexes, coeffs), names=('indexes', 'coeffs'))
        result[xs, ys] = v[idxs1] * coeffs1 + v[idxs2] * coeffs2 + v[idxs3] * coeffs3 + v[idxs4] * coeffs4
        return result, intertable


_Interpolator.register(ScipyNearest)
_Interpolator.register(ScipyInvdist)
_Interpolator.register(GribNearest)
_Interpolator.register(GribInvdist)


class Interpolator(object):

    scipy_nearest = ScipyNearest
    scipy_invdist = ScipyInvdist
    grib_nearest = GribNearest
    grib_invdist = GribInvdist

    def __init__(self, source_lats, source_lons, source_grid_details, **kwargs):
        self.source_lons = source_lons
        self.source_lats = source_lats
        self.grid_details = source_grid_details
        self._mode = kwargs.get('mode', 'nearest')
        self._method = kwargs.get('method', 'scipy')
        self.target_mv = kwargs.get('target_mv', np.nan)
        self.source_mv = kwargs.get('source_mv', np.nan)
        self.rotated_target = kwargs.get('rotated_target', False)
        self.parallel = kwargs.get('parallel', True)
        self.gid = kwargs.get('gid', -1)  # id of grib message (comes from reader)
        self.interpolation_method = '{}_{}'.format(self._method, self._mode)
        self.intertable_filename = '{}_{}.npy'.format(self.grid_details.grid_id.replace('$', '_'),
                                                      self.interpolation_method)
        self.intertables_dir = kwargs.get('store', './')
        if not os.path.exists(self.intertables_dir):
            os.makedirs(self.intertables_dir)
        self.intertable_path = os.path.join(self.intertables_dir, self.intertable_filename)
        self._interpolator = getattr(self, self.interpolation_method)(source_lons, source_lats,
                                                                      self.grid_details,
                                                                      self.source_mv, self.target_mv,
                                                                      self.rotated_target, self.parallel)

    def interpolate(self, source_values, target_lons, target_lats):

        if not os.path.exists(self.intertable_path):
            print 'Creating intertable {}'.format(self.intertable_path)
            result, intertable = self._interpolator.interpolate(source_values, target_lons, target_lats)
            np.save(self.intertable_path, intertable)
            return result
        intertable = np.load(self.intertable_path)
        return self._interpolator.interpolate_with_table(intertable, source_values, target_lons, target_lats)

