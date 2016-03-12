"""
This software comes as Open Source and licensed via AGPL v3.
It was developed under the initiative Copernicus, EFAS operational center @ECMWF (Reading, UK).
"""

from __future__ import division

from itertools import izip
from math import radians
from sys import stdout

import numexpr as ne
import numpy as np
from scipy.spatial import cKDTree as KDTree

from grib_interpolator.utils import mask_it, progress_step_and_backchar, empty, now_string

np.seterr(all='ignore')


class InverseDistance(object):
    """
    http://docs.scipy.org/doc/scipy/reference/spatial.html
    """

    def __init__(self, sourcelons, sourcelats, grid_details, nnear, target_mv, source_mv,
                 rotated_target=False, parallel=False):
        stdout.write('Start scipy interpolation: {}\n'.format(now_string()))
        self.geodetic_info = grid_details
        self.target_grid_is_rotated = rotated_target
        self.njobs = 1 if not parallel else -1
        self.nnear = nnear
        # we receive rotated coords from GRIB_API iterator before 1.14.3
        x, y, zz = self.to_3d(sourcelons, sourcelats)
        source_locations = np.vstack((x.ravel(), y.ravel(), zz.ravel())).T
        self.source_locations = source_locations

        stdout.write('Building KDTree...\n')
        self.tree = KDTree(source_locations, leafsize=30)  # build the tree

        self._mv_target = target_mv
        self._mv_source = source_mv
        distances, indexes = self.tree.query(source_locations, k=2, n_jobs=self.njobs)
        self.min_upper_bound = np.max(distances) + np.max(distances) * 4 / self.geodetic_info.get('Nj')
        stdout.write('Skipping neighbors at distance > {}\n'.format(self.min_upper_bound))

    def interpolate(self, source_values, target_lons, target_lats):
        # Target coordinates  HAVE to be rotated coords in case GRIB grid is rotated
        # Examples of target rotated coords are COSMO lat/lon/dem PCRASTER maps
        x, y, z = self.to_3d(target_lons, target_lats, to_regular=self.target_grid_is_rotated)
        target_locations = np.vstack((x.ravel(), y.ravel(), z.ravel())).T

        stdout.write('Finding indexes for nearest neighbour k={}\n'.format(self.nnear))

        distances, indexes = self.tree.query(target_locations, k=self.nnear, n_jobs=self.njobs)

        if self.nnear == 1:
            # return distances, distances, indexes
            result, indexes = self._build_nn(source_values, distances, indexes)
            weights = distances
        else:
            # return distances, weights, indexes
            result, weights, indexes = self._build_weights(source_values, distances, indexes, self.nnear)

        stdout.write('End scipy interpolation: {}\n'.format(now_string()))
        return result, indexes, weights

    def to_3d(self, lons, lats, rotate=False, to_regular=False):

        lons = np.radians(lons)
        lats = np.radians(lats)
        x_formula = 'cos(lons) * cos(lats)'
        y_formula = 'sin(lons) * cos(lats)'
        z_formula = 'sin(lats)'

        if to_regular:
            teta = - radians((90 + self.geodetic_info.get('latitudeOfSouthernPoleInDegrees')))
            fi = - radians(self.geodetic_info.get('longitudeOfSouthernPoleInDegrees'))
            x = ne.evaluate('(cos(teta) * cos(fi) * ({x})) + (sin(fi)  * ({y})) + (sin(teta) * cos(fi) * ({z}))'.format(x=x_formula, y=y_formula, z=z_formula))
            y = ne.evaluate('(cos(teta) * sin(fi) * ({x})) + (cos(fi)  * ({y})) - (sin(teta) * sin(fi) * ({z}))'.format(x=x_formula, y=y_formula, z=z_formula))
            z = ne.evaluate('(-sin(teta) * ({x})) + (cos(teta) * ({z}))'.format(x=x_formula, z=z_formula))
        elif rotate:
            teta = radians((90 + self.geodetic_info.get('latitudeOfSouthernPoleInDegrees')))
            fi = radians(self.geodetic_info.get('longitudeOfSouthernPoleInDegrees'))
            x = ne.evaluate('(cos(teta) * cos(fi) * ({x})) + (cos(teta) * sin(fi) * ({y})) + (sin(teta) * ({z}))'.format(x=x_formula, y=y_formula, z=z_formula))
            y = ne.evaluate('(-sin(fi) * ({x})) + (cos(fi) * ({y}))'.format(x=x_formula, y=y_formula))
            z = ne.evaluate('(-sin(teta) * cos(fi) * ({x})) - (sin(teta) * sin(fi) * ({y})) + (cos(teta) * ({z}))'.format(x=x_formula, y=y_formula, z=z_formula))
        else:
            r = self.geodetic_info.get('radius')
            x = ne.evaluate('r * {x}'.format(x=x_formula))
            y = ne.evaluate('r * {y}'.format(y=y_formula))
            z = ne.evaluate('r * {z}'.format(z=z_formula))
        return x, y, z

    def _build_nn(self, z, distances, indexes):
        z = mask_it(z, self._mv_source)
        # TODO probably we don't need to mask but just an empty array
        result = mask_it(np.empty((len(distances),) + np.shape(z[0])), self._mv_target, 1)
        jinterpol = 0
        num_cells = result.size
        back_char, progress_step = progress_step_and_backchar(num_cells)

        stdout.write('{}Building coeffs: 0/{} [outs: 0] (0%)'.format(back_char, num_cells))
        stdout.flush()
        idxs = empty((len(indexes),), fill_value=z.size, dtype=np.int64)
        # wsum will be saved in intertable
        outs = 0
        for dist, ix in izip(distances, indexes):
            if jinterpol % progress_step == 0:
                stdout.write('{}Building coeffs: {}/{} [outs: {}] ({:.2f}%)'.format(back_char, jinterpol, num_cells, outs, jinterpol * 100. / num_cells))
                stdout.flush()
            if dist <= self.min_upper_bound:
                wz = z[ix]
                idxs[jinterpol] = ix
            else:
                outs += 1
                wz = self._mv_target
            result[jinterpol] = wz
            jinterpol += 1
        stdout.write('{}{:>100}'.format(back_char, ' '))
        stdout.write('{}Building coeffs: {}/{} [outs: {}] (100%)\n'.format(back_char, jinterpol, num_cells, outs))
        stdout.flush()
        return result, idxs

    def _build_weights(self, z, distances, indexes, nnear):

        # TODO CHECK: maybe we don't need to mask here
        z = mask_it(z, self._mv_source)
        # no intertable found for inverse distance nnear = 8
        # TODO CHECK if we need mask here (maybe just need an empty array)
        result = mask_it(np.empty((len(distances),) + np.shape(z[0])), self._mv_target, 1)
        jinterpol = 0
        num_cells = result.size

        back_char, progress_step = progress_step_and_backchar(num_cells)

        stdout.write('{}Building coeffs: 0/{} [outs: 0] (0%)'.format(back_char, num_cells))
        stdout.flush()
        # weights will be saved in intertable along with indexes
        weights = empty((len(distances),) + (nnear,))
        idxs = empty((len(indexes),) + (nnear,), fill_value=z.size, dtype=int)
        empty_array = empty(z[0].shape, self._mv_target)
        outs = 0
        for dist, ix in izip(distances, indexes):
            if jinterpol % progress_step == 0:
                stdout.write('{}Building coeffs: {}/{} [outs: {}] ({:.2f}%)'.format(back_char, jinterpol, num_cells, outs, jinterpol * 100. / num_cells))
                stdout.flush()
            if dist[0] <= 1e-10:
                wz = z[ix[0]]  # take exactly the point, weight = 1
                idxs[jinterpol] = ix
                weights[jinterpol] = np.array([1., 0., 0., 0.])
            elif dist[0] <= self.min_upper_bound:
                w = ne.evaluate('1 / dist ** 2')
                sums = ne.evaluate('sum(w)')
                ne.evaluate('w/sums', out=w)
                wz = np.dot(w, z[ix])  # weighted values (result)
                weights[jinterpol] = w
                idxs[jinterpol] = ix
            else:
                outs += 1
                weights[jinterpol] = np.array([1., 0., 0., 0.])
                wz = empty_array
            result[jinterpol] = wz
            jinterpol += 1
        stdout.write('{}{:>100}'.format(back_char, ' '))
        stdout.write('{}Building coeffs: {}/{} [outs: {}] (100%)\n'.format(back_char, jinterpol, num_cells, outs))
        stdout.flush()
        return result, weights, idxs
