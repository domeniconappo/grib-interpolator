import unittest
import os
import numpy as np

from grib_interpolator.base import Interpolator
from grib_interpolator.gribreader import GRIBReader

current_dir = os.path.dirname(__file__)


class TestReader(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.input_file = '/dataset/test_2013330702/EpsN320-2013063000.grb'
        cls.reader = GRIBReader(cls.input_file, True)

    @classmethod
    def tearDownClass(cls):
        cls.reader.close()

    def test_messages(self):
        args = {'shortName': '2t', 'perturbationNumber': 10}
        messages = self.reader.select_messages(**args)
        self.assertEqual(len(messages), 41)
        lat, lons = messages.latlons
        self.assertEqual(lat.shape, (21489, ))
        self.assertEqual(messages.grid_details.get('Nj'), 135)


class TestInterpolator(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.input_file = '/dataset/test_2013330702/EpsN320-2013063000.grb'
        cls.reader = GRIBReader(cls.input_file, True)

        cls.target_lats = np.load(current_dir + '/target_lats.npy')
        cls.target_lons = np.load(current_dir + '/target_lons.npy')

    @classmethod
    def tearDownClass(cls):
        cls.reader.close()

    def test_interpolator_scipy(self):
        args = {'shortName': '2t', 'perturbationNumber': 10, 'endStep': 6}
        messages = self.reader.select_messages(**args)
        lats, lons = messages.latlons
        grid_details = messages.grid_details
        aux_g, aux_v, aux_g2, aux_v2 = self.reader.get_gids_for_intertable()

        kwargs = {'source_lons': lons, 'source_lats': lats,
                  'source_grid_details': grid_details, 'gid': aux_g}

        interpolator = Interpolator(**kwargs)
        result, weights, indexes = interpolator.interpolate(aux_v, self.target_lons, self.target_lats)
        self.assertEqual(result.size, self.target_lons.size)

    def test_interpolator_gribapi(self):
        args = {'shortName': '2t', 'perturbationNumber': 10, 'endStep': 6}
        messages = self.reader.select_messages(**args)
        lats, lons = messages.latlons
        grid_details = messages.grid_details
        aux_g, aux_v, aux_g2, aux_v2 = self.reader.get_gids_for_intertable()

        kwargs = {'source_lons': lons, 'source_lats': lats,
                  'source_grid_details': grid_details, 'gid': aux_g, 'mode': 'invdist', 'method': 'grib'}

        interpolator = Interpolator(**kwargs)
        result, weights, indexes = interpolator.interpolate(aux_v, self.target_lons, self.target_lats)
        self.assertEqual(result.size, self.target_lons.size)


