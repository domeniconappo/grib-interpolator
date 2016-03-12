"""
This software comes as Open Source and licensed via AGPL v3.
It was developed under the initiative Copernicus, EFAS operational center @ECMWF (Reading, UK).
"""

import collections

import gribapi


class Step(object):
    def __init__(self, start_step_, end_step_, points_meridian_, input_step_):
        self.start_step = start_step_
        self.end_step = end_step_
        # spatial resolution
        self.resolution = points_meridian_
        # temporal resolution
        self.input_step = input_step_

    def __hash__(self):
        return hash((self.start_step, self.end_step, self.resolution, self.input_step))

    def __eq__(self, other):
        return (self.start_step, self.end_step, self.resolution, self.input_step) == (
               (other.start_step, other.end_step, other.resolution, other.input_step))

    def __str__(self):
        return 's:{} e:{} res:{} step-lenght:{}'.format(self.start_step, self.end_step,
                                                        self.resolution, self.input_step)

    def __lt__(self, other):
        return self.start_step < other.start_step

    def __le__(self, other):
        return self.start_step < other.start_step


class GribGridDetails(object):
    """
    # Managed grid types:
        * regular_gg, regular_ll
        * reduced_ll, reduced_gg (include octahedral grid)
        * rotated_ll, rotated_gg

    """

    keys = (('gridType', 'string'), ('radius', 'double'), ('numberOfValues', 'long'),
            ('Ni', 'long'), ('Nj', 'long'), ('missingValue', 'double'),
            ('longitudeOfFirstGridPointInDegrees', 'double'), ('longitudeOfLastGridPointInDegrees', 'double'),
            ('latitudeOfSouthernPoleInDegrees', 'double'), ('longitudeOfSouthernPoleInDegrees', 'double'),
            ('angleOfRotationInDegrees', 'double'))
    check_for_missing_keys = ('Ni', 'Nj',)

    def __init__(self, gid):

        self._gid = gid
        self._geo_keys = {
            key_: getattr(gribapi, 'grib_get_{}'.format(type_))(gid, key_)
            for key_, type_ in self.keys
            if gribapi.grib_is_defined(gid, key_)
        }
        self._missing_keys = {
            key_: 'MISSING'
            for key_ in self.check_for_missing_keys
            if gribapi.grib_is_missing(gid, key_)
        }
        self._grid_type = self._geo_keys.get('gridType')
        self._points_meridian = self._geo_keys.get('Nj')
        self._missing_value = self._geo_keys.get('missingValue')
        self._grid_id = self._build_id()
        # lazy computation
        self._lats = None
        self._longs = None

        self._grid_details_2nd = None
        self._change_resolution_step = None

    def _build_id(self):
        ni = 'M' if 'Ni' in self._missing_keys else self._geo_keys.get('Ni')
        nj = 'M' if 'Nj' in self._missing_keys else self._geo_keys.get('Nj')
        num_of_values = self._geo_keys.get('numberOfValues')
        long_first = ('%.4f' % (self._geo_keys.get('longitudeOfFirstGridPointInDegrees'),)).rstrip('0').rstrip('.')
        long_last = ('%.4f' % (self._geo_keys.get('longitudeOfLastGridPointInDegrees'),)).rstrip('0').rstrip('.')
        grid_id = '{}${}${}${}${}${}'.format(long_first, long_last, ni, nj, num_of_values, self._grid_type)
        return grid_id

    def set_2nd_resolution(self, grid2nd, step_range_):
        self._grid_details_2nd = grid2nd
        self._change_resolution_step = step_range_
        # change of points along meridian!
        self._points_meridian = grid2nd.num_points_along_meridian

    def get_2nd_resolution(self):
        return self._grid_details_2nd

    def get_change_res_step(self):
        return self._change_resolution_step

    @staticmethod
    def _compute_latlongs(gid):

        lats = gribapi.grib_get_double_array(gid, 'latitudes')
        lons = gribapi.grib_get_double_array(gid, 'longitudes')
        return lats, lons

    @property
    def latlons(self):
        # this method is called only for scipy interpolation
        if self._lats is None:
            self._lats, self._longs = self._compute_latlongs(self._gid)
        return self._lats, self._longs

    @property
    def grid_id(self):
        return self._grid_id

    @property
    def num_points_along_meridian(self):
        return self._points_meridian

    def get(self, geo_key):
        return self._geo_keys[geo_key]

    def __str__(self):
        return str(self._geo_keys)


class Messages(object):

    def __init__(self, values, mv, unit, type_of_level, type_of_step, grid_details, val_2nd=None):
        self.values_first_or_single_res = values
        self.values_second_res = val_2nd or {}
        self.type_of_step = type_of_step
        self.type_of_level = type_of_level
        self.unit = unit
        self.missing_value = mv

        self.grid_details = grid_details
        # order key list to get first step
        self.first_step_range = sorted(self.values_first_or_single_res.keys(), key=lambda k: (int(k.end_step)))[0]

    def append_2nd_res_messages(self, messages):
        # messages is a Messages object from second set at different resolution
        self.grid_details.set_2nd_resolution(messages.grid_details, messages.first_step_range)
        self.values_second_res = messages.first_resolution_values()

    def first_resolution_values(self):
        self.values_first_or_single_res = collections.OrderedDict(sorted(self.values_first_or_single_res.iteritems(),
                                                                         key=lambda (k, v_): (int(k.end_step), v_)))
        return self.values_first_or_single_res

    def second_resolution_values(self):
        self.values_second_res = collections.OrderedDict(sorted(self.values_second_res.iteritems(),
                                                                key=lambda (k, v_): (int(k.end_step), v_)))

        return self.values_second_res

    @property
    def grid_id(self):
        return self.grid_details.grid_id

    @property
    def grid2_id(self):
        return self.grid_details.get_2nd_resolution().grid_id

    @property
    def latlons(self):
        return self.grid_details.latlons

    @property
    def latlons_2nd(self):
        return self.grid_details.get_2nd_resolution().latlons

    def have_resolution_change(self):
        return self.grid_details.get_2nd_resolution() is not None

    def change_resolution_step(self):
        return self.grid_details.get_change_res_step()

    def __len__(self):
        return len(self.values_first_or_single_res) + len(self.values_second_res)
