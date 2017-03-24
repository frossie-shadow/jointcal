# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import unittest
import os

from astropy import units as u

import lsst.afw.coord
import lsst.afw.geom
import lsst.utils
import lsst.pex.exceptions

import jointcalTestBase

try:
    data_dir = lsst.utils.getPackageDir('testdata_jointcal')
    os.environ['ASTROMETRY_NET_DATA_DIR'] = os.path.join(data_dir, 'cfht_and_index')
except lsst.pex.exceptions.NotFoundError:
    data_dir = None

# We don't want the absolute astrometry to become significantly worse
# than the single-epoch astrometry (about 0.040").
# This value was empirically determined from the first run of jointcal on
# this data, and will likely vary from survey to survey.
absolute_error = 48e-3*u.arcsecond


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()


class JointcalTestCFHT(jointcalTestBase.JointcalTestBase, lsst.utils.tests.TestCase):

    def setUp(self):
        do_plot = False

        # center of the cfht validation_data catalog
        center = lsst.afw.coord.IcrsCoord(214.884832*lsst.afw.geom.degrees, 52.6622199*lsst.afw.geom.degrees)
        radius = 3*lsst.afw.geom.degrees

        input_dir = os.path.join(data_dir, 'cfht')
        all_visits = [849375, 850587]

        self.setUp_base(center, radius,
                        input_dir=input_dir,
                        all_visits=all_visits,
                        do_plot=do_plot)

    @unittest.skipIf(data_dir is None, "testdata_jointcal not setup")
    def test_jointcalTask_2_visits(self):
        # NOTE: The relative RMS limit was empirically determined from the
        # first run of jointcal on this data. We should always do better than
        # this in the future!
        relative_error = 25e-3*u.arcsecond
        pa1 = 0.019
        metrics = {'collected_astrometry_refStars': 825,
                   'collected_photometry_refStars': 825,
                   'selected_astrometry_refStars': 825,
                   'selected_photometry_refStars': 825,
                   'associated_astrometry_fittedStars': 2269,
                   'associated_photometry_fittedStars': 2269,
                   'selected_astrometry_fittedStars': 1239,
                   'selected_photometry_fittedStars': 1239,
                   'selected_astrometry_ccdImages': 12,
                   'selected_photometry_ccdImages': 12,
                   'astrometry_chisq': 1150.62,
                   'astrometry_ndof': 2550,
                   'photometry_chisq': 13363.6,
                   'photometry_ndof': 1089
                   }

        self._testJointcalTask(2, relative_error, absolute_error, pa1, metrics=metrics)


# TODO: the memory test cases currently fail in jointcal. Filed as DM-6626.
# class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
#     pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
