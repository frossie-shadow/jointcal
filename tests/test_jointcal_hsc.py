# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

import inspect
import unittest
import os

from astropy import units as u

import lsst.afw.coord
import lsst.afw.geom
import lsst.utils
import lsst.pex.exceptions
from lsst.meas.algorithms import LoadIndexedReferenceObjectsTask

import jointcalTestBase

try:
    data_dir = lsst.utils.getPackageDir('testdata_jointcal')
    os.environ['ASTROMETRY_NET_DATA_DIR'] = os.path.join(data_dir, 'hsc_and_index')
except lsst.pex.exceptions.NotFoundError:
    data_dir = None

# We don't want the absolute astrometry to become significantly worse
# than the single-epoch astrometry (about 0.040").
# This value was empirically determined from the first run of jointcal on
# this data, and will likely vary from survey to survey.
dist_rms_absolute = 53e-3*u.arcsecond


# for MemoryTestCase
def setup_module(module):
    lsst.utils.tests.init()


class JointcalTestHSC(jointcalTestBase.JointcalTestBase, lsst.utils.tests.TestCase):

    def setUp(self):
        do_plot = False

        # center of the hsc validation_data catalog
        center = lsst.afw.coord.IcrsCoord(320.367492*lsst.afw.geom.degrees, 0.3131554*lsst.afw.geom.degrees)
        radius = 5*lsst.afw.geom.degrees

        input_dir = os.path.join(data_dir, 'hsc')
        all_visits = [903334, 903336, 903338, 903342, 903344, 903346, 903986, 903988, 903990, 904010, 904014]

        self.setUp_base(center, radius,
                        input_dir=input_dir,
                        all_visits=all_visits,
                        do_plot=do_plot)

    @unittest.skipIf(data_dir is None, "testdata_jointcal not setup")
    def test_jointcalTask_2_visits(self):
        # NOTE: The relative RMS limit was empirically determined from the
        # first run of jointcal on this data. We should always do better than
        # this in the future!
        dist_rms_relative = 17e-3*u.arcsecond
        pa1 = 0.024
        metrics = {'collected_astrometry_refStars': 2187,
                   'collected_photometry_refStars': 2187,
                   'selected_astrometry_refStars': 2187,
                   'selected_photometry_refStars': 2187,
                   'associated_astrometry_fittedStars': 1151,
                   'associated_photometry_fittedStars': 1151,
                   'selected_astrometry_fittedStars': 770,
                   'selected_photometry_fittedStars': 770,
                   'selected_astrometry_ccdImages': 6,
                   'selected_photometry_ccdImages': 6,
                   'astrometry_chisq': 691.12,
                   'astrometry_ndof': 1858,
                   'photometry_chisq': 3753.82,
                   'photometry_ndof': 504
                   }
        self._testJointcalTask(2, dist_rms_relative, dist_rms_absolute, pa1, metrics=metrics)

    @unittest.skipIf(data_dir is None, "testdata_jointcal not setup")
    def test_jointcalTask_11_visits(self):
        # NOTE: The relative RMS limit was empirically determined from the
        # first run of jointcal on this data. We should always do better than
        # this in the future!
        dist_rms_relative = 17e-3*u.arcsecond
        pa1 = 0.134
        metrics = {'collected_astrometry_refStars': 3649,
                   'collected_photometry_refStars': 3649,
                   'selected_astrometry_refStars': 3649,
                   'selected_photometry_refStars': 3649,
                   'associated_astrometry_fittedStars': 2908,
                   'associated_photometry_fittedStars': 2908,
                   'selected_astrometry_fittedStars': 2203,
                   'selected_photometry_fittedStars': 2203,
                   'selected_astrometry_ccdImages': 33,
                   'selected_photometry_ccdImages': 33,
                   'astrometry_chisq': 7929.656,
                   'astrometry_ndof': 14262,
                   'photometry_chisq': 16773556.5,
                   'photometry_ndof': 6569
                   }
        self._testJointcalTask(11, dist_rms_relative, dist_rms_absolute, pa1, metrics=metrics)

    @unittest.skipIf(data_dir is None, "testdata_jointcal not setup")
    def testJointcalTask_2_visits_no_astrometry(self):
        """Test turning off fitting astrometry."""
        pa1 = 0.024
        metrics = {'collected_photometry_refStars': 2187,
                   'selected_photometry_refStars': 2187,
                   'associated_photometry_fittedStars': 1151,
                   'selected_photometry_fittedStars': 770,
                   'selected_photometry_ccdImages': 6,
                   'photometry_chisq': 3753.82,
                   'photometry_ndof': 504
                   }

        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.doAstrometry = False
        self.jointcalStatistics.do_astrometry = False

        caller = inspect.stack()[0][3]  # NOTE: could be inspect.stack()[0].function in py3.5
        result = self._runJointcalTask(2, caller, metrics=metrics)
        data_refs = result.resultList[0].result.dataRefs
        oldWcsList = result.resultList[0].result.oldWcsList
        rms_result = self.jointcalStatistics.compute_rms(data_refs, self.reference)

        if self.do_plot:
            self._plotJointcalTask(data_refs, oldWcsList, caller)

        self.assertIsNone(rms_result.dist_relative)
        self.assertIsNone(rms_result.dist_absolute)
        self.assertLess(rms_result.pa1, pa1)

        for data_ref in data_refs:
            wcs = data_ref.get('wcs').getWcs()
            self.assertIsNone(wcs)

    @unittest.skipIf(data_dir is None, "testdata_jointcal not setup")
    def testJointcalTask_2_visits_no_photometry(self):
        """Test turning off fitting photometry."""
        dist_rms_relative = 17e-3*u.arcsecond
        metrics = {'collected_astrometry_refStars': 2187,
                   'selected_astrometry_refStars': 2187,
                   'associated_astrometry_fittedStars': 1151,
                   'selected_astrometry_fittedStars': 770,
                   'selected_astrometry_ccdImages': 6,
                   'astrometry_chisq': 691.1210,
                   'astrometry_ndof': 1858,
                   }

        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.doPhotometry = False
        self.jointcalStatistics.do_photometry = False

        caller = inspect.stack()[0][3]  # NOTE: could be inspect.stack()[0].function in py3.5
        result = self._runJointcalTask(2, caller, metrics=metrics)
        data_refs = result.resultList[0].result.dataRefs
        oldWcsList = result.resultList[0].result.oldWcsList
        rms_result = self.jointcalStatistics.compute_rms(data_refs, self.reference)

        if self.do_plot:
            self._plotJointcalTask(data_refs, oldWcsList, caller)

        self.assertLess(rms_result.dist_relative, dist_rms_relative)
        self.assertLess(rms_result.dist_absolute, dist_rms_absolute)
        self.assertIsNone(rms_result.pa1)

        for data_ref in data_refs:
            calib = data_ref.get('wcs').getCalib()
            blank_calib = lsst.afw.image.Calib()
            self.assertEqual(calib, blank_calib)

    @unittest.skipIf(data_dir is None, "testdata_jointcal not setup")
    def test_jointcalTask_2_visits_gaia_refcat(self):
        self.config = lsst.jointcal.jointcal.JointcalConfig()
        self.config.astrometryRefObjLoader.retarget(LoadIndexedReferenceObjectsTask)

        test_config = os.path.join(lsst.utils.getPackageDir('jointcal'), 'tests/config/hsc-config.py')
        self.other_args.extend(['--configfile', test_config])
        dist_rms_relative = 17e-3*u.arcsecond
        # NOTE: PA1 is slightly different here, because the number of SDSS
        # cross-matches within 0.1" goes down after we apply the GAIA-fit WCS.
        pa1 = 0.02405
        metrics = {'collected_astrometry_refStars': 1425,
                   'collected_photometry_refStars': 2187,
                   'selected_astrometry_refStars': 1425,
                   'selected_photometry_refStars': 2187,
                   'associated_astrometry_fittedStars': 1151,
                   'associated_photometry_fittedStars': 1151,
                   'selected_astrometry_fittedStars': 645,
                   'selected_photometry_fittedStars': 770,
                   'selected_astrometry_ccdImages': 6,
                   'selected_photometry_ccdImages': 6,
                   'astrometry_chisq': 435.01995,
                   'astrometry_ndof': 1412,
                   'photometry_chisq': 3753.82,
                   'photometry_ndof': 504
                   }
        # NOTE: The astrometry/photometry tests are computed using the a.net SDSS refcat,
        # so the absolute astrometry RMS will be larger (because GAIA is better, so
        # comparing against SDSS will look "worse").
        dist_rms_absolute = 56e-3*u.arcsecond
        self._testJointcalTask(2, dist_rms_relative, dist_rms_absolute, pa1, metrics=metrics)

# TODO: the memory test cases currently fail in jointcal. Filed as DM-6626.
# class MyMemoryTestCase(lsst.utils.tests.MemoryTestCase):
#     pass

if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
