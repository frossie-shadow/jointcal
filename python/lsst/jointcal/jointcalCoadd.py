# See COPYRIGHT file at the top of the source tree.

from __future__ import division, absolute_import, print_function

from lsst.pipe.tasks.makeCoaddTempExp import MakeCoaddTempExpTask
from lsst.pipe.base import Struct


class JointcalCoaddTask(MakeCoaddTempExpTask):

    def getCalExp(self, dataRef, bgSubtracted):
        """!Return one "calexp" calibrated exposure

        @param[in] dataRef        a sensor-level data reference
        @param[in] bgSubtracted   return calexp with background subtracted? If False get the
                                  calexp's background background model and add it to the calexp.
        @return calibrated exposure

        If config.doApplyUberCal, meas_mosaic calibrations will be applied to
        the returned exposure using applyMosaicResults.
        """
        exposure = dataRef.get("calexp")

        if not bgSubtracted:
            background = dataRef.get("calexpBackground")
            mi = exposure.getMaskedImage()
            mi += background.getImage()
            del mi
        if not self.config.doApplyUberCal:
            return exposure
        # if we are here, it means that we have to apply the improved calibrations coming from jointcal
        self.log.info("doApplyUberCal is set - Using jointcal updated calibrations")
        self.applyJointcalResultsExposure(dataRef, calexp=exposure)
        return exposure

    def applyJointcalResultsExposure(self, dataRef, calexp=None):
        """Update an Exposure with the Wcs, from meas_jointcal
        (Calib and flux sacling will be also used later).
        If None, the calexp will be loaded from the dataRef.  Otherwise it is
        updated in-place.
        """
        if calexp is None:
            calexp = dataRef.get("calexp")

        wcsCont = dataRef.get("wcs")
        calexp.setWcs(wcsCont.getWcs())

        return Struct(exposure=calexp)
