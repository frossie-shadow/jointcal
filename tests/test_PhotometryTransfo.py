import numpy as np

import unittest
import lsst.utils.tests

import lsst.jointcal.photometryTransfo


class PhotometryTransfoSpatiallyInvariantTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.instFlux = 3.0
        self.transfo = lsst.jointcal.photometryTransfo.PhotometryTransfoSpatiallyInvariant(1)

    def test_apply(self):
        result = self.transfo.apply(1, 2, self.instFlux)
        self.assertEqual(result, self.instFlux)

    def test_offsetParams(self):
        """Test offsetting; note that offsetParams offsets by `-delta`."""
        delta = np.array(-1)
        self.transfo.offsetParams(delta)
        result = self.transfo.apply(1, 2, self.instFlux)
        self.assertEqual(result, self.instFlux*2)
        delta = np.array(2)
        self.transfo.offsetParams(delta)
        result = self.transfo.apply(1, 2, self.instFlux)
        self.assertEqual(result, 0)


class PhotometryTransfoChebyshevTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.transfo = lsst.jointcal.photometryTransfo.PhotometryTransfoChebyshev(2)

    def test_apply(self):
        instFlux = 3.0
        result = self.transfo.apply(1, 2, instFlux)
        self.assertEqual(result, instFlux)

    def test_offsetParams(self):
        """Test offsetting; note that offsetParams offsets by `-delta`."""
        delta = np.array((-1, -1, -1, -1))
        self.transfo.offsetParams(delta)
        result = self.transfo.apply(1, 2, self.instFlux)
        self.assertEqual(result, self.instFlux*2)
        delta = np.array(2)
        self.transfo.offsetParams(delta)
        result = self.transfo.apply(1, 2, self.instFlux)
        self.assertEqual(result, 0)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
