import numpy as np

import unittest
import lsst.utils.tests

import lsst.afw.geom
import lsst.jointcal.photometryTransfo


CHEBYSHEV_T = [
    lambda x: 1,
    lambda x: x,
    lambda x: 2*x**2 - 1,
    lambda x: (4*x**2 - 3)*x,
    lambda x: (8*x**2 - 8)*x**2 + 1,
    lambda x: ((16*x**2 - 20)*x**2 + 5)*x,
]


class PhotometryTransfoTestBase():
    def _test_offsetParams(self, delta, expect):
        self.transfo.offsetParams(delta)
        self.assertFloatsAlmostEqual(expect, self.transfo.getCoefficients())


class PhotometryTransfoSpatiallyInvariantTestCase(lsst.utils.tests.TestCase):
    def setUp(self):
        self.instFlux = 5.0
        self.transfo = lsst.jointcal.photometryTransfo.PhotometryTransfoSpatiallyInvariant()

    def test_apply(self):
        result = self.transfo.apply(1, 2, self.instFlux)
        self.assertEqual(result, self.instFlux)

    def _test_offsetParams(self, delta, expect):
        self.transfo.offsetParams(delta)
        self.assertFloatsAlmostEqual(expect, self.transfo.apply(1, 2, 1))

    def test_offsetParams(self):
        """Test offsetting; note that offsetParams offsets by `-delta`."""
        delta = np.zeros(1, dtype=float)
        self._test_offsetParams(delta, np.array(1.0))
        delta -= 1
        self._test_offsetParams(delta, np.array(2.0))

    def test_parameterDerivatives(self):
        """Test that the derivative of a spatially invariant transform is always the same."""
        result = self.transfo.parameterDerivatives(1, 2, self.instFlux)
        self.assertEqual(self.instFlux, result)
        result = self.transfo.parameterDerivatives(-5, -100, self.instFlux)
        self.assertEqual(self.instFlux, result)
        transfo = lsst.jointcal.photometryTransfo.PhotometryTransfoSpatiallyInvariant(1000.0)
        result = transfo.parameterDerivatives(1, 2, self.instFlux)
        self.assertEqual(self.instFlux, result)


class PhotometryTransfoChebyshevTestCase(PhotometryTransfoTestBase, lsst.utils.tests.TestCase):
    def setUp(self):
        # np.random.seed(100)
        # self.x = (np.random.random(50) - 0.5) * 100
        # self.y = (np.random.random(50) - 0.5) * 100
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-3, -3), lsst.afw.geom.Point2I(3, 3))
        self.instFlux = 5.0
        self.point = [2.0, 3.0]
        self.degree = 2
        self.transfo = lsst.jointcal.photometryTransfo.PhotometryTransfoChebyshev(self.degree, self.bbox)

    def test_apply(self):
        result = self.transfo.apply(1, 2, self.instFlux)
        self.assertEqual(result, self.instFlux)

    def test_offsetParams(self):
        """Test offsetting; note that offsetParams offsets by `-delta`."""
        delta = np.zeros(self.transfo.getNpar(), dtype=float)
        expect = np.zeros((self.degree+1, self.degree+1), dtype=float)
        expect[0, 0] = 1
        self._test_offsetParams(delta, expect)
        delta[0] = 1
        delta[1] = -2
        delta[2] = -3
        delta[3] = -4
        delta[4] = -5
        delta[5] = -6
        expect[0, 0] = 0
        expect[0, 1] = 2
        expect[0, 2] = 3
        expect[1, 0] = 4
        expect[1, 1] = 5
        expect[2, 0] = 6
        self._test_offsetParams(delta, expect)

    def test_parameterDerivatives(self):
        result = self.transfo.parameterDerivatives(self.point[0], self.point[1], self.instFlux)
        Tx = np.array([CHEBYSHEV_T[i](self.point[0]) for i in range(self.degree+1)], dtype=float)
        Ty = np.array([CHEBYSHEV_T[i](self.point[1]) for i in range(self.degree+1)], dtype=float)
        expect = []
        for j in range(len(Ty)):
            for i in range(0, self.degree-j+1):
                expect.append(Ty[j]*Tx[i]*self.instFlux)
        self.assertFloatsAlmostEqual(np.array(expect), result)


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
