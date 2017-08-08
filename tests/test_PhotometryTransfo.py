import numpy as np

import unittest
import lsst.utils.tests

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
        output = np.zeros(1, dtype=float)
        self.transfo.parameterDerivatives(1, 2, self.instFlux, output)
        self.assertEqual(self.instFlux, output[0])


class PhotometryTransfoChebyshevTestCase(PhotometryTransfoTestBase, lsst.utils.tests.TestCase):
    def setUp(self):
        np.random.seed(100)
        self.instFlux = 5.0
        self.point = [2.0, 3.0]
        self.x = (np.random.random(50) - 0.5) * 100
        self.y = (np.random.random(50) - 0.5) * 100
        self.order = 2
        self.transfo = lsst.jointcal.photometryTransfo.PhotometryTransfoChebyshev(self.order)

    def test_apply(self):
        result = self.transfo.apply(1, 2, self.instFlux)
        self.assertEqual(result, self.instFlux)

    def test_offsetParams(self):
        """Test offsetting; note that offsetParams offsets by `-delta`."""
        delta = np.zeros(self.order**2, dtype=float)
        self._test_offsetParams(delta, np.array([[1, 0], [0, 0]]))
        delta = np.array([1, -2, -3, -4], dtype=float)
        self._test_offsetParams(delta, np.array([[0, 2], [3, 4]]))
        # self.transfo.offsetParams(delta)
        # result = self.transfo.apply(self.point[0], self.point[1], self.instFlux)
        # self.assertEqual(result, self.instFlux)
        # delta = np.array([-1, -2, -3, -4], dtype=float)
        # coefficients = np.array([1, 0, 0, 0]) - delta
        # self.transfo.offsetParams(delta.T)
        # result = self.transfo.apply(self.point[0], self.point[1], self.instFlux)
        # import os; print(os.getpid()); import ipdb; ipdb.set_trace();
        # Tx = np.array([CHEBYSHEV_T[i](self.point[0]) for i in range(self.order)], dtype=float)
        # Ty = np.array([CHEBYSHEV_T[i](self.point[1]) for i in range(self.order)], dtype=float)
        # expect = np.dot(Ty[], np.dot(coefficients, Tx[i]))
        # self.assertEqual(result, expect)

    def test_parameterDerivatives(self):
        pass


class MemoryTester(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
