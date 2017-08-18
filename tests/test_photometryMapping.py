import numpy as np

import unittest
import lsst.utils.tests

import lsst.afw.geom
import lsst.jointcal.photometryMappings
import lsst.jointcal.photometryTransfo
import lsst.jointcal.star


CHEBYSHEV_T = [
    lambda x: 1,
    lambda x: x,
    lambda x: 2*x**2 - 1,
    lambda x: (4*x**2 - 3)*x,
    lambda x: (8*x**2 - 8)*x**2 + 1,
    lambda x: ((16*x**2 - 20)*x**2 + 5)*x,
]


class PhotometryMappingTestBase():
    def setUp(self):
        self.instFlux = 5.0
        self.instFluxErr = 2.0

        baseStar0 = lsst.jointcal.star.BaseStar(0, 0, 1, 2)
        self.star0 = lsst.jointcal.star.MeasuredStar(baseStar0)
        baseStar1 = lsst.jointcal.star.BaseStar(1, 2, 3, 4)
        self.star1 = lsst.jointcal.star.MeasuredStar(baseStar1)
        # self.star.setInstFlux(instFlux)
        # self.star.setInstFluxErr(instFluxErr)

    def _test_offsetParams(self, delta, expect):
        self.mapping.offsetParams(delta)
        self.assertFloatsAlmostEqual(expect, self.mapping.getTransfo().getParameters())


class PhotometryMappingTestCase(PhotometryMappingTestBase, lsst.utils.tests.TestCase):
    def setUp(self):
        super(PhotometryMappingTestCase, self).setUp()
        transfo = lsst.jointcal.photometryTransfo.PhotometryTransfoSpatiallyInvariant()
        self.mapping = lsst.jointcal.photometryMappings.PhotometryMapping(transfo)

    def test_transformFlux(self):
        result = self.mapping.transformFlux(self.star0, self.instFlux)
        self.assertEqual(result, self.instFlux)

    def test_offsetParams(self):
        """Test offsetting; note that offsetParams offsets by `-delta`."""
        delta = np.zeros(1, dtype=float)
        self._test_offsetParams(delta, np.array(1.0))
        delta -= 1
        self._test_offsetParams(delta, np.array(2.0))

    def test_computeParameterDerivatives(self):
        """Test that the derivative of a spatially invariant transform is always the same."""
        result = self.mapping.computeParameterDerivatives(self.star0, self.instFlux)
        self.assertEqual(self.instFlux, result)
        result = self.mapping.computeParameterDerivatives(self.star1, self.instFlux)
        self.assertEqual(self.instFlux, result)
        transfo = lsst.jointcal.photometryTransfo.PhotometryTransfoSpatiallyInvariant(1000.0)
        mapping = lsst.jointcal.photometryMappings.PhotometryMapping(transfo)
        result = mapping.computeParameterDerivatives(self.star0, self.instFlux)
        self.assertEqual(self.instFlux, result)


@unittest.skip('djfkd;asdf')
class ChipVisitPhotometryMappingTestCase(PhotometryMappingTestBase, lsst.utils.tests.TestCase):
    def setUp(self):
        # np.random.seed(100)
        # self.x = (np.random.random(50) - 0.5) * 100
        # self.y = (np.random.random(50) - 0.5) * 100
        super(PhotometryMappingTestCase, self).setUp()
        self.bbox = lsst.afw.geom.Box2I(lsst.afw.geom.Point2I(-3, -3), lsst.afw.geom.Point2I(3, 3))
        self.degree = 2
        chipTransfo = lsst.jointcal.photometryTransfo.PhotometryTransfoSpatiallyInvariant()
        visitTransfo = lsst.jointcal.photometryTransfo.PhotometryTransfoSpatiallyInvariant()
        self.mappingInvariants = lsst.jointcal.photometryMappings.ChipVisitPhotometryMapping(chipTransfo,
                                                                                             visitTransfo)
        visitTransfo = lsst.jointcal.photometryMappings.PhotometryMappingChebyshev(self.degree, self.bbox)
        self.mappingConstrained = lsst.jointcal.photometryMappings.ChipVisitPhotometryMapping(chipTransfo,
                                                                                              visitTransfo)

    def test_transformFlux(self):
        result = self.mapping.transformFlux(1, 2, self.instFlux)
        self.assertEqual(result, self.instFlux)

    def test_offsetParams(self):
        """Test offsetting; note that offsetParams offsets by `-delta`."""
        delta = np.zeros(self.mapping.getNpar(), dtype=float)
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

    def test_computeParameterDerivatives(self):
        result = self.mapping.computeParameterDerivatives(self.point[0], self.point[1], self.instFlux)
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
