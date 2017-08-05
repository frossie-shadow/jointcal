#include "ndarray.h"

#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/PhotometryTransfo.h"

namespace lsst {
namespace jointcal {

// = PhotometryTransfoChebyshev

void PhotometryTransfo::computeDerivative(const Point &where, PhotometryTransfoSpatiallyInvariant &derivative,
                                          const double step) const {
    // double result = 0;
    double flux1 = apply(where, 1);
    derivative.setValue((1 - flux1) / step);
}

// = PhotometryTransfoChebyshev
PhotometryTransfoChebyshev::PhotometryTransfoChebyshev(dize_t order) {
    coeffs = ndarray::allocate(makeVector(2, order));
    afw::math::ChebyshevBoundedField(afw::geom::Box2I(), coeffs);
}

void PhotometryTransfoChebyshev::offsetParams(Eigen::VectorXd const &delta) override {
    for (size_t i = 0; i < delta.size(); i++) {
        ndarray::flatten<1>(_coefficients)[i] -= delta[i];
    }
}

}  // namespace jointcal
}  // namespace lsst
