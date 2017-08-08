#include "ndarray.h"
#include "Eigen/Core"

#include "lsst/afw/math/detail/TrapezoidalPacker.h"

#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/PhotometryTransfo.h"

namespace lsst {
namespace jointcal {

void PhotometryTransfo::computeDerivative(const Point &where, PhotometryTransfoSpatiallyInvariant &derivative,
                                          const double step) const {
    // double result = 0;
    double flux1 = apply(where, 1);
    derivative.setValue((1 - flux1) / step);
}

/**
 * PhotometryTransfoChebyshev
 */

// ------------------ PhotometryTransfoChebyshev helpers ---------------------------------------------------

namespace {

// To evaluate a 1-d Chebyshev function without needing to have workspace, we use the
// Clenshaw algorith, which is like going through the recurrence relation in reverse.
// The CoeffGetter argument g is something that behaves like an array, providing access
// to the coefficients.
template <typename CoeffGetter>
double evaluateFunction1d(CoeffGetter g, double x, int size) {
    double b_kp2 = 0.0, b_kp1 = 0.0;
    for (int k = (size - 1); k > 0; --k) {
        double b_k = g[k] + 2 * x * b_kp1 - b_kp2;
        b_kp2 = b_kp1;
        b_kp1 = b_k;
    }
    return g[0] + x * b_kp1 - b_kp2;
}

// This class imitates a 1-d array, by running evaluateFunction1d on a nested dimension;
// this lets us reuse the logic in evaluateFunction1d for both dimensions.  Essentially,
// we run evaluateFunction1d on a column of coefficients to evaluate T_i(x), then pass
// the result of that to evaluateFunction1d with the results as the "coefficients" associated
// with the T_j(y) functions.
struct RecursionArrayImitator {
    double operator[](int i) const {
        return evaluateFunction1d(coefficients[i], x, coefficients.getSize<1>());
    }

    RecursionArrayImitator(ndarray::Array<double const, 2, 2> const &coefficients_, double x_)
            : coefficients(coefficients_), x(x_) {}

    ndarray::Array<double const, 2, 2> coefficients;
    double x;
};

// Compute an affine transform that maps an arbitrary box to [-1,1]x[-1,1]
afw::geom::AffineTransform makeChebyshevRangeTransform(afw::geom::Box2D const bbox) {
    return afw::geom::AffineTransform(
            afw::geom::LinearTransform::makeScaling(2.0 / bbox.getWidth(), 2.0 / bbox.getHeight()),
            afw::geom::Extent2D(-(2.0 * bbox.getCenterX()) / bbox.getWidth(),
                                -(2.0 * bbox.getCenterY()) / bbox.getHeight()));
}

// Initialize a "unit" Chebyshev
ndarray::Array<double, 2, 2> _identityChebyshev(size_t order) {
    ndarray::Array<double, 2, 2> coeffs = ndarray::allocate(ndarray::makeVector(order, order));
    coeffs.deep() = 0.0;
    coeffs[0][0] = 1;
    return coeffs;
}
}  // namespace

PhotometryTransfoChebyshev::PhotometryTransfoChebyshev(size_t order)
        : _toChebyshevRange(makeChebyshevRangeTransform(afw::geom::Box2D())),
          _coefficients(_identityChebyshev(order)) {}

double PhotometryTransfoChebyshev::apply(double x, double y, double instFlux) const {
    return instFlux *
           evaluateFunction1d(RecursionArrayImitator(_coefficients, x), y, _coefficients.getSize<0>());
}

void PhotometryTransfoChebyshev::offsetParams(Eigen::VectorXd const &delta) {
    // NOTE: the indexing in this method and parameterDerivatives must be kept consistent
    Eigen::VectorXd::Index k = 0;
    for (ndarray::Size j = 0; j < _coefficients.getSize<0>(); ++j) {
        for (ndarray::Size i = 0; i < _coefficients.getSize<1>(); ++i, ++k) {
            _coefficients[j][i] -= delta[k];
        }
    }
}

void PhotometryTransfoChebyshev::parameterDerivatives(double x, double y, double instFlux,
                                                      Eigen::VectorXd &derivatives) const override {
    // NOTE: the indexing in this method and offsetParams must be kept consistent
    Eigen::VectorXd::Index k = 0;
    for (ndarray::Size j = 0; j < _coefficients.getSize<0>(); ++j) {
        for (ndarray::Size i = 0; i < _coefficients.getSize<1>(); ++i, ++k) {
            derivatives[k] = ;
        }
    }
}
}  // namespace jointcal
}  // namespace lsst
