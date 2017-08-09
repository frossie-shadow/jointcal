// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_PHOTOMETRY_TRANSFO_H
#define LSST_JOINTCAL_PHOTOMETRY_TRANSFO_H

#include <iostream>
#include <sstream>
#include <memory>

#include "Eigen/Core"
#include "ndarray.h"

#include "lsst/afw/geom/AffineTransform.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/afw/geom/Point.h"
#include "lsst/jointcal/Point.h"

namespace lsst {
namespace jointcal {

class Point;

class PhotometryTransfoSpatiallyInvariant;

/*
 * A photometric transform, defined as a scale factor of the input calibration.
 *
 *     initialCalibFlux (Maggies) * transfo(x,y) -> correctedFlux (Maggies)
 *
 * @todo Eventually will be defined the same as PhotoCalib:
 *     instFlux (ADU) -> flux (maggies)
 *
 * @seealso lsst::afw::image::PhotoCalib
 */
class PhotometryTransfo {
public:
    /// Apply the transform to instFlux at (x,y), put result in flux
    virtual double apply(double x, double y, double instFlux) const = 0;

    double apply(Point const &in, double instFlux) const { return apply(in.x, in.y, instFlux); }

    /// dumps the transfo coefficients to stream.
    virtual void dump(std::ostream &stream = std::cout) const = 0;

    /// Return a string describing this transfo. For the pybind11/python layer.
    std::string __str__() const {
        std::stringstream s;
        dump(s);
        return s.str();
    }

    /// Return the number of parameters (used to compute chisq)
    virtual int getNpar() const { return 0; }

    /**
     * Offset the parameters by some (negative) amount during fitting.
     *
     * Equivalent to `flatten(parameters) -= delta`
     *
     * Ordering of delta is the same as the ordering of the derivatives returned from `parameterDerivatives`.
     */
    virtual void offsetParams(Eigen::VectorXd const &delta) = 0;

    /// return a copy (allocated by new) of the transformation.
    virtual std::unique_ptr<PhotometryTransfo> clone() const = 0;

    /**
     * Compute the derivatives with respect to the parameters (i.e. the coefficients).
     *
     * @param[in]  x        The x coordinate to compute at (in the appropriate units for this transfo).
     * @param[in]  y        The y coordinate to compute at (in the appropriate units for this transfo).
     * @param[in]  instFlux     The instrument flux to compute the derivative at.
     * @param[out] derivatives  The computed derivatives, in the same order as the deltas in offsetParams.
     */
    virtual void parameterDerivatives(double x, double y, double instFlux,
                                      Eigen::VectorXd &derivatives) const = 0;

    void computeDerivative(Point const &where, PhotometryTransfoSpatiallyInvariant &derivative,
                           const double step = 0.01) const;
};

/*
 * Photometric offset independent of position.
 *
 * initialCalibFlux (Maggies) * SpatiallyInvariantTransfo -> correctedFlux (Maggies)
 *
 * @todo Eventually to be defined as:
 *     instFlux / value = flux
 */
class PhotometryTransfoSpatiallyInvariant : public PhotometryTransfo {
public:
    PhotometryTransfoSpatiallyInvariant(double value = 1) : _value(value) {}

    /// @copydoc PhotometryTransfo::apply
    double apply(double x, double y, double instFlux) const override { return instFlux * _value; }

    /// @copydoc PhotometryTransfo::dump
    void dump(std::ostream &stream = std::cout) const override { stream << _value; }

    /// @copydoc PhotometryTransfo::getNpar
    int getNpar() const override { return 1; }

    /// @copydoc PhotometryTransfo::offsetParams
    void offsetParams(Eigen::VectorXd const &delta) override { _value -= delta[0]; };

    /// @copydoc PhotometryTransfo::clone
    std::unique_ptr<PhotometryTransfo> clone() const override {
        return std::unique_ptr<PhotometryTransfo>(new PhotometryTransfoSpatiallyInvariant(_value));
    }

    /// The spatial derivative of a constant zeropoint is 1.
    void computeDerivative(Point const &where, PhotometryTransfoSpatiallyInvariant &derivative,
                           const double step = 0.01) const {
        derivative.setValue(1);
    }

    /// @copydoc PhotometryTransfo::parameterDerivatives
    void parameterDerivatives(double x, double y, double instFlux,
                              Eigen::VectorXd &derivatives) const override {
        // the derivative of a spatially constant transfo w.r.t. that value is just the instFlux.
        derivatives[0] = instFlux;
    }

protected:
    void setValue(double value) { _value = value; }

    friend class PhotometryTransfo;

private:
    /// value of this transform at all locations.
    double _value;
};

/**
 * nth-order 2d Chebyshev photometry transfo.
 *
 * The 2-d Chebyshev polynomial used here is defined as:
 *
 *  @f[
 *  f(x,y) = \sum_i \sum_j a_{i,j} T_i(x) T_j(y)
 *  @f]
 *
 * where @f$T_n(x)@f$ is the n-th order Chebyshev polynomial of @f$x@f$ and
 * @f$a_{i,j}@f$ is the corresponding coefficient of the (i,j) polynomial term.
 *
 * Note that the polynomial order=n means that the highest terms will be of the form:
 *   @f[
 *   a_{0,n}*x^n*y^0, a_{n-1,1}*x^(n-1)*y^1, ..., a_{1,n-1}*x^1*y^(n-1), a_{n,0}*x^0*y^n
 *   @f]
 */
class PhotometryTransfoChebyshev : public PhotometryTransfo {
public:
    /**
     * Create a Chebyshev transfo with terms up to order in (x*y)
     *
     * @param[in]  order  The maximum order in (x*y).
     */
    PhotometryTransfoChebyshev(size_t order, afw::geom::Box2I const &bbox);

    /// @copydoc PhotometryTransfo::apply
    double apply(double x, double y, double instFlux) const override;

    /// @copydoc PhotometryTransfo::dump
    void dump(std::ostream &stream = std::cout) const override {
        stream << "PhotometryTransfoChebyshev (" << _coefficients.getShape() << " coefficients in y,x)";
    }

    /// @copydoc PhotometryTransfo::getNpar
    int getNpar() const override { return _nParameters; }

    /// @copydoc PhotometryTransfo::offsetParams
    void offsetParams(Eigen::VectorXd const &delta) override;

    /// @copydoc PhotometryTransfo::clone
    std::unique_ptr<PhotometryTransfo> clone() const override {
        return nullptr;
        // return std::unique_ptr<PhotometryTransfo>(new PhotometryTransfoChebyshev());
    }

    /// @copydoc PhotometryTransfo::parameterDerivatives
    void parameterDerivatives(double x, double y, double instFlux,
                              Eigen::VectorXd &derivatives) const override;

    /// Get a copy of the coefficients of the polynomials, as a 2d array (NOTE: order is [y][x])
    ndarray::Array<double, 2, 2> getCoefficients() { return ndarray::copy(_coefficients); }

private:
    afw::geom::AffineTransform _toChebyshevRange;  // maps points from the bbox to [-1,1]x[-1,1]

    ndarray::Array<double, 2, 2> _coefficients;  // shape=(order+1, order+1)
    ndarray::Size _order;
    ndarray::Size _nParameters;
};

}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_PHOTOMETRY_TRANSFO_H
