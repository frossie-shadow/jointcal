// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_PHOTOMETRY_MAPPING_H
#define LSST_JOINTCAL_PHOTOMETRY_MAPPING_H

#include <memory>

#include "lsst/afw/image/PhotoCalib.h"

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/MeasuredStar.h"
#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/PhotometryTransfo.h"

namespace lsst {
namespace jointcal {

class PhotometryMapping {
public:
    explicit PhotometryMapping(PhotometryTransfo const &_transfo) : index(-1), transfo(_transfo.clone()) {}

    /// No copy or move: there is only ever one instance of a given mapping (i.e. per ccd+visit)
    PhotometryMapping(PhotometryMapping const &) = delete;
    PhotometryMapping(PhotometryMapping &&) = delete;
    PhotometryMapping &operator=(PhotometryMapping const &) = delete;
    PhotometryMapping &operator=(PhotometryMapping &&) = delete;

    /// Number of total parameters in this mapping
    unsigned getNpar() const { return transfo->getNpar(); }

    /*
     * Sets how this set of parameters (of length getNpar()) map into the "grand" fit.
     * Expects that indices has enough space reserved.
     */
    void setMappingIndices(std::vector<unsigned> &indices) const {
        indices.reserve(getNpar());
        for (unsigned k = 0; k < getNpar(); ++k) {
            indices[k] = index + k;
        }
    }

    /**
     * Return the on-sky transformed flux at (x,y).
     *
     * @param[in]  x         The x coordinate to compute at (in the appropriate units for this transfo).
     * @param[in]  y         The y coordinate to compute at (in the appropriate units for this transfo).
     * @param[in]  instFlux  The instrument flux to transform.
     *
     * @return     The on-sky flux transformed from instFlux at measuredStar's position.
     */
    double transformFlux(double x, double y, double instFlux) const { return transfo->apply(x, y, instFlux); }

    /**
     * Compute the derivatives with respect to the parameters (i.e. the coefficients).
     *
     * @param[in]  x        The x coordinate to compute at (in the appropriate units for this transfo).
     * @param[in]  y        The y coordinate to compute at (in the appropriate units for this transfo).
     * @param[in]  instFlux     The instrument flux to compute the derivative at.
     * @param[out] derivatives  The computed derivatives, in the same order as the deltas in offsetParams.
     */
    void computeParameterDerivatives(double x, double y, double instFlux,
                                     Eigen::VectorXd &derivatives) const {
        transfo->computeParameterDerivatives(x, y, instFlux, derivatives);
    }

    void offsetParams(Eigen::VectorXd const &delta) { transfo->offsetParams(delta); }

    /// Get the index of this mapping in the grand fit.
    unsigned getIndex() { return index; }

    /// Set the index of this mapping in the grand fit.
    void setIndex(unsigned i) { index = i; }

    /**
     * Return this mapping represented as a PhotoCalib.
     */
    std::unique_ptr<afw::image::PhotoCalib> toPhotoCalib() const {
        return std::unique_ptr<afw::image::PhotoCalib>(new afw::image::PhotoCalib(0.0));
    }

    std::shared_ptr<PhotometryTransfo> getTransfo() { return transfo; }

protected:
    // Start index of this mapping in the "grand" fit
    unsigned index;

    // the actual transformation to be fit
    std::shared_ptr<PhotometryTransfo> transfo;
};

/**
 * Class representing a two-level photometric transform: one for the ccd and one for the visit.
 */
class TwoTransfoPhotometryMapping : public PhotometryTransfo {
    std::shared_ptr<afw::image::PhotoCalib> toPhotoCalib() const {
        return std::unique_ptr<afw::image::PhotoCalib>(new afw::image::PhotoCalib(0.0));
    }

protected:
    // Start index of this mapping in the "grand" fit
    unsigned index;

    // the actual transformation to be fit
    std::shared_ptr<PhotometryTransfo> transfo;
};

}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_PHOTOMETRY_MAPPING_H
