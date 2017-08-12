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

/**
 * Relates transfo(s) to their position in the fitting matrix and allows interaction with the transfo(s).
 */
class PhotometryMappingBase {
public:
    PhotometryMappingBase() : index(-1) {}

    /// No copy or move: there is only ever one instance of a given mapping (i.e. per ccd+visit)
    PhotometryMappingBase(PhotometryMappingBase const &) = delete;
    PhotometryMappingBase(PhotometryMappingBase &&) = delete;
    PhotometryMappingBase &operator=(PhotometryMappingBase const &) = delete;
    PhotometryMappingBase &operator=(PhotometryMappingBase &&) = delete;

    /// Number of total parameters in this mapping
    virtual unsigned getNpar() const = 0;

    /**
     * Return the on-sky transformed flux at (x,y).
     *
     * @param[in]  measuredStar  The measured star position to transform.
     * @param[in]  instFlux      The instrument flux to transform.
     *
     * @return     The on-sky flux transformed from instFlux at measuredStar's
     *             position.
     */
    virtual double transformFlux(MeasuredStar const &measuredStar, double instFlux) const = 0;

    /**
     * Compute the derivatives with respect to the parameters (i.e. the coefficients).
     *
     * @param[in]  measuredStar The measured star position to transform.
     * @param[in]  instFlux     The instrument flux to compute the derivative at.
     * @param[out] derivatives  The computed derivatives, in the same order as the deltas in offsetParams.
     */
    virtual void computeParameterDerivatives(MeasuredStar const &measuredStar, double instFlux,
                                             Eigen::Ref<Eigen::VectorXd> derivatives) const = 0;

    /**
     * Offset the transfo parameters by delta.
     *
     * @param[in]   delta vector to offset transfo parameters. Same ordering as derivatives in
     *              computeParameterDerivatives();
     */
    virtual void offsetParams(Eigen::VectorXd const &delta) = 0;

    virtual Eigen::VectorXd getParameters() = 0;

    /**
     * Gets how this set of parameters (of length getNpar()) map into the "grand" fit.
     * Expects that indices has enough space reserved.
     */
    virtual void getMappingIndices(std::vector<unsigned> &indices) const = 0;

    /// Dump the contents of the transfos, for debugging.
    virtual void dump(std::ostream &stream = std::cout) const = 0;

    /// Get the index of this mapping in the grand fit.
    unsigned getIndex() { return index; }

    /// Set the index of this mapping in the grand fit.
    void setIndex(unsigned i) { index = i; }

protected:
    // Start index of this mapping in the "grand" fit
    unsigned index;
};

/**
 * A single-transfo mapping.
 */
class PhotometryMapping : public PhotometryMappingBase {
public:
    explicit PhotometryMapping(std::shared_ptr<PhotometryTransfo> transfo)
            : PhotometryMappingBase(), _transfo(std::move(transfo)) {}

    /// @copydoc PhotometryMappingBase::getNpar
    unsigned getNpar() const override { return _transfo->getNpar(); }

    /// @copydoc PhotometryMappingBase::transformFlux
    double transformFlux(MeasuredStar const &measuredStar, double instFlux) const override {
        return _transfo->apply(measuredStar.x, measuredStar.y, instFlux);
    }

    /// @copydoc PhotometryMappingBase::computeParameterDerivatives
    void computeParameterDerivatives(MeasuredStar const &measuredStar, double instFlux,
                                     Eigen::Ref<Eigen::VectorXd> derivatives) const override {
        _transfo->computeParameterDerivatives(measuredStar.x, measuredStar.y, instFlux, derivatives);
    }

    /// @copydoc PhotometryMappingBase::offsetParams
    void offsetParams(Eigen::VectorXd const &delta) override { _transfo->offsetParams(delta); }

    /// @copydoc PhotometryMappingBase::getParameters
    Eigen::VectorXd getParameters() override { return _transfo->getParameters(); }

    /// @copydoc PhotometryMappingBase::getMappingIndices
    void getMappingIndices(std::vector<unsigned> &indices) const override {
        indices.reserve(getNpar());
        for (unsigned k = 0; k < getNpar(); ++k) {
            indices[k] = index + k;
        }
    }

    /// @copydoc PhotometryMappingBase::dump
    void dump(std::ostream &stream = std::cout) const override {
        stream << "index: " << index << " params: ";
        _transfo->dump(stream);
    }

    std::shared_ptr<PhotometryTransfo> getTransfo() const { return _transfo; }

private:
    // the actual transformation to be fit
    std::shared_ptr<PhotometryTransfo> _transfo;
};

/**
 * A two-level photometric transform: one for the ccd and one for the visit.
 */
class ChipVisitPhotometryMapping : public PhotometryMappingBase {
public:
    ChipVisitPhotometryMapping(std::shared_ptr<PhotometryMapping> chipMapping,
                               std::shared_ptr<PhotometryMapping> visitMapping)
            : PhotometryMappingBase(),
              _chipMapping(std::move(chipMapping)),
              _visitMapping(std::move(visitMapping)) {}

    // std::shared_ptr<afw::image::PhotoCalib> toPhotoCalib() const {
    //     return std::shared_ptr<afw::image::PhotoCalib>(new afw::image::PhotoCalib(0.0));
    // }

    /// @copydoc PhotometryMappingBase::getNpar
    unsigned getNpar() const override { return _chipMapping->getNpar() + _visitMapping->getNpar(); }

    /// @copydoc PhotometryMappingBase::transformFlux
    double transformFlux(MeasuredStar const &measuredStar, double instFlux) const override {
        double tempFlux = _chipMapping->getTransfo()->apply(measuredStar.x, measuredStar.y, instFlux);
        return _visitMapping->getTransfo()->apply(measuredStar.getXFocal(), measuredStar.getYFocal(),
                                                  tempFlux);
    }

    /// @copydoc PhotometryMappingBase::computeParameterDerivatives
    void computeParameterDerivatives(MeasuredStar const &measuredStar, double instFlux,
                                     Eigen::Ref<Eigen::VectorXd> derivatives) const override;

    /// @copydoc PhotometryMappingBase::offsetParams
    void offsetParams(Eigen::VectorXd const &delta) override {
        _chipMapping->offsetParams(delta.segment(0, _chipMapping->getNpar()));
        _visitMapping->offsetParams(delta.segment(_chipMapping->getNpar(), _visitMapping->getNpar()));
    }

    /// @copydoc PhotometryMappingBase::getParameters
    Eigen::VectorXd getParameters() override {
        Eigen::VectorXd joined(getNpar());
        joined << _chipMapping->getParameters(), _visitMapping->getParameters();
        return joined;
    }

    /// @copydoc PhotometryMappingBase::getMappingIndices
    void getMappingIndices(std::vector<unsigned> &indices) const override;

    /// @copydoc PhotometryMappingBase::dump
    void dump(std::ostream &stream = std::cout) const override {
        stream << "index: " << index << " chipMapping: ";
        _chipMapping->dump(stream);
        stream << "visitMapping: ";
        _visitMapping->dump(stream);
    }

private:
    // the actual transformation to be fit
    std::shared_ptr<PhotometryMapping> _chipMapping;
    std::shared_ptr<PhotometryMapping> _visitMapping;
};

}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_PHOTOMETRY_MAPPING_H
