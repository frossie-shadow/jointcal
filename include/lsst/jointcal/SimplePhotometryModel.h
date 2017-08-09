// -*- LSST-C++ -*-
#ifndef LSST_JOINTCAL_SIMPLE_PHOTOMETRY_MODEL_H
#define LSST_JOINTCAL_SIMPLE_PHOTOMETRY_MODEL_H

#include <map>

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/Eigenstuff.h"
#include "lsst/jointcal/PhotometryMapping.h"
#include "lsst/jointcal/PhotometryModel.h"
#include "lsst/jointcal/Point.h"

namespace lsst {
namespace jointcal {

class CcdImage;
class Point;

//! Photometric response model which has a single photometric factor per CcdImage.
/*! It considers a full exposure as reference. */
class SimplePhotometryModel : public PhotometryModel {
public:
    SimplePhotometryModel(CcdImageList const &ccdImageList);

    /// No copy or move: there is only ever one instance of a given model (i.e. per ccd+visit)
    SimplePhotometryModel(SimplePhotometryModel const &) = delete;
    SimplePhotometryModel(SimplePhotometryModel &&) = delete;
    SimplePhotometryModel &operator=(SimplePhotometryModel const &) = delete;
    SimplePhotometryModel &operator=(SimplePhotometryModel &&) = delete;

    /**
     * Assign indices to parameters involved in mappings, starting at firstIndex.
     *
     * @param[in]  whatToFit   Ignored.
     * @param[in]  firstIndex  Index to start assigning at.
     *
     * @return     The highest assigned index.
     */
    unsigned assignIndices(std::string const &whatToFit, unsigned firstIndex) override;

    /// @copydoc PhotometryModel::offsetParams
    void offsetParams(Eigen::VectorXd const &delta) override;

    /// @copydoc PhotometryModel::transformFlux
    double transformFlux(CcdImage const &ccdImage, MeasuredStar const &star, double instFlux) const override;

    /// @copydoc PhotometryModel::getMappingIndices
    void getMappingIndices(CcdImage const &ccdImage, std::vector<unsigned> &indices) override;

    /// @copydoc PhotometryModel::computeParameterDerivatives
    void computeParameterDerivatives(MeasuredStar const &measuredStar, CcdImage const &ccdImage,
                                     Eigen::VectorXd &derivatives) override;

private:
    typedef std::map<CcdImage const *, std::unique_ptr<PhotometryMapping>> MapType;
    MapType _myMap;

    /// Return the mapping associated with this ccdImage. name is a descriptor for error messages.
    PhotometryMapping *findMapping(CcdImage const &ccdImage, std::string name) const override;
};

}  // namespace jointcal
}  // namespace lsst

#endif  // LSST_JOINTCAL_SIMPLE_PHOTOMETRY_MODEL_H
