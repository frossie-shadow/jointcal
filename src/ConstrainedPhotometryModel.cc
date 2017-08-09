#include <map>

#include "lsst/log/Log.h"

#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/ConstrainedPhotometryModel.h"
#include "lsst/jointcal/PhotometryMapping.h"

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.ConstrainedPhotometryModel");
}

namespace lsst {
namespace jointcal {

unsigned ConstrainedPhotometryModel::assignIndices(std::string const &whatToFit, unsigned firstIndex) {
    return 0;
}

void ConstrainedPhotometryModel::offsetParams(Eigen::VectorXd const &delta) {
    for (auto &i : _myMap) {
        auto mapping = i.second.get();
        mapping->offsetParams(delta.segment(mapping->getIndex(), mapping->getNpar()));
    }
}

double ConstrainedPhotometryModel::transformFlux(CcdImage const &ccdImage, MeasuredStar const &star,
                                                 double instFlux) const {
    return 0;
}

void ConstrainedPhotometryModel::getMappingIndices(CcdImage const &ccdImage, std::vector<unsigned> &indices) {
}

void ConstrainedPhotometryModel::computeParameterDerivatives(MeasuredStar const &measuredStar,
                                                             CcdImage const &ccdImage,
                                                             Eigen::VectorXd &derivatives) {}
}  // namespace jointcal
}  // namespace lsst
