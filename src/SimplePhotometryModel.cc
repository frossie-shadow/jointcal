#include <iostream>

#include "lsst/log/Log.h"
#include "lsst/jointcal/PhotometryMapping.h"
#include "lsst/jointcal/PhotometryTransfo.h"
#include "lsst/jointcal/SimplePhotometryModel.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/MeasuredStar.h"

namespace {
LOG_LOGGER _log = LOG_GET("jointcal.SimplePhotometryModel");
}

namespace lsst {
namespace jointcal {

SimplePhotometryModel::SimplePhotometryModel(CcdImageList const &ccdImageList) {
    for (auto const &ccdImage : ccdImageList) {
        auto photoCalib = ccdImage->getPhotoCalib();
        // Use (fluxMag0)^-1 from the PhotoCalib as the default.
        _myMap[ccdImage.get()] = std::unique_ptr<PhotometryMapping>(new PhotometryMapping(
                PhotometryTransfoSpatiallyInvariant(1.0 / photoCalib->getInstFluxMag0())));
    }
    LOGLS_INFO(_log, "SimplePhotometryModel got " << _myMap.size() << " ccdImage mappings.");
}

unsigned SimplePhotometryModel::assignIndices(std::string const &whatToFit, unsigned firstIndex) {
    unsigned ipar = firstIndex;
    for (auto const &i : _myMap) {
        auto mapping = i.second.get();
        mapping->setIndex(ipar);
        ipar += mapping->getNpar();
    }
    return ipar;
}

void SimplePhotometryModel::offsetParams(Eigen::VectorXd const &delta) {
    for (auto &i : _myMap) {
        auto mapping = i.second.get();
        mapping->offsetParams(delta.segment(mapping->getIndex(), mapping->getNpar()));
    }
}

double SimplePhotometryModel::transformFlux(CcdImage const &ccdImage, MeasuredStar const &star,
                                            double instFlux) const {
    auto mapping = this->findMapping(ccdImage, "transformFlux");
    return mapping->transformFlux(star.x, star.y, instFlux);
}

void SimplePhotometryModel::getMappingIndices(CcdImage const &ccdImage, std::vector<unsigned> &indices) {
    auto mapping = this->findMapping(ccdImage, "getMappingIndices");
    if (indices.size() < mapping->getNpar()) indices.resize(mapping->getNpar());
    indices[0] = mapping->getIndex();
}

void SimplePhotometryModel::computeParameterDerivatives(MeasuredStar const &measuredStar,
                                                        CcdImage const &ccdImage,
                                                        Eigen::VectorXd &derivatives) {
    auto mapping = this->findMapping(ccdImage, "computeParameterDerivatives");
    mapping->computeParameterDerivatives(measuredStar.x, measuredStar.y, measuredStar.getInstFlux(),
                                         derivatives);
}

PhotometryMapping *SimplePhotometryModel::findMapping(CcdImage const &ccdImage, std::string name) const {
    auto i = _myMap.find(&ccdImage);
    if (i == _myMap.end())
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterError,
                          "SimplePolyModel::" + name + ", cannot find CcdImage " + ccdImage.getName());
    return i->second.get();
}

}  // namespace jointcal
}  // namespace lsst
