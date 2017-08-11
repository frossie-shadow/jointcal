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

ConstrainedPhotometryModel::ConstrainedPhotometryModel(CcdImageList const &ccdImageList, int degree) {
    // First initialize all visit and ccd transfos, before we make the ccdImage mappings.
    for (auto const &ccdImage : ccdImageList) {
        auto visit = ccdImage->getVisit();
        auto chip = ccdImage->getCcdId();
        auto visitPair = _visitMap.find(visit);
        auto chipPair = _chipMap.find(chip);
        // If the chip is not in the map, add it, otherwise continue.
        if (chipPair == _chipMap.end()) {
            auto photoCalib = ccdImage->getPhotoCalib();
            // Use (fluxMag0)^-1 from the PhotoCalib as the default.
            auto transfo = std::make_shared<PhotometryTransfoSpatiallyInvariant>(
                    1.0 / photoCalib->getInstFluxMag0());
            _chipMap[chip] = std::unique_ptr<PhotometryMapping>(new PhotometryMapping(transfo));
        }
        // If the visit is not in the map, add it, otherwise continue.
        if (visitPair == _visitMap.end()) {
            // if (_visitMap.size() == 0) {
            //     _visitMap[visit] =
            //             std::unique_ptr<SimpleGtransfoMapping>(new
            //             SimpleGtransfoMapping(GtransfoIdentity()));
            // } else {
            auto transfo = std::make_shared<PhotometryTransfoChebyshev>(degree);
            _visitMap[visit] = std::unique_ptr<PhotometryMapping>(new PhotometryMapping(transfo));
            // }
        }
    }
    // Now create the ccdImage mappings, which are combinations of the chip/visit mappings above.
    _myMap.reserve(ccdImageList.size());  // we know how big it will be, so pre-allocate space.
    for (auto const &ccdImage : ccdImageList) {
        auto visit = ccdImage->getVisit();
        auto chip = ccdImage->getCcdId();
        _myMap.emplace(ccdImage->getHashKey(),
                       std::unique_ptr<ChipVisitPhotometryMapping>(new ChipVisitPhotometryMapping(
                               _chipMap[chip]->getTransfo(), _visitMap[visit]->getTransfo())));
    }
    LOGLS_INFO(_log, "Constructor got " << _chipMap.size() << " chip mappings and " << _visitMap.size()
                                        << " visit mappings.");
    LOGLS_DEBUG(_log, "CcdImage map has " << _myMap.size() << " mappings, with " << _myMap.bucket_count()
                                          << " buckets and a load factor of " << _myMap.load_factor());
}

unsigned ConstrainedPhotometryModel::assignIndices(std::string const &whatToFit, unsigned firstIndex) {
    // TODO: currently ignoring whatToFit. Should we even care? Maybe it helps to initialize with fitting just
    // visits or chips first before fitting both simultaneously?
    unsigned index = firstIndex;
    for (auto &chip : _chipMap) {
        chip.second->setIndex(index);
        index += chip.second->getNpar();
    }
    for (auto &visit : _visitMap) {
        visit.second->setIndex(index);
        index += visit.second->getNpar();
    }
    return index;
}

void ConstrainedPhotometryModel::offsetParams(Eigen::VectorXd const &delta) {
    for (auto &i : _chipMap) {
        auto mapping = i.second.get();
        mapping->offsetParams(delta.segment(mapping->getIndex(), mapping->getNpar()));
    }
    for (auto &i : _visitMap) {
        auto mapping = i.second.get();
        mapping->offsetParams(delta.segment(mapping->getIndex(), mapping->getNpar()));
    }
}

double ConstrainedPhotometryModel::transformFlux(CcdImage const &ccdImage, MeasuredStar const &measuredStar,
                                                 double instFlux) const {
    auto mapping = this->findMapping(ccdImage, "transformFlux");
    return mapping->transformFlux(measuredStar, instFlux);
}

void ConstrainedPhotometryModel::getMappingIndices(CcdImage const &ccdImage, std::vector<unsigned> &indices) {
    auto mapping = this->findMapping(ccdImage, "getMappingIndices");
    if (indices.size() < mapping->getNpar()) indices.resize(mapping->getNpar());
    indices[0] = mapping->getIndex();
    // TODO: I think I need a for loop here, from the above value to that +mapping->getNpar()?
}

void ConstrainedPhotometryModel::computeParameterDerivatives(MeasuredStar const &measuredStar,
                                                             CcdImage const &ccdImage,
                                                             Eigen::VectorXd &derivatives) {}

PhotometryMappingBase *ConstrainedPhotometryModel::findMapping(CcdImage const &ccdImage,
                                                               std::string name) const {
    auto i = _myMap.find(ccdImage.getHashKey());
    if (i == _myMap.end())
        throw LSST_EXCEPT(
                pex::exceptions::InvalidParameterError,
                "ConstrainedPhotometryModel::" + name + ", cannot find CcdImage " + ccdImage.getName());
    return i->second.get();
}

}  // namespace jointcal
}  // namespace lsst
