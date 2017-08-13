/*
 * LSST Data Management System
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 * See the COPYRIGHT file
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the LSST License Statement and
 * the GNU General Public License along with this program.  If not,
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */

#include "pybind11/pybind11.h"
#include "numpy/arrayobject.h"
#include "ndarray/pybind11.h"
#include "ndarray/eigen.h"
#include "Eigen/Core"

#include "lsst/jointcal/PhotometryTransfo.h"
#include "lsst/jointcal/PhotometryMapping.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace jointcal {
namespace {

void declarePhotometryMappingBase(py::module &mod) {
    py::class_<PhotometryMappingBase, std::shared_ptr<PhotometryMappingBase>> cls(mod,
                                                                                  "PhotometryMappingBase");

    cls.def("getNpar", &PhotometryMappingBase::getNpar);
    cls.def("transformFlux", &PhotometryMappingBase::transformFlux);
    cls.def("computeParameterDerivatives", &PhotometryMappingBase::computeParameterDerivatives);
    cls.def("offsetParams", &PhotometryMappingBase::offsetParams);
    cls.def("getParameters", &PhotometryMappingBase::getParameters);
    cls.def("getMappingIndices", &PhotometryMappingBase::getMappingIndices);

    cls.def("getIndex", &PhotometryMappingBase::getIndex);
    cls.def("setIndex", &PhotometryMappingBase::setIndex);
}

void declarePhotometryMapping(py::module &mod) {
    py::class_<PhotometryMapping, std::shared_ptr<PhotometryMapping>, PhotometryMappingBase> cls(
            mod, "PhotometryMapping");
    cls.def(py::init<std::shared_ptr<PhotometryTransfo>>(), "transfo"_a);
}

void declareChipVisitPhotometryMapping(py::module &mod) {
    py::class_<ChipVisitPhotometryMapping, std::shared_ptr<ChipVisitPhotometryMapping>, PhotometryMappingBase>
            cls(mod, "ChipVisitPhotometryMapping");

    cls.def(py::init<std::shared_ptr<PhotometryMapping>, std::shared_ptr<PhotometryMapping>>(),
            "chipMapping"_a, "visitMapping"_a);

    // cls.def("getTransfo1", &ChipVisitPhotometryMapping::getTransfo1,
    //         py::return_value_policy::reference_internal);
    // cls.def_property_readonly("transfo1", &ChipVisitPhotometryMapping::getTransfo1,
    //                           py::return_value_policy::reference_internal);

    // cls.def("getTransfo2", &ChipVisitPhotometryMapping::getTransfo2,
    //         py::return_value_policy::reference_internal);
    // cls.def_property_readonly("transfo2", &ChipVisitPhotometryMapping::getTransfo2,
    //                           py::return_value_policy::reference_internal);
}

PYBIND11_PLUGIN(PhotometryMappings) {
    py::module::import("lsst.jointcal.star");
    py::module::import("lsst.jointcal.photometryTransfo");
    py::module mod("PhotometryMappings");

    if (_import_array() < 0) {
        PyErr_SetString(PyExc_ImportError, "numpy.core.multiarray failed to import");
        return nullptr;
    }

    declarePhotometryMappingBase(mod);
    declareChipVisitPhotometryMapping(mod);
    declarePhotometryMapping(mod);

    return mod.ptr();
}
}  // namespace
}  // namespace jointcal
}  // namespace lsst
