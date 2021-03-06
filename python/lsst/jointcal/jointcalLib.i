// -*- lsst-c++ -*-
%define jointcalLib_DOCSTRING
"
Python interface to lsst::jointcal classes
"
%enddef

%feature("autodoc", "1");

%module(package="lsst.jointcal", docstring=jointcalLib_DOCSTRING) jointcalLib

// Suppress swig complaints; see afw/image/imageLib.i for more
#pragma SWIG nowarn=314                 // print is a python keyword (--> _print)
#pragma SWIG nowarn=362                 // operator=  ignored

%{
#include <exception>
#include <list>
#include "lsst/afw/table.h"
#include "lsst/afw/detection.h"
#include "lsst/afw/cameraGeom.h"
#include "lsst/afw/math.h"
#include "lsst/afw/image.h"
#include "lsst/afw/coord.h"
#include "lsst/afw/geom.h"
#include "lsst/afw/image/VisitInfo.h"
#include "lsst/afw/geom/Box.h"
#include "lsst/daf/base/PropertySet.h"
#include "lsst/jointcal/BaseStar.h"
#include "lsst/jointcal/StarList.h"
#include "lsst/jointcal/FittedStar.h"
#include "lsst/jointcal/RefStar.h"
#include "lsst/jointcal/Jointcal.h"
#include "lsst/jointcal/CcdImage.h"
#include "lsst/jointcal/AstromFit.h"
#include "lsst/jointcal/Associations.h"
#include "lsst/jointcal/DistortionModel.h"
#include "lsst/jointcal/SimplePolyModel.h"
#include "lsst/jointcal/Projectionhandler.h"
#include "lsst/jointcal/CountedRef.h"
#include "lsst/jointcal/SipToGtransfo.h"
#include "lsst/jointcal/Gtransfo.h"
#include "lsst/jointcal/FatPoint.h"
#include "lsst/jointcal/Point.h"
#include "lsst/jointcal/ConstrainedPolyModel.h"
#include "lsst/jointcal/PhotomFit.h"
#include "lsst/jointcal/SimplePhotomModel.h"
%}

%include "lsst/p_lsstSwig.i"
%initializeNumPy(jointcal)
%{
#include "ndarray/swig.h"
#include "ndarray/converter/eigen.h"
%}

%import "lsst/afw/table/tableLib.i"
%import "lsst/afw/table/io/ioLib.i"
%import "lsst/afw/image/imageLib.i"

%import "lsst/pex/config.h"

%shared_ptr(lsst::daf::base::PropertySet);
%shared_ptr(lsst::jointcal::JointcalControl);

%template(SourceList) std::vector<lsst::afw::table::SortedCatalogT< lsst::afw::table::SourceRecord> >;
%template(WcsList) std::vector<std::shared_ptr<lsst::afw::image::TanWcs> >;
%template(PropertySetList) std::vector<std::shared_ptr<lsst::daf::base::PropertySet> >;
%template(CalibList) std::vector<std::shared_ptr< lsst::afw::image::Calib > >;
%template(BboxList) std::vector<lsst::afw::geom::Box2I>;
%template(StringList) std::vector<std::string>;
%template(IntList) std::vector<int>;

%include "lsst/jointcal/Jointcal.h"
namespace lsst {
namespace jointcal {
class DistortionModel;
class TripletList;
class CcdImage;
class CcdImageList;
class MeasuredStar;
class FittedStarList;
class RefStarList;
class MeasuredStarList;
class BaseStar;
class RefStar;
class StarList;
class FittedStar;
class Point;
}}
%include "lsst/jointcal/Chi2.h"
%include "lsst/jointcal/AstromFit.h"
%include "lsst/jointcal/Associations.h"
namespace lsst {
namespace jointcal {
  class Gtransfo;
  class Mapping;
class SimplePolyMapping;
 class ProjectionHandler;
}}
%include "lsst/jointcal/DistortionModel.h"
namespace lsst {
namespace jointcal {
  class TanSipPix2RaDec;
}}
%include "lsst/jointcal/Projectionhandler.h"

namespace lsst {
namespace jointcal {
  class MeasuredStarList;
  class Frame;
}}

%shared_ptr(lsst::jointcal::CcdImage);
// TODO: NOTE: this may not be needed? Unclear.
// %template(CcdImageRef) std::shared_ptr<lsst::jointcal::CcdImage>;
%template(CcdImageRefList) std::list<std::shared_ptr<lsst::jointcal::CcdImage> >;
%template(TanSipPix2RaDecRef) std::shared_ptr<lsst::jointcal::TanSipPix2RaDec>;

// TODO: SWIG doesn't like boost::intrusive_ptr so we can't wrap these.
// Once DM-4043 is solved, we should clean this up to better wrap the lists.
// %shared_ptr(lsst::jointcal::BaseStar);
// %shared_ptr(lsst::jointcal::RefStar);
// %shared_ptr(lsst::jointcal::FittedStar);
// %template(RefStarRef) std::shared_ptr<lsst::jointcal::RefStar>;
// %intrusive_ptr(lsst::jointcal::RefStar);
// %template(StarListOfRefStar) std::list<std::intrusive_ptr<lsst::jointcal::RefStar> >;
%include "lsst/jointcal/StarList.h"
// %template(StarListMonkey) lsst::jointcal::StarList<lsst::jointcal::RefStar>;
// typedef StarList std::list<lsst::jointcal::CountedRef>;
// %template(RefStarRefList) std::list<lsst::jointcal::RefStar>;

%shared_ptr(lsst::jointcal::TanSipPix2RaDec);
%include "lsst/jointcal/Point.h"
%include "lsst/jointcal/FatPoint.h"
%include "lsst/jointcal/Gtransfo.h"

%include "lsst/jointcal/CcdImage.h"
%include "lsst/jointcal/BaseStar.h"
%include "lsst/jointcal/FittedStar.h"
%include "lsst/jointcal/RefStar.h"
%include "lsst/jointcal/SimplePolyModel.h"
%include "lsst/jointcal/ConstrainedPolyModel.h"

%include "lsst/jointcal/SipToGtransfo.h"
%include "lsst/jointcal/PhotomModel.h"
%include "lsst/jointcal/PhotomFit.h"
%include "lsst/jointcal/SimplePhotomModel.h"

